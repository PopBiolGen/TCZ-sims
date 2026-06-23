# dea-drying-probs.R
#
# Queries the Digital Earth Australia (DEA) Waterbodies API to derive
# empirical drying probabilities for waterpoints in the living.waters dataset.
#
# Input:  living.waters  (sf POINT, WGS84) — must be in session
# Output: living.waters.dea  (sf, same points + DEA-derived columns)
#         annual_metrics      (tibble, per-waterbody per-year dry-season summary)
#
# Methodology:
#   1. Query DEA WFS for waterbody polygons within the living.waters bbox.
#   2. Match each living.waters point to the nearest DEA polygon (≤ MAX_DIST_M).
#   3. Download the ~16-day time series CSV for each matched waterbody.
#   4. Filter to the dry season (April–October) to avoid monsoon cloud gaps.
#   5. A year is classified "dry" when mean dry-season pc_wet == 0 across all
#      valid observations (requires ≥ MIN_OBS_PER_YEAR valid obs).
#   6. p_dry_empirical = fraction of years classified dry (1986–present).
#
# Caveats:
#   - DEA Waterbodies only maps features > 2700 m² present > 10% of the time.
#     Small soaks below this threshold will be unmatched (p_dry = NA).
#   - Dense phreatophytic vegetation can mask water signal, potentially
#     inflating dry-year counts at spring sites.
#   - These are marginal probabilities; see annual_metrics for year-level data
#     to build conditional probabilities on rainfall.

library(sf)
library(httr2)
library(dplyr)
library(lubridate)

# ── Parameters ────────────────────────────────────────────────────────────────
MAX_DIST_M      <- 1000   # max distance (m) to nearest DEA polygon to be "matched"
MIN_OBS_PER_YEAR <- 3     # minimum valid dry-season observations to classify a year
DRY_SEASON_MONTHS <- 6:12 # April–October (avoids monsoon cloud cover)

# ── 1. Query DEA WFS for waterbody polygons ───────────────────────────────────
bbox <- st_bbox(living.waters)
buf  <- 0.05  # degrees buffer around bbox

wfs_url <- paste0(
  "https://geoserver.dea.ga.gov.au/geoserver/wfs?",
  "service=WFS&version=1.1.0&request=GetFeature",
  "&typeName=dea:DigitalEarthAustraliaWaterbodies_v3",
  "&bbox=", paste(bbox["xmin"] - buf, bbox["ymin"] - buf,
                  bbox["xmax"] + buf, bbox["ymax"] + buf,
                  "EPSG:4326", sep = ","),
  "&outputFormat=application/json",
  "&maxFeatures=500"
)

dea_wb <- request(wfs_url) |>
  req_perform() |>
  resp_body_string() |>
  st_read(quiet = TRUE)

message(nrow(dea_wb), " DEA waterbodies found in region")

# ── 2. Spatial matching ───────────────────────────────────────────────────────
# DEA polygons are in GDA94/Albers (EPSG:3577); transform living.waters to match
lw_albers    <- living.waters |> st_transform(3577)
nearest_idx  <- st_nearest_feature(lw_albers, dea_wb)
nearest_dist <- st_distance(lw_albers, dea_wb[nearest_idx, ], by_element = TRUE)

lw_matched <- lw_albers |>
  mutate(
    dea_uid      = dea_wb$uid[nearest_idx],
    dea_ts_url   = dea_wb$timeseries[nearest_idx],
    dea_area_m2  = dea_wb$area_m2[nearest_idx],
    dist_to_wb_m = as.numeric(nearest_dist),
    dea_matched  = dist_to_wb_m <= MAX_DIST_M
  )

message(sum(lw_matched$dea_matched),  " points matched within ", MAX_DIST_M, "m")
message(sum(!lw_matched$dea_matched), " points unmatched (likely below DEA detection threshold)")

# ── 3. Download time series CSVs ──────────────────────────────────────────────
# Deduplicate: multiple living.waters points may share a DEA waterbody
urls_to_fetch <- lw_matched |>
  st_drop_geometry() |>
  filter(dea_matched) |>
  distinct(dea_uid, dea_ts_url)

fetch_ts <- function(uid, csv_url) {
  tryCatch({
    df        <- read.csv(url(csv_url), stringsAsFactors = FALSE)
    df$date   <- as.POSIXct(df$date, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    df$dea_uid <- uid
    df
  }, error = function(e) {
    warning("Failed to fetch: ", uid, " — ", conditionMessage(e))
    NULL
  })
}

message("Downloading ", nrow(urls_to_fetch), " time series...")
ts_list <- Map(fetch_ts, urls_to_fetch$dea_uid, urls_to_fetch$dea_ts_url)
ts_all  <- do.call(bind_rows, ts_list)

message(nrow(ts_all), " observations downloaded (", n_distinct(ts_all$dea_uid), " waterbodies)")

# ── 4. Compute annual dry-season metrics ──────────────────────────────────────
annual_metrics <- ts_all |>
  filter(!is.na(pc_wet)) |>
  mutate(
    year  = year(date),
    month = month(date)
  ) |>
  filter(month %in% DRY_SEASON_MONTHS) |>
  group_by(dea_uid, year) |>
  summarise(
    n_obs       = n(),
    mean_pc_wet = mean(pc_wet),
    max_pc_wet  = max(pc_wet),
    dry = any(pc_wet == 0),
    .groups = "drop"
  ) |>
  filter(n_obs >= MIN_OBS_PER_YEAR)

# ── 5. Overall drying probability per waterbody ───────────────────────────────
drying_probs <- annual_metrics |>
  group_by(dea_uid) |>
  summarise(
    n_years_dea     = n(),
    p_dry_empirical = mean(dry),
    mean_pc_wet_dry = mean(mean_pc_wet),
    .groups = "drop"
  )

# ── 6. Join back to living.waters ─────────────────────────────────────────────
living.waters.dea <- lw_matched |>
  st_transform(4326) |>
  left_join(drying_probs, by = "dea_uid")

# Clear unmatched points' DEA columns for clarity
living.waters.dea <- living.waters.dea |>
  mutate(across(c(p_dry_empirical, n_years_dea, mean_pc_wet_dry),
                ~ ifelse(dea_matched, ., NA_real_)))

message("Done. ", sum(!is.na(living.waters.dea$p_dry_empirical)),
        " waterpoints have empirical drying probabilities.")
message(sum(is.na(living.waters.dea$p_dry_empirical)),
        " waterpoints have no DEA match — p_dry set to NA.")
