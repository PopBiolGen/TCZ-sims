# src/rainfall-drying-model.R
#
# Builds an empirical relationship between wet-season rainfall and dry-season
# drying probability for DEA-matched waterpoints.
#
# Wet season:  December (year - 1) through May (year)
# Response:    dry (logical) from annual_metrics
# Model:       glmer(dry ~ wet_mm + (1 | dea_uid), family = binomial)
#
# Inputs (must be in session):
#   annual_metrics    — from dea-drying-probs.R
#   living.waters.dea — from dea-drying-probs.R
#
# Outputs:
#   silo_daily      — daily SILO rainfall per site (tbl_df)
#   silo_monthly    — monthly totals per site (tbl_df)
#   wet_season_rain — wet-season totals per site-year (tbl_df)
#   model_data      — annual_metrics joined with wet_season_rain
#   fit_drying      — glmer model object

library(sf)
library(tidyverse)
library(lme4)
library(weatherOz)

SILO_KEY <- Sys.getenv("SILO_API_KEY")

# ── 1. Unique sites with coordinates ─────────────────────────────────────────
sites <- living.waters.dea |>
  filter(dea_matched) |>
  st_transform(4326) |>
  mutate(
    lon = st_coordinates(geometry)[, 1],
    lat = st_coordinates(geometry)[, 2]
  ) |>
  st_drop_geometry() |>
  group_by(dea_uid) |>
  slice(1) |>
  ungroup() |>
  select(dea_uid, lat, lon)

message(nrow(sites), " unique DEA sites to fetch")

# ── 2. Fetch daily rainfall from SILO for each site ──────────────────────────
# Dec 1986 – May 2025 spans all wet seasons preceding dry-season years 1987–2025
fetch_silo_rain <- function(uid, lat, lon) {
  message("  Fetching ", uid)
  tryCatch(
    get_data_drill(
      longitude  = lon,
      latitude   = lat,
      start_date = "1986-12-01",
      end_date   = "2025-05-31",
      values     = "rain",
      api_key    = SILO_KEY
    ) |>
      select(year, month, day, date, rainfall) |>
      mutate(dea_uid = uid),
    error = function(e) {
      warning("SILO fetch failed for ", uid, ": ", conditionMessage(e))
      NULL
    }
  )
}

message("Fetching SILO rainfall for ", nrow(sites), " sites ...")
silo_daily <- pmap(
  list(uid = sites$dea_uid, lat = sites$lat, lon = sites$lon),
  fetch_silo_rain
) |>
  list_rbind()

message(
  nrow(silo_daily), " daily observations fetched across ",
  n_distinct(silo_daily$dea_uid), " sites"
)

# ── 3. Monthly totals ─────────────────────────────────────────────────────────
silo_monthly <- silo_daily |>
  group_by(dea_uid, year, month) |>
  summarise(rain_mm = sum(rainfall, na.rm = TRUE), .groups = "drop")

# ── 4. Wet-season totals ──────────────────────────────────────────────────────
# Dec(Y-1) + Jan–May(Y) → assign Dec to the following dry-season year
wet_season_rain <- silo_monthly |>
  filter(month %in% c(12L, 1:5)) |>
  mutate(dry_year = if_else(month == 12L, year + 1L, as.integer(year))) |>
  group_by(dea_uid, dry_year) |>
  filter(n_distinct(month) == 6L) |>    # require all 6 months to be present
  summarise(wet_mm = sum(rain_mm), .groups = "drop") |>
  rename(year = dry_year)

# ── 5. Join to annual_metrics ─────────────────────────────────────────────────
model_data <- annual_metrics |>
  inner_join(wet_season_rain, by = c("dea_uid", "year")) |>
  mutate(
    wet_mm_sc  = as.numeric(scale(wet_mm)),
    log_wet_sc = as.numeric(scale(log(wet_mm + 1)))
  )

message(nrow(model_data), " site-years with both DEA and SILO data")
message(n_distinct(model_data$dea_uid), " sites represented")

# ── 6. Fit mixed-effects logistic regression ──────────────────────────────────
# wet_mm_sc: standardised linear rainfall
# log_wet_sc: standardised log-rainfall (often more linear on logit scale)

fit_drying <- glmer(
  dry ~ wet_mm_sc + (1 | dea_uid),
  data   = model_data,
  family = binomial(link = "logit")
)

fit_drying_log <- glmer(
  dry ~ log_wet_sc + (1 | dea_uid),
  data   = model_data,
  family = binomial(link = "logit")
)

summary(fit_drying)
summary(fit_drying_log)

# AIC comparison
AIC(fit_drying, fit_drying_log)
