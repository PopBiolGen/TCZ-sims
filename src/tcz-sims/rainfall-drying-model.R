# src/rainfall-drying-model.R
#
# Builds an empirical relationship between wet-season rainfall and dry-season
# drying probability for DEA-matched waterpoints.
#
# Wet season:  December (year - 1) through May (year)
# Response:    dry (logical) from annual_metrics
# Model:       glmer(dry ~ wet_mm_sc + wet_mm_lag1_sc + (1 | dea_uid),
#                    family = binomial)
#
# The random intercept per waterpoint captures baseline drying tendency
# and enables site-specific predictions via predict(..., re.form = NULL).
# Scaling parameters are stored in `rain_scale` for use with new data.
#
# Inputs (must be in session):
#   annual_metrics    — from dea-drying-probs.R
#   living.waters.dea — from dea-drying-probs.R
#
# Outputs:
#   silo_daily      — daily SILO rainfall per site (tbl_df)
#   silo_monthly    — monthly totals per site (tbl_df)
#   wet_season_rain — wet-season totals + lag-1 per site-year (tbl_df)
#   model_data      — annual_metrics joined with wet_season_rain
#   rain_scale      — list(wet_mm, wet_mm_lag1) of scale() attributes
#   fit_drying      — glmer model object (primary model)

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
# Start Dec 1985 so the lag-1 wet season is available for dry-season year 1987
fetch_silo_rain <- function(uid, lat, lon) {
  message("  Fetching ", uid)
  tryCatch(
    get_data_drill(
      longitude  = lon,
      latitude   = lat,
      start_date = "1985-12-01",
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

# ── 4. Wet-season totals + lag-1 ──────────────────────────────────────────────
# Dec(Y-1) + Jan–May(Y) → dry_year Y
# wet_mm_lag1: the wet season from the prior year (dry_year Y-1)
wet_season_rain <- silo_monthly |>
  filter(month %in% c(12L, 1:5)) |>
  mutate(dry_year = if_else(month == 12L, year + 1L, as.integer(year))) |>
  group_by(dea_uid, dry_year) |>
  filter(n_distinct(month) == 6L) |>    # require all 6 months present
  summarise(wet_mm = sum(rain_mm), .groups = "drop") |>
  rename(year = dry_year) |>
  arrange(dea_uid, year) |>
  group_by(dea_uid) |>
  mutate(wet_mm_lag1 = lag(wet_mm)) |>  # previous wet season
  ungroup()

# ── 5. Join to annual_metrics and scale predictors ────────────────────────────
# Scaling is applied globally (not per-site) so the intercept variance
# remains interpretable. Attributes are stored for prediction on new data.
model_data_raw <- annual_metrics |>
  inner_join(wet_season_rain, by = c("dea_uid", "year")) |>
  # Drop rows missing the lag (typically the first observed year per site)
  filter(!is.na(wet_mm_lag1))

sc_wet      <- scale(model_data_raw$wet_mm)
sc_wet_lag1 <- scale(model_data_raw$wet_mm_lag1)

rain_scale <- list(
  wet_mm      = list(center = attr(sc_wet, "scaled:center"),
                     scale  = attr(sc_wet, "scaled:scale")),
  wet_mm_lag1 = list(center = attr(sc_wet_lag1, "scaled:center"),
                     scale  = attr(sc_wet_lag1, "scaled:scale"))
)

model_data <- model_data_raw |>
  mutate(
    wet_mm_sc      = as.numeric(sc_wet),
    wet_mm_lag1_sc = as.numeric(sc_wet_lag1)
  )

message(nrow(model_data), " site-years with both DEA and SILO data")
message(n_distinct(model_data$dea_uid), " sites represented")

# ── 6. Fit mixed-effects logistic regression ──────────────────────────────────
# Random intercept per waterpoint captures baseline drying tendency.
# Predictions for known sites use predict(..., re.form = NULL) to include BLUPs.
# Predictions for new sites use predict(..., re.form = NA) for the population mean.

# Current wet season only (baseline)
fit_drying_current <- glmer(
  dry ~ wet_mm_sc + (1 | dea_uid),
  data   = model_data,
  family = binomial(link = "logit")
)

# Current + previous wet season
fit_drying <- glmer(
  dry ~ wet_mm_sc + wet_mm_lag1_sc + (1 | dea_uid),
  data   = model_data,
  family = binomial(link = "logit")
)

summary(fit_drying_current)
summary(fit_drying)

# AIC comparison: does the lag-1 term improve fit?
AIC(fit_drying_current, fit_drying)

# ── 7. Predict drying probability for each site-year ─────────────────────────
# re.form = NULL → include site-specific random intercept (best for known sites)
drying_probs <- model_data |>
  mutate(
    p_dry = predict(fit_drying, newdata = model_data,
                    type = "response", re.form = NULL)
  ) |>
  select(dea_uid, year, wet_mm, wet_mm_lag1, dry, p_dry)

# Helper: predict p(dry) for a new rainfall record at a known site
# Supply wet_mm and wet_mm_lag1 in mm; uses stored scaling params.
predict_drying <- function(dea_uid, wet_mm, wet_mm_lag1, model = fit_drying) {
  newdata <- tibble(
    dea_uid        = dea_uid,
    wet_mm_sc      = (wet_mm      - rain_scale$wet_mm$center)      / rain_scale$wet_mm$scale,
    wet_mm_lag1_sc = (wet_mm_lag1 - rain_scale$wet_mm_lag1$center) / rain_scale$wet_mm_lag1$scale
  )
  predict(model, newdata = newdata, type = "response", re.form = NULL)
}
