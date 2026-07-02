# prep_app_data.R
#
# Run this script once after each annual forecast update to prepare data
# for the interactive Shiny app.
#
# Requires:
#   - out/forecast.RData  (produced by save_outputs())
#   - tcz.boundary        (sf object in session, loaded via src/load-data.R)
#
# Output:
#   - shiny/data/app_data.rds

library(dplyr)
library(sf)

start.year <- 2025

cat("Loading simulation output...\n")
load("out/forecast.RData")
cat("Processing", length(output), "simulation replicates.\n")

# ---- Extract per-simulation arrival years ----
arrivals_list <- lapply(seq_along(output), function(i) {
  mat <- as.data.frame(output[[i]]$popmatrix)
  mat <- mat[mat$age > 0, ]
  if (nrow(mat) == 0) return(NULL)
  max.age <- max(mat$age)
  data.frame(
    ID      = mat$ID,
    X       = mat$X,
    Y       = mat$Y,
    arrival = round(start.year + max.age - mat$age)
  )
})
all_arrivals <- bind_rows(arrivals_list)

years <- start.year:max(all_arrivals$arrival)
cat("Years covered:", min(years), "to", max(years), "\n")

# ---- Per-point, per-year colonisation probability ----
cat("Computing per-year probabilities...\n")
prob_data <- all_arrivals |>
  group_by(ID, X, Y) |>
  reframe(
    year = years,
    prob  = sapply(years, function(y) mean(arrival <= y))
  )

# ---- Convert Albers (3577) to WGS84 for leaflet ----
cat("Converting coordinates to WGS84...\n")
coords <- all_arrivals |>
  select(ID, X, Y) |>
  distinct() |>
  st_as_sf(coords = c("X", "Y"), crs = 3577) |>
  st_transform(4326) |>
  mutate(
    lon = st_coordinates(geometry)[, 1],
    lat = st_coordinates(geometry)[, 2]
  ) |>
  st_drop_geometry() |>
  select(ID, lon, lat)

# ---- Mean arrival year per point (for map popups) ----
mean_arrival <- all_arrivals |>
  group_by(ID) |>
  summarise(mean_arrival = round(mean(arrival)), .groups = "drop")

prob_data <- prob_data |>
  left_join(coords, by = "ID") |>
  left_join(mean_arrival, by = "ID") |>
  select(ID, lon, lat, year, prob, mean_arrival)

# ---- All waterpoints for background display (including never-colonised) ----
all_points <- as.data.frame(output[[1]]$popmatrix) |>
  select(ID, X, Y) |>
  st_as_sf(coords = c("X", "Y"), crs = 3577) |>
  st_transform(4326) |>
  mutate(
    lon = st_coordinates(geometry)[, 1],
    lat = st_coordinates(geometry)[, 2]
  ) |>
  st_drop_geometry() |>
  select(ID, lon, lat)

# ---- TCZ boundary in WGS84 ----
tcz_wgs84 <- st_transform(tcz.boundary, 4326)

# ---- Pastoral boundaries in WGS84 ----
pastoral_wgs84 <- st_transform(pastoral.boundaries, 4326)

# ---- Initial map extent: span from TCZ boundary to colonised waterpoints ----
# Using colonised waterpoints rather than inv.front, which can have geometry
# extending well beyond actual toad locations.
colonised_bbox <- all_arrivals |>
  filter(arrival <= start.year) |>
  select(X, Y) |>
  distinct() |>
  st_as_sf(coords = c("X", "Y"), crs = 3577) |>
  st_transform(4326) |>
  st_bbox()

tcz_bbox <- st_bbox(tcz_wgs84)

initial_bounds <- list(
  lng1 = unname(tcz_bbox["xmin"])      - 0.5,
  lat1 = unname(tcz_bbox["ymin"])      - 0.5,
  lng2 = unname(colonised_bbox["xmax"]) + 0.5,
  lat2 = unname(colonised_bbox["ymax"]) + 0.5
)

# ---- Save ----
saveRDS(
  list(
    prob_data           = prob_data,
    all_points          = all_points,
    tcz_boundary        = tcz_wgs84,
    pastoral_boundaries = pastoral_wgs84,
    initial_bounds      = initial_bounds,
    years               = years,
    n_sims              = length(output),
    updated             = Sys.Date()
  ),
  "shiny/data/app_data.rds"
)

cat(
  "Saved shiny/data/app_data.rds\n",
  " Points :", n_distinct(prob_data$ID), "\n",
  " Years  :", length(years), "\n",
  " Rows   :", nrow(prob_data), "\n"
)
