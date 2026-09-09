# prep_app_data.R
#
# Run this script once after each annual forecast update to prepare data
# for the interactive Shiny app.
#
# Requires:
#   - out/forecast.RData         (produced by save_outputs())
#   - tcz.boundary               (sf; in session from src/tcz-sims/load-data.R, else
#                                 read from DATA_PATH below)
#   - pastoral.boundaries        (sf; in session from load-data.R, else read from
#                                 out/pastoral-boundaries.shp which load-data.R writes)
#
# Output:
#   - shiny/data/app_data.rds
#
# It is sourced at the end of src/tcz-sims/forecast.R (boundaries already in session).
# It can also be run standalone (e.g. to rebuild the shinylive site without re-running
# the forecast) as long as out/forecast.RData and out/pastoral-boundaries.shp exist.

library(dplyr)
library(sf)

start.year <- start_year

cat("Loading simulation output...\n")
load("out/forecast.RData")
cat("Processing", length(output), "simulation replicates.\n")

# ---- Boundaries: reuse session objects if present, otherwise load from disk ----
if (!exists("tcz.boundary")) {
  spatial.dir <- file.path(Sys.getenv("DATA_PATH"), "GIS - General",
                           "GIS_layers_read_only/")
  tcz.boundary <- st_read(
    file.path(spatial.dir, "TCZ_boundary/Toad_Containment_Zone_Boundary_July 26.shp"),
    quiet = TRUE
  )
}
if (!exists("pastoral.boundaries")) {
  pastoral.boundaries <- st_read("out/pastoral-boundaries.shp", quiet = TRUE)
}

# ---- Extract per-simulation arrival years ----
# Uses first.col.gen (generation a point was first colonised, relative to start.year) rather
# than back-solving from age/max.age: this is the same value setup()/spread.pilb() already
# maintain, and -- unlike the age-based calc -- it correctly carries a real historical
# colonisation year for points that were already colonised at forecast start (see
# setup()'s col.year.id/start.year args in pprocess_functions.R), instead of collapsing them
# all to start.year.
arrivals_list <- lapply(seq_along(output), function(i) {
  mat <- as.data.frame(output[[i]]$popmatrix)
  mat <- mat[mat$age > 0, ]
  if (nrow(mat) == 0) return(NULL)
  data.frame(
    ID      = mat$ID,
    X       = mat$X,
    Y       = mat$Y,
    arrival = round(start.year + mat$first.col.gen)
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

# ---- Waterpoints for background display ----
# Only points colonised in at least one replicate (i.e. present in all_arrivals).
# Points that are never colonised -- everything south of the fastest replicate's
# front, which stops at the first TCZ breach -- carry no information but, under
# shinylive/webR, each is a live browser marker. Dropping them cuts the app's
# memory and render cost. prob_data is already colonised-only (built above from
# all_arrivals), so all_points is the last place these points survive.
n_all_points <- nrow(as.data.frame(output[[1]]$popmatrix))
all_points <- as.data.frame(output[[1]]$popmatrix) |>
  select(ID, X, Y) |>
  filter(ID %in% unique(all_arrivals$ID)) |>
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

# ---- Pastoral boundaries in WGS84, trimmed to the TCZ latitude and north ----
# Pastoral polygons lying entirely south of the TCZ are off-screen context that
# only adds vertices (memory + SVG paths under webR). Keep any polygon that reaches
# the TCZ's southern edge or further north; drop the rest.
pastoral_wgs84 <- st_transform(pastoral.boundaries, 4326)
n_pastoral_all <- nrow(pastoral_wgs84)
pastoral_bbox <- st_bbox(pastoral_wgs84)
pastoral_keep_box <- st_as_sfc(st_bbox(c(
  xmin = unname(pastoral_bbox["xmin"]) - 1,
  xmax = unname(pastoral_bbox["xmax"]) + 1,
  ymin = unname(st_bbox(tcz_wgs84)["ymin"]) - 0.1,
  ymax = unname(pastoral_bbox["ymax"]) + 1
), crs = 4326))
pastoral_wgs84 <- st_filter(pastoral_wgs84, pastoral_keep_box, .predicate = st_intersects)

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
    start_year          = start.year,
    n_sims              = length(output),
    updated             = Sys.Date()
  ),
  "shiny/data/app_data.rds"
)

cat(
  "Saved shiny/data/app_data.rds\n",
  " Points         :", n_distinct(prob_data$ID), "\n",
  " Background pts  :", nrow(all_points), "of", n_all_points,
  "(dropped", n_all_points - nrow(all_points), "never-colonised)\n",
  " Pastoral polys  :", nrow(pastoral_wgs84), "of", n_pastoral_all,
  "(dropped", n_pastoral_all - nrow(pastoral_wgs84), "south of TCZ)\n",
  " Years          :", length(years), "\n",
  " Rows           :", nrow(prob_data), "\n"
)
