# Script to load all background data required for TCZ sims

######## load functions and libraries ########
source("src/pprocess_functions.R")

######## load posteriors ########
load("dat/Posteriors.RData")

# define the data directories
data.dir <- file.path(Sys.getenv("DATA_PATH"), "Toads/TCZ/infrastructure")
spatial.dir <- file.path(Sys.getenv("DATA_PATH"), "GIS - General", "GIS_layers_read_only/")


# load TCZ infrastructure data 
data.dump.id <- "483d4936-3b85-41c1-bf75-eb0f0fda3c27"
tcz.sites <- st_read(file.path(data.dir, data.dump.id, "water_point_audit.geojson"))
tcz.infrastructure <- st_read(file.path(data.dir, data.dump.id, "water_point_audit_infrastructure_item.geojson"))

# load background points (from Southwell et al)
old.lagrange.points <- st_read(file.path(spatial.dir, "Edited-layers/merged-points_rainfall_LaGrange.shp")) |> 
  select(fcsubtype_, full_name, perennia_1, origin_des, watercou_1, area_m, st_perimet)
  # read in additional points from aerial imagery ()
d_extra <- st_read("dat/tims_points.kml") |> 
  st_transform(d_extra, crs = 3577) |>  # convert to Albers to match old.lagrange.points
  mutate(origin_des = "Manmade") 
  # add additional records onto old.lagrange.points
old.lagrange.points <- bind_rows(old.lagrange.points, d_extra)



# load rainfall data (Number of days where at least 1mm of rain falls)
rainfall.raster <- rast(file.path(spatial.dir, "annual-rainfall/rdann-1.asc"))
