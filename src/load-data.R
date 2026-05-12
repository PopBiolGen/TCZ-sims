# Script to load all background data required for TCZ sims

######## load functions and libraries ########
source("src/pprocess_functions.R") # functions for point process model
source("src/get-invasion-front.R") # function for bringing in current invasion front

######## load posteriors ########
load("dat/Posteriors.RData")

######## define the data directories ########
data.dir <- file.path(Sys.getenv("DATA_PATH"), "Toads/TCZ/infrastructure")
spatial.dir <- file.path(Sys.getenv("DATA_PATH"), "GIS - General", "GIS_layers_read_only/")

######## Load point data ########
# load TCZ infrastructure data (downloaded as geojson)
data.dump.id <- "483d4936-3b85-41c1-bf75-eb0f0fda3c27"
tcz.sites <- st_read(file.path(data.dir, data.dump.id, "water_point_audit.geojson")) |> 
  mutate(origin_des = "Manmade", inside_tcz = TRUE) |> 
  select(-observer, -station_name, -traditional_owner_country, -whole_site_photo)
tcz.infrastructure <- st_read(file.path(data.dir, data.dump.id, "water_point_audit_infrastructure_item.geojson"))

# load background points (from Southwell et al)
old.lagrange.points <- st_read(file.path(spatial.dir, "Edited-layers/merged-points_rainfall_LaGrange.shp")) |> 
  select(fcsubtype_, full_name, perennia_1, origin_des, watercou_1, area_m, st_perimet) |> 
  st_transform(4326) # switch to WGS84 to match other data sources

  # read in additional points from aerial imagery ()
d_extra <- st_read("dat/tims_points.kml") |> 
  mutate(origin_des = "Manmade") 
  # add additional records onto old.lagrange.points
old.lagrange.points <- bind_rows(old.lagrange.points, d_extra)
rm(d_extra)

######## Load rasters ########
# load rainfall data (Number of days where at least 1mm of rain falls)
rainfall.raster <- terra::rast(file.path(spatial.dir, "annual-rainfall/rdann-1.asc"))

######## Load polygons ########
# load TCZ poly
tcz.boundary <- st_read(file.path(spatial.dir, "TCZ_boundary/toad_containment_zone.shp"))
# load wa coastline poly
wa.coast <- st_read(file.path(spatial.dir, "coast-poly_wa.shp")) |> 
  st_transform(crs = 4326)

####### Merge and filter datasets #######
# 1. Flag which lagrange points fall inside the TCZ boundary
inside_tcz <- st_within(old.lagrange.points, tcz.boundary, sparse = FALSE)[, 1]
old.lagrange.points <- old.lagrange.points |>
  mutate(inside_tcz = inside_tcz)

# 2. Remove Manmade points that are inside the TCZ
lagrange.filtered <- old.lagrange.points |>
  filter(!(origin_des == "Manmade" & inside_tcz))

# 3. Merge tcz.sites with filtered lagrange points
# (unmatched columns will fill with NA)
all.points <- bind_rows(lagrange.filtered, tcz.sites)

# 4. Extract rainfall values to the merged point set
rainfall.vals <- terra::extract(rainfall.raster, terra::vect(all.points))
all.points <- all.points |>
  mutate(rainfall = rainfall.vals[[2]])

# Bring in estimated invasion front and score east of there as colonised
colnsd <- score_colonised(all.points)
all.points <- colnsd$scored.points
inv.front <- colnsd$front; rm(colnsd)
  

####### Plot it #######
bbox <- st_bbox(all.points)

ggplot() +
  geom_sf(data = wa.coast, fill = NA, color = "grey30") +
  geom_sf(data = tcz.boundary, fill = NA, color = "red") +
  geom_sf(data = inv.front, color = "red", lty = 2) +
  geom_sf(data = all.points, aes(color = colonised)) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
           ylim = c(bbox["ymin"], bbox["ymax"]))
