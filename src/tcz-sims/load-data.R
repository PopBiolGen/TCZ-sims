# Script to load all background data required for TCZ sims

######## load functions and libraries ########
source("src/pprocess_functions.R") # functions for point process model
source("src/get-invasion-front.R") # function for bringing in current invasion front

######## load posteriors ########
load("dat/Posteriors_2026.RData")

######## define the data directories ########
data.dir <- file.path(Sys.getenv("DATA_PATH"), "Toads/TCZ/infrastructure")
spatial.dir <- file.path(Sys.getenv("DATA_PATH"), "GIS - General", "GIS_layers_read_only/")

######## Load point data ########
# load TCZ infrastructure data (downloaded as geojson)
data.dump.id <- "089e7bef-6ea8-45d7-a819-a1c964d552fb"
tcz.sites <- st_read(file.path(data.dir, data.dump.id, "water_point_audit.geojson")) |> 
  mutate(origin_des = "Manmade", inside_tcz = TRUE) |> 
  select(-observer, -station_name, -traditional_owner_country, -whole_site_photo)
tcz.infrastructure <- st_read(
  file.path(data.dir, 
            data.dump.id, 
            "water_point_audit_infrastructure_item.geojson"))

# load costing data
## Total costs
total.cost <- readxl::read_xlsx(path = file.path(data.dir, "Activity Steps Toad project - 4 crew - 9Jul26.xlsx"),
                                sheet = "Analysis") |> 
  select('Total Cost') |> 
  slice(1) |> 
  unlist()
## Site-level costs
tcz.costs <- readxl::read_xlsx(path = file.path(data.dir, "Activity Steps Toad project - 4 crew - 9Jul26.xlsx"),
                               sheet = "Activity Steps Toad project",
                               skip = 1) |> 
  filter(!is.na(Activity))
names(tcz.costs) <- tolower(make.names(names(tcz.costs)))
tcz.costs <- tcz.costs |> 
  group_by(activity) |> 
  summarise(site.id = first(activity.id),
            site = first(activity),
            effort = first(site.duration)) |> 
  mutate(cost = effort/sum(effort) * total.cost)
#### To Do ####
# merge costs to infrastructure
normalize_name <- function(x) {
  x |>
    gsub("\u2019|\u2018", "'", x = _) |>   # curly single quotes → straight
    gsub("\u201C|\u201D", '"', x = _) |>   # curly double quotes → straight
    trimws() |>
    tolower()
}

costs_for_join <- tcz.costs |>
  mutate(activity = recode(activity,
    "Anna Plains Homestead Perimeter Fence and Gates" = "Anna Plains Homestead",
    "Eeganah well king brown well"                    = "Eeganah well \"king brown well\"",
    "Gilbert 2"                                       = "Gilberts 2",
    "Herbs 1"                                         = "Herb's 1",
    "No.5"                                            = "No. 5",
    "Wotans Bore"                                     = "Wotens bore",
    "Yard bore - Nita"                                = "Yard bore"
  )) |>
  mutate(.key = normalize_name(activity)) |>
  select(.key, cost)

tcz.sites <- tcz.sites |>
  mutate(.key = normalize_name(name_of_site)) |>
  left_join(costs_for_join, by = ".key") |>
  select(-.key)

# merge the two sets to apply ruleset for toad colonisable points
temp <- left_join(tcz.sites, st_drop_geometry(tcz.infrastructure), by = "X_record_id") |> 
  select(X_record_id, 
         name_of_site, 
         brief_description, 
         infrastructure_item, 
         overall_site_comments, 
         fs_control_device_comment, 
         proposed_works_at_water_point,
         item_type,
         controls_required,
         new_items,
         origin_des,
         inside_tcz,
         geometry)
# identify sites with works required
something.ss <- !vapply(temp$proposed_works_at_water_point, function(x){"A - no works required" %in% x}, FUN.VALUE = logical(1))
# take only infrastructure with something to do and collapse back to site-level
tcz.sites <- temp[something.ss,] |>
  group_by(X_record_id) |>
  summarise(
    # Columns with consistent values — take first
    name_of_site                  = first(name_of_site),
    brief_description             = first(brief_description),
    overall_site_comments         = first(overall_site_comments),
    fs_control_device_comment     = first(fs_control_device_comment),
    proposed_works_at_water_point = list(unique(proposed_works_at_water_point)),
    # Columns with multiple values — collapse into list
    infrastructure_item           = list(infrastructure_item),
    item_type                     = list(item_type),
    controls_required             = list(controls_required),
    new_items                     = list(new_items),
    origin_des                    = "Manmade" # all TCZ infrastructure points are Manmade
  )
rm(temp)

# items that contain fences
#fence.ss <- vapply(temp$item_type, function(x){any(c("Fence", "Gate") %in% x)}, FUN.VALUE = logical(1))

# load background points (from Southwell et al)
old.lagrange.points <- st_read(file.path(spatial.dir, "Edited-layers/merged-points_rainfall_LaGrange.shp")) |> 
  select(fcsubtype_, full_name, perennia_1, origin_des, watercou_1, area_m, st_perimet) |> 
  st_transform(4326) # switch to WGS84 to match other data sources

## Note the manmade points in the TCZ here ^ have been replaced by the audit in tcz.sites
# these old.lagrange.points also do a bad job of identifying natural waterpoints in TCZ.  Need to replace with something else
# Living waters points make the most sense
living.waters <- st_read(file.path(spatial.dir, "living-waters/living-waters-points.shp")) |>
  mutate(origin_des = "Natural") |> 
  st_transform(4326) # switch to WGS84 to match other data sources


# read in additional points from aerial imagery ()
d_extra <- st_read(file.path(spatial.dir, "additional-waterpoints/additional-waterpoints.kml")) |> 
  mutate(
    origin_des = case_when(
      # "point" in name → Manmade (catches watering point, waterpoint, etc.)
      str_detect(Name, regex("point", ignore_case = TRUE)) ~ "Manmade",
      # Natural water features
      str_detect(Name, regex("lagoon|lgoon|lake|geegully|creek|seep|spring|soak|salt", ignore_case = TRUE)) ~ "Natural",
      # Remaining manmade infrastructure
      str_detect(Name, regex("water|dam|tank|trough|irrigat|borrow|pond|ditch|resort|homestead|campsite|station|park", ignore_case = TRUE)) ~ "Manmade",
      TRUE ~ NA_character_
    )
  )
  # add additional records onto old.lagrange.points
old.lagrange.points <- bind_rows(old.lagrange.points, d_extra)
rm(d_extra)

######## Load rasters ########
# load rainfall data (Number of days where at least 1mm of rain falls)
rainfall.raster <- terra::rast(file.path(spatial.dir, "annual-rainfall/rdann-1.asc"))

######## Load polygons ########
# load TCZ poly
tcz.boundary <- st_read(file.path(spatial.dir, "TCZ_boundary/Toad_Containment_Zone_Boundary_July 26.shp"))
# load wa coastline poly
wa.coast <- st_read(file.path(spatial.dir, "coast-poly_wa.shp")) |> 
  st_transform(crs = 4326) |>
  st_make_valid() |>
  st_intersection(old.lagrange.points |> st_transform(4326) |> st_bbox() |> st_as_sfc())
# load pastoral boundaries
pastoral.boundaries <- st_read(file.path(spatial.dir, "Pastoral_Stations_DPLH_083.gdb")) |> 
  st_transform(crs = 4326) |>
  st_make_valid() |>
  st_intersection(old.lagrange.points |> st_transform(4326) |> st_bbox() |> st_as_sfc())

####### Merge and filter datasets #######
# 1. Flag which lagrange points fall inside the TCZ boundary
inside_tcz <- st_within(old.lagrange.points, tcz.boundary, sparse = FALSE)[, 1]
old.lagrange.points <- old.lagrange.points |>
  mutate(inside_tcz = inside_tcz)

# 2. Remove Manmade points that are inside the TCZ
lagrange.filtered <- old.lagrange.points |>
  filter(!inside_tcz) # remove all these points, natural and manmade

# 3. Merge tcz.sites with filtered lagrange points
# (unmatched columns will fill with NA)
all.points <- bind_rows(lagrange.filtered, tcz.sites, living.waters) 
all.points$inside_tcz <- st_within(all.points, tcz.boundary, sparse = FALSE)[, 1]

# 4. Extract rainfall values to the merged point set
rainfall.vals <- terra::extract(rainfall.raster, terra::vect(all.points))
all.points <- all.points |>
  mutate(rainfall = rainfall.vals[[2]])

####### Bring in estimated invasion front and scored colonised points #######
colnsd <- score_colonised(all.points)
all.points <- colnsd$scored.points
inv.front <- colnsd$front
fp <- colnsd$fp; rm(colnsd)

# get coordinates in albers
albers <- all.points |> 
  st_transform(crs = 3577) |> 
  st_coordinates()
all.points <- bind_cols(all.points, albers)

####### Plot it #######
bbox <- st_bbox(all.points)

time.to.tcz.fig <-ggplot() +
  geom_sf(data = wa.coast, fill = NA, color = "grey30") +
  geom_sf(data = tcz.boundary, fill = NA, color = "red") +
  geom_sf(data = inv.front, color = "red", lty = 2) +
  geom_sf(data = all.points, aes(color = factor(colonised))) +
  scale_color_manual(values = c("0" = "#132B43", "1" = "#56B1F7")) +
  geom_sf(data = fp, color = "red") +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
           ylim = c(bbox["ymin"], bbox["ymax"]))
ggsave(filename = "out/time-to-tcz.png", plot = time.to.tcz.fig)

bbox <- st_bbox(tcz.boundary)
tcz.waterpoints.fig <- ggplot() +
  geom_sf(data = wa.coast, fill = NA, color = "grey30") +
  geom_sf(data = tcz.boundary, fill = NA, color = "red") +
  geom_sf(data = all.points, aes(color = origin_des)) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
           ylim = c(bbox["ymin"], bbox["ymax"]))
ggsave(filename = "out/tcz-waterpoints.pdf", plot = tcz.waterpoints.fig)
