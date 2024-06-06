# Some map outputs built with advice from this nice tutorial: https://www.littlemissdata.com/blog/maps

library(lubridate)
library(ggplot2)
library(ggrepel)
library(ggmap)
library(sf)
library(tidyverse)

# read in the point data
d <- read.csv(file = "out/spread-TCZ-timing.csv")
# cast to sf
d <- st_as_sf(d, coords = c("X", "Y"))
# set the CRS (Australian Albers)
d <- st_set_crs(d, 3577)

register_google(key = Sys.getenv("GOOGLE_MAP_API_KEY"))

basemap <- get_googlemap(
    center = c(lon = 121.5, lat = -19.123426),
    zoom = 7,
    size = c(640, 640),
    maptype = "hybrid"
    )

# function to hack the ggmap bounding box to map google map default
# from https://stackoverflow.com/questions/47749078/how-to-put-a-geom-sf-produced-map-on-top-of-a-ggmap-produced-raster
ggmap_bbox <- function(map) {
  if (!inherits(map, "ggmap")) stop("map must be a ggmap object")
  # Extract the bounding box (in lat/lon) from the ggmap to a numeric vector, 
  # and set the names to what sf::st_bbox expects:
  map_bbox <- setNames(unlist(attr(map, "bb")), 
                       c("ymin", "xmin", "ymax", "xmax"))
  
  # Coonvert the bbox to an sf polygon, transform it to 3857, 
  # and convert back to a bbox (convoluted, but it works)
  bbox_3857 <- st_bbox(st_transform(st_as_sfc(st_bbox(map_bbox, crs = 4326)), 3857))
  
  # Overwrite the bbox of the ggmap object with the transformed coordinates 
  attr(map, "bb")$ll.lat <- bbox_3857["ymin"]
  attr(map, "bb")$ll.lon <- bbox_3857["xmin"]
  attr(map, "bb")$ur.lat <- bbox_3857["ymax"]
  attr(map, "bb")$ur.lon <- bbox_3857["xmax"]
  map
}

basemap <- ggmap_bbox(basemap)

p <- ggmap(basemap) + 
  coord_sf(crs = st_crs(3857)) + # force the ggplot2 map to be in 3857
  geom_sf(data = st_transform(d, 3857), mapping = aes(colour = arrival), inherit.aes = FALSE)

p


