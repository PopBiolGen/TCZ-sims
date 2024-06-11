# useful resources:
# https://cran.r-project.org/web/packages/tmap/vignettes/tmap-getstarted.html
# https://www.r-bloggers.com/2019/10/map-coloring-the-color-scale-styles-available-in-the-tmap-package/

library(ggplot2)
library(sf)
library(tmap)
library(tmaptools)
library(dplyr)
library(magick)

######### Make a static map of estimated arrival time #########
# function to read point data (in Albers) and cast to sf with a CRS
read.point.data <- function(fname) {
  # read in the point data
  d <- read.csv(file = fname)
  # cast to sf
  d <- st_as_sf(d, coords = c("X", "Y"))
  # set the CRS (Australian Albers)
  d <- st_set_crs(d, 3577)
  d <- st_transform(d, 3857) # switch to web map default CRS
}

d <- read.point.data("out/spread-TCZ-timing.csv")



# read a satellite image basemap from ESRI
bm <- read_osm(
  d,
  type = "osm", # for satellite image, "esri-imagery",
  zoom = 8,
  ext = 1
)

p <- tm_shape(bm,
         unit = "km") +
  tm_rgb() +
  tm_shape(d) +
  tm_dots(size = 0.2,
          col = "arrival",
          breaks = 2024:2042, 
          legend.format = list(big.mark = ""),
          title = "Predicted year of toad arrival") 

tmap_save(p, filename = "out/year-of-arrival.pdf")

# other options for plotting
# ggplot() +
#   stars::geom_stars()
# 
# terra::rast(bm)
# 
# ggplot() +
#   geom_spatraster

######### Make a dynamic map of estimated arrival time #########
fpath <- "out/dynamic_maps/"
input.flist <- list.files(path = fpath, pattern = ".csv")
for (yy in input.flist){
  temp <- read.point.data(fname = paste0(fpath, yy))
  year.name <- gsub(".csv", "", yy)
  fname <- paste0(fpath, year.name, ".png")
  p <- tm_shape(bm,
                unit = "km") +
    tm_rgb() +
    tm_shape(temp) +
    tm_dots(size = 0.2,
            col = "prob.colonised",
            breaks = seq(0, 1, length.out = 5),
            title = "Probability of colonisation") +
    tm_layout(title = year.name)
  tmap_save(p, filename = fname)
}

flist <- list.files(path = fpath, pattern = ".png")
images <- image_read(paste0(fpath, flist))
animation <- image_animate(images, fps = 1)
image_write(animation, path = "out/dynamic_maps/animated_map.gif")

# clean up
file.remove(paste0(fpath, input.flist))
file.remove(paste0(fpath, flist))