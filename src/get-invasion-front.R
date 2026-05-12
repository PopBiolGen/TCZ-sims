# A function to get the invasion front from the invasion front inference tool
# The invasion front is represented as a line y = ax + b in Australian Albers
# This script gets the line and then makes all points east of the line colonised

score_colonised <- function(pts) {
  # load parameters
  load(file.path(Sys.getenv("DATA_PATH"), "invasion-front-parameters.Rdata"))
  a <- mod1[[1]]["a",1] # inferred a
  b <- mod1[[1]]["b",1]+mean.coord["Y"] # inferred b (intercept)
  
  init.crs <- st_crs(pts) # get original crs
  pts <- st_transform(pts, crs = 3577) # transform points to  Albers
  bb <- st_bbox(pts)
  
  # make a polyline matching a, b, within the bounding box of the points
  x_vals <- c(bb["xmin"], bb["xmax"])
  y_vals <- a * (x_vals-mean.coord["X"]) + b
  
  line <- st_sf(
    geometry = st_sfc(
      st_linestring(cbind(x_vals, y_vals)),
      crs = st_crs(pts)
    )
  )
  
  # score whether each point is east or west of line
  coords <- st_coordinates(pts)
  
  pts$colonised <- ifelse(
    coords[, "Y"] - a * (coords[, "X"]-mean.coord["X"]) - b > 0,
    1,  # east
    0   # west
  )
  
  # 
  # ggplot() +
  #   geom_sf(data = pts, aes(color = colonised)) +
  #   geom_sf(data = line, col = "red") 

  
  
  # transform back to original crs
  pts <- st_transform(pts, crs = init.crs)
  line <- st_transform(line, crs = init.crs)
  list(scored.points = pts, front = line)
}


