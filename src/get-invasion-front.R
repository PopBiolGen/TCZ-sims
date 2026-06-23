# A function to get the invasion front from the invasion front inference tool
# The invasion front is represented as a line y = ax + b in Australian Albers
# This script gets the line and then makes all points east of the line colonised

score_colonised <- function(pts) {
  # load parameters
  load(file.path(Sys.getenv("DATA_PATH"), "invasion-front-parameters.Rdata"))
  a <- mod1[[1]]["a",1] # inferred a
  b <- mod1[[1]]["b",1]+mean.coord["Y"] # inferred b (intercept)
  
  ### Find a northern boundary W of the Fitzroy
  # load a point representing the mouth of the Fitzroy that the front will have to go around
  fp <- data.frame(lat = -17.633614, lon = 123.444807)
  # cast this to a spatial dataframe with initial web mercator CRS, transform to albers
  fp <- st_as_sf(fp, coords = c("lon", "lat"), crs = 4326) |> 
    st_transform(crs = 3577)
  # get perpendicular distance from this point to the invasion front
  coords.fp <- st_coordinates(fp)
  num_fp <- a * (coords.fp[, "X"] - mean.coord["X"]) - coords.fp[, "Y"] + b
  dist.to.fp <- abs(num_fp) / sqrt(a^2 + 1)
  # northing of the foot of the perpendicular from fp to the invasion front
  Y_intersect <- coords.fp[, "Y"] + num_fp / (a^2 + 1)
  
  ### Get a southern boundary 
  sp <- data.frame(lat = -18.514959, lon = 123.069856) |> 
    st_as_sf(coords = c("lon", "lat"), crs = 4326) |> 
    st_transform(crs = 3577) |> 
    st_coordinates()
  
  ### Get our points
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
  
  # score whether each point is colonised
  coords <- st_coordinates(pts)

  east.of.line <- coords[, "Y"] - a * (coords[, "X"] - mean.coord["X"]) - b > 0
  north.of.sp <- coords[, "Y"] > sp[1, "Y"]
  east.of.fp <- coords[, "X"] > coords.fp[1, "X"]
  north.of.intersect <- coords[, "Y"] > Y_intersect
  dist.to.pt      <- sqrt((coords[, "X"] - coords.fp[1, "X"])^2 + (coords[, "Y"] - coords.fp[1, "Y"])^2)
  
  pts$colonised <- as.integer(
    north.of.sp & east.of.line & (!north.of.intersect | (dist.to.pt <= dist.to.fp) | (east.of.fp & north.of.intersect))
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


