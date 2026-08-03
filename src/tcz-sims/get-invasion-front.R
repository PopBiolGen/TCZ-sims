# A function to get the invasion front from the invasion front inference tool
# The invasion front is represented as a line y = ax + b in Australian Albers
# This script gets the line and then makes all points east of the line colonised

score_colonised <- function(pts) {
  # load parameters
  load(file.path(Sys.getenv("DATA_PATH"), "invasion-front-parameters.Rdata"))
  a <- mod1[[1]]["a",1] # inferred a (slope; dimensionless, unaffected by equal scaling of X and Y)
  b_raw <- mod1[[1]]["b",1] # inferred intercept in centred-km model space
  
  # mean.coord is in km (coordinates were divided by scale = 1000 before centring).
  # Convert everything to Albers metres for use with sf geometries.
  # In model space: y_km_centred = a * x_km_centred + b_raw
  # In Albers metres: Y = a * (X - mean.coord_m["X"]) + b
  mean.coord_m <- mean.coord * scale
  b <- (b_raw + mean.coord["Y"]) * scale  # effective y-intercept in Albers metres
  
  ### Find a northern boundary W of the Fitzroy
  # load a point representing the mouth of the Fitzroy that the front will have to go around
  fp <- data.frame(lat = -17.633614, lon = 123.444807)
  # cast this to a spatial dataframe with initial web mercator CRS, transform to albers
  fp <- st_as_sf(fp, coords = c("lon", "lat"), crs = 4326) |> 
    st_transform(crs = 3577)
  # get perpendicular distance from this point to the invasion front
  coords.fp <- st_coordinates(fp)
  num_fp <- a * (coords.fp[, "X"] - mean.coord_m["X"]) - coords.fp[, "Y"] + b
  dist.to.fp <- abs(num_fp) / sqrt(a^2 + 1)
  # northing of the foot of the perpendicular from fp to the invasion front
  Y_intersect <- coords.fp[, "Y"] + num_fp / (a^2 + 1)
  
  ### Get a southern boundary 
  sp <- data.frame(lat = -19.06, lon = 123.069856) |> 
    st_as_sf(coords = c("lon", "lat"), crs = 4326) |> 
    st_transform(crs = 3577) |> 
    st_coordinates()
  
  ### Get our points
  init.crs <- st_crs(pts) # get original crs
  pts <- st_transform(pts, crs = 3577) # transform points to  Albers
  bb <- st_bbox(pts)
  
  # make a polyline matching a, b, within the bounding box of the points (Albers metres)
  x_vals <- c(bb["xmin"], bb["xmax"])
  y_vals <- a * (x_vals - mean.coord_m["X"]) + b
  
  line <- st_sf(
    geometry = st_sfc(
      st_linestring(cbind(x_vals, y_vals)),
      crs = 3577
    )
  )
  
  # score whether each point is colonised
  coords <- st_coordinates(pts)

  east.of.line <- coords[, "Y"] - a * (coords[, "X"] - mean.coord_m["X"]) - b > 0
  north.of.sp <- coords[, "Y"] > sp[1, "Y"]
  east.of.fp <- coords[, "X"] > coords.fp[1, "X"]
  north.of.intersect <- coords[, "Y"] > Y_intersect
  dist.to.pt      <- sqrt((coords[, "X"] - coords.fp[1, "X"])^2 + (coords[, "Y"] - coords.fp[1, "Y"])^2)
  # has the front itself already passed fp? if not, the "go around the river
  # mouth" correction below is premature and shouldn't be applied yet
  fp.reached <- num_fp < 0

  pts$colonised <- as.integer(
    north.of.sp & east.of.line &
      (!fp.reached | !north.of.intersect | (dist.to.pt <= dist.to.fp) | (east.of.fp & north.of.intersect))
  )
  
  # 
  # ggplot() +
  #   geom_sf(data = pts, aes(color = colonised)) +
  #   geom_sf(data = line, col = "red") 

  
  
  # transform back to original crs
  # densify line before transforming so the straight Albers line isn't
  # approximated by only 2 vertices in geographic space
  pts <- st_transform(pts, crs = init.crs)
  line <- st_segmentize(line, dfMaxLength = 10000) |>  # 10 km intervals
    st_transform(crs = init.crs)
  fp <- st_transform(fp, crs = init.crs)
  list(scored.points = pts, front = line, fp = fp)
}


# Multi-year version: uses the shared-slope/year-specific-intercept model
# (invasion-front-parameters-multi-year.Rdata) to derive, for each point, the
# earliest modelled year whose front line it's already east of -- i.e. an
# estimated colonisation year rather than just a current colonised/not flag.
# Assumes the yearly fronts are nested/advancing (a point colonised under an
# earlier year's line stays colonised under every later year's line); this is
# checked but not enforced -- posterior noise can occasionally make a line
# retreat slightly, in which case the earliest-colonising year still wins.
#
# Returns colonisation_year = NA for points not yet colonised even under the
# latest modelled year's front. Points colonised before the earliest modelled
# year (years.all[1]) are left-censored: they get stamped with years.all[1]
# even though the true colonisation year may be earlier -- there's no data to
# distinguish further back than that.
score_colonised_multi_year <- function(pts, diagnostic_plot = FALSE) {
  # load multi-year parameters: mod.multi, mean.coord, scale, years.all, run.info
  load(file.path(Sys.getenv("DATA_PATH"), "invasion-front-parameters-multi-year.Rdata"))
  a <- mod.multi[[1]]["a", 1] # shared slope across years

  mean.coord_m <- mean.coord * scale

  ### Find a northern boundary W of the Fitzroy (same anchor as score_colonised())
  fp <- data.frame(lat = -17.633614, lon = 123.444807)
  fp <- st_as_sf(fp, coords = c("lon", "lat"), crs = 4326) |>
    st_transform(crs = 3577)
  coords.fp <- st_coordinates(fp)

  ### Get a southern boundary
  sp <- data.frame(lat = -19.06, lon = 123.069856) |>
    st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
    st_transform(crs = 3577) |>
    st_coordinates()

  ### Get our points
  init.crs <- st_crs(pts) # get original crs
  pts <- st_transform(pts, crs = 3577) # transform points to Albers
  bb <- st_bbox(pts)
  coords <- st_coordinates(pts)
  x_vals <- c(bb["xmin"], bb["xmax"]) # bounding box, for drawing each year's line

  n.years <- length(years.all)
  col_year <- rep(NA_real_, nrow(pts))
  lines.list <- vector("list", n.years)

  # step through years oldest -> newest; the first year a point classifies as
  # colonised is taken as its estimated colonisation year
  for (k in seq_len(n.years)) {
    b_raw_k <- mod.multi[[1]][paste0("b[", k, "]"), 1]
    b_k <- (b_raw_k + mean.coord["Y"]) * scale # effective y-intercept in Albers metres, year k

    num_fp <- a * (coords.fp[, "X"] - mean.coord_m["X"]) - coords.fp[, "Y"] + b_k
    dist.to.fp <- abs(num_fp) / sqrt(a^2 + 1)
    Y_intersect <- coords.fp[, "Y"] + num_fp / (a^2 + 1)

    east.of.line <- coords[, "Y"] - a * (coords[, "X"] - mean.coord_m["X"]) - b_k > 0
    north.of.sp <- coords[, "Y"] > sp[1, "Y"]
    east.of.fp <- coords[, "X"] > coords.fp[1, "X"]
    north.of.intersect <- coords[, "Y"] > Y_intersect
    dist.to.pt <- sqrt((coords[, "X"] - coords.fp[1, "X"])^2 + (coords[, "Y"] - coords.fp[1, "Y"])^2)
    # has this year's front already passed fp? if not, the "go around the
    # river mouth" correction below is premature and shouldn't be applied yet
    fp.reached <- num_fp < 0

    colonised.k <- as.logical(
      north.of.sp & east.of.line &
        (!fp.reached | !north.of.intersect | (dist.to.pt <= dist.to.fp) | (east.of.fp & north.of.intersect))
    )

    newly.colonised <- colonised.k & is.na(col_year)
    col_year[newly.colonised] <- years.all[k]

    y_vals <- a * (x_vals - mean.coord_m["X"]) + b_k
    lines.list[[k]] <- st_sf(
      year = years.all[k],
      geometry = st_sfc(st_linestring(cbind(x_vals, y_vals)), crs = 3577)
    )
  }

  pts$colonisation_year <- col_year
  pts$colonised <- as.integer(!is.na(col_year)) # = colonised under the latest year's line

  fronts <- do.call(rbind, lines.list)

  if (diagnostic_plot) {
    print(
      ggplot() +
        geom_sf(data = fronts, aes(colour = factor(year)), linetype = "dashed") +
        geom_sf(data = pts, aes(colour = factor(colonisation_year)), size = 1) +
        labs(colour = "Year") +
        theme_bw()
    )
  }

  # transform back to original crs; densify lines first (see score_colonised())
  pts <- st_transform(pts, crs = init.crs)
  fronts <- st_segmentize(fronts, dfMaxLength = 10000) |>
    st_transform(crs = init.crs)
  fp <- st_transform(fp, crs = init.crs)

  list(scored.points = pts, fronts = fronts, fp = fp)
}


