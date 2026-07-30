library(sf)
library(tmap)
library(terra)
library(tmaptools)
library(magick)
library(tidyverse)
library(Matrix)

# Control-type codes used by the (optional) control step in spread.pilb().
# CONTROL_NONE is the default for every point unless a scenario opts in via setup()'s
# control.type.id argument -- this keeps the control step a strict no-op for existing scripts.
CONTROL_NONE  <- 0L
CONTROL_TANK  <- 1L
CONTROL_FENCE <- 2L

# maps a character control-type column ("none"/"tank"/"fence", case-insensitive) onto the
# numeric codes above; unmatched or NA values fall back to CONTROL_NONE
encode_control_type <- function(x){
  code <- c(none = CONTROL_NONE, tank = CONTROL_TANK, fence = CONTROL_FENCE)
  out <- unname(code[tolower(as.character(x))])
  out[is.na(out)] <- CONTROL_NONE
  out
}

#The kernel
dcncross<-function(x, u, v) {  #cauchy-normal distribution in 2D
  (u^v*v*sqrt(v^v*(u^2*v+x^2)^(-2-v)))/(2*pi)
} 

# returns probability density for a truncated kernel
# Assumes trunc.area has been worked out during kernel fitting
dcncross.trunc <- function(x, u, v, trunc.dist, trunc.area){
  if (is.null(trunc.dist)) return(dcncross(x, u, v))
  density <- dcncross(x, u, v)/trunc.area
  density[x > trunc.dist] <- 0
  density
}

# given a point in the pdist list, takes out that point and its n nearest neighbours from the
#   pairwise distance list.  Returns modified list.
knock.out.nn<-function(pdist.list, point, n.n.neighb, natural){
  if (length(point)>1) warning("Multiple point removal not allowed")
  temp<-pdist.list[[point]] # get ID of points to knock out
  temp<-temp[order(temp[,'dists'])[1:(n.n.neighb+1)],"snk.ID"]
  temp<-temp[!temp%in%natural] #can't knock out natural points
  lapply(pdist.list, function(x) {subset(x, !x[,1]%in%temp)}) # and remove them from everywhere in the list
}

# given a point not in the spread table, takes out that point's n nearest (artificial) neighbours from the
#   spread table.  Returns modified spread.table and pdist.list.
knock.out.nn.xy<-function(X, Y, spread.table, n, natural){
  if (length(X)>1 | length(Y)>1) {warning("Multiple point removal not allowed"); return(NULL)}
  pdists<-sqrt((spread.table[,"X"]-X)^2+(spread.table[,"Y"]-Y)^2)
  top<-order(pdists)[!order(pdists)%in%natural]
  spread.table<-spread.table[-top[1:n],]
  pairs.mod<-pdist.fast(spread.table[,"X"], spread.table[,"Y"], maximum=500000, space.size=500000)
  spread.table[, "n.pairs"]<-do.call("c",lapply(pairs.mod,nrow))
  spread.table[, "ID"]<-1:nrow(spread.table)
  list(spread.table=spread.table, pairs.mod=pairs.mod)
}

# given a point, takes out that point and its n sized shortest path from the
#   pairwise distance list.  Returns modified list.
knock.out.path<-function(pdist.list, point, n.points=2, natural){
  temp<-pathway(pdist.list, point, n.points=2, natural)
  if (is.null(temp)) return(NULL)
  lapply(pdist.list, function(x) {subset(x, !x[,1]%in%temp)}) # and remove them from everywhere in the list
}

# function to make plots based on a scenario name
make_plots <- function(scenario.name, plot.year = TRUE, tcz.boundary = NULL) {
  in.name <- paste0("out/", scenario.name) # get filename for scenario
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
  
  # just to define mapping extent (varies among scenarios, otherwise)
  d <- read.point.data("out/basemap_points.csv")
  
  # read a basemap in from ESRI
  bm <- read_osm(
    d,
    type = "esri-topo", # for satellite image, "esri-imagery",
    zoom = 8,
    ext = 1
  )
  
  # replace with data to be mapped
  d <- read.point.data(paste0(in.name, ".csv"))
  
  # optionally transform tcz.boundary to Web Mercator to match basemap
  if (!is.null(tcz.boundary)) tcz.bm <- st_transform(tcz.boundary, 3857)
  
  p <- tm_shape(bm,
                unit = "km") +
    tm_rgb() +
    tm_shape(d) +
    tm_dots(size = 0.2,
            fill = "arrival",
            fill.scale = tm_scale_intervals(breaks = 2024:2042),
            fill.legend = tm_legend(title = "Predicted year of toad arrival"))
  
  if (!is.null(tcz.boundary)) p <- p + tm_shape(tcz.bm) + tm_borders(col = "red", lwd = 1.5)
  
  tmap_save(p, filename = paste0("out/year-of-arrival_", scenario.name, ".pdf"))
  
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
  input.flist <- list.files(path = fpath, pattern = "\\.csv$")
  for (yy in input.flist){
    temp <- read.point.data(fname = paste0(fpath, yy))
    year.name <- gsub("\\.csv$", "", yy)
    fname <- paste0(fpath, year.name, ".png")
    if (plot.year) y.name.plot <- year.name else y.name.plot <- ""
    p <- tm_shape(bm,
                  unit = "km") +
      tm_rgb() +
      tm_shape(temp) +
      tm_dots(size = 0.2,
              fill = "prob.colonised",
              fill.scale = tm_scale_intervals(breaks = seq(0, 1, length.out = 5)),
              fill.legend = tm_legend(title = "Probability of colonisation")) +
      tm_title(gsub("[^0-9]", "", y.name.plot))
    if (!is.null(tcz.boundary)) p <- p + tm_shape(tcz.bm) + tm_borders(col = "red", lwd = 1.5)
    tmap_save(p, filename = fname)
  }
  
  flist <- list.files(path = fpath, pattern = "\\.png$")
  
  if (length(flist) == 0) {
    warning("No PNG frames found in ", fpath, " — skipping animation. Were dynamic map CSVs written by save_outputs()?")
    return(invisible(NULL))
  }
  
  images <- image_read(paste0(fpath, flist))
  animation <- image_animate(images, fps = 1)
  image_write(animation, path = paste0("out/dynamic_maps/animated-map_", scenario.name, ".gif"))
  
  # clean up
  file.remove(paste0(fpath, input.flist))
  file.remove(paste0(fpath, flist))
}


# Approximate correction for pi*r2 calculation by the proportional overlap between waterbodies that are less than 2r distant from one another
neigh.corr<-function(pairs, r){
  #collect pairs less than 2r distant
  corrn<-function(x){
    ss <- x[,"dists"]<=2*r & x[,"dists"]>0
    temp<-matrix(x[ss,], ncol=2)
    if (length(temp)==0) return(1)
    theta<-2*acos(temp[,2]/(2*r))
    out<-1-(theta-sin(theta))/(2*pi)
    out<-prod(out)
    out
  }
  unlist(lapply(pairs, corrn))
}

obs_pred_cf<-function(preds, obs){
  preds[,5]<-preds[,5]+2007
  preds<-as.data.frame(preds)
  temp<-merge(preds, obs, by.x=c(2, 3, 5), by.y=c(9, 10, 15))
  #browser()
  temp<-temp[,"Pres"]==temp[,"OCCUPIED"]
  temp<-na.exclude(temp)
  sum(temp)
}

output_val<-function(pop, gens){
  pop<-pop[,1:4]
  pop<-cbind(pop,rep(gens, nrow(pop)))
  colnames(pop)[5]<-"Generation"	
  pop
}

# given a point, identifies n points making the shortest overall pathway between points
pathway<-function(pdist.list, point, n.points=2, natural){
  if (length(point)>1) {warning("Multiple point specification not allowed"); return(NULL)}
  IDs<-point
  temp<-pdist.list[[point]] # get first matrix
  for (ii in 1:(n.points-1)){
    temp<-temp[!temp[,"snk.ID"]%in%IDs,] #remove rows already identified in IDs
    next.id<-temp[which(temp[,"dists"]==min(temp[,"dists"]))[1], "snk.ID"]
    temp<-rbind(temp, pdist.list[[next.id]])
    IDs<-c(IDs, next.id)
  }
  if (sum(IDs%in%natural)>0) return(NULL)
  IDs
}


# Computes full distance matrix given vectors of X, and Y coordinates
pdist <- function(X, Y){
  if (length(X) < 10000){
    X<-X-min(X) # start coordinates at zero
    Y<-Y-min(Y)
    # function for calculating squared distances along each axis
    sq.dist <- function(x_1, x_2){
      (x_1-x_2)^2
    }
    # returns euclidean distance given squared distances along x and y (thanks Pythagoras)
    euc.dist <- function(sq.dist.x, sq.dist.y){
      sqrt(sq.dist.x + sq.dist.y)
    }
    s.d.x <- outer(X, X, FUN = sq.dist)
    s.d.y <- outer(Y, Y, FUN = sq.dist)
    p.dist <- euc.dist(s.d.x, s.d.y)
    return(p.dist)
  } else {
    return(as.matrix(dist(cbind(X, Y))))
  }
}

#plots opportunities and colonised populations
plotter<-function(popmatrix, file.name="temp.png", gen){
  png(filename=file.name, width=7, height=7, units="cm", res=150, pointsize=6)
  plot(popmatrix[,2], popmatrix[,3], xlab="False easting (kms)", ylab="False northing (kms)", pch=19)
  occp<-subset(popmatrix, popmatrix[,"Pres"]==1)
  points(occp[,2], occp[,3], pch=19, col="red")
  legend('topleft', legend=paste("Time =", gen), bty="n", pch=NA, cex=1.5)
  dev.off()
}

# function mapping days of rain to days of movement
rain_to_days <- function(rain.days, days.per.rain = 4){
  days <- function(r.d){
    denom <- 365:(365-(r.d-1))
    p.no.move <- 1-days.per.rain/denom
    p.no.move <- prod(p.no.move)
    p.move <- 1 - p.no.move
    days <- p.move * 365
    round(days)
  }
  sapply(rain.days, days)
}

# Finds points closer than threshold distance apart and removes one of the points
# returns a filtered population table, and a filtered pairwise distance table
remove_spatial_duplicates <- function(pop.mat, threshold){
  X <- pop.mat[,"X"]
  Y <- pop.mat[, "Y"]
  gc.switch <- ifelse(length(X) > 10000, TRUE, FALSE) # do we need to manage memory?
  outList <- vector(mode = "list", length = 2) #list to take outputs
  cat("Calculating pairwise distance matrix...\n")
  outList[[1]] <- pdist(X, Y) # get pairwise distances
  if (gc.switch) gc() # free up memory
  cat("Finding spatial duplicates...\n")
  cat("\t Thresholding...\n")
  pd <- outList[[1]] < threshold # logical matrix for what comes next..
  cat("\t Removing symmetry...\n")
  pd[lower.tri(pd, diag = TRUE)] <- FALSE # set lower triangle to FALSE
  if (gc.switch) gc() # free up memory from lower.tri
  cat("\t Finding duplicates...\n")
  tooClose <- apply(pd, 2, FUN = sum) == 0 # colsums == 0
  if (gc.switch) rm(pd); gc() # free up memory from lower.tri
  cat("Removing spatial duplicates...\n")
  outList[[1]] <- outList[[1]][tooClose, tooClose] #subset pairwise matrix
  outList[[2]] <- pop.mat[tooClose,]
  names(outList) <- c("pairs", "spread.table")
  outList
}

# function that will run sims on whatever has been thrown into environment by setup()
run_sims <- function(n.sims = 100, gens, plot = FALSE, rollup,
                      control.step = FALSE, stop.on.target = TRUE, drying.matrix = NULL) {
  if (control.step && rollup) stop("rollup and control.step cannot both be TRUE -- control-step scenarios must process every site every generation")
  output<-vector("list", length=n.sims) # vector to take outputs

  for (rr in 1:n.sims){ # for reps
    cat("Rep ", rr, " of ", n.sims, "\n")
    lambda.samp<-10^rnorm(1, mean=sample.lambda, sd=sample.lambda.sd)
    r.samp<-10^2
    temp<-spread.pilb(pop=spread.table, gens=gens, pairs=pairs, delta=lambda.samp, r=r.samp, plot = plot, rollup = rollup,
                       control.step = control.step, stop.on.target = stop.on.target, drying.matrix = drying.matrix)
    temp<-c(temp, list(pars=cbind(lambda=lambda.samp, r=r.samp)))
    output[[rr]]<-temp
    gc() # cleanup memory
  }
output
}

save_outputs <- function(output, path, scenario.name, start.year, plot.time = FALSE, ABC = FALSE, extinction.aware = FALSE){
  ######## save output and generate summaries ########
  # time to arrive at target
  out.name <- paste0(path, "/", scenario.name) # make filename for scenario
  
  save(output, file = paste0(out.name, ".RData"))
  if (ABC) {
    return()
  }
  write.csv(output[[1]]$popmatrix, file = "out/basemap_points.csv", row.names = FALSE)
  
  time.vec <- unlist(lapply(output, FUN = function(x){c(x$gen)}))
  if (plot.time & sum(is.finite(time.vec))>0) {
    pdf(file = paste0(out.name, ".pdf"))
      boxplot(time.vec, xlab = "Time to target (y)", horizontal = TRUE)
    dev.off()
  }
  
  # Work out mean time to and arrival year for each point
  
  add_index_column <- function(x, index) { # Function to add index column to each matrix
    mat <- x$popmatrix[x$popmatrix[, "age"] > 0, ] # remove points never colonised
    index_col <- rep(index, n = nrow(mat))
    if (extinction.aware) {
      # first.col.gen is a direct record of first-colonisation generation, immune to the gaps in
      # `age` that occur once a point can go extinct and later recolonise
      arrival <- round(start.year + mat[, "first.col.gen"])
    } else {
      arrival <- round(start.year+max(mat[, "age"])-mat[, "age"]) # calculate arrival year
    }
    cbind(index_col, mat, arrival)
  }
  # Apply the function to each element of the list using lapply
  modified_matrices <- lapply(seq_along(output), function(i) {
    add_index_column(output[[i]], i)
  })
  # bind the lot together into single matrix
  pop.out <- do.call("rbind", modified_matrices)
  
  # get mean arrival time for each point, for making a static map
  pop.summary <- pop.out %>% 
    as.data.frame() %>%
    group_by(ID) %>%
    summarise_all(mean)
  write.csv(pop.summary, file = paste0(out.name, ".csv"), row.names = FALSE) 
  
  # to make a dynamic map...
  # for each year, 2024 to max(mean arrival time), generate a csv to plot, that reports probability of colonisation at that time for each waterpoint
  max.time <- max(pop.summary$arrival) + 2
  for (yy in start.year:max.time){
    fname <- paste0(path, "/dynamic_maps/", scenario.name, "_", yy, ".csv")
    pop.summary <- pop.out %>%
      as.data.frame() %>%
      group_by(ID) %>%
      summarise(X = mean(X), Y = mean(Y), prob.colonised = mean(arrival <= yy)) %>%
      filter(prob.colonised > 0.1)
    write.csv(pop.summary, file = fname, row.names = FALSE) 
  }
}


# Sets up the spread table and pulls parameters ready for simulations
  # returns the spread table
setup <- function(point.data = "dat/art_nat_clp.csv", 
                  X.id = "POINT_X", 
                  Y.id = "POINT_Y", 
                  present.id = "ARRIVE_MCP",
                  artificial.natural.id = "art_nat",
                  rain.id = "rain_1mm",
                  observations = NULL,
                  remove_duplicates = TRUE,
                  threshold = 100, # metres within which to filter out duplicates
                  constant.rain = NULL, # else number of days you want across whole area
                  trunc.dist = TRUE, # false for full kernel
                  TCZ = FALSE, # implement the TCZ, or not?
                  control.type.id = NULL, # name of a character column ("none"/"tank"/"fence") -- enables the control step
                  fence.area.id = NULL, # name of a column giving fenced-polygon area (m^2, Albers)
                  fence.units.id = NULL, # name of a column giving km-of-fenceline + gate count ("g" in the qmd)
                  fence.breach.rate = 0.05, # assumed per-km-of-fenceline/gate breach rate
                  nat.surv.prob.id = NULL, # name of a column giving per-point natural survival probability (1 - drying prob)
                  control.fail.prob.id = NULL # name of a column overriding the default annual control-defeat probability
){
  cat("Loading kernel parameters...\n")
  load("dat/Kernel-fits_truncated.RData")
  
  cat("Loading point data...\n")
  if (is.object(point.data)) { 
    pData <- point.data
  } else {
    pData <- read.csv(point.data)
  }
  
  if (TCZ) pData <- subset(pData, !(pData[["TCZ"]] == 1 & pData[[artificial.natural.id]] == 0))
  
  cat("Loading posterior estimates...\n")
  load("dat/Posteriors.RData")
  
  
  cat("Assigning kernel parameters to waterpoints...\n")
  # assign kernel values to waterpoints
  if (is.null(constant.rain)){
    u <- rain_to_days(pData[[rain.id]])
  }else {u <- rep(constant.rain, nrow(pData))}
  #browser()
  u<-fits[u, c("u", "v", "max.dist", "area")]
  if (!trunc.dist) u[, "max.dist"] <- NULL # to switch to infinite positive bounds on kernel
  
  cat("Building spread table...\n")
  
  pData[[present.id]][is.na(pData[[present.id]])] <- 0 # set NAs to 0
  tg <- pData[[present.id]]==2 # identify target sites
  pData[[present.id]][pData[[present.id]]==2] <- 0 # re-set targetted sites to 0

  nats_vec <- as.numeric(pData[[artificial.natural.id]]==0)

  cat("Assigning control-step parameters to waterpoints...\n")
  n.pts <- nrow(pData)
  # control type: defaults to CONTROL_NONE for everyone unless control.type.id is supplied;
  # natural points can never carry an artificial control type
  control.type <- if (!is.null(control.type.id)) encode_control_type(pData[[control.type.id]]) else rep(CONTROL_NONE, n.pts)
  control.type[!is.na(nats_vec) & nats_vec==1] <- CONTROL_NONE

  # fenced-area effective detection radius (circle of equal area to the fenced polygon)
  fence.area <- if (!is.null(fence.area.id)) pData[[fence.area.id]] else rep(NA_real_, n.pts)
  r.eff <- ifelse(control.type==CONTROL_FENCE & !is.na(fence.area), sqrt(fence.area/pi), NA_real_)

  # fence breach probability, from km-of-fenceline + gate count ("g" in the qmd)
  fence.units <- if (!is.null(fence.units.id)) pData[[fence.units.id]] else rep(0, n.pts)
  fence.units[is.na(fence.units)] <- 0
  p.breach <- ifelse(control.type==CONTROL_FENCE, 1-(1-fence.breach.rate)^(fence.units/4), 0)

  # natural-point survival probability (1 - drying probability); static default of 1 (always survives)
  # until a real per-point/per-year drying model is wired in
  surv.prob <- if (!is.null(nat.surv.prob.id)) pData[[nat.surv.prob.id]] else rep(1, n.pts)
  surv.prob[is.na(surv.prob)] <- 1

  # annual control-defeat probability for tanks/fences, overridable per-row
  control.fail.prob.default <- ifelse(control.type==CONTROL_TANK, 0.01,
                                ifelse(control.type==CONTROL_FENCE, 0.05, 0))
  control.fail.prob <- if (!is.null(control.fail.prob.id)) {
    cfp <- pData[[control.fail.prob.id]]
    ifelse(is.na(cfp), control.fail.prob.default, cfp)
  } else control.fail.prob.default

  # bookkeeping columns, always present regardless of control-step usage
  occ.years <- pData[[present.id]] # mirrors age's initial value: 1 for already-colonised points, 0 otherwise
  first.col.gen <- ifelse(pData[[present.id]]==1, 0, NA_real_)

  spread.table <-  cbind(ID = 1:nrow(pData),
                         X = pData[[X.id]],
                         Y = pData[[Y.id]],
                         target = tg,
                         u = u,
                         Pres = pData[[present.id]], # replace 2 from target with 0
                         age = pData[[present.id]], # set already colonised to age = 1
                         nats = nats_vec,
                         control.type = control.type,
                         r.eff = r.eff,
                         fence.units = fence.units,
                         p.breach = p.breach,
                         surv.prob = surv.prob,
                         control.fail.prob = control.fail.prob,
                         occ.years = occ.years,
                         first.col.gen = first.col.gen)
  
  if (!is.null(observations)){ # does nothing if NULL, else vector of column names
    spread.table <- cbind(spread.table, obs = as.matrix(pData[observations]))
  }
  
  if (remove_duplicates) {
    outList <- remove_spatial_duplicates(spread.table, threshold)
  }else {
    cat("Calculating pairwise distance matrix...\n")
    pairs_pdist<-pdist(X = spread.table[, "X"],Y = spread.table[, "Y"])
    
    outList <- list(spread.table = spread.table, pairs = pairs_pdist)
  }
  gc() # free up whatever memory we can before embarking on this next step
  cat("Calculating sparse dispersal matrices... \n")
  d_mat <- dcncross(outList$pairs, 
                            u = outList$spread.table[,"u"], 
                            v = outList$spread.table[,"v"])/outList$spread.table[, "area"]
  d_mat[outList$pairs > outList$spread.table[, "max.dist"]] <- 0L
  outList$pairs <- Matrix(d_mat, sparse = TRUE) 
  
  cat("Placing spread table and dispersal matrix in: ")
  list2env(outList, envir = globalenv())
}



# spreads the population over gens generations or until target sites are reached
# returns number of generations
# target is a vector of rows of pop that contain targets
spread.pilb<-function(pop, gens, pairs, delta, r, plot=FALSE, rollup,
                       control.step=FALSE, # opt-in within-timestep control/extinction step
                       stop.on.target=TRUE, # FALSE to keep running the full horizon after a target breach (for point-years accounting)
                       drying.matrix=NULL # optional gens x n matrix of natural drying probabilities, refreshed into surv.prob each generation
                       ){ #pairs is a list from pdist.fast
  if (control.step && rollup) stop("rollup and control.step cannot both be TRUE -- control-step scenarios must process every site every generation")
  if (control.step){
    required.cols <- c("control.type", "r.eff", "fence.units", "p.breach", "surv.prob", "control.fail.prob", "occ.years", "first.col.gen")
    missing.cols <- setdiff(required.cols, colnames(pop))
    if (length(missing.cols) > 0) stop("control.step=TRUE requires spread.table columns built via setup()'s control-step arguments; missing: ", paste(missing.cols, collapse=", "))
  }
  trigger <- TRUE # to catch time to first arrival in pilbara
  time.to.target <- NA
  # A progress bar
  pb <- txtProgressBar(min = 0, max = gens, style = 3)

  # point-specific detection radius: fenced areas use r.eff (inflated to the fenced polygon's
  # effective radius), everything else falls back to the global scalar r
  r_vec <- if ("r.eff" %in% colnames(pop)) ifelse(is.na(pop[,"r.eff"]), r, pop[,"r.eff"]) else r

  if (control.step){ # static per-point masks, computed once (control.type/nats never change during the run)
    is.fence <- pop[,"control.type"]==CONTROL_FENCE
    is.tank.or.fence <- pop[,"control.type"] %in% c(CONTROL_TANK, CONTROL_FENCE) # tank and fence resolve identically at the control step (see below)
    is.nat   <- pop[,"nats"]==1
  }

  # make full dispersal matrix under normal conditions
  #d_mat <- dcncross(pairs, u = pop[,"u"], v = pop[,"v"])/pop[, "area"]
  #d_mat[pairs > pop[, "max.dist"]] <- 0
  #d_mat <- Matrix(d_mat, sparse = TRUE) # cast across to sparse matrix

  for (i in 1:gens){
    pres.before <- pop[,"Pres"]==1 # snapshot at generation start, used by the control step below
    #which sites are occupied
    if (rollup) {
      occp <- pop[,"Pres"]==1 & pop[, "age"] < 4 # for large simulations, can stop processing points 4+ y colonised
    }else {
      occp <- pop[,"Pres"]==1
    }
    lambda_t_x <- rpois(sum(occp), delta) # stochastic propagules from occupied site x time t
    gma <- sum(lambda_t_x) # total propagules at this time step
    pairs_t<-pairs[occp, , drop = FALSE] #collect relevant rows of pairwise dispersal matrix
    expected_n <- drop(crossprod(pairs_t, lambda_t_x))*(pi*r_vec^2) # sum density contributions from all sources and convert to expected n
    if (gma == 0){ # no propagules produced this generation (rmultinom errors on an all-zero probability vector)
      colonised <- rep(FALSE, nrow(pop))
    } else if (gma < .Machine$integer.max){ # when we go over machine tolerance, skip draw from multinom (realised likely to be very close to expected)
      failures <- gma-sum(expected_n)
      if(failures<0) failures<-0 # catches the approximately statement (primarily happens at large r)
      expected_n <- c(expected_n, failures) # add failures
      realised_n <- rmultinom(1, gma, expected_n)[-length(expected_n)] # draw propagules
      if (control.step){ # fence's own stopping power thins arrivals at not-yet-colonised fenced points
        fence.uncol <- is.fence & !pres.before
        if (any(fence.uncol)) realised_n[fence.uncol] <- rbinom(sum(fence.uncol), size=realised_n[fence.uncol], prob=1-pop[fence.uncol,"p.breach"])
      }
      colonised <- realised_n > 2
    }else {
      if (control.step){
        fence.uncol <- is.fence & !pres.before
        if (any(fence.uncol)) expected_n[fence.uncol] <- expected_n[fence.uncol]*(1-pop[fence.uncol,"p.breach"])
      }
      colonised <- expected_n > 2
    }
    pop[colonised, "Pres"] <- 1 #set to colonised

    if (control.step){ # control step: colonised points of a controlled type may be knocked back to
      # uncolonised this generation; every occupied point of a controlled type is re-tested fresh
      # each generation, regardless of when it was first colonised
      occ.now <- pop[,"Pres"]==1
      nat.idx <- is.nat & occ.now
      if (any(nat.idx)){
        if (!is.null(drying.matrix)) pop[nat.idx, "surv.prob"] <- 1 - drying.matrix[i, nat.idx]
        survives <- rbinom(sum(nat.idx), 1, pop[nat.idx, "surv.prob"])
        pop[which(nat.idx)[survives==0], "Pres"] <- 0
      }
      # control.fail.prob is P(control fails, i.e. toads persist) -- so P(extinguished) = 1-control.fail.prob.
      # Tank and fence resolve identically here (only their route into colonisation differs: fences
      # are breach-thinned at the colonisation step above), so both are drawn in one rbinom() call
      # over their per-row control.fail.prob values.
      control.idx <- is.tank.or.fence & occ.now
      if (any(control.idx)){
        persists <- rbinom(sum(control.idx), 1, pop[control.idx, "control.fail.prob"])
        pop[which(control.idx)[persists==0], "Pres"] <- 0
      }
    }

    occp <- pop[,"Pres"]==1 #which sites are occupied now (post-control, when applicable)
    pop[occp, "age"] <- pop[occp, "age"] + 1 # age each of the colonised populations
    # occ.years/first.col.gen are updated unconditionally (not gated on control.step): with no
    # extinction ever occurring, occ.years stays numerically identical to age, and first.col.gen
    # simply records each point's one-and-only colonisation generation
    pop[occp, "occ.years"] <- pop[occp, "occ.years"] + 1 # cumulative point-years occupied, gap-tolerant across extinction/recolonisation
    first.col.new <- occp & !pres.before & is.na(pop[,"first.col.gen"])
    pop[first.col.new, "first.col.gen"] <- i # stamped once, never overwritten on later recolonisation
    if (plot) plotter(pop, file.name=paste("out/", i,".png", sep=""), gen=i)
    test.condition <- sum(pop[pop[, "target"]==1,"Pres"]) # number of target sites occupied
    if (test.condition > 0 && trigger) {
      time.to.target <- i #record time of arrival
      trigger <- FALSE
      if (stop.on.target) break # stop as soon as any target sites are hit
    }
    if (stop.on.target && test.condition > 0 & test.condition == sum(pop[, "target"]==1)) break # stop if all target points colonised
    setTxtProgressBar(pb, i)
  }
  close(pb) # close progress bar
  list(gen=time.to.target, popmatrix=pop)
}
