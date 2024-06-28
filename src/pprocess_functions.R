library(sf)
library(tmap)
library(tmaptools)
library(magick)
library(dplyr)

#The kernel
dcncross<-function(x, u, v) {  #cauchy-normal distribution in 2D
  (x*u^v*v*sqrt(v^v*(u^2*v+x^2)^(-2-v)))/(2*pi*x)
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
make_plots <- function(scenario.name) {
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
    type = "osm", # for satellite image, "esri-imagery",
    zoom = 8,
    ext = 1
  )
  
  # replace with data to be mapped
  d <- read.point.data(paste0(in.name, ".csv"))
  
  p <- tm_shape(bm,
                unit = "km") +
    tm_rgb() +
    tm_shape(d) +
    tm_dots(size = 0.2,
            col = "arrival",
            breaks = 2024:2042, 
            legend.format = list(big.mark = ""),
            title = "Predicted year of toad arrival") 
  
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
      tm_layout(title = gsub("[^0-9]", "", year.name))
    tmap_save(p, filename = fname)
  }
  
  flist <- list.files(path = fpath, pattern = ".png")
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
  if (gc.switch) gc()
  names(outList) <- c("pairs_pdist", "spread.table")
  outList
}

# function that will run sims on whatever has been thrown into environment by setup()
run_sims <- function(n.sims = 100, gens, plot = FALSE, rollup) { 
  output<-vector("list", length=n.sims) # vector to take outputs
  
  for (rr in 1:n.sims){ # for reps
    cat("Rep ", rr, "\n")
    lambda.samp<-10^rnorm(1, mean=sample.lambda, sd=sample.lambda.sd)
    r.samp<-10^2
    temp<-spread.pilb(pop=spread.table, gens=gens, pairs=pairs_pdist, delta=lambda.samp, r=r.samp, plot = plot, rollup = rollup)
    temp<-c(temp, list(pars=cbind(lambda=lambda.samp, r=r.samp)))
    output[[rr]]<-temp
  }
output
}

save_outputs <- function(output, path, scenario.name, start.year, plot.time = FALSE, ABC = FALSE){
  ######## save output and generate summaries ########
  # time to arrive at target
  out.name <- paste0(path, "/", scenario.name) # make filename for scenario
  
  save(output, file = paste0(out.name, ".RData"))
  if (ABC) {
    return()
  }
  
  time.vec <- unlist(lapply(output, FUN = function(x){c(x$gen)}))
  if (plot.time & sum(is.finite(time.vec))>0) {
    pdf(file = paste0(out.name, ".pdf"))
      hist(time.vec, xlab = "Time to the Pilbara (y)")
    dev.off()
  }
  
  # Work out mean time to and arrival year for each point
  
  add_index_column <- function(x, index) { # Function to add index column to each matrix
    mat <- x$popmatrix[x$popmatrix[, "age"] > 0, ] # remove points never colonised
    index_col <- rep(index, n = nrow(mat))
    arrival <- round(start.year+max(mat[, "age"])-mat[, "age"]) # calculate arrival year
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
                  TCZ = FALSE # implement the TCZ, or not?
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
  
  u<-fits[u, c("u", "v", "max.dist", "area")]
  if (!trunc.dist) u[, "max.dist"] <- NULL # to switch to infinite positive bounds on kernel
  
  cat("Building spread table...\n")
  
  pData[[present.id]][is.na(pData[[present.id]])] <- 0 # set NAs to 0
  tg <- pData[[present.id]]==2 # identify target sites
  pData[[present.id]][pData[[present.id]]==2] <- 0 # re-set targetted sites to 0
  
  spread.table <-  cbind(ID = 1:nrow(pData),
                         X = pData[[X.id]],
                         Y = pData[[Y.id]],
                         target = tg,
                         u = u,
                         Pres = pData[[present.id]], # replace 2 from target with 0
                         age = pData[[present.id]], # set already colonised to age = 1
                         nats = as.numeric(pData[[artificial.natural.id]]==0),
                         obs = as.matrix(pData[observations])) # does nothing if NULL, else vector of column names
  
  
  if (remove_duplicates) {
    outList <- remove_spatial_duplicates(spread.table, threshold)
  }else {
    cat("Calculating pairwise distance matrix...\n")
    pairs_pdist<-pdist(X = spread.table[, "X"],Y = spread.table[, "Y"])
    
    outList <- list(spread.table = spread.table, pairs = pairs_pdist)
  }
  write.csv(spread.table, file = "out/basemap_points.csv") # for mapping later
  cat("Placing spread table and pairwise distance matrix in: ")
  list2env(outList, envir = globalenv())
}


# spreads the population over gens generations or until target sites are reached
# returns number of generations
# target is a vector of rows of pop that contain targets
spread.pilb<-function(pop, gens, pairs, delta, r, plot=FALSE, rollup){ #pairs is a list from pdist.fast  
  trigger <- TRUE # to catch time to first arrival in pilbara
  time.to.target <- NA
  # A progress bar
  pb <- txtProgressBar(min = 0, max = gens, style = 3)
  for (i in 1:gens){
		#which sites are occupied
		if (rollup) {
		  occp <- pop[,"Pres"]==1 & pop[, "age"] < 4 # for large simulations, can stop processing points 4+ y colonised
		}else {
		  occp <- pop[,"Pres"]==1 
		}
    #browser()
		lambda_t_x <- rpois(sum(occp), delta) # stochastic propagules from occupied site x time t
		gma <- sum(lambda_t_x) # total propagules at this time step
		pairs_t<-pairs[occp, , drop = FALSE] #collect relevant rows of pair matrix
		marg_dens <- apply(pairs_t, # get probability density accruing from each source
		                   MARGIN = 2, 
		                   FUN = dcncross.trunc, 
		                   u = pop[occp,"u"], 
		                   v = pop[occp,"v"],
		                   trunc.area = pop[occp,"area"],
		                   trunc.dist = pop[occp,"max.dist"])
		marg_dens <- sweep(marg_dens, # make into toad density
		                   MARGIN = 1, 
		                   STATS = lambda_t_x, 
		                   FUN = "*")
		marg_expected_n <- sweep(marg_dens, # make into expected count
		                         MARGIN = c(1,2), 
		                         STATS = pi*r^2, 
		                         FUN = "*")
		gc() # tidy up memory
		expected_n <- colSums(marg_expected_n, na.rm = TRUE) # sum contributions from all sources
		if (gma < .Machine$integer.max){ # when we go over machine tolerance, skip draw from multinom (realised likely to be very close to expected)
		  failures <- gma-sum(expected_n)
		  if(failures<0) failures<-0 # catches the approximately statement (primarily happens at large r)
		  expected_n <- c(expected_n, failures) # add failures
		  realised_n <- rmultinom(1, gma, expected_n)[-length(expected_n)] # draw propagules
		  colonised <- realised_n > 2
		}else colonised <- expected_n > 2
		pop[colonised, "Pres"] <- 1 #set to colonised
		occp <- pop[,"Pres"]==1 #which sites are occupied now
		pop[occp, "age"] <- pop[occp, "age"] + 1 # age each of the colonised populations
		if (plot) plotter(pop, file.name=paste("out/", i,".png", sep=""), gen=i)
		test.condition <- sum(pop[pop[, "target"]==1,"Pres"]) # number of target sites occupied
    if (test.condition > 0 && trigger) {
      time.to.target <- i #record time of arrival
      trigger <- FALSE
    }
		if (test.condition > 0 & test.condition == sum(pop[, "target"]==1)) break # stop if all target points colonised
		setTxtProgressBar(pb, i)
  }
  close(pb) # close progress bar
	list(gen=time.to.target, popmatrix=pop)	
}
