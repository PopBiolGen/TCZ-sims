#The kernel
dcncross<-function(x, u, v) {  #stuart's cauchy-normal distribution in 2D
  (x*u^v*v*sqrt(v^v*(u^2*v+x^2)^(-2-v)))/(2*pi*x)
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
  p.dist
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

# Finds points closer than threshold distance apart and removes one of the points
# returns a filtered population table
# to be used before creation of the spread table.
remove_spatial_duplicates <- function(pop.mat, threshold){
  X <- pop.mat[,"X"]
  Y <- pop.mat[, "Y"]
  pDists <- pdist(X, Y) # get pairwise distances
  pDists[lower.tri(pDists, diag = TRUE)] <- threshold # set relevant parts to >threshold
  below_threshold <- function(x){x < threshold} 
  tooClose <- apply(pDists, MARGIN = 1, FUN = below_threshold) # matrix
  tooClose <- apply(tooClose, 2, FUN = sum) == 0 # colsums == 0
  pop.mat[tooClose,]
}

# spreads the population over gen generations, and compares predictions to observed spread
spread<-function(pop, gens, pairs, delta, r, obs){ #pairs is a list from pdist.fast
  preds<-NULL	
for (i in 1:gens){
		#if (i%%5==0) output(pop, i, K)
		occp<-subset(pop, pop[,"Pres"]==1) #collect occupied sites
		occp<-cbind(occp, lambda=rpois(nrow(occp), delta))
		gma<-sum(occp[,"lambda"])
		if (length(occp[,1])==length(pop[,1])) {
			print(paste("No more vacant opportunities at generation", i))
			break()	
		}
		potl<-pairs[occp[,"ID"]] #collect relevant parts of pair list
		potl<-do.call("rbind", potl)
		src.ID<-rep(occp[,"ID"], times=occp[,"n.pairs"])
		potl<-cbind(src.ID, potl)
		lambda<-rep(occp[,"lambda"], times=occp[,"n.pairs"])
		U<-rep(occp[,"u"], times=occp[,"n.pairs"]) #expand source specific kernel parameters
		V<-rep(occp[,"v"], times=occp[,"n.pairs"])
		recruits<-lambda*dcncross(potl[,"dists"]+0.05, U, V) #calculate densities attributable to each pair
		recruits<-(pi*r^2*neigh.corr(pairs, r))/gma*tapply(recruits, potl[,"snk.ID"], sum) #sum densities from colonised waterbodies over all waterbodies and convert to proportion
		failures<-1-sum(recruits)
		if(failures<0) failures<-0 # catches the approximately statement (primarily happens at large r)
		recruits<-c(recruits, failures) #add on the failures
		recruit.ID<-as.integer(names(recruits)[-length(recruits)])
		recruits<-rmultinom(1, gma, recruits)[-length(recruits)]
		recruits<-cbind(recruit.ID, recruits)
		recruits<-subset(recruits, recruits[,"recruits"]>2)
		pop[match(recruits[,"recruit.ID"], pop[,"ID"]), "Pres"]<-1
		pop[which(pop[,"Pres"]==1), "age"]<-1+pop[which(pop[,"Pres"]==1), "age"]
		preds<-rbind(preds,output_val(pop, i))
	#plotter(pop, file.name=paste("gen",i,".png", sep=""))
	#print(obs_pred_cf(preds,obs))
	}
		obs_pred_cf(preds,obs)
}


# Sets up the spread table and pulls parameters ready for simulations
  # returns the spread table
setup <- function(point.data = "dat/art_nat_clp.csv", 
                  X.id = "POINT_X", 
                  Y.id = "POINT_Y", 
                  present.id = "ARRIVE_MCP",
                  artificial.natural.id = "art_nat",
                  rain.id = "rain_1mm",
                  remove_duplicates = TRUE,
                  threshold = 100, # metres within which to filter out duplicates
                  constant.rain = NULL # else number of days you want across whole area 
                  ){
  load("dat/Kernel_fits.RData")
  if (is.object(point.data)) { 
    pData <- point.data
  } else {
      pData<-read.csv(point.data)
    }
  
  load("dat/Posteriors.RData")
  #get matrix for the 'spread' function
  # need matrix containing:
  # "ID, X, Y, Pres (0s), n.pairs, u (rainy days*85.35[which is estimate of u]), 
  # age (0s)"
  
  ID <- 1:nrow(pData)
  X <- pData[[X.id]]
  Y <- pData[[Y.id]]
  Pres <- pData[[present.id]]
  target <- Pres==2
  Pres[Pres==2] <- 0
  age <- Pres # set already colonised to age = 1
  nats <- pData[[artificial.natural.id]]==0
  
  if (is.null(constant.rain)){
    u <- (pData[[rain.id]]-1)/364
    u <- 3*(u-u^2) + u^3
    u <- pData[[rain.id]]+3*pData[[rain.id]]*(1-u)
    u <- floor(u)
  }else {u <- rep(constant.rain, nrow(pData))}
  
  u<-fits[u,1:2]
  
  spread.table <-  cbind(ID, X, Y, Pres, target, u, age, nats)
  
  if (remove_duplicates) {
    spread.table <- remove_spatial_duplicates(spread.table, threshold)
  }
  
  pairs_pdist<-pdist(X = spread.table[, "X"],Y = spread.table[, "Y"])
  
  outList <- list(spread.table = spread.table, pairs = pairs_pdist)
  list2env(outList, envir = globalenv())
}


# spreads the population over gens generations or until target sites are reached
# returns number of generations
# target is a vector of rows of pop that contain targets
spread.pilb<-function(pop, gens, pairs, delta, r, plot=FALSE){ #pairs is a list from pdist.fast  
  trigger <- TRUE # to catch time to first arrival in pilbara
  for (i in 1:gens){
		occp <- pop[,"Pres"]==1 #which sites are occupied
		lambda_t_x <- rpois(sum(occp), delta) # stochastic propagules from occupied site x time t
		gma <- sum(lambda_t_x) # total propagules at this time step
		pairs_t<-pairs[occp, , drop = FALSE] #collect relevant rows of pair matrix
		marg_dens <- apply(pairs_t, # get probability density accruing from each source
		                   MARGIN = 2, 
		                   FUN = dcncross, 
		                   u = pop[occp,"u"], 
		                   v = pop[occp,"v"])
		marg_dens <- sweep(marg_dens, # make into toad density
		                   MARGIN = 1, 
		                   STATS = lambda_t_x, 
		                   FUN = "*")
		marg_expected_n <- sweep(marg_dens, # make into expected count
		                         MARGIN = c(1,2), 
		                         STATS = pi*r^2, 
		                         FUN = "*")
		expected_n <- colSums(marg_expected_n, na.rm = TRUE) # sum contributions from all sources
		failures <- gma-sum(expected_n)
		if(failures<0) failures<-0 # catches the approximately statement (primarily happens at large r)
		expected_n <- c(expected_n, failures) # add failures
		realised_n <- rmultinom(1, gma, expected_n)[-length(expected_n)] # draw propagules
		colonised <- realised_n > 2
		pop[colonised, "Pres"] <- 1 #set to colonised
		occp <- pop[,"Pres"]==1 #which sites are occupied now
		pop[occp, "age"] <- pop[occp, "age"] + 1 # age each of the colonised populations
		if (plot==TRUE) plotter(pop, file.name=paste(i,".png", sep=""), gen=i)
		test.condition <- sum(pop[pop[, "target"]==1,"Pres"]) # number of target sites occupied
    if (test.condition > 0 && trigger) {
      time.to.pilbara <- i #record time of arrival
      trigger <- FALSE
    }
		if (test.condition == sum(pop[, "target"]==1)) break # stop if all target points colonised
	}
	list(gen=time.to.pilbara, popmatrix=pop)	
}