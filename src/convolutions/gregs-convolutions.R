source("src/convolutions/convolution-functions.R")

#library(adehabitat)
library(amt)
library(dplyr)

ndays<-250 #days over which to convolve

#load data and convert date to POSIXct
d<-read.table("dat/gregsData.txt", header=T, sep="\t")
d$date<-as.POSIXct(d$date*24*60*60, origin=as.Date("1904-1-1"))

#Grab only data from Jan-March 2005/2006
d<-subset(d, d$month<=3)

# make into track object
d.track <- d %>% make_track(.x = EW.UTM,  # make into track object
                           .y = NS.UTM, 
                           .t = date, 
                           id = Toad.ID, 
                           crs = 32752)

# store individuals in column list format
d.track.cl <- d.track %>% nest(data = -"id") 

# create new data column, cast to steps, unlist and calculate step length per day
d.steps <- d.track.cl %>% mutate(steps = map(data, steps)) %>% 
                          unnest(steps) %>%
                          mutate(sl.dt = sl_/as.numeric(dt_)*86400)


# subset data and resample for each individual with more than 5 observations
resamps<-c()
ids <- unique(d.steps$id)
for (i in 1:length(ids)){
  cat(i, "\n")
	temp <- subset(d.steps, d.steps$id==ids[i])
	if (sum(is.finite(temp$ta_))<5) next
	sclr <- temp$sl.dt[is.finite(temp$sl.dt)]
	ta_ <- temp$ta_[is.finite(temp$ta_)]
	#simulate 1000 random walks by resampling dist and turn angle
	temp2 <- nday(sclr, ndays, ta_, 10000)
	resamps <- rbind(resamps, temp2)
}

nrow(resamps)
max.dist <- apply(resamps, 2, max) # maximum from each resampling of number of steps
max.dist <- max.dist * 1.1 # increase by 10% on resampled limit

# fit kernels 
fits<-nwise(resamps+0.05, init.v=1.5)

# cbind with truncation distance
fits <- cbind(fits[, !grepl(pattern = "LL", colnames(fits))], max.dist = max.dist)

# add area correction for truncation
fits <- cbind(fits, area = kernel.truncation(u = fits[,"u"], v = fits[,"v"], trunc.dist = fits[,"max.dist"]))



save(fits, file="dat/Kernel-fits_truncated.RData")

