setwd("~/Dropbox/Papers/Submitted/Artificial waterbodies/Data")

source("convolution functions.R")

library(MASS)
library(adehabitat)

ndays<-200 #days over which to convolve

#load data and convert date to POSIXct
d<-read.table("Gregs data.txt", header=T, sep="\t")
d$date<-as.POSIXct(d$date*24*60*60, origin=as.Date("1904-1-1"))

#convert to ltraj

#Grab only data from Jan-March 2005/2006
d<-subset(d, d$month<=3)
#d<-subset(d, d$year<=2006)

#convert to ltraj
temp<-as.ltraj(d[,c(9, 8)], d$date, d$Toad.ID)
temp2<-ltraj2traj(temp)

# summ<-summary(temp)
# mean(summ$date.end-summ$date.begin)/(60*60*24)
# mean(summ$nb.reloc)

#grab dists and turn angles for simulating random walks
out<-data.frame(as.character(temp2$id), temp2$dist, temp2$rel.angle)
colnames(out)<-c("ID", "dist", "rel.angle")

# subset data and resample for each individual with more than 5 observations
resamps<-c()
for (i in 1: length(levels(out$ID))){
	temp<-subset(out, out$ID==levels(out$ID)[i])
	if (length(which(is.na(temp$rel.angle)==F))<5) next
	sclr<-temp$dist[which(is.na(temp$dist)==F)]
	rel.angle<-temp$rel.angle[which(is.na(temp$rel.angle)==F)]
	#simulate 1000 random walks by resampling dist and turn angle
	temp2<-nday(sclr, ndays, rel.angle, 1000)
	resamps<-rbind(resamps, temp2)
}

## Alternative resampling procedure that ignores individual differences
#   out<-na.omit(out)
#   sclr<-out$dist
#   rel.angle<-out$rel.angle
#   #simulate 1000 random walks by resampling dist and turn angle
#   resamps<-nday(sclr, ndays, rel.angle, 100000)

#look at resulting kernel for ndays days of movement 
hist(resamps[,ndays])

#save samples for ndays=1:208 (big file!)
save(resamps, file="Greg's convolution resamples.RData")

#fit kernel
fits<-nwise(resamps+0.05, init.v=1.5)

save(fits, file="Kernel_fits.RData")

