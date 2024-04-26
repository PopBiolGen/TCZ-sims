args=(commandArgs(TRUE))

#evaluate the arguments
# input arguments will be File.ID

for(i in 1:length(args)) {
   eval(parse(text=args[[i]]))
}

#set the working directory
#setwd("/Users/ben/Desktop/Spread model stuff from Reid")
setwd("~/final_natural_water_files")
#source("/Users/ben/evo-dispersal/art_wbdies/VRD/pprocess_functions.R")
source("~/evo-dispersal/art_wbdies/VRD/pprocess_functions.R")

vrd<-read.csv("vrd_val_done_edit_csv_Events_WGS_rain_done_UTM_UNIQUE.csv")
#max(vrd$POINT_X)-min(vrd$POINT_X) # 221702.5
#max(vrd$POINT_Y)-min(vrd$POINT_Y) # 196911.4
pairs_pdist<-pdist.fast(X=vrd$POINT_X,Y=vrd$POINT_Y,maximum=300000,space.size= 300000)

##get matrix for Ben's 'spread' function

# need matrix containing:
# "ID, X, Y, Pres (0s), n.pairs, u (rainy days*85.35[which is estimate of u]), 
# age (0s)"

ID<-as.numeric(rownames(vrd))
X<-vrd$POINT_X
Y<-vrd$POINT_Y
Pres<-vrd$ARRIVE_MCP
age<-rep(0,length(X))

# calculate n.pairs using pdist
n.pairs<-do.call("c",lapply(pairs_pdist,nrow))

u<-(vrd$rain_1mm-1)/364
u<-3*(u-u^2) + u^3
u<-vrd$rain_1mm+3*vrd$rain_1mm*(1-u)
u<-floor(u)
load("Kernel_fits.RData")

# 
u<-fits[u,1:2]

spread.table<-as.matrix(cbind(ID,X,Y,Pres,n.pairs,u,age),
nrow=length(age),ncol=7)

##Load up a table of observations (= Letnic's observed data)
ob<-read.csv("vrd_val_done_edit_csv_Events_WGS_rain_done_UTM_no_06_07.csv")

out<-c()
for (ii in 1:5000){
  ## draw r and delta from uniform priors
  r<- runif(1, 1e1, 1e3)
  delta<-runif(1, 1e3, 1e6)
  ##Feed them all to the spread function for validation
  temp<-spread(spread.table, 2, pairs_pdist, delta, r, ob)
  temp<-cbind(delta, r, temp)
  out<-rbind(out, temp)
}

save(out, file=paste("~/final_natural_water_files/ABC_output/", File.ID, ".rdata", sep=""))
