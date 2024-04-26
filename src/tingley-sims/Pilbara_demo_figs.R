# This script runs the model and makes demo figures (without a base map)

source("src/pprocess_functions.R")

load("dat/Posteriors.RData")
load(file="dat/low density points.RData")
plb<-read.csv("dat/art_nat_clp.csv")

#max(plb$POINT_X)-min(plb$POINT_X) # 436252.5
#max(plb$POINT_Y)-min(plb$POINT_Y) # 318785.0
pairs_pdist<-pdist.fast(X=plb$POINT_X,Y=plb$POINT_Y,maximum=500000,space.size= 500000)

##get matrix for Ben's 'spread' function

# need matrix containing:
# "ID, X, Y, Pres (0s), n.pairs, u (rainy days*85.35[which is estimate of u]), 
# age (0s)"

ID<-as.numeric(rownames(plb))
X<-plb$POINT_X
Y<-plb$POINT_Y
Pres<-plb$ARRIVE_MCP
age<-rep(0,length(X))
target<-which(Pres==2)
Pres[Pres==2]<-0
nats<-which(plb$art_nat==0)
arts<-which(plb$art_nat==1)


# calculate n.pairs using pdist
n.pairs<-do.call("c",lapply(pairs_pdist,nrow))

u<-(plb$rain_1mm-1)/364
u<-3*(u-u^2) + u^3
u<-plb$rain_1mm+3*plb$rain_1mm*(1-u)
u<-floor(u)
load("dat/Kernel_fits.RData")

u<-fits[u,1:2]

spread.table<-as.matrix(cbind(ID,X,Y,Pres,n.pairs,u,age),
nrow=length(age),ncol=7)

plotter(spread.table, file.name="out/0.png", gen=0)
#################################################
setwd("out")
  lambda.samp<-10^sample.lambda
  r.samp<-10^2
  temp<-spread.pilb(pop=spread.table, gens=100, pairs=pairs_pdist, target=target, delta=lambda.samp, r=r.samp, plot=TRUE)
  
  temp$popmatrix[,"Pres"]<-1
  plotter(temp$popmatrix, file.name=paste(temp$gen+1, ".png", sep=""), gen=temp$gen+1)

################################################

for(ii in 1:nrow(ld.points)){ # for each point
  	mod<-knock.out.nn.xy(X=ld.points[ii,"x"], Y=ld.points[ii,"y"], spread.table=spread.table, n=100, natural=nats)
  	target<-which(mod$spread.table[,"Pres"]==2)
	  mod$spread.table[mod$spread.table[,"Pres"]==2,"Pres"]<-0
	  plotter(mod$spread.table, file.name="0.png", gen=0)
  	temp<-spread.pilb(pop=mod$spread.table, gens=100, pairs=mod$pairs.mod, target=target, delta=lambda.samp, r=r.samp, plot=TRUE)
}
  
setwd("..")