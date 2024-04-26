args=(commandArgs(TRUE))

#evaluate the arguments
# input arguments will be File.ID

for(i in 1:length(args)) {
   eval(parse(text=args[[i]]))
}

setwd("~/final_natural_water_files/Pilbara")
source("~/evo-dispersal/art_wbdies/VRD/pprocess_functions.R")

load("../Posteriors.RData")

plb<-read.csv("../art_nat_clp.csv")
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
load("../Kernel_fits.RData")

u<-fits[u,1:2]

spread.table<-as.matrix(cbind(ID,X,Y,Pres,n.pairs,u,age),
nrow=length(age),ncol=7)



#################################################
#can toads spread just using natural  and artificial waterbodies?
if (File.ID==1) save(spread.table, file="Natural and artificial waterbodies table.RData")
times<-c()
points.colonised<-c()
for (ii in 1:10){
  lambda.samp<-10^rnorm(1, mean=sample.lambda, sd=sample.lambda.sd)
  r.samp<-10^2
  temp<-spread.pilb(pop=spread.table, gens=100, pairs=pairs_pdist, target=target, delta=lambda.samp, r=r.samp)
  times<-c(times, temp[[1]])
  points.colonised<-rbind(points.colonised, cbind(repl=rep(ii, nrow(temp[[2]])), temp[[2]]))
}
save(times, points.colonised, file=paste("times colonised natural and artificial", File.ID, ".RData", sep=""))
################################################