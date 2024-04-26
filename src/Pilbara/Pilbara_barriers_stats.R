#Local script

#set the working directory
setwd("~/Dropbox/Papers/Submitted/Artificial waterbodies/Figures/Data")
source("~/evo-dispersal/art_wbdies/VRD/pprocess_functions.R")
load(file="low density points.RData")
load("Kernel_fits.RData")
plb<-read.csv("art_nat_clp.csv")
load("Posteriors.RData")

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
nats<-which(plb$art_nat==0)
arts<-which(plb$art_nat==1)


# calculate n.pairs using pdist
n.pairs<-do.call("c",lapply(pairs_pdist,nrow))

u<-(plb$rain_1mm-1)/364
u<-3*(u-u^2) + u^3
u<-plb$rain_1mm+3*plb$rain_1mm*(1-u)
u<-floor(u)
u<-fits[u,1:2]
artificial<-rep(0, length(X))
artificial[arts]<-1

spread.table<-as.matrix(cbind(ID,X,Y,Pres,n.pairs,u,age, artificial),
nrow=length(age),ncol=8)



###################################################
n<-100
art.points<-spread.table[arts, ]
included<-rep(0, nrow(art.points)) ; included[1:n]<-1
for(ii in 1:nrow(ld.points)){
  pdists<-sqrt((art.points[,"X"]-ld.points[ii, "x"])^2+(art.points[,"Y"]-ld.points[ii,"y"])^2)
  art.points<-cbind(art.points[order(pdists),], included)
}
nat.points<-cbind(spread.table[nats,], matrix(0, nrow=length(nats), ncol=3))


out.table<-rbind(art.points, nat.points)

pdf(file="../Barriersmap.pdf")
cls<-c("black", "grey40")
barr<-10
  for(ii in 1:nrow(ld.points)){
    plot(out.table[,"X"], out.table[,"Y"], pch=19, col=cls[out.table[,barr]+1])
    nats2<-!out.table[,"artificial"]
    points(out.table[nats2, "X"], out.table[nats2,"Y"], pch=19, col="darkorange")
    barr<-barr+1
  }
dev.off()

save(out.table, file="Data for Barriers map.RData")


