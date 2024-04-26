setwd("/Users/ben/Documents/Papers/Artificial Waterbodies/VRD")

a<-read.csv("all_watercourses_springs_lakes_holes_bores_CLIP_GDA94_RAIN_UTM.csv")

##get matrix for Ben's 'spread' function

# need matrix containing:
#"ID, X, Y, Pres (0s), n.pairs, u (rainy days*85.35), age (0s)"

ID<-as.numeric(rownames(a))
X<-a$POINT_X
Y<-a$POINT_Y
Pres<-rep(0,length(X))
age<-rep(0,length(X))

# calculate n.pairs using pdist
load("vrd_pairs.RData")
n.pairs<-do.call("c",lapply(pairs.art,nrow))

u<-(a$RASTERVALU-1)/364
u<-3*(u-u^2) + u^3
u<-a$RASTERVALU+3*a$RASTERVALU*(1-u)
u<-floor(u)
load("Kernel fits.RData")

# ndays> 100 in all cases, so only gaussian kernel necessary...
u<-fits[u,1:2]

spread.table<-as.matrix(cbind(ID,X,Y,Pres,n.pairs,u,age),nrow=length(age),ncol=7)

save(spread.table,file="spread.table.RData")



