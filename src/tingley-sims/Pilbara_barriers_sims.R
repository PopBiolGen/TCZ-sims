File.ID <- "test"


#set the working directory

source("src/ABC/pprocess_functions.R")
load(file="dat/low density points.RData")
load("dat/Kernel_fits.RData")
plb<-read.csv("dat/art_nat_clp.csv")
load("dat/Posteriors.RData")

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

spread.table<-as.matrix(cbind(ID,X,Y,Pres,n.pairs,u,age),
nrow=length(age),ncol=7)



###################################################
# Does removing n NNs around the points stop toads?
nn<-seq(40, 110, 5)
reps<-5
output<-vector("list", length=length(nn)*nrow(ld.points)*reps) # vector to take outputs

for (kk in 1:length(nn)){ # for each number of NNs
  for(ii in 1:nrow(ld.points)){ # for each point
  	mod<-knock.out.nn.xy(X=ld.points[ii,"x"], Y=ld.points[ii,"y"], spread.table=spread.table, n=nn[kk], natural=nats)
  	target<-which(mod$spread.table[,"Pres"]==2)
	mod$spread.table[mod$spread.table[,"Pres"]==2,"Pres"]<-0
  	for (jj in 1:reps){ # for ten reps
      lambda.samp<-10^rnorm(1, mean=sample.lambda, sd=sample.lambda.sd)
  	  r.samp<-10^2
  	  temp<-spread.pilb(pop=mod$spread.table, gens=100, pairs=mod$pairs.mod, target=target, delta=lambda.samp, r=r.samp)
      temp<-c(temp, list(pars=cbind(lambda=lambda.samp, r=r.samp)), list(nns=nn[kk]), list(point=ii))
      output[[nrow(ld.points)*reps*(kk-1)+reps*(ii-1)+jj]]<-temp
    }
  }
}
save(output, file=paste("out/nn", File.ID, ".RData", sep=""))



