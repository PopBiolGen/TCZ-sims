#local script

source("src/ABC/pprocess_functions.R")


alpha<-2.113e+10 #best fit value of K from VRD analysis
alpha.sd<-4016222910
plb<-read.csv("dat/art_nat_clp.csv")
#max(plb$POINT_X)-min(plb$POINT_X) # 436252.5
#max(plb$POINT_Y)-min(plb$POINT_Y) # 318785.0
pairs_pdist<-pdist.fast(X=plb$POINT_X,Y=plb$POINT_Y,maximum=100000,space.size= 500000)

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


###############################################################
## Heuristic - density plot.
# put a smoother through, calculate densities and make plots
png("out/heuristic.png", width=14, height=24, units="cm", res=300)
par(mfrow=c(2,1))
sm<-lowess(spread.table[,"X"], spread.table[,"Y"], f=0.2)
bw.adj<-0.2
dns<-density.line(X1=sm$x, Y1=sm$y, X2=spread.table[arts,"X"], Y2=spread.table[arts,"Y"], adjust=bw.adj)
dns.nat<-density.line(X1=sm$x, Y1=sm$y, X2=spread.table[nats,"X"], Y2=spread.table[nats,"Y"], adjust=bw.adj)
dns.comb<-density.line(X1=sm$x, Y1=sm$y, X2=spread.table[,"X"], Y2=spread.table[,"Y"], adjust=bw.adj)
plot(dns, type="l", bty="l", ylab="Density of waterbodies (per transect metre)", xlab="Easting (m)", ylim=c(0, max(c(dns$y, dns.nat$y, dns.comb$y))))
lines(dns.nat, type="l", bty="l", col="darkorange")
lines(dns.comb, type="l", bty="l", col="green")
legend("topright", legend=c("Natural", "Artificial", "Combined"), col=c("darkorange", "black", "green"), lty=1, bty="n", title="Waterbody type")
root1<-which(dns.nat$y==min(dns.nat$y[which(dns$x>-1300000 & dns$x<(-1200000))]))
root2<-which(dns.nat$y==min(dns.nat$y[which(dns$x>-1200000 & dns$x<(-1100000))]))
root3<-which(dns.nat$y==min(dns.nat$y[which(dns$x>-1100000 & dns$x<(-1000000))]))
arrows(dns$x[root1], 0, dns$x[root1], max(dns$y), length=0, col="red")
arrows(dns$x[root2], 0, dns$x[root2], max(dns$y), length=0, col="red")
arrows(dns$x[root3], 0, dns$x[root3], max(dns$y), length=0, col="red")
mtext(side=3, line=1, adj=0, text="A)")

plot(spread.table[,"X"], spread.table[,"Y"], pch=19, xlab="Easting (m)", ylab="Northing (m)", bty="l")
points(spread.table[nats,"X"], spread.table[nats,"Y"], pch=19, col="darkorange")
lines(sm, col="red")
points(sm$x[c(root1, root2, root3)], sm$y[c(root1, root2, root3)], col="red", pch=19, cex=2)
text(sm$x[c(root1, root2, root3)], sm$y[c(root1, root2, root3)], labels=c(1:3))
legend("topleft", legend=c("Natural", "Artificial", "Low density points"), col=c("darkorange", "black", "red"), pch=19, bty="n")
mtext(side=3, line=1, adj=0, text="B)")
dev.off()

out.data<-cbind(X=spread.table[arts,"X"], Y=spread.table[arts,"Y"], Category=rep("Artificial", length(arts)))
out.data<-rbind(out.data, cbind(spread.table[nats,"X"], spread.table[nats,"Y"], Category=rep("Natural", length(nats))))
out.data<-rbind(out.data, cbind(sm$x[c(root1, root2, root3)], sm$y[c(root1, root2, root3)], c("Eastern", "Central", "Western")))
save(out.data, file="out/HeurB.RData")

##############################################################
ld.points<-cbind(x=sm$x[c(root1, root2, root3)], y=sm$y[c(root1, root2, root3)])
save(ld.points, file="out/low density points.RData")
###################################################
# Does removing n NNs around the points stop toads?
nn<-seq(50, 75, 5)
for (kk in 1:length(nn)){
  print(nn[kk])
  results<-c()
  for (jj in 1:100){
    rep<-c()
    for(ii in 1:nrow(ld.points)){
      PseudoK<-rnorm(1, mean=21753929728, sd=0.5*9669058589)
      mod<-knock.out.nn.xy(X=ld.points[ii,"x"], Y=ld.points[ii,"y"], spread.table=spread.table, n=nn[kk], natural=nats)
      temp<-spread.pilb(pop = mod$spread.table,
                        gens = 200,
                        pairs = mod$pairs.mod,
                        target = target,
                        delta = 600,
                        r = 2000)
      rep<-c(rep, temp[[1]])
    }
    results<-rbind(results, rep)
  }
  save(ld.points, results, file=paste("out/nn", nn[kk], ".RData", sep=""))
}

summ<-c()
for(ii in 1:length(nn)){
  load(file=paste("nn", nn[ii], ".RData", sep=""))
  summ<-rbind(summ, apply(results, 2, function(x){sum(x!=200)/100}))
}

png(file="../Figures/n_remove.png", height=12, width=12, res=300, units="cm")
plot(nn, summ[,2], pch=19, xlab="Number of artificial waterbodies removed", ylab="Probability of stopping toads", bty="l")
lines(nn, summ[,2])
points(nn, summ[,1], pch=19, col="green")
lines(nn, summ[,1], col="green")
legend("topright", legend=c("Barrier 1", "Barrier 2"), pch=19, col=c("black", "green"), bty="n")
dev.off()
