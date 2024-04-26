setwd("~/final_natural_water_files/ABC_output")

#load(file="ABC Output.RData")
concat<-c()
for (ii in 1:100){
  load(paste(ii, ".rdata", sep=""))
  concat<-rbind(concat, out)
  rm(out)
}

concat<-cbind(concat, log.r=log(concat[, "r"], 10))
concat<-cbind(concat, log.lambda=log(concat[, "delta"], 10))
lm.set<-concat[concat[,"temp"]>=29, ]
mod.r<-lm(lm.set[,"log.r"]~lm.set[,"temp"])
mod.lambda<-lm(lm.set[,"log.lambda"]~lm.set[,"temp"])
best<-concat[concat[,"temp"]>=35, ]
zero.r<-sum(mod.r$coefficients*c(1, 36))
zero.lambda<-sum(mod.lambda$coefficients*c(1, 36))
best.diff<-36-best[,"temp"]
adj.r<-mod.r$coefficients[2]*best.diff
adj.lambda<-mod.lambda$coefficients[2]*best.diff
post.r<-best[,"log.r"]+adj.r
post.lambda<-best[,"log.lambda"]+adj.lambda

library(emdbook)
cvr<-cov(post.r, post.lambda)
M<-matrix(c(var(post.r),cvr,cvr,var(post.lambda)),nrow=2)
X<-seq(min(post.r), max(post.r), length.out=100)
Y<-seq(min(post.lambda), max(post.lambda), length.out=100)
XY<-as.matrix(expand.grid(x=X, y=Y))


png(file="Posterior bivariate.png", width=16, height=16, res=300, units="cm")
par(mar=c(5,5,4,2))
#fig.mat<-matrix(c(1,1,4,3,3,2,3,3,2), nrow=3)
#layout(fig.mat)
#dens<-density(post.lambda)
#plot(dens$y, dens$x, bty="l", ylab=expression(log[10](italic(hat(C)))), xlab="Density", main="", type="l", ylim=c(3.5,6)) #"log10 mean.lambda"
#dens<-density(post.r)
#plot(dens$x, dens$y, bty="l", ylab="Density", xlab="log10 r", main="", type="l", xlim=c(1.5,3))
plot(post.r, post.lambda, xlab=expression(log[10](italic(r))), ylab=expression(log[10](italic(hat(C)))),
	bty="l", xlim=c(1.5,3), ylim=c(3.5,6))
#dens<-dmvnorm(XY, mu=c(mean(post.r), mean(post.lambda)), Sigma=M)
#contour(X, Y, matrix(dens, nrow=100), add=TRUE)
mod1<-lm(post.lambda~post.r)
abline(reg=mod1)
dev.off()


sample.lambda<-predict(mod1, data.frame(post.r=2))
sample.lambda.sd<-sd(mod1$residuals)
mean.r<-summary(post.r)
mean.lambda<-summary(post.lambda)


save(concat, best, post.r, post.lambda, file="ABC Output.RData")
save(M, mean.r, mean.lambda, sample.lambda, sample.lambda.sd, file="../Posteriors.RData")