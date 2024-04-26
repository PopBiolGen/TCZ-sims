setwd("~/Dropbox/Papers/Submitted/Artificial waterbodies/Data")

load(file="Greg's convolution resamples.RData")
load(file="Kernel_fits.RData")

source("convolution functions.R")

################################################
# representative fits #
xmax<-10000
samps<-c(1, 80, 160)
png(file="../Figures/Kernel fits.png", res=300, height=36, width=12, units="cm")
par(mfrow=c(3,1), mar=c(5,5,4,2), cex.lab=1.5)
labels=c("A)", "B)", "C)")
for(ii in 1:length(samps)){
  ttl<-paste(samps[ii], ifelse(ii==1, "day of movement", "days of movement"))
  hist(resamps[,samps[ii]], breaks=80, xlim=c(0, xmax), freq=F, col="grey80", border="grey80", ylab="Probability density", 
       xlab="Distance (m)", main=ttl, ylim=c(0, max(collapse(1:xmax, fits[samps[ii],1], fits[samps[ii],2]))))
  lines(1:xmax, collapse(1:xmax, fits[samps[ii],1], fits[samps[ii],2]))
  mtext(side=3, line=1, adj=0, text=labels[ii])
  if(ii==1) legend("center", legend=c("Resampled data", "Fitted kernel"), fill=c("grey80", NA), border=NA, lty=c(NA, 1), 
    bty="n", cex=1.3, merge=T)
}
dev.off()
################################################

################################################
# Trends in u and v over n
png(file="../Figures/Kernel parameters over n.png", res=300, height=16, width=12, units="cm")
par(mar=c(5,4,4,5))
plot(fits[,"u"], xlab="Number of convolutions", ylab=expression("Fitted scale ("*italic("u")*")"), type="l")
par(new=T)
plot(fits[,"v"], axes="F", xlab="", ylab="", type="l", col="green")
axis(side=4)
mtext(text=expression("Fitted shape ("*italic("v")*")"), side=4, line=3)
legend("bottomright", legend=c(expression(italic("u")), expression(italic("v"))), lty=1, col=c("black", "green"), bty="n")
dev.off()
################################################