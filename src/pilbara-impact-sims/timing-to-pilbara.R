######## load functions ########
source("src/pprocess_functions.R")

######## load data ########
#load("dat/Kernel_fits.RData")
load("dat/Posteriors.RData")

tczPoints <-read.csv("dat/merged_clipped_for_simulation.csv")
tczPoints$origin_des <- as.numeric(as.factor(tczPoints$origin_des))-1


######## setup the environment ########
setup(point.data = tczPoints, 
      X.id = "X", 
      Y.id = "Y",
      present.id = "colonised",
      artificial.natural.id = "origin_des",
      rain.id = "daysRain_1",
      maximum=100000)


######## how long to the pilbara (in wet seasons from 2023/4) ########
reps <- 100
output<-vector("list", length=reps) # vector to take outputs

for (rr in 1:reps){ # for reps
  lambda.samp<-10^rnorm(1, mean=sample.lambda, sd=sample.lambda.sd)
  r.samp<-10^2
  temp<-spread.pilb(pop=spread.table, gens=100, pairs=pairs, target=target, delta=lambda.samp, r=r.samp)
  temp<-c(temp, list(pars=cbind(lambda=lambda.samp, r=r.samp)))
  output[[rr]]<-temp
}

save(output, file = "out/timing-to-pilbara.Rdata")

time.vec <- unlist(lapply(output, FUN = function(x){c(x$gen)}))

pdf(file = "out/time-to-pilbara.pdf")
  hist(time.vec, xlab = "Time to the Pilbara (y)")
dev.off()

summary(time.vec)