library(dplyr)

######## load functions ########
source("src/pprocess_functions.R")

######## load data ########
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
      threshold = 50,
      constant.rain = 180,
      trunc.dist = NULL)


######## how long to the pilbara (in wet seasons from dry season of 2024) ########
reps <- 100
output<-vector("list", length=reps) # vector to take outputs

for (rr in 1:reps){ # for reps
  cat("Rep ", rr, "\n")
  lambda.samp<-10^rnorm(1, mean=sample.lambda, sd=sample.lambda.sd)
  r.samp<-10^2
  temp<-spread.pilb(pop=spread.table, gens=100, pairs=pairs, delta=lambda.samp, r=r.samp)
  temp<-c(temp, list(pars=cbind(lambda=lambda.samp, r=r.samp)))
  output[[rr]]<-temp
}

######## save output and generate summaries ########
# time to arrive in Pilbara
save(output, file = "out/timing-to-pilbara.Rdata")

time.vec <- unlist(lapply(output, FUN = function(x){c(x$gen)}))

pdf(file = "out/time-to-pilbara.pdf")
  hist(time.vec, xlab = "Time to the Pilbara (y)")
dev.off()

summary(time.vec)

# Work out mean time to and arrival year for each point

add_index_column <- function(x, index) { # Function to add index column to each matrix
  mat <- x$popmatrix[x$popmatrix[, "age"] > 0, ] # remove points never colonised
  index_col <- rep(index, n = nrow(mat))
  arrival <- round(2024+max(mat[, "age"])-mat[, "age"]) # calculate arrival year
  cbind(index_col, mat, arrival)
}
# Apply the function to each element of the list using lapply
modified_matrices <- lapply(seq_along(output), function(i) {
  add_index_column(output[[i]], i)
})
# bind the lot together into single matrix
pop.out <- do.call("rbind", modified_matrices)

# get mean arrival time for each point, or making a static map
pop.summary <- pop.out %>% 
  as.data.frame() %>%
  group_by(ID) %>%
  summarise_all(mean)
write.csv(pop.summary, file = "out/spread-TCZ-timing.csv", row.names = FALSE) 

# to make a dynamic map...
# for each year, 2024 to max(mean arrival time), generate a csv to plot, that reports probability of colonisation at that time for each waterpoint
max.time <- max(pop.summary$arrival) + 2
for (yy in 2024:max.time){
  fname <- paste0("out/dynamic_maps/", yy, ".csv")
  pop.summary <- pop.out %>%
    as.data.frame() %>%
    group_by(ID) %>%
    summarise(X = mean(X), Y = mean(Y), prob.colonised = mean(arrival <= yy)) %>%
    filter(prob.colonised > 0.1)
  write.csv(pop.summary, file = fname, row.names = FALSE) 
}

# make figures
source("src/figures/output-maps.R")
