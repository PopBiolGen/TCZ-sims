######## load functions and libraries ########
source("src/pprocess_functions.R")

######## note ########
# large number of points = big memory demand.  Useful to set
# usethis::edit_r_environ()
# R_MAX_VSIZE=100Gb

######## load data ########
load("dat/Posteriors_2023.RData")

k.points <- read.csv("dat/waterpoint-data_Kimberley_trimmed.csv")
k.points$origin_des <- as.numeric(as.factor(k.points$origin_des))-1
k.points$colonised[is.na(k.points$colonised)] <- 0
k.points$colonised[k.points$colonised==2009] <- 1

setup(point.data = k.points, 
      X.id = "X", 
      Y.id = "Y",
      present.id = "colonised",
      artificial.natural.id = "origin_des",
      rain.id = "ndays_1",
      threshold = 100,
      constant.rain = NULL,
      trunc.dist = TRUE,
      TCZ = FALSE)

rm(k.points) #free up some memory
gc()

# make a target
spread.table[spread.table[,"X"] == min(spread.table[,"X"]), "target"] <- 1



system.time( temp <- run_sims(n.sims = 1, gens = 14, plot = TRUE, rollup = TRUE))
