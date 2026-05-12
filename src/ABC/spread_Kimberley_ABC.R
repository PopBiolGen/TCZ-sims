######## load functions and libraries ########
source("src/pprocess_functions.R")

######## note ########
# large number of points = big memory demand.  Useful to set
# usethis::edit_r_environ()
# R_MAX_VSIZE=100Gb

######## load data ########
load("dat/Posteriors.RData")

# k.points <- read.csv("dat/waterpoint-data_Kimberley_trimmed.csv")
# k.points$origin_des <- as.numeric(as.factor(k.points$origin_des))-1
# k.points$colonised[is.na(k.points$colonised)] <- 0
# k.points$colonised[k.points$colonised==2009] <- 1
# 
# setup(point.data = k.points,
#       X.id = "X",
#       Y.id = "Y",
#       present.id = "colonised",
#       artificial.natural.id = "origin_des",
#       rain.id = "ndays_1",
#       observations = c("X2021_pres", "X2022_pres", "X2023_pres"),
#       threshold = 100,
#       constant.rain = NULL,
#       trunc.dist = TRUE,
#       TCZ = FALSE)
# save(spread.table, pairs, file = "out/setup_kimberley.RData")
# rm(k.points) #free up some memory
# gc()

load(file = "out/setup_kimberley.RData")

# Place some more variance on the old priors...
sample.lambda.sd <- 2.5 * sample.lambda.sd


# run
ii <- 1
  sim.set <- run_sims(n.sims = 2000, gens = 14, plot = FALSE, rollup = TRUE)
  save_outputs(sim.set, path = "out/ABC", scenario.name = paste0("ABC_", ii), start.year = 2009, ABC = TRUE)

