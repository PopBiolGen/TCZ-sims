######## load functions and libraries ########
source("src/pprocess_functions.R")

######## load data ########
load("dat/Posteriors.RData")

tczPoints <-read.csv("dat/merged_clipped_for_simulation.csv")
tczPoints$origin_des <- as.numeric(as.factor(tczPoints$origin_des))-1


######## worst case scenario ########
scen.name <- "worst-case"
setup(point.data = tczPoints, 
      X.id = "X", 
      Y.id = "Y",
      present.id = "colonised",
      artificial.natural.id = "origin_des",
      rain.id = "daysRain_1",
      threshold = 50,
      constant.rain = 180,
      trunc.dist = NULL)

run_sims(scen.name)
make_plots(scen.name)

######## likely case scenario ########
scen.name <- "likely-case"
setup(point.data = tczPoints, 
      X.id = "X", 
      Y.id = "Y",
      present.id = "colonised",
      artificial.natural.id = "origin_des",
      rain.id = "daysRain_1",
      threshold = 50,
      constant.rain = NULL,
      trunc.dist = 80000)

run_sims(scen.name)
make_plots(scen.name)