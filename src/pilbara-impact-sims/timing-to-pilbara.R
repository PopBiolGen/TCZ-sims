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
      trunc.dist = TRUE,
      TCZ = FALSE)

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
      trunc.dist = TRUE,
      TCZ = FALSE)

run_sims(scen.name)
make_plots(scen.name)

######## worst case scenario, but with TCZ ########
scen.name <- "worst-case_TCZ"
setup(point.data = tczPoints, 
      X.id = "X", 
      Y.id = "Y",
      present.id = "colonised",
      artificial.natural.id = "origin_des",
      rain.id = "daysRain_1",
      threshold = 50,
      constant.rain = 180,
      trunc.dist = TRUE,
      TCZ = TRUE)

run_sims(scen.name)
make_plots(scen.name)

######## likely scenario, but with TCZ ########
scen.name <- "likely-case_TCZ"
setup(point.data = tczPoints, 
      X.id = "X", 
      Y.id = "Y",
      present.id = "colonised",
      artificial.natural.id = "origin_des",
      rain.id = "daysRain_1",
      threshold = 50,
      constant.rain = NULL,
      trunc.dist = TRUE,
      TCZ = TRUE)

run_sims(scen.name)
make_plots(scen.name)