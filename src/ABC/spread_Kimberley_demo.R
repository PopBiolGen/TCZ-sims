######## load functions and libraries ########
source("src/pprocess_functions.R")

######## note ########
# large number of points = big memory demand.  Useful to set
# usethis::edit_r_environ()
# R_MAX_VSIZE=100Gb

######## load data ########
load("dat/Posteriors_2023.RData")
scen.name <- "kimberley-demo"

k.points <- read.csv("dat/waterpoint-data_Kimberley.csv")
# fix a few idiosyncratic things
k.points$origin_des <- as.numeric(as.factor(k.points$origin_des))-1
k.points$colonised[is.na(k.points$colonised)] <- 0
k.points$colonised[k.points$colonised==2009] <- 1
# remove ABC data
k.point <- k.points[, !(names(k.points) %in% c("X2021_pres", "X2022_pres", "X2023_pres"))]

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
# run a simulation to make a map
k.spread <- run_sims(n.sims = 1, gens = 18, plot = FALSE, rollup = TRUE)
save_outputs(output = k.spread, path = "out", scenario.name = scen.name, start.year=2009)
make_plots(scen.name)
