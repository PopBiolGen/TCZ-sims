######## load functions and libraries ########
source("src/pprocess_functions.R")

######## load data ########
load("dat/Posteriors_2023.RData")

tczPoints <-read.csv("dat/waterpoint-data_LaGrange.csv")
tczPoints$origin_des <- as.numeric(as.factor(tczPoints$origin_des))-1


######## do-nothing scenario ########
scen.name <- "do-nothing"
setup(point.data = tczPoints, 
      X.id = "X", 
      Y.id = "Y",
      present.id = "colonised",
      artificial.natural.id = "origin_des",
      rain.id = "ndays_1",
      threshold = 100,
      constant.rain = NULL,
      trunc.dist = TRUE,
      TCZ = FALSE)
# write out points for basemapping
write.csv(spread.table, file = "out/basemap_points.csv", row.names = FALSE)
# run sims..
sim_out <- run_sims(gens = 100, plot = FALSE, rollup = FALSE)
save_outputs(output = sim_out, path = "out", scenario.name = scen.name, start.year=2023)
make_plots(scen.name)

######## TCZ scenario ########
scen.name <- "TCZ"
setup(point.data = tczPoints, 
      X.id = "X", 
      Y.id = "Y",
      present.id = "colonised",
      artificial.natural.id = "origin_des",
      rain.id = "ndays_1",
      threshold = 100,
      constant.rain = NULL,
      trunc.dist = TRUE,
      TCZ = TRUE)

sim_out <- run_sims(gens = 50, plot = FALSE, rollup = FALSE)
save_outputs(output = sim_out, path = "out", scenario.name = scen.name, start.year=2023)
make_plots(scen.name)