######## load functions and libraries ########
source("src/pprocess_functions.R")


######## load data ########
load("dat/Posteriors.RData")
scen.name <- "paruku-arrival"

pk.points <- read.csv("dat/waterpoint-data_Paruku.csv")
# fix a few idiosyncratic things
pk.points$origin_des <- as.numeric(as.factor(pk.points$origin_des))-1
pk.points$colonised[is.na(pk.points$colonised)] <- 0


setup(point.data = pk.points,
      X.id = "X",
      Y.id = "Y",
      present.id = "colonised",
      artificial.natural.id = "origin_des",
      rain.id = "ndays_1",
      threshold = 100,
      constant.rain = NULL,
      trunc.dist = TRUE,
      TCZ = FALSE)

save(spread.table, pairs_pdist, file = "out/setup_paruku_complete.RData")


# run a simulation to make a map
pk.spread <- run_sims(n.sims = 100, gens = 50, plot = FALSE, rollup = FALSE)
save_outputs(output = pk.spread, path = "out", scenario.name = scen.name, start.year=2011)
make_plots(scen.name)
