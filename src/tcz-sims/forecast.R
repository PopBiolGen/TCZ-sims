# Script to run sims from the current invasion front, until the tcz is reached

######## load functions and libraries ########
source("src/tcz-sims/load-data.R")

######## define target points ########
all.points <- all.points |>
  mutate(colonised = ifelse(st_within(geometry, tcz.boundary, sparse = FALSE)[, 1], 2L, colonised))

######## load simulation posteriors ########
load("dat/Posteriors_2026.RData")

start_year <- 2026

######## tcz-arrival scenario ########
scen.name <- "forecast"
setup(point.data = all.points,
      X.id = "X",
      Y.id = "Y",
      present.id = "colonised",
      artificial.natural.id = "origin_des",
      rain.id = "rainfall",
      threshold = 100,
      constant.rain = NULL,
      trunc.dist = TRUE,
      TCZ = FALSE,
      col.year.id = "colonisation_year",
      start.year = start_year)
# write out points for basemapping
write.csv(spread.table, file = "out/basemap_points.csv", row.names = FALSE)
# run sims..
sim_out <- run_sims(n.sims = 100, gens = 20, plot = FALSE, rollup = FALSE)
save_outputs(output = sim_out,
             path = "out",
             scenario.name = scen.name,
             start.year=start_year,
             plot.time = TRUE,
             extinction.aware = TRUE)
make_plots(scen.name, plot.year = TRUE, tcz.boundary = tcz.boundary)

## Output data to shiny app
source("shiny/prep_app_data.R")
