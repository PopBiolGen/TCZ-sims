######## load functions and libraries ########
source("src/pprocess_functions.R")

######## load data ########
load("dat/Posteriors.RData")

tczPoints <-read.csv("dat/waterpoint-data_LaGrange.csv")
tczPoints$origin_des <- as.numeric(as.factor(tczPoints$origin_des))-1
tczPoints <- select(tczPoints, X, Y, origin_des, ndays_1, colonised, TCZ)
tczPoints$colonised[is.na(tczPoints$colonised)] <- 0

# read in additional points from Tim
d_extra <- st_read("dat/tims_points.kml")
d_extra <- st_transform(d_extra, crs = 3577) # convert to Albers
coords <- st_coordinates(d_extra) %>% 
  as.data.frame() %>%
  mutate(origin_des = 0, ndays_1 = 34, colonised = 0, TCZ = 0) %>%
  select(-Z)

tczPoints <- rbind(tczPoints, coords) 
            

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
sim_out <- run_sims(n.sims = 50, gens = 50, plot = FALSE, rollup = FALSE)
save_outputs(output = sim_out, path = "out", scenario.name = scen.name, start.year=2023)
make_plots(scen.name, plot.year = FALSE)

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

sim_out <- run_sims(n.sims = 50, gens = 50, plot = FALSE, rollup = FALSE)
save_outputs(output = sim_out, path = "out", scenario.name = scen.name, start.year=2023)
make_plots(scen.name)