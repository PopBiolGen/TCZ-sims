# Script to run sims over the TCZ
# to test various control scenarios for efficacy and optimise investment in control interventions

######## load functions and libraries ########
source("src/pprocess_functions.R")

######## load data ########
load("dat/Posteriors.RData")

# define the data directory
data.dir <- file.path(Sys.getenv("DATA_PATH"), "Toads/TCZ/infrastructure")
data.dump.id <- "483d4936-3b85-41c1-bf75-eb0f0fda3c27"
# load TCZ data 
tcz.sites <- st_read(file.path(data.dir, data.dump.id, "water_point_audit.geojson"))
tcz.infrastructure <- st_read(file.path(data.dir, data.dump.id, "water_point_audit_infrastructure_item.geojson"))
