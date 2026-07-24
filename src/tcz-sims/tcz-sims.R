# Script to run sims over the TCZ
# to test various control scenarios for efficacy and optimise investment in control interventions

######## load functions and libraries ########
source("src/tcz-sims/load-data.R")

######## filter to TCZ region ########
# buffer TCZ by 80km (project to Albers for metric distances, then back)
tcz_buffer_80 <- tcz.boundary |>
  st_transform(3577) |>
  st_buffer(80000) |>
  st_union() |>
  st_transform(st_crs(all.points))

all.points <- all.points[st_within(all.points, tcz_buffer_80, sparse = FALSE)[, 1], ]

######## define colonised and target points ########
# Use TCZ bounding box to define north/south: 
#   north of TCZ → already colonised (colonised = 1)
#   south of TCZ → target points    (colonised = 2)
tcz_bbox   <- st_bbox(st_transform(tcz.boundary, st_crs(all.points)))
pt_y       <- st_coordinates(all.points)[, "Y"]
within_tcz <- st_within(all.points, tcz.boundary, sparse = FALSE)[, 1]

all.points <- all.points |>
  mutate(
    colonised = case_when(
      within_tcz                 ~ 0L,  # inside TCZ → uncolonised simulation space
      pt_y > tcz_bbox[["ymax"]] ~ 1L,  # north of TCZ → already colonised
      pt_y < tcz_bbox[["ymin"]] ~ 2L,  # south of TCZ → target points
      TRUE ~ colonised
    ),
    origin_des = (origin_des == "Natural") # to match logical required of setup()
  ) |> 
  rename(TCZ = inside_tcz) # to match required name in setup()

# Plot to make sure all is as it should be
bbox <- st_bbox(tcz_buffer_80)

ggplot() +
  geom_sf(data = wa.coast, fill = NA, color = "grey30") +
  geom_sf(data = tcz_buffer_80, fill = NA, color = "blue", lty = 2) +
  geom_sf(data = tcz.boundary, fill = NA, color = "red") +
  geom_sf(data = all.points,
          aes(color = factor(colonised,
                             levels = c(0, 1, 2),
                             labels = c("Uncolonised (TCZ)", "Colonised (N)", "Target (S)")),
              shape = origin_des),
          size = 0.8) +
  scale_color_manual(values = c("Uncolonised (TCZ)" = "grey50",
                                "Colonised (N)"     = "firebrick",
                                "Target (S)"        = "steelblue"),
                     name = "Status") +
  scale_shape_manual(values = c(`TRUE` = 8, `FALSE` = 16),
                     labels = c("Natural", "Manmade"),
                     name = "Origin") +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
           ylim = c(bbox["ymin"], bbox["ymax"]))


######## load simulation posteriors ########
load("dat/Posteriors_2026.RData")


######## tcznfull-imlementation scenario ########
scen.name <- "tcz-full"
setup(point.data = all.points, 
      X.id = "X", 
      Y.id = "Y",
      present.id = "colonised",
      artificial.natural.id = "origin_des",
      rain.id = "rainfall",
      threshold = 100,
      constant.rain = NULL,
      trunc.dist = TRUE,
      TCZ = TRUE)

# write out points for basemapping
write.csv(spread.table, file = "out/basemap_points.csv", row.names = FALSE)
# run sims..
sim_out <- run_sims(n.sims = 100, gens = 20, plot = FALSE, rollup = FALSE)
save_outputs(output = sim_out, 
             path = "out", 
             scenario.name = scen.name, 
             start.year=2032,
             plot.time = TRUE)
make_plots(scen.name, plot.year = TRUE, tcz.boundary = tcz.boundary)
