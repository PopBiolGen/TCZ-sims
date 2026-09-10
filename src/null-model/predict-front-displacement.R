# Null model for invasion speed: for each year with an empirical front-position
# estimate, start the point-process spread model from that front and run it
# forward one generation (~one wet season = one calendar year -- see rain_to_days()
# in pprocess_functions.R), many times, then compare the simulated front
# displacement against the empirical year-over-year displacement estimated in the
# invasion-front-monitoring repo (time-to-detection-multi-year_JAGS.R).
#
# Requires invasion-front-parameters-multi-year.Rdata (mod.multi, post.samples,
# bbox.fit, mean.coord, scale, years.all) to already exist in $DATA_PATH --
# produced by (re-)running that script in invasion-front-monitoring.
#
# Notes:
# - The starting front for year k is fixed at its posterior mean (a, b_k) -- not
#   resampled per rep -- so the resulting spread of b_sim across reps reflects only
#   the stochastic spread process (propagule draws, colonisation draws, lambda
#   resampling in run_sims()), not front-location uncertainty. That's compared
#   against the empirical delta[k] posterior, which instead reflects detection-data
#   uncertainty -- these are two different kinds of uncertainty.
# - b_sim is recovered per rep via a deterministic misclassification-minimising
#   threshold on the simulated Pres pattern (fit_b_threshold()), not a smooth
#   logistic fit -- this matches the deterministic step-function front definition
#   used by score_colonised()/score_colonised_multi_year(), and avoids
#   (quasi-)separation problems a logistic fit would hit when only a narrow ring of
#   points near the front changes status in a single generation.
# - The threshold fit is restricted to points inside bbox.fit (the Albers-metres
#   bounding box of the data the empirical front line was actually fit to), even
#   though the simulation itself runs on the full point network -- the front is
#   only approximately straight over that shorter span.

######## load functions, libraries, and point data ########
source("src/tcz-sims/load-data.R") # produces all.points (colonisation_year, Albers X/Y already attached)

######## load empirical front parameters (produced by invasion-front-monitoring) ########
# Moved up from its original position (after all.points/survey.point.bbox construction) --
# survey.point.bbox below needs a/mean.coord_m/b_metres() to build a front-relative cutoff.
load(file.path(Sys.getenv("DATA_PATH"), "invasion-front-parameters-multi-year.Rdata"))
# objects brought into scope: mod.multi, post.samples, bbox.fit, mean.coord, scale, years.all, run.info

a <- mod.multi[[1]]["a", 1]        # posterior-mean slope (shared across years, fixed across reps)
mean.coord_m <- mean.coord * scale # mean.coord (km, model space) -> Albers metres

# posterior-mean b (Albers metres) for years.all[tt] -- same transform score_colonised_multi_year() uses
b_metres <- function(tt) {
  b_raw <- mod.multi[[1]][paste0("b[", tt, "]"), 1]
  (b_raw + mean.coord["Y"]) * scale
}
b1 <- b_metres(1) # earliest modelled year's intercept, used to build survey.point.bbox below

# replace old.lagrange.points with points from a different geographic area, to cover the area
# covered by toad monitoring surveys plus an 80km buffer -- except on the already-colonised
# (east) side, where the boundary follows a line parallel to the invasion front's axis, 80km
# (perpendicular) behind the earliest modelled front (years.all[1]). A plain axis-aligned bbox
# buffer can't guarantee uniform perpendicular clearance from a diagonal front (slope `a`); this
# does. West/north/south sides aren't reported as problematic and keep the original 80km buffer.
# Note, invasion front analysis needs to be run first, and its repo needs to be sitting in same
# parent directory as this repo.
load(file = "../invasion-front-monitoring/out/merged-visual-surveys.RData")
survey.bb <- df |>
  st_transform(crs = 3577) |>
  st_bbox()

# west/north/south: original 80km axis-aligned buffer. East side left generously unconstrained
# (1000km) -- irrelevant once intersected with front.cutoff.poly below (that's always the
# binding edge in practice), so its exact value doesn't need tuning.
survey.box <- st_bbox(c(xmin = survey.bb[["xmin"]] - 80000,
                         ymin = survey.bb[["ymin"]] - 80000,
                         xmax = survey.bb[["xmax"]] + 1e6,
                         ymax = survey.bb[["ymax"]] + 80000),
                       crs = 3577) |>
  st_as_sfc()

# front-relative cutoff: half-plane D <= b1.offset, D = Y - a*(X - meanX_m) (score_colonised()'s
# convention: colonised iff D > b). b1.offset shifts b1 80km into colonised (increasing-D)
# territory: a perpendicular offset `dist` in D-space is `dist * sqrt(1+a^2)` (same conversion
# used later in this script for disp.sim.km).
b1.offset <- b1 + 80000 * sqrt(1 + a^2)

# Offset line D = b1.offset, built generously far along the front's axis (u = X - meanX_m).
# singleSide buffers give the LEFT side of the line's direction of travel for positive dist,
# RIGHT for negative (?st_buffer). For a line traversed in the direction of increasing u (tangent
# (1,a)), the left-hand normal is (-a,1) == grad(D) (D's increasing direction) -- so the *right*
# side (negative dist) is the decreasing-D/uncolonised-ish side we want, for any sign of `a`.
u.vals <- c(-2e6, 2e6) # +/-2000km along the front's axis -- generously covers WA
front.offset.line <- st_linestring(cbind(mean.coord_m[["X"]] + u.vals,
                                          a * u.vals + b1.offset))
front.cutoff.poly <- st_sfc(front.offset.line, crs = 3577) |>
  st_buffer(dist = -2e6, singleSide = TRUE, endCapStyle = "FLAT")

survey.point.bbox <- st_intersection(survey.box, front.cutoff.poly)
stopifnot(length(survey.point.bbox) > 0, !st_is_empty(survey.point.bbox))
# essentially all df survey points should fall inside the final box -- if not, something is
# badly wrong with the cutoff direction/parameters, not just a few points on a tight edge
n.survey.outside <- sum(!st_within(st_transform(df, crs = 3577), survey.point.bbox, sparse = FALSE))
if (n.survey.outside > 0) {
  warning(n.survey.outside, " of ", nrow(df), " survey points fall outside survey.point.bbox")
}

# replace all.points with the relevant set of waterpoints from combined waterpoint data
all.points <- st_read(file.path(spatial.dir, "Edited-layers/merged-points.shp")) |> 
  select(fcsubtype_, full_name, perennia_1, origin_des, watercou_1, area_m, st_perimet) |> 
  st_transform(crs = 4326) |>
  bind_rows(d_extra) |> # add additional points from aerial imagery
  bind_rows(df |> select(full_name = display_label, geometry)) |>  # add additional points from on-ground surveys
  st_transform(crs = 3577) |>
  st_filter(survey.point.bbox, .predicate = st_within)
  
rainfall.vals <- terra::extract(rainfall.raster, terra::vect(all.points))
all.points <- all.points |>
  mutate(rainfall = rainfall.vals[[2]])

####### Bring in estimated invasion front(s) and scored colonised points #######
colnsd <- score_colonised_multi_year(all.points)
all.points <- colnsd$scored.points # includes colonised (0/1) and colonisation_year
inv.front <- colnsd$fronts # one front line per modelled year
fp <- colnsd$fp; rm(colnsd)

# get coordinates in albers
albers <- all.points |> 
  st_transform(crs = 3577) |> 
  st_coordinates()
all.points <- bind_cols(all.points, albers)


######## load simulation posteriors ########
load("dat/Posteriors_2026.RData")

## Make a static figure here to show selected points, and the invasion front lines on a map.
dir.create("out/null-model", showWarnings = FALSE)
bbox <- st_bbox(all.points)
selected.points.fig <- ggplot() +
  geom_sf(data = wa.coast, fill = NA, color = "grey30") +
  geom_sf(data = tcz.boundary, fill = NA, color = "grey60", linetype = "dashed") +
  geom_sf(data = inv.front, aes(group = factor(year)), color = "red", lty = 2) +
  geom_sf(data = all.points, aes(color = factor(colonisation_year)), size = 0.7) +
  geom_sf(data = fp, color = "red") +
  labs(color = "Colonisation year") +
  # all.points/inv.front/fp are Albers (3577) here, unlike wa.coast/tcz.boundary (4326) -- pin the
  # plot CRS explicitly so the Albers-metre bbox isn't misread as lon/lat by coord_sf's default (first-layer) CRS
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
           ylim = c(bbox["ymin"], bbox["ymax"]),
           crs = st_crs(all.points))
ggsave(filename = "out/null-model/selected-points-and-fronts.pdf", plot = selected.points.fig,
       width = 180, height = 160, units = "mm")


######## helper: recover b from a colonisation pattern ########
# Given perpendicular-axis position D = Y - a*(X - meanX_m) and colonisation status Pres
# (0/1) for a set of points, find the b that best separates colonised from uncolonised
# points under the rule "colonised if D > b" (score_colonised()'s own rule) -- i.e. the
# threshold minimising misclassifications against the observed Pres pattern. Ties broken
# by taking the middle of the tied set. Returns NA if every point falls on one side (no
# separating threshold exists within the data).
fit_b_threshold <- function(D, Pres) {
  ord <- order(D)
  D <- D[ord]; Pres <- Pres[ord]
  n <- length(D)
  if (n < 2) return(NA_real_)
  total1 <- sum(Pres)
  cum1 <- c(0, cumsum(Pres))          # cum1[i+1] = colonised count among first i points, i = 0..n
  i.seq <- 0:n
  errors <- 2 * cum1 - i.seq - total1 + n  # errors[i+1] = misclassifications if 1..i called uncolonised
  best <- which(errors == min(errors)) - 1L # back to i-indexing (0..n)
  i.star <- best[ceiling(length(best) / 2)]
  if (i.star <= 0L || i.star >= n) return(NA_real_) # degenerate: whole window on one side
  mean(D[i.star:(i.star + 1L)])
}

######## quick self-test of fit_b_threshold() on synthetic data ########
# (independent of real data/geometry -- just checks the fitting logic itself)
{
  set.seed(1)
  D.test <- sort(runif(500, 0, 1000))
  b.true <- 550
  Pres.test <- as.integer(D.test > b.true)
  b.rec <- fit_b_threshold(D.test, Pres.test)
  stopifnot(abs(b.rec - b.true) < diff(range(D.test)) / length(D.test) * 2)
}

######## build spread table + dispersal matrix once (year-invariant part of setup()) ########
# Everything setup() computes other than Pres/age/occ.years/first.col.gen -- point geometry,
# the deduplicated point set, per-point kernel params (u/v/max.dist/area, from rainfall only),
# and above all the sparse pairwise dispersal matrix `pairs` (the expensive O(n^2)-ish step) --
# depends only on point location and rainfall, never on colonisation state. So it only needs
# building once; each year below just overwrites the colonisation-dependent columns of the
# global spread.table in place before calling run_sims(), instead of rebuilding from scratch.
# present_k is a placeholder here (overwritten every iteration below); TCZ/control-step args are
# left at their defaults (no targets, no control step) so target/control columns stay constant.
all.points$present_k <- 0L
setup(point.data = all.points,
      X.id = "X",
      Y.id = "Y",
      present.id = "present_k",
      artificial.natural.id = "origin_des",
      rain.id = "rainfall",
      threshold = 100,
      constant.rain = NULL,
      trunc.dist = TRUE,
      TCZ = FALSE)

# original all.points row index for each surviving (post-dedup) spread.table row, so a given
# year's present_k (indexed by all.points row) can be aligned to spread.table's row order
kept.ids <- spread.table[, "ID"]

######## years to simulate forward from (every year except the last has a "next year") ########
sim.years <- years.all[-length(years.all)]
n.sims <- 100 # sim reps per year -- raise for a smoother null distribution (runtime scales linearly)

null.results <- vector("list", length(sim.years))

for (tt in seq_along(sim.years)) {
  yr <- sim.years[tt]
  cat("\n==== Year", yr, "====\n")

  ######## starting colonisation state for this year, from the multi-year front ########
  # written directly into spread.table's Pres/age/occ.years/first.col.gen columns -- this is
  # exactly what setup() would compute for these columns given col.year.id = NULL (see
  # pprocess_functions.R), just without rebuilding the dispersal matrix each time
  present_k_full <- ifelse(!is.na(all.points$colonisation_year) &
                              all.points$colonisation_year <= yr, 1L, 0L)
  present_k <- present_k_full[kept.ids]
  spread.table[, "Pres"] <- present_k
  spread.table[, "age"] <- present_k
  spread.table[, "occ.years"] <- present_k
  spread.table[, "first.col.gen"] <- ifelse(present_k == 1, 0, NA_real_)

  ######## sanity check: threshold-fit should recover ~b_k from the (deterministic) starting state ########
  b_k <- b_metres(tt)
  in.bbox0 <- spread.table[, "X"] >= bbox.fit["xmin"] & spread.table[, "X"] <= bbox.fit["xmax"] &
    spread.table[, "Y"] >= bbox.fit["ymin"] & spread.table[, "Y"] <= bbox.fit["ymax"]
  D0 <- spread.table[, "Y"] - a * (spread.table[, "X"] - mean.coord_m["X"])
  b0.check <- fit_b_threshold(D0[in.bbox0], spread.table[in.bbox0, "Pres"])
  cat("  Sanity check: b recovered from starting state =", round(b0.check),
      "m vs fitted b_k =", round(b_k), "m (diff =", round(b0.check - b_k), "m)\n")

  ######## run one generation (= one year) forward, many reps ########
  sim_out <- run_sims(n.sims = n.sims, gens = 1, plot = FALSE, rollup = FALSE,
                       stop.on.target = FALSE)

  ######## fit b_sim per rep, restricted to the empirical fitting bbox ########
  disp.sim.km <- vapply(sim_out, function(rep) {
    pm <- rep$popmatrix
    in.bbox <- pm[, "X"] >= bbox.fit["xmin"] & pm[, "X"] <= bbox.fit["xmax"] &
      pm[, "Y"] >= bbox.fit["ymin"] & pm[, "Y"] <= bbox.fit["ymax"]
    pm <- pm[in.bbox, , drop = FALSE]
    D <- pm[, "Y"] - a * (pm[, "X"] - mean.coord_m["X"])
    b_sim <- fit_b_threshold(D, pm[, "Pres"])
    (b_sim - b_k) / sqrt(1 + a^2) / 1000 # metres -> km; signed, same convention as empirical delta
  }, numeric(1))

  n.degenerate <- sum(is.na(disp.sim.km))
  if (n.degenerate > 0) {
    cat("  Warning:", n.degenerate, "of", n.sims, "reps had no separating threshold in bbox.fit (dropped)\n")
  }

  ######## empirical comparison: delta[tt] = displacement from years.all[tt] to years.all[tt+1] ########
  disp.emp.km <- post.samples[, paste0("delta[", tt, "]")]

  null.results[[tt]] <- list(year = yr,
                              sim = disp.sim.km,
                              empirical = disp.emp.km,
                              n.degenerate = n.degenerate,
                              n.sims = n.sims)
}
names(null.results) <- sim.years

######## assemble comparison data and plot ########
comparison.df <- do.call(rbind, lapply(null.results, function(x) {
  rbind(
    data.frame(year = x$year, source = "empirical", distance_km = abs(x$empirical)),
    data.frame(year = x$year, source = "simulated (null)", distance_km = abs(x$sim[!is.na(x$sim)]))
  )
}))

dir.create("out/null-model", showWarnings = FALSE)
save(null.results, comparison.df, file = "out/null-model/front-displacement-null-vs-empirical.RData")

print(
  ggplot(comparison.df, aes(x = factor(year), y = distance_km, fill = source)) +
    geom_violin(alpha = 0.6, colour = NA, position = position_dodge(width = 0.8)) +
    geom_boxplot(width = 0.1, outlier.shape = NA, position = position_dodge(width = 0.8)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    labs(x = "Starting year", y = "Perpendicular front movement (km)", fill = NULL) +
    theme_bw()
)
ggsave("out/null-model/front-displacement-null-vs-empirical.pdf", width = 160, height = 120, units = "mm")

######## degenerate-rep summary across all years ########
degen <- vapply(null.results, `[[`, numeric(1), "n.degenerate")
if (any(degen > 0)) {
  cat("\nDegenerate (whole-bbox-one-sided) reps by year:\n")
  print(setNames(degen, sim.years))
}

######################################################################
######## ACTUAL-RAINFALL SCENARIO: how much of the empirical/null ########
######## discrepancy is explained by real interannual rainfall?   ########
######################################################################
# Everything below is a second, fully independent set of sims that swaps the static
# climatological rain-days-per-year value used above (rainfall.raster / all.points$rainfall,
# a single long-run-average value per point, loaded in load-data.R and reused unchanged
# for every simulated year) for the ACTUAL number of >=1mm rain-days recorded at each
# point during its relevant wet season (1 Jul of the starting year -- 30 Jun of the
# following year -- the same "one generation = one wet season" convention documented at
# the top of this script), sourced from SILO. Nothing above this point is touched or
# reused for global state (spread.table/pairs get rebuilt fresh below); the two scenarios
# are combined only at the very end, into one three-way comparison plot.
library(weatherOz)

######## unique SILO DataDrill grid cells covering all.points (dedup to keep API calls tractable) ########
# all.points covers thousands of waterpoints, but SILO's DataDrill product is a 0.05-degree
# (~5km) interpolated grid -- querying every point individually would be thousands of
# redundant calls for the same cell, so points are rounded onto the grid and only unique
# cells are queried.
silo.cell.res <- 0.05
point.lonlat <- all.points |>
  st_transform(crs = 4326) |>
  st_coordinates()
all.points$silo_lon <- round(point.lonlat[, "X"] / silo.cell.res) * silo.cell.res
all.points$silo_lat <- round(point.lonlat[, "Y"] / silo.cell.res) * silo.cell.res
# fixed-width key: round()*0.05 can leave floating-point crumbs (122.05000000000001), which
# would make the cache key representation-dependent across sessions/joins
all.points$silo_cell <- sprintf("%.2f_%.2f", all.points$silo_lon, all.points$silo_lat)

unique.cells <- all.points |>
  st_drop_geometry() |>
  distinct(silo_cell, silo_lon, silo_lat)
cat("\n", nrow(unique.cells), "unique SILO grid cells covering", nrow(all.points), "points\n")

######## wet-season windows needed: one per sim.years entry (Jul yr -- Jun yr+1) ########
wet.season.start <- function(yr) as.Date(paste0(yr, "-07-01"))
wet.season.end   <- function(yr) as.Date(paste0(yr + 1, "-06-30"))
fetch.start <- wet.season.start(min(sim.years))
fetch.end   <- wet.season.end(max(sim.years))

######## fetch (with an on-disk cache) daily rainfall per unique grid cell ########
SILO_KEY <- Sys.getenv("SILO_API_KEY")
silo.cache.path <- "dat/silo-actual-rain-days-cache.RData"
if (file.exists(silo.cache.path)) {
  # brings in `silo.daily.cache` (one row per silo_cell/day already fetched) and
  # `silo.cache.range` (the date span those rows were fetched over)
  load(silo.cache.path)
  # a cached cell is only reusable if its rows actually span the range we now need -- otherwise
  # (e.g. years.all has gained a year since the cache was written) cells would look "already
  # fetched" while silently missing days, truncating this year's rain-day counts
  if (!exists("silo.cache.range") || silo.cache.range["start"] > fetch.start ||
      silo.cache.range["end"] < fetch.end) {
    cat("SILO cache does not cover", as.character(fetch.start), "--", as.character(fetch.end),
        "; discarding and refetching\n")
    silo.daily.cache <- data.frame()
  }
} else {
  silo.daily.cache <- data.frame()
}

cached.cells <- if (nrow(silo.daily.cache) > 0) unique(silo.daily.cache$silo_cell) else character(0)
to.fetch <- unique.cells |> filter(!silo_cell %in% cached.cells)
cat(nrow(to.fetch), "grid cells need fetching from SILO (", length(cached.cells), "already cached)\n")

fetch_silo_cell <- function(cell, lon, lat) {
  cat("  Fetching", cell, "\n")
  tryCatch(
    get_data_drill(longitude = lon, latitude = lat,
                    start_date = fetch.start, end_date = fetch.end,
                    values = "rain", api_key = SILO_KEY) |>
      as.data.frame() |>
      # date rebuilt from year/month/day rather than taken from SILO's own `date` column,
      # whose returned format (Date vs YYYYMMDD string) is ambiguous and would silently
      # break as.Date() downstream
      select(year, month, day, rainfall) |>
      mutate(silo_cell = cell),
    error = function(e) {
      warning("SILO fetch failed for ", cell, ": ", conditionMessage(e))
      NULL
    }
  )
}

if (nrow(to.fetch) > 0) {
  new.daily <- purrr::pmap(list(to.fetch$silo_cell, to.fetch$silo_lon, to.fetch$silo_lat), fetch_silo_cell) |>
    bind_rows()
  silo.daily.cache <- bind_rows(silo.daily.cache, new.daily)
  silo.cache.range <- c(start = fetch.start, end = fetch.end)
  save(silo.daily.cache, silo.cache.range, file = silo.cache.path)
}
silo.daily.cache <- silo.daily.cache |>
  mutate(date = as.Date(sprintf("%04d-%02d-%02d", year, month, day)))

######## aggregate daily rainfall into per-wet-season rain-day counts (>=1mm), per cell ########
# ">=1mm" matches the threshold already used by the static climatological raster
# (rdann-1.asc, "days where at least 1mm of rain falls" -- see load-data.R) so the two
# scenarios are otherwise like-for-like.
rain.days.by.cell.year <- purrr::map_dfr(sim.years, function(yr) {
  silo.daily.cache |>
    filter(date >= wet.season.start(yr), date <= wet.season.end(yr)) |>
    group_by(silo_cell) |>
    summarise(rain.days = sum(rainfall >= 1, na.rm = TRUE), .groups = "drop") |>
    mutate(year = yr)
})

######## clip to the dispersal-kernel-fits table's valid range ########
# setup() uses rain_to_days()'s output as a direct row index into `fits`
# (dat/Kernel-fits_truncated.RData); rain-days values above a certain point map to a
# row beyond the table (fits[u, ] silently returns NA kernel params for that point).
# The static climatological raster never seems to reach this range, but a real wet
# wet-season plausibly could, so any point-year's rain-days input is capped at the
# largest value that stays in range, and how often this bites is logged per year below.
load("dat/Kernel-fits_truncated.RData") # brings `fits` back into scope (setup() only loads it into its own local frame)
max.kernel.index <- nrow(fits)
max.rain.days <- max((0:365)[rain_to_days(0:365) <= max.kernel.index])
cat("Actual rain-days values above", max.rain.days, "/season will be clipped (kernel-fits table has",
    max.kernel.index, "rows)\n")

######## second simulation loop: actual (per-year, per-point) rainfall scenario ########
# Mirrors the static-rainfall loop above, but must rebuild spread.table/pairs (via a
# fresh setup() call) every year, since the kernel now varies by year -- unlike the
# static-rainfall loop, which safely reuses one spread.table/pairs built once outside
# the loop because rainfall never changed there.
null.results.actual.rain <- vector("list", length(sim.years))

for (tt in seq_along(sim.years)) {
  yr <- sim.years[tt]
  cat("\n==== [actual rainfall] Year", yr, "====\n")

  ######## this year's actual rain-days per point, clipped, joined onto all.points ########
  yr.rain.days <- rain.days.by.cell.year |> filter(year == yr) |> select(silo_cell, rain.days)
  all.points.yr <- all.points |>
    left_join(yr.rain.days, by = "silo_cell") |>
    mutate(rainfall = pmin(rain.days, max.rain.days))
  # the join must not change row count/order -- present_k below is indexed by all.points row,
  # so a duplicated silo_cell key would silently misalign colonisation state against geometry
  stopifnot(nrow(all.points.yr) == nrow(all.points))

  n.missing <- sum(is.na(all.points.yr$rainfall))
  if (n.missing > 0) {
    cat("  Warning:", n.missing, "points had no SILO data for this wet season; falling back to climatological rainfall\n")
    all.points.yr$rainfall[is.na(all.points.yr$rainfall)] <- all.points$rainfall[is.na(all.points.yr$rainfall)]
  }
  n.clipped <- sum(yr.rain.days$rain.days > max.rain.days, na.rm = TRUE)
  if (n.clipped > 0) cat("  Note:", n.clipped, "point-years clipped to", max.rain.days, "rain-days (kernel-fits table range)\n")

  ######## rebuild spread table + dispersal matrix for this year's actual rainfall ########
  setup(point.data = all.points.yr,
        X.id = "X",
        Y.id = "Y",
        present.id = "present_k",
        artificial.natural.id = "origin_des",
        rain.id = "rainfall",
        threshold = 100,
        constant.rain = NULL,
        trunc.dist = TRUE,
        TCZ = FALSE)

  # this year's own kept.ids -- not reused from the static-rainfall loop above, since
  # this is a fresh setup() call (dedup only depends on X/Y/threshold so in practice
  # it will match, but recomputing costs nothing and removes the need to assume that)
  kept.ids.yr <- spread.table[, "ID"]

  ######## starting colonisation state for this year (identical logic to the static loop) ########
  present_k_full <- ifelse(!is.na(all.points$colonisation_year) &
                              all.points$colonisation_year <= yr, 1L, 0L)
  present_k <- present_k_full[kept.ids.yr]
  spread.table[, "Pres"] <- present_k
  spread.table[, "age"] <- present_k
  spread.table[, "occ.years"] <- present_k
  spread.table[, "first.col.gen"] <- ifelse(present_k == 1, 0, NA_real_)

  ######## sanity check: threshold-fit should recover ~b_k from the (deterministic) starting state ########
  b_k <- b_metres(tt)
  in.bbox0 <- spread.table[, "X"] >= bbox.fit["xmin"] & spread.table[, "X"] <= bbox.fit["xmax"] &
    spread.table[, "Y"] >= bbox.fit["ymin"] & spread.table[, "Y"] <= bbox.fit["ymax"]
  D0 <- spread.table[, "Y"] - a * (spread.table[, "X"] - mean.coord_m["X"])
  b0.check <- fit_b_threshold(D0[in.bbox0], spread.table[in.bbox0, "Pres"])
  cat("  Sanity check: b recovered from starting state =", round(b0.check),
      "m vs fitted b_k =", round(b_k), "m (diff =", round(b0.check - b_k), "m)\n")

  ######## run one generation (= one wet season) forward, many reps ########
  sim_out <- run_sims(n.sims = n.sims, gens = 1, plot = FALSE, rollup = FALSE,
                       stop.on.target = FALSE)

  ######## fit b_sim per rep, restricted to the empirical fitting bbox ########
  disp.sim.km <- vapply(sim_out, function(rep) {
    pm <- rep$popmatrix
    in.bbox <- pm[, "X"] >= bbox.fit["xmin"] & pm[, "X"] <= bbox.fit["xmax"] &
      pm[, "Y"] >= bbox.fit["ymin"] & pm[, "Y"] <= bbox.fit["ymax"]
    pm <- pm[in.bbox, , drop = FALSE]
    D <- pm[, "Y"] - a * (pm[, "X"] - mean.coord_m["X"])
    b_sim <- fit_b_threshold(D, pm[, "Pres"])
    (b_sim - b_k) / sqrt(1 + a^2) / 1000 # metres -> km; signed, same convention as empirical delta
  }, numeric(1))

  n.degenerate <- sum(is.na(disp.sim.km))
  if (n.degenerate > 0) {
    cat("  Warning:", n.degenerate, "of", n.sims, "reps had no separating threshold in bbox.fit (dropped)\n")
  }

  disp.emp.km <- post.samples[, paste0("delta[", tt, "]")]

  null.results.actual.rain[[tt]] <- list(year = yr,
                                          sim = disp.sim.km,
                                          empirical = disp.emp.km,
                                          n.degenerate = n.degenerate,
                                          n.sims = n.sims,
                                          n.clipped = n.clipped,
                                          n.missing.silo = n.missing)
}
names(null.results.actual.rain) <- sim.years

######## combine both scenarios with empirical into one three-way comparison ########
comparison.df.actual <- do.call(rbind, lapply(null.results.actual.rain, function(x) {
  data.frame(year = x$year, source = "simulated (actual rainfall)",
             distance_km = abs(x$sim[!is.na(x$sim)]))
}))

comparison.df.all <- rbind(comparison.df, comparison.df.actual)
comparison.df.all$source <- factor(comparison.df.all$source,
                                    levels = c("empirical", "simulated (null)", "simulated (actual rainfall)"))

dir.create("out/null-model", showWarnings = FALSE)
save(null.results, null.results.actual.rain, comparison.df.all,
     file = "out/null-model/front-displacement-null-vs-empirical-vs-actual-rainfall.RData")

print(
  ggplot(comparison.df.all, aes(x = factor(year), y = distance_km, fill = source)) +
    geom_violin(alpha = 0.6, colour = NA, position = position_dodge(width = 0.8)) +
    geom_boxplot(width = 0.1, outlier.shape = NA, position = position_dodge(width = 0.8)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    labs(x = "Starting year", y = "Perpendicular front movement (km)", fill = NULL) +
    theme_bw()
)
ggsave("out/null-model/front-displacement-three-way.pdf", width = 180, height = 120, units = "mm")

######## degenerate-rep and clipping summary across all years (actual-rainfall scenario) ########
degen.actual <- vapply(null.results.actual.rain, `[[`, numeric(1), "n.degenerate")
if (any(degen.actual > 0)) {
  cat("\n[actual rainfall] Degenerate (whole-bbox-one-sided) reps by year:\n")
  print(setNames(degen.actual, sim.years))
}
clipped.actual <- vapply(null.results.actual.rain, `[[`, numeric(1), "n.clipped")
cat("\n[actual rainfall] Point-years clipped to kernel-fits range, by year:\n")
print(setNames(clipped.actual, sim.years))

