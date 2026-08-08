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

######## load simulation posteriors ########
load("dat/Posteriors_2026.RData")

######## load empirical front parameters (produced by invasion-front-monitoring) ########
load(file.path(Sys.getenv("DATA_PATH"), "invasion-front-parameters-multi-year.Rdata"))
# objects brought into scope: mod.multi, post.samples, bbox.fit, mean.coord, scale, years.all, run.info

a <- mod.multi[[1]]["a", 1]        # posterior-mean slope (shared across years, fixed across reps)
mean.coord_m <- mean.coord * scale # mean.coord (km, model space) -> Albers metres

# posterior-mean b (Albers metres) for years.all[tt] -- same transform score_colonised_multi_year() uses
b_metres <- function(tt) {
  b_raw <- mod.multi[[1]][paste0("b[", tt, "]"), 1]
  (b_raw + mean.coord["Y"]) * scale
}

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
n.sims <- 200 # sim reps per year -- raise for a smoother null distribution (runtime scales linearly)

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
