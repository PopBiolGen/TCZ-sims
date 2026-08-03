# Synthetic demo/regression script for the control-step (management-driven extinction) mechanics
# added to src/pprocess_functions.R. Uses made-up point data (no DATA_PATH dependency) so it can be
# run standalone to sanity-check the new setup()/spread.pilb() arguments before wiring real TCZ
# infrastructure-audit data through load-data.R.
#
# Three checks:
#   A. Regression: control.step=FALSE reproduces today's permanent-occupation behavior exactly,
#      even though the new columns are present in spread.table.
#   B. Attrition: with control.step=TRUE, tank/fence points go extinct at their expected annual
#      rates (1%/5%) and natural points (default surv.prob=1) never do.
#   C. Fence breach-thinning: fenced target points (r.eff == r, isolating p.breach's effect) get
#      colonised less often than otherwise-identical unfenced targets at the same distance.

source("src/pprocess_functions.R")
set.seed(20260730)

r <- 300 # detection radius (m), matching the qmd's assumed value

######## Dataset 1: attrition test (tank / fence / natural points, all pre-colonised) ########
n.per.type <- 800
make_grid <- function(n, spacing = 1000, offset = 0) {
  side <- ceiling(sqrt(n))
  g <- expand.grid(gx = seq_len(side), gy = seq_len(side))[seq_len(n), ]
  data.frame(X = g$gx * spacing + offset, Y = g$gy * spacing)
}

tank.pts    <- cbind(make_grid(n.per.type, offset = 0),       control_type = "tank",  art_nat = 1)
fence.pts   <- cbind(make_grid(n.per.type, offset = 1e6),     control_type = "fence", art_nat = 1)
natural.pts <- cbind(make_grid(n.per.type, offset = 2e6),     control_type = "none",  art_nat = 0)

attrition.data <- rbind(tank.pts, fence.pts, natural.pts)
attrition.data$present     <- 1 # all pre-colonised
attrition.data$fence_area  <- ifelse(attrition.data$control_type=="fence", pi*r^2, NA) # r.eff == r
attrition.data$fence_units <- ifelse(attrition.data$control_type=="fence", 54, 0)      # p.breach ~ 0.5 (unused here, no new colonisation)

setup(point.data = attrition.data,
      X.id = "X", Y.id = "Y", present.id = "present", artificial.natural.id = "art_nat",
      remove_duplicates = FALSE, constant.rain = 100, TCZ = FALSE,
      control.type.id = "control_type", fence.area.id = "fence_area", fence.units.id = "fence_units")
attrition.table <- spread.table
attrition.pairs <- pairs

cat("\n== A. Regression check (control.step=FALSE) ==\n")
out.A <- spread.pilb(pop = attrition.table, gens = 10, pairs = attrition.pairs,
                      delta = 0, r = r, rollup = FALSE, control.step = FALSE)
cat("All points still colonised:", all(out.A$popmatrix[, "Pres"]==1), "\n")
cat("occ.years == age for every point:", all(out.A$popmatrix[, "occ.years"]==out.A$popmatrix[, "age"]), "\n")
cat("age advanced by 10 (from its initial value of 1) for every point:", all(out.A$popmatrix[, "age"]==11), "\n")

cat("\n== B. Attrition rates (control.step=TRUE, delta=0, 1 generation) ==\n")
# control.fail.prob is P(control fails, i.e. toads persist), so extinction probability is
# 1-control.fail.prob: ~99% for tanks (1% failure rate per the qmd), ~95% for fences (5% eradication
# failure rate) -- both are deliberately HIGH per-generation extinction probabilities; it's the fence's
# breach-thinning at the colonisation step (tested separately below) that limits how often a fence
# point gets colonised in the first place.
out.B <- spread.pilb(pop = attrition.table, gens = 1, pairs = attrition.pairs,
                      delta = 0, r = r, rollup = FALSE, control.step = TRUE)
pm <- out.B$popmatrix
rate <- function(type) 1 - mean(pm[pm[,"control.type"]==type, "Pres"])
cat(sprintf("Tank extinction rate:    %.4f (expected ~0.99)\n", rate(CONTROL_TANK)))
cat(sprintf("Fence extinction rate:   %.4f (expected ~0.95)\n", rate(CONTROL_FENCE)))
cat(sprintf("Natural extinction rate: %.4f (expected 0, surv.prob defaults to 1)\n", rate(CONTROL_NONE)))

cat("\n== occ.years / first.col.gen sanity (control.step=TRUE, 10 generations) ==\n")
out.C <- spread.pilb(pop = attrition.table, gens = 10, pairs = attrition.pairs,
                      delta = 0, r = r, rollup = FALSE, control.step = TRUE)
pm <- out.C$popmatrix
cat("occ.years <= age for every point:", all(pm[, "occ.years"] <= pm[, "age"]), "\n")
cat("first.col.gen is 0 for all (all pre-colonised at t=0):", all(pm[, "first.col.gen"]==0), "\n")

cat("\n== rollup + control.step incompatibility is rejected ==\n")
rejected <- tryCatch({
  spread.pilb(pop = attrition.table, gens = 1, pairs = attrition.pairs, delta = 0, r = r,
              rollup = TRUE, control.step = TRUE)
  FALSE
}, error = function(e) TRUE)
cat("stop() raised when rollup=TRUE & control.step=TRUE:", rejected, "\n")

######## Dataset 2: fence breach-thinning test ########
n.targets <- 30
angles <- seq(0, 2*pi, length.out = n.targets*2 + 1)[seq_len(n.targets*2)]
dist.to.seed <- 1000

seed.pt      <- data.frame(X = 0, Y = 0, control_type = "none", art_nat = 1, present = 1,
                           fence_area = NA, fence_units = 0)
unfenced.pts <- data.frame(X = dist.to.seed*cos(angles[1:n.targets]),
                            Y = dist.to.seed*sin(angles[1:n.targets]),
                            control_type = "none", art_nat = 1, present = 0,
                            fence_area = NA, fence_units = 0)
fenced.pts   <- data.frame(X = dist.to.seed*cos(angles[(n.targets+1):(2*n.targets)]),
                            Y = dist.to.seed*sin(angles[(n.targets+1):(2*n.targets)]),
                            control_type = "fence", art_nat = 1, present = 0,
                            fence_area = pi*r^2, fence_units = 54) # r.eff == r, p.breach ~ 0.5

breach.data <- rbind(seed.pt, unfenced.pts, fenced.pts)

setup(point.data = breach.data,
      X.id = "X", Y.id = "Y", present.id = "present", artificial.natural.id = "art_nat",
      remove_duplicates = FALSE, constant.rain = 100, TCZ = FALSE,
      control.type.id = "control_type", fence.area.id = "fence_area", fence.units.id = "fence_units")
breach.table <- spread.table
breach.pairs <- pairs
p.breach.val <- unique(breach.table[breach.table[,"control.type"]==CONTROL_FENCE, "p.breach"])
cat(sprintf("\n== C. Fence breach-thinning test (p.breach = %.3f, r.eff == r == %d) ==\n", p.breach.val, r))

# rows are in input order (remove_duplicates=FALSE): 1 = seed, 2:(n.targets+1) = unfenced targets,
# (n.targets+2):(2*n.targets+1) = fenced targets
unfenced.idx <- breach.table[, "ID"] %in% 2:(n.targets+1)
fenced.idx   <- breach.table[, "ID"] %in% (n.targets+2):(2*n.targets+1)
stopifnot(all(breach.table[unfenced.idx, "control.type"]==CONTROL_NONE),
          all(breach.table[fenced.idx, "control.type"]==CONTROL_FENCE))

n.reps <- 300
unfenced.col <- numeric(n.reps); fenced.col <- numeric(n.reps)
for (rr in seq_len(n.reps)) {
  out <- spread.pilb(pop = breach.table, gens = 1, pairs = breach.pairs,
                      delta = 200, r = r, rollup = FALSE, control.step = TRUE)
  pm <- out$popmatrix
  unfenced.col[rr] <- mean(pm[unfenced.idx, "Pres"])
  fenced.col[rr]   <- mean(pm[fenced.idx, "Pres"])
}
cat(sprintf("Unfenced target colonisation rate: %.3f\n", mean(unfenced.col)))
cat(sprintf("Fenced target colonisation rate:   %.3f\n", mean(fenced.col)))
cat(sprintf("(fenced rate reflects two compounding effects, both measured within this same generation:\n"))
cat(sprintf(" breach-thinning reduces arrivals by ~(1 - p.breach) = %.3f, and any point that does get\n", 1 - p.breach.val))
cat(sprintf(" colonised still faces the recurring control.fail.prob=0.05 eradication test immediately;\n"))
cat(sprintf(" combined, the fenced rate should be roughly (1-p.breach)*0.05 = %.4f of the unfenced rate)\n", (1-p.breach.val)*0.05))

# age==0 iff first.col.gen is NA holds for any point never yet colonised (using the last rep's result,
# which always has some never-colonised target points since colonisation isn't certain)
cat("age==0 exactly where first.col.gen is NA (last rep):",
    all((pm[, "age"]==0) == is.na(pm[, "first.col.gen"])), "\n")

