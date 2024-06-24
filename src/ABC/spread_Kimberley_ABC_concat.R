# Script to concatenate the ABC outputs and generate priors for lambda

flist <- list.files("out/ABC", pattern = "ABC_", full.names = TRUE)

# retrieve and concatenate all outputs
concat <- c()
for (ff in 1:length(flist)){
  load(file = flist[ff])
  concat <- c(concat, output)
}

# get lambda values
l.vals <- unlist(lapply(concat, function(x) {x$pars[, "lambda"]}))

# get maximum score for a fit
x1 <- concat[[1]]$popmatrix
max.fit.val <- sum(is.finite(x1[, "X2023_pres"])) + sum(is.finite(x1[, "X2022_pres"]))

# function to score each simulation for match to 2022, 2023 data
score_sim <- function(x) {
  pm <- x$popmatrix
  test.2023 <- as.numeric(pm[, "age"] > 0)
  m.2023 <- sum(test.2023 == pm[, "X2023_pres"], na.rm = TRUE)
  test.2022 <- as.numeric((pm[, "age"]-1) > 0)
  m.2022 <- sum(test.2022 == pm[, "X2022_pres"], na.rm = TRUE)
  (m.2023+m.2022)/max.fit.val
}

fit.vals <- unlist(lapply(concat, score_sim))

plot(fit.vals~l.vals)

# take top 5% of sims
cut.q <- quantile(fit.vals, probs = 0.95)
samp.l <- l.vals[fit.vals>=cut.q]
samp.l.log <- log10(samp.l)
sample.lambda <- mean(samp.l.log)
sample.lambda.sd <- sd(samp.l.log)

save(sample.lambda, sample.lambda.sd, file = "dat/Posteriors_2023.RData")
rm(list=ls())
