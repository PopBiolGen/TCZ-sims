load("out/do-nothing.RData")

time.vec <- unlist(lapply(output, FUN = function(x){c(x$gen)}))
time.vec <- time.vec + 1 # bring to 2023 start date
arrival.year <- 2023+time.vec

other.ests <- c(2035, 2041) # arrival estimates from extrapolation methods (see appendix 1)

pdf(file = "out/timing-to-pilbara_likely-case_boxplot.pdf", width = 10, height = 5)
  boxplot(arrival.year, 
          horizontal = TRUE, 
          xlab = "Arrival year",
          ylim = c(2030, max(arrival.year)))
  points(other.ests, rep(1, 2), pch = c(17, 19))
  text(x = c(other.ests), y = rep(1,2), labels = c("Constant speed", "Rain-dependent speed"), pos = 1)
  text(x = median(arrival.year), y = 0.75, labels = c("Stochastic model"))
dev.off()
