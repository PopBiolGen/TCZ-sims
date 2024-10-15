load("out/tcz-arrival.RData")

time.vec <- unlist(lapply(output, FUN = function(x){c(x$gen)}))
time.vec <- time.vec + 1 # bring to 2023 start date
arrival.year <- 2023+time.vec

source("src/pilbara-impact-sims/spread-rate-rainfall.R") # estimate arrival time taking into account only variation in rainfall

other.ests <- c(153/43+2024, tt.tcz+2024) # arrival estimates from extrapolation methods (see pilbara impact paper Appendix 1 for details)

pdf(file = "out/timing-to-tcz_boxplot.pdf", width = 10, height = 5)
  boxplot(arrival.year, 
          horizontal = TRUE, 
          xlab = "Arrival year",
          ylim = c(2024, max(arrival.year)))
  points(other.ests, rep(1, 2), pch = c(17, 19))
  text(x = c(other.ests), y = rep(1,2), labels = c("Constant speed", "Rain-dependent speed"), pos = 1)
  text(x = median(arrival.year), y = 0.75, labels = c("Stochastic model"))
dev.off()
