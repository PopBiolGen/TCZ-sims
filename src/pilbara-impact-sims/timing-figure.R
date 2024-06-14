load("out/timing-to-pilbara_likely-case.RData")

time.vec <- unlist(lapply(output, FUN = function(x){c(x$gen)}))
time.vec <- time.vec + 1 # bring to 2023 start date
arrival.year <- 2023+time.vec

pdf(file = "out/timing-to-pilbara_likely-case_boxplot.pdf")
  boxplot(arrival.year, horizontal = TRUE, xlab = "Arrival year")
dev.off()
