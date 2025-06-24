# This is spread rate method 2 in the Dunlop et al paper
# Analysis to unearth relationship between ndays of rain and spread rate.
# using invasion of the Kimberley as the observation
# then applies this relationship to forecast time to arrive in the Pilbara.

library(dplyr)
library(sf)

#load data of n rainy days along transect between Lake argyle and Willare = 14 years to invade this distance
d <- read.csv(file = "dat/Invasion_pathway_rainfall.csv")
d <- st_as_sf(x = d, coords = c("X", "Y"), crs = 4326) # cast to sf

nyears <- 14

k.path <- d %>% filter(Pathway == "Kimberley")
k.path$distance <- as.numeric(cumsum(c(0,diag(st_distance(k.path[-nrow(k.path), ], k.path[-1, ])))))
total.dist.k <- max(k.path$distance)/1000 # km from Argyle to Willare

# make an interpolation function to extract n rainy days at a set of distances
ndays.k <- approxfun(x = k.path$distance, y = k.path$nDays_rain_1)
# get breakpoints for ndays
dist.vec <- seq(0, max(k.path$distance), length.out = nyears)
ndays.vec <- ndays.k(dist.vec)
# this gives us (approximately) the number of rainy days on the invasion front each year

# we know they spread total.dist in nyears.  We assume spread distance in a given year is proportional to ndays
# we wish to know what this proportionality is, so:
prop <- total.dist.k/sum(ndays.vec)

# we now apply this to the LaGrange pathway..
l.path <- d %>% filter(Pathway == "La Grange")
consecutive.distances <- diag(st_distance(l.path[-nrow(l.path), ], l.path[-1, ]))
l.path$distance <- as.numeric(c(0, cumsum(consecutive.distances))/1000)

ndays.l <- approxfun(x = l.path$distance, y = l.path$nDays_rain_1)

#time to arrive in the Pilbara
test.d <- max(l.path$distance)
spread.dist <- 0
year <- 1
tt <- 1
while (test.d > 0){
  this.yr <- prop * ndays.l(spread.dist[tt])
  spread.dist <- c(spread.dist, this.yr)
  year <- c(year, year[tt]+1)
  test.d <- test.d - this.yr
  tt <- tt+1
}
tt.pilb <- tt

cbind(year, spread.dist, total.spread.dist = cumsum(spread.dist))

#time to arrive at head of TCZ
test.d <- 153 # distance in km between invasion front and Shamrock station as of Oct 24
spread.dist <- 0
year <- 1
tt <- 1
while (test.d > 0){
  this.yr <- prop * ndays.l(spread.dist[tt])
  spread.dist <- c(spread.dist, this.yr)
  year <- c(year, year[tt]+1)
  test.d <- test.d - this.yr
  tt <- tt+1
}
tt.tcz <- tt


