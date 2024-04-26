source("/Users/ben/evo-dispersal/art_wbdies/VRD/pprocess_functions.R")

setwd("/Users/ben/Documents/Papers/Artificial Waterbodies/VRD")
d<-read.csv("all_watercourses_springs_lakes_holes_bores_CLIP_GDA94_RAIN_UTM.csv")
pairs.art<-pdist.fast(d$POINT_X, d$POINT_Y, 100000, 8300000)
save(pairs.art, file="vrd_pairs.RData")
