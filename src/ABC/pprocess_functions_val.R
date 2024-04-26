#takes x, y, indices of a matrix and places them into a vector
flatmat<-function(x, y, size){ #size is size of one side of matrix
	(y-1)*size+x
}


neighbours.init<-function(space.size, cell.size){
	lth<-space.size/cell.size
	neigh<-array(0, dim=c(9,2,lth,lth))	
	for (x in 1:lth){
		for (y in 1:lth){
			out<-matrix(nrow=9, ncol=2)
			ifelse((x-1)==0, X<-NA, X<-x-1) #correct for indices going to zero
			ifelse((y-1)==0, Y<-NA, Y<-y-1)
			ifelse((y+1)>lth, Yp<-NA, Yp<-y+1) # correct for indices going to > space
			ifelse((x+1)>lth, Xp<-NA, Xp<-x+1)
			# Assign neighbours
			out[1,]<-c(x, y) #record target cell
			out[2,]<-c(Xp, y)
			out[3,]<-c(X, y)
			out[4,]<-c(x, Y)
			out[5,]<-c(Xp, Y)	
			out[6,]<-c(X, Y)
			out[7,]<-c(x, Yp)
			out[8,]<-c(Xp, Yp)
			out[9,]<-c(X, Yp)
			neigh[,,x,y]<-out
		}	
	}
	neigh	#return the array
}


pdist.fast<-function(X, Y, maximum, space.size){
	if (min(X)<0) X<-X-min(X) # move negative coordinates to positive land
	if (min(Y)<0) Y<-Y-min(Y)
	lth<-space.size/maximum
	if (lth%%1!=0) {
		print("Error: maximum must divide into space.size perfectly")
		return(NULL)	
	}
	out<-vector("list", length=length(X)) #list to take neightbour.ID and distance
	points<-vector("list", length=lth^2)#matrix of lists of point IDs
	gX<-X%/%maximum+1 #collapse to grid refs
	gY<-Y%/%maximum+1
	neigh<-neighbours.init(space.size, maximum)
	for (i in 1:length(X)){ #throw point IDs into grid cell list
		temp<-flatmat(gX[i], gY[i], lth)
		points[[temp]]<-c(points[[temp]], i)
	}
	for (i in 1:length(X)){
		nb<-neigh[,,gX[i], gY[i]] #find neighbouring cells
		nb<-subset(nb, is.na(apply(nb,1,sum))==F)
		temp<-flatmat(nb[,1], nb[,2], lth) #find relevant points
		snk.ID<-unlist(points[temp])
		dists<-sqrt((X[i]-X[snk.ID])^2+(Y[i]-Y[snk.ID])^2)
		temp<-cbind(snk.ID, dists)
		temp<-subset(temp, temp[,"dists"]<=maximum)
		out[[i]]<-temp
	}	
	out
}




#The kernel
dcncross<-function(x, u, v) {  #stuart's cauchy-normal distribution in 2D
	(x*u^v*v*sqrt(v^v*(u^2*v+x^2)^(-2-v)))/(2*pi*x)
} 

# outputs data.  With or without a plot as well.

output<-function(pop, gen, id, plot=T){
	dname<-paste(id, "pop", gen, ".txt", sep="")
	write.table(pop, dname, sep="\t", row.names=F)
	if (plot==T){
		pname<-paste(id, "plot", gen, ".png", sep="")
		plotter(pop, pname)
	}	
}

# spreads the population over gen generations
spread<-function(pop, gens, pairs, K, obs){ #pairs is a list from pdist.fast
	preds<-NULL	
for (i in 1:gens){
		#if (i%%5==0) output(pop, i, K)
		occp<-subset(pop, pop[,"Pres"]==1 & pop[,"age"]<51) #collect occupied sites
		if (length(occp[,1])<1) {
			print(paste("No more vacant opportunities at generation", i))
			break()	
		}
		potl<-pairs[occp[,"ID"]] #collect relevant parts of pair list
		potl<-do.call("rbind", potl)
		src.ID<-rep(occp[,"ID"], times=occp[,"n.pairs"])
		potl<-cbind(src.ID, potl)
		U<-rep(occp[,"u"], times=occp[,"n.pairs"]) #expand source specific kernel parameters
		V<-rep(occp[,"v"], times=occp[,"n.pairs"])
		recruits<-dcncross(potl[,"dists"]+0.05, U, V) #calculate probabilities of arrival for each pair
    n.inds<-sum(recruits)*K #Calculates number of individuals to draw
    no.recruits<-1-recruits
		recruits<-1-tapply(no.recruits, potl[,"snk.ID"], vec.prod) #vector of probabilities for the multinomial
		recruit.ID<-as.integer(row.names(recruits))
		recruits<-rmultinom(1, n.inds, recruits)
		recruits<-cbind(recruit.ID, recruits)
		colnames(recruits)[2]<-"recruits"
		recruits<-subset(recruits, recruits[,"recruits"]>2)
		pop[match(recruits[,"recruit.ID"], pop[,"ID"]), "Pres"]<-1
		pop[which(pop[,"Pres"]==1), "age"]<-1+pop[which(pop[,"Pres"]==1), "age"]
	preds<-rbind(preds,output_val(pop, i))
	}
		obs_pred_cf(preds,obs)
}

# calculates the product of all elements in a vector
vec.prod<-function(vec){
  last<-length(vec)
  cumprod(vec)[last]
}

#plots opportunities and colonised populations
plotter<-function(popmatrix, file.name){
	png(filename=file.name, width=7, height=7, units="cm", res=600, pointsize=6)
	plot(popmatrix[,2], popmatrix[,3], xlab="False easting (kms)", ylab="False northing (kms)", pch=19)
	occp<-subset(popmatrix, popmatrix[,"Pres"]==1)
	points(occp[,2], occp[,3], pch=19, col="red")
	dev.off()
}


output_val<-function(pop, gen){
		pop<-pop[,1:4]
		pop<-cbind(pop,rep(gen, nrow(pop)))	
		pop
}

obs_pred_cf<-function(preds, obs){
	preds[,5]<-preds[,5]+2007
	temp<-merge(preds, obs, by.x=c(2, 3, 5), by.y=c(4, 5, 10))
	#browser()
	temp2<-ifelse(temp[,"Pres"]==temp[,"OCCUPIED"],1,0)
	sum(temp2)
}
