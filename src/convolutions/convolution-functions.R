require(MASS)

# generates angles based on a correlated random walk
# initial angle is random, following angles change preceding angle by draw from ang_disps
cor_angle<-function(ang_disps, n_days, samp_n=10000){
	init<-runif(samp_n, min=0, max=2*pi)
	delta<-matrix(sample(ang_disps, samp_n*(n_days-1), replace=T), nrow=samp_n)
	delta<-cbind(init, delta)
	delta<-t(apply(delta,1,cumsum))%%(2*pi)
	delta
}
# resamples scalar distances samp_n times, creates vectors based on random initial
# 	direction and subsequent correlated walk (defined by ang_disps data) over n_days.
#	Calculates displacements over n_days
nday<-function(sc_dists, n_days, ang_disps, samp_n=10000){
	r.scdists<-sample(sc_dists, samp_n*n_days, replace=T)
	r.scdists<-matrix(r.scdists, nrow=samp_n)
	r.angles<-cor_angle(ang_disps, n_days, samp_n)
	X<-r.scdists*cos(r.angles)
	Y<-r.scdists*sin(r.angles)
	X<-t(apply(X, 1, cumsum)) #cumulant across days
	Y<-t(apply(Y, 1, cumsum))
	r.scdists<-sqrt(X^2+Y^2) # cumulative displacement over days
	r.scdists
}

#Stuart's 2D -> 1D collapse of a 2D pdf that changes between cauchy and normal distributions dependent on the shape parameter.
collapse<-function(x, u, v) {
	x*u^v*v*sqrt(v^v*(u^2*v+x^2)^(-2-v))
}

collapse_gaus<-function(x, u) {  
	(exp(-x^2/(2*u^2))/(2*pi*u^2))*(2*pi*x)
}

# Function to return the area under the curve at a given truncation distance for a vector of u, v
kernel.truncation<-function(u, v, trunc.dist) {
  area <- rep(NA,length(u))
  for (ii in 1:length(u)) {
    u1 <- u[ii]
    v1 <- v[ii]
    td1 <- trunc.dist[ii]
    integrand <- function(x) {(x*u1^v1*v1*sqrt(v1^v1*(u1^2*v1+x^2)^(-2-v1)))} # 1D kernel
    area[ii] <- integrate(integrand, lower = 0, upper = td1)$value
  }
  return(area)
} 

# fits collapse to resamples of data from nday
nwise<-function(d.dists, init.v){
	init.u<-mean(d.dists[,1])
	#init.v<-0.5
	out<-c()
	for (i in 1: length(d.dists[1,])){ #step through so as to adaptively optimize starting values
		if (init.v<35) {
			temp<-fitdistr(d.dists[,i], collapse, start=list(u=init.u, v=init.v), control=list(maxit=1000))
			out<-rbind(out, c(temp[[1]], temp[[4]]))
			init.u<-temp[[1]][1]
			init.v<-temp[[1]][2]
			print(c(i, out[i,]))
		}
		else {
			temp<-fitdistr(d.dists[,i], collapse_gaus, start=list(u=init.u), control=list(maxit=1000))
			out<-rbind(out, c(temp[[1]], NA, temp[[4]]))
			init.u<-temp[[1]][1]
			print(c(i, out[i,]))
		}
	}
	colnames(out)<-c("u", "v", "LL")
	out	
}

