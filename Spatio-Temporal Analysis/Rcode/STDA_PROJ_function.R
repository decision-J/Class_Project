library(mapproj)
library(geoR)
library(spBayes)
library(stats)
library(cubature)
library(ggplot2)
library(scales)

##############################################################################################
gam.integrand = function(t,j,s0,u,ut,sigmasq,phi,coord) {
  ##############################################################################################
  delta = unlist(s0) - coord[j,] + u%*%t
  delta.mag = sqrt(delta[1,]^2 + delta[2,]^2)
  gradK = -sigmasq*phi^2*exp(-phi*delta.mag) * t(delta)
  gradK %*% ut
}

##############################################################################################
K.integrand = function(t,s0,u,ut,sigmasq,phi,coord) {
  ##############################################################################################
  
  if(diff(t) == 0){return(0)}
  delta1 = u[1] * (t[2] - t[1])
  delta2 = u[2] * (t[2] - t[1])
  delta.mag = sqrt(delta1^2 + delta2^2)
  #if(delta.mag == 0){ return(phi*ut[1]^2 + phi*ut[2]^2) } # limit of H as distance -> 0
  H11 = -sigmasq*phi^2*exp(-phi*delta.mag)*(1-phi*delta1^2/delta.mag)
  H22 = -sigmasq*phi^2*exp(-phi*delta.mag)*(1-phi*delta2^2/delta.mag)
  -(H11*ut[1]^2 + H22*ut[2]^2)
}



##############################################################################################
flux.line = function(s0,u,ut,tstar,jump.data) {
  ##############################################################################################
  with(jump.data, {
    beta = params[1:n.beta]
    sigmasq = params[n.beta+1]
    tausq = params[n.beta+2]
    phi = params[n.beta+3] # make sure parameterization agrees w/ Banerjee
    
    u = matrix(u); ut = matrix(ut);
    mu.flux = tstar*t(ut)%*%beta[2:3] 	# correct for beta0+beta1*x+beta2*y
    
    if(dim(covariates)[1]!=0){ covar.mat = as.matrix(cbind(1,covariates)) }
    else{ covar.mat = 1 }
    mu <- covar.mat %*% beta
    
    # integrate gamma for each j
    n = length(dates)
    gamma.flux = rep(NA,n)
    for(j in 1:n) {
      gamma.flux[j] = integrate(gam.integrand,lower=0,upper=tstar,j=j,s0=s0,u=u,ut=ut,
                                sigmasq=sigmasq,phi=phi,coord=coord)$value
    }
    # integrate to find K only once
    K.flux = adaptIntegrate(K.integrand, lowerLimit=c(0,0), upperLimit=c(tstar,tstar),
                            s0=s0,u=u,ut=ut,sigmasq=sigmasq,phi=phi,coord=coord)$integral
    # get mean and variance of the total flux across the line segment
    #Kinv = solve( sigmasq * (1+phi*dmat) * exp(-phi*dmat ) + tausq*sigmasq*diag(n) )
    Kinv = solve( sigmasq * (1+phi*dmat) * exp(-phi*dmat ) + tausq*diag(n) )
    
    mean.flux = mu.flux + as.vector( -t(gamma.flux)%*%Kinv%*%(dates - mu) )
    var.flux = K.flux - t(gamma.flux)%*%Kinv%*%gamma.flux
    return(c(mean.flux,var.flux))
  })
}


##############################################################################################
flux.box = function(center,jump.data) {
  ##############################################################################################
  
  with(jump.data, {
    flux.up = flux.line(s0 = c(center[1] - r, center[2] + r), u = c(1,0), ut = c(0,1), tstar = 2*r, jump.data)
    flux.right = flux.line(c(center[1] + r, center[2] + r), c(0,-1), c(1,0), 2*r, jump.data)
    flux.down = flux.line(c(center[1] + r, center[2] - r), c(-1,0), c(0,-1), 2*r, jump.data)
    flux.left = flux.line(c(center[1] - r, center[2] - r), c(0,1), c(-1,0), 2*r, jump.data)
    
    # get 10000 posterior samples of the average flux out of the box
    #samp.up = rnorm(10000,mean=flux.up[1],sd=sqrt(flux.up[2]))
    #samp.right = rnorm(10000,mean=flux.right[1],sd=sqrt(flux.right[2]))
    #samp.down = rnorm(10000,mean=flux.down[1],sd=sqrt(flux.down[2]))
    #samp.left = rnorm(10000,mean=flux.left[1],sd=sqrt(flux.left[2]))
    #samps = list(samp.up/(2*r),samp.down/(2*r),samp.left/(2*r),samp.right/(2*r))
    
    #samp.flux = (samp.up+samp.right+samp.down+samp.left)/(8*r)
    # use xtable to print a nice table showing the mean + stdev of the average flux along each side,
    # in addition to the overall average
    
    #meangrad = sapply(samps,mean)
    #sdgrad = sapply(samps,sd)
    
    meangrad = c(flux.up[1]/(2*r),flux.down[1]/(2*r),flux.left[1]/(2*r),flux.right[1]/(2*r))
    sdgrad = c(sqrt(flux.up[2])/(2*r),sqrt(flux.down[2])/(2*r),sqrt(flux.left[2])/(2*r),sqrt(flux.right[2])/(2*r))
    #sides = c("Top","Bottom","Left","Right","Total")
    #out = data.frame(meangrad,sdgrad)
    #rownames(out) = sides
    #colnames(out) = c("Avg Gradient","StDev")
    #print(xtable(out))
    
    # return whether each side has a significantly negative gradient (-1), significantly positive gradient (+1) or non-significant gradient (0)
    sig.sides = c(0,0,0,0)
    sig.sides = sig.sides - ((meangrad + 1.96*sdgrad) < 0) # significantly small gradient (long range jump)
    sig.sides = sig.sides + ((meangrad - 1.96*sdgrad) > 0) # significantly large gradient
    return( sig.sides )
  })
}

##############################################################################################
jump.scan = function(dates, coord.longlat, newcoord.longlat=NULL, r=NULL, params, Albers = NULL, covar = NULL) {
  ##############################################################################################
  
  # without Albers equal area conic projection
  if(is.null(Albers)){
    coord = coord.longlat
    conv.factor = 1
    # Default to predicting gradient at observed locations
    if(is.null(newcoord.longlat)){
      newcoord.longlat <- coord.longlat
    }
    newcoord = newcoord.longlat
    
    # Albers equal area conic projection
  }else{
    proj = mapproject(coord.longlat[,1], coord.longlat[,2], projection="albers", parameters=Albers, orientation=NULL)
    coord = cbind(proj$x,proj$y)
    # scale to km (minimize distortion)
    conv.factor = conversion(coord,coord.longlat)
    
    # Default to predicting gradient at observed locations
    if(is.null(newcoord.longlat)){
      newcoord.longlat <- coord.longlat
    }
    
    # Albers equal area conic projection for newcoordinate
    proj.n = mapproject(newcoord.longlat[,1], newcoord.longlat[,2], projection="albers", parameters=Albers, orientation=NULL)
    newcoord = cbind(proj.n$x,proj.n$y)
  }
  
  # Construct covariates
  covariates = cbind(coord,covar)
  n.beta = dim(covariates)[2]+1
  all.sides = matrix(NA,ncol=4,nrow=dim(newcoord)[1])
  
  # compute pairwise distance matrix 'dmat' of coord
  dmat <- sqrt(outer(coord[,1], coord[,1], "-")^2 + outer(coord[,2], coord[,2], "-")^2)
  if(is.null(r)){   r <- quantile(dmat,.025)    }
  jump.data = list(newcoord=newcoord, r=r, params=params, dates=dates, coord=coord, n.beta=n.beta, covariates=covariates, dmat=dmat)
  
  for(i in 1:dim(newcoord)[1]) {
    all.sides[i,] = flux.box(newcoord[i,], jump.data)
  }
  res <- list(coord.longlat=coord.longlat,newcoord.longlat=newcoord.longlat, conv.factor=conv.factor, r=r, all.sides=all.sides)
  class(res) <- "jump.scan"
  res
}

##############################################################################################
plotboxes = function(r = NULL, object, num.sides, method = c("box", "arrow"),point=FALSE,...) {
  ##############################################################################################
  
  coord.longlat <- object[[1]]
  newcoord.longlat <- object[[2]]
  conv.factor <- object[[3]]
  jump <- object[[4]]
  all.sides <- object[[5]]
  if(is.null(r)){   r <- 1    }
  
  jump.centers = NULL
  jump.sides = NULL
  
  x0 = NULL; x1 = NULL; y0 = NULL; y1 = NULL;
  color.sides = NULL
  for(i in 1:dim(newcoord.longlat)[1]) {
    if(max(all.sides[i,]) <= 0 && sort(all.sides[i,])[num.sides] == -1) {
      jump.centers = rbind(jump.centers,newcoord.longlat[i,])
      jump.sides = rbind(jump.sides,all.sides[i,])
    }
  }
  if(length(jump.centers)==0) {
    print("No sites of potential long-range jumps.")
    return(0)
  }
  
  for(i in 1:dim(jump.centers)[1]) {
    center = jump.centers[i,]
    color.sides = c(color.sides,-jump.sides[i,] + 1) # black if nonsignificant, red if significant
    upper_x = center[1] + r
    lower_x = center[1] - r
    upper_y = center[2] + r
    lower_y = center[2] - r
    x0 = c(x0,lower_x,upper_x,lower_x,upper_x)
    x1 = c(x1,upper_x,lower_x,lower_x,upper_x)
    y0 = c(y0,upper_y,lower_y,lower_y,upper_y)
    y1 = c(y1,upper_y,lower_y,upper_y,lower_y)
    # plot in order top, bottom, left, right
  }
  color.sides[color.sides==1] = NA # remove nonsignificant segments  plotgrad(out.grad,cex=1,pch=".",database="state")
  
  if(method=="box"){               # Draw box
    segments(unlist(x0), unlist(y0), unlist(x1), unlist(y1), col=color.sides,...)}else{
      # Draw arrow
      points(jump.centers,cex=1,col=2,pch=7,lwd=2)
      colorsig <- matrix(color.sides,ncol=4,byrow=T)
      arrows(jump.centers[,1],jump.centers[,2],jump.centers[,1],y1[seq(1,length(color.sides),by=4)],col=colorsig[,1],length=0.1,...)  #1. (x,y+r)
      arrows(jump.centers[,1],jump.centers[,2],jump.centers[,1],y1[seq(2,length(color.sides),by=4)],col=colorsig[,2],length=0.1,...)  #2. (x,y-r)
      arrows(jump.centers[,1],jump.centers[,2],x1[seq(2,length(color.sides),by=4)],jump.centers[,2],col=colorsig[,3],length=0.1,...)  #3. (x-r,y)
      arrows(jump.centers[,1],jump.centers[,2],x1[seq(1,length(color.sides),by=4)],jump.centers[,2],col=colorsig[,4],length=0.1,...)  #4. (x+r,y)
    }
  
  # print
  if(conv.factor==1){
    cat("Long range jumps tested for" ,jump,".\n")
  }else{
    cat("Long range jumps tested for ",conv.factor*jump,"km.\n") # scale to km #
  }
}

##############################################################################################
localgrad <- function(dates, covar = NULL, coord.longlat, newcoord.longlat = NULL, Albers = NULL, n.samp = 10000, thin = 10, knots = NA, amcmc = NA, verbose = FALSE){
  ##############################################################################################
  
  # without Albers equal area conic projection
  if(is.null(Albers)){
    coord = coord.longlat
    conv.factor = 1
    # Default to predicting gradient at observed locations
    if(is.null(newcoord.longlat)){
      newcoord.longlat <- coord.longlat
    }
    newcoord = newcoord.longlat
    
    # Albers equal area conic projection
  }else{
    proj = mapproject(coord.longlat[,1], coord.longlat[,2], projection="albers", parameters=Albers, orientation=NULL)
    coord = cbind(proj$x,proj$y)
    # scale to km (minimize distortion)
    conv.factor = conversion(coord,coord.longlat)
    
    # Default to predicting gradient at observed locations
    if(is.null(newcoord.longlat)){
      newcoord.longlat <- coord.longlat
    }
    
    # Albers equal area conic projection for newcoordinate
    proj.n = mapproject(newcoord.longlat[,1], newcoord.longlat[,2], projection="albers", parameters=Albers, orientation=NULL)
    newcoord = cbind(proj.n$x,proj.n$y)
  }
  
  # Construct covariates
  n.obs <- length(dates)
  n.pred <- dim(newcoord)[1]
  covariates = cbind(coord,covar)
  n.beta = dim(covariates)[2]+1
  formula = dates~covariates
  
  # matrix of pairwise distances
  dmat <- sqrt(outer(coord[,1], coord[,1], "-")^2 + outer(coord[,2], coord[,2], "-")^2)
  
  # get prior estimates of covariance parameters from variogram (trend is only using coordinates)
  vario.out <- variofit(variog(coords = coord, data=dates, trend = ~ coord[,1] + coord[,2], messages=verbose), fix.kappa = TRUE, kappa = 1.5, weights="cressie",max.dist=as.numeric(quantile(dmat[dmat!=0],.5)), messages=verbose)
  tausq.vario <- vario.out$nugget
  sigmasq.vario <- vario.out$cov.pars[1]
  phi.vario <- 1/vario.out$cov.pars[2] # transform to agree w/ spBayes parameterization
  
  # provision for bad variogram fit
  if(phi.vario == 0) { phi.vario = min(dmat[dmat!=0]) }
  if(sigmasq.vario == 0) { sigmasq.vario = min(dmat[dmat!=0]) }
  if(tausq.vario == 0) { tausq.vario = sigmasq.vario/100 }
  
  # fit a Gaussian process, Matern nu=3/2 using spBayes package
  # starting value for the regression parameter
  start.beta <- c(mean(dates),rep(1,n.beta-1))
  priors <- list("phi.Unif"=c(1/max(dmat),1/min(dmat[dmat!=0])), "sigma.sq.IG"=c(2, sigmasq.vario),
                 "tau.sq.IG"=c(2, tausq.vario), "nu.unif"=c(1,2))
  starting <- c(start.beta,"phi"=phi.vario, "sigma.sq"=sigmasq.vario, "tau.sq"=tausq.vario, "nu"=1.5)
  tuning <- lapply(starting,function(x){x/100})
  tuning$nu <- 0 # fix smoothness parameter at 1.5
  
  if(is.na(knots) && is.na(amcmc)) {
    krige.out <- spLM(formula, coords=coord, starting=starting, tuning=tuning, priors=priors,
                      cov.model = "matern",  n.samples = n.samp, n.report = 500, verbose = verbose)
  } else if(!is.na(knots) && is.na(amcmc)) {
    krige.out <- spLM(formula, coords=coord, starting=starting, tuning=tuning, priors=priors,
                      cov.model = "matern",  n.samples = n.samp, n.report = 500, verbose = verbose, knots=knots)
  } else if(is.na(knots) && !is.na(amcmc)) {
    krige.out <- spLM(formula, coords=coord, starting=starting, tuning=tuning, priors=priors,
                      cov.model = "matern",  n.samples = n.samp, n.report = 500, verbose = verbose, amcmc=amcmc)
  } else {
    krige.out <- spLM(formula, coords=coord, starting=starting, tuning=tuning, priors=priors,
                      cov.model = "matern",  n.samples = n.samp, n.report = 500, verbose = verbose, knots=knots, amcmc=amcmc)
  }
  burn.in <- 0.5*n.samp
  # recover beta and spatial random effects
  krige.samp <- spRecover(krige.out, start=burn.in, verbose=verbose)
  posterior.samples <- cbind(krige.samp$p.beta.recover.samples,krige.samp$p.theta.recover.samples)
  
  
  # thin posterior samples to numsamp
  numsamp <- n.samp/thin
  thin.ind <- seq(1,dim(posterior.samples)[1],length=floor(numsamp))
  posterior.samples <- posterior.samples[thin.ind,]
  
  # draw posterior gradients from each location in newcoord following (Banerjee 2003)
  posterior.out <- list()
  for(i in 1:n.pred) {
    posterior.out[[i]] <- matrix(nrow=0,ncol=2)
  }
  if(verbose==TRUE){cat("Computing gradients for",numsamp,"samples:\n")}
  
  if(dim(covariates)[1]!=0){ covar.mat = as.matrix(cbind(1,covariates)) }else { covar.mat = 1 }
  
  for(j in 1:dim(posterior.samples)[1]) {
    beta <- posterior.samples[j,1:(n.beta)]
    sigmasq <- posterior.samples[j,n.beta+1]
    tausq <- posterior.samples[j,n.beta+2]
    phi <- posterior.samples[j,n.beta+3]
    
    Kinv <- solve( sigmasq * (1+phi*dmat) * exp(-phi*dmat) + tausq*diag(n.obs) )
    mu <- covar.mat %*% beta
    grad.mu <- beta[2:3]
    grad.mu = matrix(rep(grad.mu,n.pred),nrow=n.pred,byrow=TRUE)
    
    # now get 10 gradient samples at each location (total samples 10*numsamp)
    for(i in 1:n.pred) {
      s0 <- newcoord[i,]
      delta <- cbind(s0[1] - coord[,1], s0[2] - coord[,2])
      delta.mag <- sqrt( (s0[2]-coord[,2])^2 + (s0[1]-coord[,1])^2 )
      
      gam <- -sigmasq*phi^2*exp(-phi*delta.mag)*delta
      mean.grad <- grad.mu[i,] - as.vector( t(gam)%*%Kinv%*%(dates-mu) )
      var.grad <- sigmasq*phi^2*diag(2) - t(gam)%*%Kinv%*%gam
      
      # multivariate normal sample of size 10
      z <- matrix(rnorm(10*2),10,2) %*% chol(var.grad)
      y <- t(mean.grad + t(z))
      posterior.out[[i]] <- rbind(posterior.out[[i]], y)
    }
    if(verbose==TRUE){
      cat(10*j,",",sep="")
      if(j%%20==0) cat("\n")
    }
  }
  
  # test for significant gradients
  normal.density <- rep(NA,n.pred)
  for(i in 1:n.pred) {
    pts <- posterior.out[[i]]
    rho <- cor(pts[,1],pts[,2])
    mu.x <- mean(pts[,1]); mu.y <- mean(pts[,2])
    sd.x <- sd(pts[,1]); sd.y <- sd(pts[,2])
    normal.density[i] <- 1/(1-rho^2) * ( mu.x^2/sd.x^2 + mu.x^2/sd.x^2 - 2*rho*mu.x*mu.y/(sd.x*sd.y) )
  }
  is.sig <- normal.density > pchisq(.95,2)
  
  # transform gradients into speeds of spread (polar transformation)
  spread.means <- t(sapply( posterior.out, function(x){apply(x,2,mean)} ))
  x <- spread.means[,1]
  y <- spread.means[,2]
  r <- sqrt(x^2+y^2)
  theta <- atan(y/x) + (x<=0)*pi
  r <- 1/r
  theta <- theta + pi
  speed.means <- cbind(r*cos(theta),r*sin(theta))
  
  res <- list(call=deparse(match.call(),width.cutoff = 500), dates = dates, coord.longlat = coord.longlat, newcoord.longlat = newcoord.longlat,
              coord = coord, newcoord = newcoord, conv.factor=conv.factor, samples = posterior.out, means = speed.means, sig = is.sig, post.samp = posterior.samples)
  class(res) <- "localgrad"
  res
}



##############################################################################################
summary.localgrad<-function(object, show=FALSE, ...){
  ##############################################################################################
  means.km <- object$conv.factor*object$means # scale to km
  mag.speed <- sqrt( means.km[,1]^2 + means.km[,2]^2 )
  
  mean.speed <- mean(mag.speed)
  median.speed <- median(mag.speed)
  
  n.pred <- length(object$samples)
  n.sig <- sum(object$sig)
  
  cat("\nCall:\n",object$call)
  cat("\n\n")
  cat("Speed of spread estimated at",n.pred,"locations.\n","Of these,",n.sig,"are significantly nonzero.\n")
  cat("\n\n")
  if(!object$conv.factor==1){
    cat("Mean speed of spread is estimated as",mean.speed,"km.\n","Median speed of spread is estimated as",median.speed,"km.\n")
  }
  if (show==TRUE){
    x11()
    par(mfrow=c(4,1))
    plot(object$post.samp[,1], type = "l", ylab="beta value", main="beta 1")
    plot(object$post.samp[,2], type = "l", ylab="beta value", main="beta 2")
    plot(object$post.samp[,3], type = "l", ylab="beta value", main="beta 3")
    plot(object$post.samp[,4], type = "l", ylab="beta value", main="beta 4")
    
    x11()
    par(mfrow=c(3,1))
    plot(object$post.samp[,4], type = "l", ylab="sigma value", main="Sigma Squared")
    plot(object$post.samp[,5], type = "l", ylab="tau value", main="Tau Squared")
    plot(object$post.samp[,6], type = "l", ylab="phi value", main="Phi")
  }
}

##############################################################################################
conversion <- function(albers,coord){
  ##############################################################################################
  dmat1 <- sqrt(outer(albers[,1], albers[,1], "-")^2 + outer(albers[,2], albers[,2], "-")^2)
  dmat2 <- matrix(0,dim(coord)[1],dim(coord)[1])
  for(i in 1:dim(coord)[1]){
    for(j in 1:i){
      dmat2[i,j] <- getDistanceFromLatLonInKm(coord[i,2],coord[i,1],coord[j,2],coord[j,1])
    }
  }
  dmat2 <- t(dmat2) + dmat2
  result <- mean(dmat2/dmat1,na.rm=T)
  return(result)
}


##############################################################################################
getDistanceFromLatLonInKm = function(lat1,lon1,lat2,lon2) {
  ##############################################################################################
  R = 6371;                      # Radius of the earth in km
  dLat = (lat2-lat1)* (pi/180);  # deg2rad below
  dLon = (lon2-lon1)* (pi/180);
  a =
    sin(dLat/2) * sin(dLat/2) +
    cos((lat1)* (pi/180)) * cos((lat2)* (pi/180)) *
    sin(dLon/2) * sin(dLon/2)
  
  c = 2 * atan2(sqrt(a), sqrt(1-a));
  d = R * c; # Distance in km
  return(d)
}

##############################################################################################
plotgrad = function(object, map=NULL,database=NULL,merge.data=NULL,
                    xlab=expression(paste(degree,"Longitude")),ylab=expression(paste(degree,"Latitude")),xlim=NA, ylim=NA,...){
  ##############################################################################################
  dates <- object$dates
  means <- object$means
  sig <- object$sig
  plotcoord_dat <- subset(merge.data, select = c("longitude","latitude","confirmed"))
  plotcoord <- subset(plotcoord_dat, select = c("longitude","latitude"))
  plotmeans <- means
  
  if(is.na(xlim[1])) { xlim=range(plotcoord[,1]) }
  if(is.na(ylim[1])) { ylim=range(plotcoord[,2]) }
  
  vect = c(seq(min(dates),max(dates),length=12)[2:11],max(dates)+1)
  colramp = heat.colors(length(vect))
  plotcol = colramp[findInterval(dates,vect,all.inside=TRUE)]

  r = sqrt(plotmeans[,1]^2 + plotmeans[,2]^2)
  theta = atan2(plotmeans[,2],plotmeans[,1])
  scale = as.numeric(scale(plotcoord_dat$confirmed))
  # r.scale = r/max(r) * scale
  r.scale = r/max(r) * 0.01

  plot_dat <- as.data.frame(cbind(plotcoord, findInterval(dates,vect,all.inside=TRUE)))
  x11()
  ggplot(data = database,add=TRUE,
         mapping = aes(x = long,
                       y = lat,
                       group = group)) +
    geom_polygon(fill = 'lightgrey',
                 color = 'black') +
    geom_point(data = plot_dat,
               mapping = aes(x = plot_dat[,1],
                             y = plot_dat[,2],
                             group = plot_dat[,3])) +
    scale_color_gradientn(colours = rev(colramp[findInterval(dates,vect,all.inside=TRUE)[order(findInterval(dates,vect,all.inside=TRUE))]]), name = "First date") +
    geom_segment(data = plot_dat,
                 aes(x = plot_dat[,1], y = plot_dat[,2], xend = plot_dat[,1] + r.scale*cos(theta), yend = plot_dat[,2] + r.scale*sin(theta),group = plot_dat[,3]), 
                 arrow=arrow(length=unit(0.3,'cm')), color=rev(colramp[findInterval(dates,vect,all.inside=TRUE)[order(findInterval(dates,vect,all.inside=TRUE))]]),
                 size=1)
}

##############################################################################################
plotgrad_plus = function(object, map=NULL,database=NULL,merge.data=NULL,
                    xlab=expression(paste(degree,"Longitude")),ylab=expression(paste(degree,"Latitude")),xlim=NA, ylim=NA,...){
  ##############################################################################################
  dates <- object$dates
  means <- object$means
  sig <- object$sig
  plotcoord_dat <- subset(merge.data, select = c("longitude","latitude","confirmed"))
  plotcoord <- subset(plotcoord_dat, select = c("longitude","latitude"))
  plotmeans <- means
  
  if(is.na(xlim[1])) { xlim=range(plotcoord[,1]) }
  if(is.na(ylim[1])) { ylim=range(plotcoord[,2]) }
  
  vect = c(seq(min(dates),max(dates),length=12)[2:11],max(dates)+1)
  colramp = heat.colors(length(vect))
  plotcol = colramp[findInterval(dates,vect,all.inside=TRUE)]
  
  r = sqrt(plotmeans[,1]^2 + plotmeans[,2]^2)
  theta = atan2(plotmeans[,2],plotmeans[,1])
  scale = as.numeric(scale(plotcoord_dat$confirmed))/50
  # r.scale = r/max(r)*scale
  r.scale = r/max(r)*0.01
  
  plot_dat <- as.data.frame(cbind(plotcoord, findInterval(dates,vect,all.inside=TRUE)))
  if(!is.null(map)){ 
    x11()
    map +
      geom_point(data = plot_dat,
                 mapping = aes(x = plot_dat[,1],
                               y = plot_dat[,2],
                               group = plot_dat[,3])) +
      # scale_color_gradientn(colours = rev(colramp[findInterval(dates,vect,all.inside=TRUE)[order(findInterval(dates,vect,all.inside=TRUE))]]), name = "First date") +
      geom_segment(data = plot_dat,
                   aes(x = plot_dat[,1], y = plot_dat[,2], xend = plot_dat[,1] + r.scale*cos(theta), yend = plot_dat[,2] + r.scale*sin(theta),group = plot_dat[,3]), 
                   arrow=arrow(length=unit(0.3,'cm')), color=rev(colramp[findInterval(dates,vect,all.inside=TRUE)[order(findInterval(dates,vect,all.inside=TRUE))]]),
                   size=1)  
  }
}


##############################################################################################
plotgrad2 = function(object, map=NULL,database=NULL,merge.data=NULL,
                    xlab=expression(paste(degree,"Longitude")),ylab=expression(paste(degree,"Latitude")),xlim=NA, ylim=NA,...){
  ##############################################################################################
  dates <- data.kor$date
  means <- object$means
  sig <- object$sig
  plotcoord_dat <- subset(merge.data, select = c("longitude","latitude","confirmed"))
  plotcoord <- subset(plotcoord_dat, select = c("longitude","latitude"))
  plotmeans <- means
  
  if(is.na(xlim[1])) { xlim=range(plotcoord[,1]) }
  if(is.na(ylim[1])) { ylim=range(plotcoord[,2]) }
  
  vect = c(seq(min(dates),max(dates),length=12)[2:11],max(dates)+1)
  colramp = heat.colors(length(vect))
  plotcol = colramp[findInterval(dates,vect,all.inside=TRUE)]
  
  r = sqrt(plotmeans[,1]^2 + plotmeans[,2]^2)
  theta = atan2(plotmeans[,2],plotmeans[,1])
  scale = as.numeric(scale(plotcoord_dat$confirmed))
  # r.scale = r/max(r) * scale
  r.scale = r/max(r) * 0.01
  
  plot_dat <- as.data.frame(cbind(plotcoord, findInterval(dates,vect,all.inside=TRUE)))
  x11()
  ggplot(data = database,add=TRUE,
         mapping = aes(x = long,
                       y = lat,
                       group = group)) +
    geom_polygon(fill = 'lightgrey',
                 color = 'black') +
    geom_point(data = plot_dat,
               mapping = aes(x = plot_dat[,1],
                             y = plot_dat[,2],
                             group = plot_dat[,3])) +
    scale_color_gradientn(colours = rev(colramp[findInterval(dates,vect,all.inside=TRUE)[order(findInterval(dates,vect,all.inside=TRUE))]]), name = "First date") +
    geom_segment(data = plot_dat,
                 aes(x = plot_dat[,1], y = plot_dat[,2], xend = plot_dat[,1] + r.scale*cos(theta), yend = plot_dat[,2] + r.scale*sin(theta),group = plot_dat[,3]), 
                 arrow=arrow(length=unit(0.3,'cm')), color=rev(colramp[findInterval(dates,vect,all.inside=TRUE)[order(findInterval(dates,vect,all.inside=TRUE))]]),
                 size=1)
}