
  plot.freq.distribution.prior.posterior = function( prior, posterior, legendadd=FALSE, br=75, xrange ) {
 
    qnt =  c(0.005, 0.995) # quants
    addr = c(0.9, 1.1)  # additional range

    # prior = NULL; prior$scale=1; prior$shape=1; prior$class="gamma"
    if ( prior$class=="none" ) {
      if (missing(xrange)) xrange = range( posterior, prior$max, prior$min, na.rm=T )
      xval = seq( xrange[1], xrange[2], length=500 )
      postdat = hist( posterior, breaks=br, plot=FALSE )
      yrange = range( 0, postdat$density, na.rm=T ) * 1.02
      hist( posterior, freq=FALSE, breaks=br, xlim=xrange, ylim=yrange, main="", xlab="", ylab="Density", col="lightgray", border="gray")  
      if (legendadd) legend( "topright", bty="n", legend=paste( c( "Prior  ")), lwd=c(1), col=c( "red" ), lty=c("solid") )
      return()
    }

    if ( prior$class == "uniform" ) {
      if (missing(xrange)) xrange = range( posterior, prior$max, prior$min, na.rm=T )
      xval = seq( xrange[1], xrange[2], length=500 )
      dprior = dunif( xval, min=prior$min, max=prior$max )
    }

    if ( prior$class == "normal" ) {
      if (missing(xrange)) xrange = quantile( posterior, qnt, na.rm=T )* addr
      xval = seq( xrange[1], xrange[2], length=500 )
      dprior = dnorm( xval, mean=prior$mean, sd=prior$sd )
    }

    if ( prior$class == "lognormal" ) {
      if (missing(xrange)) xrange = quantile( posterior, qnt, na.rm=T ) * addr
      xrange[1] = min(posterior)
      xval =  exp( seq( log(xrange[1]), log(xrange[2]), length=500 ) )
      xrange[1] = 0
      dprior =  dlnorm( xval, meanlog=prior$meanlog, sdlog=prior$sdlog ) 
    }
   
    if ( prior$class == "beta" ) {
      if (missing(xrange)) xrange = c( -0.1, 1.1 )
      xval = seq( 0, 1, length=500 )
      dprior = dbeta( xval, prior$a, prior$b )
    }
  
    if ( prior$class == "cauchy" ) {
      if (missing(xrange)) xrange = c( -0.1, 1.1 )
      xval = seq( xrange[1], xrange[2], length=500 )
      dprior = dcauchy( xval, prior$a, prior$b )
    }

    if ( prior$class == "gamma" ) {
      if (missing(xrange)) xrange = range( posterior, na.rm=T ) * addr
      xval = seq( xrange[1], xrange[2], length=500 )
      if ( exists( "rate", prior ) ) { prior$scale = 1/prior$rate }
      dprior = dgamma( xval, shape=prior$shape, scale=prior$scale )
    }
  
    postdat = hist( posterior, breaks=br, plot=FALSE )
    #dprior = dprior / sum( dprior) * sum(postdat$density)
    yrange = range( 0, dprior, postdat$density, na.rm=T ) * 1.02
    
    hist( posterior, freq=FALSE, breaks=br, xlim=xrange, ylim=yrange, main="", xlab="", ylab="Density", col="lightgray", border="gray")  
    lines( xval, dprior, col="red", lwd=1, lty="solid" )
    if (legendadd) legend( "topright", bty="n", legend=paste( c( "Prior  ")), lwd=c(1), col=c( "red" ), lty=c("solid") )
    
  }



