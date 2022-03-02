
plot_x__y_density = function( x=NULL, y_mean=NULL, y_sd=NULL, ns=1000, dist="gaussian", xlab="", ylab="" ) {

  if ( is.null(x)) {
    # test
    x =1:10
    y_mean =runif(10)
    y_sd = runif(10) / 10
    xlab="xtest"
    ylab="ytest"
    ns = 1000
  }

  if (dist == "gaussian" ) {
    nd = length(x)

    y = matrix( NA, nrow=ns, ncol=nd )
    for (i in 1:nd) y[,i] = rnorm( ns, mean=y_mean[i], sd=y_sd[i])

    prs = seq( from=0.025, to=0.975, length.out=500)
    Bq =  apply( y, 2, quantile, probs=prs, na.rm=T  )

    yran = range( pretty( Bq ) )
    plot( x, Bq[1,], type="n", ylim=yran, xlim=range(x), xlab=xlab, ylab=ylab ) #change xlim to yrs0 to remove 3 yr projection
    cols = gray.colors( floor(length( prs)/2) )
    cols2 = c(cols[length(cols):1], cols )
    for ( j in 1:length(prs) ) lines ( x, Bq[j,], lwd=4, col=cols2[j] )
  }
}
