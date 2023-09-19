
identify_modes = function( Z, sigdigits=2, lowpassfilter=1, lowpassfilter2=0, n_min=30, n_cutoff=0, grad_method="simple",
    plot_solutions=FALSE, bw=1, kernel="gaussian", kexp=1) {
    # find modes and troughs from data vector of histograms
    
    # kernel density based approach
    require(numDeriv)
    require(sf)

    # simple determination from histogram
    Z = Z[ which(is.finite(Z) ) ]
    Zt = sum(Z)
    v = data.table(Z=Z)
    vt = v[, .N, by=Z]
  
    vt$frac = vt$N / sum(vt$N)
    todrop = vt[ N <= n_cutoff, Z] # drop low values
    if (length(todrop) > 0) {
        tdp = which(Z %in% todrop)
        if (length(tdp) > 0 ) Z = Z[- tdp ]
    }
    if (length(Z) <  n_min) return(NA)
 
    u = try( density(Z, bw=bw, kernel=kernel) , silent = TRUE )
    if ("try-error" %in% class(u) )  return(NA)

    u$y0 = u$y
    if (kexp != 1) u$y = u$y ^ kexp  # kexp > 1 augments peaks ..
    
    u$y = u$y / sum(u$y) # normalize
    u$y = zapsmall(u$y)
    u$y = u$y - lowpassfilter
    u$y[u$y < 0] = 0 
 
    dZ  = smooth_data(u$x, numDeriv::grad( approxfun( u$x, u$y, rule=2 ), u$x , method=grad_method ) )
    dZ = dZ / sum( abs(dZ)) # normalize to make scale a probability
    dZ = zapsmall(dZ)

    ddZ = smooth_data(u$x, numDeriv::grad( approxfun( u$x,  dZ, rule=2 ), u$x , method=grad_method ) )
    ddZ = ddZ / sum( abs(ddZ))  
    ddZ = zapsmall(ddZ)
    #ID crossing 0 line
    zero = st_linestring( cbind( range(u$x) *c(0.9, 1.1) , rep(0, 2 )) )

    dZ0 =  as.matrix( st_intersection( st_linestring( cbind( u$x, dZ ) ), zero)  )[,1]
    
    # find closest match
    oo = abs( outer( u$x, dZ0, FUN="-" ) )
    xi = sort( unique( apply(oo, 2, which.min) ) )
    
    ip = xi[ which( ddZ[xi] < -lowpassfilter2 ) ]  # index of peaks
    pky = u$y0[ ip ] 
    pks = u$x[ ip ] 

    trs = u$x[ xi[ which( ddZ[xi] >  lowpassfilter2 ) ] ]

    if (length(pks) > 0) pks = sort( unique( round( pks, sigdigits) ))
    if (length(trs) > 0) trs = sort( unique( round( trs, sigdigits) ))

    if (plot_solutions) {
        plot(u)
        abline( v=pks, col="blue", lty="dotted" )
        abline( v=trs, col="red", lty ="dotted" )
    }
 
    return(list(peaks=pks, troughs=trs, pkvalue=pks, N=Zt))
}


