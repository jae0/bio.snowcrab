
identify_modes = function( Z, eps=0, lower_filter=0,  n_min=30, plot_solutions=FALSE, bw=1, kernel="gaussian") {
    # find modes and troughs from data vector of histograms
    
    # kernel density based approach
    require(numDeriv)
    require(sf)

    # simple determination from histogram
    Z = Z[ which(is.finite(Z) ) ]
    v = data.table(Z=Z)
    vt = v[, .N, by=Z]
    todrop = vt[N <= lower_filter, Z] # drop low values
    if (length(todrop) > 0) {
        tdp = which(Z %in% todrop)
        if (length(tdp) > 0 ) Z = Z[- tdp ]
    }
    if (length(Z) <  n_min) return(NA)
 
    u = try( density(Z, bw=bw, kernel=kernel) , silent = TRUE )
    if ("try-error" %in% class(u) )  return(NA)

    u$y = u$y / sum(u$y) # normalize
    u$y = zapsmall(u$y)
 
    dZ  = smooth_data(u$x, numDeriv::grad( approxfun( u$x, u$y, rule=2 ), u$x  ) )
    dZ = zapsmall(dZ)
    dZ = dZ / sum( abs(dZ)) # normalize to make scale a probability

    ddZ = smooth_data(u$x, numDeriv::grad( approxfun( u$x,  dZ, rule=2 ), u$x  ) )
    ddZ = zapsmall(ddZ)
    ddZ = ddZ / sum( abs(ddZ)) # normalize to make eps be a probability

    zero = st_linestring( cbind( range(u$x), rep(0, 2 )) )

    dZ0 =  as.matrix( st_intersection( st_linestring( cbind( u$x, dZ ) ), zero)  )[,1]
    
    # find closest match
    oo = abs( outer( u$x, dZ0, FUN="-" ) )
    xi = sort( unique( apply(oo, 2, which.min) ) )
 
    pks = u$x[ xi[ which( ddZ[xi] < -eps ) ] ]
    trs = u$x[ xi[ which( ddZ[xi] >  eps ) ] ]

    if (length(pks) > 0) pks = sort( unique( round( pks, 3) ))
    if (length(trs) > 0) trs = sort( unique( round( trs, 3) ))

    if (plot_solutions) {
        plot(u)
        abline( v=pks, col="blue" )
        abline( v=trs, col="red" )
    }
 
    return(list(peaks=pks, troughs=trs))
}


