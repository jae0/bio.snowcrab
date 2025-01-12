
identify_modes = function( Z, W=NULL, T=NULL, V=NULL, X=NULL, 
    sigdigits=2, lowpassfilter=0, lowpassfilter2=0, dx=NULL, xvals=NULL,
    n_min=30, grad_method="simple", decompose_distributions=FALSE, 
    plot_solutions=FALSE, bw=1, kernel="gaussian", kexp=1, override_range=NULL ) {
    # find modes and troughs from data vector of histograms
    
    require(numDeriv)
    require(sf)
 
    if (0){
        W=NULL
        T=NULL
        V=NULL
        X=NULL        
        n_min=30
        grad_method="simple"
        decompose_distributions=FALSE
        plot_solutions=TRUE
        kernel="gaussian"
        kexp=1
        override_range=NULL
        Z = unlist(MO[ sex=="0" & mat=="1" , cw])
        lowpassfilter2=0.0005
        dx=ldx
        bw=0.02
        sigdigits=3
        plot=TRUE
    }

    if (!is.null(X)) {
        # treat Z and X as density 

        u = data.table(x=X, y=Z )
        class(u) = "density"
        u$bw = bw
        u$n = 1L

    } else {
   
     # kernel density based approach
        # simple determination from histogram
        hasdata = which( is.finite(Z) )
        Z = Z[ hasdata ]
        if (length(Z) <  n_min) return(list(peaks=NA, troughs=NA, peak_values=NA, N=NA))

        if (!is.null(W)) {
            W = W[hasdata]
            hasdata2 = which( is.finite(W) )
            W = W[hasdata2]
            Z = Z[hasdata2]
            W = W/sum(W)
        } else {
            W = rep(1/length(Z), length(Z) ) 
        } 

        u = try( density(Z, bw=bw, kernel=kernel, weights=W), silent = TRUE )
    }
    
    if ("try-error" %in% class(u) )  return(list(peaks=NA, troughs=NA, peak_values=NA, N=NA, u=NA, res=NA))

    if (!is.null(T)) {
        # troughs also given ... multiply with peaks to adjust PDFs
 
        DZ = u
        #DZ$x = trunc( DZ$x/dx ) * dx 
        DZ = data.table(x=DZ$x, Py=DZ$y) 
        #DZ = DZ[, .(Py=mean(Py, na.rm=TRUE)), by=x]

        hasdata = which(is.finite(T) )
        if (!is.null(V)) {
            hasdata = intersect( hasdata, which(is.finite(V) ) )
            V = V[hasdata]
        }
        T = T[ hasdata ]

        if (!is.null(V))  V = V/sum(V)
        DT = density( T, bw=bw, kernel=kernel,  weights=V )
        #DT$x = trunc( DT$x/dx ) * dx 
        DT = data.table(x=DT$x, Ty=DT$y) 
        #DT = DT[, .(Ty=mean(Ty, na.rm=TRUE)), by=x]

        if (!is.null(override_range)) {
            i = which( DT$x >= override_range[1] & DT$x <= override_range[2] )
            if (length(i) > 0) DT$Ty[i] = 0 
        }

        DC = data.table( x = xvals )
#        DC = DZ[ DC, on="x" ]  
#        DC = DT[ DC, on="x" ]  
 
        DZ = unique(DZ)
        DC$Py = smooth_data( DC$x, approxfun( DZ$x, DZ$Py, rule=2 )(DC$x) ) 
        DC$Py = DC$Py / sum( DC$Py) / dx # normalize to density

        DT = unique(DT)
        DC$Ty = smooth_data( DC$x, approxfun( DT$x, DT$Ty, rule=2 )(DC$x) ) 
        DC$Ty = DC$Ty / sum( DC$Ty) / dx # normalize to density

        DC[, e:=Py*(1-Ty)]  
#        DC[ e<0, "e" ]= 0
#        DC$e = zapsmall(DC$e)

        # plot( Py ~ x, DC, type="l", col="gray" )
        # lines( e ~ x, DC, col="red")
    
        u$x = DC$x 
        u$y = DC$e

    }  

    u$y0 = u$y
    if (kexp != 1) u$y = u$y ^ kexp  # kexp > 1 augments peaks ..
    
    u$y = zapsmall(u$y)
    u$y = u$y / dx
    u$y = u$y / sum(u$y)  # normalize
    u$y = u$y - lowpassfilter
    u$y[u$y < 0] = 0 

    dZ  = smooth_data(u$x, numDeriv::grad( approxfun( u$x, u$y, rule=2 ), u$x , method=grad_method ) )
    dZ = zapsmall(dZ)
    dZ = dZ / sum( abs(dZ)) # normalize to make scale a probability

    ddZ = smooth_data(u$x, numDeriv::grad( approxfun( u$x,  dZ, rule=2 ), u$x , method=grad_method ) )
    ddZ = zapsmall(ddZ)
    ddZ = ddZ / sum( abs(ddZ))  
    #ID crossing 0 line
    zero = st_linestring( cbind( range(u$x) *c(0.9, 1.1) , rep(0, 2 )) )

    dZ0 =  as.matrix( st_intersection( st_linestring( cbind( u$x, dZ ) ), zero)  )[,1]
    
    # find closest match
    oo = abs( outer( u$x, dZ0, FUN="-" ) )
    xi = sort( unique( apply(oo, 2, which.min) ) )
    
    ip = xi[ which( ddZ[xi] < -lowpassfilter2 ) ]  # index of peaks
    pks = u$x[ ip ] 
    pky = u$y0[ ip ] 

    it = xi[ which( ddZ[xi] >  lowpassfilter2 ) ]
    trs = u$x[ it ]
    try = u$y0[ it ]

    if (length(pks) > 0) pks = sort( unique( round( pks, sigdigits) ))
    if (length(trs) > 0) trs = sort( unique( round( trs, sigdigits) ))

    if (plot_solutions) {
        plot(u)
        abline( v=pks, col="blue", lty="dotted" )
        abline( v=trs, col="red", lty ="dotted" )
    } 

    res = list()
    
    if (decompose_distributions) {
        # individiual weights exist ... flag to decompose the distributions
        max_classes = 14
        finished = FALSE
        xp = sort(pks, decreasing=TRUE)
        count = length(xp)
        xd = diff(u$x)[1]
        xv = u$x
        yv = u$y*xd
        yv = yv / sum(yv)
        # lines( yv~xv)
        for (h in 1:count) {
            rg =  xp[h] + c(-bw, bw)
            xi = which( xv>xp[h]-3*bw & xv<xp[h]+3*bw) 
            x = sample( xv[ xi ] , length(xv), replace=TRUE, prob=yv[xi]/sum(yv[xi]) )
            n <- length(x)
            xvariance <- crossprod(x - xp[h])[1L] / n  # (a bised estimator, though)
            xsd = sqrt(  xvariance )
            oo = dnorm( xv, mean=xp[h], sd=xsd )
            # not yet resolved
            yv_sum = sum(yv)
##             oo = oo * xd / (yv() oo[which.max(oo)]
            # sum(oo) # 1
            yv = yv - oo  # residual
            xtozero = xp[h] + c(-1, 1) * xsd
            yv[ which(xv > xtozero[1] & xv < xtozero[2]) ] = 0  
            res[[h]] = c( xp[h], xsd, n )
            if (plot_plot_solutionssolutions)  {
                points( oo ~ u$x, col="blue", pch=19 )
                points( yv ~ u$x, col="red", pch=19 )
            }
        }
    }
    
    out = list(peaks=pks, troughs=trs, peak_values=pky, trough_values=try, 
      N=sum(Z), u=u, res=res )

    return(out)
}
 