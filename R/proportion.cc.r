
  proportion.cc = function(odb, region, year) {
    # sex codes
    male = 0
    female = 1
    sex.unknown = 2

    out= NULL

    # Remove CW's outside norms and remove production (pre-sorted) samples
    i = which( odb$sex==male & odb$fishyr==year )  # --- fishing year is used and not the actual year caught
    r = emgis::polygon_inside(x=odb, region=emgis::polygon_internal_code(region), planar=F)

    c.r.y = unique( intersect (r, i) )

    # soft = which( odb$durometer < 68 ) # unreliable

    ccs = c(1:5)
    for (cc in ccs) {
      icc = which( odb$shell == cc )
     # if (cc==1) icc = unique( union( icc, soft ) )
     # if (cc>1)  icc = unique( setdiff( icc, soft ) )
      out = c(out, length( intersect( c.r.y, icc ) ) )
    }

    ntot = sum(out, na.rm=T)

    out = c( round(out/ntot*100,2), ntot)

    return(out)
  }



