timeseries_simple = function( dat, regions, yrs, vn, lookup.table=NULL, sdci=FALSE ) {

    setDT(dat)
    
    tsdata = CJ( region=regions, year=as.character(yrs), variable=vn )
    tsdata[, uid := paste(region, year, variable, sep="_") ]

    tsdata$mean = NA
    tsdata$se = NA
    tsdata$sd = NA
    tsdata$n = NA
    tsdata$ub = NA
    tsdata$lb = NA
 
    nv = length(vn)

    for (vi in 1:nv ) {
      v = vn[vi]
      if ( !is.numeric( dat[[v]] ) ) next()
      if (v %in% c("uid", "year" ) ) next()

      message( vi, " of ", nv, " -- ", v ) 

      if (!is.null(lookup.table)) {
        dat$XX = 0
        dat$XX = bio.snowcrab::variable.recode(
            dat[[v]],
            v,
            direction="forward",
            lookup.table=lookup.table
        ) # transform variables where necessary
      } else {
        dat$XX = dat[[v]]
      }


      for (r in regions) {
        ri = which( dat[[r]] == 1 )
        if (length(ri)==0) next()
  
        res = dat[ri, .(
          mean = mean( XX, na.rm=TRUE), 
          n = length(which(is.finite(XX))),
          sd = sd(XX, na.rm=TRUE)
          ), 
          by=.(year) 
        ]
        res[["se"]] = res[["sd"]] / sqrt(res[["n"]] - 1)
        
        if(sdci) {
          res[["lb"]] = res[["mean"]] - res[["sd"]]
          res[["ub"]] = res[["mean"]] + res[["sd"]]
        } else {
          res[["lb"]] = res[["mean"]] - res[["se"]]
          res[["ub"]] = res[["mean"]] + res[["se"]]
        }

        res[["region"]] = r
        res[["variable"]] = v
        res[["uid"]] = paste( r, res[["year"]], v, sep="_" )  

        tsi = match( res$uid, tsdata$uid )

        if (!is.null(lookup.table)) {
            res[["mean"]] = bio.snowcrab::variable.recode ( res[["mean"]], v, direction="backward", lookup.table=lookup.table )
            res[["n"]] = res[["n"]]
            res[["se"]] = bio.snowcrab::variable.recode (res[["se"]], v, direction="backward" , lookup.table=lookup.table, is.sd=TRUE)
            res[["sd"]] = bio.snowcrab::variable.recode (res[["sd"]], v, direction="backward" , lookup.table=lookup.table, is.sd=TRUE)
            res[["lb"]] = bio.snowcrab::variable.recode (res[["lb"]], v, direction="backward" , lookup.table=lookup.table)
            res[["ub"]] = bio.snowcrab::variable.recode (res[["ub"]], v, direction="backward" , lookup.table=lookup.table)
        }

        tsdata[["mean"]][ tsi ] =  res[["mean"]]
        tsdata[["n"]][ tsi ] = res[["n"]]
        tsdata[["se"]][ tsi ] = res[["se"]]
        tsdata[["sd"]][ tsi ] = res[["sd"]]
        tsdata[["lb"]][ tsi ] = res[["lb"]]
        tsdata[["ub"]][ tsi ] = res[["ub"]]
        
      }
    }
    
    tsdata$year = as.numeric( tsdata$year)

    return(tsdata)
}
