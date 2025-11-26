
  map.logbook.locations.by.subarea = function(p, basedir, years=NULL, map.method="lattice"  ) {

    x = logbook.db( DS="logbook" )
    x = x[polygon_inside(x, region="isobath1000m"),]

    if (is.null(years)) years = sort( unique( x$yr ) )

    x = x[, c("yr", "lon", "lat")]
    x = x[ is.finite( rowSums(x) ) ,]
  

    if (map.method=="lattice" ) {
      dir.create (basedir, showWarnings=FALSE, recursive =TRUE)
    
      for (y in years) {

        ii =  which(x$yr==y)
        if ( length(ii)  < 10 ) next()
        toplot = x[ ii, c("lon", "lat")]
        fn = file.path( basedir, paste("logbook.locations", y, "png", sep=".") )
        print(fn)
        png( filename=fn, width=3072, height=2304, pointsize=40, res=300 )
        lp = aegis_map( toplot, depthcontours=TRUE, annot=y, annot.cex=2.8, corners=p$corners, plotlines="cfa.regions", pt.cex=1.5   )
        print(lp)
        dev.off()
      }
    }
  }
