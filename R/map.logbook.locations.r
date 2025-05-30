
  map.logbook.locations = function(p, basedir, years=NULL, map.method="lattice"  ) {

    x = logbook.db( DS="logbook" )
    x = x[polygon_inside(x, region="isobath1000m"),]

    if (is.null(years)) years = sort( unique( x$yr ) )

    x = x[, c("yr", "lon", "lat")]
    x = x[ is.finite( rowSums(x) ) ,]

    if (map.method=="lattice" ) {
      for (y in years) {
        ii =  which(x$yr==y)
        if ( length(ii)  < 10 ) next()
        toplot = x[ ii, c("lon", "lat")]
        # annot = paste ("Logbook", y)
        annot = y
        outfn = paste("logbook.locations", y, sep=".")
        print(outfn)
        dir.create (basedir, showWarnings=FALSE, recursive =TRUE)
        fn = file.path( basedir, paste(outfn, "png", sep="." ) )
        png( filename=fn, width=3072, height=2304, pointsize=40, res=300 )
        lp = aegis_map( toplot, depthcontours=TRUE, annot=annot, annot.cex=2.8, corners=p$corners, plotlines="cfa.regions", pt.cex=1.5   )
        print(lp)
        dev.off()
      }
    }
  }
