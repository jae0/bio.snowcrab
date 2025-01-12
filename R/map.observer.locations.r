

  map.observer.locations = function(p, basedir,  map.method="lattice", years=NULL ) {

    odb = observer.db( DS="odb")
    odb$yr = odb$fishyr  # use fishyr and not the real year ###################

    if (is.null(years)) years = sort( unique( odb$yr ) )

    odb = odb[, c("yr", "lon", "lat")]
    odb = odb[ is.finite( rowSums(odb) ) ,]

    if (map.method=="lattice" ) {

      corners = data.frame(rbind( cbind( plon=c(220, 990), plat=c(4750, 5270) )))
      for (y in years) {
        ii =  which(odb$yr==y)
        if ( length(ii)  < 10 ) next()
        toplot = odb[ ii, c("lon", "lat")]
        # annot = paste ("Observer", y)
        annot =  y
        outfn = paste("observer.locations", y, sep=".")
        print(outfn)
        dir.create (basedir, showWarnings=FALSE, recursive =TRUE)
        fn = file.path( basedir, paste(outfn, "png", sep="." ) )
        png( filename=fn, width=3072, height=2304, pointsize=40, res=300 )
        lp = aegis_map( xyz=toplot,  depthcontours=TRUE, annot=annot, annot.cex=2.8, corners=corners, plotlines="cfa.regions", pt.cex=1.5 )
        print(lp)
        dev.off()
      }
    }
    return ("Done" )
  }
