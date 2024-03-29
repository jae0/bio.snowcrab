 map.survey.locations = function(p, basedir, newyear=T, map.method="lattice" ) {

    set = snowcrab.db( DS="set.clean")
    years = sort( unique( set$yr ) )
    if (newyear) years = p$year.assessment

    if (map.method=="lattice" ) {

      set = set[, c("yr", "plon", "plat")]
      set = set[ is.finite( rowSums(set) ) ,]
      corners = data.frame(rbind( cbind( plon=c(220, 990), plat=c(4750, 5270) )))

      for (y in years) {
        toplot = set[ which(set$yr==y), c("plon", "plat")]
        annot = paste (y)
        fn = file.path(basedir, paste( "survey.locations", y, "png", sep="." ) ) 
        print(fn)
        dir.create (basedir, showWarnings=FALSE, recursive =TRUE)
        png( filename=fn, width=3072, height=2304, pointsize=40, res=300 )
        print(
          aegis_map( xyz=toplot, depthcontours=TRUE, annot=annot, annot.cex=2.8, corners=corners, plotlines="cfa.regions", pt.cex=1.5 )
        )
        dev.off()
      }
    }

    return ("Done" )
  }
