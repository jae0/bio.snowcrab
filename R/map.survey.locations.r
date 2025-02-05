 map.survey.locations = function(p, 
  basedir=project.datadirectory("bio.snowcrab", "output", "maps", "survey.locations" ),
  years=NULL, map.method="lattice", set=NULL ) {

    if (is.null(set)) set = snowcrab.db( DS="set.clean")
    if (is.null(years)) years = sort( unique( set$yr ) )
 
    if (map.method=="lattice" ) {

      set = set[, c("yr", "plon", "plat")]
      set = set[ is.finite( rowSums(set) ) ,]
      corners = data.frame(rbind( cbind( plon=c(220, 990), plat=c(4750, 5270) )))

      for (y in years) {
        toplot = set[ which(set$yr==y), c("plon", "plat")]
        if (nrow(toplot) == 0) next()
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
