 
#TODO BC add functionality for pdf&kml outputs

map.set.information = function(p, outdir, variables, mapyears,
  plot_crs=st_crs("EPSG:32620"),  # UTM20N,
  interpolate.method='tps', theta=p$pres*25, ptheta=theta/2.3,
  idp=2, log.variable=TRUE, predlocs=NULL, positive_only=TRUE,
  minN=10, probs=c(0.025, 0.975) ) {

    set = snowcrab.db( p=p, DS="set.biologicals")

    if(missing(variables)){
      variables = bio.snowcrab::snowcrab.variablelist("all.data")
      variables = intersect( variables, names(set) )
    }

    # define compact list of variable year combinations for parallel processing
    if (missing(mapyears)) mapyears = sort( unique(set$yr) )
    if (exists( "libs", p)) RLibrary( p$libs )

    nplon = length( seq(min(p$corners$plon), max(p$corners$plon), by = p$pres) )
    nplat = length( seq(min(p$corners$plat), max(p$corners$plat), by = p$pres) )

    if (is.null(predlocs)) {
      o = get_predlocs(p)
      predlocs = o[["predlocs"]]  
      aoi = o[["aoi"]]
      o = NULL
    } 

    for ( v in variables ) {

      ratio=FALSE
      if (grepl('ratio', v)) ratio=TRUE

      sv = set[[v]]
      sv0 = sv[ which( sv>0 ) ]
      
      # offset for log transformation
      offset = quantile( sv0, probs=0, na.rm=TRUE ) 

      # range of all years
      er =  quantile( sv0, probs=probs, na.rm=TRUE ) 

      for ( y in mapyears ) {
        
          outfn = paste( v,y, sep=".")
          outloc = file.path( outdir,v) 
 
          vns = c("plon","plat",v)
          set_xyz = set[ yr==y, ..vns ]
          setnames( set_xyz, v, "z" )
          set_xyz = set_xyz[ , .( z=mean(z, na.rm=TRUE) ), by=.(plon,plat) ]
          set_xyz = set_xyz[is.finite(z),]

          if(nrow(set_xyz)<minN)next() #skip to next variable if not enough data
 
          if (ratio) {
            er=c(0,1)
            withdata = which(is.finite( set_xyz$z ))
          } else {
            if (positive_only) withdata = which(set_xyz$z > 0)
          }

          ler = er
          
          if (length(withdata) < 3) {
            print(paste("skipped",v, y, "<3 data points to create map", sep="."))
            next()
          }


          distances =  rdist( 
            predlocs[ aoi, .(plon, plat)], 
            set_xyz[ withdata, .(plon, plat) ]
          )

          distances[ which(distances < ptheta) ] = NA
          shortrange = which( !is.finite( rowSums(distances) ) )
          ips = aoi[ shortrange ]
          
          if (log.variable){
            set_xyz$z = log10(set_xyz$z+offset)
            ler=log10(er+offset) 
          }

          datarange = seq( ler[1], ler[2], length.out=50)
          set_xyz = na.omit(set_xyz)

          if (nrow(set_xyz) < minN || is.na(er[1])) next() #skip to next variable if not enough data
  

          if(interpolate.method=='mba'){
            # seems suspect .. 
            u= MBA::mba.surf(x=set_xyz[,.(plon,plat,z)], nplon, nplat, sp=TRUE, extend=TRUE   )
            pred_xyz = cbind( predlocs[ips, 1:2], u$xyz.est@data$z[ips] )
          }
        
          if(interpolate.method=='tps'){
            
            u= fastTps(x=set_xyz[,.(plon,plat)] , Y=set_xyz[["z"]], theta=theta )
            xyz = cbind( predlocs[ips, 1:2], predict(u, xnew=predlocs[ips,1:2]))
          }

          if(interpolate.method=='idw'){
            # broken?
            require(gstat)
            u = gstat(id = "z", formula = z ~ 1, locations = ~ plon + plat, data = set_xyz, set = list(idp = idp))
            xyz = predict(u, predlocs[ips, 1:2])[,1:3]
          }
          
          setDT(xyz)
          names( xyz) = c("plon", "plat", "z")
          
          cols = colorRampPalette(c("darkblue","cyan","green", "yellow", "orange","darkred", "black"), space = "Lab")

          xyz$z[xyz$z>ler[2]] = ler[2]
          if (ratio) xyz$z[xyz$z<ler[1]] = ler[1]

          ckey=NULL
          if (log.variable){
            # create labels for legend on the real scale
            labs=as.vector(c(1,2,5)%o%10^(-4:5))
            labs=labs[which(labs>er[1]&labs<er[2])]
            ckey=list(labels=list(at=log10(labs+offset),labels=labs,cex=2))
          }

          dir.create (outloc, showWarnings=FALSE, recursive =TRUE)

          filename=file.path(outloc, paste(outfn, "png", sep="."))
          print(filename)
          png( filename=filename, width=3072, height=2304, pointsize=40, res=300 )
          lp = aegis_map( xyz, xyz.coords="planar", depthcontours=TRUE, 
            pts=set_xyz[,.(plon,plat)],
            annot=y, annot.cex=4, at=datarange, col.regions=cols(length(datarange)+1),
            colpts=FALSE, corners=p$corners, display=FALSE, colorkey=ckey, plotlines="cfa.regions" )
          print(lp)

          dev.off()

      }
    }

    return("Done")
  }
