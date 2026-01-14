 
#TODO BC add functionality for pdf&kml outputs

map.set.information = function(
  p, 
  outdir, 
  variables, 
  mapyears,
  plotmethod = "ggplot",
  plot_crs=st_crs( projection_proj4string("utm20N") ) ,  # not used yet
  interpolate.method='tps', 
  theta=p$pres*25, 
  ptheta=theta/2.3,
  idp=2, 
  log.variable=TRUE, 
  predlocs=NULL, 
  positive_only=TRUE,
  minN=10, 
  probs=c(0.025, 0.975) 
) {

  x = snowcrab.db( p=p, DS="set.biologicals")
  x = st_as_sf( x, coords= c("plon", "plat"), crs=st_crs( projection_proj4string("utm20") ) )  # note in km

  if (missing(variables)) {
    variables = bio.snowcrab::snowcrab.variablelist("all.data")
    variables = intersect( variables, names(x) )
  }

  # define compact list of variable year combinations for parallel processing
  if (missing(mapyears)) mapyears = sort( unique(x$yr) )

  if (is.null(predlocs)) {
    # identify locations of interest based upon depth and position
    o = get_predlocs(p)
    predlocs = o[["predlocs"]]  
    aoi = o[["aoi"]]
    o = NULL
  } 

  for ( v in variables ) {

    ratio=FALSE
    if (grepl('ratio', v)) ratio=TRUE

    sv = x[[v]]
    sv0 = sv[ which( sv>0 ) ]
      
    # data_offset for log transformation (~ magnitude of observation error) 
    data_offset = quantile( sv0, probs=0, na.rm=TRUE ) 

    # range of all years
    ler =  quantile( sv0, probs=probs, na.rm=TRUE ) 

    for ( y in mapyears ) {
      
      outfn = paste( v,y, sep=".")
      outloc = file.path( outdir, v ) 

      vns = c("plon","plat",v)
      xyz = x[ yr==y, ..vns ]
      setnames( xyz, v, "z" )
      xyz = xyz[ , .( z=mean(z, na.rm=TRUE) ), by=.(plon, plat) ]
      xyz = xyz[is.finite(z),]

      if (nrow(xyz) < minN) next() #skip to next variable if not enough data

      withdata = 1:nrow(xyz)

      if (ratio) {
        ler = c(0,1)
        withdata = which(is.finite( xyz$z ))
      } else {
        if (positive_only) {
          withdata = which(xyz$z > 0)
        }
      } 
      if (length(withdata) < 3) {
        print( paste( "skipped",v, y, "<3 data points to create map", sep="." ) )
        next()
      }
  
      distances =  rdist( 
        predlocs[ aoi, .(plon, plat)], 
        xyz[ withdata, .(plon, plat) ]
      )

      distances[ which(distances < ptheta) ] = NA
      shortrange = which( !is.finite( rowSums(distances) ) )
      ips = aoi[ shortrange ]
  
      if (log.variable){
        xyz$z = log10( xyz$z + data_offset )
        ler = log10( ler + data_offset ) 
      }

      datarange = seq( ler[1], ler[2], length.out=50)
      xyz = na.omit(xyz)

      if (nrow(xyz) < minN || is.na(ler[1])) next() #skip to next variable if not enough data


      if(interpolate.method=='mba'){
        # seems suspect .. deprecated
        nplon = length( seq(min(p$corners$plon), max(p$corners$plon), by = p$pres) )
        nplat = length( seq(min(p$corners$plat), max(p$corners$plat), by = p$pres) )

        u= MBA::mba.surf(x=xyz[,.(plon,plat,z)], nplon, nplat, sp=TRUE, extend=TRUE   )
        pred_xyz = cbind( predlocs[ips, 1:2], u$xyz.est@data$z[ips] )
      }
    
      if(interpolate.method=='tps'){
        u= fastTps(x=xyz[,.(plon,plat)] , Y=xyz[["z"]], theta=theta )
        pred_xyz = cbind( predlocs[ips, 1:2], predict(u, xnew=predlocs[ips,1:2]))
      }

      if(interpolate.method=='idw'){
        # broken?
        require(gstat)
        u = gstat(id = "z", formula = z ~ 1, locations = ~ plon + plat, data = xyz, set = list(idp = idp))
        pred_xyz = predict(u, predlocs[ips, 1:2])[,1:3]
      }
      
      setDT( pred_xyz )
      names( pred_xyz ) = c("plon", "plat", "z")
      
      cols = colorRampPalette(c("darkblue","cyan","green", "yellow", "orange","darkred", "black"), space = "Lab")

      pred_xyz$z[ pred_xyz$z > ler[2] ] = ler[2]
      pred_xyz$z[ pred_xyz$z < ler[1] ] = ler[1]

      col_key = NULL
      if (log.variable){
        # create labels for legend on the real scale
        labs = as.vector( c(1,2,5) %o% 10^(-4:5) )
        labs = labs[ which(labs>ler[1] & labs<ler[2]) ]
        col_key = list(
          labels=list( at=log10(labs+data_offset), labels=labs, cex=2) 
        )
      }

      dir.create (outloc, showWarnings=FALSE, recursive =TRUE)

      filename=file.path(outloc, paste(outfn, "png", sep="."))
      print(filename)
      
          
          stop("not working yet")

          
      if (plotmethod=="ggplot") {
        
        require(ggplot2)

        o = ggplot() +
          geom_sf( data = sppoly, aes(fill=.data[["z"]]), lwd=0, alpha=1  )  +
          coord_sf(xlim = bb$x, ylim =bb$y, expand = FALSE, crs=plot_crs ) +
          scale_fill_gradientn(
            name = yrs[i], 
            limits=range(datarange),
            colors=colors, 
            na.value=NA ) +
          guides(
            fill = guide_colorbar(
              title.position = "bottom",
              # title.theme = element_blank(), 
              # title.theme = element_text(size = 20),
              label.theme = element_text(size = 14) ) )  

        if (!is.null(additional_features)) o = o + additional_features

        o = o + 
          theme(
            axis.line=element_blank(),
            # axis.text.x=element_blank(),
            # axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(), 
            legend.position = "inside",
            legend.position.inside=legend.position.inside,
            legend.title = element_text( paste0("\n", yrs[i]), size=20, vjust = -2 ) ,
            # panel.background=element_blank(),
            panel.background = element_rect(fill =NA),
            panel.border=element_blank(),
            # panel.grid.major=element_blank(),
            panel.grid.major = element_line(color = "grey"),
            panel.grid.minor=element_blank(),
            plot.background = element_rect(fill="white"), 
            plot.caption = element_text(hjust = 0, size = 12)
          )

        ggsave( fn, o )
        print(fn)
      }

      if (plotmethod=="aegis_map") {
        png( filename=filename, width=3072, height=2304, pointsize=40, res=300 )
          lp = aegis_map( pred_xyz, xyz.coords="planar", depthcontours=TRUE, 
            pts=pred_xyz[,.(plon,plat)],
            annot=y, annot.cex=4, at=datarange, col.regions=cols(length(datarange)+1),
            colpts=FALSE, corners=p$corners, display=FALSE, colorkey=col_key, plotlines="cfa.regions" )
          print(lp)
        dev.off()
      }

    } # end for y
  } # end for v
  
  return("Done")
}
