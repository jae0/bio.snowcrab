
map.set.information = function(
  p, 
  outdir, 
  variables, 
  mapyears,
  plot_method = "aegis_map",
  plot_crs=st_crs( projection_proj4string("utm20N") ) ,  
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



  if (0) {
    plot_crs=st_crs( projection_proj4string("utm20N") )
    theta=p$pres*25
    ptheta=theta/2.3
    probs=c(0.025, 0.975)
    positive_only=TRUE
    log.variable=TRUE
    predlocs = NULL
    mapyears = y = 2025
  }



  # colors = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
  cols = colorRampPalette(c("darkblue","cyan","green", "yellow", "orange","darkred", "black"), space = "Lab")

  if (plot_method == "ggplot") {
    stop("not working yet .. colour map is off")
    require(ggplot2)
    bb = point_to_bbox( p$corners, plot_crs=plot_crs )
    additional_features = snowcrab_mapping_features(p, redo=FALSE )
  }

  set = snowcrab.db( p=p, DS="set.biologicals")

  if(missing(variables)){
    variables = bio.snowcrab::snowcrab.variablelist("all.data")
    variables = intersect( variables, names(set) )
  }

  # define compact list of variable year combinations for parallel processing
  if (missing(mapyears)) mapyears = sort( unique(set$yr) )
  if (exists( "libs", p)) RLibrary( p$libs )

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
        # seems suspect .. do not use
        nplon = length( seq(min(p$corners$plon), max(p$corners$plon), by = p$pres) )
        nplat = length( seq(min(p$corners$plat), max(p$corners$plat), by = p$pres) )
        u= MBA::mba.surf(x=set_xyz[,.(plon,plat,z)], nplon, nplat, sp=TRUE, extend=TRUE   )
        pred_xyz = cbind( predlocs[ips, 1:2], u$xyz.est@data$z[ips] )
      }
    
      if(interpolate.method=='tps'){
        u= fastTps(x=set_xyz[,.(plon,plat)] , Y=set_xyz[["z"]], theta=theta )
        pred_xyz = cbind( predlocs[ips, 1:2], predict(u, xnew=predlocs[ips,1:2]))
      }

      if(interpolate.method=='idw'){
        # broken?
        require(gstat)
        u = gstat(id = "z", formula = z ~ 1, locations = ~ plon + plat, data = set_xyz, set = list(idp = idp))
        pred_xyz = predict(u, predlocs[ips, 1:2])[,1:3]
      }
      
      setDT(pred_xyz)
      names( pred_xyz) = c("plon", "plat", "z")

      pred_xyz$z[ pred_xyz$z > ler[2] ] = ler[2]
      if (ratio) pred_xyz$z[ pred_xyz$z < ler[1] ] = ler[1]

      ckey=NULL
      if (log.variable){
        # create labels for legend on the real scale
        labs = as.vector(c(1,2,5) %o% 10^(-4:5))
        labs = labs[ which( labs > er[1] & labs < er[2] ) ]
        ckey = list( 
          labels = list( at = log10(labs + offset), labels = labs, cex = 2 ) 
        )
      }

      dir.create (outloc, showWarnings=FALSE, recursive =TRUE)

      filename=file.path(outloc, paste(outfn, "png", sep="."))
      print(filename)
          
      if (plot_method == "aegis_map") {
        png( filename=filename, width=3072, height=2304, pointsize=40, res=300 )
        lp = aegis_map( pred_xyz, xyz.coords="planar", depthcontours=TRUE, 
          pts=set_xyz[,.(plon,plat)],
          annot=y, annot.cex=4, at=datarange, col.regions=cols(length(datarange)+1),
          colpts=FALSE, corners=p$corners, display=FALSE, colorkey=ckey, plotlines="cfa.regions" )
        print(lp)
        dev.off()
      }

      if (plot_method == "ggplot") {
        
        Z = st_as_sf( pred_xyz, coords= c("plon", "plat") )
        st_crs(Z) = st_crs( projection_proj4string("utm20") )    # note in km
        Z = stars::st_rasterize( Z["z"], dx=p$pres, dy=p$pres )
        Z = st_transform(Z, plot_crs )  #now im m

        o = ggplot() +
          geom_stars(data=Z, aes(fill=z) ) +  ## <<-- issue is here or scale_fill*
          scale_fill_gradientn(
            name = y,
            limits=range(datarange),
            colours=cols(length(datarange)+1),
            #values = datarange,
            na.value=NA ) +
          coord_sf(xlim = bb$x, ylim =bb$y, expand = FALSE, crs=plot_crs ) +
          guides(
            fill = guide_colorbar(
              title.position = "bottom",
              # title.theme = element_blank(),
              # title.theme = element_text(size = 20),
              label.theme = element_text(size = 14) )
          ) + 
          additional_features +
          theme(
            axis.line=element_blank(),
            # axis.text.x=element_blank(),
            # axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            legend.position = "inside",
            legend.position.inside=c( 0.925, 0.15 ),
            legend.title = element_text( paste0("\n", y), size=20, vjust = -2 ) ,
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

    }
  }

  return("Done")
}
