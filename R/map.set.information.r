
map.set.information = function(
  p, 
  outdir, 
  variables, 
  mapyears,
  plot_method = "ggplot",
  plot_crs = st_crs( projection_proj4string("utm20N") ) ,  
  interpolate.method = "tps", 
  theta=p$pres*30, 
  ptheta=theta/2,
  idp=2, 
  predlocs=NULL,   
  minN=1, 
  probs=c(0.005, 0.995) 
) { 

  if (0) {
    plot_crs=st_crs( projection_proj4string("utm20N") )
    theta=p$pres*25
    ptheta=theta/2.3
    probs=c(0.025, 0.975)
    log.variable=TRUE
    predlocs = NULL
    mapyears = y = 2025
    plot_method = "ggplot"
    variables = NULL
    minN=3
 
  }

  # colors = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
  cols = colorRampPalette(c("darkblue","cyan","green", "yellow", "orange","darkred", "black"), space = "Lab")

  if (plot_method == "ggplot") {
    
    require(ggplot2)
    bb = point_to_bbox( p$corners, plot_crs=plot_crs )
    additional_features = snowcrab_mapping_features(p, redo=FALSE )
  }

  set = snowcrab.db( p=p, DS="set.biologicals")

  if(missing(variables)){

    variables = bio.snowcrab::snowcrab.variablelist("all.data")
    variables = intersect( variables, names(set) )

    ratio_vars = c("sexratio.all", "sexratio.mat", "sexratio.imm")
    nolog.variables = c("t", "z", "julian", variables[grep("cw", variables)])
    log.variables = variables[ !variables %in% nolog.variables ]

    mass_vars = log.variables[ grep('mass', log.variables)]
    no_vars = log.variables[ grep('no', log.variables)]

  }

  # define compact list of variable year combinations for parallel processing
  if (missing(mapyears)) mapyears = sort( unique(set$yr) )

  if (is.null(predlocs)) {
    o = get_predlocs(p)
    predlocs = o[["predlocs"]]  
    aoi = o[["aoi"]]
    o = NULL
  } 

  for ( v in variables ) {
 
    sv = set[[v]]
    sv0 = sv[ which( sv>0 ) ]

    # range of all years
    if (grepl('ratio', v)) {
      data_range = c(0, 1)
    } else {
      data_range = quantile( sv0, probs=probs, na.rm=TRUE ) 
    }

    if ( v %in% log.variables ){
      # data_offset for log transformation
      data_offset = quantile( sv0, probs=0, na.rm=TRUE ) 
      data_range = log10( data_range + data_offset ) 
      theta = 40
      ptheta = theta/2
    }

    if ( any(is.na(data_range))) {
      message(v, " has no valid data?" )
      next()
    }

    nd = 50

    color_range = seq( data_range[1], data_range[2], length.out=nd )
  
    if (plot_method == "ggplot") {
      # interval for each color value -- ggplot wants the actual intervals
      ggvalues = rep(NA, 2*nd)
      uu = 1:nd *2
      ggvalues[uu-1] = color_range 
      ggvalues[uu] = color_range + min(diff(color_range)) / (nd*10) 
      ggvalues = scales::rescale(ggvalues)
    }
    
    for ( y in mapyears ) {
        
      outfn = paste( v,y, sep=".")
      outloc = file.path( outdir,v) 

      vns = c("plon","plat",v)
      set_xyz = set[ yr==y, ..vns ]
      setnames( set_xyz, v, "z" )
      set_xyz = set_xyz[ , .(z=mean(z, na.rm=TRUE) ), by=.(plon,plat) ]
      set_xyz$z[ !is.finite(set_xyz$z) ] = NA

      set_xyz$plon = jitter(set_xyz$plon)
      set_xyz$plat = jitter(set_xyz$plat)
      set_xyz$z = jitter(set_xyz$z)
      
      if (nrow(set_xyz) < minN) next() #skip to next variable if not enough data

      if (grepl('ratio', v)) {
        withdata = which(is.finite( set_xyz$z ))
      } else {
        withdata = 1:nrow(set_xyz)
      }
      
      if (length(withdata) < 1) {
        message("Skipping", v, y, "... < 1 data points \n")
        next()
      }

      distances =  rdist( 
        predlocs[ aoi, .(plon, plat)], 
        set_xyz[ withdata, .(plon, plat) ]
      )

      distances[ which(distances < ptheta) ] = NA
      shortrange = which( !is.finite( rowSums(distances) ) )
      ips = aoi[ shortrange ]
      
      if ( v %in% log.variables ){
        set_xyz$z = log10(set_xyz$z + data_offset)
      }
      
      set_xyz$z[ !is.finite(set_xyz$z) ] = NA

      if ( nrow(set_xyz) < 1 ) {
        message( "Skipping", v, y, "... unexpected data/Inf \n")
        next() 
      }

      if (interpolate.method=='mba'){
        # seems suspect .. do not use
        nplon = length( seq(min(p$corners$plon), max(p$corners$plon), by = p$pres) )
        nplat = length( seq(min(p$corners$plat), max(p$corners$plat), by = p$pres) )
        u= MBA::mba.surf(x=set_xyz[,.(plon,plat,z)], nplon, nplat, sp=TRUE, extend=TRUE   )
        pred_xyz = cbind( predlocs[ips, 1:2], u$xyz.est@data$z[ips] )
      }
      
      if (interpolate.method=='tps'){
        u = fastTps(x=set_xyz[,.(plon,plat)] , Y=set_xyz[["z"]], theta=theta )
        pred_xyz = cbind( predlocs[ips, 1:2], predict(u, xnew=predlocs[ips,1:2]))
      }

      if (interpolate.method=='idw'){
        # broken?
        require(gstat)
        u = gstat(id = "z", formula = z ~ 1, locations = ~ plon + plat, data = set_xyz, set = list(idp = idp))
        pred_xyz = predict(u, predlocs[ips, 1:2])[,1:3]
      }
      
      setDT(pred_xyz)
      names( pred_xyz) = c("plon", "plat", "z")

      pred_xyz$z[ pred_xyz$z > data_range[2] ] = data_range[2]
      if (grepl('ratio', v)) pred_xyz$z[ pred_xyz$z < data_range[1] ] = data_range[1]

      dir.create (outloc, showWarnings=FALSE, recursive =TRUE)

      filename = file.path(outloc, paste(outfn, "png", sep="."))
      print(filename)
          
      if (plot_method == "aegis_map") {
        ckey=NULL
        if ( v %in% log.variables ){
          # create labels for legend on the real scale
          labs = as.vector(c(1,2,5) %o% 10^(-4:5))
          labs = labs[ which( labs > data_range[1] & labs < data_range[2] ) ]
          ckey = list( 
            labels = list( at = log10(labs + data_offset), labels = labs, cex = 2 ) 
          )
        }

        png( filename=filename, width=3072, height=2304, pointsize=40, res=300 )
        lp = aegis_map( pred_xyz, xyz.coords="planar", depthcontours=TRUE, 
          pts=set_xyz[,.(plon,plat)],
          annot=y, annot.cex=4, at=color_range, col.regions=cols(length(color_range)+1),
          colpts=FALSE, corners=p$corners, display=FALSE, colorkey=ckey, plotlines="cfa.regions" )
        print(lp)
        dev.off()
      }

      if (plot_method == "ggplot") {
        
        # create labels for legend on the real scale

        if ( v %in% log.variables ){
          dr = data_range - log10(data_offset)
          dx = 1
          labs = seq( -10, 10, by=dx) 
          labs = labs[ which( labs >= (dr[1]-dx*1.1) & labs <= (dr[2] + dx*1.1) ) ]
          brks = labs
          labs = signif(10^labs, 1)
        } else {
          labs = pretty(color_range, n=3)
          brks = labs
        }

        Z = st_as_sf( pred_xyz, coords= c("plon", "plat") )
        st_crs(Z) = st_crs( projection_proj4string("utm20") )    # note in km
        Z = stars::st_rasterize( Z["z"], dx=p$pres, dy=p$pres )
        Z = st_transform(Z, plot_crs )  #now im m
        Z = as.data.table(Z)

        set_xyz = st_as_sf( set_xyz, coords= c("plon", "plat") )
        set_xyz$z = NULL

        st_crs(set_xyz) = st_crs( projection_proj4string("utm20") )    # note in km
        set_xyz = st_transform( set_xyz, plot_crs )  #now im m

        if ( v %in% log.variables ){
          Z$z = Z$z - log10( data_offset)
          Z$z[ Z$z < 0 ] = 0
        }

        o = ggplot() +
          geom_raster( data=Z, aes(x=x, y=y, fill=z), alpha=0.99 ) +  ## <<-- issue is here or scale_fill*
          scale_fill_gradientn(
            name = y,
            limits = range(brks),
            colors = cols(length(color_range)),
            values = ggvalues,  # interval for each color
            labels = labs ,
            breaks = brks,
            na.value=NA ) +
          labs(caption = v ) +
          coord_sf( xlim = bb$x, ylim =bb$y, expand = FALSE, crs=plot_crs ) +
          geom_sf( data=set_xyz, size=0.75, alpha=0.25)   +
          guides(
            fill = guide_colorbar(
              title.position = "top",
              label.theme = element_text(size = 10) )
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
            legend.position.inside=c( 0.07, 0.8 ),
            legend.title = element_text( as.character(y), size=16, vjust = 2 ) ,
            # panel.background=element_blank(),
            panel.background = element_rect(fill =NA),
            panel.border=element_blank(),
            # panel.grid.major=element_blank(),
            panel.grid.major = element_line(color = "grey"),
            panel.grid.minor=element_blank(),
            plot.background = element_rect(fill="white"),
            plot.caption = element_text(hjust = 0, size = 10)
          )

        ggsave( filename, o )
 
      }

    }
  }

  return("Done")
}
