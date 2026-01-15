map.fisheries.data = function(
  outdir, 
  yrs, 
  variables=c('effort', 'cpue', 'landings') ,
  FUN = c(sum, mean, sum),
  probs=c(0.025, 0.975), 
  pres = 10, 
  plot_crs = st_crs( projection_proj4string("utm20N"))  ,
  outformat = "png",
  plotmethod = "ggplot",
  ...
  ) {

  # colsfn = colorRampPalette(c("darkblue","cyan","green", "yellow", "orange","darkred", "black"), space = "Lab")
  colors_brew = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))

  require(sf)
  
  additional_features = snowcrab_mapping_features(p, redo=FALSE ) 

  x = logbook.db( DS="logbook" )
  x = x [polygon_inside( x, region="isobath1000m"),]
  x = x[ which(x$effort <= 300) ,]
  x = x[ which(x$cpue < 500),]
  x$year = x$yr #this creates proper offset for 4X, 2017-18 season =2017
  x$landings = x$landings/1000 
 
  x = st_as_sf( x, coords= c("lon", "lat"), crs=st_crs( projection_proj4string("lonlat_wgs84") ) )
  x = st_transform(x, st_crs( projection_proj4string("utm20") ) )
  
  if (missing(yrs)) yrs=sort(unique(x$year))
 
  polygrid = st_as_sf( st_make_grid( x, cellsize=pres, what="polygons", square=TRUE ) )
  polygrid$AUID = as.character( 1:nrow(polygrid) )  # row index
  row.names(polygrid) = polygrid$AUID

  x$AUID = st_points_in_polygons( pts=x, polys=polygrid[, "AUID"], varname="AUID" )

  
  ellps = list(...)

  if ( exists("legend.position.inside", ellps)) {
    legend.position.inside = ellps[["legend.position.inside"]] 
  } else {
    legend.position.inside = c( 0.925, 0.15 ) 
  }
 
  polygrid = st_transform(polygrid, plot_crs )
  
  bb = point_to_bbox( 
    data.table( lon=c(-66.1, -56.6), lat = c(42.8,  47.4) ), 
    plot_crs=plot_crs 
  )

  for (v in 1: length(variables)) {
    vn =variables[v]
    er = quantile( x[[vn]], probs=probs, na.rm=TRUE )
    datarange = pretty( er  )  
    
    outloc = file.path( outdir, vn )
    dir.create (outloc, showWarnings=FALSE, recursive =TRUE)

    for (i in 1:length(yrs)){

      outfn = paste( vn, yrs[i], sep=".")
      iy = which(x$year == yrs[i] )

      toplot = x[[vn]] [iy]  

      oo = tapply( toplot, x[["AUID"]] [iy], FUN[[v]], na.rm=TRUE )
 
      polygrid$z = NA
      polygrid$z[ match( names(oo) , polygrid$AUID ) ] = oo

      polygrid$z[ polygrid$z < er[1] ] = er[1] # set levels higher than max datarange to max datarange
      polygrid$z[ polygrid$z > er[2] ] = er[2] # set levels higher than max datarange to max datarange

      # do the map
      fn = file.path( outloc, paste(outfn, outformat, sep="." ) )

      if (plotmethod=="ggplot") {
        
        require(ggplot2)

        o = ggplot() +
          geom_sf( data = polygrid, aes(fill=.data[["z"]]), lwd=0, alpha=0.95  )  +
          coord_sf(xlim = bb$x, ylim =bb$y, expand = FALSE, crs=plot_crs ) +
          scale_fill_gradientn(
            name = yrs[i], 
            limits=range(datarange),
            # colors=colsfn(length(datarange)), 
            colors = colors_brew,  
            na.value=NA ) +
          guides(
            fill = guide_colorbar(
              title.position = "bottom",
              # title.theme = element_blank(), 
              # title.theme = element_text(size = 20),
              label.theme = element_text(size = 14) ) )  +
          additional_features + 
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

      if (plotmethod=="tmap") {

        require(tmap)

        if (outformat=="pdf") pdf( file=fn, width=8, height=6, bg='white', pointsize=20 )
        if (outformat=="svg") svg( filename=fn, width=8, height=6, bg='white', pointsize=20   )
        if (outformat=="png") png( filename=fn, width=1200, height=800, pointsize=20, res=100 )

          o = tm_shape( polygrid, crs=plot_crs ) +
            tm_polygons(
              fill="z",
              fill.scale = tm_scale_continuous(
                  ticks = datarange,  
                  values = "brewer.yl_or_rd",
                  value.na = NA
              ) ,
              fill.legend = tm_legend(
                position=c("left", "top"),  
                frame=TRUE, 
                scale = 1 , 
                title.size=1.5, 
                text.size=0.80, 
                legend.width=0.75,
                na.show=FALSE,
                orientation = "landscape",
                title = yrs[i]
              ),
              col = NULL,
  #            constrast=c(0,0.6),
              lwd = 0.5, 
              col_alpha = 0.5,
              fill_alpha = 0.95 ) +
          tm_shape( coastline, crs=plot_crs ) +
            tm_polygons( col="grey80" ) +
          tm_shape( isobaths, crs=plot_crs ) +
            tm_lines( col="lightgray", col_alpha=0.5) +
          tm_shape( managementlines, crs=plot_crs ) +
            tm_lines( col="grey40", col_alpha=0.6, lwd=2) +
          tm_compass( position=c( "right", "top")) + 
          tm_scalebar( position=c("right", "bottom" ), text.size=0.7, width=20) +
          tm_layout( frame=FALSE ) 
          print(o)
        dev.off()
        print(fn)
      } #tmap
    }  # yr
  }  #var
}


 