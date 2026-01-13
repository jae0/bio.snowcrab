

map.observer.locations = function(p, basedir, years=NULL,
  plot_crs=st_crs( projection_proj4string("utm20")  ) ) {

  x = observer.db( DS="odb")
  x = x[ polygon_inside(x, region="isobath1000m"), ]

  x$yr = x$fishyr  # use fishyr and not the calendar year ###################
  setDT(x)

  x = x[, .(yr, lon, lat)]
  x = x[ is.finite( rowSums(x) ) ,]

  x = st_as_sf( x, coords= c("lon", "lat") )
  st_crs(x) =  st_crs( projection_proj4string("lonlat_wgs84") ) 
  
  x = st_transform(x, plot_crs )  # redundant .. in case input data is another projection

  if (is.null(years)) years = sort( unique( x$yr ) )

  if (!file.exists(basedir)) dir.create (basedir, showWarnings=FALSE, recursive =TRUE)
 
  bb = c(p$corners$lon, p$corners$lat)
 
  additional_features = snowcrab_mapping_features(p, redo=FALSE ) 

  local_theme = theme(
      axis.line=element_blank(),
      axis.ticks=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(), 
      legend.position = element_blank(),
      #legend.title = element_blank(),
      panel.background = element_rect(fill =NA),
      panel.border=element_blank(),
      panel.grid.major = element_line(color = "grey"),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(), 
      plot.caption = element_text(hjust = 0, size = 14)
  )

  for (y in years) {

    ii =  which(x$yr==y)
    if ( length(ii)  < 10 ) next()
    xy = x[ ii, ]

    plt = ggplot( ) +
        geom_sf(data=xy, aes(), col="darkgray", lwd=0, cex=3, alpha=0.95) +  
        additional_features +
        labs(caption = paste("Observer locations: ", y)) +
        coord_sf(xlim =bb[ 1:2 ], ylim =bb[3:4], expand = FALSE, crs=plot_crs) +  #
        local_theme 
 
    fn = file.path( basedir, paste("observer.locations", y, "png", sep="." ) )
    print(fn)
   
    ggsave(filename=fn, plot=plt,  width=12, height = 8)
 
  }

  return ("Done" )
}
