map.survey.locations = function(
  p, 
  basedir=project.datadirectory("bio.snowcrab", "output", "maps", "survey.locations" ),
  years=NULL, 
  plot_crs=st_crs( projection_proj4string("utm20N") ) 
) {

    set = snowcrab.db( DS="set.clean")
    setDT(set)

    if (is.null(years)) years = sort( unique( set$yr ) )
    x = set[ yr %in% years, .(lon, lat, yr ) ]
    x = x[ is.finite( rowSums(x) ) ,]
    x = st_as_sf( x, coords= c("lon", "lat") )
    st_crs(x) =  st_crs( projection_proj4string("lonlat_wgs84") ) 
    
    x = st_transform(x, plot_crs )  # redundant .. in case input data is another projection
  
    if (!file.exists(basedir)) dir.create (basedir, showWarnings=FALSE, recursive =TRUE)
  
    bb = point_to_bbox( p$corners, plot_crs=plot_crs )
  
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
      if ( length(ii)  < 1 ) next()
      xy = x[ ii, ]

      plt = ggplot( ) +
          geom_sf(data=xy, aes(), col="darkgray", lwd=0, cex=3, alpha=0.95) +  
          additional_features +
          labs(caption = paste("Survey locations: ", y)) +
          coord_sf(xlim =bb$x, ylim =bb$y, expand = FALSE, crs=plot_crs) +  #
          local_theme 
      
      fn = file.path(basedir, paste( "survey.locations", y, "png", sep="." ) ) 
      print(fn)

      ggsave(filename=fn, plot=plt,  width=12, height = 8)
 
    }
 
  return ("Done" )
}
