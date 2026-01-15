
map.fisheries.data.all.locations = function(p, plot_crs = st_crs( projection_proj4string("utm20N")) ) { 
    # maps all fishing locations

    lgbk = logbook.db( DS="fishing.grounds.global",  p=p )  # already gridded
    setDT(lgbk)
    
    lgbk$log_effort = log(lgbk$total.effort)
    i = which( is.finite( lgbk$log_effort ))
    lgbk = lgbk[i, .(plon, plat, log_effort) ]

    lgbk = st_as_sf(lgbk, coords=c("plon", "plat"), crs=p$aegis_proj4string_planar_km )
    lgbk = st_transform(lgbk, crs=plot_crs)
    lgbk = cbind(lgbk,  st_coordinates(lgbk) )  #X, Y
 
    bb = point_to_bbox( p$corners, plot_crs=plot_crs )
 
    o = ggplot() +
      geom_point( data = lgbk, aes(x=X, y=Y, colour=log_effort), alpha=0.99, size=0.5 ) +  
      scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu")) ) +
      coord_sf(xlim =bb$x, ylim =bb$y, expand = FALSE, crs=plot_crs, default_crs=plot_crs) +
      snowcrab_mapping_features(p, redo=FALSE )  + 
      theme(
        axis.line=element_blank(),
        # axis.text.x=element_blank(),
        # axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), 
        legend.position = "inside",
        legend.position.inside= c( 0.925, 0.15 ) ,
        legend.title = element_text("Effort", vjust = +2),
        # panel.background=element_blank(),
        panel.background = element_rect(fill ="white"),
        panel.border=element_blank(),
        # panel.grid.major=element_blank(),
        panel.grid.major = element_line(color = "grey"),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(), 
        plot.caption = element_text(hjust = 0, size = 12)
    )

    media_loc = file.path( p$datadir, "media" )
    fn = file.path( media_loc, "snowcrab_fishing.png" )
    print(fn)
    ggsave( fn, o )
}
