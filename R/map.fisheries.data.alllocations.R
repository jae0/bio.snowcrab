
map.fisheries.data.alllocations = function(p){ 
    # maps all fishing locations
    plot_crs  = st_crs(p$aegis_proj4string_planar_km )
    coastline = aegis.coastline::coastline_db( DS="eastcoast_gadm", project_to=plot_crs ) 
    isobaths = aegis.bathymetry::isobath_db( depths=c(50, 100, 200, 400, 800), project_to=plot_crs  )
    managementlines = aegis.polygons::area_lines.db( DS="cfa.regions", returntype="sf", project_to=plot_crs )

    lgbk = logbook.db( DS="fishing.grounds.global",  p=p )
    lgbk$log_effort = log(lgbk$total.effort)
    i = which( is.finite( lgbk$log_effort ))
    lgbk = lgbk[i, ]

    lgbk = st_as_sf(lgbk, coords=c("plon", "plat"), crs=plot_crs )
    lgbk = cbind(lgbk,  st_coordinates(lgbk) )
    xr = p$corners$plon 
    yr = p$corners$plat

    ggplot( ) +
      geom_sf(data = isobaths, col = "grey90", fill = NA)  +
      geom_point( data = lgbk, aes(x=X, y=Y, colour=log_effort), alpha=0.8, size=0.5 ) +  
      scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(5, "RdYlBu")) ) +
      geom_sf(data = managementlines, col = "grey45", fill = NA)  +
      geom_sf(data = coastline, col = "grey45", fill = NA)  +
      coord_sf(xlim =xr, ylim =yr, expand = FALSE) +
      theme(
        axis.line=element_blank(),
        # axis.text.x=element_blank(),
        # axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), 
        legend.position = "inside",
        legend.position.inside= c( 0.925, 0.15 ) ,
        legend.title = element_blank(),
        # panel.background=element_blank(),
        panel.background = element_rect(fill =NA),
        panel.border=element_blank(),
        # panel.grid.major=element_blank(),
        panel.grid.major = element_line(color = "grey"),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(), 
        plot.caption = element_text(hjust = 0, size = 14)
    )

    media_loc = file.path( p$datadir, "media" )
    fn = file.path( media_loc, "snowcrab_fishing.png" )
    ggsave( fn )
}
