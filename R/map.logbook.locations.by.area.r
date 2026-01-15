
map.logbook.locations.by.area = function(p, basedir, years=NULL, mau="subarea", 
    plot_crs=st_crs( projection_proj4string("utm20N")  )) {

    maus = management_areal_units( mau=mau ) 

    additional_features = snowcrab_mapping_features(p, redo=FALSE ) 

    local_theme = theme(
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), 
        legend.position = "bottom",
        #legend.title = element_blank(),
        panel.background = element_rect(fill =NA),
        panel.border=element_blank(),
        panel.grid.major = element_line(color = "grey"),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(), 
        plot.caption = element_text(hjust = 0, size = 14)
    )

    x = logbook.db( DS="logbook" )
    x = x[ polygon_inside(x, region="isobath1000m"), ]
    x = x[, .(lon, lat, yr)] # note: "yr" is fishing year, in 4x: 1999-2000 is yr=1999
    x = x[ is.finite( rowSums(x) ) ,]

    x = st_as_sf( x, coords= c("lon", "lat"))
    st_crs(x) = st_crs( projection_proj4string("lonlat_wgs84") )  
    
    if (!file.exists(basedir)) dir.create (basedir, showWarnings=FALSE, recursive =TRUE)

    bbx = list(
      cfanorth = data.table( lon=c(-60.8, -59), lat=c(47.4, 45.8) ),
      cfa23 = data.table( lon=c( -60.6, -56.8), lat=c(43.6,  46.2) ),
      cfa24 = data.table( lon=c(-63.7, -59), lat=c(43.4, 45.8) ),
      cfasouth = data.table( lon=c(-63.7, -56.8), lat=c(43.3, 46.3) ),
      cfa4x = data.table( lon=c(-65.8, -63.2), lat=c(43.2, 44.7) )
    )
 
    for (r in maus[["internal"]] ) {

        j = polygon_inside(x, region=r )
        y = x[j,]

        y = st_transform(y, plot_crs )  #  in case input data is another projection
        bb = point_to_bbox( bbx[[r]], plot_crs=plot_crs )
    
        y = y[ which(y$yr %in% years ) , ] # note: "yr" is fishing year, in 4x: 1999-2000 is yr=1999
        y$year = as.factor(y$yr) 

        plt = ggplot( ) +
            geom_sf(data=y, aes( fill=year , col=year ), lwd=0, cex=6, alpha=0.3) +  
            scale_color_manual(values=c("yellow", "red")) +
            additional_features +
            labs(caption = "Logbook locations") +
            coord_sf(xlim =bb$x, ylim =bb$y, expand = FALSE, crs=plot_crs ) +  #
            local_theme 
              
        fn = file.path( basedir, paste("logbook_locations_recent_", r, ".png", sep="") )
        print(fn)

        ggsave(filename=fn, plot=plt,  width=12, height = 8)

    } 

  }


 
