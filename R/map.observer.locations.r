

map.observer.locations = function(p, basedir,  map.method="lattice", years=NULL ) {

  odb = observer.db( DS="odb")
  odb$yr = odb$fishyr  # use fishyr and not the real year ###################

  if (is.null(years)) years = sort( unique( odb$yr ) )

  odb = odb[, c("yr", "lon", "lat")]
  odb = odb[ is.finite( rowSums(odb) ) ,]
  
  dir.create (basedir, showWarnings=FALSE, recursive =TRUE)
  corners = data.frame(rbind( cbind( plon=c(220, 990), plat=c(4750, 5270) )))

  if (map.method=="ggplot") additional_features = snowcrab_mapping_features(p, plot_crs=projection_proj4string("lonlat_wgs84")) 
  
  for (y in years) {

    ii =  which(odb$yr==y)
    if ( length(ii)  < 10 ) next()
    toplot = odb[ ii, c("lon", "lat")]
    fn = file.path( basedir, paste("observer.locations", y, "png", sep="." ) )
    print(fn)

    if (map.method=="lattice" ) {
      png( filename=fn, width=3072, height=2304, pointsize=40, res=300 )
        lp = aegis_map( xyz=toplot, spatial_domain=p$spatial_domain, depthcontours=TRUE, annot=y, annot.cex=2.8, corners=corners, plotlines="cfa.regions", pt.cex=1.5 )
        print(lp)
      dev.off()
    }

    if (map.method=="ggplot" ) {
      
      # plt = ggplot() +
      #   coord_sf(xlim =p$corners$lon, ylim =p$corners$lat, expand = FALSE, crs=st_crs(p$plot_crs) ) +
      #   additional_features +
      #   theme(
      #     axis.line=element_blank(),
      #     # axis.text.x=element_blank(),
      #     # axis.text.y=element_blank(),
      #     axis.ticks=element_blank(),
      #     axis.title.x=element_blank(),
      #     axis.title.y=element_blank(), 
      #     legend.position = "inside",
      #     legend.position.inside=legend.position.inside,
      #     legend.title = element_blank(),
      #     # panel.background=element_blank(),
      #     panel.background = element_rect(fill =NA),
      #     panel.border=element_blank(),
      #     # panel.grid.major=element_blank(),
      #     panel.grid.major = element_line(color = "grey"),
      #     panel.grid.minor=element_blank(),
      #     plot.background=element_blank(), 
      #     plot.caption = element_text(hjust = 0, size = 14)
      #   )

    }

  }
  return ("Done" )
}
