
snowcrab_mapping_features = function( p, 
  area_lines=NULL, 
  isobaths=c( 100, 200, 300, 400, 500 ), 
  coastline=c("Canada", "United States of America"), #ignored
  plot_crs=projection_proj4string("lonlat_wgs84"),
  xlim=c(-85,-35), 
  ylim=c(35, 65),
  redo=FALSE, 
  target="ggplot" 
)   {


    if (0) {
      plot_crs=projection_proj4string("lonlat_wgs84")
      xlim=c(-85,-35)
      ylim=c(35, 65)
    }

    # same as carstm::features-to_add, but with different defaults

    fn = file.path( p$project.outputdir, paste0("snowcrab_mapping_features_", target, ".rdz" ) )
    if (!redo){
      O = NULL
      if (file.exists(fn)) O = aegis::read_write_fast(fn)
      if (!is.null(O)) return(O)
    }
 
    rg = NULL
    if (!is.null(area_lines )) {
      rg = area_lines.db( DS=area_lines, returntype="sf", project_to=plot_crs )
    }

    z = NULL
    if (!is.null(isobaths)) {
      z = aegis.bathymetry::isobath_db( depths=isobaths, project_to=plot_crs )
    }

    rg = area_lines.db( DS="cfa.regions", returntype="sf", project_to=plot_crs )
    
    cl = st_transform( polygons_rnaturalearth(xlim=xlim, ylim=ylim), st_crs(plot_crs) )
 
    O = list()

    if ("ggplot" %in% target) {
      require(ggplot2)

      O =  ggplot() +
        geom_sf( data=z,  fill=NA, col = "slategray",  lwd=0.25) +
        geom_sf( data=rg, fill=NA, col = "slategray",  lwd=2.0) + 
        geom_sf( data=cl, fill="gray80", col = "slategray", lwd=0.4)
      O = O[["layers"]]
    }

    if ("tmap" %in% target) {
      require(tmap)

      O =  
        tm_shape( z,  crs=plot_crs ) + tm_lines( col="slategray", col_alpha=0.5, lwd=0.2) +
        tm_shape( rg, crs=plot_crs ) + tm_lines( col="slategray", col_alpha=0.6, lwd=1.5)   + 
        tm_shape( cl, crs=plot_crs ) + tm_borders( col = "slategray", alpha=0.5, lwd=0.5)
    }

    dir.create( p$project.outputdir, showWarnings = FALSE, recursive = TRUE )
    read_write_fast( data=O, fn=fn )
    return(O)
}
