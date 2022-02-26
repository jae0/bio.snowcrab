
snowcrab_features_tmap = function( p ) {
  require(tmap)
  plot_crs = p$aegis_proj4string_planar_km
  O =  
    tm_shape( aegis.bathymetry::isobath_db(  depths=c(  10, 50, 100, 200, 300, 400, 800 ), project_to=plot_crs  ), projection=plot_crs ) +
      tm_lines( col="slategray", alpha=0.5, lwd=0.5) +
    tm_shape( aegis.polygons::area_lines.db( DS="cfa.regions", returntype="sf", project_to=plot_crs ), projection=plot_crs ) + 
      tm_lines( col="slategray", alpha=0.75, lwd=2)   + 
    tm_shape( st_transform( polygons_rnaturalearth(), st_crs(plot_crs) ), projection=plot_crs ) +
      tm_borders( col = "slategray", alpha=0.75, lwd=2)
  return(O)
}
