
snowcrab_features_tmap = function( p ) {
  require(tmap)
    map_centre = c( (p$lon0+p$lon1)/2  , (p$lat0+p$lat1)/2  )
    map_zoom = 7.5
    plot_crs = p$aegis_proj4string_planar_km
    additional_features =  
      tm_shape( aegis.polygons::area_lines.db( DS="cfa.regions", returntype="sf", project_to=plot_crs ), projection=plot_crs ) + 
        tm_lines( col="slategray", alpha=0.75, lwd=2)   + 
      tm_shape( aegis.bathymetry::isobath_db(  depths=c( seq(0, 400, by=50), 1000), project_to=plot_crs  ), projection=plot_crs ) +
        tm_lines( col="slategray", alpha=0.5, lwd=0.5) 
      # tm_shape( aegis.coastline::coastline_db( DS="eastcoast_gadm", project_to=plot_crs ), projection=plot_crs ) +
        # tm_polygons( col="lightgray", alpha=0.5 , border.alpha =0.5)
    return(additional_features)
}
