
snowcrab_mapping_features = function( p, redo=FALSE ) {

    fn = file.path( p$project.outputdir, "snowcrab_mapping_features.RDS" )
    if (!redo){
      O = NULL
      if (file.exists(fn)) O = readRDS(fn)
      if (!is.null(O)) return(O)
    }

    plot_crs = p$aegis_proj4string_planar_km
    z = aegis.bathymetry::isobath_db(  depths=c( 100, 200, 300, 400, 500 ), project_to=plot_crs  )
    rg = aegis.polygons::area_lines.db( DS="cfa.regions", returntype="sf", project_to=plot_crs )
    cl = st_transform( polygons_rnaturalearth(), st_crs(plot_crs) )

    require(ggplot2)
    require(tmap)

    O = list()
    O[["tmap"]] =  
      tm_shape( z,  projection=plot_crs ) + tm_lines( col="slategray", alpha=0.5, lwd=0.2) +
      tm_shape( rg, projection=plot_crs ) + tm_lines( col="slategray", alpha=0.75, lwd=2)   + 
      tm_shape( cl, projection=plot_crs ) + tm_borders( col = "slategray", alpha=0.5, lwd=0.5)

    O[["ggplot2"]] =  ggplot() +
      geom_sf( data=z,  fill=NA, col = "slategray",  lwd=0.25) +
      geom_sf( data=rg, fill=NA, col = "slategray",  lwd=2.0) + 
      geom_sf( data=cl, fill=NA, col = "slategray", lwd=0.5)


    dir.create( p$project.outputdir, showWarnings = FALSE, recursive = TRUE )
    saveRDS( O, file=fn, compress=TRUE )
    return(O)
}
