
snowcrab_mapping_features = function( 
  p, 
  area_lines=NULL, 
  isobaths=c( 100, 200, 300, 400, 500 ),  
  plot_crs=st_crs( projection_proj4string("utm20") ),  
  redo=FALSE, 
  target="ggplot" 
)   {


    if (0) {
      plot_crs= st_crs( projection_proj4string("utm20") )
      isobaths=c( 100, 200, 300, 400, 500 )
    }

    # same as carstm::features-to_add, but with different defaults

    fn = file.path( p$project.outputdir, paste0("snowcrab_mapping_features_", target, ".rdz" ) )
    if (!redo){
      O = NULL
      if (file.exists(fn)) O = aegis::read_write_fast(fn)
      if (!is.null(O)) return(O)
    }
 
    crs_lonlat = st_crs( projection_proj4string("lonlat_wgs84") )
 
    area_lines = "cfa.regions"
    rg = area_lines.db( DS=area_lines, returntype="sf", project_to=crs_lonlat )
    rg = st_transform( rg, plot_crs )

    z = isobath_db( depths=isobaths, project_to=crs_lonlat, spatial_domain="canada.east" )  # aegis.bathymetry
    z = st_transform(z,  plot_crs )

    cl = polygons_rnaturalearth() 
    cl = st_transform( cl, plot_crs )
     
    require(ggplot2)

    O =  ggplot() +
      geom_sf( data=z,  fill=NA, col = "slategray",  lwd=0.2, alpha=0.4 ) +
      geom_sf( data=rg, fill=NA, col = "slategray",  lwd=1.0, alpha=0.8) + 
      geom_sf( data=cl, fill="gray90", col = "slategray", lwd=0.4, alpha=0.8) 
      
    O = O[["layers"]]
  
    dir.create( p$project.outputdir, showWarnings = FALSE, recursive = TRUE )
    read_write_fast( data=O, fn=fn )
    return(O)
}
