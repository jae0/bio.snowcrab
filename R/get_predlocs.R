get_predlocs = function(p, redo=FALSE) {
      # depth look up can be slow so pre-construct lookup area
      predlocs =NULL
      fn_predlocs = file.path(p$project.outputdir, paste("predlocs_", p$spatial_domain, ".rdz", sep=""))
      if (!redo){
        if (file.exists(fn_predlocs)) {
          predlocs = read_write_fast(fn_predlocs)
          return( predlocs )
        }
      }
 
      predlocs = spatial_grid(p) 
      predlocs = planar2lonlat( predlocs,  proj.type=p$aegis_proj4string_planar_km   )

      pL = aegis.bathymetry::bathymetry_parameters( project_class="stmv"  )
      LUT= aegis_survey_lookuptable( aegis_project="bathymetry", 
        project_class="core", DS="aggregated_data", pL=pL )

      predlocs$z = aegis_lookup( pL=pL, LOCS=predlocs[, c("lon", "lat")],  LUT=LUT,
        output_format="points",  
        variable_name="z.mean", space_resolution=p$pres ) # core=="rawdata"
      setDT(predlocs)
      
      aoi = filter_by_spatial_domain( spatial_domain=p$spatial_domain, Z=predlocs   )
      # aoi = geo_subset( spatial_domain=p$spatial_domain, Z=predlocs   )

      predlocs = list(predlocs=predlocs, aoi=aoi)

      read_write_fast( predlocs, fn=fn_predlocs )

      return( predlocs )
      
      if (0) {
        predlocs = predlocs[ aoi, ]
        lp = levelplot( z ~ plon+plat, data=predlocs[aoi,])
      }
}