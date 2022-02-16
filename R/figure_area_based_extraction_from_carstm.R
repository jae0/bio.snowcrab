


figure_area_based_extraction_from_carstm = function(DS="temperature")  {

    year.assessment = 2021
    require(bio.snowcrab)   # loadfunctions("bio.snowcrab") 


    # params for number
    pN = snowcrab_parameters(
      project_class="carstm",
      yrs=1999:year.assessment,   
      areal_units_type="tesselation",
      family="poisson",
      carstm_model_label = "1999_present",  # 1999_present is the default anything else and you are on your own
      selection = list(type = "number")
    )


    # params for mean size .. mostly the same as pN
    pW = snowcrab_parameters(
      project_class="carstm",
      yrs=1999:year.assessment,   
      areal_units_type="tesselation",
      carstm_model_label = "1999_present",  # 1999_present is the default anything else and you are on your own
    #   carstm_model_label = paste(   carstm_model_label,   variabletomodel, sep="_")  
      family =  "gaussian" ,  
      selection = list(type = "meansize")
    )


   # first define  polygons to use for extraction:
    year.assessment = 2021
 
    areal_units_proj4string_planar_km = aegis::projection_proj4string("utm20")
    subareas =  c("cfanorth",   "cfa23", "cfa24", "cfa4x" ) 
    ns = length(subareas)

    sc_sppoly = st_union( areal_units( p=pN ) )
    polys = NULL
    for ( i in 1:ns ) {
      subarea = subareas[i]
      print(subarea) # # snowcrab areal units
      au = polygon_managementareas( species="maritimes", area=subarea )
      au = st_transform( st_as_sf(au), st_crs( areal_units_proj4string_planar_km ) )
      au = st_intersection( au, sc_sppoly )
      au$AUID=subarea
      polys = rbind( polys, au )
    }
    # plot(polys)


    if (DS=="number") {
      # number:
      yrs =1999:year.assessment
      tss = aegis_lookup(  
        parameters=pN, 
        LOCS=expand.grid( AUID=polys$AUID, timestamp= yrs + 0.75 ), LOCS_AU=polys, 
        project_class="carstm", output_format="areal_units", 
        variable_name=list( "predictions" ), statvars=c("mean", "sd"), space_resolution=pN$pres,
        returntype = "data.table"
      ) 

      plot( rep(0, length(yrs)) ~ yrs, type="n", ylim=c(2, 9) )
      for (i in 1:ns ){
        lines(predictions_mean~yr, tss[which(tss$AUID==subareas[i]),], type="b", col=i, pch=i)
      }
      legend("topleft", legend=subareas, col=1:ns, pch=1:ns )
    }


    if (DS== "mean_size" ) {
      # mean size:
      yrs =1999:year.assessment
      tss = aegis_lookup(  
        parameters=pW, 
        LOCS=expand.grid( AUID=polys$AUID, timestamp= yrs + 0.75 ), LOCS_AU=polys, 
        project_class="carstm", output_format="areal_units", 
        variable_name=list( "predictions" ), statvars=c("mean", "sd"), space_resolution=pW$pres,
        returntype = "data.table"
      ) 

      plot( rep(0, length(yrs)) ~ yrs, type="n", ylim=c(2, 9) )
      for (i in 1:ns ){
        lines(predictions_mean~yr, tss[which(tss$AUID==subareas[i]),], type="b", col=i, pch=i)
      }
      legend("topleft", legend=subareas, col=1:ns, pch=1:ns )

    }


    if (DS=="temperature") {
      # extract from area-based estimates: temperature
      
      yrs =1970:year.assessment
  
      # default paramerters (copied from 03_temperature_carstm.R )
      require(aegis.temperature)
      params = list( 
        temperature = temperature_parameters( 
          project_class="carstm", 
          yrs=yrs, 
          carstm_model_label="1970_present"
        ) 
      )

      tss = aegis_lookup(  
        parameters=params["temperature"], 
        LOCS=expand.grid( AUID=polys$AUID, timestamp= yrs + 0.75 ), LOCS_AU=polys, 
        project_class="carstm", output_format="areal_units", 
        variable_name=list( "predictions" ), statvars=c("mean", "sd"), space_resolution=pN$pres,
        returntype = "data.table"
      ) 

      basedir = file.path( project.datadirectory("bio.snowcrab" ), "assessments", year.assessment, "timeseries" )
      dir.create (basedir, showWarnings=FALSE, recursive =TRUE)
      
      fn = file.path(basedir,   "temperature_bottom.pdf"  ) 
      print(fn)
      pdf( file=fn )

        plot( rep(0, length(yrs)) ~ yrs, type="n", ylim=c(2, 9) )
        for (i in 1:ns ){
          lines(predictions_mean~yr, tss[which(tss$AUID==subareas[i]),], type="b", col=i, pch=i)
        }
        legend("topleft", legend=subareas, col=1:ns, pch=1:ns )
      dev.off()
    }
 }

