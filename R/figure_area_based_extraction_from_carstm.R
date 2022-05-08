


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


   # first define  polygons to use for extraction:
    year.assessment = 2021
 
    areal_units_proj4string_planar_km = aegis::projection_proj4string("utm20")
    subareas =  c("cfanorth", "cfasouth", "cfa4x" ) 
    subareas_label = c( "NENS", "SENS", "4X" )
    
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
      legend("topleft", legend=subareas_label, col=1:ns, pch=1:ns )
    }


    if (DS== "mean_size" ) {
      # mean size:
      yrs =1999:year.assessment
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
      legend("topleft", legend=subareas_label, col=1:ns, pch=1:ns )

    }

    if (DS== "presence_absence" ) {
      # mean size:
      yrs =1999:year.assessment
      # params for mean size .. mostly the same as pN
      pH = snowcrab_parameters(
        project_class="carstm",
        yrs=1999:year.assessment,   
        areal_units_type="tesselation",
        carstm_model_label = "1999_present",  # 1999_present is the default anything else and you are on your own
      #   carstm_model_label = paste(   carstm_model_label,   variabletomodel, sep="_")  
        family =  "binomial" ,  
        selection = list(type = "presence_absence")
      )

      tss = aegis_lookup(  
        parameters=pH, 
        LOCS=expand.grid( AUID=polys$AUID, timestamp= yrs + 0.75 ), LOCS_AU=polys, 
        project_class="carstm", output_format="areal_units", 
        variable_name=list( "predictions" ), statvars=c("mean", "sd"), space_resolution=pW$pres,
        returntype = "data.table"
      ) 

      plot( rep(0, length(yrs)) ~ yrs, type="n", ylim=c(2, 9) )
      for (i in 1:ns ){
        lines(predictions_mean~yr, tss[which(tss$AUID==subareas[i]),], type="b", col=i, pch=i)
      }
      legend("topleft", legend=subareas_label, col=1:ns, pch=1:ns )

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
      
      fn = file.path(basedir, "temperature_bottom.pdf"  ) 
      print(fn)
      pdf( file=fn, width=6, height=8 )

        yran = c( 1, 10 )
        yrs = params$temperature$yrs
        nrep = 5000
        nd = length(yrs)
        y = matrix( NA, nrow=nrep, ncol=nd )
        ncolors = 600
        prs = seq( from=0.025, to=0.975, length.out=ncolors*2 )  # 95% CI
        alpha = seq( from=0.85, to=0.95, length.out=ncolors )

        cols_plot = c("Blues", "Purples",  "Reds", "Dark Mint", "Greens", "Light Grays")
        ncc = length(cols_plot)
        
        cols = matrix( NA, ncol=ncc, nrow=ncolors*2 )
        for (i in 1:ncc ) {
          cols[,i] = c( rev(hcl.colors(ncolors, cols_plot[i], alpha=rev(alpha))), hcl.colors(ncolors, cols_plot[i], alpha=(alpha) ) ) 
        }
 
        layout( matrix( c(1,2,3), ncol=1, nrow=3 ) )
        par(oma=c(6, 6, 6, 1)) # outer margins default:  c(0, 1, 0, 1)'c(bottom, left, top, right)'
        par(mar=c(0, 0, 0.4, 0))

        for (j in 1:ns ) {
          if (j==ns) {
            plot( 0, 0, type="n", ylim=yran, xlim=range(yrs), axes=FALSE, xlab="Year", ylab="Bottom temperature (Celsius)" ) #change xlim to yrs0 to remove 3 yr projection
          } else {
            plot( 0, 0, type="n", ylim=yran, xlim=range(yrs), axes=FALSE, xlab=NULL, ylab=NULL )             
          }

          X = tss[which(tss$AUID==subareas[j]),]
          for (i in 1:nd) y[,i] = rnorm( nrep, mean=X$predictions_mean[i], sd=X$predictions_sd[i])
          Bq =  apply( y, 2, quantile, probs=prs, na.rm=T  )
          for ( k in 1:length(prs) ) lines ( yrs, Bq[k,], lwd=3, col=cols[k,6] )
          lines ( yrs, Bq[ncolors,], lwd=3, col=cols[ncolors, 6]  ) # median
          text( yrs[3] , yran[2]*0.95, subareas_label[j])
          lines( yrs, rep(7, length(yrs)), col="darkred", lwd=2)
          axis(2)
        }
        # legend("topleft", legend=toupper(subareas), lwd=5, col=cols[ncolors,1:ns]  )
        axis(1 )

      dev.off()
 
    }
  
  }
 



