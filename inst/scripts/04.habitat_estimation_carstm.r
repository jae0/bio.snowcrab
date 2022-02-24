

# Snow crab --- Areal unit modelling of habitat   

# -------------------------------------------------
# Part 1 -- construct basic parameter list defining the main characteristics of the study
  year.assessment = 2021
  require(bio.snowcrab)   # loadfunctions("bio.snowcrab") 
  
  runtypes = list()

  runtypes[["R0"]] = list( 
    carstm_label= "1999_present_male_R0", 
    theta=c(-0.790, 3.663, 5.293, 1.211, -1.382, 3.752, 4.514, -1.619, 1.632, -0.959, 2.544, 2.728 ) 
  )

  runtypes[["R1"]] = list( 
    carstm_label= "1999_present_male_R1", 
    theta=c( -0.896, 2.567, -0.420, 0.142, 1.870, 4.258, 5.024, -0.203, 3.068, -0.073, 2.852, 1.503 )
  )
  
  runtypes[["imm"]] = list( 
    carstm_label= "1999_present_imm", 
    theta=c( 1.070, 1.125, 2.441, 2.433, -1.111, 3.215, 4.590, -1.624, 2.283, -0.353, 0.089, 2.646 )
  )
 
  runtypes[["f.mat"]] = list( 
    carstm_label= "1999_present_female_mature", 
    theta=c( 1.822, 3.625, -0.537, 0.218, 2.792, 4.227, 4.881, -0.184, 3.138, -0.082, 2.903, 1.623 )
  )

  # runtypes[["f.berried"]] = list( 
  #   carstm_label= "1999_present_female_berried", 
  #   theta=c( 1.171, 1.573, -0.035, 3.908, 1.883, 4.993, 2.954, -0.882, 0.464, -1.343, 1.763, 2.565  )
  # )
 
  runtypes[["primiparous"]] = list( 
    carstm_label= "1999_present_female_primiparous", 
    theta=c( 1.822, 3.625, -0.537, 0.218, 2.792, 4.227, 4.881, -0.184, 3.138, -0.082, 2.903, 1.623 )
  )

  runtypes[["multiparous"]] = list( 
    carstm_label= "1999_present_female_multiparous", 
    theta=c( 1.822, 3.625, -0.537, 0.218, 2.792, 4.227, 4.881, -0.184, 3.138, -0.082, 2.903, 1.623 )
  )
 

  for (i in 1: length(runtypes) ) {

      pH = snowcrab_parameters( 
        project_class="carstm", 
        yrs=1999:year.assessment,  
        areal_units_type="tesselation", 
        carstm_model_label = runtypes[[i]]$carstm_label,
        family = "binomial",  # "nbinomial", "betabinomial", "zeroinflatedbinomial0" , "zeroinflatednbinomial0"
        carstm_model_inla_control_familiy = list(control.link=list(model='logit')),
        selection = list(
          type = "presence_absence",
          biologicals=list(
            spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 )
          ),
          biologicals_using_snowcrab_filter_class=names(runtypes)[i]
        )
      )

      # these area shared across all categories of crab as they are survey station constrained
      sppoly = areal_units( p=pH )  # to reload
      plot( sppoly[, "au_sa_km2"]  )

      # could use the same as pN as it is also caluclated there and so avoid doing another lookup but keeping it separate gives more flexibility
      M = snowcrab.db( p=pH, DS="carstm_inputs", redo=TRUE  )  # will redo if not found

      # model
      fit = carstm_model( p=pH, data=M, sppoly=sppoly, posterior_simulations_to_retain="predictions",
        theta = runtypes[i]$theta,
        control.inla = list( int.strategy="eb"  ), 
      ) 
  
      fit = NULL
      gc()

      if (0) {
        # very large files .. slow 
        fit = carstm_model( p=pH, DS="carstm_modelled_fit" )  # extract currently saved model fit
        plot(fit)
        plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
      }

      # extract results
      res = carstm_model( p=pH, DS="carstm_modelled_summary"  ) # to load currently saved results
      res$summary$dic$dic
      res$summary$dic$p.eff


      # a few plots
      carstm_plotxy( res, vn=c( "res", "random", "time" ), 
        type="b", ylim=c(0, 1), xlab="Year", ylab="Probabilty", h=0   )

      carstm_plotxy( res, vn=c( "res", "random", "cyclic" ), 
        type="b", col="slategray", pch=19, lty=1, lwd=2.5, ylim=c(0.35, 0.65),
        xlab="Season", ylab="Probabilty" )

      carstm_plotxy( res, vn=c( "res", "random", "inla.group(t, method = \"quantile\", n = 11)" ), 
        type="b", col="slategray", pch=19, lty=1, lwd=2.5, ylim=c(0, 0.8) ,
        xlab="Bottom temperature (degrees Celcius)", ylab="Probabilty" )

      carstm_plotxy( res, vn=c( "res", "random", "inla.group(z, method = \"quantile\", n = 11)" ), 
        type="b", col="slategray", pch=19, lty=1, lwd=2.5, ylim=c(0, 0.9) ,
        xlab="Depth (m)", ylab="Probabilty" )


      # a few maps
      additional_features = snowcrab_features_tmap(pH)  # for mapping below
      map_centre = c( (pH$lon0+pH$lon1)/2  , (pH$lat0+pH$lat1)/2  )
      map_zoom = 7.5
      brks = pretty(  c(0,1)  )
      
      # map all :
      outputdir = file.path( pH$modeldir, pH$carstm_model_label, "predicted.probability_of_observation" )
      if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

      vn = c( "random", "space", "combined" ) 
       tmout = carstm_map(  res=res, vn=vn, 
        sppoly = sppoly, 
        breaks = brks,
        palette="-RdYlBu",
        # palette="-RdYlBu",
        plot_elements=c(  "compass", "scale_bar", "legend" ),
        map_mode="view",
        tmap_zoom= c(map_centre, map_zoom),
        additional_features=additional_features,
        title=paste("Habitat probability",  names(runtypes)[i] )
      )  
      tmout

      fn_root = paste("Predicted_habitat_probability_persistent_spatial_effect", sep="_")
      outfilename = file.path( outputdir, paste(fn_root, "png", sep=".") )
      mapview::mapshot( tmap_leaflet(tmout), file=outfilename, vwidth = 1600, vheight = 1200 )  # very slow: consider 
      print(outfilename)

      vn = "predictions"
      brks = pretty( c(0,1)  )

      for (y in res$yrs ){

        tmatch = as.character(y)
        fn_root = paste("Predicted_habitat_probability", tmatch, sep="_")
        outfilename = file.path( outputdir, paste(fn_root, "png", sep=".") )

        tmout = carstm_map( res=res, vn=vn, tmatch=tmatch, 
          sppoly = sppoly, 
          breaks = brks,
          palette="-RdYlBu",
          alpha = 0.9,
          plot_elements=c(  "compass", "scale_bar", "legend" ),
          tmap_zoom= c(map_centre, map_zoom),
          map_mode="view",
          additional_features=additional_features,
          title=paste("Habitat probability", names(runtypes)[i], tmatch ) 
        )  
        mapview::mapshot( tmap_leaflet(tmout), file=outfilename, vwidth = 1600, vheight = 1200 )  # very slow: consider 
        print(outfilename)

      }
      
      # Aggregations across space:
      snowcrab.db(p=pH, DS="carstm_output_compute" )
      RES = snowcrab.db(p=pH, DS="carstm_output_timeseries" )
      pa = snowcrab.db(p=pH, DS="carstm_output_spacetime_pa"  )

      # plots with 95% PI
      outputdir = file.path( pH$modeldir, pH$carstm_model_label, "aggregated_biomass_timeseries" )
      if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

      (fn = file.path( outputdir, "cfa_all.png"))
      png( filename=fn, width=3072, height=2304, pointsize=12, res=300 )
        plot( cfaall ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Prob of observing snow crab", xlab="",  ylim=c(0,1))
        lines( cfaall_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
        lines( cfaall_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
      dev.off()

      (fn = file.path( outputdir, "cfa_south.png") )
      png( filename=fn, width=3072, height=2304, pointsize=12, res=300 )
        plot( cfasouth ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Prob of observing snow crab", xlab="",  ylim=c(0,1))
        lines( cfasouth_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
        lines( cfasouth_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
      dev.off()

      (fn = file.path( outputdir, "cfa_north.png"))
      png( filename=fn, width=3072, height=2304, pointsize=12, res=300 )
        plot( cfanorth ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Prob of observing snow crab", xlab="",  ylim=c(0,1))
        lines( cfanorth_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
        lines( cfanorth_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
      dev.off()

      (fn = file.path( outputdir, "cfa_4x.png"))
      png( filename=fn, width=3072, height=2304, pointsize=12, res=300 )
        plot( cfa4x ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Prob of observing snow crab", xlab="",  ylim=c(0,1))
        lines( cfa4x_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
        lines( cfa4x_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
      dev.off()

  }

# end
