

# Snow crab --- Areal unit modelling of habitat   

# -------------------------------------------------
# Part 1 -- construct basic parameter list defining the main characteristics of the study
   year.assessment = 2021
  require(bio.snowcrab)   # loadfunctions("bio.snowcrab") 
 

  pH = snowcrab_parameters( 
    project_class="carstm", 
    yrs=1999:year.assessment,  
    areal_units_type="tesselation", 
    carstm_model_label = "1999_present_male_commercial",
    selection = list(type = "presence_absence")
  )

  # these area shared across all categories of crab as they are survey station constrained
  sppoly = areal_units( p=pH )  # to reload
  plot( sppoly[, "au_sa_km2"]  )

  # could use the same as pN as it is also caluclated there and so avoid doing another lookup but keeping it separate gives more flexibility
  M = snowcrab.db( p=pH, DS="carstm_inputs"  )  # will redo if not found

  # model
  fit = carstm_model( p=pH, data=M ) # 151 configs and long optim .. 19 hrs
 
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
  res$dyear

  # a few plots
  carstm_plotxy( res, vn=c( "res", "random", "time" ), 
    type="b", ylim=c(0.1, 0.9), xlab="Year", ylab="Probabilty", h=0   )

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
  
  vn="predictions" 
  tmatch = "2021"
  tmout = carstm_map(  res=res, vn=vn, tmatch=tmatch,
    sppoly = sppoly, 
    breaks =brks,
    palette="-RdYlBu",
    # palette="-RdYlBu",
    plot_elements=c(  "compass", "scale_bar", "legend" ),
    map_mode="view",
    tmap_zoom= c(map_centre, map_zoom),
    additional_features=additional_features,
    title=paste("Habitat probability - mature male ",  tmatch )
  )  
  tmout

  if (0) {
      fn_root = paste("Predicted_habitat_probability",  tmatch, sep="_")
      outfilename = file.path( outputdir, paste(fn_root, "png", sep=".") )
      mapview::mapshot( tmap_leaflet(tmout), file=outfilename, vwidth = 1600, vheight = 1200 )  # very slow: consider 
      print(outfilename)
  }
  
  # map all :
  vn "= predictions"

  outputdir = file.path( pH$modeldir, pH$carstm_model_label, "predicted.probability_of_observation" )

  if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

  brks = pretty( c(0,1)  )

  for (y in res$yrs ){

      tmatch = as.character(y)
      fn_root = paste("Predicted_habitat_probability", tmatch, sep="_")
      outfilename = file.path( outputdir, paste(fn_root, "png", sep=".") )

      tmout = carstm_map( res=res, vn=vn, tmatch=tmatch, 
        sppoly = sppoly, 
        breaks = brks,
        palette="-RdYlBu",
        alpha = 0.8,
        plot_elements=c(  "compass", "scale_bar", "legend" ),
        tmap_zoom= c(map_centre, map_zoom),
        map_mode="view",
        additional_features=additional_features,
        title=paste("Habitat probability - mature male ", tmatch ) 
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





 
# end
