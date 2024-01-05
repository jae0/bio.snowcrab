

# Snow crab --- Areal unit modelling Delta model variation
# comination of three models via posterior simulation
# 1. Poisson on positive valued numbers offset by swept are
# 2. Meansize in space and time 
# 3  Presence-absence
# the convolution of all three after simulation is sometimes called a Hurdle or Delta model

# For any that pick this up, this works well 
# It was retracted due to reviewers not being able to understand 

# TODO::: move plotting calls to self-contained functions:


# -------------------------------------------------
# Part 1 -- construct basic parameter list defining the main characteristics of the study


  source( file.path( code_root, "bio_startup.R" )  )

  require(bio.snowcrab)   # loadfunctions("bio.snowcrab") 

  year.assessment = 2021
  yrs = 1999:year.assessment
  spec_bio = bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 )


  snowcrab_filter_class = "fb"     # fishable biomass (including soft-shelled )
  # snowcrab_filter_class = "R0"   # fishable biomass (excluding soft-shelled )
  # snowcrab_filter_class = "recruits"  # potential recruitment into fb next year (m11+, R1, R2, soft-shelled > 95mm CW )
  # snowcrab_filter_class = "pre.recruits.i6_i8" # all sexually immature crab (male and female)
  #  snowcrab_filter_class = "m.mat"
  #  snowcrab_filter_class = "f.mat"
  #  snowcrab_filter_class = "imm"

  runlabel = "1999_present_fb"

  # params for number
  pN = snowcrab_parameters(
    project_class="carstm",
    yrs=yrs,   
    areal_units_type="tesselation",
    family="poisson",
    carstm_model_label= runlabel,  
    selection = list(
      type = "number",
      biologicals=list( spec_bio=spec_bio ),
      biologicals_using_snowcrab_filter_class=snowcrab_filter_class
    )
  )


  # params for mean size .. mostly the same as pN
  pW = snowcrab_parameters(
    project_class="carstm",
    yrs=yrs,   
    areal_units_type="tesselation",
    family =  "gaussian",
    carstm_model_label= runlabel,  
    selection = list(
      type = "meansize",
      biologicals=list( spec_bio=spec_bio ),
      biologicals_using_snowcrab_filter_class=snowcrab_filter_class
    )

  )

  # params for rob of observation
  pH = snowcrab_parameters( 
    project_class="carstm", 
    yrs=yrs,  
    areal_units_type="tesselation", 
    family = "binomial",  # "binomial",  # "nbinomial", "betabinomial", "zeroinflatedbinomial0" , "zeroinflatednbinomial0"
    carstm_model_label= runlabel,  
    selection = list(
      type = "presence_absence",
      biologicals=list( spec_bio=spec_bio ),
      biologicals_using_snowcrab_filter_class=snowcrab_filter_class
    )
  )

  
  if (areal_units) {
    # polygon structure:: create if not yet made
    # for (au in c("cfanorth", "cfasouth", "cfa4x", "cfaall" )) plot(polygon_managementareas( species="snowcrab", au))
    xydata = snowcrab.db( p=pN, DS="areal_units_input", redo=TRUE )
    xydata = snowcrab.db( p=pN, DS="areal_units_input" )
    sppoly = areal_units( p=pN, xydata=xydata[ which(xydata$yr %in% pN$yrs), ], redo=TRUE, verbose=TRUE )  # create constrained polygons with neighbourhood as an attribute
  
    sppoly=areal_units( p=pN )

    plot(sppoly["npts"])

    
    additional_features = snowcrab_mapping_featuresresresres(pN)  # for mapping below

    figure_area_based_extraction_from_carstm(DS="temperature", year.assessment )  # can only do done once we have an sppoly for snow crab
  
    M = snowcrab.db( p=pN, DS="carstm_inputs", sppoly=sppoly, redo=TRUE )  # will redo if not found
  
  }

  sppoly=areal_units( p=pN )
  
  

# ------------------------------------------------
# Part 2 -- spatiotemporal statistical model

  if ( spatiotemporal_model ) {

    # total numbers
    sppoly = areal_units( p=pN )
    M = snowcrab.db( p=pN, DS="carstm_inputs", sppoly=sppoly  )  # will redo if not found

    io = which(M$tag=="observations")
    ip = which(M$tag=="predictions")
    iq = unique( c( which( M$totno > 0), ip ) )
    iw = unique( c( which( M$totno > 30), ip ) )  # need a good sample to estimate mean size

    
    # number
    fit = NULL; gc()
    fit = carstm_model( p=pN, data=M[ iq, ], sppoly=sppoly, 
      posterior_simulations_to_retain="predictions", improve.hyperparam.estimates=TRUE
    )

    # mean size
    fit = NULL; gc()
    fit = carstm_model( p=pW, data=M[ iw, ], sppoly = sppoly, 
      posterior_simulations_to_retain="predictions", improve.hyperparam.estimates=TRUE,
      control.inla = list( strategy="laplace", int.strategy="eb" )
    ) 

    # model pa using all data
    fit = NULL; gc()
    fit = carstm_model( p=pH, data=M, sppoly=sppoly, 
      posterior_simulations_to_retain="predictions", improve.hyperparam.estimates=TRUE,
      # control.family=list(control.link=list(model="logit")),  # default
      control.inla = list( strategy="laplace", int.strategy="eb" )
    )

    # choose:  
    p = pN
    p = pW
    p = pH

    if (0) {
      # extract results
      fit = carstm_model( p=p, DS="carstm_modelled_fit",  sppoly = sppoly )  # extract currently saved model fit
      fit$summary$dic$dic
      fit$summary$dic$p.eff
      plot(fit)
      plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
      plot( fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
      plot( fit$marginals.hyperpar$"Phi for space_time", type="l")  # posterior distribution of phi nonspatial dominates
      plot( fit$marginals.hyperpar$"Precision for space_time", type="l")
      plot( fit$marginals.hyperpar$"Precision for setno", type="l")
      fit = NULL
    }

    res = carstm_model( p=p, DS="carstm_modelled_summary",  sppoly = sppoly ) # to load currently saved results

    if (0) {
      p = pH
      res = carstm_model( p=p, DS="carstm_modelled_summary",  sppoly = sppoly ) # to load currently saved results

      o = carstm_2D_effects_probability( 
        res,
        xvar = "inla.group(t, method = \"quantile\", n = 11)",  
        yvar = "inla.group(z, method = \"quantile\", n = 11)", 
        xgrid = seq( -1, 10.5, by=0.5),
        ygrid = seq( 25, 350, by=25),
        xslice = 4,
        yslice = -200,
        nx=200, ny=200,
        theta = 140,
        phi = 15
      )
    
      # use a larger domain than sppoly for the following estimate:
      # sppoly is constrained to sampled locations, and so missing a lot of the inshore areas

      crs_plot = st_crs( p$aegis_proj4string_planar_km )
      domain = polygon_managementareas( species="maritimes" )
      domain = st_transform( domain, crs_plot )
      x11(); plot(domain[1])  
      
      o = carstm_optimal_habitat( 
        res = res,
        xvar = "inla.group(t, method = \"quantile\", n = 11)",  
        yvar = "inla.group(z, method = \"quantile\", n = 11)",
        domain=domain, 
        depths=c(100, 350),   # from visual analysis (carstm_2D_effects_probability)
        temperatures=c(-1.5,6) # from visual analysis (carstm_2D_effects_probability)
      ) 
      
      dev.new();
      print( o["depth_plot"] )

    }



    if (0) {

      vn=c( "random", "space", "combined" )
      vn=c( "random", "spacetime", "combined" )
      vn="predictions"  # numerical density (km^-2)

      tmatch= as.character(year.assessment)

      carstm_map(  res=res, vn=vn, tmatch=tmatch, 
          sppoly = sppoly, 
          palette="-RdYlBu",
          plot_elements=c(  "compass", "scale_bar", "legend" ),
          additional_features=additional_features,
          title =paste( vn, paste0(tmatch, collapse="-"), "no/m^2"  )
      )


      # map all :
      if ( number ) {
        outputdir = file.path( p$modeldir, p$carstm_model_label, "predicted.numerical.densities" )
        if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )
        fn_root_prefix = "Predicted_numerical_abundance"
        fn_root =  "Predicted_numerical_abundance_persistent_spatial_effect" 
        outfilename = file.path( outputdir, paste(fn_root, "png", sep=".") )
        title= paste( snowcrab_filter_class, "Predicted numerical density (no./m^2) - persistent spatial effect"  )
      }
      if ( meansize) {
        outputdir = file.path( p$modeldir, p$carstm_model_label, "predicted.meansize" )
        if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )
        fn_root_prefix = "Predicted_meansize"
        fn_root =  "Predicted_meansize_persistent_spatial_effect" 
        outfilename = file.path( outputdir, paste(fn_root, "png", sep=".") )
        title= paste( snowcrab_filter_class, "Predicted meansize (kg) - persistent spatial effect" ) 
      }
      if ( presence_absence ) {
        outputdir = file.path( p$modeldir, p$carstm_model_label, "predicted.presence_absence" )
        if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )
        fn_root_prefix = "Predicted_presence_absence"
        fn_root =  "Predicted_presence_absence_persistent_spatial_effect" 
        outfilename = file.path( outputdir, paste(fn_root, "png", sep=".") )
        title= paste( snowcrab_filter_class, "Predicted habitat probability - persistent spatial effect")  
      }

      vn = c( "random", "space", "combined" ) 
      toplot = carstm_results_unpack( res, vn )
      brks = pretty(  quantile(toplot[,"mean"], probs=c(0,0.975), na.rm=TRUE )  )

      plt = carstm_map(  res=res, vn=vn, 
        sppoly = sppoly, 
        breaks = brks,
        palette="-RdYlBu",
        plot_elements=c(  "compass", "scale_bar", "legend" ),
        additional_features=additional_features,
        outfilename=outfilename,
        title= title
      )  
      plt
    

      vn="predictions"
      toplot = carstm_results_unpack( res, vn )
      brks = pretty(  quantile(toplot[,,"mean"], probs=c(0,0.975), na.rm=TRUE )  )

      for (y in res$time ){
        tmatch = as.character(y)
        fn_root = paste(fn_root_prefix, paste0(tmatch, collapse="-"), sep="_")
        outfilename = file.path( outputdir, paste(fn_root, "png", sep=".") )

        plt = carstm_map(  res=res, vn=vn, tmatch=tmatch,
          sppoly = sppoly, 
          breaks =brks,
          palette="-RdYlBu",
          plot_elements=c(   "compass", "scale_bar", "legend" ),
          additional_features=additional_features,
          outfilename=outfilename,
          title=paste(fn_root_prefix, snowcrab_filter_class,  paste0(tmatch, collapse="-") )
        )
        plt
        print(outfilename)
      
      }

      # plots with 95% PI
      outputdir = file.path( p$modeldir, p$carstm_model_label, "effects" )
      if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

      (fn = file.path( outputdir, "time.png"))
      png( filename=fn, width=1024, height=1024, pointsize=12, res=196 )
        carstm_plotxy( res, vn=c( "res", "random", "time" ), 
          type="b", ylim=c(0, 1), xlab="Year", ylab="Probabilty", h=0, cex=1.25, cex.axis=1.25, cex.lab=1.25   )
      dev.off()

      (fn = file.path( outputdir, "cyclic.png"))
      png( filename=fn, width=1024, height=1024, pointsize=12, res=196 )
        carstm_plotxy( res, vn=c( "res", "random", "cyclic" ), 
          type="b", col="slategray", pch=19, lty=1, lwd=2.5, ylim=c(0.35, 0.65),
          xlab="Season", ylab="Probabilty", cex=1.25, cex.axis=1.25, cex.lab=1.25   )
      dev.off()


      (fn = file.path( outputdir, "temperature.png"))
      png( filename=fn, width=1024, height=1024, pointsize=12, res=196 )
        carstm_plotxy( res, vn=c( "res", "random", "inla.group(t, method = \"quantile\", n = 11)" ), 
          type="b", col="slategray", pch=19, lty=1, lwd=2.5, ylim=c(0, 0.8) ,
          xlab="Bottom temperature (degrees Celsius)", ylab="Probabilty", cex=1.25, cex.axis=1.25, cex.lab=1.25 )
      dev.off()


      (fn = file.path( outputdir, "depth.png"))
      png( filename=fn, width=1024, height=1024, pointsize=12, res=196 )
        carstm_plotxy( res, vn=c( "res", "random", "inla.group(z, method = \"quantile\", n = 11)" ), 
          type="b", col="slategray", pch=19, lty=1, lwd=2.5, ylim=c(0, 0.9) ,
          xlab="Depth (m)", ylab="Probabilty", cex=1.25, cex.axis=1.25, cex.lab=1.25 )
      dev.off()


      fit = carstm_model( p=pW, DS="carstm_modelled_fit",  sppoly = sppoly ) # to load currently saved results
    
      plot( fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
      plot( fit$marginals.hyperpar$"Phi for space_time", type="l")  # posterior distribution of phi nonspatial dominates
      plot( fit$marginals.hyperpar$"Precision for space_time", type="l")
      plot( fit$marginals.hyperpar$"Precision for setno", type="l")

    }


  }  # end spatiotemporal model


# ----------------------

  assimilate_numbers_and_size = TRUE

  if (assimilate_numbers_and_size ) {

    # wgts_max = 1.1 # kg, hard upper limit
    # N_max = NULL
    # #  quantile( M$totno[ipositive]/M$data_offset[ipositive], probs=0.95, na.rm=TRUE )  
    
    # posterior sims 
  
    sims = carstm_posterior_simulations( pN=pN, pW=pW, pH=pH, sppoly=sppoly, pa_threshold=0.05, qmax=0.99 )
    sims = sims  / 10^6 # 10^6 kg -> kt;; kt/km^2

    
    SM = aggregate_simulations( 
      sims=sims, 
      sppoly=sppoly, 
      fn=carstm_filenames( pN, returnvalue="filename", fn="aggregated_timeseries" ), 
      yrs=pN$yrs, 
      method="sum", 
      redo=TRUE 
    ) 
    
    RES= SM$RES
    # RES = aggregate_simulations( fn=carstm_filenames( pN, returnvalue="filename", fn="aggregated_timeseries" ) )$RES

    outputdir = file.path( carstm_filenames( pN, returnvalue="output_directory"), "aggregated_biomass_timeseries" )

    if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )


    ( fn = file.path( outputdir, "cfa_all.png") )
    png( filename=fn, width=3072, height=2304, pointsize=12, res=300 )
      plot( cfaall ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Biomass index (kt)", xlab="")
      lines( cfaall_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
      lines( cfaall_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    dev.off()


    ( fn = file.path( outputdir, "cfa_south.png") )
    png( filename=fn, width=3072, height=2304, pointsize=12, res=300 )
      plot( cfasouth ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Biomass index (kt)", xlab="")
      lines( cfasouth_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
      lines( cfasouth_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    dev.off()

    ( fn = file.path( outputdir, "cfa_north.png") )
    png( filename=fn, width=3072, height=2304, pointsize=12, res=300 )
      plot( cfanorth ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Biomass index (kt)", xlab="")
      lines( cfanorth_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
      lines( cfanorth_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    dev.off()

    ( fn = file.path( outputdir, "cfa_4x.png") )
    png( filename=fn, width=3072, height=2304, pointsize=12, res=300 )
      plot( cfa4x ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Biomass index (kt)", xlab="")
      lines( cfa4x_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
      lines( cfa4x_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    dev.off()


    # map it ..mean density
    sppoly = areal_units( p=pN )  # to reload

    vn = paste("biomass", "predicted", sep=".")

    outputdir = file.path( carstm_filenames( pN, returnvalue="output_directory"), "predicted_biomass_densitites" )

    if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

    B = apply( sims, c(1,2), mean ) 
    
    brks = pretty( log10( quantile( B[], probs=c(0.05, 0.95) )* 10^6)  )
  
    additional_features = snowcrab_mapping_features(pN)  # for mapping below

    for (i in 1:length(pN$yrs) ){
      y = as.character( pN$yrs[i] )
      sppoly[,vn] = log10( B[,y]* 10^6 )
      outfilename = file.path( outputdir , paste( "biomass", y, "png", sep=".") )
      plt =  carstm_map(  sppoly=sppoly, vn=vn,
          breaks=brks,
          additional_features=additional_features,
          title=paste( "log_10( Predicted biomass density; kg/km^2 )", y ),
          palette="-RdYlBu",
          plot_elements=c( "compass", "scale_bar", "legend" ), 
          outfilename=outfilename
      )
      plt
      
    }
  
  }  # end assimilate size and numbers



 


# end
