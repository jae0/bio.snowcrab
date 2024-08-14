


# -------------------------------------------------
# Snow crab --- Areal unit modelling Hurdle / Delta model  
# combination of three models via posterior simulation
# 1. Poisson on positive valued numbers offset by swept are
# 2. Meansize in space and time 
# 3  Presence-absence
# the convolution of all three after simulation is called a Hurdle or Delta model
# -------------------------------------------------
 

# TODO::: move plotting calls to self-contained functions:


# -------------------------------------------------
# Part 1 -- construct basic parameter list defining the main characteristics of the study

  source( file.path( code_root, "bio_startup.R" )  )
   
  require(bio.snowcrab)   # loadfunctions("bio.snowcrab") 

  year.assessment = 2023
  
  yrs = 1999:year.assessment
  spec_bio = bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 )
  snowcrab_filter_class = "fb"     # fishable biomass (including soft-shelled )  "m.mat" "f.mat" "imm"
  # snowcrab_filter_class = "imm"  # note poisson will not work due to var inflation .. nbinomial is a better choice 
  # snowcrab_filter_class = "f.mat"
  # snowcrab_filter_class = "m.mat"


  
  carstm_model_label= paste( "default", snowcrab_filter_class, sep="_" )

  # params for number
  pN = snowcrab_parameters(
    project_class="carstm",
    yrs=yrs,   
    areal_units_type="tesselation",
    family = switch( snowcrab_filter_class, 
      imm = "nbinomial",
      f.mat= "nbinomial",
      m.mat= "nbinomial",
      fb = "nbinomial",
      "poisson"),  
    carstm_model_label= carstm_model_label,  
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
    carstm_model_label= carstm_model_label,  
    selection = list(
      type = "meansize",
      biologicals=list( spec_bio=spec_bio ),
      biologicals_using_snowcrab_filter_class=snowcrab_filter_class
    )

  )

  # params for probability of observation
  pH = snowcrab_parameters( 
    project_class="carstm", 
    yrs=yrs,  
    areal_units_type="tesselation", 
    family = "binomial",  # "binomial",  # "nbinomial", "betabinomial", "zeroinflatedbinomial0" , "zeroinflatednbinomial0"
    carstm_model_label= carstm_model_label,  
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

    sppoly$dummyvar = ""
    xydata = st_as_sf( xydata, coords=c("lon","lat") )
    st_crs(xydata) = st_crs( projection_proj4string("lonlat_wgs84") )

    additional_features = snowcrab_mapping_features(pN)  # for mapping below
  
    tmap_mode("plot")
    
    plt = 
      tm_shape(sppoly) +
        tm_borders(col = "slategray", alpha = 0.5, lwd = 0.5) + 
        tm_shape( xydata ) + tm_sf() +
        additional_features +
        tm_compass(position = c("right", "TOP"), size = 1.5) +
        tm_scale_bar(position = c("RIGHT", "BOTTOM"), width =0.1, text.size = 0.5) +
        tm_layout(frame = FALSE, scale = 2) +
        tm_shape( st_transform(polygons_rnaturalearth(), st_crs(sppoly) )) + 
        tm_borders(col = "slategray", alpha = 0.5, lwd = 0.5)

    dev.new(width=14, height=8, pointsize=20)
    plt

  }


  if (temperature_figures) {
   
    # area-specific figures
    # /home/jae/bio.data/bio.snowcrab/assessments/2022/timeseries/temperature_bottom.pdf

    figure_area_based_extraction_from_carstm(DS="temperature", year.assessment=year.assessment )  # can only do done once we have an sppoly for snow crab
  

    # full domain:
    # default paramerters (copied from 03_temperature_carstm.R )
    require(aegis.temperature)
    params = list( 
      temperature = temperature_parameters( 
        project_class="carstm", 
        yrs=1970:year.assessment, 
        carstm_model_label="default"
      ) 
    )

    sppoly=areal_units( p=pN )
    pL = aegis.temperature::temperature_parameters( project_class="carstm", carstm_model_label="default" , yrs=p$yrs )

    LUT= aegis_survey_lookuptable( aegis_project="temperature", 
        project_class="carstm", DS="carstm_predictions", pL=pL )

    tss = aegis_lookup(  
      pl=pL, LUT=LUT,
      LOCS=expand.grid( AUID=sppoly$AUID, timestamp= yrs + 0.75 ), LOCS_AU=sppoly, 
      project_class="carstm", output_format="areal_units", 
      variable_name=list( "predictions" ), statvars=c("mean", "sd"), space_resolution=pN$pres,
      returntype = "data.table"
    ) 
   
  }

  sppoly=areal_units( p=pN )
  
  M = snowcrab.db( p=pN, DS="carstm_inputs", sppoly=sppoly, redo=TRUE )  # will redo if not found
  
  additional_features = snowcrab_mapping_featuresres(pN)  # for mapping below
  
   
  

# ------------------------------------------------
# Part 2 -- spatiotemporal statistical model

  if ( spatiotemporal_model ) {

    # total numbers
    sppoly = areal_units( p=pN )
    M = snowcrab.db( p=pN, DS="carstm_inputs", sppoly=sppoly  )  # will redo if not found
  
    io = which(M$tag=="observations")
    ip = which(M$tag=="predictions")
    iq = unique( c( which( M$totno > 0), ip ) )
    iw = unique( c( which( M$totno > 5), ip ) )  # need a good sample to estimate mean size
 
    # number 
    res = NULL; gc()
    res = carstm_model( p=pN, data=M[ iq, ], sppoly=sppoly, 
      space_id = sppoly$AUID,
      time_id =  pN$yrs,
      cyclic_id = pN$cyclic_levels,
      # redo_fit=FALSE, 
      # debug = "summary",
      theta=c(1.544, 2.744, 1.626, -2.066, 3.979, 0.683, 4.899, 4.395, 0.768, -5.015, 1.204, -2.771, 1.526),
      nposteriors=5000,
      posterior_simulations_to_retain=c( "summary", "random_spatial", "predictions"), 
      verbose=TRUE,
      # redo_fit=FALSE, 
      # debug = "summary",
      num.threads="4:3"  
    )

    # mean size
    res = NULL; gc()
    res = carstm_model( p=pW, data=M[ iw, ], sppoly = sppoly, 
      space_id = sppoly$AUID,
      time_id =  pW$yrs,
      cyclic_id = pW$cyclic_levels,
      theta=c(6.309, 7.943, 1.805, 1.646, 8.713, 3.879, 13.561, 10.775, 6.306, -0.488, 6.004, -2.065, 1.391  ),
      nposteriors=5000,
      posterior_simulations_to_retain=c( "summary", "random_spatial", "predictions"), 
      verbose=TRUE,
      num.threads="4:3", 
      # redo_fit=FALSE, 
      # debug = "summary",
      control.inla = list( strategy="laplace", int.strategy="eb" )
    ) 

    # model pa using all data
    res = NULL; gc()
    res = carstm_model( p=pH, data=M, sppoly=sppoly, 
      space_id = sppoly$AUID,
      time_id =  pH$yrs,
      cyclic_id = pH$cyclic_levels,
      theta = c(1.007, 1.793, -4.280, 0.707, -3.044, 2.929, 3.479, -1.262, -1.965, -0.510, -2.144, 2.893 ),
      nposteriors=5000,
      posterior_simulations_to_retain=c( "summary", "random_spatial", "predictions"), 
      verbose=TRUE,
      num.threads="4:3",
      #redo_fit=FALSE, 
      # debug = "summary",
      # control.family=list(control.link=list(model="logit")),  # default for binomial .. no need to specify
      control.inla = list( strategy="laplace", int.strategy="eb" )
    )

    # choose:  
    p = pN
    p = pW
    p = pH

    if (0) {
      # extract results
      fit = carstm_model( p=p, DS="modelled_fit",  sppoly = sppoly )  # extract currently saved model fit
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
   
      x11()
      crs_plot = st_crs( sppoly )
      domain = polygon_managementareas( species="maritimes" )
      domain = st_transform( domain, crs_plot )
      data_mask = st_union( sppoly[which(sppoly$filter==1),1] ) 
      # all = st_union( domain, data_mask )
      nearshore = st_cast( st_difference( domain, data_mask ), "POLYGON")[1]
      domain_new = st_union( data_mask, nearshore )
        
 
      o = carstm_optimal_habitat( 
        res = res,
        xvar = "inla.group(t, method = \"quantile\", n = 11)",  
        yvar = "inla.group(z, method = \"quantile\", n = 11)",
        depths=switch( snowcrab_filter_class, 
          fb = c(100, 350),
          imm = c( 160, 350),
          f.mat = c(100, 160),
          m.mat = c(160, 300)
        ),
        probability_limit = 0.25,
        nsims = 100,
        domain=domain_new 
      ) 
      
      dev.new();
      print( o["depth_plot"] )

      if (0) {
        u = aegis::read_write_fast('/home/jae/tmp/temp_depth_habitat.RDS')
        dev.new()
        plot( habitat~yr, u, type="b", ylim=c(0.1, 0.33))
        lines( habitat_lb~yr, u)
        lines( habitat_ub~yr, u)
        abline(v=1993)
        abline(v=2012)
      
        dev.new()
        plot( habitat_sa~yr, u, type="b" )
        lines( habitat_sa_lb~yr, u)
        lines( habitat_sa_ub~yr, u)
        abline(v=1993)
        abline(v=2012)

        ll = loess(habitat~yr, u, span=0.25 )
        pp = predict( ll, u )
        lines(pp ~ u$yr)

      }

      outputdir = file.path( p$modeldir, p$carstm_model_label )
      fn_optimal = file.path( outputdir, "optimal_habitat_temperature_depth_effect.RDS" )
      read_write_fast( data=o, file=fn_optimal )
      o = aegis::read_write_fast(fn_optimal)

      library(ggplot2)

      dev.new(width=14, height=8, pointsize=20)
      ggplot( o[["temperature_depth"]], aes(yr, habitat ) ) +
        geom_ribbon(aes(ymin=habitat_lb, max=habitat_ub), alpha=0.2, colour=NA) +
        geom_line() +
        labs(x="Year", y="Habitat probabtility", size = rel(1.5)) +
        # scale_y_continuous( limits=c(0, 300) )  
        theme_light( base_size = 22 ) 
      

      dev.new(width=14, height=8, pointsize=20)
      ggplot( o[["temperature_depth"]], aes(yr, habitat_sa ) ) +
        geom_ribbon(aes(ymin=habitat_sa_lb, max=habitat_sa_ub), alpha=0.2, colour=NA) +
        geom_line() +
        labs(x="Year", y=bquote("Habitat surface area;" ~ km^2), size = rel(1.5)) +
        # scale_y_continuous( limits=c(0, 300) )  
        theme_light( base_size = 22 ) 
        
    }



    if (0) {
      # quick plots
      vn=c( "random", "space", "re_total" )
      vn=c( "random", "spacetime", "re_total" )
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
        p=pN
        outputdir = file.path( p$modeldir, p$carstm_model_label, "predicted.numerical.densities" )
        ylab = "Number"
        if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )
        fn_root_prefix = "Predicted_numerical_abundance"
        # title= paste( snowcrab_filter_class, "Number; no./m^2"  )
      }
      if ( meansize) {
        p=pW
        ylab = "Mean weight"
        outputdir = file.path( p$modeldir, p$carstm_model_label, "predicted.meansize" )
        if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )
        fn_root_prefix = "Predicted_meansize"
        # title= paste( snowcrab_filter_class, "Mean weight; kg" ) 
      }
      if ( presence_absence ) {
        p=pH
        ylab = "Probability"
        outputdir = file.path( p$modeldir, p$carstm_model_label, "predicted.presence_absence" )
        if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )
        fn_root_prefix = "Predicted_presence_absence"
        title= paste( snowcrab_filter_class, "Probability")  
      }

      if (0) {
           # choose:  
          #p = pN
          #p = pW
          #p = pH

        # extract results
        fit = carstm_model( p=p, DS="modelled_fit",  sppoly = sppoly )  # extract currently saved model fit
        fit$summary$dic$dic
        fit$summary$dic$p.eff
        plot(fit)
        plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
        plot( fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
        plot( fit$marginals.hyperpar$"Phi for space_time", type="l")  # posterior distribution of phi nonspatial dominates
        plot( fit$marginals.hyperpar$"Precision for space_time", type="l")
        plot( fit$marginals.hyperpar$"Precision for setno", type="l")
           
          fit = carstm_model( p=pN,  DS="modelled_fit")
          fit = carstm_model( p=pH,  DS="modelled_fit")
          fit = carstm_model( p=pW,  DS="modelled_fit")

          all.hypers = INLA:::inla.all.hyper.postprocess(fit$all.hyper)
          hypers = fit$marginals.hyperpar

          carstm_prior_posterior_compare( hypers=hypers, all.hypers=all.hypers, i=2 )  
          carstm_prior_posterior_compare( hypers=hypers, all.hypers=all.hypers, i=3 )  
          carstm_prior_posterior_compare( hypers=hypers, all.hypers=all.hypers, i=5 ) 

          names(hypers)
  
        fit = NULL
  
      }

    oeffdir = file.path(p$modeldir, p$carstm_model_label, "figures")
    carstm_plot_marginaleffects( p=p, outputdir=oeffdir, fn_root_prefix=fn_root_prefix ) 

  
    additional_features = features_to_add( 
      p=p, 
#      area_lines="cfa.regions",
      isobaths=c( 100, 200, 300, 400, 500  ), 
      xlim=c(-80,-40), 
      ylim=c(38, 60) # ,redo=TRUE 
    )

  # map all bottom temps: 
  
    outputdir = file.path(p$modeldir, p$carstm_model_label, "maps" )

    fn_root_prefix = "Substrate grainsize (mm)"
     
    carstm_plot_map( p=p, outputdir=outputdir, 
      additional_features=additional_features, 
      toplot="random_spatial", probs=c(0.025, 0.975),  transf=log10, 
      colors=rev(RColorBrewer::brewer.pal(5, "RdYlBu")) ) 

    carstm_plot_map( p=p, outputdir=outputdir, additional_features=additional_features, 
      toplot="predictions", colors=rev(RColorBrewer::brewer.pal(5, "RdYlBu")) )

 
    } #end each seletion


  }  # end spatiotemporal model




# ----------------------
# Part 3: assimilation of models


  assimilate_numbers_and_size = TRUE

  if (assimilate_numbers_and_size ) {

    # wgts_max = 1.1 # kg, hard upper limit
    # N_max = NULL
    # #  quantile( M$totno[ipositive]/M$data_offset[ipositive], probs=0.95, na.rm=TRUE )  
    
    # posterior sims 
  
    sims = carstm_posterior_simulations( pN=pN, pW=pW, pH=pH, pa_threshold=0.05, qmax=0.99 )
    sims = sims  / 10^6 # units:  kg ; div (10^6) -> kt ;;
   
    SM = aggregate_simulations( 
      sims=sims, 
      sppoly=sppoly, 
      fn=carstm_filenames( pN, returnvalue="filename", fn="aggregated_timeseries" ), 
      yrs=pN$yrs, 
      method="sum", 
      redo=TRUE 
    ) 
    # units: kt/km^2

    if (0) {
      # to compute habitat prob
      sims = carstm_posterior_simulations( pH=pH, pa_threshold=0.05, qmax=0.99 )
      SM = aggregate_simulations( 
        sims=sims, 
        sppoly=sppoly, 
        fn=carstm_filenames( pN, returnvalue="filename", fn="aggregated_timeseries" ), 
        yrs=pN$yrs, 
        method="mean", 
        redo=TRUE 
      ) 
      outputdir = file.path( carstm_filenames( pN, returnvalue="output_directory"), "aggregated_habitat_timeseries" )
      RES= SM$RES  
 
    }      
    
    RES= SM$RES  # units: kt
    # RES = aggregate_simulations( fn=carstm_filenames( pN, returnvalue="filename", fn="aggregated_timeseries" ) )$RES

    # note: using pN, even though this is biomass 
    
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


    regions = c("cfanorth", "cfasouth",  "cfa4x" )
    region_label = c("N-ENS", "S-ENS", "4X")
 
    a= cbind( "cfanorth", RES[,c("yrs", "cfanorth", "cfanorth_lb", "cfanorth_ub")] )
    b= cbind( "cfasouth", RES[,c("yrs", "cfasouth", "cfasouth_lb", "cfasouth_ub")] )
    c= cbind( "cfa4x", RES[,c("yrs", "cfa4x", "cfa4x_lb", "cfa4x_ub")] )
    names(a) = names(b) = names(c) = c("region", "year", "mean", "lb", "ub")
    tdb = rbind(a, b, c)

    tdb$region = factor(tdb$region, levels=regions, labels =region_label)
    tdb = tdb[(which(!is.na(tdb$region))), ]
   
    fn = file.path( outputdir, "biomass_M0.png" )
    
    require(ggplot2)
    library(ggbreak) 

    color_map = c("#E69F00", "#56B4E9",  "#CC79A7" )  
    
    out = ggplot(tdb, aes(x=year, y=mean, fill=region, colour=region)) +
      geom_line( alpha=0.9, linewidth=1.2 ) +
      geom_point(aes(shape=region), size=3, alpha=0.7 ) +
      geom_errorbar(aes(ymin=lb,ymax=ub), linewidth=0.8, alpha=0.8, width=0.3)  +
      labs(x="Year/Année", y="Biomass index (kt) / Indice de biomasse (kt)", size = rel(1.5)) +
      scale_colour_manual(values=color_map) +
      scale_fill_manual(values=color_map) +
      scale_shape_manual(values = c(15, 17, 19)) +
      theme_light( base_size = 22) + 
      theme( legend.position="inside", legend.position.inside=c(0.75, 0.9), legend.title=element_blank()) +
      scale_y_break(c(14, 28), scales = 1)
      
      # scale_y_continuous( limits=c(0, 300) )  
      ggsave(filename=fn, plot=out, device="png", width=12, height = 8)
 


    # map it ..mean density
    sppoly = areal_units( p=pN )  # to reload

    vn = paste("biomass", "predicted", sep=".")

    outputdir = file.path( carstm_filenames( pN, returnvalue="output_directory"), "predicted_biomass_densities" )

    if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

    B = apply( sims, c(1,2), mean )  # sims units (kt);  
    B[ which(!is.finite(B)) ] = NA

    # brks = pretty( log10( quantile( B[], probs=c(0.05, 0.95), na.rm=TRUE )* 10^6)  )
    sa = units::drop_units(sppoly$au_sa_km2)
    brks = pretty( ( quantile( log(B * 10^6 / sa), probs=c(0.05, 0.95), na.rm=TRUE ))  )
  
    additional_features = snowcrab_mapping_features(pN)  # for mapping below

    for (i in 1:length(pN$yrs) ) {
      y = as.character( pN$yrs[i] )
      # u = log10( B[,y]* 10^6 )   ## Total kt->kg: log10( kg )
      u = log( B[,y]* 10^6 / sa) # ;; density  log10( kg /km^2 )
      
      sppoly[,vn] = u
      outfilename = file.path( outputdir , paste( "biomass", y, "png", sep=".") )
      plt =  carstm_map(  sppoly=sppoly, vn=vn,
          breaks=brks,
          additional_features=additional_features,
          title=y,
          # title=paste( "log_10( Predicted biomass density; kg/km^2 )", y ),
          palette="-RdYlBu",
          plot_elements=c( "compass", "scale_bar", "legend" ), 
          outfilename=outfilename
      )
      plt
      
    }
  
  }  # end assimilate size and numbers
 


# end
