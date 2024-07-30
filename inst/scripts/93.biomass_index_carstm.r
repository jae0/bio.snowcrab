


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

  require(bio.snowcrab)   # loadfunctions("bio.snowcrab") 

  year.assessment = 2023
  
  yrs = 1999:year.assessment
  spec_bio = bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 )
  snowcrab_filter_class = "fb"     # fishable biomass (including soft-shelled )  "m.mat" "f.mat" "imm"
 
  
  carstm_model_label= paste( "default", snowcrab_filter_class, sep="_" )

  # params for number
  pN = snowcrab_parameters(
    project_class="carstm",
    yrs=yrs,   
    areal_units_type="tesselation",
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
    carstm_model_label= carstm_model_label,  
    selection = list(
      type = "meansize",
      biologicals=list( spec_bio=spec_bio ),
      biologicals_using_snowcrab_filter_class=snowcrab_filter_class
    )

  )

  # params for probability of observation
  pH = snowcrab_parameters( 
    project_clastheta_inits="carstm", 
    yrs=yrs,  
    areal_units_type="tesselation", 
    carstm_model_label= carstm_model_label,  
    selection = list(
      type = "presence_absence",
      biologicals=list( spec_bio=spec_bio ),
      biologicals_using_snowcrab_filter_class=snowcrab_filter_class
    )
  )


 
  sppoly=areal_units( p=pN )
  
  pN$space_name = sppoly$AUID 
  pN$space_id = 1:nrow(sppoly)  # must match M$space

  pN$time_name = as.character(pN$yrs)
  pN$time_id =  1:pN$ny

  pN$cyclic_name = as.character(pN$cyclic_levels)
  pN$cyclic_id = 1:pN$nw

  pW$space_name = sppoly$AUID 
  pW$space_id = 1:nrow(sppoly)  # must match M$space

  pW$time_name = as.character(pW$yrs)
  pW$time_id =  1:pW$ny

  pW$cyclic_name = as.character(pW$cyclic_levels)
  pW$cyclic_id = 1:pW$nw

  pH$space_name = sppoly$AUID 
  pH$space_id = 1:nrow(sppoly)  # must match M$space

  pH$time_name = as.character(pH$yrs)
  pH$time_id =  1:pH$ny

  pH$cyclic_name = as.character(pH$cyclic_levels)
  pH$cyclic_id = 1:pH$nw

  
  M = snowcrab.db( p=pN, DS="carstm_inputs", sppoly=sppoly, redo=TRUE )  # will redo if not found
  
  additional_features = snowcrab_mapping_features(pN)  # for mapping below
  
   
  

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
      theta=c( 2.409, 1.874, 0.772, 2.092, -1.490, 5.145, 4.509, 2.178, 5.453, 0.182, 2.742, 0.525, 0.051, 0.779 ),
      nposteriors=5000,
      posterior_simulations_to_retain=c( "summary", "random_spatial", "predictions"), 
      family = "poisson",
      verbose=TRUE,
      # redo_fit=FALSE, 
      # debug = "summary",
      # debug = "predictions",
      num.threads="4:3"  
    )


    # mean size
    res = NULL; gc()
    res = carstm_model( p=pW, data=M[ iw, ], sppoly = sppoly, 
      theta=c( 6.108, 8.632, 0.883, 2.946, 9.801, 7.265, 10.726, 12.214, 11.849, 9.826, 6.556, 3.456, 5.832, 2.939, 1.625 ),
      nposteriors=5000,
      posterior_simulations_to_retain=c( "summary", "random_spatial", "predictions"), 
      family =  "gaussian",
      verbose=TRUE,
      # redo_fit=FALSE, 
      # debug = "summary",
      # control.inla = list( strategy="laplace", int.strategy="eb" ),  
      num.threads="4:3" 
    ) 

    # model pa using all data
    res = NULL; gc()
    res = carstm_model( p=pH, data=M, sppoly=sppoly, 
      theta = c( 0.926, 1.743, -0.401, 0.705, -2.574, 1.408, 2.390, 3.459, 3.321, -2.138, 3.083, -1.014, 3.558, 2.703 ),
      nposteriors=5000,
      posterior_simulations_to_retain=c( "summary", "random_spatial", "predictions"), 
      family = "binomial",  # "binomial",  # "nbinomial", "betabinomial", "zeroinflatedbinomial0" , "zeroinflatednbinomial0"
      verbose=TRUE,
      #redo_fit=FALSE, 
      # debug = "summary",
      # control.family=list(control.link=list(model="logit")),  # default for binomial .. no need to specify
      # control.inla = list( strategy="laplace", int.strategy="eb" ),
      num.threads="4:3"
    )
 
  }  # end spatiotemporal model




# ----------------------
# Part 3: assimilation of models


  assimilate_numbers_and_size = TRUE

  if (assimilate_numbers_and_size ) {

    # wgts_max = 1.1 # kg, hard upper limit
    # N_max = NULL
    # #  quantile( M$totno[ipositive]/M$data_offset[ipositive], probs=0.95, na.rm=TRUE )  
    
    # posterior sims 
  
    sims = carstm_posterior_simulations( pN=pN, pW=pW, pH=pH, pa_threshold=0.05, qmax=0.95 )
    sims = sims  / 10^6 # units:  kg ; div (10^6) -> kt ;;
#    sims[ which(!is.finite(sppoly$npts)),, ] = 0

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
      sims = carstm_posterior_simulations( pH=pH, pa_threshold=0.05, qmax=0.95 )
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
 
    print(out)

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
          colors=rev(RColorBrewer::brewer.pal(5, "RdYlBu")),
          outfilename=outfilename
      )
      plt
      
    }
  
  }  # end assimilate size and numbers


 
  
# prep data for discrete version
# Rdata files are ready load them through julia and model
# for production
fishery_model_data_inputs( year.assessment=year.assessment,  type="biomass_dynamics", for_julia=TRUE ) ## note the output directory .. this is used for the next script

# for development
# carstm_results_directory = file.path( homedir, "projects", "dynamical_model", "snowcrab", "data" )
# fishery_model_data_inputs( year.assessment=year.assessment,  type="biomass_dynamics", for_julia=TRUE, save_location=carstm_results_directory )



# end
