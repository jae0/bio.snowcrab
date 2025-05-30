

(Note: though this could be run as an automated process, it is better to run step wise in case of tweaks being needed.)

 
    

## Purpose

The snow crab biomass index is derived from the convolution of three separate Bayesian spatiotemporal model solutions via posterior simulation. This is completed through the [carstm](https://github.com/jae0/carstm) front-end to [INLA](https://www.r-inla.org/), used to perform "non-separable" spatial Conditional autocorrelation (CAR) and temporal (AR1) models. 

  - Poisson model of positive valued numbers offset by swept area in space and time 
  - Mean body size (mass; kg) in space and time 
  - Presence-absence of organisms in space and time, where presence or absence is defined on a quantile threshold of numerical density presence_absence

The convolution of all three after posterior simulation is also known as a Hurdle or Delta model.
  

## Part 1 -- construct basic parameter list defining the main characteristics of the study


```{r}
#| label: setup
#| eval: true
#| output: false

  require(aegis)
  require(bio.snowcrab)   # loadfunctions("bio.snowcrab") 

  if (0) {
    media_loc = params$media_loc
    year_assessment = params$year_assessment
    year_start = params$year_start
  }

  media_loc = file.path("~", "bio.data", "bio.snowcrab" )
  year_assessment = 2024
  year_start = 1999

  yrs = year_start:year_assessment

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
    project_class="carstm", 
    yrs=yrs,  
    areal_units_type="tesselation", 
    carstm_model_label= carstm_model_label,  
    selection = list(
      type = "presence_absence",
      biologicals=list( spec_bio=spec_bio ),
      biologicals_using_snowcrab_filter_class=snowcrab_filter_class
    )
  )
 
  # areal units upon which carstm will operate ... this is made in 01.snowcrab.r
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
  
```

Now that parameter are loaded, we can create (run manually) or reload the input data (later).

```{r}
#| label: input data
#| eval: false
#| output: false
   
  # create model data inputs and the output fields 
  M = snowcrab.db( p=pN, DS="carstm_inputs", sppoly=sppoly, redo=TRUE )  # will redo if not found

```
  

## Part 2 -- spatiotemporal statistical model

With all required parameters defined, the modelling is straightforward. Each variable is modelled with the same covariates. 

```{r}
#| label: carstm-run
#| eval: false
#| output: false

    # total numbers
     M = snowcrab.db( p=pN, DS="carstm_inputs", sppoly=sppoly  )  # will redo if not found
 
    io = which(M$tag=="observations")
    ip = which(M$tag=="predictions")
    iq = unique( c( which( M$totno > 0), ip ) )
    iw = unique( c( which( M$totno > 5), ip ) )  # need a good sample to estimate mean size
 
    # number 
    carstm_model( p=pN, data=M[ iq, ], sppoly=sppoly, 
      # theta=c( 2.7291,1.8146,2.9382,0.0132,3.8666,-0.0211,4.2673,5.5037,6.1421,0.2391,4.2522,0.7666,-0.0100,0.8763 ),
      nposteriors=5000,
      toget = c("summary", "random_spatial", "predictions"),
      posterior_simulations_to_retain = c("predictions"),
      family = "poisson",
      verbose=TRUE,
      # redo_fit=FALSE, 
      # debug = "summary",
      # debug = "predictions",
      num.threads="4:3"  
    )

    # posterior predictive check
    carstm_posterior_predictive_check(p=pN, M=M[ iq, ]  )

    # EXAMINE POSTERIORS AND PRIORS
    res = carstm_model(  p=pN, DS="carstm_summary" )  # parameters in p and summary

    outputdir = file.path(pN$modeldir, pN$carstm_model_label)

    res_vars = c( names( res$hypers), names(res$fixed) )
    for (i in 1:length(res_vars) ) {
      o = carstm_prior_posterior_compare( res, vn=res_vars[i], outputdir=outputdir )  
      dev.new(); print(o)
    }     


    # mean size 
    carstm_model( p=pW, data=M[ iw, ], sppoly = sppoly, 
      # theta=c( 6.0911, 8.6746, 0.9708, 11.4664, -0.0007, 10.6392, 6.7992, 11.4451, 12.4703, 11.5656, 6.6841, 3.4669, 5.8501, 3.1671, 1.7484 ),
      nposteriors=5000,
      toget = c("summary", "random_spatial", "predictions"),
      posterior_simulations_to_retain = c("predictions"),
      family =  "gaussian",
      verbose=TRUE,
      # redo_fit=FALSE, 
      # debug = "summary",
      # control.inla = list(  optimiser="gsl" ),
      # control.inla = list( strategy="laplace", optimiser="gsl", restart=1 ),  # gsl = gsl::bfgs2
      # control.mode = list( restart=TRUE ),
      num.threads="4:3" 
    ) 

    # posterior predictive check
    carstm_posterior_predictive_check(p=pW, M=M[ iw, ]  )

    # EXAMINE POSTERIORS AND PRIORS
    res = carstm_model(  p=pW, DS="carstm_summary" )  # parameters in p and summary
 
    outputdir = file.path(pW$modeldir, pW$carstm_model_label)

    res_vars = c( names( res$hypers), names(res$fixed) )
    for (i in 1:length(res_vars) ) {
      o = carstm_prior_posterior_compare( res, vn=res_vars[i], outputdir=outputdir )  
      dev.new(); print(o)
    }     

    # model pa using all data
    carstm_model( p=pH, data=M, sppoly=sppoly, 
      # theta = c( 0.8917, 2.0052, 4.5021, -0.0000, 1.5400, -2.4689, 1.1762, 2.6536, 2.9546, -2.1406, 3.5352, -0.7465, 3.2443, 2.4420 ),
      nposteriors=5000,
      toget = c("summary", "random_spatial", "predictions"),
      posterior_simulations_to_retain = c("predictions"),
      family = "binomial",  # "binomial",  # "nbinomial", "betabinomial", "zeroinflatedbinomial0" , "zeroinflatednbinomial0"
      verbose=TRUE,
      # redo_fit=FALSE, 
      # debug = "summary",
      
      # control.family=list(control.link=list(model="logit")),  # default for binomial .. no need to specify
      # control.inla = list( strategy="laplace", int.strategy="eb" ),
      num.threads="4:3"
    )
 
    # posterior predictive check
    carstm_posterior_predictive_check(p=pH, M=M[ , ]  )

    # EXAMINE POSTERIORS AND PRIORS
    res = carstm_model(  p=pH, DS="carstm_summary" )  # parameters in p and summary

    outputdir = file.path(pH$modeldir, pH$carstm_model_label)

    res_vars = c( names( res$hypers), names(res$fixed) )
    for (i in 1:length(res_vars) ) {
      o = carstm_prior_posterior_compare( res, vn=res_vars[i], outputdir=outputdir  )  
      dev.new(); print(o)
    }     

  
  # end spatiotemporal model

  # some maps and plots

  # for mapping below, some bathymetry and polygons
  additional_features = snowcrab_mapping_features(pN)  
  
    for (vns in c( "number", "meansize", "habitat") ) {
      #vns ="number"
      #vns ="meansize"
      #vns ="habitat"
  
      if ( vns=="number" ) {
        p=pN
        ylab = "Number"
        fn_root_prefix = "Predicted_numerical_abundance"
        fn_root = "number"
        # title= paste( snowcrab_filter_class, "Number; no./m^2"  )
      }
      if ( vns=="meansize") {
        p=pW
        ylab = "Mean weight"
        fn_root_prefix = "Predicted_meansize"
        fn_root = "weight"
        # title= paste( snowcrab_filter_class, "Mean weight; kg" ) 
      }
      if ( vns=="habitat" ) {
        p=pH
        ylab = "Probability"
        fn_root_prefix = "Predicted_presence_absence"
        fn_root = "habitat"
        # title= paste( snowcrab_filter_class, "Probability")  
      }

      outputdir = file.path( p$modeldir, p$carstm_model_label, "figures" )
      if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )


          #   # to compute habitat prob
          #   sims = carstm_posterior_simulations( pH=pH, pa_threshold=0.05, qmax=0.95 )
          #   SM = aggregate_simulations( 
          #     sims=sims, 
          #     sppoly=sppoly, 
          #     fn=carstm_filenames( pN, returnvalue="filename", fn="aggregated_timeseries" ), 
          #     yrs=pN$yrs, 
          #     method="mean", 
          #     redo=TRUE 
          #   ) 
          #   outputdir = file.path( carstm_filenames( pN, returnvalue="output_directory"), "aggregated_habitat_timeseries" )
          #   ylabel = "Habitat probability"
          #   fn_ts = "habitat_M0.png"
          #   vn = paste("habitat", "predicted", sep=".")
          #   outputdir2 = file.path( carstm_filenames( pN, returnvalue="output_directory"), "predicted_habitat" )
 
      

   
      # to load currently saved results
      res = carstm_model( p=p,  DS="carstm_summary" )  # parameters in p and direct summary
      res$direct
      res = NULL; gc()


      # plots with 95% PI 
      carstm_plot_marginaleffects( p, outputdir, fn_root ) 
   

      # maps of some of the results
      outputdirmap = file.path(p$modeldir, p$carstm_model_label, "maps" )
       
      carstm_plot_map( p=p, outputdir=outputdirmap, fn_root_prefix=fn_root_prefix , additional_features=additional_features, 
        toplot="random_spatial", probs=c(0.025, 0.975) ) 
   
      carstm_plot_map( p=p, outputdir=outputdirmap, fn_root_prefix=fn_root_prefix , additional_features=additional_features, 
        toplot="predictions", probs=c(0.1, 0.9)) 
  
  }


```

 
## Part 3: assimilation of models

Convolution is straightforward as it is operating upon joint posterior simulations. Add more maps/figures as required.

```r

  # wgts_max = 1.1 # kg, hard upper limit
  # N_max = NULL
  # #  quantile( M$totno[ipositive]/M$data_offset[ipositive], probs=0.95, na.rm=TRUE )  
  
  # posterior sims 

  require(ggplot2)
  library(ggbreak) 

  regions = c("cfanorth", "cfasouth",  "cfa4x" )
  region_label = c("N-ENS", "S-ENS", "4X")
  color_map = c("#E69F00", "#56B4E9",  "#CC79A7" )  
  
  additional_features = snowcrab_mapping_features(pN)  
  sppoly = areal_units( p=pN )  # to reload

  for ( vns in c("abundance", "habitat") ) {

    if (vns=="abundance") {
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
      # units: kt
      # note: using pN, even though this is biomass 
      outputdir = file.path( carstm_filenames( pN, returnvalue="output_directory"), "aggregated_biomass_timeseries" )
      ylabel = "Biomass index (kt)"
      fn_ts = "biomass_M0.png"
      vn = paste("biomass", "predicted", sep=".")
      outputdir2 = file.path( carstm_filenames( pN, returnvalue="output_directory"), "predicted_biomass_densities" )
    }  


    if (vns=="habitat") {
      sims = carstm_posterior_simulations( pH=pH, pa_threshold=0.05, qmax=0.95 )
      SM = aggregate_simulations( 
        sims=sims, 
        sppoly=sppoly, 
        fn=carstm_filenames( pN, returnvalue="filename", fn="aggregated_timeseries" ), 
        yrs=pN$yrs, 
        method="mean", 
        redo=TRUE 
      ) 
      # units: probability
      # note: using pN, even though this is habitat 
      outputdir = file.path( carstm_filenames( pN, returnvalue="output_directory"), "aggregated_habitat_timeseries" )
      ylabel = "Habitat probability"
      fn_ts = "habitat_M0.png"
      vn = paste("habitat", "predicted", sep=".")
      outputdir2 = file.path( carstm_filenames( pN, returnvalue="output_directory"), "predicted_habitat" )

    }  

    RES= SM$RES  
    
    if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )
    
    # plot effects

    ( fn = file.path( outputdir, "cfa_all.png") )
    png( filename=fn, width=3072, height=2304, pointsize=12, res=300 )
      plot( cfaall ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab=ylabel, xlab="")
      lines( cfaall_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
      lines( cfaall_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    dev.off()

    ( fn = file.path( outputdir, "cfa_south.png") )
    png( filename=fn, width=3072, height=2304, pointsize=12, res=300 )
      plot( cfasouth ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab=ylabel, xlab="")
      lines( cfasouth_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
      lines( cfasouth_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    dev.off()

    ( fn = file.path( outputdir, "cfa_north.png") )
    png( filename=fn, width=3072, height=2304, pointsize=12, res=300 )
      plot( cfanorth ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab=ylabel, xlab="")
      lines( cfanorth_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
      lines( cfanorth_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    dev.off()

    ( fn = file.path( outputdir, "cfa_4x.png") )
    png( filename=fn, width=3072, height=2304, pointsize=12, res=300 )
      plot( cfa4x ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab=ylabel, xlab="")
      lines( cfa4x_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
      lines( cfa4x_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    dev.off()
 
    a = cbind( "cfanorth", RES[,c("yrs", "cfanorth", "cfanorth_lb", "cfanorth_ub")] )
    b = cbind( "cfasouth", RES[,c("yrs", "cfasouth", "cfasouth_lb", "cfasouth_ub")] )
    c = cbind( "cfa4x", RES[,c("yrs", "cfa4x", "cfa4x_lb", "cfa4x_ub")] )
    names(a) = names(b) = names(c) = c("region", "year", "mean", "lb", "ub")
    tdb = rbind(a, b, c)

    tdb$region = factor(tdb$region, levels=regions, labels =region_label)
    tdb = tdb[(which(!is.na(tdb$region))), ]
  
    fn = file.path( outputdir, fn_ts )

    if (vns=="abundance") {
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
        ggsave(filename=fn, plot=out,  width=12, height = 8)

        print(out)
            # map it ..mean density
    
        if ( !file.exists(outputdir2)) dir.create( outputdir2, recursive=TRUE, showWarnings=FALSE )

        B = apply( sims, c(1,2), mean )  # sims units (kt);  
        B[ which(!is.finite(B)) ] = NA

        # brks = pretty( log10( quantile( B[], probs=c(0.05, 0.95), na.rm=TRUE )* 10^6)  )
        sa = units::drop_units(sppoly$au_sa_km2)
        brks = pretty( ( quantile( log(B * 10^6 / sa), probs=c(0.05, 0.95), na.rm=TRUE ))  )
      
        for (i in 1:length(pN$yrs) ) {
          y = as.character( pN$yrs[i] )
          # u = log10( B[,y]* 10^6 )   ## Total kt->kg: log10( kg )
          u = log( B[,y]* 10^6 / sa) # ;; density  log10( kg /km^2 )
          
          sppoly[,vn] = u
          outfilename = file.path( outputdir2 , paste( "biomass", y, "png", sep=".") )
          carstm_map(  sppoly=sppoly, vn=vn,
              breaks=brks,
              additional_features=additional_features,
              legend.position.inside=c( 0.1, 0.9 ),
              annotation=y,
              # annotation=paste( "log_10( Predicted biomass density; kg/km^2 )", y ),
              colors=rev(RColorBrewer::brewer.pal(5, "RdYlBu")),
              outfilename=outfilename
          )

        } # end year loop

    }

    if (vns=="habitat") {
      out = ggplot(tdb, aes(x=year, y=mean, fill=region, colour=region)) +
        geom_line( alpha=0.9, linewidth=1.2 ) +
        geom_point(aes(shape=region), size=3, alpha=0.7 ) +
        geom_errorbar(aes(ymin=lb,ymax=ub), linewidth=0.8, alpha=0.8, width=0.3)  +
        labs(x="Year/Année", y="Viable habitat (probability) /\nHabitat viable (probabilité)", size = rel(1.0)) +
        scale_colour_manual(values=color_map) +
        scale_fill_manual(values=color_map) +
        scale_shape_manual(values = c(15, 17, 19)) +
        theme_light( base_size = 22) + 
        theme( legend.position="inside", legend.position.inside=c(0.75, 0.9), legend.title=element_blank()) 
        # scale_y_continuous( limits=c(0, 300) )  
        ggsave(filename=fn, plot=out,  width=12, height = 8)

        print(out)
    
        if ( !file.exists(outputdir2)) dir.create( outputdir2, recursive=TRUE, showWarnings=FALSE )

        B = apply( sims, c(1,2), mean )  # sims units (kt);  
        B[ which(!is.finite(B)) ] = NA

        # brks = pretty( log10( quantile( B[], probs=c(0.05, 0.95), na.rm=TRUE )* 10^6)  )
        sa = units::drop_units(sppoly$au_sa_km2)
        brks = pretty( c(0, 1)  )
      
        for (i in 1:length(pN$yrs) ) {
          y = as.character( pN$yrs[i] )
          # u = log10( B[,y]* 10^6 )   ## Total kt->kg: log10( kg )
          u =  B[,y] 

          sppoly[,vn] = u
          outfilename = file.path( outputdir2 , paste( "habitat", y, "png", sep=".") )
          carstm_map(  sppoly=sppoly, vn=vn,
              breaks=brks,
              additional_features=additional_features,
              legend.position.inside=c( 0.1, 0.9 ),
              annotation=y,
              # annotation=paste( "log_10( Predicted biomass density; kg/km^2 )", y ),
              colors=rev(RColorBrewer::brewer.pal(5, "RdYlBu")),
              outfilename=outfilename
          )

        } # end year loop
  
    }

  
  
  }  # end vns loop

   # end assimilate size and numbers
```

      

## Part 4: prep data for integration into the fishery model (discrete version of the logistic model)

The fishery model used for stock assessment is the biomass dynamics model. It is operated through Julia. This step creates the data required.

```r

  ## note the output directory .. this is used for the next script
   # saves in carstm_directory = file.path( modeldir, carstm_model_label )

  fishery_model_data_inputs( 
    year.assessment=year.assessment,   
    type="biomass_dynamics", 
    snowcrab_filter_class="fb",
    modeldir= pN$modeldir,
    carstm_model_label=carstm_model_label,
    for_julia=TRUE,
    fishery_model_label="turing1"  
  ) 

  
```

Rdata files are ready to load through Julia. Continue to step [04_snowcrab_fishery_model_turing.md](04_snowcrab_fishery_model_turing.md) to complete the assessment.

# end
