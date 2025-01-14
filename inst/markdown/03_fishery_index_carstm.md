---
title: "Snow crab Areal unit modelling of fishery data"
author:
  - name: 
      given: Snow Crab Unit
      family: DFO Science
    # orcid: 0000-0003-3632-5723 
    # email: jae.choi@dfo-mpo.gc.ca
    # email: choi.jae.seok@gmail.com
    # corresponding: true
    affiliation: 
      - name: Bedford Institute of Oceanography, Fisheries and Oceans Canada
        city: Dartmouth
        state: NS
        # url: www.bio.gc.ca
# date: "`r format(Sys.time(), '%d %B, %Y')`"
date: last-modified
date-format: "YYYY-MM-D"
keywords: 
  - snow crab fishery performance assessment
  - areal-unit modelling of catch and effort   
abstract: |
  Snow crab modelling of snow crab catch and effort using conditional autoregressive models.   
toc: true
number-sections: true
highlight-style: pygments
# bibliography: media/references.bib  
# csl: media/canadian-journal-of-fisheries-and-aquatic-sciences.csl  # see https://www.zotero.org/styles for more
license: "CC BY"
copyright: 
  holder: Jae S. Choi
  year: 2024
citation: 
  container-title: https://github.com/jae0/bio.snowcrab/
  doi: NA
funding: "The snow crab scientific survey was funded by the snow crab fishers of Maritimes Region of Atlantic Canada."
editor:
  render-on-save: false
format:
  html: 
    code-fold: true
    html-math-method: katex
    embed-resources: true
params:
  year_assessment: 2024
  year_start: 1999
  media_loc: "media"
  sens: 1
  debugging: FALSE

---

 

   
<!-- Preamble

Could be run as an automated process but probably better to run step wise in case of tweaks being needed.

-->

 
    

## Purpose

The snow crab fishery performance is subjected to a Bayesian spatiotemporal model with the [carstm](https://github.com/jae0/carstm) front-end to [INLA](https://www.r-inla.org/), used to perform "non-separable" spatial Conditional autocorrelation (CAR) and temporal (AR1) models. 

  - Poisson model of positive valued numbers offset by swept area in space and time 
  

## Part 1 -- construct basic parameter list defining the main characteristics of the study


```{r}
#| label: setup
#| eval: true
#| output: false

  require(aegis)
  require(bio.snowcrab)   # loadfunctions("bio.snowcrab") 

  media_loc = params$media_loc
  year_assessment = params$year_assessment
  year_start = params$year_start

  yrs = year_start:year_assessment

  spec_bio = bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 )
  
  snowcrab_filter_class = "fb" # fishable biomass (including soft-shelled )  
   
  carstm_model_label= paste( "fishery", snowcrab_filter_class, sep="_" )

  # params for catch and effort
  pB = snowcrab_parameters(
    project_class="carstm",
    yrs=yrs,   
    areal_units_type="tesselation",
    carstm_model_label= carstm_model_label,  
    selection = list(
      type = "biomass",
      biologicals=list( spec_bio=spec_bio ),
      biologicals_using_snowcrab_filter_class=snowcrab_filter_class
    ),
    variabletomodel = "catch" 
  )
  
  pB$formula = update.formula( pB$formula, catch ~ . + offset( data_offset ) ) 
   
  # areal units upon which carstm will operate ... this is made in 01.snowcrab.r
  sppoly=areal_units( p=pB )
  
  pB$space_name = sppoly$AUID 
  pB$space_id = 1:nrow(sppoly)  # must match M$space

  pB$time_name = as.character(pB$yrs)
  pB$time_id =  1:pB$ny

  pB$cyclic_name = as.character(pB$cyclic_levels)
  pB$cyclic_id = 1:pB$nw
 

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
  M = logbook.db( p=pB, DS="carstm_inputs", sppoly=sppoly, redo=TRUE )  # will redo if not found

```
  

## Part 2 -- spatiotemporal statistical model

With all required parameters defined, the modelling is straightforward.  

```{r}
#| label: carstm-run
#| eval: false
#| output: false
  
    M = logbook.db( p=pB, DS="carstm_inputs", sppoly=sppoly  )

    io = which(M$tag=="observations")
    ip = which(M$tag=="predictions") 
   
    iq = unique( c( which( M$catch > 0), ip ) )
 
    # model pa using all data
    carstm_model( p=pH, data=M, sppoly=sppoly, 
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


    carstm_model( p=pB, data=M[ iq, ], sppoly=sppoly, 
      nposteriors=5000,
      toget = c("summary", "random_spatial", "predictions"),
      posterior_simulations_to_retain = c("predictions"),
      family = "lognormal", # CARSTM does log-transformation internally for "lognormal"
      verbose=TRUE,
      # redo_fit=FALSE, 
      # debug = "summary",
      # debug = "predictions",
      num.threads="4:3"  
    )

    # posterior predictive check
    carstm_posterior_predictive_check(p=pB, M=M[ iq, ]  )

    # EXAMINE POSTERIORS AND PRIORS
    res = carstm_model(  p=pB, DS="carstm_summary" )  # parameters in p and summary

    outputdir = file.path(pB$modeldir, pB$carstm_model_label)

    res_vars = c( names( res$hypers), names(res$fixed) )
    for (i in 1:length(res_vars) ) {
      o = carstm_prior_posterior_compare( res, vn=res_vars[i], outputdir=outputdir )  
      dev.new(); print(o)
    }     

 
  
  # end spatiotemporal model

  # some maps and plots

  # for mapping below, some bathymetry and polygons
  additional_features = snowcrab_mapping_features(pB)  
  
    for (vns in c( "catch", "habitat") ) {
      #vns ="catch"
      #vns ="habitat"
  
      if ( vns=="catch" ) {
        p=pB
        ylab = "catch"
        fn_root_prefix = "Predicted_catch"
        fn_root = "catch"
        # title= paste( snowcrab_filter_class, "catch; kg/trap haul"  )
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
          #     fn=carstm_filenames( pB, returnvalue="filename", fn="aggregated_timeseries" ), 
          #     yrs=pB$yrs, 
          #     method="mean", 
          #     redo=TRUE 
          #   ) 
          #   outputdir = file.path( carstm_filenames( pB, returnvalue="output_directory"), "aggregated_habitat_timeseries" )
          #   ylabel = "Habitat probability"
          #   fn_ts = "habitat_M0.png"
          #   vn = paste("habitat", "predicted", sep=".")
          #   outputdir2 = file.path( carstm_filenames( pB, returnvalue="output_directory"), "predicted_habitat" )
 
      

   
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
  
  additional_features = snowcrab_mapping_features(pB)  
  sppoly = areal_units( p=pB )  # to reload

  for ( vns in c("catch", "habitat") ) {

    if (vns=="catch") {
      sims = carstm_posterior_simulations( pB=pB, pH=pH, pa_threshold=0.05, qmax=0.95 ) 
      SM = aggregate_simulations( 
        sims=sims, 
        sppoly=sppoly, 
        fn=carstm_filenames( pB, returnvalue="filename", fn="aggregated_timeseries" ), 
        yrs=pB$yrs, 
        method="sum", 
        redo=TRUE 
      )  
      # note: using pB, even though this is catch 
      outputdir = file.path( carstm_filenames( pB, returnvalue="output_directory"), "aggregated_timeseries" )
      ylabel = "CPUE (kg/trap haul)"
      fn_ts = "cpue.png"
      vn = paste("catch", "predicted", sep=".")
      outputdir2 = file.path( carstm_filenames( pB, returnvalue="output_directory"), "predicted_cpue" )
    }  


    if (vns=="habitat") {
      sims = carstm_posterior_simulations( pH=pH, pa_threshold=0.05, qmax=0.95 )
      SM = aggregate_simulations( 
        sims=sims, 
        sppoly=sppoly, 
        fn=carstm_filenames( pB, returnvalue="filename", fn="aggregated_timeseries" ), 
        yrs=pB$yrs, 
        method="mean", 
        redo=TRUE 
      ) 
      # units: probability
      # note: using pB, even though this is habitat 
      outputdir = file.path( carstm_filenames( pB, returnvalue="output_directory"), "aggregated_habitat_timeseries" )
      ylabel = "Habitat probability"
      fn_ts = "habitat_M0.png"
      vn = paste("habitat", "predicted", sep=".")
      outputdir2 = file.path( carstm_filenames( pB, returnvalue="output_directory"), "predicted_habitat" )

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

    if (vns=="catch") {
      out = ggplot(tdb, aes(x=year, y=mean, fill=region, colour=region)) +
        geom_line( alpha=0.9, linewidth=1.2 ) +
        geom_point(aes(shape=region), size=3, alpha=0.7 ) +
        geom_errorbar(aes(ymin=lb,ymax=ub), linewidth=0.8, alpha=0.8, width=0.3)  +
        labs(x="Year/Année", y="CPUE (kt/trap haul)", size = rel(1.5)) +
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
      
        for (i in 1:length(pB$yrs) ) {
          y = as.character( pB$yrs[i] )
          # u = log10( B[,y]* 10^6 )   ## Total kt->kg: log10( kg )
          u = log( B[,y]* 10^6 / sa) # ;; density  log10( kg /km^2 )
          
          sppoly[,vn] = u
          outfilename = file.path( outputdir2 , paste( "cpue", y, "png", sep=".") )
          carstm_map(  sppoly=sppoly, vn=vn,
              breaks=brks,
              additional_features=additional_features,
              legend.position.inside=c( 0.1, 0.9 ),
              annotation=y,
              # annotation=paste( "log_10( Predicted cpue; kg/trap haul )", y ),
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

        B = apply( sims, c(1,2), mean )  # sims units (kg/th );  
        B[ which(!is.finite(B)) ] = NA

        # brks = pretty( log10( quantile( B[], probs=c(0.05, 0.95), na.rm=TRUE )* 10^6)  )
        sa = units::drop_units(sppoly$au_sa_km2)
        brks = pretty( c(0, 1)  )
      
        for (i in 1:length(pB$yrs) ) {
          y = as.character( pB$yrs[i] )
          u =  B[,y] 
          sppoly[,vn] = u
          outfilename = file.path( outputdir2 , paste( "habitat", y, "png", sep=".") )
          carstm_map(  sppoly=sppoly, vn=vn,
              breaks=brks,
              additional_features=additional_features,
              legend.position.inside=c( 0.1, 0.9 ),
              annotation=y,
              # annotation=paste( "log_10( Predicted cpue; kg/trap haul )", y ),
              colors=rev(RColorBrewer::brewer.pal(5, "RdYlBu")),
              outfilename=outfilename
          )

        } # end year loop
  
    }

  
  
  }  # end vns loop

   # end assimilate size and numbers
```



# end
