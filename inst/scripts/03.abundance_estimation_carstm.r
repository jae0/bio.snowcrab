

# Snow crab --- Areal unit modelling of habitat  -- no reliance upon stmv fields


# -------------------------------------------------
# Part 1 -- construct basic parameter list defining the main characteristics of the study
# require(aegis)

  year.assessment = 2020

  p = bio.snowcrab::snowcrab_parameters( 
    project_class="carstm", 
    assessment.years=2000:year.assessment, 
    areal_units_type="lattice",
    selection = list(type = "number")
 )

# ------------------------------------------------
# Part 2 -- spatiotemporal statistical model 

  if ( spataiotemporal_model ) {

    if (0) {
      # polygon structure:: create if not yet made
      # adjust based upon RAM requirements and ncores
      
      p$areal_units_constraint_nmin = 3
      p$areal_units_constraint_ntarget = 10

      for (au in c("cfanorth", "cfasouth", "cfa4x", "cfaall" )) plot(polygon_managementareas( species="snowcrab", au))
      xydata = snowcrab.db( p=p, DS="areal_units_input", redo=TRUE )

      sppoly = areal_units( p=p, redo=TRUE )  # create constrained polygons with neighbourhood as an attribute
      MS = NULL
    }

    sppoly = areal_units( p=p )  # to reload
    plot( sppoly[, "au_sa_km2"]  )


    # -------------------------------------------------
    M = snowcrab.db( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found
    M = NULL; gc()

    fit = carstm_model( p=p, M='snowcrab.db( p=p, DS="carstm_inputs" )' ) # 151 configs and long optim .. 19 hrs
    # fit = carstm_model( p=p, DS="carstm_modelled_fit")

      # extract results
      if (0) {
        # very large files .. slow 
        fit = carstm_model( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit
        plot(fit)
        plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
      }


    res = carstm_model( p=p, DS="carstm_modelled_summary"  ) # to load currently saved results
    res$summary$dic$dic
    res$summary$dic$p.eff
    res$dyear


    plot_crs = p$aegis_proj4string_planar_km
    coastline=aegis.coastline::coastline_db( DS="eastcoast_gadm", project_to=plot_crs )
    isobaths=aegis.bathymetry::isobath_db( depths=c(50, 100, 200, 400 ), project_to=plot_crs )
    managementlines = aegis.polygons::area_lines.db( DS="cfa.regions", returntype="sf", project_to=plot_crs )
  
    time_match = list( year=as.character(2020)  )
    carstm_map(  res=res, 
        vn=paste(p$variabletomodel, "predicted", sep="."), 
        time_match=time_match, 
        coastline=coastline,
        managementlines=managementlines,
        isobaths=isobaths,
        main=paste("Predicted abundance", paste0(time_match, collapse="-") )  
    )
      

    # map all :
    vn = paste(p$variabletomodel, "predicted", sep=".")

    outputdir = file.path( p$modeldir, p$carstm_model_label, "predicted.numerical.densitites" )
    if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

    brks = pretty(  quantile(res[[vn]], probs=c(0,0.975))  )

    for (y in res$year ){

        time_match = list( year=as.character(y)  )
        fn_root = paste("Predicted_abundance", paste0(time_match, collapse="-"), sep="_")
        fn = file.path( outputdir, paste(fn_root, "png", sep=".") )

          carstm_map(  res=res, vn=vn, time_match=time_match, 
            breaks =brks,
            coastline=coastline,
            isobaths=isobaths,
            managementlines=managementlines,
            main=paste("Predicted numerial abundance", paste0(time_match, collapse="-") ),
            outfilename=fn
          )  

    }
    


    snowcrab.db(p=p, DS="carstm_output_compute" )
    
    RES = snowcrab.db(p=p, DS="carstm_output_timeseries" )

    bio = snowcrab.db(p=p, DS="carstm_output_spacetime_biomass" )
    num = snowcrab.db(p=p, DS="carstm_output_spacetime_number" )

    outputdir = file.path( p$modeldir, p$carstm_model_label, "aggregated_biomass_timeseries" )

    if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )


    png( filename=file.path( outputdir, "cfa_all.png"), width=3072, height=2304, pointsize=12, res=300 )
      plot( cfaall ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Biomass index (kt)", xlab="")
      lines( cfaall_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
      lines( cfaall_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    dev.off()


    png( filename=file.path( outputdir, "cfa_south.png"), width=3072, height=2304, pointsize=12, res=300 )
      plot( cfasouth ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Biomass index (kt)", xlab="")
      lines( cfasouth_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
      lines( cfasouth_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    dev.off()


    png( filename=file.path( outputdir, "cfa_north.png"), width=3072, height=2304, pointsize=12, res=300 )
      plot( cfanorth ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Biomass index (kt)", xlab="")
      lines( cfanorth_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
      lines( cfanorth_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    dev.off()


    png( filename=file.path( outputdir, "cfa_4x.png"), width=3072, height=2304, pointsize=12, res=300 )
      plot( cfa4x ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Biomass index (kt)", xlab="")
      lines( cfa4x_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
      lines( cfa4x_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    dev.off()





    # map it ..mean density

    sppoly = areal_units( p=p )  # to reload

    plot_crs = p$aegis_proj4string_planar_km
    coastline=aegis.coastline::coastline_db( DS="eastcoast_gadm", project_to=plot_crs )
    isobaths=aegis.bathymetry::isobath_db( depths=c(50, 100, 200, 400, 800), project_to=plot_crs  )
    managementlines = aegis.polygons::area_lines.db( DS="cfa.regions", returntype="sf", project_to=plot_crs )
  
    vn = paste("biomass", "predicted", sep=".")

    outputdir = file.path( p$modeldir, p$carstm_model_label, "predicted.biomass.densitites" )

    if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

    brks = pretty(  quantile( bio[], probs=c(0,0.975) )* 10^6 )

    for (i in 1:length(p$yrs) ){
      y = as.character( p$yrs[i] )
      sppoly[,vn] = bio[,y] * 10^6
      fn = file.path( outputdir , paste( "biomass", y, "png", sep=".") )

        carstm_map(  sppoly=sppoly, vn=vn,    
          breaks=brks, 
          coastline=coastline,
          isobaths=isobaths,
          managementlines=managementlines,
          main=paste("Predicted biomass density", y ),  
          outfilename=fn
        )
    }

    plot( fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
    plot( fit$marginals.hyperpar$"Phi for space", type="l")  # posterior distribution of phi nonspatial dominates
    plot( fit$marginals.hyperpar$"Precision for space", type="l")
    plot( fit$marginals.hyperpar$"Precision for setno", type="l")



    (res$summary)
 

Time used:
    Pre = 13.7, Running = 11277, Post = 6.19, Total = 11297 
Random effects:
  Name	  Model
    dyri AR1 model
   yr AR1 model
   inla.group(t, method = "quantile", n = 11) RW2 model
   inla.group(z, method = "quantile", n = 11) RW2 model
   inla.group(substrate.grainsize, method = "quantile", n = 11) RW2 model
   inla.group(pca1, method = "quantile", n = 11) RW2 model
   inla.group(pca2, method = "quantile", n = 11) RW2 model
   space BYM2 model
   space_time BYM2 model

Model hyperparameters:
                                                                             mean    sd 0.025quant 0.5quant 0.975quant
Precision for dyri                                                         52.326 0.080     52.163   52.328     52.476
Rho for dyri                                                                0.765 0.000      0.764    0.765      0.765
Precision for yr                                                           55.148 0.102     54.956   55.145     55.355
Rho for yr                                                                  0.761 0.000      0.761    0.761      0.762
Precision for inla.group(t, method = "quantile", n = 11)                   52.842 0.264     52.402   52.833     53.290
Precision for inla.group(z, method = "quantile", n = 11)                   54.802 0.209     54.489   54.771     55.260
Precision for inla.group(substrate.grainsize, method = "quantile", n = 11) 53.714 0.159     53.488   53.689     54.061
Precision for inla.group(pca1, method = "quantile", n = 11)                54.755 0.140     54.556   54.732     55.068
Precision for inla.group(pca2, method = "quantile", n = 11)                54.875 0.227     54.529   54.850     55.319
Precision for space                                                    54.589 0.049     54.488   54.588     54.700
Phi for space                                                           0.049 0.000      0.049    0.049      0.049
Precision for space_time                                                         56.525 0.036     56.452   56.525     56.595
Phi for space_time                                                                0.050 0.000      0.049    0.050      0.050
GroupRho for space_time                                                           0.773 0.000      0.773    0.773      0.773
                                                                             mode
Precision for dyri                                                         52.338
Rho for dyri                                                                0.765
Precision for yr                                                           55.134
Rho for yr                                                                  0.762
Precision for inla.group(t, method = "quantile", n = 11)                   52.521
Precision for inla.group(z, method = "quantile", n = 11)                   54.601
Precision for inla.group(substrate.grainsize, method = "quantile", n = 11) 53.527
Precision for inla.group(pca1, method = "quantile", n = 11)                54.601
Precision for inla.group(pca2, method = "quantile", n = 11)                54.602
Precision for space                                                    54.590
Phi for space                                                           0.049
Precision for space_time                                                         56.527
Phi for space_time                                                                0.049
GroupRho for space_time                                                           0.773

Expected number of effective parameters(stdev): 480.66(15.72)
Number of equivalent replicates : 15.10 

Deviance Information Criterion (DIC) ...............: 51941.31
Deviance Information Criterion (DIC, saturated) ....: 94121.92
Effective number of parameters .....................: -680.14

Watanabe-Akaike information criterion (WAIC) ...: 57902.25
Effective number of parameters .................: 3659.12

Marginal log-Likelihood:  -29248.16 
Posterior marginals for the linear predictor and
 the fitted values are computed



  }


  ########## 


  if (fishery_model) {

      p$fishery_model = fishery_model( DS = "logistic_parameters", p=p, tag=p$areal_units_type )
      p$fishery_model$stancode = stan_initialize( stan_code=fishery_model( p=p, DS="stan_surplus_production" ) )

      str( p$fishery_model)

      p$fishery_model$stancode$compile()

      # res = fishery_model( p=p, DS="samples", tag=p$areal_units_type )  # to get samples
      
      if (0) {
      
        fit = fishery_model( p=p,   DS="fit", tag=p$areal_units_type )  # to get samples
      
        print( fit, max_rows=30 )
        # fit$summary("K", "r", "q")
        
        fit$cmdstan_diagnose()
        fit$cmdstan_summary()
  
          # (penalized) maximum likelihood estimate (MLE) 
        fit_mle =  p$fishery_model$stancode$optimize(data =p$fishery_model$standata, seed = 123)
        fit_mle$summary( c("K", "r", "q") )

        u = stan_extract( as_draws_df(fit_mle$draws() ) )

        mcmc_hist(fit$draws("K")) +
          vline_at(fit_mle$mle(), size = 1.5)

        # Variational Bayes  
        fit_vb = p$fishery_model$stancode$variational( data =p$fishery_model$standata, seed = 123, output_samples = 4000)
        fit_vb$summary(c("K", "r", "q"))

        u = stan_extract( as_draws_df(fit_vb$draws() ) )

        bayesplot_grid(
          mcmc_hist(fit$draws("K"), binwidth = 0.025),
          mcmc_hist(fit_vb$draws("K"), binwidth = 0.025),
          titles = c("Posterior distribution from MCMC", "Approximate posterior from VB")
        )

        color_scheme_set("gray")
        mcmc_dens(fit$draws("K"), facet_args = list(nrow = 3, labeller = ggplot2::label_parsed ) ) + facet_text(size = 14 )   
        # mcmc_hist( fit$draws("K"))

        res = fishery_model( 
          DS="logistic_model", 
          p=p, 
          tag=p$areal_units_type,
          fit = fit_vb
          # from here down are params for cmdstanr::sample()
        )
        
      }



      fit = p$fishery_model$stancode$sample(         
        data=p$fishery_model$standata, 
        iter_warmup = 4000,
        iter_sampling = 2000,
        seed = 123,
        chains = 3,
        parallel_chains = 3,  # The maximum number of MCMC chains to run in parallel.
        max_treedepth = 18,
        adapt_delta = 0.99,
        refresh = 500
      )

      # save fit and get draws
      res = fishery_model( 
        DS="logistic_model", 
        p=p, 
        tag=p$areal_units_type,
        fit = fit
        # from here down are params for cmdstanr::sample()
      )



      color_scheme_set("gray")
      mcmc_dens(res$mcmc, regex_pars="K",  facet_args = list(nrow = 3, labeller = ggplot2::label_parsed ) ) + facet_text(size = 14 )   
      # mcmc_hist( fit$draws("K"))


      # frequency density of key parameters
      fishery_model( DS="plot", vname="K", res=res )
      fishery_model( DS="plot", vname="r", res=res )
      fishery_model( DS="plot", vname="q", res=res, xrange=c(0.5, 3))
      fishery_model( DS="plot", vname="FMSY", res=res  )
      # fishery_model( DS="plot", vname="bosd", res=res  )
      # fishery_model( DS="plot", vname="bpsd", res=res  )

      # timeseries
      fishery_model( DS="plot", type="timeseries", vname="biomass", res=res  )
      fishery_model( DS="plot", type="timeseries", vname="fishingmortality", res=res)

      # Summary table of mean values for inclusion in document
      biomass.summary.table(x)

      # Harvest control rules
      fishery_model( DS="plot", type="hcr", vname="default", res=res  )
      fishery_model( DS="plot", type="hcr", vname="simple", res=res  )

      # diagnostics
      # fishery_model( DS="plot", type="diagnostic.errors", res=res )
      # fishery_model( DS="plot", type="diagnostic.phase", res=res  )


      NN = res$p$fishery_model$standata$N

      # K
      plot.new()
      layout( matrix(c(1,2,3), 3, 1 ))
      par(mar = c(4.4, 4.4, 0.65, 0.75))
      for (i in 1:3) plot(density(res$mcmc$K[,i] ), main="")
      ( qs = apply(  res$mcmc$K[,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) )

      # R
      plot.new()
      layout( matrix(c(1,2,3), 3, 1 ))
      par(mar = c(4.4, 4.4, 0.65, 0.75))
      for (i in 1:3) plot(density(res$mcmc$r[,i] ), main="")
      ( qs = apply(  res$mcmc$r[,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) )

      # q
      plot.new()
      layout( matrix(c(1,2,3), 3, 1 ))
      par(mar = c(4.4, 4.4, 0.65, 0.75))
      for (i in 1:3) plot(density(res$mcmc$q[,i] ), main="")
      ( qs = apply(  res$mcmc$q[,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) )

      # FMSY
      plot.new()
      layout( matrix(c(1,2,3), 3, 1 ))
      par(mar = c(4.4, 4.4, 0.65, 0.75))
      for (i in 1:3) plot(density(res$mcmc$FMSY[,i] ), main="")
      ( qs = apply(  res$mcmc$FMSY[,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) )


      # densities of biomass estimates for the year.assessment
      plot.new()
      layout( matrix(c(1,2,3), 3, 1 ))
      par(mar = c(4.4, 4.4, 0.65, 0.75))
      for (i in 1:3) plot(density(res$mcmc$B[,NN,i] ), main="")
      ( qs = apply(  res$mcmc$B[,NN,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) )

      # densities of biomass estimates for the previous year
      plot.new()
      layout( matrix(c(1,2,3), 3, 1 ))
      par(mar = c(4.4, 4.4, 0.65, 0.75))
      for (i in 1:3) plot(density( res$mcmc$B[,NN-1,i] ), main="")
      ( qs = apply(  res$mcmc$B[,NN-1,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) )

      # densities of F in assessment year
      plot.new()
      layout( matrix(c(1,2,3), 3, 1 ))
      par(mar = c(4.4, 4.4, 0.65, 0.75))
      for (i in 1:3) plot(density(  res$mcmc$F[,NN,i] ), xlim=c(0.01, 0.6), main="")
      ( qs = apply(  res$mcmc$F[,NN,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) )
      ( qs = apply(  res$mcmc$F[,NN,], 2, mean ) )

      # densities of F in previous year
      plot.new()
      layout( matrix(c(1,2,3), 3, 1 ))
      par(mar = c(4.4, 4.4, 0.65, 0.75))
      for (i in 1:3) plot(density(  res$mcmc$F[,NN-1,i] ), xlim=c(0.01, 0.6), main="")
      ( qs = apply(  res$mcmc$F[,NN-1,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) )
      ( qs = apply(  res$mcmc$F[,NN-1,], 2, mean ) )

      # F for table ---
      summary( res$mcmc$F, median)
  }


# end
