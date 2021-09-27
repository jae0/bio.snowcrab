

# Snow crab --- Areal unit modelling of habitat  -- no reliance upon stmv fields


# -------------------------------------------------
# Part 1 -- construct basic parameter list defining the main characteristics of the study
# require(aegis)

  year.assessment = 2020
  require(bio.snowcrab)   # loadfunctions("bio.snowcrab") 

  p = snowcrab_parameters(
    project_class="carstm",
    assessment.years=2000:year.assessment,
    areal_units_type="tesselation",
#    areal_units_constraint_ntarget = 20,
#    areal_units_constraint_nmin = 5,
    carstm_model_label = "tesselation",   # default is the name of areal_units_type
    selection = list(type = "number")
 )



# ------------------------------------------------
# Part 2 -- spatiotemporal statistical model

  if ( spatiotemporal_model ) {

      if (0) {
        # polygon structure:: create if not yet made
        # adjust based upon RAM requirements and ncores
          # create if not yet made
        for (au in c("cfanorth", "cfasouth", "cfa4x", "cfaall" )) plot(polygon_managementareas( species="snowcrab", au))
        xydata = snowcrab.db( p=p, DS="areal_units_input", redo=TRUE )

        sppoly = areal_units( p=p, hull_alpha=15, redo=TRUE, verbose=TRUE )  # create constrained polygons with neighbourhood as an attribute
        plot( sppoly[, "npts"]  )

        MS = NULL

        # p$carstm_model_label = "tesselation_overdispersed"   # default is the name of areal_units_type
        # p$family  = "zeroinflatedpoisson0" #  "binomial",  # "nbinomial", "betabinomial", "zeroinflatedbinomial0" , "zeroinflatednbinomial0"
        # p$carstm_model_inla_control_familiy = NULL

      }

      sppoly = areal_units( p=p )  # to reload


      # -------------------------------------------------
      M = snowcrab.db( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found
      M = NULL; gc()

      fit = carstm_model( 
        p=p, 
        data='snowcrab.db( p=p, DS="carstm_inputs" )', 
        posterior_simulations_to_retain="predictions" ,
        scale_offsets = TRUE,  # required to stabilize  as offsets are so small, required for : inla.mode = "experimental" 
        # redo_fit = FALSE,  # only to redo sims and extractions 
        # toget="predictions",  # this updates a specific subset of calc
        control.inla = list( strategy='adaptive' ),  # strategy='laplace', "adaptive" int.strategy="eb" 
        num.threads="4:2"
      )

      if (0) {
        # control.compute=list(smtp="default", dic=TRUE, waic=TRUE, cpo=FALSE, config=TRUE, return.marginals.predictor=TRUE)
        # control.fixed=list(prec=1,prec.intercept=1)  
        # 151 configs and long optim .. 19 hrs
        # fit = carstm_model( p=p, DS="carstm_modelled_fit")
  
        # extract results
        # very large files .. slow
          fit = carstm_model( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit
          fit$summary$dic$dic
          fit$summary$dic$p.eff

          plot(fit)
          plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )

      }


      res = carstm_model( p=p, DS="carstm_modelled_summary"  ) # to load currently saved results


      map_centre = c( (p$lon0+p$lon1)/2 - 0.5, (p$lat0+p$lat1)/2 -0.8 )
      map_zoom = 5

      plot_crs = p$aegis_proj4string_planar_km

      
      require(tmap)
      
      additional_features =  
        tm_shape( aegis.polygons::area_lines.db( DS="cfa.regions", returntype="sf", project_to=plot_crs ), projection=plot_crs ) + 
          tm_lines( col="slategray", alpha=0.75, lwd=2)   + 
        tm_shape( aegis.bathymetry::isobath_db(  depths=c( seq(0, 400, by=50), 1000), project_to=plot_crs  ), projection=plot_crs ) +
          tm_lines( col="slategray", alpha=0.5, lwd=0.5) +
        tm_shape( aegis.coastline::coastline_db( DS="eastcoast_gadm", project_to=plot_crs ), projection=plot_crs ) +
          tm_polygons( col="lightgray", alpha=0.5 , border.alpha =0.5)

      (additional_features)

      vn=c( "random", "space", "combined" )
      vn=c( "random", "spacetime", "combined" )
      vn="predictions"  # numerical density (km^-2)

      tmatch="2015"

      carstm_map(  res=res, vn=vn, tmatch=tmatch, 
          palette="-RdYlBu",
          plot_elements=c(  "compass", "scale_bar", "legend" ),
          additional_features=additional_features,
          tmap_zoom= c(map_centre, map_zoom),
          title =paste( vn, paste0(tmatch, collapse="-") )
      )


      # map all :
      outputdir = file.path( p$modeldir, p$carstm_model_label, "predicted.numerical.densitites" )
      if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

      vn="predictions"
      toplot = carstm_results_unpack( res, vn )
      brks = pretty(  quantile(toplot[,,"mean"], probs=c(0,0.975), na.rm=TRUE )  )

     
      for (y in res$time ){
          tmatch = as.character(y)
          fn_root = paste("Predicted_numerical_abundance", paste0(tmatch, collapse="-"), sep="_")
          fn = file.path( outputdir, paste(fn_root, "png", sep=".") )

            o = carstm_map(  res=res, vn=vn, tmatch=tmatch,
              breaks =brks,
              palette="-RdYlBu",
              plot_elements=c(    "compass", "scale_bar", "legend" ),
              additional_features=additional_features,
              title=paste("Predicted numerical density (per km^2) ", paste0(tmatch, collapse="-") ),
              map_mode="plot",
              scale=0.75,
              outformat="tmap",
              outfilename=fn
            )

      }


      fit = meanweights_by_arealunit_modelled( p=p, redo=TRUE, returntype="carstm_modelled_fit" )  ## used in carstm_output_compute

         if (0) {
          fit = meanweights_by_arealunit_modelled( p=p, returntype="carstm_modelled_fit" )  
          
          res = meanweights_by_arealunit_modelled( p=p, returntype="carstm_modelled_summary" )  ## used in carstm_output_compute

          map_centre = c( (p$lon0+p$lon1)/2 - 0.5, (p$lat0+p$lat1)/2 -0.8 )
          map_zoom = 5

          plot_crs = p$aegis_proj4string_planar_km


          additional_features =  
            tm_shape( aegis.polygons::area_lines.db( DS="cfa.regions", returntype="sf", project_to=plot_crs ), projection=plot_crs ) + 
              tm_lines( col="slategray", alpha=0.75, lwd=2)   + 
            tm_shape( aegis.bathymetry::isobath_db( depths=c( seq(50, 400, by=50), 500), project_to=plot_crs  ), projection=plot_crs ) +
              tm_lines( col="slategray", alpha=0.5, lwd=0.5) +
            tm_shape( aegis.coastline::coastline_db( DS="eastcoast_gadm", project_to=plot_crs ), projection=plot_crs ) +
              tm_polygons( col="lightgray", alpha=0.5 , border.alpha =0.5)

          (additional_features)

          vn="predictions"
          tmatch = "2020"

          outputdir = file.path( p$modeldir, p$carstm_model_label, "predicted.mean.weight" )
          if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

          fn_root = paste("Predicted_mean_size", paste0(tmatch, collapse="-"), sep="_")

          fn = file.path( outputdir, paste(fn_root, "png", sep=".") )
   
          carstm_map(  res=res, vn=vn, tmatch=tmatch, 
              palette="-RdYlBu",
              plot_elements=c(  "compass", "scale_bar", "legend" ),
              additional_features=additional_features,
              tmap_zoom= c(map_centre, map_zoom),
              title =paste("Predicted  mean weight of individual (kg)", paste0(tmatch, collapse="-") )
          )

       
        }

      snowcrab.db(p=p, DS="carstm_output_compute" )

      RES = snowcrab.db(p=p, DS="carstm_output_timeseries" )

      bio = snowcrab.db(p=p, DS="carstm_output_spacetime_biomass" )
      num = snowcrab.db(p=p, DS="carstm_output_spacetime_number" )


      outputdir = file.path( p$modeldir, p$carstm_model_label, "aggregated_biomass_timeseries" )

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

      sppoly = areal_units( p=p )  # to reload



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
            additional_features=additional_features,
            title=paste("Predicted biomass density", y ),
            outfilename=fn
          )
      }

      plot( fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
      plot( fit$marginals.hyperpar$"Phi for space_time", type="l")  # posterior distribution of phi nonspatial dominates
      plot( fit$marginals.hyperpar$"Precision for space_time", type="l")
      plot( fit$marginals.hyperpar$"Precision for setno", type="l")


    (res$summary)

Time used:
    Pre = 6.84, Running = 12106, Post = 14.3, Total = 12127
Fixed effects:
             mean   sd 0.025quant 0.5quant 0.975quant  mode kld
(Intercept) 7.333 0.13      7.079    7.333      7.587 7.333   0

Random effects:
  Name	  Model
    dyri AR1 model
   yr AR1 model
   inla.group(t, method = "quantile", n = 11) RW2 model
   inla.group(z, method = "quantile", n = 11) RW2 model
   inla.group(substrate.grainsize, method = "quantile", n = 11) RW2 model
   inla.group(pca1, method = "quantile", n = 11) RW2 model
   inla.group(pca2, method = "quantile", n = 11) RW2 model
   space Besags ICAR model
   space_time BYM2 model

Model hyperparameters:
                                                                             mean    sd 0.025quant 0.5quant 0.975quant
zero-probability parameter for zero-inflated poisson_0                      0.259 0.000      0.259    0.259      0.260
Precision for dyri                                                         53.979 0.463     53.315   53.901     55.031
Rho for dyri                                                                0.772 0.001      0.771    0.772      0.773
Precision for yr                                                           53.530 0.130     53.349   53.508     53.818
Rho for yr                                                                  0.753 0.001      0.752    0.753      0.754
Precision for inla.group(t, method = "quantile", n = 11)                   53.695 0.158     53.405   53.697     53.948
Precision for inla.group(z, method = "quantile", n = 11)                   55.764 0.176     55.425   55.770     56.041
Precision for inla.group(substrate.grainsize, method = "quantile", n = 11) 53.831 0.179     53.542   53.826     54.131
Precision for inla.group(pca1, method = "quantile", n = 11)                53.898 0.239     53.543   53.870     54.359
Precision for inla.group(pca2, method = "quantile", n = 11)                55.105 0.235     54.733   55.090     55.520
Precision for space                                                    52.718 0.035     52.655   52.715     52.792
Precision for space_time                                                         56.221 0.035     56.153   56.222     56.288
Phi for space_time                                                                0.048 0.000      0.048    0.048      0.048
GroupRho for space_time                                                           0.772 0.000      0.772    0.772      0.773
                                                                             mode
zero-probability parameter for zero-inflated poisson_0                      0.259
Precision for dyri                                                         53.503
Rho for dyri                                                                0.771
Precision for yr                                                           53.374
Rho for yr                                                                  0.752
Precision for inla.group(t, method = "quantile", n = 11)                   53.560
Precision for inla.group(z, method = "quantile", n = 11)                   55.658
Precision for inla.group(substrate.grainsize, method = "quantile", n = 11) 53.562
Precision for inla.group(pca1, method = "quantile", n = 11)                53.564
Precision for inla.group(pca2, method = "quantile", n = 11)                54.779
Precision for space                                                    52.707
Precision for space_time                                                         56.223
Phi for space_time                                                                0.048
GroupRho for space_time                                                           0.772

Expected number of effective parameters(stdev): 658.29(3.85)
Number of equivalent replicates : 10.98

Deviance Information Criterion (DIC) ...............: 46260.72
Deviance Information Criterion (DIC, saturated) ....: 105353.22
Effective number of parameters .....................: 39.98

Watanabe-Akaike information criterion (WAIC) ...: 51091.12
Effective number of parameters .................: 3606.46

Marginal log-Likelihood:  -22329.74
Posterior marginals for the linear predictor and
 the fitted values are computed


  }


  ##########


  if (fishery_model) {

    # you need a stan installation on your system as well (outside of R)
    # install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
    
    require(cmdstanr)

      p$fishery_model = fishery_model( DS = "logistic_parameters", p=p, tag=p$areal_units_type )
      p$fishery_model$stancode = stan_initialize( stan_code=fishery_model( p=p, DS="stan_surplus_production" ) )

      str( p$fishery_model)

      p$fishery_model$stancode$compile()


      # res = fishery_model( p=p, DS="samples", tag=p$areal_units_type )  # to get samples

      if (0) {

        fit = fishery_model( p=p,   DS="fit", tag=p$areal_units_type )  # to get samples
        fit$summary(c("K", "r", "q", "qc"))

        print( fit, max_rows=30 )
        # fit$summary("K", "r", "q")

          # (penalized) maximum likelihood estimate (MLE)
        fit_mle =  p$fishery_model$stancode$optimize(data =p$fishery_model$standata, seed = 123)
        fit_mle$summary( c("K", "r", "q") )

        # fit$cmdstan_diagnose()
        # fit$cmdstan_summary()


        u = stan_extract( as_draws_df(fit_mle$draws() ) )

        mcmc_hist(fit$draws("K")) +
          vline_at(fit_mle$mle(), size = 1.5)

        # Variational Bayes
        fit_vb = p$fishery_model$stancode$variational( data =p$fishery_model$standata, seed = 123, output_samples = 4000)
        fit_vb$summary(c("K", "r", "q", "qc"))

        res = fishery_model(
          DS="logistic_model",
          p=p,
          tag=p$areal_units_type,
          fit = fit_vb
          # from here down are params for cmdstanr::sample()
        )

        u = stan_extract( as_draws_df(fit_vb$draws() ) )

        bayesplot_grid(
          mcmc_hist(fit$draws("K"), binwidth = 0.025),
          mcmc_hist(fit_vb$draws("K"), binwidth = 0.025),
          titles = c("Posterior distribution from MCMC", "Approximate posterior from VB")
        )

        color_scheme_set("gray")
        mcmc_dens(fit$draws("K"), facet_args = list(nrow = 3, labeller = ggplot2::label_parsed ) ) + facet_text(size = 14 )
        # mcmc_hist( fit$draws("K"))


      }



      fit = p$fishery_model$stancode$sample(
        data=p$fishery_model$standata,
        iter_warmup = 15000,
        iter_sampling = 8000,
        seed = 123,
        chains = 3,
        parallel_chains = 3,  # The maximum number of MCMC chains to run in parallel.
        max_treedepth = 16,
        adapt_delta = 0.995,
        refresh = 1000
      )

      fit$summary(c("K", "r", "q", "qc"))

      # save fit and get draws
      res = fishery_model( p=p, DS="logistic_model", tag=p$areal_units_type, fit=fit )       # from here down are params for cmdstanr::sample()

      # frequency density of key parameters
      fishery_model( DS="plot", vname="K", res=res )
      fishery_model( DS="plot", vname="r", res=res )
      fishery_model( DS="plot", vname="q", res=res, xrange=c(0.5, 2.5))
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

      # bosd
      plot.new()
      layout( matrix(c(1,2,3), 3, 1 ))
      par(mar = c(4.4, 4.4, 0.65, 0.75))
      for (i in 1:3) plot(density(res$mcmc$bosd[,i] ), main="")
      ( qs = apply(  res$mcmc$bosd[,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) )


      # bpsd
      plot.new()
      layout( matrix(c(1,2,3), 3, 1 ))
      par(mar = c(4.4, 4.4, 0.65, 0.75))
      for (i in 1:3) plot(density(res$mcmc$bpsd[,i] ), main="")
      ( qs = apply(  res$mcmc$bpsd[,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) )

      # rem_sd
      plot.new()
      layout( matrix(c(1,2,3), 3, 1 ))
      par(mar = c(4.4, 4.4, 0.65, 0.75))
      for (i in 1:3) plot(density(res$mcmc$rem_sd[,i] ), main="")
      ( qs = apply(  res$mcmc$rem_sd[,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) )


    # qc
      plot.new()
      layout( matrix(c(1,2,3), 3, 1 ))
      par(mar = c(4.4, 4.4, 0.65, 0.75))
      for (i in 1:3) plot(density(res$mcmc$qc[,i] ), main="")
      ( qs = apply(  res$mcmc$qc[,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) )

     # b0
      plot.new()
      layout( matrix(c(1,2,3), 3, 1 ))
      par(mar = c(4.4, 4.4, 0.65, 0.75))
      for (i in 1:3) plot(density(res$mcmc$b0[,i] ), main="")
      ( qs = apply(  res$mcmc$b0[,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) )


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
