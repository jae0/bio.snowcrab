

# Snow crab size structured model (via stan --- retired)

  require(bio.snowcrab)   # loadfunctions("bio.snowcrab") 

  year.assessment = 2021
  yrs = 1999:year.assessment
  spec_bio = bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 )
  

  loadfunctions("bio.snowcrab")

  pS = snowcrab_parameters( project_class="size_structured", yrs=yrs,
    selection = list( type = "number", biologicals=list( spec_bio=spec_bio ) ) 
  )



  M = snowcrab.db( p=pS, DS="set.complete" )  
  setDT(M)

  # must check why we have dups here
  i = which(duplicated(M[,1:10]))
  M = M[-i,]

  s_male = c("mi123.no", "mi4.no", "mi5.no", "mi6.no", "mi7.no", "mi8.no", "mi9.no", "mi10.no", "mi11.no", "mi12.no", "ma9.no", "ma10.no", "ma11.no", "ma12.no", "ma13.no" )
  i_male = c( 3:12, 9:13 )  # instar
  m_male = c( rep(0, 10), rep(1, 5) ) # maturity

  # convert to numbers
  males = M[ , c("uid",  s_male ) ]
  males[, s_male] = males[, s_male]  * M$sa
  m = data.table::melt( setDT(males), id.vars="uid", measure.vars=s_male, variable.factor=FALSE  )
  setnames(m, "value", "n")
  setnames(m, "variable", "stage")
  # plot(xtabs(n~stage, m)[s_male])

  m$instar = i_male[ match( m$stage, s_male ) ]
  m$maturity = m_male[ match( m$stage, s_male ) ]

  tot = m[, .(ntot=sum(n)), by=.(uid, instar) ]
  mat = m[ maturity==1, .(nmat=sum(n)),  by=.(uid, instar)]
  cnts = mat[ tot, on=.(uid, instar) ]
  cnts$fraction_mature = cnts$nmat / cnts$ntot 
 
  cnts = cnts[ M[,c("uid", "yr")], on=.(uid) ] 

  mx = dcast( cnts[,.(ntot=sum(ntot, na.rm=TRUE)), by=.( instar, yr) ], instar ~ yr, value.var="ntot" )
 
  mtm = dcast( cnts[,.(tm=sum(fraction_mature*ntot, na.rm=TRUE)), by=.( instar, yr) ], instar ~ yr, value.var="tm" )

  mm = as.matrix(mtm / mx)[ , 2:26]

  image( as.matrix( round( mm, 3) ) )

  plot(rowMeans(mm))



  s_female = c("fi1234.no", "fi5.no", "fi6.no", "fi7.no", "fi8.no", "fi9.no", "fi10.no", "fa7.no", "fa8.no", "fa9.no", "fa10.no" )

  females = M[ , c("uid",  s_female ) ]
  females[, s_female] = females[, s_female]  * M$sa
  f = data.table::melt( setDT(females), id.vars="uid", measure.vars=s_female, variable.factor=FALSE  )
  setnames(f, "value", "n")
  setnames(f, "variable", "stage")
  plot(xtabs(n~stage, f)[s_female])



  # some survey timestamps extend into January (e.g., 2020) force them to be part of the correct "survey year", i.e., "yr"
  i = which(lubridate::month(M$timestamp)==1)
  if (length(i) > 0) M$timestamp[i] = M$timestamp[i] - lubridate::duration(month=1)

  M$tiyr = lubridate::decimal_date(M$timestamp)
  M$dyear[ M$dyear > 1] = 0.99  # a survey year can run into the next year, cap the seasonal compenent at the calendar year for modellng 

  # reduce size
  M = M[ which( M$lon > p$corners$lon[1] & M$lon < p$corners$lon[2]  & M$lat > p$corners$lat[1] & M$lat < p$corners$lat[2] ), ]
  # levelplot(z.mean~plon+plat, data=M, aspect="iso")


















    if (areal_units) {
      # polygon structure:: create if not yet made
      # for (au in c("cfanorth", "cfasouth", "cfa4x", "cfaall" )) plot(polygon_managementareas( species="snowcrab", au))
      xydata = snowcrab.db( p=pN, DS="areal_units_input", redo=TRUE )
      xydata = snowcrab.db( p=pN, DS="areal_units_input" )
      sppoly = areal_units( p=pN, xydata=xydata[ which(xydata$yr %in% pN$yrs), ], redo=TRUE, verbose=TRUE )  # create constrained polygons with neighbourhood as an attribute
    
      sppoly=areal_units( p=pN )

      plot(sppoly["npts"])

     
      additional_features = snowcrab_mapping_features(pN)  # for mapping below

      figure_area_based_extraction_from_carstm(DS="temperature", year.assessment  )  # can only do done once we have an sppoly for snow crab
    
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
 
 
      # subset to positive definite data (for number and size)
      isubset =1:nrow(M)
      ipositive = unique( c( which( M$totno > 0), ip ) )
      if (grepl("hurdle", pN$carstm_model_label)) {
        isubset = ipositive
        hist(M$data_offset[ipositive])
      }



      # number
      fit = NULL; gc()
      fit = carstm_model( p=pN, data=M[ isubset, ], sppoly=sppoly, 
        posterior_simulations_to_retain="predictions", improve.hyperparam.estimates=TRUE
      )

      # mean size
      fit = NULL; gc()
      fit = carstm_model( p=pW, data=M[ isubset, ], sppoly = sppoly, 
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


        if ( number ) {
          p = pN
          outputdir = file.path( p$modeldir, p$carstm_model_label, "predicted.numerical.densities" )
          if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )
          fn_root_prefix = "Predicted_numerical_abundance"
          title= paste( snowcrab_filter_class, "Predicted numerical density (no./m^2) - persistent spatial effect"  )
        }
        if ( meansize) {
          p = pW
          outputdir = file.path( p$modeldir, p$carstm_model_label, "predicted.meansize" )
          if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )
          fn_root_prefix = "Predicted_meansize"
          title= paste( snowcrab_filter_class, "Predicted meansize (kg) - persistent spatial effect" ) 
        }
        if ( presence_absence ) {
          p = pH
          outputdir = file.path( p$modeldir, p$carstm_model_label, "predicted.presence_absence" )
          if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )
          fn_root_prefix = "Predicted_presence_absence"
          title= paste( snowcrab_filter_class, "Predicted habitat probability - persistent spatial effect")  
        }
  
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
    
        
        oeffdir = file.path(p$data_root, p$carstm_model_label, "figures")
        fn_root_prefix = "Predicted_numerical_abundance"
        carstm_plot_marginaleffects( p, oeffdir, fn_root_prefix ) 
    

        # maps of some of the results
        outputdir = file.path(p$data_root, p$carstm_model_label, "maps" )
        carstm_plot_map( p, outputdir, fn_root_prefix , additional_features, toplot="random_spatial", probs=c(0.025, 0.975) ) 
        carstm_plot_map( p, outputdir, fn_root_prefix , additional_features, toplot="predictions", probs=c(0.1, 0.9)) 

    }  # end spatiotemporal model


    # ----------------------

    assimilate_numbers_and_size = TRUE

    if (assimilate_numbers_and_size ) {

      # wgts_max = 1.1 # kg, hard upper limit
      # N_max = NULL
      # #  quantile( M$totno[ipositive]/M$data_offset[ipositive], probs=0.95, na.rm=TRUE )  
      
      # posterior sims 
      sims = carstm_posterior_simulations( pN=pN, pW=pW, pH=pH, pa_threshold=0.05, qmax=0.99 )
      sims = sims / 10^6 # 10^6 kg -> kt;; kt/km^2
      
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
    
      additional_features = snowcrab_mapping_featuresresresres(pN)  # for mapping below

      for (i in 1:length(pN$yrs) ){
        y = as.character( pN$yrs[i] )
        sppoly[,vn] = log10( B[,y]* 10^6 )
        outfilename = file.path( outputdir , paste( "biomass", y, "png", sep=".") )
        plt =  carstm_map(  
            sppoly=sppoly, 
            vn=vn,
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





##########

# this part is only relevent for R0 

fishery_model = FALSE

if (fishery_model) {

  # you need a stan installation on your system as well (outside of R), and the R-interface "cmdstanr":
  # install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
  
  require(cmdstanr)
  
  loadfunctions("bio.snowcrab")


  # choose:
  pN$fishery_model_label = "stan_surplus_production_2022_model_qc_uniform"
  # pN$fishery_model_label = "stan_surplus_production_2022_model_variation1_wider_qc_uniform"
  # pN$fishery_model_label = "stan_surplus_production_2022_model_variation1_wider_qc_normal"
  # pN$fishery_model_label = "stan_surplus_production_2022_model_qc_cauchy_wider"

  # pN$fishery_model_label = "stan_surplus_production_2022_model_qc_beta"  # qc_beta postive definite
  # pN$fishery_model_label = "stan_surplus_production_2022_model_qc_cauchy"
  # pN$fishery_model_label = "stan_surplus_production_2019_model"


  pN$fishery_model = fishery_model( DS = "logistic_parameters", p=pN, tag=pN$fishery_model_label )
  to_look = c("K", "r", "q", "qc", "log_lik" )

      if ( model_version=="framework_2019" ) {
        # this is to create results for reviewers in 2022 that wanted a comparison with previous methods .. can be deleted in future 
        # bring in unscaled abundance index
        a = fishery_model( DS="data_aggregated_timeseries", p=pN  )
        a$IOA[ !is.finite(a$IOA) ] = 0
        pN$fishery_model$fmdata$IOA = a$IOA
        to_look = c("K", "r", "q", "log_lik" )

      }  


  #  str( pN$fishery_model)

  pN$fishery_model$stancode = stan_initialize( stan_code=fishery_model( p=pN, DS=pN$fishery_model_label ) )
  pN$fishery_model$stancode$compile()
  

  fit = pN$fishery_model$stancode$sample(
    data=pN$fishery_model$fmdata,
    iter_warmup = 4000,
    iter_sampling = 4000,
    seed = 45678,
    chains = 3,
    parallel_chains = 3,  # The maximum number of MCMC chains to run in parallel.
    max_treedepth = 18,
    adapt_delta = 0.99,
    refresh = 1000
  )

  fit$summary(to_look)

  require(loo)
  waic(fit$draws("log_lik"))
  loo(fit$draws("log_lik"))
   
   
  # save fit and get draws
  res = fishery_model( p=pN, DS="logistic_model", tag=pN$fishery_model_label, fit=fit )       # from here down are params for cmdstanr::sample()

  if (0) {
    # reload saved fit and results
    load(pN$fishery_model$fnres)
    fit = aegis::read_write_fast(pN$fishery_model$fnfit)

  }



  # frequency density of key parameters
  fishery_model( DS="plot", vname="K", res=res )
  fishery_model( DS="plot", vname="r", res=res )
  fishery_model( DS="plot", vname="q", res=res, xrange=c(0.5, 2.5))
  fishery_model( DS="plot", vname="qc", res=res, xrange=c(-1, 1))
  fishery_model( DS="plot", vname="FMSY", res=res  )

  # timeseries
  fishery_model( DS="plot", type="timeseries", vname="biomass", res=res  )
  fishery_model( DS="plot", type="timeseries", vname="fishingmortality", res=res)

  # Harvest control rules
  fishery_model( DS="plot", type="hcr", vname="default", res=res  )

  # Summary table of mean values for inclusion in document
  
  ( qs = apply(  res$mcmc$K[,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) )  # carrying capactiy

  ( qs = apply(  res$mcmc$FMSY[,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) ) # FMSY


  biomass = as.data.table( fit$summary("B") )
  np = year.assessment+c(1:pN$fishery_model$fmdata$M)
  biomass$yr = rep( c(pN$yrs, np ), 3)
  nt = pN$fishery_model$fmdata$N +pN$fishery_model$fmdata$M
  biomass$region = c( rep("cfanorth", nt), rep("cfasouth", nt), rep("cfa4x", nt) )
  (biomass)

  NN = res$pN$fishery_model$fmdata$N

  # densities of biomass estimates for the year.assessment
  ( qs = apply(  res$mcmc$B[,NN,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) )

  # densities of biomass estimates for the previous year
  ( qs = apply(  res$mcmc$B[,NN-1,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) )

  # densities of F in assessment year
  ( qs = apply(  res$mcmc$F[,NN,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) )
  ( qs = apply(  res$mcmc$F[,NN,], 2, mean ) )

  # densities of F in previous year
  ( qs = apply(  res$mcmc$F[,NN-1,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) )
  ( qs = apply(  res$mcmc$F[,NN-1,], 2, mean ) )



  if (0) {
      # obsolete:

      fishery_model( DS="plot", vname="bosd", res=res  )
      fishery_model( DS="plot", vname="bpsd", res=res  )
      fishery_model( DS="plot", type="hcr", vname="simple", res=res  )
    
      # biomass.summary.table()

      fit = fishery_model( p=pN,   DS="fit", tag=pN$fishery_model_label )  # to load samples (results)
      fit$summary(c("K", "r", "q", "qc"))
      print( fit, max_rows=30 )
      fit$cmdstan_diagnose()
      fit$cmdstan_summary()

      # testing other samplers and optimizsers ... faster , good for debugging

      # (penalized) maximum likelihood estimate (MLE)
      fit_mle =  pN$fishery_model$stancode$optimize(data =pN$fishery_model$fmdata, seed = 123)
      fit_mle$summary( to_look )
      u = stan_extract( as_draws_df(fit_mle$draws() ) )

      # mcmc_hist(fit$draws("K")) + vline_at(fit_mle$mle(), size = 1.5)

      # # Variational Bayes
      fit_vb = pN$fishery_model$stancode$variational( data =pN$fishery_model$fmdata, seed = 123, output_samples = 4000)
      fit_vb$summary(to_look)
      fit_vb$cmdstan_diagnose()
      fit_vb$cmdstan_summary()


      u = stan_extract( as_draws_df(fit_vb$draws() ) )

      # bayesplot_grid(
      #   mcmc_hist(fit$draws("K"), binwidth = 0.025),
      #   mcmc_hist(fit_vb$draws("K"), binwidth = 0.025),
      #   titles = c("Posterior distribution from MCMC", "Approximate posterior from VB")
      # )

      # color_scheme_set("gray")
      # mcmc_dens(fit$draws("K"), facet_args = list(nrow = 3, labeller = ggplot2::label_parsed ) ) + facet_text(size = 14 )
      # mcmc_hist( fit$draws("K"))

      # obtain mcmc samples from vb solution
      res_vb = fishery_model(
        DS="logistic_model",
        p=pN,
        tag=pN$fishery_model_label,
        fit = fit_vb
      )

      names(res_vb$mcmc)

      # other diagnostics
      # fishery_model( DS="plot", type="diagnostic.errors", res=res )
      # fishery_model( DS="plot", type="diagnostic.phase", res=res  )

      NN = res$pN$fishery_model$fmdata$N

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

    }  # end skip

}  # end fishery model





# end
