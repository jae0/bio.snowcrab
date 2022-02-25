

# Snow crab --- Areal unit modelling of habitat  -- no reliance upon stmv fields


# TODO::: move plotting calls to self-contained functions:


NOTE :::: #######################################################################33
NOTE :::: For this to run, you must run three other projects that are dependencies 
NOTE ::::   (actually more if this is your first time through to get static fields - step 0)
NOTE ::::   0. aegis.bathymetry and aegis.substrate ( 01_ and 03_carstm.. ) .. static so only once
NOTE ::::   1. aegis.temperature  (01_temperature_data.R, 03_temperature_carstm_1970_present.R) 
NOTE ::::      greater time depth is required to get more spatial resolution 
NOTE ::::      using options for the shorter period: yrs=1999:year.assessment and 
NOTE ::::   2. aegis.survey (01_survey_data.R )
NOTE ::::   3. aegis.speciescomposition (01_speciescomposition_carstm_1999_to_present.R)
NOTE :::: #######################################################################33


# -------------------------------------------------
# Part 1 -- construct basic parameter list defining the main characteristics of the study
# require(aegis)

if( basic_parameters) {
  year.assessment = 2021
  require(bio.snowcrab)   # loadfunctions("bio.snowcrab") 


  # params for number
  pN = snowcrab_parameters(
    project_class="carstm",
    yrs=1999:year.assessment,   
    areal_units_type="tesselation",
    family="poisson",
    carstm_model_label = "1999_present",  # 1999_present is the default anything else and you are on your own
    # offset_shift=10^3,  # multiplier for data_offset and totno to bring data_offset closer in magnitude to 1
    selection = list(type = "number")
  )


  # params for mean size .. mostly the same as pN
  pW = snowcrab_parameters(
    project_class="carstm",
    yrs=1999:year.assessment,   
    areal_units_type="tesselation",
    carstm_model_label = "1999_present",  # 1999_present is the default anything else and you are on your own
    family =  "gaussian" ,  
    selection = list(type = "meansize")
  )

}



if (areal_units) {
  # polygon structure:: create if not yet made
  # for (au in c("cfanorth", "cfasouth", "cfa4x", "cfaall" )) plot(polygon_managementareas( species="snowcrab", au))
  xydata = snowcrab.db( p=pN, DS="areal_units_input", redo=TRUE )
  xydata = snowcrab.db( p=pN, DS="areal_units_input" )
  sppoly = areal_units( p=pN, xydata=xydata[ which(xydata$yr %in% pN$yrs), ], redo=TRUE, verbose=TRUE )  # create constrained polygons with neighbourhood as an attribute
   
 
  sppoly=areal_units( p=pN )
  plot(sppoly["npts"])
 
  M = snowcrab.db( p=pN, DS="carstm_inputs", sppoly=sppoly, redo=TRUE )  # will redo if not found
  # M = snowcrab.db( p=pN, DS="carstm_inputs", sppoly=sppoly  )  # will redo if not found
  
}
 
figure_area_based_extraction_from_carstm(DS="temperature", sppoly=sppoly)  # can only do done once we have an sppoly for snow crab


# ------------------------------------------------
# Part 2 -- spatiotemporal statistical model

if ( spatiotemporal_model ) {

  # total numbers
  sppoly=areal_units( p=pN )
  M = snowcrab.db( p=pN, DS="carstm_inputs", sppoly=sppoly  )  # will redo if not found
  
  io = which(M$tag=="observations")
  ip = which(M$tag=="predictions")

  mo = 1/ median(M$data_offset[io] ) 
  mo = 10^6
  M$data_offset[io] = M$data_offset[io] * mo  # ( number / km^2 ) * mo ==> this forces everyting to be expressed as  no. / m^2  .. including predictions (offset==1)
 

  hist(M$data_offset)

  fit = carstm_model( 
    p=pN, 
    data=M, 
    sppoly = sppoly, 
    posterior_simulations_to_retain="predictions" ,
    # control.inla = list( strategy="laplace"  ), 
    # redo_fit = FALSE,  # only to redo sims and extractions 
    redo_fit=TRUE,  
    theta = c( -0.813, 3.341, -2.447, -1.534, -2.768, 1.182, -0.274, -3.195, 1.095, -2.752, -0.812, 0.678 ),
    # theta=c( 1.226, 1.463, 1.528, -2.071, -3.106, -1.857, 2.405, -1.172, -0.003, 1.016  ), # to start optim from a solution close to the final in 2021 ...
    # debug = TRUE,
    # inla.mode="classic",
    num.threads="4:2"
  )

# List of hyperparameters: 
# 		theta[0] = [Log precision for time]
# 		theta[1] = [Rho_intern for time]
# 		theta[2] = [Log precision for cyclic]
# 		theta[3] = [Log precision for inla.group(t, method = "quantile", n = 11)]
# 		theta[4] = [Log precision for inla.group(z, method = "quantile", n = 11)]
# 		theta[5] = [Log precision for inla.group(pca1, method = "quantile", n = 11)]
# 		theta[6] = [Log precision for inla.group(pca2, method = "quantile", n = 11)]
# 		theta[7] = [Log precision for space]
# 		theta[8] = [Logit phi for space]
# 		theta[9] = [Log precision for space_time]
# 		theta[10] = [Logit phi for space_time]
# 		theta[11] = [Group rho_intern for space_time]

  
  
  if (0) {
    # extract results
    fit = carstm_model( p=pN, DS="carstm_modelled_fit",  sppoly = sppoly )  # extract currently saved model fit
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

  res = carstm_model( p=pN, DS="carstm_modelled_summary",  sppoly = sppoly ) # to load currently saved results


  additional_features = snowcrab_features_tmap(pN)  # for mapping below
  map_centre = c( (pN$lon0+pN$lon1)/2  , (pN$lat0+pN$lat1)/2  )
  map_zoom = 7.5

  if (quick_view) {

    vn=c( "random", "space", "combined" )
    vn=c( "random", "spacetime", "combined" )
    vn="predictions"  # numerical density (km^-2)
    tmatch="2015"

    carstm_map(  res=res, vn=vn, tmatch=tmatch, 
        sppoly = sppoly, 
        palette="-RdYlBu",
        plot_elements=c(  "compass", "scale_bar", "legend" ),
        additional_features=additional_features,
        tmap_zoom= c(map_centre, map_zoom),
        title =paste( vn, paste0(tmatch, collapse="-"), "no/m^2"  )
    )

  }  

  # map all :
   
  outputdir = file.path( pN$modeldir, pN$carstm_model_label, "predicted.numerical.densitites" )
  if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

  vn="predictions"
  toplot = carstm_results_unpack( res, vn )
  brks = pretty(  quantile(toplot[,,"mean"], probs=c(0,0.975), na.rm=TRUE )  )

  for (y in res$time ){
    tmatch = as.character(y)
    fn_root = paste("Predicted_numerical_abundance", paste0(tmatch, collapse="-"), sep="_")
    outfilename = file.path( outputdir, paste(fn_root, "png", sep=".") )

    tmout = carstm_map(  res=res, vn=vn, tmatch=tmatch,
      sppoly = sppoly, 
      breaks =brks,
      palette="-RdYlBu",
      plot_elements=c(   "compass", "scale_bar", "legend" ),
      additional_features=additional_features,
      map_mode="view",
      tmap_zoom= c(map_centre, map_zoom),
      title=paste("Predicted numerical density (no./m^2) ", paste0(tmatch, collapse="-") )
    )
  
    mapview::mapshot( tmap_leaflet(tmout), file=outfilename, vwidth = 1600, vheight = 1200 )  # very slow: consider 
    print(outfilename)
  
  }

  carstm_plotxy( res, vn=c( "res", "random", "time" ), 
    type="b", ylim=c(0,4), xlab="Year", ylab="No km^-2 x 10^4", h=0.5, v=1992   )

  carstm_plotxy( res, vn=c( "res", "random", "cyclic" ), 
    type="b", col="slategray", pch=19, lty=1, lwd=2.5, ylim=c(0.0, 2.5),
    xlab="Season", ylab="No km^-2 x 10^4" )

  carstm_plotxy( res, vn=c( "res", "random", "inla.group(t, method = \"quantile\", n = 11)" ), 
    type="b", col="slategray", pch=19, lty=1, lwd=2.5, ylim=c(0, 3) ,
    xlab="Bottom temperature (degrees Celcius)", ylab="No km^-2 x 10^4" )

  carstm_plotxy( res, vn=c( "res", "random", "inla.group(z, method = \"quantile\", n = 11)" ), 
    type="b", col="slategray", pch=19, lty=1, lwd=2.5, ylim=c(0, 5) ,
    xlab="Depth (m)", ylab="No km^-2 x 10^4" )

  # carstm_plotxy( res, vn=c( "res", "random", "inla.group(log.substrate.grainsize, method = \"quantile\", n = 11)" ), 
  #   type="b", col="slategray", pch=19, lty=1, lwd=2.5,
  #   xlab="Ln(grain size; mm)", ylab="No km^-2" )

  
  # -------------------------------------
  # mean size

  fit = carstm_model( 
    p=pW, 
    data=M, 
    sppoly = sppoly, 
    posterior_simulations_to_retain="predictions", 
    # control.inla = list( int.strategy='eb'  ), 
    # control.inla =#  list( strategy='adaptive' )
    # control.inla = list( strategy="laplace",  int.strategy='eb'  ), 
    # theta = c( 5.104, 9.052, 0.718, 9.880, 11.394, 8.104, 9.275, 2.938, 5.509, 4.319, 2.611 ),
    redo_fit=TRUE, # to start optim from a solution close to the final in 2021 ... 
    # redo_fit=FALSE, # to start optim from a solution close to the final in 2021 ... 
    # debug = TRUE,
    # inla.mode="classic", 
    num.threads="4:2"
  ) 
        

  if (0) {
    fit = carstm_model( p=pW, DS="carstm_modelled_fit",  sppoly = sppoly ) # to load currently saved results
  
    plot( fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
    plot( fit$marginals.hyperpar$"Phi for space_time", type="l")  # posterior distribution of phi nonspatial dominates
    plot( fit$marginals.hyperpar$"Precision for space_time", type="l")
    plot( fit$marginals.hyperpar$"Precision for setno", type="l")

  }


  res = carstm_model( p=pW, DS="carstm_modelled_summary",  sppoly = sppoly ) # to load currently saved results

  additional_features = snowcrab_features_tmap(pW)  # for mapping below
  map_centre = c( (p$lon0+p$lon1)/2  , (p$lat0+p$lat1)/2  )
  map_zoom = 7.5

  if (quick_view) {

    vn=c( "random", "space", "combined" )
    vn=c( "random", "spacetime", "combined" )
    vn="predictions"  # numerical density (km^-2)
    tmatch="2015"

    tmatch="2021"


    tmout = carstm_map(  res=res, vn=vn, tmatch=tmatch, 
        sppoly = sppoly,
        palette="-RdYlBu",
        plot_elements=c(  "compass", "scale_bar", "legend" ),
        additional_features=additional_features,
        tmap_zoom= c(map_centre, map_zoom),
        title =paste("Predicted  mean weight of individual (kg)", paste0(tmatch, collapse="-") )
    )

  }  

  # map all :
  outputdir = file.path( pW$modeldir, pW$carstm_model_label, "predicted.mean.weight" )
  if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )
  fn_root = paste("Predicted_mean_size", paste0(tmatch, collapse="-"), sep="_")
  fn = file.path( outputdir, paste(fn_root, "png", sep=".") )

  vn="predictions"
  toplot = carstm_results_unpack( res, vn )
  brks = pretty(  quantile(toplot[,,"mean"], probs=c(0,0.975), na.rm=TRUE )  )

  for (y in res$time ){
    tmatch = as.character(y)
    fn_root = paste("Predicted_numerical_abundance", paste0(tmatch, collapse="-"), sep="_")
    outfilename = file.path( outputdir, paste(fn_root, "png", sep=".") )

    tmout = carstm_map(  res=res, vn=vn, tmatch=tmatch,
      sppoly = sppoly, 
      breaks =brks,
      palette="-RdYlBu",
      plot_elements=c(   "compass", "scale_bar", "legend" ),
      additional_features=additional_features,
      map_mode="view",
      tmap_zoom= c(map_centre, map_zoom),
      title=paste("Predicted numerical density (no./km^2) ", paste0(tmatch, collapse="-") )
    )
    mapview::mapshot( tmap_leaflet(tmout), file=outfilename, vwidth = 1600, vheight = 1200 )  # very slow: consider 
    print(outfilename)
  }



  carstm_plotxy( res, vn=c( "res", "random", "time" ), 
    type="b", ylim=c(-0.04, 0.04), xlab="Year", ylab="Mean weight (kg)", h=0   )

  carstm_plotxy( res, vn=c( "res", "random", "cyclic" ), 
    type="b", col="slategray", pch=19, lty=1, lwd=2.5, ylim=c(-0.04, 0.04),
    xlab="Season", ylab="Mean weight (kg)" )

  carstm_plotxy( res, vn=c( "res", "random", "inla.group(t, method = \"quantile\", n = 11)" ), 
    type="b", col="slategray", pch=19, lty=1, lwd=2.5, ylim=c(-0.02, 0.02) ,
    xlab="Bottom temperature (degrees Celcius)", ylab="Mean weight (kg)" )

  carstm_plotxy( res, vn=c( "res", "random", "inla.group(z, method = \"quantile\", n = 11)" ), 
    type="b", col="slategray", pch=19, lty=1, lwd=2.5, ylim=c(-0.04, 0.04) ,
    xlab="Depth (m)", ylab="Mean weight (kg)" )


}



if (assimilate_numbers_and_size ) {

  p = pN
  p$pW = pW   # copy W params into pN

  # assimimilate no and meansize
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

  B = apply( bio, c(1,2), mean ) 
  
  brks = pretty( log10( quantile( B[], probs=c(0.05, 0.95) )* 10^6)  )
 
  additional_features = snowcrab_features_tmap(p)  # for mapping below
  map_centre = c( (p$lon0+p$lon1)/2  , (p$lat0+p$lat1)/2  )
  map_zoom = 7.5


  for (i in 1:length(p$yrs) ){
    y = as.character( p$yrs[i] )
    sppoly[,vn] = log10( B[,y]* 10^6 )
    outfilename = file.path( outputdir , paste( "biomass", y, "png", sep=".") )
    tmout =  carstm_map(  sppoly=sppoly, vn=vn,
        breaks=brks,
        additional_features=additional_features,
        title=paste("log_10( Predicted biomass density; kg/km^2 )", y ),
        palette="-RdYlBu",
        legend.text.size=3,
        plot_elements=c( "compass", "scale_bar", "legend" ),
        map_mode="view",
        tmap_zoom= c(map_centre, map_zoom) 
    )
    mapview::mapshot( tmap_leaflet(tmout), file=outfilename, vwidth = 1600, vheight = 1200 )  # very slow: consider 
    print(outfilename)
  }


}


##########


if (fishery_model) {

  # you need a stan installation on your system as well (outside of R), and the R-interface "cmdstanr":
  # install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
  
  require(cmdstanr)
  
  p=pN

  p$fishery_model = fishery_model( DS = "logistic_parameters", p=p, tag=p$areal_units_type )
  p$fishery_model$stancode = stan_initialize( 
    stan_code=fishery_model( p=p, DS="stan_surplus_production" ) 
  )
  #  str( p$fishery_model)
  p$fishery_model$stancode$compile()
  to_look = c("K", "r", "q", "qc" )

  fit = p$fishery_model$stancode$sample(
    data=p$fishery_model$standata,
    iter_warmup = 5000,
    iter_sampling = 10000,
    seed = 456,
    chains = 3,
    parallel_chains = 3,  # The maximum number of MCMC chains to run in parallel.
    max_treedepth = 16,
    adapt_delta = 0.99,
    refresh = 1000
  )

  fit$summary(to_look)

  # save fit and get draws
  res = fishery_model( p=p, DS="logistic_model", tag=p$areal_units_type, fit=fit )       # from here down are params for cmdstanr::sample()

  # load(p$fishery_model$fnres)
  # fit = readRDS(p$fishery_model$fnfit)

  # frequency density of key parameters
  fishery_model( DS="plot", vname="K", res=res )
  fishery_model( DS="plot", vname="r", res=res )
  fishery_model( DS="plot", vname="q", res=res, xrange=c(0.5, 2.5))
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
  np = year.assessment+c(1:p$fishery_model$standata$M)
  biomass$yr = rep( c(p$yrs, np ), 3)
  nt = p$fishery_model$standata$N +p$fishery_model$standata$M
  biomass$region = c( rep("cfanorth", nt), rep("cfasouth", nt), rep("cfa4x", nt) )
  (biomass)

  NN = res$p$fishery_model$standata$N

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

      fit = fishery_model( p=p,   DS="fit", tag=p$areal_units_type )  # to load samples (results)
      fit$summary(c("K", "r", "q", "qc"))
      print( fit, max_rows=30 )
      fit$cmdstan_diagnose()
      fit$cmdstan_summary()

      # testing other samplers and optimizsers ... faster , good for debugging

      # (penalized) maximum likelihood estimate (MLE)
      fit_mle =  p$fishery_model$stancode$optimize(data =p$fishery_model$standata, seed = 123)
      fit_mle$summary( to_look )
      u = stan_extract( as_draws_df(fit_mle$draws() ) )

      # mcmc_hist(fit$draws("K")) + vline_at(fit_mle$mle(), size = 1.5)

      # # Variational Bayes
      fit_vb = p$fishery_model$stancode$variational( data =p$fishery_model$standata, seed = 123, output_samples = 4000)
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
        p=p,
        tag=p$areal_units_type,
        fit = fit_vb
      )

      names(res_vb$mcmc)

      # other diagnostics
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

}





# end
