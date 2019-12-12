

# Snow crab --- Areal unit modelling of habitat  -- no reliance upon stmv fields


# -------------------------------------------------
# Part 1 -- construct basic parameter list defining the main characteristics of the study
  p = bio.snowcrab::snowcrab_carstm( DS="parameters", assessment.years=1999:2018 )

  # misc run params adjustments here:
  p$inla_num.threads = 4
  p$inla_blas.num.threads = 4


# -------------------------------------------------
# Part 2 -- polygon structure
  sppoly = areal_units( p=p )  # to reload
  plot(sppoly)
  spplot( sppoly, "au_sa_km2", main="AUID", sp.layout=p$coastLayout )
  if (0) {
    # create if not yet made
    MS = snowcrab.db( p=p, DS="biological_data"  )  # domain is  sse     # the underlying observations/data
    # ensure if polys exist and create if required
    for (au in c("cfanorth", "cfasouth", "cfa4x", "cfaall" )) plot(polygons_managementarea( species="snowcrab", au))
    sppoly = areal_units( p=p, areal_units_constraint=MS[, c("lon", "lat")], redo=TRUE )  # create constrained polygons with neighbourhood as an attribute
    coastLayout = aegis.coastline::coastline_layout( p=p, redo=TRUE )
    MS = NULL
  }


# -------------------------------------------------
# Part 3 -- create covariate field for bathymetry
# bathymetry -- ensure the data assimilation in bathymetry is first completed :: 01.bathymetry_data.R
# about 4.4 hrs to redo; 15 configs @ 0.5 hrs each

  pB = bathymetry_carstm( p=p, DS="parameters_override" )
  M = bathymetry.db( p=pB, DS="aggregated_data", redo=TRUE )
  M = bathymetry_carstm( p=pB, DS="carstm_inputs", redo=TRUE  ) # will redo if not found
  M = NULL; gc()
  res = carstm_model( p=pB, M='bathymetry_carstm( p=pB, DS="carstm_inputs" )', DS="redo", carstm_model_label="production"  ) # run model and obtain predictions

  if (0) {
    # to use a saved instance
    res = carstm_model( p=pB, DS="carstm_modelled", carstm_model_label="production" ) # run model and obtain predictions
    fit = carstm_model( p=pB, DS="carstm_modelled_fit", carstm_model_label="production" )  # extract currently saved model fit
    # maps of some of the results
    vn = paste(pB$variabletomodel, "predicted", sep=".")
    carstm_plot( p=pB, res=res, vn=vn )

    vn = paste(pB$variabletomodel, "predicted_se", sep=".")
    carstm_plot( p=pB, res=res, vn=vn )

    vn = paste(pB$variabletomodel, "random_auid_nonspatial", sep=".")
    carstm_plot( p=pB, res=res, vn=vn )

    vn = paste(pB$variabletomodel, "random_auid_spatial", sep=".")
    carstm_plot( p=pB, res=res, vn=vn )

    plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE, single=TRUE )

    # Time used:
    #     Pre = 4.71, Running = 10928, Post = 4.51, Total = 10938
    # Fixed effects:
    #              mean    sd 0.025quant 0.5quant 0.975quant  mode kld
    # (Intercept) 7.883 0.001      7.882    7.883      7.885 7.883   0

    # Random effects:
    #   Name	  Model
    #     auid BYM2 model

    # Model hyperparameters:
    #                                             mean     sd 0.025quant 0.5quant 0.975quant    mode
    # Precision for the lognormal observations 752.050  3.198    745.879  752.003    758.450 751.849
    # Precision for auid                     286.200 23.004    239.367  287.143    329.336 291.721
    # Phi for auid                             0.979  0.021      0.922    0.985      0.999   0.996

    # Expected number of effective parameters(stdev): 698.31(1.78)
    # Number of equivalent replicates : 161.11

    # Deviance Information Criterion (DIC) ...............: 1348750.89
    # Deviance Information Criterion (DIC, saturated) ....: 113171.52
    # Effective number of parameters .....................: 698.85

    # Marginal log-Likelihood:  -674832.61
    # Posterior marginals for the linear predictor and
    #  the fitted values are computed

  }




# -------------------------------------------------
# Part 4 -- create covariate field for  substrate
# ensure the data assimilation in substrate is first completed :: 01.substrate_data.R
# 25 configs @ 2 hrs each, total time 32 hrs
  pS = substrate_carstm(p=p, DS="parameters_override" )
  M = substrate.db( p=pS, DS="aggregated_data", redo=TRUE )  # used for data matching/lookup in other aegis projects that use substrate
  M = substrate_carstm( p=pS, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  M = NULL; gc()
  res = carstm_model( p=pS, M='substrate_carstm( p=pS, DS="carstm_inputs")', DS="redo", carstm_model_label="production"  )  # run model and obtain predictions

  if(0) {
    res = carstm_model( p=pS, DS="carstm_modelled", carstm_model_label="production"   ) # run model and obtain predictions
    fit = carstm_model( p=pS, DS="carstm_modelled_fit", carstm_model_label="production" )  # extract currently saved model fit
    summary(fit)
    # Model hyperparameters:
    #                                           mean    sd 0.025quant 0.5quant 0.975quant  mode
    # Precision for the lognormal observations 1.405 0.006      1.392    1.405      1.417 1.405
    # Precision for zi                         4.318 2.467      1.105    3.822     10.495 2.771
    # Precision for auid                       0.840 0.124      0.652    0.820      1.134 0.770
    # Phi for auid                             0.959 0.034      0.867    0.969      0.995 0.986

    # Expected number of effective parameters(stdev): 191.31(0.21)
    # Number of equivalent replicates : 502.93

    # Deviance Information Criterion (DIC) ...............: 41922.78
    # Deviance Information Criterion (DIC, saturated) ....: 96422.26
    # Effective number of parameters .....................: 192.34

    # Marginal log-Likelihood:  -21357.76

    vn = paste(pS$variabletomodel, "predicted", sep=".")
    carstm_plot( p=pS, res=res, vn=vn ) # maps of some of the results

    vn = paste(pS$variabletomodel, "predicted_se", sep=".")
    carstm_plot( p=pS, res=res, vn=vn )

    vn = paste(pS$variabletomodel, "random_auid_nonspatial", sep=".")
    carstm_plot( p=pS, res=res, vn=vn )

    vn = paste(pS$variabletomodel, "random_auid_spatial", sep=".")
    carstm_plot( p=pS, res=res, vn=vn )



    plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE, single=TRUE )

    #------------
    # 6567.727 sec = 1.8hrs

    pS$carstm_model_label = "production_inla.group_quantile_25"
    # DIC:
    #   Mean of Deviance ................. 33840.1
    #   Deviance at Mean ................. 33636.8
    #   Effective number of parameters ... 203.318
    #   DIC .............................. 34043.4
    # DIC (Saturated):
    #   Mean of Deviance ................. 95897.9
    #   Deviance at Mean ................. 95694.6
    #   Effective number of parameters ... 203.318
    #   DIC .............................. 96101.3




    pS$carstm_model_label = "production_inla.group_quantile_20"
    # DIC:
    # 	Mean of Deviance ................. 33995.4
    # 	Deviance at Mean ................. 33794.6
    # 	Effective number of parameters ... 200.881
    # 	DIC .............................. 34196.3
    # DIC (Saturated):
    # 	Mean of Deviance ................. 96187.3
    # 	Deviance at Mean ................. 95986.4
    # 	Effective number of parameters ... 200.881
    # 	DIC .............................. 96388.2


    pS$carstm_modelcall = paste('
      inla(
        formula =', pS$variabletomodel, ' ~ 1
          + f( inla.group(z, method="quantile", n=20) ,  model="rw2", scale.model=TRUE, hyper=H$rw2)
          + f(auid, model="bym2", graph=sppoly@nb, scale.model=TRUE, constr=TRUE, hyper=H$bym2),
        family = "lognormal",
        data= M,
        control.compute=list(dic=TRUE, config=TRUE),
        control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
        control.predictor=list(compute=FALSE, link=1 ),
        control.fixed=H$fixed,  # priors for fixed effects, generic is ok
        # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),  # extra work to get tails
        # control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
        # control.mode = list( restart=TRUE, result=RES ), # restart from previous estimates
        # control.inla = list(cmin = 0 ),
        # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
        # control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
        # control.inla = list(h=1e-3, tolerance=1e-9, cmin=0), # restart=3), # restart a few times in case posteriors are poorly defined
        # control.mode = list( restart=TRUE, result=RES ), # restart from previous estimates
        verbose=TRUE
      ) ' )
    res = carstm_model( p=pS, M='substrate_carstm( p=pS, DS="carstm_inputs")', DS="redo", carstm_model_label=pS$carstm_model_label  )  # run model and obtain predictions
    fit = carstm_model( p=pS, DS="carstm_modelled_fit", carstm_model_label=pS$carstm_model_label )  # extract currently saved model fit

  }




# -------------------------------------------------
# Part 5 -- create covariate field for temperature
# ensure the data assimilation in temperature is first completed :: 01.temperature_data.R
# total: 30 min, 80 configs .. fast
  pT = temperature_carstm(p=p, DS="parameters_override" )
  M = temperature.db( p=pT, DS="aggregated_data", redo=TRUE )  #  used for data matching/lookup in other aegis projects that use temperature
  M = temperature_carstm( p=pT, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  M = NULL; gc()
  res = carstm_model( p=pT, M='temperature_carstm( p=pT, DS="carstm_inputs")', DS="redo", carstm_model_label="production"  ) # run model and obtain predictions
  if (0) {
    res = carstm_model( p=pT, DS="carstm_modelled", carstm_model_label="production" ) # run model and obtain predictions
    fit = carstm_model(  p=pT, DS="carstm_modelled_fit", carstm_model_label="production" )  # extract currently saved model fit
    summary(fit)

# Time used:
#     Pre = 3.23, Running = 3680, Post = 16.4, Total = 3700
# Fixed effects:
#              mean    sd 0.025quant 0.5quant 0.975quant  mode kld
# (Intercept) 4.454 0.388      3.682    4.453      5.228 4.451   0

# Random effects:
#   Name	  Model
#     year_factor AR1 model
#    dyri AR1 model
#    inla.group(z, method = "quantile", n = 25) RW2 model
#    auid BYM2 model

# Model hyperparameters:
#                                                           mean    sd 0.025quant 0.5quant 0.975quant  mode
# Precision for the Gaussian observations                  0.429 0.003      0.422    0.429      0.436 0.429
# Precision for year_factor                                2.878 0.841      1.549    2.772      4.821 2.571
# Rho for year_factor                                      0.386 0.155      0.059    0.395      0.663 0.411
# Precision for dyri                                       3.092 1.126      1.383    2.934      5.749 2.622
# Rho for dyri                                             0.601 0.147      0.268    0.617      0.838 0.653
# Precision for inla.group(z, method = "quantile", n = 25) 0.128 0.028      0.079    0.126      0.190 0.122
# Precision for auid                                       0.410 0.033      0.352    0.408      0.482 0.401
# Phi for auid                                             0.996 0.005      0.983    0.998      1.000 1.000
# GroupRho for auid                                        0.723 0.025      0.670    0.725      0.767 0.731

# Expected number of effective parameters(stdev): 1563.56(31.92)
# Number of equivalent replicates : 22.46

# Deviance Information Criterion (DIC) ...............: 130925.13
# Deviance Information Criterion (DIC, saturated) ....: 36676.07
# Effective number of parameters .....................: 1565.21

# Marginal log-Likelihood:  -64612.59
# Posterior marginals for the linear predictor and
#  the fitted values are computed


  # Time used:
  #     Pre = 2.67, Running = 559, Post = 2.69, Total = 564
  # Fixed effects:
  #             mean    sd 0.025quant 0.5quant 0.975quant  mode kld
  # (Intercept) 5.127 0.309       4.52    5.127      5.734 5.127   0

  # Random effects:
  #   Name	  Model
  #     year_factor AR1 model
  #   dyri RW2 model
  #   zi RW2 model
  #   auid BYM2 model

  # Model hyperparameters:
  #                                         mean    sd 0.025quant 0.5quant 0.975quant  mode
  # Precision for the Gaussian observations 0.421 0.000      0.421    0.421   4.52e-01 0.421
  # Precision for year_factor               1.325 0.057      1.181    1.308   2.29e+05 1.344
  # Rho for year_factor                     0.458 0.022      0.398    0.452   1.00e+00 0.466
  # Precision for dyri                      1.131 0.061      1.005    1.121   1.30e+00 1.085
  # Precision for zi                        0.001 0.000      0.001    0.001   1.06e-01 0.001
  # Precision for auid                    0.220 0.001      0.214    0.220   2.23e-01 0.220
  # Phi for auid                          0.993 0.000      0.992    0.993   1.00e+00 0.994

  # Expected number of effective parameters(stdev): 215.49(0.00)
  # Number of equivalent replicates : 146.64

  # Deviance Information Criterion (DIC) ...............: 128672.37
  # Deviance Information Criterion (DIC, saturated) ....: 43266.47
  # Effective number of parameters .....................: 215.48

  # Marginal log-Likelihood:  -64703.07
  # Posterior marginals for the linear predictor and
  # the fitted values are computed



  # Time used:
  #     Pre = 1.5, Running = 848, Post = 2.71, Total = 852
  # Fixed effects:
  #             mean    sd 0.025quant 0.5quant 0.975quant  mode kld
  # (Intercept) 5.101 0.206      4.698    5.101      5.504 5.101   0

  # Random effects:
  #   Name	  Model
  #     year_factor AR1 model
  #   dyri RW2 model
  #   zi RW2 model
  #   auid BYM2 model

  # Model hyperparameters:
  #                                         mean    sd 0.025quant 0.5quant 0.975quant  mode
  # Precision for the Gaussian observations 0.355 0.003      0.349    0.355      0.361 0.355
  # Precision for year_factor               2.417 0.850      1.036    2.335      4.306 2.137
  # Rho for year_factor                     0.440 0.176      0.090    0.444      0.761 0.442
  # Precision for dyri                      0.556 0.305      0.167    0.492      1.331 0.375
  # Precision for zi                        0.001 0.000      0.000    0.001      0.001 0.001
  # Precision for auid                    0.328 0.027      0.280    0.325      0.387 0.319
  # Phi for auid                          0.992 0.010      0.966    0.995      1.000 1.000
  # GroupRho for auid                     0.843 0.015      0.810    0.844      0.868 0.848

  # Expected number of effective parameters(stdev): 1216.66(0.00)
  # Number of equivalent replicates : 25.97

  # Deviance Information Criterion (DIC) ...............: 123610.43
  # Deviance Information Criterion (DIC, saturated) ....: 32822.21
  # Effective number of parameters .....................: 1216.63

  # Marginal log-Likelihood:  -60756.52



    vn = paste(pT$variabletomodel, "predicted", sep=".")
    carstm_plot( p=pT, res=res, vn=vn, time_match=list(year="2000", dyear="0.85" ) )       # maps of some of the results

    vn = paste(pT$variabletomodel, "predicted_se", sep=".")
    carstm_plot( p=pT, res=res, vn=vn, time_match=list(year="2000", dyear="0.85" ) )       # maps of some of the results

    vn = paste(pT$variabletomodel, "random_auid_nonspatial", sep=".")
    carstm_plot( p=pT, res=res, vn=vn, time_match=list(year="2000"  ) )       # maps of some of the results
    vn = paste(pT$variabletomodel, "random_auid_spatial", sep=".")
    carstm_plot( p=pT, res=res, vn=vn, time_match=list(year="2000"  ) )       # maps of some of the results



    plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE, single=TRUE )


        pT$carstm_model_label = "production"
        pT$carstm_modelcall = paste('
          inla(
            formula = ', pT$variabletomodel, ' ~ 1
              + f( year_factor, model="ar1", hyper=H$ar1 )
              + f( dyri, model="ar1", scale.model=TRUE, hyper=H$ar1 )
              + f( inla.group( z, method="quantile", n=25 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
              + f( auid, model="bym2", graph=sppoly@nb, group=year_factor, scale.model=TRUE, constr=TRUE, hyper=H$bym2),
            family = "normal",
            data= M,
            control.compute=list(dic=TRUE, config=TRUE),
            control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
            control.predictor=list(compute=FALSE, link=1 ),
            control.fixed=H$fixed,  # priors for fixed effects, generic is ok
            # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
            # control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
            # control.mode = list( restart=TRUE, result=RES ), # restart from previous estimates
            # control.inla = list(cmin = 0 ),
            # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
            # control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
            # control.inla = list(h=1e-3, tolerance=1e-9, cmin=0), # restart=3), # restart a few times in case posteriors are poorly defined
            # control.mode = list( restart=TRUE, result=RES ), # restart from previous estimates
            verbose=TRUE
          ) ' )

        #  + f(tiyr2, model="seasonal", season.length=10 )
        #  + f(dyear, model="ar1", hyper=H$ar1 )
        #  + f(seasonal, model="seasonal", season.length=', pT$n.season, ', scale.model=TRUE )  # using seasonal effect is not recommended as it is not smoothed well .. rw2 is better



  }





# -------------------------------------------------
# Part 6 -- create covariate field for species composition 1
# ensure that survey data is assimilated : bio.snowcrab::01snowcb_data.R, aegis.survey::01.surveys.data.R , etc.
# 30 min, 150 configs
  pPC1 = speciescomposition_carstm(p=p, DS="parameters_override", varnametomodel="pca1" )
  M = speciescomposition_carstm( p=pPC1, DS="carstm_inputs", varnametomodel="pca1", redo=TRUE )  # will redo if not found
  M = NULL; gc()
  res = carstm_model( p=pPC1, M='speciescomposition_carstm( p=p, DS="carstm_inputs", varnametomodel="pca1" )', DS="redo", carstm_model_label="production"   ) # run model and obtain predictions
  if (0) {
    res = carstm_model( p=pPC1, DS="carstm_modelled", carstm_model_label="production"  ) # run model and obtain predictions
    fit = carstm_model( p=pPC1, DS="carstm_modelled_fit", carstm_model_label="production" )  # extract currently saved model fit
    summary(fit)
    vn = paste(pPC1$variabletomodel, "predicted", sep=".")
    carstm_plot( p=pPC1, res=res, vn=vn, time_match=list(year="2000" ), dyear="0.85" ) # maps of some of the results
    vn = paste(pPC1$variabletomodel, "random_auid_nonspatial", sep=".")
    carstm_plot( p=pPC1, res=res, vn=vn, time_match=list(year="2000" ), dyear="0.85" )       # maps of some of the results , dyear="0.85"
    vn = paste(pPC1$variabletomodel, "random_auid_spatial", sep=".")
    carstm_plot( p=pPC1, res=res, vn=vn, time_match=list(year="2000" ), dyear="0.85" )       # maps of some of the results , dyear="0.85"
  }





# -------------------------------------------------
# Part 7 -- create covariate field for species composition 2
# ensure that survey data is assimilated : bio.snowcrab::01snowcb_data.R, aegis.survey::01.surveys.data.R ,
# etc.30 min, 30 min
  pPC2 = speciescomposition_carstm(p=p, DS="parameters_override", varnametomodel="pca2" )
  M = speciescomposition_carstm( p=pPC2, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  M = NULL; gc()
  res = carstm_model( p=pPC2, M='speciescomposition_carstm( p=p, DS="carstm_inputs", varnametomodel="pca2" )', DS="redo" , carstm_model_label="production"  ) # run model and obtain predictions
  if (0) {
    res = carstm_model( p=pPC2, DS="carstm_modelled", carstm_model_label="production"   ) # run model and obtain predictions
    fit = carstm_model( p=pPC2, DS="carstm_modelled_fit", carstm_model_label="production"  )  # extract currently saved model fit
    summary(fit)
    vn = paste(pPC2$variabletomodel, "predicted", sep=".")
    carstm_plot( p=pPC2, res=res, vn=vn, time_match=list(year="2000" ), dyear="0.85" )       # maps of some of the results
    vn = paste(pPC2$variabletomodel, "random_auid_nonspatial", sep=".")
    carstm_plot( p=pPC2, res=res, vn=vn, time_match=list(year="2000" ), dyear="0.85" )       # maps of some of the results , dyear="0.85"
    vn = paste(pPC2$variabletomodel, "random_auid_spatial", sep=".")
    carstm_plot( p=pPC2, res=res, vn=vn, time_match=list(year="2000" ), dyear="0.85" )       # maps of some of the results , dyear="0.85"
  }



# finished covariates ... move onto abundance index estimation



# -------------------------------------------------
# Part 8 -- Snow crab anbundance -- main mode used for production
  M = snowcrab_carstm( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  M = NULL; gc()
  res = carstm_model( p=p, M='snowcrab_carstm( p=p, DS="carstm_inputs" )' ) # 151 configs and long optim .. 19 hrs
  m = get("inla.models", INLA:::inla.get.inlaEnv())
	m$latent$rw2$min.diff = NULL
	assign("inla.models", m, INLA:::inla.get.inlaEnv())

  if (0) {

    # extract results
    res = carstm_model( p=p, DS="carstm_modelled" ) # to load currently saved res
    fit = carstm_model( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit
    summary(fit)
    vn = paste(p$variabletomodel, "predicted", sep=".")
    carstm_plot( p=p, res=res, vn=vn, time_match=list(year="2000" ) )     # maps of some of the results

    plot(fit)
    plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE, single=TRUE )
    s = summary(fit)
    s$dic$dic
# [1] 47779.408  # poisson

    s$dic$p.eff
# [1] 1962.805 # poisson

    # maps of some of the results

    vn = paste(p$variabletomodel, "random_sample_iid", sep=".")
    if (exists(vn, res)) carstm_plot( p=p, res=res, vn=vn, time_match=list(year="1950", dyear="0.05") )

    vn = paste(p$variabletomodel, "random_auid_nonspatial", sep=".")
    if (exists(vn, res)) {
      res_dim = dim( res[[vn]] )
      if (res_dim == 1 ) time_match = NULL
      if (res_dim == 2 ) time_match = list(year="2000")
      if (res_dim == 3 ) time_match = list(year="2000", dyear="0.85" )
      carstm_plot( p=p, res=res, vn=vn, time_match=time_match )
    }

    vn = paste(p$variabletomodel, "random_auid_spatial", sep=".")
    if (exists(vn, res)) {
      res_dim = dim( res[[vn]] )
      if (res_dim == 1 ) time_match = NULL
      if (res_dim == 2 ) time_match = list(year="2000")
      if (res_dim == 3 ) time_match = list(year="2000", dyear="0.85" )
      carstm_plot( p=p, res=res, vn=vn, time_match=time_match )
    }

  }


  snowcrab_abundance_index( p=p, operation="compute", carstm_model_label=p$carstm_model_label )

  RES = snowcrab_abundance_index(p=p, operation="load_timeseries", carstm_model_label=p$carstm_model_label  )

  bio = snowcrab_abundance_index(p=p, operation="load_spacetime_biomass", carstm_model_label=p$carstm_model_label  )
  num = snowcrab_abundance_index(p=p, operation="load_spacetime_number", carstm_model_label=p$carstm_model_label  )
  wt = snowcrab_abundance_index(p=p, operation="load_spacetime_weights", carstm_model_label=p$carstm_model_label  )


  plot( cfaall ~ yrs, data=RES, lty=1, lwd=2.5, col="red", type="b")
  plot( cfasouth ~ yrs, data=RES, lty=1, lwd=2.5, col="red", type="b")
  plot( cfanorth ~ yrs, data=RES, lty=1, lwd=2.5, col="red", type="b")
  plot( cfa4x ~ yrs, data=RES, lty=1, lwd=2.5, col="red", type="b")

  p$boundingbox = list( xlim=p$corners$lon, ylim=p$corners$lat) # bounding box for plots using spplot

  p$coastLayout = aegis.coastline::coastline_layout(p=p)
  p$mypalette=RColorBrewer::brewer.pal(9, "YlOrRd")
    sppoly = areal_units( p=p )  # to reload


  # map it ..mean density
  vn = "pred"
  sppoly@data[,vn] = bio[,"2017"]
  brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
  spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )


  plot( fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
  plot( fit$marginals.hyperpar$"Phi for auid", type="l")  # posterior distribution of phi nonspatial dominates
  plot( fit$marginals.hyperpar$"Precision for auid", type="l")
  plot( fit$marginals.hyperpar$"Precision for setno", type="l")





# end


# other model variations for testing follows:

if (0) {

  # 7 hrs
	# Mean of Deviance ................. 45495.7
	# Deviance at Mean ................. 43774.3
	# Effective number of parameters ... 1721.47
	# DIC .............................. 47217.2

  p$carstm_model_label = "zeroinflatedpoisson0"
  p$carstm_modelcall = paste(
    'inla( formula =', p$variabletomodel,
    ' ~ 1
      + offset( log(data_offset))
      + f(year_factor, model="ar1", hyper=H$ar1 )
      + f(dyri, model="rw2", scale.model=TRUE, hyper=H$rw2 )
      + f(ti, model="rw2", scale.model=TRUE, hyper=H$rw2)
      + f(zi, model="rw2", scale.model=TRUE, hyper=H$rw2)
      + f(gsi, model="rw2", scale.model=TRUE, hyper=H$rw2)
      + f(pca1i, model="rw2", scale.model=TRUE, hyper=H$rw2)
      + f(pca2i, model="rw2", scale.model=TRUE, hyper=H$rw2)
      + f(auid, model="bym2", graph=sppoly@nb, group=year_factor, scale.model=TRUE, constr=TRUE, hyper=H$bym2),
      family = "zeroinflatedpoisson0",
      data= M,
      control.compute=list(dic=TRUE, config=TRUE),
      control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
      control.predictor=list(compute=FALSE, link=1 ),
      control.fixed=H$fixed,  # priors for fixed effects, generic is ok
      # control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
      # control.mode = list( restart=TRUE, result=RES ), # restart from previous estimates
      # control.inla = list(cmin = 0 ),
      # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
      # control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
      # control.inla = list(h=1e-3, tolerance=1e-9, cmin=0), # restart=3), # restart a few times in case posteriors are poorly defined
      # control.mode = list( restart=TRUE, result=RES ), # restart from previous estimates
    verbose=TRUE
    )'
  )


}

# -------------------------------------
# simpler models:




# -------------------------------------
# model 1 - simple glm
M = snowcrab_carstm( p=p, DS="carstm_inputs"  )  # will redo if not found
sppoly = areal_units( p=p )

fit = glm(
  formula = totno ~ 1 + offset( log( data_offset) ) + AUID + year_factor,
  family = "poisson", # "zeroinflatedpoisson0",
  data= M[ which(M$tag=="observations"), ]
)

s = summary(fit)
AIC(fit)  # 77326

# reformat predictions into matrix form
ii = which(
  M$tag=="predictions" &
  M$AUID %in% M[ which(M$tag=="observations"), "AUID"] &
  M$year_factor %in% M[ which(M$tag=="observations"), "year_factor"]
)

preds = predict( fit, newdata=M[ii,], type="response", na.action=na.omit, se.fit=TRUE )  # no/km2

out = reformat_to_array(
  input = preds$fit,
  matchfrom = list( AUID=as.character(M$AUID[ii]),  year_factor=as.character(M$year[ii]) ),
  matchto   = list( AUID=as.character(sppoly$AUID), year_factor=as.character(factor(p$yrs)) )
)
# out[ out>1e10] = NA
# convert numbers/km to biomass/auid (kg)..


# / 10^6  # 10^6 kg -> kt .. kg/km * km
RES$poisson_glm = list(
  yrs = p$yrs,
  cfaall    = colSums( out * weight_year * sppoly$au_sa_km2, na.rm=TRUE ) / 10^6 ,
  cfanorth  = colSums( out * weight_year * sppoly$cfanorth_surfacearea, na.rm=TRUE )/ 10^6 ,
  cfasouth  = colSums( out * weight_year * sppoly$cfasouth_surfacearea, na.rm=TRUE )/ 10^6 ,
  cfa23     = colSums( out * weight_year * sppoly$cfa23_surfacearea, na.rm=TRUE )/ 10^6 ,
  cfa24     = colSums( out * weight_year * sppoly$cfa24_surfacearea, na.rm=TRUE )/ 10^6 ,
  cfa4x     = colSums( out * weight_year * sppoly$cfa4x_surfacearea, na.rm=TRUE )/ 10^6
)


plot( cfaall ~ yrs, data=RES$poisson_glm, lty=1, lwd=2.5, col="blue", type="b")
plot( cfanorth ~ yrs, data=RES$poisson_glm, lty=1, lwd=2.5, col="green", type="b")
plot( cfasouth ~ yrs, data=RES$poisson_glm, lty=1, lwd=2.5, col="green", type="b")
plot( cfa4x ~ yrs, data=RES$poisson_glm, lty=1, lwd=2.5, col="green", type="b")



# map it
vn = "pred"
yr = "2018"
sppoly@data[,vn] = out[,yr] * weight_year[,yr]  # biomass density
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )





# -------------------------------------
# simple gam
require(mgcv)
fit = gam(
  formula = Y ~ 1 + offset( log( data_offset) ) + AUID + yr_factor + s(auid, year, bs="ts"),
  family = "poisson", # "zeroinflatedpoisson0",
  data= M[ which(M$tag=="observations"), ]
)

s = summary(fit)
AIC(fit)  # 76752

# reformat predictions into matrix form
ii = which(
  M$tag=="predictions" &
  M$AUID %in% M[ which(M$tag=="observations"), "AUID"] &
  M$yr_factor %in% M[ which(M$tag=="observations"), "yr_factor"]
)

preds = predict( fit, newdata=M[ii,], type="response", na.action=na.omit, se.fit=TRUE )  # no/km2

out = reformat_to_array(
  input = preds$fit,
  matchfrom = list( AUID=M$AUID[ii], yr_factor=M$yr_factor[ii]),
  matchto   = list( AUID=sppoly$AUID, yr_factor=factor(p$yrs) )
)

RES$poisson_gam = list(
  yrs = p$yrs,
  cfaall    = colSums( out * weight_year * sppoly$au_sa_km2, na.rm=TRUE ) / 10^6 ,
  cfanorth  = colSums( out * weight_year * sppoly$cfanorth_surfacearea, na.rm=TRUE )/ 10^6 ,
  cfasouth  = colSums( out * weight_year * sppoly$cfasouth_surfacearea, na.rm=TRUE )/ 10^6 ,
  cfa23     = colSums( out * weight_year * sppoly$cfa23_surfacearea, na.rm=TRUE )/ 10^6 ,
  cfa24     = colSums( out * weight_year * sppoly$cfa24_surfacearea, na.rm=TRUE )/ 10^6 ,
  cfa4x     = colSums( out * weight_year * sppoly$cfa4x_surfacearea, na.rm=TRUE )/ 10^6
)


plot( poisson_gam ~ yr, data=RES$poisson_gam, lty=1, lwd=2.5, col="blue", type="b")
plot( poisson_gam_cfanorth ~ yr, data=RES$poisson_gam, lty=1, lwd=2.5, col="green", type="b")
plot( poisson_gam_cfasouth ~ yr, data=RES$poisson_gam, lty=1, lwd=2.5, col="green", type="b")
plot( poisson_gam_cfa4x ~ yr, data=RES$poisson_gam, lty=1, lwd=2.5, col="green", type="b")

# map it
vn = "pred"
yr = "2018"
sppoly@data[,vn] = out[,yr] * weight_year[,yr]  # biomass density
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )





# - -----------------------------
# simple with default priors
fit = inla(
  formula = Y ~ 1 + offset( log( data_offset) ) + AUID + yr_factor ,
  family = "poisson", # "zeroinflatedpoisson0",
  data= M,
  control.compute=list(dic=TRUE, config=TRUE),
  control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
  control.predictor=list(compute=FALSE, link=1 ),
  control.fixed=H$fixed,  # priors for fixed effects, generic is ok
  control.inla=list(int.strategy="eb") ,# to get empirical Bayes results much faster.
 # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
  # num.threads=4,
  # blas.num.threads=4,
  verbose=TRUE
)

fn_test ="~/tmp/snowcrab_20km_bym.Rdata"
# save( fit, file=fn_test)
# load(fn_test)



s = summary(fit)
s$dic$dic  # 32072
s$dic$p.eff # 4865


# reformat predictions into matrix form
ii = which(M$tag=="predictions")
out = reformat_to_array(
  input = fit$summary.fitted.values[ ii, "mean" ],
  matchfrom = list( AUID=M$AUID[ii], yr_factor=M$yr_factor[ii]),
  matchto   = list( AUID=sppoly$AUID, yr_factor=factor(p$yrs) )
)
# out[ out>1e10] = NA
RES$poisson_basic = colSums( {out * weight_year * sppoly$au_sa_km2}, na.rm=TRUE )  / 10^6 # kt
RES$poisson_basic_cfanorth = colSums( {out * weight_year * sppoly$cfanorth_surfacearea}, na.rm=TRUE ) / 10^6  # 10^6 kg -> kt # kg/km * km
RES$poisson_basic_cfasouth = colSums( {out * weight_year * sppoly$cfasouth_surfacearea}, na.rm=TRUE ) / 10^6  # 10^6 kg -> kt # kg/km * km
RES$poisson_basic_cfa4x = colSums( {out * weight_year * sppoly$cfa4x_surfacearea}, na.rm=TRUE ) / 10^6  # 10^6 kg -> kt # kg/km * km

lines( poisson_basic ~ yr, data=RES, lty=1, lwd=2.5, col="blue", type="b")
lines( poisson_basic_cfanorth ~ yr, data=RES, lty=1, lwd=2.5, col="green", type="b")
lines( poisson_basic_cfasouth ~ yr, data=RES, lty=1, lwd=2.5, col="green", type="b")
lines( poisson_basic_cfa4x ~ yr, data=RES, lty=1, lwd=2.5, col="green", type="b")


# --------------------------------------
# simple with priors
fit = inla(
  formula =
    Y ~ 1 + offset( log( data_offset) )
      + f(auid, model="iid", hyper=H$iid)
      + f(year, model="iid", hyper=H$iid )
   ,
  family = "poisson", # "zeroinflatedpoisson0",
  data= M,
  control.compute=list(dic=TRUE, config=TRUE),
  control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
  control.predictor=list(compute=FALSE, link=1 ),
  control.fixed=H$fixed,  # priors for fixed effects, generic is ok
  control.inla=list(int.strategy="eb") ,# to get empirical Bayes results much faster.
  #  # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
  num.threads=4,
  blas.num.threads=4,
  verbose=TRUE
)

fn3 ="~/tmp/snowcrab_30km_iid_space_time.Rdata"
# save( fit, file=fn3)
# load(fn3)

# Fixed effects:/
#               mean   sd 0.025quant 0.5quant 0.975quant   mode kld
# (Intercept) -6.269 1.12     -8.468   -6.269     -4.071 -6.269   0

# Random effects:
#   Name	  Model
#     auid BYM2 model
#    year IID model
#    iid_error IID model
#    ti RW2 model
#    tidday RW2 model
#    zi RW2 model
#    zid RW2 model
#    zidd RW2 model
#    si RW2 model

# Model hyperparameters:
#                           mean    sd 0.025quant 0.5quant 0.975quant   mode
# Precision for auid     0.738 0.017      0.707    0.738      0.772  0.736
# Phi for auid           0.328 0.007      0.314    0.328      0.343  0.328
# Precision for year      18.798 0.425     17.974   18.794     19.645 18.788
# Precision for iid_error  1.140 0.021      1.099    1.140      1.182  1.140
# Precision for ti         3.352 0.070      3.216    3.351      3.493  3.349
# Precision for tidday     3.524 0.102      3.329    3.523      3.728  3.520
# Precision for zi         3.465 0.076      3.322    3.462      3.621  3.453
# Precision for zid        4.847 0.111      4.622    4.850      5.058  4.864
# Precision for zidd      54.613 1.207     52.279   54.599     57.025 54.572
# Precision for si         0.054 0.002      0.051    0.054      0.057  0.054

# Expected number of effective parameters(stdev): 5008.73(2.07)
# Number of equivalent replicates : 1.52

# Deviance Information Criterion (DIC) ...............: 32072.22
# Deviance Information Criterion (DIC, saturated) ....: 13129.59
# Effective number of parameters .....................: 4865.00

# Marginal log-Likelihood:  -23435.97


s = summary(fit)
s$dic$dic  # 32072
s$dic$p.eff # 4865


# reformat predictions into matrix form
ii = which(M$tag=="predictions")
out = reformat_to_array(
  input = fit$summary.fitted.values[ ii, "mean" ],
  matchfrom = list( AUID=M$AUID[ii], yr_factor=M$yr_factor[ii]),
  matchto   = list( AUID=sppoly$AUID, yr_factor=factor(p$yrs) )
)
# out[ out>1e10] = NA
RES$poisson_basic2 = colSums( {out * weight_year * sppoly$au_sa_km2}, na.rm=TRUE ) /10^6  # km
lines( poisson_basic2 ~ yr, data=RES, lty=1, lwd=2.5, col="blue", type="b")


# map it
vn = "pred"
yr = "2017"
sppoly@data[,vn] = out[,yr] * weight_year[,yr]  # biomass density
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )




# -----------------------------------------
# car simple each posterior config takes @30km:: 10sec .. x 25 configs = 4 min  //  @20km :: 40 sec x 25 config = 20 min // @ 10 km 123.90s tot 45 min
fit = inla(
  formula =
    Y ~ 1 + offset( log( data_offset) )
      + f(auid, model="bym2", graph=sppoly@nb, scale.model=TRUE, constr=TRUE, hyper=H$bym2)
      + f(year, model="iid", hyper=H$iid )
    ,
  family = "poisson", # "zeroinflatedpoisson0",
  data= M,
  control.compute=list(dic=TRUE, config=TRUE),
  control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
  control.predictor=list(compute=FALSE, link=1 ),
  # control.fixed=H$fixed,  # priors for fixed effects, generic is ok
  # control.inla=list(int.strategy="eb") ,# to get empirical Bayes results much faster.
   # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
  num.threads=2,
  blas.num.threads=2,
  verbose=TRUE
)


fn80 ="~/tmp/snowcrab_30km_bym_iid.Rdata"
# save( fit, file=fn80)
# load(fn80)

s = summary(fit)
s$dic$dic  # 31225
s$dic$p.eff # 5200

plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )

# reformat predictions into matrix form
ii = which(M$tag=="predictions")
out = reformat_to_array(
  input = fit$summary.fitted.values[ ii, "mean" ],
  matchfrom = list( AUID=M$AUID[ii], yr_factor=M$yr_factor[ii]),
  matchto   = list( AUID=sppoly$AUID, yr_factor=factor(p$yrs) )
)
# out[ out>1e10] = NA
RES$poisson_car_simple = colSums( {out * weight_year * sppoly$au_sa_km2}, na.rm=TRUE ) / 1e6

lines( poisson_car_simple ~ yr, data=RES, lty=1, lwd=2.5, col="blue", type="b")

# map it
vn = "pred"
yr = "2017"
sppoly@data[,vn] = out[,yr] * weight_year[,yr]  # biomass density
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )



# -----------------------------------------
# simple car grouped by year  27 configs x 100s each = 60 min
fit = inla(
  formula =
    Y ~ 1 + offset( log( data_offset) )
      + f(auid, model="bym2", graph=sppoly@nb, group=year, scale.model=TRUE, constr=TRUE, hyper=H$bym2)
      + f(year, model="iid", hyper=H$iid )
      # + f(ti, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      # + f(tisd, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      # + f(timin, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      # + f(timax, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      # + f(tidday, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      # + f(zi, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      # + f(zid, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      # + f(zidd, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      # + f(si, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2),
      ,
  family = "poisson", # "zeroinflatedpoisson0",
  data= M,
  control.compute=list(dic=TRUE, config=TRUE),
  control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
  control.predictor=list(compute=FALSE, link=1 ),
  # control.fixed=H$fixed,  # priors for fixed effects, generic is ok
  # control.inla=list(int.strategy="eb") ,# to get empirical Bayes results much faster.
   # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
  verbose=FALSE
)

fn90 ="~/tmp/snowcrab_30km_bym|yr.Rdata"
# save( fit, file=fn90)
# load(fn90)

s = summary(fit)
s$dic$dic  # 31107
s$dic$p.eff # 5124

# Fixed effects:
#              mean     sd 0.025quant 0.5quant 0.975quant  mode kld
# (Intercept) 5.564 0.1096      5.348    5.564       5.78 5.565   0

# Random effects:
# Name	  Model
#  auid   BYM2 model
# year   IID model
# iid_error   IID model

# Model hyperparameters:
#                           mean     sd 0.025quant 0.5quant 0.975quant   mode
# Precision for auid    0.3016 0.0408     0.2274   0.2996     0.3879 0.2963
# Phi for auid          0.9970 0.0051     0.9827   0.9991     1.0000 0.9999
# GroupRho for auid     0.9176 0.0146     0.8864   0.9185     0.9434 0.9201
# Precision for year      6.3591 2.4395     2.9280   5.9114    12.3501 5.1296
# Precision for iid_error 0.7329 0.0203     0.6955   0.7320     0.7751 0.7291

# Expected number of effective parameters(std dev): 5278.80(16.21)
# Number of equivalent replicates : 1.399

# Deviance Information Criterion (DIC) ...............: 31107.11
# Deviance Information Criterion (DIC, saturated) ....: 12918.07
# Effective number of parameters .....................: 5124.12

# Marginal log-Likelihood:  -18060.75

# reformat predictions into matrix form
ii = which(M$tag=="predictions")
out = reformat_to_array(
  input = fit$summary.fitted.values[ ii, "mean" ],
  matchfrom = list( AUID=M$AUID[ii], yr_factor=M$yr_factor[ii]),
  matchto   = list( AUID=sppoly$AUID, yr_factor=factor(p$yrs) )
)
# out[ out>1e10] = NA
RES$poisson_car.year_iid_yr = colSums( {out * weight_year * sppoly$au_sa_km2}, na.rm=TRUE ) / 10^6
plot( poisson_car.year_iid_yr ~ yr, data=RES, lty=1, lwd=2.5, col="blue", type="b")

# map it
vn = "pred"
yr = "2017"
sppoly@data[,vn] = out[,yr] * weight_year[,yr]  # biomass density
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )








# -----------------------------------------
# simple car grouped by year  45 configs x 22s each = 60 min // 25 min @ 20 km
fit = inla(
  formula =
    Y ~ 1 + offset( log( data_offset) )
      + f(auid, model="bym2", graph=sppoly@nb, scale.model=TRUE, constr=TRUE, hyper=H$bym2)
      + f(year, model="iid", hyper=H$iid )
      + f(ti, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      # + f(tisd, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      # + f(timin, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      # + f(timax, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      # + f(tidday, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
       + f(zi, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      # + f(zid, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      # + f(zidd, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      # + f(si, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      ,
  family = "poisson", # "zeroinflatedpoisson0",
  data= M,
  control.compute=list(dic=TRUE, config=TRUE),
  control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
  control.predictor=list(compute=FALSE, link=1 ),
  # control.fixed=H$fixed,  # priors for fixed effects, generic is ok
  # control.inla=list(int.strategy="eb") , # to get empirical Bayes results much faster.
   # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
  verbose=TRUE
)

fn100 ="~/tmp/snowcrab_30km_bym_envir_temp.Rdata"
# save( fit, file=fn100)
# load(fn100)

s = summary(fit)
s$dic$dic  # 30288
s$dic$p.eff # 4942

# Fixed effects:
#              mean     sd 0.025quant 0.5quant 0.975quant  mode kld
# (Intercept) 1.324 0.0679      1.191    1.324      1.458 1.324   0

# Random effects:
# Name	  Model
#  auid   BYM2 model
# year   IID model
# iid_error   IID model
# ti   RW2 model
# zi   RW2 model

# Model hyperparameters:
#                            mean     sd 0.025quant 0.5quant 0.975quant   mode
# Precision for auid     0.3665 0.0582     0.2686   0.3606     0.4970 0.3479
# Phi for auid           0.4514 0.1155     0.2228   0.4554     0.6655 0.4796
# Precision for year       0.5738 0.1993     0.2793   0.5421     1.0536 0.4841
# Precision for iid_error  0.6356 0.0168     0.6032   0.6354     0.6693 0.6349
# Precision for ti        12.5146 9.8082     2.0243   9.9645    38.1867 5.4997
# Precision for zi         0.1644 0.1040     0.0375   0.1412     0.4284 0.0958

# Expected number of effective parameters(std dev): 5348.15(21.85)
# Number of equivalent replicates : 1.381

# Deviance Information Criterion (DIC) ...............: 30288.36
# Deviance Information Criterion (DIC, saturated) ....: 12099.32
# Effective number of parameters .....................: 4941.85

# Marginal log-Likelihood:  -23563.46
# Posterior marginals for linear predictor and fitted values computed


# reformat predictions into matrix form
ii = which(M$tag=="predictions")
out = reformat_to_array(
  input = fit$summary.fitted.values[ ii, "mean" ],
  matchfrom = list( AUID=M$AUID[ii], yr_factor=M$yr_factor[ii]),
  matchto   = list( AUID=sppoly$AUID, yr_factor=factor(p$yrs) )
)
# out[ out>1e10] = NA
RES$poisson_car.year_iid_yr = colSums( {out * weight_year * sppoly$au_sa_km2}, na.rm=TRUE ) / 10^6
plot( poisson_car.year_iid_yr ~ yr, data=RES, lty=1, lwd=2.5, col="blue", type="b")

# map it
vn = "pred"
yr = "2017"
sppoly@data[,vn] = out[,yr] * weight_year[,yr]  # biomass density
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )






# -----------------------------------------
# simple car grouped by year nn configs x nn s each = nn min
fit = inla(
  formula =
    Y ~ 1 + offset( log( data_offset) )
      + f(auid, model="bym2", graph=sppoly@nb, scale.model=TRUE, constr=TRUE, hyper=H$bym2)
      + f(year, model="iid", hyper=H$iid )
      + f(ti, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      + f(tisd, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      + f(timin, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      + f(timax, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      + f(tidday, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      + f(zi, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      + f(zid, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      + f(zidd, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      + f(si, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      ,
  family = "poisson", # "zeroinflatedpoisson0",
  data= M,
  control.compute=list(dic=TRUE, config=TRUE),
  control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
  control.predictor=list(compute=FALSE, link=1 ),
  # control.fixed=H$fixed,  # priors for fixed effects, generic is ok
  # control.inla=list(int.strategy="eb") , # to get empirical Bayes results much faster.
   # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
  verbose=TRUE
)

fn101 ="~/tmp/snowcrab_30km_bym_envir_all.Rdata"
# save( fit, file=fn101)
# load(fn101)

s = summary(fit)
s$dic$dic  # 30288
s$dic$p.eff # 4942

# Fixed effects:
#              mean     sd 0.025quant 0.5quant 0.975quant  mode kld
# (Intercept) 1.324 0.0679      1.191    1.324      1.458 1.324   0

# Random effects:
# Name	  Model
#  auid   BYM2 model
# year   IID model
# iid_error   IID model
# ti   RW2 model
# zi   RW2 model

# Model hyperparameters:
#                            mean     sd 0.025quant 0.5quant 0.975quant   mode
# Precision for auid     0.3665 0.0582     0.2686   0.3606     0.4970 0.3479
# Phi for auid           0.4514 0.1155     0.2228   0.4554     0.6655 0.4796
# Precision for year       0.5738 0.1993     0.2793   0.5421     1.0536 0.4841
# Precision for iid_error  0.6356 0.0168     0.6032   0.6354     0.6693 0.6349
# Precision for ti        12.5146 9.8082     2.0243   9.9645    38.1867 5.4997
# Precision for zi         0.1644 0.1040     0.0375   0.1412     0.4284 0.0958

# Expected number of effective parameters(std dev): 5348.15(21.85)
# Number of equivalent replicates : 1.381

# Deviance Information Criterion (DIC) ...............: 30288.36
# Deviance Information Criterion (DIC, saturated) ....: 12099.32
# Effective number of parameters .....................: 4941.85

# Marginal log-Likelihood:  -23563.46
# Posterior marginals for linear predictor and fitted values computed


# reformat predictions into matrix form
ii = which(M$tag=="predictions")
out = reformat_to_array(
  input = fit$summary.fitted.values[ ii, "mean" ],
  matchfrom = list( AUID=M$AUID[ii], yr_factor=M$yr_factor[ii]),
  matchto   = list( AUID=sppoly$AUID, yr_factor=factor(p$yrs) )
)
# out[ out>1e10] = NA
RES$poisson_car.year_iid_yr = colSums( {out * weight_year * sppoly$au_sa_km2}, na.rm=TRUE ) / 10^6
plot( poisson_car.year_iid_yr ~ yr, data=RES, lty=1, lwd=2.5, col="blue", type="b")

# map it
vn = "pred"
yr = "2017"
sppoly@data[,vn] = out[,yr] * weight_year[,yr]  # biomass density
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )






# ------------------------------------------------
# Model 1: a binomial (presence-absence) aks, habitat probability model with linear covariate effects

set$Y = trunc(set$Y)

fit = glm(
  formula = Y ~  offset(data_offset) + 1 + z + dZ + ddZ + t + tsd, #+ tmin + tmax + degreedays + log(substrate.grainsize),
  family=poisson(link="log"),
  data=set[ok,]
)


fit = glm(
  formula = Y ~ 1 + z + dZ + ddZ + t + tsd + tmin + tmax + degreedays + log(substrate.grainsize),
  family=binomial(link="logit"),
  data=set[ok,]
)

s = summary(fit)
AIC(fit)  # 20686688

toplot = expand.grid( t=seq(-1,20), z=c(5,10,20,40,80,160,320,640),degreedays=seq(0, 5000, by=100) )
toplot$predictions = predict(fit, newdata=toplot, type="response", se.fit=FALSE  )

plot( predictions ~ z, toplot[ which( {toplot$t==min(toplot$t)} & {toplot$degreedays==min(toplot$degreedays)} ), ], type="b" )
plot( predictions ~ t, toplot[ which( {toplot$z==min(toplot$z)} & {toplot$degreedays==min(toplot$degreedays)} ), ], type="b" )
plot( predictions ~ degreedays, toplot[ which( {toplot$z==min(toplot$z)} & {toplot$t==min(toplot$t)} ), ], type="b")


# predicted probabilities of observing cod, given covariates (temperature, depth, etc)

APS = aegis_prediction_surface( aegis_data=res$means )
APS$data_offset=1
APS$yr = APS$year
APS$yr_factor = factor( APS$year, levels=p$yrs)
APS$iyr = match(APS$yr_factor, p$yrs)
APS$iauid = match( APS$AUID, sppoly$AUID )

predictions = predict(fit, APS, type="response", se.fit=TRUE  )
APS$predictions = predictions$fit
APS$predictions.se = predictions$se.fit

# reformat predictions into matrix form
out = matrix(NA, nrow=length(sppoly$AUID), ncol=length(p$yrs), dimnames=list( sppoly$AUID, p$yrs) )
out[ cbind(APS$iauid, APS$iyr) ] = APS$predictions
RES$habitat_glm = colSums( {out * sppoly$au_sa_km2 }, na.rm=TRUE ) /sum(sppoly$au_sa_km2) # sa weighted average prob habitat


# map it onto auid means of temperature and depth
aps = APS[ APS$year==2017,  ]
iy = match( as.character(sppoly$AUID), aps$AUID )
vn = "pred"
sppoly@data[,vn] = APS$predictions[iy]
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )




# ------------------------------------------------
# Model 5b: habitat model with a smoothed covariate effect
fit = mgcv::gam(
  formula = pa ~  1 + s(t, k=3, bs="ts")  +  s(z, k=3, bs="ts")  +  s(degreedays, k=3, bs="ts"),
  family=binomial(link="logit"),
  data=set[ok,]
)

fit = mgcv::gam(
  formula = Y ~ offset(log(data_offset)) + 1 + AUID:yr_factor + AUID + yr_factor + s(t, bs="tp", k=3)  + s(z, bs="tp", k=3),
  family=poisson(link="log"),
  data=set[ok, ]
)

s = summa

inverse.logit = function( x ) {
  # x should be the log odds ratio
  oddsratio = exp(x)
  prob = oddsratio / (1 + oddsratio )
  return (prob)
}


plot(fit, all.terms=TRUE, trans=inverse.logit, seWithMean=TRUE, jit=TRUE, rug=TRUE )
s = summary(fit)
AIC(fit)  #  10579 .. a little better than the glm

toplot = expand.grid( t=seq(-1,20), z=c(5,10,20,40,80,160,320,640),degreedays=seq(0, 5000, by=100) )
toplot$predictions = predict(fit, newdata=toplot, type="response", se.fit=FALSE  )

plot( predictions ~ z, toplot[ which( {toplot$t==min(toplot$t)} & {toplot$degreedays==min(toplot$degreedays)} ), ], type="b" )
plot( predictions ~ t, toplot[ which( {toplot$z==min(toplot$z)} & {toplot$degreedays==min(toplot$degreedays)} ), ], type="b" )
plot( predictions ~ degreedays, toplot[ which( {toplot$z==min(toplot$z)} & {toplot$t==min(toplot$t)} ), ], type="b")

# predicted probabilities of observing cod, given temperature and depth

APS = aegis_prediction_surface( aegis_data=res$means )
APS$data_offset=1
APS$yr = APS$year
APS$yr_factor = factor( APS$year, levels=p$yrs)
APS$iyr = match(APS$yr_factor, p$yrs)
APS$iauid = match( APS$AUID, sppoly$AUID )

predictions = predict(fit, APS, type="response", se.fit=TRUE  )
APS$predictions = predictions$fit
APS$predictions.se = predictions$se.fit

# reformat predictions into matrix form
out = matrix(NA, nrow=length(sppoly$AUID), ncol=length(p$yrs), dimnames=list( sppoly$AUID, p$yrs) )
out[ cbind(APS$iauid, APS$iyr) ] = APS$predictions

RES$habitat_gam = colSums( {out * sppoly$au_sa_km2 }, na.rm=TRUE ) /sum(sppoly$au_sa_km2) # sa weighted average prob habitat


# map it
iy = which( APS$year=="2017")
it = match( as.character(sppoly$AUID), APS$AUID[iy] )
vn = "pred"
sppoly@data[,vn] = APS$predictions[iy][it]
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )


# ------------------------------------------------
# Model 5c: as above but via INLA  ... very slow
## similar to GAM

APS = aegis_prediction_surface( aegis_data=res$means )
APS$data_offset=1
APS$pa = NA  # what we are trying to predict
APS$tag = "predictions"
APS$yr = APS$year
APS$yr_factor = factor( APS$year, levels=p$yrs)

basic_vars = unique(c( "pa", "t", "z", "degreedays", "data_offset", "tag", "yr", "AUID"))

M = rbind( set[, basic_vars], APS[,basic_vars] )

M$t[!is.finite(M$t)] = median(M$t, na.rm=TRUE )  # missing data .. quick fix .. do something better for
M$z[!is.finite(M$z)] = median(M$z, na.rm=TRUE )  # missing data .. quick fix .. do something better for

M$yr_factor = factor( M$yr, levels=p$yrs )
M$AUID  = factor( M$AUID, levels=levels(sppoly$AUID ))


M$ti = discretize_data( M$t, p$discretization$t )
M$di = discretize_data( M$t, p$discretization$degreedays )
M$zi = discretize_data( M$t, p$discretization$z )



fit = inla(
  formula = pa ~ 1
      + f(ti, model="rw2", scale.model=TRUE, diagonal=1e-5, hyper=H$prec)
      + f(zi, model="rw2", scale.model=TRUE, diagonal=1e-5, hyper=H$prec)
      + f(di, model="rw2", scale.model=TRUE, diagonal=1e-5, hyper=H$prec),
  family="binomial",  # alternates family="zeroinflatedbinomial0", family="zeroinflatedbinomial1",
  data=M,
  control.family=list(control.link=list(model="logit")),
  control.compute=list(dic=TRUE, config=TRUE),
  control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
  control.predictor=list(compute=TRUE, link=1 ), # compute=TRUE on each data location
  control.fixed=H$fixed,  # priors for fixed effects
  control.inla=list(  correct=TRUE, correct.verbose=FALSE ), # strategy="laplace", cutoff=1e-6,
  verbose=TRUE
)
plot(fit )
s = summary(fit)
s$dic$dic  # 10585   .. not sure why ..
s$dic$p.eff # 17.73

APS = cbind( APS, fit$summary.fitted.values[ which(M$tag=="predictions"), ] )
APS$iyr = match(APS$yr_factor, p$yrs)
APS$iauid = match( APS$AUID, sppoly$AUID )

# reformat predictions into matrix form
out = matrix(NA, nrow=length(sppoly$AUID), ncol=length(p$yrs), dimnames=list( sppoly$AUID, p$yrs) )
out[ cbind(APS$iauid, APS$iyr) ] = APS$mean
RES$habitat_inla = colSums( {out * sppoly$au_sa_km2 }, na.rm=TRUE ) /sum(sppoly$au_sa_km2) # sa weighted average prob habitat

# map it
vn = "pred"
sppoly@data[,vn] = out[,"2017"]
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )



# NOTE most of the decline occured when "habitat" was stable !
# NOTE habitat has become more variable since 1995
# NOTE habitat has declined in 2017
# to do : compute SE's and add to the graph



## -------------------------------------------------------
# bym, iid_year
# Model XX:


APS = aegis_prediction_surface( aegis_data=res$means )
APS$data_offset=1
APS$pa = NA  # what we are trying to predict
APS$tag = "predictions"
APS$yr = APS$year
APS$yr_factor = factor( APS$year, levels=p$yrs)

basic_vars = unique(c( "pa", "t", "z", "degreedays", "data_offset", "tag", "yr", "AUID"))

M = rbind( set[, basic_vars], APS[,basic_vars] )

M$t[!is.finite(M$t)] = median(M$t, na.rm=TRUE )  # missing data .. quick fix .. do something better for
M$z[!is.finite(M$z)] = median(M$z, na.rm=TRUE )  # missing data .. quick fix .. do something better for

M$yr_factor = factor( as.character(M$yr) )
M$AUID  = factor( M$AUID, levels=levels(sppoly$AUID ))
M$auid  = as.numeric( factor(M$AUID))
M$year  = as.numeric( M$yr_factor)


M$ti = discretize_data( M$t, p$discretization$t )
M$di = discretize_data( M$t, p$discretization$degreedays )
M$zi = discretize_data( M$t, p$discretization$z )


fit = inla(
  formula = pa ~ 1
  + f(ti, model="rw2", scale.model=TRUE, hyper=H$rw2)
  + f(zi, model="rw2", scale.model=TRUE, hyper=H$rw2)
  + f(di, model="rw2", scale.model=TRUE, hyper=H$rw2)
  + f(year, model="iid", hyper=H$iid)
  + f(auid, model="bym2", graph=sppoly@nb, scale.model=TRUE, constr=TRUE, hyper=H$bym2),
  family="binomial",  # alternates family="zeroinflatedbinomial0", family="zeroinflatedbinomial1",
  data=M,
  control.family=list(control.link=list(model="logit")),
  control.compute=list(dic=TRUE, config=TRUE),
  control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
  control.predictor=list(compute=TRUE, link=1 ), # compute=TRUE on each data location
  control.fixed=H$fixed,  # priors for fixed effects
  control.inla=list(  correct=TRUE, correct.verbose=FALSE ), # strategy="laplace", cutoff=1e-6,
  verbose=TRUE
)
plot(fit )
s = summary(fit)
s$dic$dic  #8885
s$dic$p.eff #163.8

APS = cbind( APS, fit$summary.fitted.values[ which(M$tag=="predictions"), ] )

APS$iyr = match(APS$yr_factor, p$yrs)
APS$iauid = match( APS$AUID, sppoly$AUID )

# reformat predictions into matrix form
out = matrix(NA, nrow=length(sppoly$AUID), ncol=length(p$yrs), dimnames=list( sppoly$AUID, p$yrs) )
out[ cbind(APS$iauid, APS$iyr) ] = APS$mean
RES$habitat_bym_yriid = colSums( {out * sppoly$au_sa_km2 }, na.rm=TRUE ) /sum(sppoly$au_sa_km2) # sa weighted average prob habitat

# map it
vn = "pred"
sppoly@data[,vn] = out[,"2017"]
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )



if (0) {
 fn = file.path( "~", "tmp", "RES.rdata" )
 # save(RES, file=fn)
 # load(fn)
}



dev.new(width=11, height=7)
# col = c("slategray", "turquoise", "darkorange", "green", "blue", "darkred", "cyan", "darkgreen", "purple" )
col = c( "darkorange", "green", "blue",  "cyan" )
pch = c(20, 21, 22, 23)#, 24, 25, 26, 27, 20)
lty = c(1, 3, 4, 5)# , 6, 7, 1, 3, 4 )
lwd = c(4, 4, 4, 4)# , 4, 4, 4, 4, 4 )
type =c("l", "l", "l", "l")#, "l", "l", "l", "l", "l")
legend=c("GLM", "GAM", "INLA", "INLA CAR")#, "INLA Envir AR1 CAR", "INLA Envir AR1 CAR|year", "INLA Envir AR1|auid CAR", "INLA Envir AR1|auid CAR|year", "INLA Envir CAR|year")

plot( habitat_glm  ~ yr, data=RES, lty=lty[1], lwd=lwd[1], col=col[1], pch=pch[1], type=type[1], ylim=c(0.375,0.825), xlab="Year", ylab="kg")
lines( habitat_gam ~ yr, data=RES, lty=lty[2], lwd=lwd[2], col=col[2], pch=pch[2], type=type[2])
lines( habitat_inla ~ yr, data=RES, lty=lty[3], lwd=lwd[3], col=col[3], pch=pch[3], type=type[3])
lines( habitat_bym_yriid ~ yr, data=RES, lty=lty[4], lwd=lwd[4], col=col[4], pch=pch[4], type=type[4])  # yr_iid
# lines( model1 ~ yr, data=RES, lty=lty[5], lwd=lwd[5], col=col[5], pch=pch[5], type=type[5])
# lines( INLA.Envir.AR1.CAR_year ~ yr, data=RES, lty=lty[6], lwd=lwd[6], col=col[6], pch=pch[6], type=type[6])
# lines( INLA.Envir.AR1_auid.CAR ~ yr, data=RES, lty=lty[7], lwd=lwd[7], col=col[7], pch=pch[7], type=type[7])
# lines( INLA.Envir.AR1_auid.CAR_year ~ yr, data=RES, lty=lty[8], lwd=lwd[8], col=col[8], pch=pch[8], type=type[8])
# lines( INLA.Envir.yr_iid.CAR_year ~ yr, data=RES, lty=lty[9], lwd=lwd[9], col=col[9], pch=pch[9], type=type[9])



legend("topright", legend=legend, lty=lty, col=col, lwd=lwd )


# end
