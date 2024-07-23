

# ----------------------------------------------
# apply fishery model to biomass indices
# NOTE::: require 03.snowcrab_carstm.r to be completed 
# (i.e.,spatiotemporal model and assimilate_numbers_and_size to have been completed 
# ----------------------------------------------
 

# TODO::: move plotting calls to self-contained functions:


# -------------------------------------------------
# Part 1 -- construct basic parameter list defining the main characteristics of the study

  source( file.path( code_root, "bio_startup.R" )  )
 

  require(bio.snowcrab)   # loadfunctions("bio.snowcrab") 
  # loadfunctions("bio.snowcrab")


  year.assessment = 2021
  yrs = 1999:year.assessment
  spec_bio = bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 )
  snowcrab_filter_class = "fb"     # fishable biomass (including soft-shelled ) 
  
  carstm_model_label= paste( "default", snowcrab_filter_class, sep="_" )

  # params for number
  p = snowcrab_parameters(
    project_class="carstm",
    yrs=yrs,   
    areal_units_type="tesselation",
    family =  "poisson",  
    carstm_model_label= carstm_model_label,  
    selection = list(
      type = "number",
      biologicals=list( spec_bio=spec_bio ),
      biologicals_using_snowcrab_filter_class=snowcrab_filter_class
    ),
    # fishery_model_label = "stan_surplus_production_2022_model_qc_uniform"
    # fishery_model_label = "stan_surplus_production_2022_model_variation1_wider_qc_uniform"
    # fishery_model_label = "stan_surplus_production_2022_model_variation1_wider_qc_normal"
    # fishery_model_label = "stan_surplus_production_2022_model_qc_cauchy_wider"
    fishery_model_label = "stan_surplus_production_2022_model_qc_beta"  # qc_beta postive definite
    # fishery_model_label = "stan_surplus_production_2022_model_qc_cauchy"
    # fishery_model_label = "stan_surplus_production_2019_model"
  )
 

  p$fishery_model = fishery_model( DS = "logistic_parameters", p=p, tag=p$fishery_model_label )

  to_look = c("K", "r", "q", "qc", "log_lik" )

      if ( p$fishery_model_label=="framework_2019" ) {
        # this is to create results for reviewers in 2022 that wanted a comparison with previous methods .. can be deleted in future 
        # bring in unscaled abundance index
        a = fishery_model( DS="data_aggregated_timeseries", p=p  )
        a$IOA[ !is.finite(a$IOA) ] = 0
        p$fishery_model$fmdata$IOA = a$IOA
        to_look = c("K", "r", "q", "log_lik" )
      }  


  #  str( p$fishery_model)

  p$fishery_model$stancode = stan_initialize( stan_code=fishery_model( p=p, DS=p$fishery_model_label ) )
  p$fishery_model$stancode$compile()
  

  fit = p$fishery_model$stancode$sample(
    data=p$fishery_model$fmdata,
    iter_warmup = 4000,
    iter_sampling = 4000,
    seed = 1,
    chains = 3,
    parallel_chains = 3,  # The maximum number of MCMC chains to run in parallel.
    max_treedepth = 18,
    adapt_delta = 0.99,
    refresh = 1000
  )

  fit$summary(to_look)

  # require(loo)
  # waic(fit$draws("log_lik"))
  # loo(fit$draws("log_lik"))
   
   
  # save fit and get draws
  res = fishery_model( p=p, DS="logistic_model", tag=p$fishery_model_label, fit=fit )       # from here down are params for cmdstanr::sample()

  if (0) {
    # reload saved fit and results
    res = aegis::read_write_fast(p$fishery_model$fnres)
    fit = aegis::read_write_fast(p$fishery_model$fnfit)
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
  np = year.assessment+c(1:p$fishery_model$fmdata$M)
  biomass$yr = rep( c(p$yrs, np ), 3)
  nt = p$fishery_model$fmdata$N +p$fishery_model$fmdata$M
  biomass$region = c( rep("cfanorth", nt), rep("cfasouth", nt), rep("cfa4x", nt) )
  (biomass)

  NN = res$p$fishery_model$fmdata$N

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

 
