

# Snow crab --- Areal unit modelling of habitat   

# -------------------------------------------------
# Part 1 -- construct basic parameter list defining the main characteristics of the study
   year.assessment = 2021
  require(bio.snowcrab)   # loadfunctions("bio.snowcrab") 
 
  pH = snowcrab_parameters( 
    project_class="carstm", 
    yrs=1999:year.assessment,  
    areal_units_type="tesselation", 
    carstm_model_label = "1999_present",
    selection = list(type = "presence_absence")
  )


  if (0) {
    # testing:
      pH$selection$type = "presence_absence"
      pH$selection$biologicals=list(
        spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=pH$groundfish_species_code ),
        sex=0, # male
        mat=1, # do not use maturity status in groundfish data as it is suspect ..
        len= c( 95, 200 )/10, #  mm -> cm ; aegis_db in cm
        ranged_data="len"
      )
      pH$selection$survey=list(
        data.source = c("snowcrab", "groundfish", "logbook"),
        yr = pH$yrs,      # time frame for comparison specified above
        settype = 1, # same as geartype in groundfish_survey_db
        polygon_enforce=TRUE,  # make sure mis-classified stations or incorrectly entered positions get filtered out
        strata_toremove = NULL #,  # emphasize that all data enters analysis initially ..
        # ranged_data = c("dyear")  # not used .. just to show how to use range_data
      )
      pH$variabletomodel = "pa"
      
      pH$carstm_model_label = "nonseparable_space-time_pa_fishable_binomial"
      pH$carstm_modelengine = "inla"
      pH$formula = as.formula( paste(
        pH$variabletomodel, ' ~ 1 ',
          ' + f( dyri, model="ar1", hyper=H$ar1 ) ',
          ' + f( inla.group( t, method="quantile", n=11 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
          ' + f( inla.group( z, method="quantile", n=11 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
          ' + f( inla.group( substrate.grainsize, method="quantile", n=11 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
          ' + f( inla.group( pca1, method="quantile", n=11 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
          ' + f( inla.group( pca2, method="quantile", n=11 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
          ' + f( space, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, constr=TRUE, hyper=H$bym2 ) '
          ' + f( space_time, model="bym2", graph=slot(sppoly, "nb"), group=time_space, scale.model=TRUE, constr=TRUE, hyper=H$bym2, control.group=list(model="ar1", hyper=H$ar1_group)) '
      ) )

      pH$family = "binomial"  # "nbinomial", "betabinomial", "zeroinflatedbinomial0" , "zeroinflatednbinomial0"
      pH$carstm_model_inla_control_familiy = list(control.link=list(model='logit'))

    #  pH$family  = "zeroinflatedbinomial1", #  "binomial",  # "nbinomial", "betabinomial", "zeroinflatedbinomial0" , "zeroinflatednbinomial0"
    #  pH$carstm_model_inla_control_familiy = NULL

  }

  # could create new polygons but simply re-using those for number:
  # params for number
  pN = snowcrab_parameters(
    project_class="carstm",
    yrs=1999:year.assessment,   
    areal_units_type="tesselation",
    family="poisson",
    carstm_model_label = "1999_present",  # 1999_present is the default anything else and you are on your own
    offset_shift=10^3,
    selection = list(type = "number")
  )

  sppoly = areal_units( p=pN )  # to reload
  plot( sppoly[, "au_sa_km2"]  )

  M = snowcrab.db( p=pN, DS="carstm_inputs"  )  # will redo if not found

  fit = carstm_model( p=pH, data=M ) # 151 configs and long optim .. 19 hrs
 

  # fit = carstm_model( p=pH, DS="carstm_modelled_fit")

  # extract results
  if (0) {
    # very large files .. slow 
    fit = carstm_model( p=pH, DS="carstm_modelled_fit" )  # extract currently saved model fit
    plot(fit)
    plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
  }

  res = carstm_model( p=pH, DS="carstm_modelled_summary"  ) # to load currently saved results
  res$summary$dic$dic
  res$summary$dic$p.eff
  res$dyear


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



  additional_features = snowcrab_features_tmap(pH)  # for mapping below

  time_match = list( year=as.character(2020)  )
  tmout = carstm_map(  res=res, 
      vn="predictions", 
      time_match=time_match , 
      # palette="-RdYlBu",
      plot_elements=c(  "compass", "scale_bar", "legend" ),
      additional_features=additional_features,
      main=paste("Habitat probability - mature male ", paste0(time_match, collapse="-") )  
  tmout
  )
    

  # map all :
  vn "= predictions"

  outputdir = file.path( pH$modeldir, pH$carstm_model_label, "predicted.probability.observation" )

  if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

  brks = pretty(  res[[vn]]  )

  for (y in res$year ){

      time_match = list( year=y  )
      fn_root = paste("Predicted_abundance", paste0(time_match, collapse="-"), sep="_")
      fn = file.path( outputdir, paste(fn_root, "png", sep=".") )

      tmout = carstm_map(  
          res=res, 
          vn=vn, 
          time_match=time_match, 
          breaks =brks,
          # palette="-RdYlBu",
          plot_elements=c(  "compass", "scale_bar", "legend" ),
          additional_features=additional_features,
          main=paste("Habitat probability - mature male ", paste0(time_match, collapse="-") ),
          outfilename=fn
        )  

  }
 
  snowcrab.db(p=pH, DS="carstm_output_compute" )
  RES = snowcrab.db(p=pH, DS="carstm_output_timeseries" )
  pa = snowcrab.db(p=pH, DS="carstm_output_spacetime_pa"  )

# plots with 95% PI

  outputdir = file.path( pH$modeldir, pH$carstm_model_label, "aggregated_biomass_timeseries" )

  if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

  (fn = file.path( outputdir, "cfa_all.png"))
  png( filename=fn, width=3072, height=2304, pointsize=12, res=300 )
    plot( cfaall ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Prob of observing snow crab", xlab="",  ylim=c(0,1))
    lines( cfaall_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    lines( cfaall_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
  dev.off()

  (fn = file.path( outputdir, "cfa_south.png") )
  png( filename=fn, width=3072, height=2304, pointsize=12, res=300 )
    plot( cfasouth ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Prob of observing snow crab", xlab="",  ylim=c(0,1))
    lines( cfasouth_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    lines( cfasouth_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
  dev.off()

  (fn = file.path( outputdir, "cfa_north.png"))
  png( filename=fn, width=3072, height=2304, pointsize=12, res=300 )
    plot( cfanorth ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Prob of observing snow crab", xlab="",  ylim=c(0,1))
    lines( cfanorth_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    lines( cfanorth_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
  dev.off()

  (fn = file.path( outputdir, "cfa_4x.png"))
  png( filename=fn, width=3072, height=2304, pointsize=12, res=300 )
    plot( cfa4x ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Prob of observing snow crab", xlab="",  ylim=c(0,1))
    lines( cfa4x_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    lines( cfa4x_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
  dev.off()





 
# end
