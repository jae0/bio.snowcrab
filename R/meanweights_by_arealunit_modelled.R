

meanweights_by_arealunit_modelled = function( p=NULL, redo=FALSE, returntype="predictions_mean" ) {
  # find mean weight for each stratum and year .. fill where possible if required
 
  if (redo) {
    message( "Re-creating mean size estimates from data .. this will take about 2 hrs ...")
    # must be fetched beofre modfying the model_label and variable name
    M = snowcrab.db( p=p, DS="carstm_inputs" ) 
    M$meansize  = M$totwgt / M$totno  # note, these are constrained by filters in size, sex, mat, etc. .. in the initial call
  }

  # alter model formula / storage location
  p_mw = p
  
  p_mw$variabletomodel = "meansize"
  p_mw$carstm_model_label = paste( p_mw$variabletomodel, p_mw$areal_units_type, p_mw$selection$type, sep="_")  

  p_mw$formula = as.formula( paste(
      p_mw$variabletomodel, ' ~ 1',
            ' + f( time, model="ar1",  hyper=H$ar1 ) ',
            ' + f( cyclic, model="seasonal", scale.model=TRUE, season.length=10, hyper=H$iid  )',
            ' + f( inla.group( t, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
            ' + f( inla.group( z, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
            ' + f( inla.group( substrate.grainsize, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
            ' + f( inla.group( pca1, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
            ' + f( inla.group( pca2, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
            ' + f( space, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, constr=TRUE, hyper=H$bym2 ) ',
            ' + f( space_time, model="bym2", graph=slot(sppoly, "nb"), group=time_space, scale.model=TRUE, constr=TRUE, hyper=H$bym2, control.group=list(model="ar1", hyper=H$ar1_group)) '
  ) )

  p_mw$family =  "gaussian"   

  if (redo) {
    fit = carstm_model( p=p_mw, data=M, redo_fit = TRUE, posterior_simulations_to_retain="predictions", control.inla = list( strategy='adaptive' ) ) # 151 configs and long optim .. 19 hrs
    # fit = carstm_model( p=p_mw, data=M, redo_fit = TRUE, posterior_simulations_to_retain="predictions", control.inla = list( strategy='laplace', int.strategy="eb" ) ) # 151 configs and long optim .. 19 hrs
  } 

  if (returntype=="modelled_fit" ) {
    fit = carstm_model( p=p_mw, DS="modelled_fit" )  # extract currently saved model fit
    return(fit)
  }
  
  if (returntype=="predictions_mean" ) {
     res = carstm_model( p=p_mw, DS="carstm_predictions"  ) # to load currently saved results
     return( res[["predictions"]][,,"mean"] )
  }
  if (returntype=="predictions" ) {
     res = carstm_model( p=p_mw, DS="carstm_predictions"  ) # to load currently saved results
     return( res[["predictions"]] )
  }
  if (returntype=="predictions_posterior_simulations" ) {
     res = carstm_model( p=p_mw, DS="carstm_samples"  ) # to load currently saved results
     return( res[["predictions"]] )
  }
  if (returntype=="summary" )  {
    res = carstm_model( p=p_mw, DS="carstm_summary"  ) # to load currently saved results
    return( res )
  }
   
  return( "Need to add a new returntype option?" )

    if (0) {

          plot_crs = p_mw$aegis_proj4string_planar_km
          managementlines = aegis.polygons::area_lines.db( DS="cfa.regions", returntype="sf", project_to=plot_crs )
        
          vn="predictions"
          tmatch = as.character(2020)

          carstm_map(  res=res,  vn=vn, tmatch=tmatch,
              palette="RdYlBu",
              additional_polygons=managementlines,
              title=paste("Predicted mean weight of indivudual (kg)", paste0(tmatch, collapse="-") )
          )
            
    }

}

