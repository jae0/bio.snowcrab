
  parameters_numerical_dynamics = function( yrs=yrs, snowcrab_filter_class = "M0", spec_bio=2526, save_location="" ) {
    # we do this here so that the same parameters are accessible when creating the data wrap up for julia

    carstm_model_label= paste( "default", snowcrab_filter_class, sep="_" )

    # poisson works too but variance is not exactly poisson (higher than mean)
    Nfamily = switch( snowcrab_filter_class, 
      M0 = "nbinomial",   
      M1 = "nbinomial",
      M2 = "nbinomial",
      M3 = "nbinomial",
      M4 = "nbinomial",
      f.mat = "nbinomial"
    )
  
    # params for number
    pN = snowcrab_parameters(
      project_class="carstm",
      yrs=yrs,   
      areal_units_type="tesselation",
      family = Nfamily,  
      carstm_model_label= carstm_model_label,  
      carstm_directory = save_location,
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
      family =  "gaussian",
      carstm_model_label= carstm_model_label,  
      carstm_directory = save_location,
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
      family = "binomial",  # "binomial",  # "nbinomial", "betabinomial", "zeroinflatedbinomial0" , "zeroinflatednbinomial0"
      carstm_model_label= carstm_model_label,  
      carstm_directory = save_location,
      selection = list(
        type = "presence_absence",
        biologicals=list( spec_bio=spec_bio ),
        biologicals_using_snowcrab_filter_class=snowcrab_filter_class
      )
    )

    sppoly_tweaks = list(
    # vary params by variable as data densities vary for these size/age/sex groups 
        areal_units_constraint_ntarget= list( M0=8, M1=10, M2=14, M3=14, M4=14, f.mat=8 ),
        n_iter_drop=list( M0=0, M1=1, M2=1, M3=1, M4=1, f.mat=0  )
    )
 
    pN$areal_units_constraint_ntarget = sppoly_tweaks[["areal_units_constraint_ntarget"]][[snowcrab_filter_class]]
    pN$n_iter_drop = sppoly_tweaks[["n_iter_drop"]][[snowcrab_filter_class]]

    pW$areal_units_constraint_ntarget = sppoly_tweaks[["areal_units_constraint_ntarget"]][[snowcrab_filter_class]]
    pW$n_iter_drop = sppoly_tweaks[["n_iter_drop"]][[snowcrab_filter_class]]

    pH$areal_units_constraint_ntarget = sppoly_tweaks[["areal_units_constraint_ntarget"]][[snowcrab_filter_class]]
    pH$n_iter_drop = sppoly_tweaks[["n_iter_drop"]][[snowcrab_filter_class]]

    return( list(pN=pN, pW=pW, pH=pH ) )
  }
