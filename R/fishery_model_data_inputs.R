
fishery_model_data_inputs = function( year.assessment=2021, snowcrab_filter_class = "fb", type="biomass_dynamics", fishery_model_label = "turing1" ) {

  yrs = 1999:year.assessment
  spec_bio = bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 )

  runlabel= paste( "1999_present", snowcrab_filter_class, sep="_" )

  p = snowcrab_parameters(
    project_class="carstm",
    yrs=yrs,   
    areal_units_type="tesselation",
    family =  "poisson",  
    carstm_model_label= runlabel,  
    selection = list(
      type = "number",
      biologicals=list( spec_bio=spec_bio ),
      biologicals_using_snowcrab_filter_class=snowcrab_filter_class
    ),
    fishery_model_label = fishery_model_label
  )

  if (type=="biomass_dynamics") {
  
    p$fishery_model = fishery_model( DS = "logistic_parameters", p=p, tag=p$fishery_model_label )

    # N  # no. years
    # U  # no. regions
    # M  # no. years to project
    # ty  # index of year at which transistion of spring -> fall surveys
    # er    # target exploitation rate
    # eps   # small real number
 
    # missing 
    # missing_n 
    # missing_ntot 

    # return NA's 
    
    Y = p$fishery_model$fmdata$IOA
    Y[ which(p$fishery_model$fmdata$missing==1)] = NA  
   
    Kmu = p$fishery_model$fmdata$Kmu
    Ksd = p$fishery_model$fmdata$Ksd
    CAT = p$fishery_model$fmdata$CAT
    ty  = p$fishery_model$fmdata$ty 


    # estimating:
     # eps> [U] K;
     # 0.25, upper=2.0> [U] r;  # biologically should be ~ 1 
     # eps, upper=2.5> [U] q;  # multiplicative factor, unlikely to be >200%
     # -1.0, upper=1.0> [U] qc;  #  offset .. unlikely to be off by > 50%
     # eps, upper=0.5> [U] bosd;  # observation error
     # eps, upper=0.5> [U] bpsd;  # process error
     # eps, upper=0.5> [U] rem_sd;  # catch error
     # eps> [U] b0;
     # eps> [missing_ntot] IOAmissing;
     # lower=eps> [M+N,U] bm;  # force bm max 1 
     # Y  

    fnout = file.path(p$fishery_model$outdir, "biodyn.RData")
    save( Y, Kmu, Ksd, CAT, ty, file=fnout ) 
    message("Data for biomass dynamics model saved to the following location:")
    
    return( fnout )
  
  }

}
