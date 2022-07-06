

fishery_model_data_inputs = function( year.assessment=2021,  
  type="biomass_dynamics", fishery_model_label = "turing1" ) {


  yrs = 1999:year.assessment
  spec_bio = bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 )
  snowcrab_filter_class = "fb"

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


    cfanorth =  1 # column index
    cfasouth =  2 # column index
    cfa4x =  3 # column index



  if (type=="biomass_dynamics") {
  
    # observations

    landings = bio.snowcrab::snowcrab_landings_db()
      # NOTE:: message( "Fishing 'yr' for CFA 4X has been set to starting year:: 2001-2002 -> 2001, etc.")
      # year is year of capture
      # yr is "fishing year" relative to the assessment cycle
    landings = landings[ which (landings$cfa %in% c( "cfanorth", "cfasouth", "cfa4x" ) ) , ]
    L = tapply( landings$landings, INDEX=landings[,c("yr", "cfa")], FUN=sum, na.rm=T )
    nL = nrow(L)

    cfaall = tapply( landings$landings, INDEX=landings[,c("yr")], FUN=sum, na.rm=T )
    L = cbind( L, cfaall )
    L = L / 1000/1000  # convert to kt pN$fishery_model_label = "stan_surplus_production_2022_model_qc_cauchy"
    L[ !is.finite(L)] = 0

    L = as.data.frame( L[ match( p$yrs, rownames(L) ),  ] )
    L = as.matrix(L)  # catches  , assume 20% handling mortality and illegal landings
    L[ which(!is.finite(L)) ] = eps  # remove NA's


    # biomass data: post-fishery biomass are determined by survey B)
    B = aggregate_biomass_from_simulations( fn=carstm_filenames( p, returnvalue="filename", fn="aggregated_timeseries" ) )$RES

    rownames(B) = B$yrs
    B = as.data.frame( B[ match( p$yrs, B$yrs ),  ] )

    # cfa4x have had no estimates prior to 2004

    cfanorth.baddata = which( p$yrs <= 2004 )
#    B[ cfanorth.baddata, cfanorth ] = NA

    cfasouth.baddata = which( p$yrs <= 2004 )
#    B[ cfasouth.baddata, cfasouth ] = NA

    cfa.nodata =   which( p$yrs <= 2004 )
    B[ cfa.nodata , cfa4x ] = NA

    Y = as.matrix(B) # observed index of abundance
    
    oo = apply( Y, 2, range, na.rm=TRUE )
    for (i in 1:ncol(oo)) {
      Y[,i] = (Y[,i] - oo[1,i] )/ diff(oo[,i])  # force median 0.5 with most data inside 0,1
    }

    er = 0.2  # target exploitation rate
    U = 3  # number of regions
    N = length(p$yrs)  # no years with data
    M = 3 # no years for projections
    ty = which(p$yrs == 2004)  # index of the transition year (2004) between spring and fall surveys
    cfa4x = 3 # column index of cfa4x
    eps = 1e-9  # small non-zero number

    missing = ifelse( is.finite(Y), 0, 1)
    missing_n = colSums(missing)
    missing_ntot = sum(missing_n)

    # this must be done last
    Y[ which(!is.finite(Y)) ] = 0 # reset NAs to 0 as stan does not take NAs
    Y[ which(missing==1)] = NA  

    # priors
    Kmu =  c( 5.5, 65.0, 2.0 )   ## based upon prior historical analyses (when stmv and kriging were attempted)
    rmu =  c( 1.0, 1.0, 1.0 )    ## biological constraint 
    qmu =  c( 1.0, 1.0, 1.0 )    ## based upon video observations q is close to 1 .. but sampling locations can of course cause bias (avoiding rocks and bedrock)

    Ksd =  c( 0.25, 0.25, 0.25 ) * Kmu   
    rsd =  c( 0.1, 0.1, 0.1 ) * rmu  # smaller SD's to encourage solutions closer to prior means
    qsd =  c( 0.1, 0.1, 0.1 ) * qmu   
    
   
    odir = file.path( p$modeldir, p$carstm_model_label, "fishery_model_results", p$fishery_model_label )
    fnout = file.path(odir, "biodyn_biomass.RData")
    save( Y, Kmu, Ksd, L, ty, file=fnout ) 
    message("Data for biomass dynamics model saved to the following location:")
    
    return( fnout )
  
  }
  
  
  if (type=="numerical_dynamics") {
  
    # observations
    eps = 1e-9  # small non-zero number
  
    # data: post-fishery  are determined by survey B)
  
    yrs = 1999:p$year.assessment
    spec_bio = bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 )
    snowcrab_filter_class = "M0"     # fishable biomass (including soft-shelled )  "m.mat" "f.mat" "imm"
    
    runlabel= paste( "1999_present", snowcrab_filter_class, sep="_" )
  
    # params for number
    pN = snowcrab_parameters(
      project_class="carstm",
      yrs=yrs,   
      areal_units_type="tesselation",
      family = "nbinomial" ,
      carstm_model_label= runlabel,  
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
      carstm_model_label= runlabel,  
      selection = list(
        type = "meansize",
        biologicals=list( spec_bio=spec_bio ),
        biologicals_using_snowcrab_filter_class=snowcrab_filter_class
      )
    )
  
    sppoly=areal_units( p=pN )
    
    simsN = carstm_posterior_simulations( pN=pN)
    SN = aggregate_biomass_from_simulations( 
      sims=simsN, 
      sppoly=sppoly, 
      yrs=pN$yrs, 
      method="sum", 
      redo=TRUE 
    ) 
    RESN = SN$RES

    simsW = carstm_posterior_simulations( pW=pW)
    SW = aggregate_biomass_from_simulations( 
      sims=simsW, 
      sppoly=sppoly, 
      yrs=pW$yrs, 
      method="mean", 
      redo=TRUE 
    ) 
    RESW = SW$RES


    landings = bio.snowcrab::snowcrab_landings_db()
      # NOTE:: message( "Fishing 'yr' for CFA 4X has been set to starting year:: 2001-2002 -> 2001, etc.")
      # year is year of capture
      # yr is "fishing year" relative to the assessment cycle
    landings = landings[ which (landings$cfa %in% c( "cfanorth", "cfasouth", "cfa4x" ) ) , ]
    landings$timestamp = landings$date.fished
    i = which(is.na( landings$timestamp )) 
    if (length(i) > 0) landings$timestamp[ i ] = landings$date.landed[ i ]

    landings$timestamp = lubridate::as_datetime(landings$date.fished)
    i = which(is.na( landings$timestamp )) 
    if (length(i) > 0) {
      # missing time of year / season .. setting to summer 
      landings$timestamp[ i ] = lubridate::ymd( paste( landings$year[i], "07", "01", sep="-" ) )
    }
    
    landings$dyear = lubridate::decimal_date( landings$timestamp ) - lubridate::year(landings$timestamp )
    time_resolution = 0.1
    landings$dyear = trunc(landings$dyear / time_resolution + 1 ) * time_resolution
    
    landings$ts = landings$year + landings$dyear
   
    L = tapply( landings$landings, INDEX=landings[,c("ts", "cfa")], FUN=sum, na.rm=T )
    nL = nrow(L)

    cfaall = tapply( landings$landings, INDEX=landings[,c("ts")], FUN=sum, na.rm=T )
    L = cbind( L, cfaall )
    L = L / 1000/1000  # convert to kt pN$fishery_model_label = "stan_surplus_production_2022_model_qc_cauchy"
    L[ !is.finite(L)] = 0
    L = as.data.frame(L)
    L$ts = as.numeric( rownames(L) )

    # landings to number:
    L$yrs = floor( L$ts )
    names(RESW) = paste( "mw", names(RESW), sep="_")
    L = merge(L, RESW, by.x="yrs", by.y="mw_yrs" )
    
    L[, "cfa4x"] = floor( L[, "cfa4x"] / L[, "mw_cfa4x"] * 1000 * 1000 ) ## kt /kg *1000 *1000 --> number
    L[, "cfasouth"] = floor( L[, "cfasouth"] / L[, "mw_cfasouth"] * 1000 * 1000 ) ## kt /kg *1000 *1000 --> number
    L[, "cfanorth"] = floor( L[, "cfanorth"] / L[, "mw_cfanorth"] * 1000 * 1000 ) ## kt /kg *1000 *1000 --> number
    L[, "cfaall"] = floor( L[, "cfaall"] / L[, "mw_cfaall"] * 1000 * 1000 ) ## kt /kg *1000 *1000 --> number

    L = L[ , c("ts", "yrs", "cfaall", "cfanorth", "cfasouth", "cfa4x") ]
    # L = as.matrix(L)  # catches  , assume 20% handling mortality and illegal landings
    
    for (i in c( "cfaall", "cfanorth", "cfasouth", "cfa4x") ) {
      j = which(!is.finite(L[,i]) )
      if (length(j) > 0) L[ j, i ] = eps  # remove NA's
    }
    

    # observed index of abundance
    Y = RESN
    rownames(Y) = Y$yrs
    Y = as.data.frame( Y )

    cfa4x = which(names(Y)=="cfa4x") # column index of cfa4x
    cfanorth =  which(names(Y)=="cfanorth")
    cfasouth =  which(names(Y)=="cfasouth")
    cfaall =  which(names(Y)=="cfaall")
    
    cfanorth.baddata = which( Y$yrs <= 2004 )
#    Y[ cfanorth.baddata, cfanorth ] = NA

    cfasouth.baddata = which( Y$yrs <= 2004 )
#    Y[ cfasouth.baddata, cfasouth ] = NA

    # cfa4x have had no estimates prior to 2004
    cfa.nodata =   which( Y$yrs <= 2004 )
    Y[ cfa.nodata , "cfa4x" ] = NA

    Y = Y[, c("yrs", "cfaall", "cfanorth", "cfasouth", "cfa4x") ]
    for (i in 2:ncol(Y)) {
      Y[,i] = scale(Y[,i])  # force mean=0 sd=1
    }

    er = 0.2  # target exploitation rate
    U = 3  # number of regions
    N = length(p$yrs)  # no years with data
    M = 3 # no years for projections
    ty = which(p$yrs == 2004)  # index of the transition year (2004) between spring and fall surveys
    
    missing = ifelse( is.finite(as.matrix(Y) ), 0, 1)
    missing_n = colSums(missing)
    missing_ntot = sum(missing_n)

    # this must be done last
    i = which(Y$yrs < 2004); Y$yrs[i] = Y$yrs[i] + 0.4  #"spring"
    i = which(Y$yrs >= 2004); Y$yrs[i] = Y$yrs[i] + 0.8  # "fall"
    Y = as.matrix(Y)
    Y[ which(!is.finite(Y)) ] = 0 # reset NAs to 0 as stan does not take NAs


    # priors
    Kmu =  c( 5.5, 65.0, 2.0 )   ## based upon prior historical analyses (when stmv and kriging were attempted)
    rmu =  c( 1.0, 1.0, 1.0 )    ## biological constraint 
    qmu =  c( 1.0, 1.0, 1.0 )    ## based upon video observations q is close to 1 .. but sampling locations can of course cause bias (avoiding rocks and bedrock)

    Ksd =  c( 0.25, 0.25, 0.25 ) * Kmu   
    rsd =  c( 0.1, 0.1, 0.1 ) * rmu  # smaller SD's to encourage solutions closer to prior means
    qsd =  c( 0.1, 0.1, 0.1 ) * qmu   
    
    Y[ which(missing==1)] = NA  
    Y = as.data.frame(Y)
    
    odir = file.path( p$modeldir, p$carstm_model_label, "fishery_model_results", p$fishery_model_label )
    fnout = file.path(odir, "biodyn_number.RData")
    save( Y, Kmu, Ksd, L, ty, file=fnout ) 
    message("Data for size structured numerical dynamics model saved to the following location:")
    
    return( fnout )
  
  }

  
  if (type=="size_structured_numerical_dynamics") {
  
    # observations
    eps = 1e-9  # small non-zero number
    er = 0.2  # target exploitation rate
    U = 3  # number of regions
    N = length(p$yrs)  # no years with data
    M = 3 # no years for projections
    ty = which(p$yrs == 2004)  # index of the transition year (2004) between spring and fall surveys
    
    # data: post-fishery  are determined by survey B)
  
    yrs = 1999:p$year.assessment
    spec_bio = bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 )
    
    M0_W = NULL
    YALL = data.frame( yr = p$yrs )
    
    for ( snowcrab_filter_class in c("M0", "M1", "M2")) {     # fishable biomass (including soft-shelled )  "m.mat" "f.mat" "imm"
    
        runlabel= paste( "1999_present", snowcrab_filter_class, sep="_" )
      
        # params for number
        pN = snowcrab_parameters(
          project_class="carstm",
          yrs=yrs,   
          areal_units_type="tesselation",
          family =  "nbinomial", 
          carstm_model_label= runlabel,  
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
          carstm_model_label= runlabel,  
          selection = list(
            type = "meansize",
            biologicals=list( spec_bio=spec_bio ),
            biologicals_using_snowcrab_filter_class=snowcrab_filter_class
          )
        )
      
        sppoly=areal_units( p=pN )
        
        simsN = carstm_posterior_simulations( pN=pN)
        SN = aggregate_biomass_from_simulations( 
          sims=simsN, 
          sppoly=sppoly, 
          yrs=pN$yrs, 
          method="sum", 
          redo=TRUE 
        ) 
        RESN = SN$RES
    
        simsW = carstm_posterior_simulations( pW=pW)
        SW = aggregate_biomass_from_simulations( 
          sims=simsW, 
          sppoly=sppoly, 
          yrs=pW$yrs, 
          method="mean", 
          redo=TRUE 
        ) 
        RESW = SW$RES
    
        # observed index of abundance
        Y = RESN
        rownames(Y) = Y$yrs
        Y = as.data.frame( Y )
    
        cfa4x = which(names(Y)=="cfa4x") # column index of cfa4x
        cfanorth =  which(names(Y)=="cfanorth")
        cfasouth =  which(names(Y)=="cfasouth")
        cfaall =  which(names(Y)=="cfaall")
        
        cfanorth.baddata = which( Y$yrs <= 2004 )
        cfasouth.baddata = which( Y$yrs <= 2004 )
        cfa.nodata =   which( Y$yrs <= 2004 )
 
        Y[ cfa.nodata , "cfa4x" ] = NA
    
        Y = Y[, c("yrs", "cfaall", "cfanorth", "cfasouth", "cfa4x") ]
        for (i in 2:ncol(Y)) {
          Y[,i] = scale(Y[,i])  # force mean=0 sd=1
        }
      
        # this must be done last
        i = which(Y$yrs < 2004); Y$yrs[i] = Y$yrs[i] + 0.4  #"spring"
        i = which(Y$yrs >= 2004); Y$yrs[i] = Y$yrs[i] + 0.8  # "fall"
        Y = as.matrix(Y)
        Y[ which(!is.finite(Y)) ] = 0 # reset NAs to 0 as stan does not take NAs
        
        missing = ifelse( is.finite(Y ), 0, 1)
        missing_n = colSums(missing)
        missing_ntot = sum(missing_n)
    
        Y[ which(missing==1)] = NA  
        Y = as.data.frame(Y)
        Y$yrs = NULL
        names(Y) = paste( names(Y), snowcrab_filter_class, sep="_")
        YALL = cbind( YALL, Y )

        if (snowcrab_filter_class =="M0") {
            M0_W = RESW
            names(M0_W) = paste( "mw", names(M0_W), sep="_")

        }
        message ( snowcrab_filter_class )
    }


    landings = bio.snowcrab::snowcrab_landings_db()
      # NOTE:: message( "Fishing 'yr' for CFA 4X has been set to starting year:: 2001-2002 -> 2001, etc.")
      # year is year of capture
      # yr is "fishing year" relative to the assessment cycle
    landings = landings[ which (landings$cfa %in% c( "cfanorth", "cfasouth", "cfa4x" ) ) , ]
    landings$timestamp = landings$date.fished
    i = which(is.na( landings$timestamp )) 
    if (length(i) > 0) landings$timestamp[ i ] = landings$date.landed[ i ]

    landings$timestamp = lubridate::as_datetime(landings$date.fished)
    i = which(is.na( landings$timestamp )) 
    if (length(i) > 0) {
      # missing time of year / season .. setting to summer 
      landings$timestamp[ i ] = lubridate::ymd( paste( landings$year[i], "07", "01", sep="-" ) )
    }
    
    landings$dyear = lubridate::decimal_date( landings$timestamp ) - lubridate::year(landings$timestamp )
    time_resolution = 0.1
    landings$dyear = trunc(landings$dyear / time_resolution + 1 ) * time_resolution
    
    landings$ts = landings$year + landings$dyear
   
    L = tapply( landings$landings, INDEX=landings[,c("ts", "cfa")], FUN=sum, na.rm=T )
    nL = nrow(L)

    cfaall = tapply( landings$landings, INDEX=landings[,c("ts")], FUN=sum, na.rm=T )
    L = cbind( L, cfaall )
    L = L / 1000/1000  # convert to kt pN$fishery_model_label = "stan_surplus_production_2022_model_qc_cauchy"
    L[ !is.finite(L)] = 0
    L = as.data.frame(L)
    L$ts = as.numeric( rownames(L) )


    # landings to number:
    L$yrs = floor( L$ts )
    
    L = merge(L, M0_W, by.x="yrs", by.y="mw_yrs" )
    
    L[, "cfa4x"] = floor( L[, "cfa4x"] / L[, "mw_cfa4x"] * 1000 * 1000 ) ## kt /kg *1000 *1000 --> number
    L[, "cfasouth"] = floor( L[, "cfasouth"] / L[, "mw_cfasouth"] * 1000 * 1000 ) ## kt /kg *1000 *1000 --> number
    L[, "cfanorth"] = floor( L[, "cfanorth"] / L[, "mw_cfanorth"] * 1000 * 1000 ) ## kt /kg *1000 *1000 --> number
    L[, "cfaall"] = floor( L[, "cfaall"] / L[, "mw_cfaall"] * 1000 * 1000 ) ## kt /kg *1000 *1000 --> number

    L = L[ , c("ts", "yrs", "cfaall", "cfanorth", "cfasouth", "cfa4x") ]
    # L = as.matrix(L)  # catches  , assume 20% handling mortality and illegal landings
    
    for (i in c( "cfaall", "cfanorth", "cfasouth", "cfa4x") ) {
      j = which(!is.finite(L[,i]) )
      if (length(j) > 0) L[ j, i ] = eps  # remove NA's
    }
     

    # priors
    Kmu =  c( 5.5, 65.0, 2.0 )   ## based upon prior historical analyses (when stmv and kriging were attempted)
    rmu =  c( 1.0, 1.0, 1.0 )    ## biological constraint 
    qmu =  c( 1.0, 1.0, 1.0 )    ## based upon video observations q is close to 1 .. but sampling locations can of course cause bias (avoiding rocks and bedrock)

    Ksd =  c( 0.25, 0.25, 0.25 ) * Kmu   
    rsd =  c( 0.1, 0.1, 0.1 ) * rmu  # smaller SD's to encourage solutions closer to prior means
    qsd =  c( 0.1, 0.1, 0.1 ) * qmu   
    
    Y = YALL
   
    odir = file.path( p$modeldir, p$carstm_model_label, "fishery_model_results", p$fishery_model_label )
    fnout = file.path(odir, "biodyn_number_size_struct.RData")
    save( Y, Kmu, Ksd, L, ty, file=fnout ) 
    message("Data for biomass dynamics model saved to the following location:")
    
    return( fnout )
  
  }


}
