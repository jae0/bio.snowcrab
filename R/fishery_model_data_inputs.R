

fishery_model_data_inputs = function( year.assessment=2021,  save_location=NULL,
  type="biomass_dynamics", fishery_model_label = "turing1", for_julia=FALSE, time_resolution=1/12, 
  sppoly_tweaks=NULL
) {

  if (0) {
    require(bio.snowcrab)   
    source( file.path( code_root, "bio_startup.R" )  )
    # loadfunctions("bio.snowcrab")
    year.assessment=2023
    # type="biomass_dynamics"
    type="size_structured_numerical_dynamics"
    for_julia=TRUE
    time_resolution=2/52
    fishery_model_label = "turing1" 
    save_location = file.path( homedir, "projects", "dynamical_model", "snowcrab", "data" )
  }

  yrs = 1999:year.assessment
  spec_bio = bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 )


  if (type=="biomass_dynamics") {

    snowcrab_filter_class = "fb"

    carstm_model_label= paste( "default", snowcrab_filter_class, sep="_" )

    p = snowcrab_parameters(
      project_class="carstm",
      yrs=yrs,   
      areal_units_type="tesselation",
      family =  "poisson",  
      carstm_model_label= carstm_model_label,  
      selection = list(
        type = "number",  # number is to get started
        biologicals=list( spec_bio=spec_bio ),
        biologicals_using_snowcrab_filter_class=snowcrab_filter_class
      ),
      fishery_model_label = fishery_model_label
    )
  

    # observations
    eps = 1e-9
    
    landings = bio.snowcrab::snowcrab_landings_db()
      # NOTE:: message( "Fishing 'yr' for CFA 4X has been set to starting year:: 2001-2002 -> 2001, etc.")
      # year is year of capture
      # yr is "fishing year" relative to the assessment cycle
    landings = landings[ which (landings$cfa %in% c( "cfanorth", "cfasouth", "cfa4x" ) ) , ]
    L = tapply( landings$landings, INDEX=landings[,c("yr", "cfa")], FUN=sum, na.rm=T )
    nL = nrow(L)

    cfaall = tapply( landings$landings, INDEX=landings[,c("yr")], FUN=sum, na.rm=T )
    L = as.data.frame( cbind( L, cfaall ) )
    L = L / 1000/1000  # convert to kt 
    
    L$yrs = as.numeric( rownames(L) )
    L =  L[ match( p$yrs, rownames(L) ),  ] 
 
    L = L[ , c(  "yrs", "cfaall", "cfanorth", "cfasouth", "cfa4x") ]
    L$cfa4x[ which(L$yrs==2018)] = 0
    for (i in c( "cfaall", "cfanorth", "cfasouth", "cfa4x") ) {
      j = which(!is.finite(L[,i]) )
      if (length(j) > 0) L[ j, i ] = eps  # remove NA's
    }
      
    # biomass data: post-fishery biomass are determined by survey B)
    B = aggregate_simulations( fn=carstm_filenames( p, returnvalue="filename", fn="aggregated_timeseries" ) )$RES

    rownames(B) = B$yrs
    B = as.data.frame( B[ match( p$yrs, B$yrs ),  ] )

    # cfa4x have had no estimates prior to 2004

    cfanorth.baddata = which( B$yrs <= 2004 )
#    B[ cfanorth.baddata, which(colnames(B)=="cfanorth") ] = NA

    cfasouth.baddata = which( B$yrs <= 2004 )
#    B[ cfasouth.baddata, which(colnames(B)=="cfasouth") ] = NA

    cfa.nodata =   which( B$yrs <= 2004 )
    B[ cfa.nodata , which(colnames(B)=="cfa4x") ] = NA

    Y = B #  index of abundance
   
    Y = as.data.frame(as.matrix(Y))
    er = 0.2  # target exploitation rate
    U = 3  # number of regions
    N = length(p$yrs)  # no years with data
    M = 3 # no years for projections
    ty = which(p$yrs == 2004)  # index of the transition year (2004) between spring and fall surveys
      
  
    # priors
    Kmu =  c( 5.0, 60.0, 1.25  )   ## based upon prior historical analyses (when stmv and kriging were attempted)
    rmu =  c( 1.0, 1.0, 1.0 )    ## biological constraint 
    qmu =  c( 1.0, 1.0, 1.0 )    ## based upon video observations q is close to 1 .. but sampling locations can of course cause bias (avoiding rocks and bedrock)

    Ksd =  c( 0.25, 0.25, 0.25 ) * Kmu   
    rsd =  c( 0.1, 0.1, 0.1 ) * rmu  # smaller SD's to encourage solutions closer to prior means
    qsd =  c( 0.1, 0.1, 0.1 ) * qmu   
    
    if (is.null(save_location)) {
      save_location = file.path( p$modeldir, p$carstm_model_label, "fishery_model_results", p$fishery_model_label )
    }  
    dir.create( save_location , showWarnings=FALSE,  recursive =TRUE)

    fnout = file.path(save_location, "biodyn_biomass.RData")

    if (for_julia) {
      Y= as.data.frame(Y)
      L= as.data.frame(L)
      L$ts = as.numeric(rownames(L))
    }


    save( Y, Kmu, Ksd, L, ty, file=fnout ) 
    message("Data for biomass dynamics model saved to the following location:")
    
    return( fnout )
  
  }
   

  if (type=="numerical_dynamics") {
   
    eps = 1e-9  # small non-zero number
    
    pA = parameters_numerical_dynamics( yrs=yrs, snowcrab_filter_class=snowcrab_filter_class,  spec_bio=spec_bio, save_location=save_location ) 
    pN = pA$pN
    pW = pA$pW
    pH = pA$pH
    pA = NULL
     
    sppoly=areal_units( p=pN )
 
  # sims = carstm_posterior_simulations( pN=pN, pW=pW, pH=pH,  pa_threshold=0.05, qmax=0.99 )
  # sims = sims  / 10^6 # 10^6 kg -> kt;; kt/km^2
 
    carstm_directory = file.path( save_location,  paste( "default", snowcrab_filter_class, sep="_" ) )

  # Hurdle model .. req Hurdle correction
    simsN = carstm_posterior_simulations( pN=pN, pH=pH, pa_threshold=0.05, qmax=0.99, carstm_directory=carstm_directory )
    SN = aggregate_simulations( 
      sims=simsN, 
      sppoly=sppoly, 
      yrs=pN$yrs, 
      method="sum", 
      redo=TRUE 
    ) 
    RESN = SN$RES
    SN = NULL

    simsW = carstm_posterior_simulations( pW=pW, carstm_directory=carstm_directory )
    SW = aggregate_simulations( 
      sims=simsW, 
      sppoly=sppoly, 
      yrs=pW$yrs, 
      method="mean", 
      redo=TRUE 
    ) 
    RESW = SW$RES
    SW = NULL

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
    # L = as.matrix(L)  # catches  
    
    for (i in c( "cfaall", "cfanorth", "cfasouth", "cfa4x") ) {
      j = which(!is.finite(L[,i]) )
      if (length(j) > 0) L[ j, i ] = eps  # remove NA's
    }
    

    # observed index of abundance
    rownames(RESN) = RESN$yrs
    RESN = as.data.frame( RESN )

    cfa4x = which(names(RESN)=="cfa4x") # column index of cfa4x
    cfanorth =  which(names(RESN)=="cfanorth")
    cfasouth =  which(names(RESN)=="cfasouth")
    cfaall =  which(names(RESN)=="cfaall")
    
    cfanorth.baddata = which( RESN$yrs <= 2004 )
#    RESN[ cfanorth.baddata, cfanorth ] = NA

    cfasouth.baddata = which( RESN$yrs <= 2004 )
#    RESN[ cfasouth.baddata, cfasouth ] = NA

    # cfa4x have had no estimates prior to 2004
    cfa.nodata =   which( RESN$yrs <= 2004 )
    RESN[ cfa.nodata , "cfa4x" ] = NA

    for (i in c("cfaall", "cfanorth", "cfasouth", "cfa4x") ) {
      RESN[,paste(i, "cv", sep="_")] = RESN[,paste(i, "sd", sep="_")] / RESN[,i]
    }

    RESN = RESN[, c("yrs", "cfaall", "cfanorth", "cfasouth", "cfa4x", "cfaall_cv", "cfanorth_cv", "cfasouth_cv", "cfa4x_cv") ]

    for (i in c("cfaall", "cfanorth", "cfasouth", "cfa4x") ) {
      RESN[,i] = RESN[,i] / max(RESN[,i], na.rm=T )  # force mean=0 sd=1
    }

    er = 0.2  # target exploitation rate
    U = 3  # number of regions
    N = length(p$yrs)  # no years with data
    M = 3 # no years for projections
    ty = which(p$yrs == 2004)  # index of the transition year (2004) between spring and fall surveys
    
    missing = ifelse( is.finite(as.matrix(RESN) ), 0, 1)
 
    # this must be done last
    i = which(RESN$yrs < 2004); RESN$yrs[i] = RESN$yrs[i] + 0.4  #"spring"
    i = which(RESN$yrs >= 2004); RESN$yrs[i] = Y$yrs[i] + 0.8  # "fall"
    
    RESN = as.matrix(RESN)

    # priors
    Kmu =  c( 5.0, 60.0, 1.25  )   ## based upon prior historical analyses (when stmv and kriging were attempted)
    rmu =  c( 1.0, 1.0, 1.0 )    ## biological constraint 
    qmu =  c( 1.0, 1.0, 1.0 )    ## based upon video observations q is close to 1 .. but sampling locations can of course cause bias (avoiding rocks and bedrock)

    Ksd =  c( 0.25, 0.25, 0.25 ) * Kmu   
    rsd =  c( 0.1, 0.1, 0.1 ) * rmu  # smaller SD's to encourage solutions closer to prior means
    qsd =  c( 0.1, 0.1, 0.1 ) * qmu   
    
    RESN[ which(missing==1)] = NA  
    
    Y = as.data.frame(RESN)
    
    if (is.null(save_location)) {
      save_location = file.path( p$modeldir, p$carstm_model_label, "fishery_model_results", p$fishery_model_label )
    }  

    dir.create( save_location , showWarnings=FALSE, recursive =TRUE)

    fnout = file.path(save_location, "biodyn_number.RData")
    save( Y, Kmu, Ksd, L, ty, file=fnout ) 
    message("Data for numerical dynamics model saved to the following location:")
    
    return( fnout )
  
  }

  
  if (type=="size_structured_numerical_dynamics") {
    
    p = bio.snowcrab::load.environment( year.assessment=year.assessment )

    # observations
    eps = 1e-9  # small non-zero number
    er = 0.2  # target exploitation rate
    U = 3  # number of regions
    N = length(p$yrs)  # no years with data
    M = 3 # no years for projections
    ty = which(p$yrs == 2004)  # index of the transition year (2004) between spring and fall surveys
    
    # data: post-fishery  are determined by survey B)
    
    M0_W = NULL
    Y = data.frame( yrs = p$yrs )
    i = which(Y$yrs < 2004); Y$yrs[i] = Y$yrs[i] + 0.4  #"spring"
    i = which(Y$yrs >= 2004); Y$yrs[i] = Y$yrs[i] + 0.8  # "fall"

    for ( snowcrab_filter_class in c("M0", "M1", "M2", "M3", "M4", "f.mat")) {     
    
        carstm_model_label= paste( "default", snowcrab_filter_class, sep="_" )
  
        pA = parameters_numerical_dynamics( yrs=yrs, snowcrab_filter_class=snowcrab_filter_class,  spec_bio=spec_bio, save_location=save_location ) 
        pN = pA$pN
        pW = pA$pW
        pH = pA$pH
        pA = NULL

        sppoly = areal_units( p=pN )
 
        # sims = carstm_posterior_simulations( pN=pN, pW=pW, pH=pH, pa_threshold=0.05, qmax=0.99 )
        # sims = sims  / 10^6 # 10^6 kg -> kt;; kt/km^2

        if (is.null(save_location)) {
          save_location = file.path( pN$modeldir, pN$carstm_model_label, "fishery_model_results", fishery_model_label )
        }  
         
        dir.create( save_location ,showWarnings=FALSE,  recursive =TRUE)

        fnout = file.path(save_location, "biodyn_number_size_struct.RData")
      
        carstm_directory = file.path( save_location,  paste( "default", snowcrab_filter_class, sep="_" ) )

  #    if (snowcrab_filter_class=="M1") browser()

      # Hurdle model .. req Hurdle correction
        simsN = carstm_posterior_simulations( pN=pN, pH=pH, pa_threshold=0.05, qmax=0.99, carstm_directory=carstm_directory  )
        SN = aggregate_simulations( 
          sims=simsN, 
          sppoly=sppoly, 
          yrs=pN$yrs, 
          method="sum", 
          redo=TRUE 
        ) 
        RESN = SN$RES
        SN = NULL
    
        simsW = carstm_posterior_simulations( pW=pW, carstm_directory=carstm_directory )
        SW = aggregate_simulations( 
          sims=simsW, 
          sppoly=sppoly, 
          yrs=pW$yrs, 
          method="mean", 
          redo=TRUE 
        ) 
        RESW = SW$RES
        SW = NULL

        # index of abundance
        rownames(RESN) = RESN$yrs
        RESN = as.data.frame( RESN )
        
        cfanorth.baddata = which( RESN$yrs <= 2004 )
        cfasouth.baddata = which( RESN$yrs <= 2004 )
        cfa.nodata =   which( RESN$yrs <= 2004 )
 
        RESN[ cfa.nodata , "cfa4x" ] = NA
    
        RESN = RESN[, c("yrs", "cfaall", "cfanorth", "cfasouth", "cfa4x", "cfaall_sd", "cfanorth_sd", "cfasouth_sd", "cfa4x_sd") ]
        # for (i in 2:ncol(RESN)) {
        #   RESN[,i] = RESN[,i] / max(RESN[,i], na.rm=T )  # force (0,1) 
        # }
      
        # this must be done last
        RESN = as.matrix(RESN)
        
        missing = ifelse( is.finite(as.matrix(RESN) ), 0, 1)
    
        RESN[ which(missing==1)] = NA  
        RESN = as.data.frame(RESN)
        RESN$yrs = NULL
        names(RESN) = paste( names(RESN), snowcrab_filter_class, sep="_")
        Y = cbind( Y, RESN )

        # estimated mean size  
        RESW = as.data.frame( RESW )
        rownames(RESW) = RESW$yrs
        colnames(RESW) = paste( "mw", colnames(RESW), sep="_")

        if (snowcrab_filter_class =="M0") {
            M0_W = RESW
        }

        # cfa4x = which(names(RESW)=="cfa4x") # column index of cfa4x
        # cfanorth =  which(names(RESW)=="cfanorth")
        # cfasouth =  which(names(RESW)=="cfasouth")
        # cfaall =  which(names(RESW)=="cfaall")
        
        cfanorth.baddata = which( RESW$yrs <= 2004 )
        cfasouth.baddata = which( RESW$yrs <= 2004 )
        cfa.nodata =   which( RESW$yrs <= 2004 )
 
        RESW[ cfa.nodata , "mw_cfa4x" ] = NA
        RESW$yrs = NULL
        names(RESW) = paste( names(RESW), snowcrab_filter_class, sep="_")
        Y = cbind( Y, RESW )

        # habitat        
        simsH = carstm_posterior_simulations( pH=pH, carstm_directory=carstm_directory  )
        simsH = ifelse( simsH >= p$habitat.threshold.quantile, 1, 0 )
        
        H = data.frame(
          Hcfaall    = rowMeans( colSums( simsH * sppoly$au_sa_km2, na.rm=TRUE ), na.rm=TRUE ),
          Hcfanorth  = rowMeans( colSums( simsH * sppoly$cfanorth_surfacearea, na.rm=TRUE ), na.rm=TRUE),
          Hcfasouth  = rowMeans( colSums( simsH * sppoly$cfasouth_surfacearea, na.rm=TRUE ), na.rm=TRUE),
          Hcfa23     = rowMeans( colSums( simsH * sppoly$cfa23_surfacearea, na.rm=TRUE ), na.rm=TRUE),
          Hcfa24     = rowMeans( colSums( simsH * sppoly$cfa24_surfacearea, na.rm=TRUE ), na.rm=TRUE),
          Hcfa4x     = rowMeans( colSums( simsH * sppoly$cfa4x_surfacearea, na.rm=TRUE ), na.rm=TRUE)
        )
        names(H) = paste(names(H), snowcrab_filter_class, sep="_")
        Y = cbind(Y, H) 
  
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

    # year correction for 4X
    # ie:: fishery from 1999-2000 in 4X coded as 1999 -- ie. assume all captured by 31 Dec
    Lts = lubridate::date_decimal( L$ts )
    to.offset = which( lubridate::month(Lts) >= 1 & lubridate::month(Lts) <= 6 )
    L$yrs_4x = L$yrs
    L$yrs_4x[to.offset] = L$yrs_4x[to.offset] - 1
    
    L = merge(L, M0_W, by.x="yrs", by.y="mw_yrs" )
    
    L[, "cfa4x"] = floor( L[, "cfa4x"] / L[, "mw_cfa4x"] * 1000 * 1000 ) ## kt /kg *1000 *1000 --> number
    L[, "cfasouth"] = floor( L[, "cfasouth"] / L[, "mw_cfasouth"] * 1000 * 1000 ) ## kt /kg *1000 *1000 --> number
    L[, "cfanorth"] = floor( L[, "cfanorth"] / L[, "mw_cfanorth"] * 1000 * 1000 ) ## kt /kg *1000 *1000 --> number
    L[, "cfaall"] = floor( L[, "cfaall"] / L[, "mw_cfaall"] * 1000 * 1000 ) ## kt /kg *1000 *1000 --> number

    L = L[ , c("ts", "yrs", "yrs_4x", "cfaall", "cfanorth", "cfasouth", "cfa4x") ]
    # L = as.matrix(L)  # catches 
    
    for (i in c( "cfaall", "cfanorth", "cfasouth", "cfa4x") ) {
      j = which(!is.finite(L[,i]) )
      if (length(j) > 0) L[ j, i ] = eps  # remove NA's
    }
     
    # priors
    Kmu =  c( 5.0, 60.0, 1.25 )   ## based upon prior historical analyses (when stmv and kriging were attempted)
    rmu =  c( 1.0, 1.0, 1.0 )    ## biological constraint 
    qmu =  c( 1.0, 1.0, 1.0 )    ## based upon video observations q is close to 1 .. but sampling locations can of course cause bias (avoiding rocks and bedrock)

    Ksd =  c( 0.25, 0.25, 0.25 ) * Kmu   
    rsd =  c( 0.1, 0.1, 0.1 ) * rmu  # smaller SD's to encourage solutions closer to prior means
    qsd =  c( 0.1, 0.1, 0.1 ) * qmu   
    
    save( Y, Kmu, Ksd, L, M0_W, file=fnout ) 
    message("Data for stage-structred numerical dynamics model saved to the following location:")
    
    return( fnout )
  
  }


}
