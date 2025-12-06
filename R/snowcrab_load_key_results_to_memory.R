
snowcrab_load_key_results_to_memory = function( 
  p, 
  data_loc =  file.path(data_root, "bio.snowcrab"),
  todo = "fishery_results", # "ecosystem", "fishery_model"
  years_model = NULL, 
  envir = parent.frame(),
  mau="region",
  debugging=FALSE,  
  return_as_list=TRUE, 
  redo=FALSE 
) {
  
  # function to bring in key fishery stats and assessment results and make available in memory 
  # primary usage is for Quarto/Rmarkdown documents

  out =NULL
  
  fn_out = file.path( data_loc, "assessments", p$year.assessment,  
    paste("sc_key_", todo, "_", mau, ".rdz", sep="") 
  )
  
  if (!redo) {
    if (file.exists(fn_out)) {
      out = read_write_fast(fn_out)
      if (return_as_list) {
        return( invisible( out ) )
      } else {
        return( invisible( list2env( out, envir) ) )
      }
    }
  }


  year_assessment = p$year.assessment
  year_previous = year_assessment - 1
 
  p$corners = data.frame(plon=c(220, 990), plat=c(4750, 5270) )
  p$mapyears = year_assessment + c(-5:0 )   # default in case not specified
  years = as.character(1996: year_assessment)
  yrs_observer = year_assessment + c(0:-4)
  years_to_show = c(year_assessment - c(0:4))


  # mau is variable name in logbook.db() , snowcrab.db() and logbook.db()
 
  maus = management_areal_units( mau=mau )  


  if ( "fishery_results" %in% todo ) {
  
    FD = fishery_data( mau=mau )  # mass in tonnes

    # fishery summaries
    dtvars = c( mau, "yr", "Licenses", "TAC", "landings", "effort", "cpue")
    dt = FD$summary_annual[ yr %in% years_to_show, ..dtvars ]  
    names(dt) = c(Capitalize(mau), "Year", "Licenses", "TAC", "Landings", "Effort", "CPUE") 
      

    # observer data
    odb0 = setDT(observer.db("odb"))

    # bycatch summaries
    BC = list()
    for ( reg in c(maus[["internal"]], "cfaall")) {
        oo = observer.db( DS="bycatch_summary", p=p,  yrs=p$yrs, region=reg  )  # using polygon test

        if (is.null(oo)) message( "bycatch data is likely broken .. check" ) 

        oo$bycatch_table_effort[ oo$bycatch_table_effort==0 ] = NA
        oo$bycatch_table_effort[ is.na(oo$bycatch_table_effort) ] = "."
        oo$bycatch_table_catch[ oo$bycatch_table_catch==0 ] = NA
        oo$bycatch_table_catch[ is.na(oo$bycatch_table_catch) ] = "."

        BC[[reg]] = oo
    }
  

    # landings
    ly = lp = ldt = list()   
    for ( reg in maus[["internal"]] ) {
      ly[[reg]]  = round(FD$summary_annual[ get(mau)==reg & yr==year_assessment, ] [["landings"]], 1)
      lp[[reg]]  = round(FD$summary_annual[ get(mau)==reg & yr==year_previous,   ] [["landings"]], 1)
      ldt[[reg]] = round((ly[[reg]] - lp[[reg]])  / lp[[reg]] * 100, 1 )
    }

    # effort
    ep = ey = edt = list()  
    for ( reg in maus[["internal"]] ) {
      ey[[reg]]  = round(FD$summary_annual[ get(mau)==reg & yr==year_assessment, ] [["effort"]], 3)
      ep[[reg]]  = round(FD$summary_annual[ get(mau)==reg & yr==year_previous,   ] [["effort"]], 3)
      edt[[reg]] = round((ey[[reg]] - ep[[reg]])  / ep[[reg]] * 100, 1 )
    }

    # cpue
    cp = cy = cdt = list()   
    for ( reg in maus[["internal"]] ) {
      cy[[reg]]  = round(FD$summary_annual[ get(mau)==reg & yr==year_assessment, ] [["cpue"]], 2)
      cp[[reg]]  = round(FD$summary_annual[ get(mau)==reg & yr==year_previous,   ] [["cpue"]], 2)
      cdt[[reg]] = round((cy[[reg]] - cp[[reg]])  / cp[[reg]] * 100, 1 )
    }


    # tac
    tacy = tacp = list()   
    for ( reg in maus[["internal"]] ) {
      tacy[[reg]] = FD$summary_annual[ 
        get(mau)==reg & yr==year_assessment,] [["TAC"]]

      tacp[[reg]] = FD$summary_annual[ 
        get(mau)==reg & yr==year_previous  ,] [["TAC"]]
    }
      

    # shell condition (% soft)
    ccsy = ccsp = list()   
    for ( reg in maus[["internal"]] ) {
      ccsy[[reg]] = FD$shell_condition[ 
        get(mau)==reg & fishyr==year_assessment & shell %in% c(1,2), 
        sum(percent)
      ]

      ccsp[[reg]] = FD$shell_condition[ 
        get(mau)==reg & fishyr==year_previous   & shell %in% c(1,2), 
        sum(percent)
      ]
    }
      
    # fraction observed
    # here mean is used to force result as a scalar
    foby = fobp = list()   
    for ( reg in maus[["internal"]] ) {
      foby[[reg]] = FD$fraction_observed[ 
        get(mau)==reg & yr==year_assessment, 
        mean(observed_landings_pct, na.rm=TRUE)
      ]

      fobp[[reg]] = FD$fraction_observed[ 
        get(mau)==reg & yr==year_previous  , 
        mean(observed_landings_pct, na.rm=TRUE)
      ]
    }
  }


  if ( "survey" %in% todo ) {

    # recode region to selection above:
    set0 = snowcrab.db(p=p, DS="set.biologicals")
    setDT(set0)
    # check towquality .. this should always == 1
    if (length( unique( set0$towquality) ) != 1 ) print("error -- not good tows")
 
    det0 = snowcrab.db( p=p, DS="det.georeferenced" )
    setDT(det0)
    det0$fishyr = det0$yr  ## the counting routine expects this variable
  
  }


  if ( "ecosystem" %in% todo ) {
    # predator diet data
    diet_data_dir = file.path( data_loc, "data", "diets" )
    
    # assimilate the CSV data tables:
    # diet = get_feeding_data( diet_data_dir, redo=TRUE )  # if there is a data update
    diet = get_feeding_data( diet_data_dir, redo=FALSE )
    tx = taxa_to_code("snow crab")  
    # matching codes are 
    #  spec    tsn                  tx                   vern tx_index
    #1  528 172379        BENTHODESMUS           BENTHODESMUS     1659
    #2 2522  98427        CHIONOECETES SPIDER QUEEN SNOW UNID      728
    #3 2526  98428 CHIONOECETES OPILIO        SNOW CRAB QUEEN      729
    # 2 and 3 are correct

    snowcrab_predators = diet[ preyspeccd %in% c(2522, 2526), ]  # n=159 oservations out of a total of 58287 observations in db (=0.28% of all data)
    snowcrab_predators$Species = code_to_taxa(snowcrab_predators$spec)$vern
    snowcrab_predators$Predator = factor(snowcrab_predators$Species)
    
    counts = snowcrab_predators[ , .(Frequency=.N), by=.(Species)]
    setorderv(counts, "Frequency", order=-1)
    
    # species composition
    psp = speciescomposition_parameters( yrs=p$yrs, carstm_model_label="default" )
    pca = speciescomposition_db( DS="pca", p=psp )  

    pcadata = as.data.frame( pca$loadings )
    pcadata$vern = stringr::str_to_title( taxonomy.recode( from="spec", to="taxa", tolookup=rownames( pcadata ) )$vern )
    

  }


  if ( "fishery_model" %in% todo ) {

    fishery_model_results = file.path( data_loc, "fishery_model" )
        
    fm_loc = file.path( data_loc, 'fishery_model', year_assessment, model_variation )
    
    # as modelled years in fishery model can differ from iput data years, make sure  "years_model" is correct
    if (is.null(years_model)) years_model = p$fishery_model_years = 2000:year_assessment
        
    method = "logistic_discrete_historical"
    loc = file.path(data_loc, "fishery_model", year_assessment, method )

    b1north = fread( file.path(loc, "results_turing_cfanorth_bio_fishing.csv"), header=TRUE, sep=";" )
    b1south = fread( file.path(loc, "results_turing_cfasouth_bio_fishing.csv"), header=TRUE, sep=";" )
    b14x = fread( file.path(loc, "results_turing_cfa4x_bio_fishing.csv"), header=TRUE, sep=";" )

    t1 = which(years_model == year_assessment -1 )
    t0 = which(years_model == year_assessment )

    B_north = rowMeans(b1north, na.rm=TRUE )
    B_south = rowMeans(b1south, na.rm=TRUE )
    B_4x = rowMeans(b14x, na.rm=TRUE )

    B_north_sd = apply(b1north, 1, sd, na.rm=TRUE )
    B_south_sd = apply(b1south, 1, sd, na.rm=TRUE )
    B_4x_sd = apply(b14x, 1, sd, na.rm=TRUE )

    fmnorth = fread( file.path(loc, "results_turing_cfanorth_fm.csv"), header=TRUE, sep=";" )
    fmsouth = fread( file.path(loc, "results_turing_cfasouth_fm.csv"), header=TRUE, sep=";" )
    fm4x = fread( file.path(loc, "results_turing_cfa4x_fm.csv"), header=TRUE, sep=";" )

    FM_north = rowMeans(fmnorth, na.rm=TRUE )
    FM_south = rowMeans(fmsouth, na.rm=TRUE )
    FM_4x = rowMeans(fm4x, na.rm=TRUE )

    FM_north_sd = apply(fmnorth, 1, sd, na.rm=TRUE )
    FM_south_sd = apply(fmsouth, 1, sd, na.rm=TRUE )
    FM_4x_sd = apply(fm4x, 1, sd, na.rm=TRUE )


    fsnorth = fread( file.path(loc, "results_turing_cfanorth_summary.csv"), header=TRUE, sep=";" )
    fssouth = fread( file.path(loc, "results_turing_cfasouth_summary.csv"), header=TRUE, sep=";" )
    fs4x = fread( file.path(loc, "results_turing_cfa4x_summary.csv"), header=TRUE, sep=";" )

    Knorth = fsnorth[which(fsnorth$parameters=="K"),]
    Ksouth = fssouth[which(fssouth$parameters=="K"),]
    K4x = fs4x[which(fs4x$parameters=="K"),]
    
    K_north = round(Knorth[["mean"]], 2 )
    K_south = round(Ksouth[["mean"]], 2 )
    K_4x = round(K4x[["mean"]], 2 )

    K_north_sd = round(Knorth[["std"]], 2 )
    K_south_sd = round(Ksouth[["std"]], 2 )
    K_4x_sd = round(K4x[["std"]], 2 )

    rnorth = fsnorth[which(fsnorth$parameters=="r"),]
    rsouth = fssouth[which(fssouth$parameters=="r"),]
    r4x = fs4x[which(fs4x$parameters=="r"),]

    r_north = round(rnorth[["mean"]], 2 )
    r_south = round(rsouth[["mean"]], 2 )
    r_4x = round(r4x[["mean"]], 2 )

    r_north_sd = round(rnorth[["std"]], 2 )
    r_south_sd = round(rsouth[["std"]], 2 )
    r_4x_sd = round(r4x[["std"]], 2 )


    qnorth = fsnorth[which(fsnorth$parameters=="q1"),]
    qsouth = fssouth[which(fssouth$parameters=="q1"),]
    q4x = fs4x[which(fs4x$parameters=="q1"),]

    q_north = round(qnorth[["mean"]], 2 )
    q_south = round(qsouth[["mean"]], 2 )
    q_4x = round(q4x[["mean"]], 2 )

    q_north_sd = round(qnorth[["std"]], 2 )
    q_south_sd = round(qsouth[["std"]], 2 )
    q_4x_sd = round(q4x[["std"]], 2 )

  }
 
  out = as.list( environment() ) 

  read_write_fast( out, fn = fn_out )
  
  print(fn_out)

  if (return_as_list) {
    return( invisible( out ) )
  } else {
    return( invisible( list2env( out, envir) ) )
  }

}
