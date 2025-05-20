
  observer.db = function( DS, p=NULL, yrs=NULL, fn_root=project.datadirectory("bio.snowcrab"), 
    region="cfaall", yrs_show=5  ) {

    # sex codes
    male = 0
    female = 1
    sex.unknown = 2

    # maturity codes
    immature = 0
    mature = 1
    mat.unknown = 2


    if (DS %in% c("rawdata.redo", "rawdata") ) {

      fn.loc =  file.path( fn_root, "data", "observer", "datadump" )
			dir.create( fn.loc, recursive = TRUE, showWarnings = FALSE )

			if (DS=="rawdata") {
				out = NULL
				for ( YR in yrs ) {
					fny = file.path( fn.loc, paste( YR, "rdata", sep="."))
					if (file.exists(fny)) {
						odb = NULL
            load (fny)
            setDT(odb)
            if (!is.null(odb)) {
              out = try( rbind( out, odb ) )
              if (inherits(out, "try-error" )) {
                stop( "problem: inconsistent number of variables", fny)
              }
            }
          }
				}
				return (out)
			}

      # for the full list of tables:
      # tbls = sqlTables(connect)
      # gs.tables = tbls[ which(tbls[,2] == "GROUNDFISH"),]
      # print(gs.tables)
      con = ROracle::dbConnect(DBI::dbDriver("Oracle"),dbname=oracle.snowcrab.server , username=oracle.snowcrab.user, password=oracle.snowcrab.password, believeNRows=F)

      for ( YR in yrs ) {
        fny = file.path( fn.loc, paste( YR, "rdata", sep="."))
        odbq = paste(
          "SELECT s.LATITUDE, s.LONGITUDE, s.LANDING_DATE, s.SET_NO, s.PRODCD_ID, s.SETCD_ID, s.EST_CATCH, s.EST_KEPT_WT," ,
          "s.NUM_HOOK_HAUL, d.BOARD_DATE, d.FISH_NO, d.SEXCD_ID, d.FISH_LENGTH, " ,
          "d.FEMALE_ABDOMEN, d.CHELA_HEIGHT, d.SHELLCOND_CD, d.DUROMETRE, d.TRIP_ID, d.TRIP  " ,
          "FROM SNOWCRAB.SNCRABDETAILS_OBS d, SNOWCRAB.SNCRABSETS_OBS s " ,
          "WHERE d.TRIP_ID = s.TRIP_ID  " ,
          "AND d.SET_NO = s.SET_NO  " ,
          "AND s.SETCD_ID = 1", ## to keep trips for actual fishery activity (and not CP activity etc) 
          "AND d.FISH_NO Is Not Null" ,
          "AND EXTRACT(YEAR from d.BOARD_DATE) = ", YR )
        odb = NULL
        odb = ROracle::dbGetQuery(con, odbq )
        save( odb, file=fny, compress=T)
        gc()  # garbage collection
        print(YR)
      }
      ROracle::dbDisconnect(con)
      return (yrs)

    }



    if (DS %in% c("bycatch.redo", "bycatch") ) {

      fn.loc =  file.path( fn_root, "data", "observer", "bycatch" )
			dir.create( fn.loc, recursive = TRUE, showWarnings = FALSE )

			if (DS=="bycatch") {
        out = NULL
				for ( YR in yrs ) {
					fny = file.path( fn.loc, paste( YR, "rdata", sep="."))
					if (file.exists(fny)) {
						odb = NULL
            load (fny)
            setDT(odb)
            if (!is.null(odb)) {
              out = try( rbind( out, odb ) )
              if (inherits(out, "try-error" )) {
                stop( "problem: inconsistent number of variables", fny)
              }
            }
          }
				}
				return (out)
			}

      con = ROracle::dbConnect(DBI::dbDriver("Oracle"),dbname=oracle.snowcrab.server , username=oracle.snowcrab.user, password=oracle.snowcrab.password, believeNRows=F)

      for ( YR in yrs ) {
        fny = file.path( fn.loc, paste( YR, "rdata", sep="."))
        odbq = paste(
          "SELECT trip.trip_id, trip.trip, trip.board_date, trip.landing_date, st.set_no, vess.vessel_name, vess.license_no, vess.cfv,",
          "  isp.LATITUDE, isp.LONGITUDE, isp.DEPTH, isp.WATER_TEMPERATURE, ", 
          "  sc.common, st.comarea_id, st.nafarea_id,  ",
          "  ca.speccd_id, ca.est_num_caught,  ca.est_kept_wt, ca.est_discard_wt, st.NUM_HOOK_HAUL,",
          "  fish.fish_no , fish.fish_length, fish.fish_weight ", 
          "FROM istrips trip, isgears gr, isfishsets st, iscatches ca, isfish fish, issetprofile isp, isvessels vess,  isobservercodes o, isspeciescodes sc ", 
          "WHERE trip.tripcd_id = 2509",
          "AND vess.vess_id = trip.vess_id",
          "AND vess.license_no = trip.license_no",
          "AND o.obscd_id = trip.obscd_id",
          "AND trip.trip_id = gr.trip_Id",
          "AND st.fishset_id = isp.fishset_id ",
          "AND isp.fishset_id = ca.fishset_id ",
          "AND ca.speccd_id = sc.speccd_id ",
          "AND trip.tripcd_id=2509 ",
          "AND st.specscd_id=2526 ",
          "AND isp.pntcd_id=4 ",
          "AND (trip.trip_id = st.trip_id AND gr.gear_id = st.gear_id)",
          "AND ca.catch_id = fish.catch_id(+)",
          "AND EXTRACT(YEAR from trip.board_date) = ", YR )

        odb = NULL
        odb = ROracle::dbGetQuery(con, odbq )
        save( odb, file=fny, compress=T)
        gc()  # garbage collection
        print(YR)
      }
      ROracle::dbDisconnect(con)
      return (yrs)

    } 
     

    if (DS %in% c("bycatch_clean_data", "bycatch_clean_data.redo")) {
      
      fn = file.path( project.datadirectory("bio.snowcrab"), "data", "observer", "odb_bycatch.rdata" )
      if (DS=="bycatch_clean_data") {
        load( fn )
        return(obs)
      }

      obs = observer.db( DS="bycatch", p=p, yrs=yrs )  # 711,413 in 2023
      names(obs) = tolower(names(obs))
      setDT(obs)

      # QA/QC 
      obs$lon = -obs$longitude
      obs$lat =  obs$latitude
      
      obs$longitude = NULL
      obs$latitude = NULL

      obs[ !is.finite(est_kept_wt), "est_kept_wt" ] = 0
      obs[ !is.finite(est_discard_wt), "est_discard_wt" ] = 0
      # obs[ !is.finite(num_hook_haul), "num_hook_haul" ] = 5 # this is protocal (assume when not recorded it is 5)
      obs[ , wgt:=est_kept_wt + est_discard_wt ]
      
      obs[ num_hook_haul > 10, "num_hook_haul"] = NA  # errors likely 
      obs[ est_discard_wt > 1000, "est_discard_wt"] = NA  # errors likely .. these are estimated 

      obs[ , cpue:= est_discard_wt /num_hook_haul ]

      obs[ cpue > 150, "cpue"] = NA  # errors likely as capture higher than 150kg / trap unlikely .. note these are cpu of all positive valued catches .. zero-valued are not records

      obs[ !is.finite(board_date), "board_date" ] = NA
      
      obs[ , yr:=year(board_date) ]
 
      # cfa 4X has a fishing season that spans two years recode "yr" to "fishyr" to accomodate this
      cfa4x = polygon_inside(obs, aegis.polygons::polygon_internal_code("cfa4x"))
      mnth =  lubridate::month(obs$board_date)

      to.offset = which(mnth >= 1 & mnth <= 7 )

      to.offset = sort(intersect(cfa4x, to.offset))
      obs$fishyr = obs$yr
      obs$fishyr[to.offset]  = obs$fishyr[to.offset] - 1

      obs = obs[obs$fishyr >= 1996 ,] # years for which observer database are good
      obs[ , uid:=paste(trip, set_no, sep="_") ]    # length(unique(obs$id))  32406

      save(obs, file=fn, compress=TRUE)

      return(fn)
    }



    if (DS %in% c("bycatch_summary")) {
  
      # fishery rates broken down by year (kg, trap hauls)
      lgbk = logbook.db( DS="logbook", p=p, yrs=yrs ) # 71,194    
      names(lgbk) = tolower(names(lgbk))
      setDT(lgbk)
      lgbk$fishyr = lgbk$yr  # fishyr is coded as "yr" in logbook.db
      lgbk = lgbk[lgbk$fishyr %in% c(yrs, min(yrs)-1, max(yrs)+1) ,]  
      j = polygon_inside( lgbk[, c("lon", "lat")], region=region )
      lgbk = lgbk[j,]
      lgyr = lgbk[, .(
          totaleffort = sum(effort, na.rm=TRUE),
          totallandings = sum(landings, na.rm=TRUE),
          cpue = sum(landings, na.rm=TRUE)/ sum(effort, na.rm=TRUE)
        ), 
        by=.(fishyr)
      ] 
      lgbk = NULL

      # observer data .. extra care needed as there are duplicated records, etc
      obs = observer.db( DS="bycatch_clean_data", p=p,  yrs=yrs )  # Prepare at sea observed data
      i = polygon_inside( obs[,  c("lon", "lat")], region=region )
      
      if (length(i) == 0) return(NULL)

      oss = obs = obs[i,]
      
      # drop unid fish, "STONES AND ROCKS", "SEAWEED ALGAE KELP", "SNOW CRAB QUEEN", "CORAL", "SPONGES", "LEATHERBACK SEA TURTLE",  "BASKING SHARK", "SEALS", "WHALES"
      # large animals will make etmates meaningless if using weight ... so drop whales, leatherback etc.
      obs = obs[ !(speccd_id %in% c(90, 233, 900, 920, 8332, 8600, 9200, 9300, 9435 )), ]    # 711,412 
      obs = obs[!grep("NA", uid),]  # remove data with NA's in uid
      uid0 = unique( obs[, .(uid, fishyr, cfv, wk=week(board_date))] )

      # sampling effort and catch .. unique as data entered seems sometimes not aggregated
      # also the SQL pull should be roken down to get the correct parts without the outer join
      observed_effort = obs[, .(no_traps=unique(num_hook_haul)), by=.(uid)][, .(no_traps=sum(no_traps, na.rm=TRUE)), by=.(uid)]
      observed_effort[ no_traps==0, "no_traps"] = NA  # 5  ?override with standard sampling procotol expectation
      observed_catch = obs[, .(wgt=unique(wgt)), by=.(uid)][, .(wgt=sum(wgt, na.rm=TRUE)), by=.(uid)] 
      observed = observed_effort[observed_catch, on=.(uid) ]
  
      # same but broken down by uid and species
      observed_catch2 = obs[, .(wgt=unique(wgt)), by=.(uid, speccd_id)][, .(wgt=sum(wgt, na.rm=TRUE)), by=.(uid, speccd_id)] 
      observed_discard2 = obs[, .(est_discard_wt=unique(est_discard_wt)), by=.(uid, speccd_id)][, .(est_discard_wt=sum(est_discard_wt, na.rm=TRUE)), by=.(uid, speccd_id)] 
      observed_kept2 = obs[, .(est_kept_wt=unique(est_kept_wt)), by=.(uid, speccd_id)][, .(est_kept_wt=sum(est_kept_wt, na.rm=TRUE)), by=.(uid, speccd_id)] 
      # kept weight is mostly snow crab as this is for the snow crab fishery 
      # non-zero (buy very low) kept species include: 2520, 2523, 2525, 2527, 10, 51  (probably for personal/direct use or sampling activity)
      # .. can otherwise be ignored as totals are mostly discarded
      observed2 = observed_catch2[observed_discard2, on=.(uid, speccd_id) ]
      observed2 = observed2[observed_kept2, on=.(uid, speccd_id) ]

      # add set/trip stats
      observed2 = observed[ observed2, on=.(uid) ]
      setnames( observed2, "i.wgt", "est_catch_wt" )
      setnames( observed2, "wgt", "est_total_catch_wt" )

      observed2[ , discard_rate := est_discard_wt/est_catch_wt ]
      observed2[ , fraction_of_catch := est_catch_wt/est_total_catch_wt ]
      observed2[ , cpue_total := est_catch_wt/no_traps ]
      observed2[ , cpue_kept := est_catch_wt/no_traps ]
 
      # efficiency of snow crab capture
      eff = observed2[ speccd_id==2526, ] 
      eff = uid0[eff, on="uid"]
      eff_summ = eff[, .(discard_rate=mean(discard_rate, na.rm=TRUE), discard_rate_sd=sd(discard_rate, na.rm=TRUE)), by=.(fishyr) ]

      # catch rates kg per trap
      # number of sampling events
      cpue = dcast( observed2, uid~speccd_id, value.var="cpue_total", 
         fun.aggregate=mean, fill=0, drop=FALSE, na.rm=TRUE)      
       # cpue = uid0[catch, on="uid"]
      
      cpue_uid = cpue$uid
   
      cpue = cpue[ , lapply(.SD, function(x) {ifelse(is.nan(x), NA, x) } ), 
        .SDcols=patterns("[[:digit:]]+") ] # reset NaNs

      bct  = colMeans(cpue, na.rm=TRUE)
      specid = colnames(cpue)
      tx = bio.taxonomy::taxonomy.recode( from="spec", to="taxa", tolookup=specid )
      species = tx$vern
      
      toshow = which(specid != "2526")  #drop as snow crab rates are very high
      spec = as.numeric(as.factor(specid[toshow]) )


      # direct computation from effort by year/spec
      bct_effort = cbind(data.table(uid=cpue_uid), cpue)
      bct_effort = bct_effort[uid0, on=.(uid)]
      bct_effort = bct_effort[ order(fishyr), lapply(.SD, function(x) {mean(x[is.finite(x)])}), .SDcols=patterns("[[:digit:]]+"), by=.(fishyr) ]  
      bct_effort = bct_effort[ , lapply(.SD, function(x) {ifelse(is.nan(x), NA, x) } ) ]  # reset NaNs
      bct_effort = lgyr[bct_effort, on=.(fishyr)]
      bct_effort = bct_effort[ , lapply(.SD, function(x) {x*totaleffort}), .SDcols=patterns("[[:digit:]]+"), by=.(fishyr) ]  
      
      bct_effort_catchmean = colMeans(bct_effort[, .SD, .SDcols=patterns("[[:digit:]]+") ], na.rm=TRUE)
      
      specs = names(bct_effort) [-1]
      tx = bio.taxonomy::taxonomy.recode( from="spec", to="taxa", tolookup=specs )
      
      bct_effort = bct_effort[ data.table(fishyr=yrs) , on=.(fishyr)]
      bct_effort$fishyr = NULL

      bct_effort = as.data.table( t(bct_effort) )

      names( bct_effort) = as.character( yrs )
      bct_effort = zapsmall( bct_effort, digits=9)
      bct_effort$spec = specs
      bct_effort$species = stringr::str_to_title(tx$vern)
      bct_effort$taxa = tx$tx

      # snow crab kept to this point as a sosuble check on computations
      yrss = p$year.assessment - (yrs_show-1):0
      yrsa = 1999:p$year.assessment 

      to_show = c( "species", as.character( yrss ) )
      to_use = c( "species", as.character( yrsa ) )
      to_keep = setdiff( which(bct_effort_catchmean > 1), which(names(bct_effort_catchmean) =="2526" ))

      # check if years missing
      missing_years = setdiff( to_use, names(bct_effort))
      j = length(missing_years)
      if ( j > 0) {
        for (i in 1:j) {
          vn = paste( "missing_years[", i, "]", sep="" )
          bct_effort[, eval(parse(text=vn)) ] = NA 
        }
      } 

      bctabe = bct_effort[to_keep, ..to_use] 
      i = 2:ncol(bctabe)
      sums = data.table( species="Total", t(colSums( bct_effort[to_keep, ..to_use][,..i] ) ) )
      bctabe  = rbind(bctabe, sums )

      lds = lgyr[ data.table(fishyr=yrsa) , on=.(fishyr)]
      landings = round(t(lds[["totallandings"]]))
      landings = data.table( "Landings/Débarquements (kg)", landings )
      names(landings ) = c("species", yrsa )
      bctabe  = rbind(bctabe, landings )

      effort2 = t(lds[["totaleffort"]])
      effort2 = data.table( "Effort (trap hauls/casiers levés)", effort2 ) 
      names(effort2 ) = c("species", yrsa )
      bctabe  = rbind(bctabe, effort2 )
  
      obs_summ = observed[ uid0, on=.(uid) ]

      coverage = obs_summ[order(fishyr), .(wgt=sum(wgt, na.rm=TRUE)), by=.(fishyr)]
      coverage = coverage[ data.table(fishyr=yrsa) , on=.(fishyr)][["wgt"]]
      coverage = data.table( "At sea observed catch \n/Débarquements observés en mer (kg)", t(coverage) ) 
      names(coverage ) = c("species", yrsa )
      bctabe  = rbind(bctabe, coverage) 

      efforts = obs_summ[order(fishyr), .(no_traps=sum(no_traps, na.rm=TRUE)), by=.(fishyr)]
      efforts = efforts[ data.table(fishyr=yrsa) , on=.(fishyr)][["no_traps"]]
      efforts = data.table( "At sea observed effort (trap hauls) \n/Efforts observés en mer (casiers levés)", t(efforts) ) 
      names(efforts ) = c("species", yrsa )
      bctabe  = rbind(bctabe, efforts) 
 
      bctabe$"Average/Moyen" = round(rowMeans(bctabe[,.SD,.SDcols=patterns("[[:digit:]]+")], na.rm=TRUE) )

      bctabe[,2:ncol(bctabe)] = round(bctabe[,2:ncol(bctabe)])

      toprint = c("species", yrss, "Average/Moyen")
      bctabe = bctabe[,..toprint]



      # return classifying variables to cpue kg/trap , standaradized to snow crab total catch rates 
      obs2 = cbind( uid=cpue_uid, cpue/cpue[["2526"]]  )  # cpu is total (kept+discard)
      obs2 = uid0[ obs2, on="uid"]
      # aggregate to year
      bct_catch = obs2[order(fishyr), lapply(.SD, function(x) {mean(x[is.finite(x)])}), .SDcols=patterns("[[:digit:]]+"), by=.(fishyr) ]  
      bct_catch = bct_catch[ , lapply(.SD, function(x) {ifelse(is.nan(x), NA, x) } ) ]  # reset NaNs
      bct_catch = lgyr[ bct_catch, on=.(fishyr)]
      bct_catch = eff_summ[ bct_catch, on=.(fishyr)]
      bct_catch = bct_catch[ data.table(fishyr=yrs) , on=.(fishyr)]
      yrs = bct_catch$fishyr
      bct_catch$fishyr = NULL
      cpue_fraction = colMeans(bct_catch[, .SD, .SDcols=patterns("[[:digit:]]+") ], na.rm=TRUE)
      # compute by catch as a fraction of snow crab landings
      bct_catch = bct_catch[, lapply(.SD, function(x) {x*totallandings* ( 1 / (1-discard_rate) ) }), .SDcols=patterns("[[:digit:]]+") ]
      bct_catch = bct_catch[, lapply(.SD, function(x) ifelse(is.nan(x), NA, x))]  #NaN's behave differently ..

      specs = names(bct_catch)  
      tx = bio.taxonomy::taxonomy.recode( from="spec", to="taxa", tolookup=specs )
      bct_catch = as.data.table( t(bct_catch) )
      names( bct_catch) = as.character( yrs )
      bct_catch = zapsmall( bct_catch, digits=9)
      bct_catch$spec = specs
      bct_catch$cpue_fraction = cpue_fraction # mean fraction of landing per year (weight/weight)
      bct_catch$species = stringr::str_to_title(tx$vern)
      bct_catch$taxa = tx$tx

      # snow crab kept to this point as a sosuble check on computations
      yrss = p$year.assessment - (yrs_show-1):0
      yrsa = 1999:p$year.assessment 

      to_show = c( "species", as.character( yrss ) )
      to_use = c( "species", as.character( yrsa ) )
      to_keep = setdiff( which(cpue_fraction >1e-9), which(names(cpue_fraction) =="2526" ))

      # check if years missing
      missing_years = setdiff( to_use, names(bct_catch))
      j = length(missing_years)
      if ( j > 0) {
        for (i in 1:j) {
          vn = paste( "missing_years[", i, "]", sep="" )
          bct_catch[, eval(parse(text=vn)) ] = NA 
        }
      }
      
      bctabc = bct_catch[to_keep, ..to_use] 
      i = 2:ncol(bctabc)
      sums = data.table( species="Total", t(colSums( bct_catch[to_keep, ..to_use][,..i] ) ) )
      bctabc  = rbind(bctabc, sums )

      lds = lgyr[ data.table(fishyr=yrsa) , on=.(fishyr)]
      landings = round(t(lds[["totallandings"]]))
      landings = data.table( "Landings/Débarquements (kg)", landings )
      names(landings ) = c("species", yrsa )
      bctabc  = rbind(bctabc, landings )

      effort2 = t(lds[["totaleffort"]])
      effort2 = data.table( "Effort (trap hauls/casiers levés)", effort2 ) 
      names(effort2 ) = c("species", yrsa )
      bctabc  = rbind(bctabc, effort2 )
  
      obs_summ = observed[ uid0, on=.(uid) ]

      coverage = obs_summ[order(fishyr), .(wgt=sum(wgt, na.rm=TRUE)), by=.(fishyr)]
      coverage = coverage[ data.table(fishyr=yrsa) , on=.(fishyr)][["wgt"]]
      coverage = data.table( "At sea observed catch \n/Débarquements observés en mer (kg)", t(coverage) ) 
      names(coverage ) = c("species", yrsa )
      bctabc  = rbind(bctabc, coverage) 

      efforts = obs_summ[order(fishyr), .(no_traps=sum(no_traps, na.rm=TRUE)), by=.(fishyr)]
      efforts = efforts[ data.table(fishyr=yrsa) , on=.(fishyr)][["no_traps"]]
      efforts = data.table( "At sea observed effort (trap hauls) \n/Efforts observés en mer (casiers levés)", t(efforts) ) 
      names(efforts ) = c("species", yrsa )
      bctabc  = rbind(bctabc, efforts) 
 
      bctabc$"Average/Moyen" = round(rowMeans(bctabc[,.SD,.SDcols=patterns("[[:digit:]]+")], na.rm=TRUE))
 
      bctabc[,2:ncol(bctabc)] = round(bctabc[,2:ncol(bctabc)])

      toprint = c("species", yrss, "Average/Moyen")
      bctabc = bctabc[,..toprint]

 
      out = list( 
        oss=oss,
        eff_summ=eff_summ,
        bct = bct[toshow],
        bct_sc = round(bct[["2526"]], 1 ),
        spec = spec,  # as factor
        specid = specid,
        species = stringr::str_to_title(species[toshow]),
        cpue_fraction= cpue_fraction, 
        bycatch_table_effort = bctabe,
        bycatch_table_catch = bctabc
      )

      return(out)
    }


    # ---------------------

    if (DS %in% c("odb", "odb.redo")) {
      fn = file.path( project.datadirectory("bio.snowcrab"), "data", "observer", "odb.rdata" )
      if (DS=="odb") {
        load( fn )
        return(odb)
      }

      mod1 = allometry.snowcrab ( "cw.mass", "male")
      mod2 = allometry.snowcrab ( "chela.mass", "male")
      mod3 = allometry.snowcrab ( "cw.chela.mat", "male")

      odb = observer.db( DS="rawdata", yrs=1996:(p$year.assessment + 1)) # add one to get new year fraction for 4x
      names(odb) = rename.bio.snowcrab.variables( names(odb) )

      i.m = which( odb$sex==1)
      i.f = which( odb$sex==2)
      i.o = setdiff( 1:nrow( odb ), c(i.m, i.f) )
      odb$sex[ i.m ] = male
      odb$sex[ i.f ] = female
      odb$sex[ i.o ] = sex.unknown

      odb$cw[ which(odb$cw < 40 | odb$cw > 185)  ] = NA
      odb$log.cw = log( odb$cw )
      odb$mass =NA
      odb$log.mass = predict( mod1, odb )
      kk = intersect( which( !is.finite( odb$mass ) ), i.m )
      odb$mass[kk] = exp( odb$log.mass[kk] )

      kk = intersect( which( !is.finite( odb$mass ) ), i.m )
      odb$log.chela = log( odb$chela )
      log.mass2 = predict( mod2, odb )
      odb$mass[kk] = exp( log.mass2[kk] )
      rm (mod1, mod2)


      odb$mat_tmp = predict( mod3, odb, type="response" )
      odb$mat = NA
      odb$mat [ which(odb$mat_tmp <0.5 & odb$sex==male) ] = immature
      odb$mat [ which(odb$mat_tmp >=0.5 & odb$sex==male) ] = mature

      odb$log.cw = NULL
      odb$log.chela = NULL
      odb$log.mass = NULL
      odb$mat_tmp = NULL

      odb$lon = -odb$lon
      odb$timestamp =  odb$sdate
      odb$yr = lubridate::year(odb$timestamp)

      # cfa 4X has a fishing season that spans two years recode "yr" to "fishyr" to accomodate this
      cfa4x = polygon_inside(odb, aegis.polygons::polygon_internal_code("cfa4x"))
      to.offset = which( lubridate::month(odb$timestamp) >= 1 & lubridate::month(odb$timestamp) <= 7 )
      to.offset = sort(intersect(cfa4x, to.offset))
      odb$fishyr = odb$yr
      odb$fishyr[to.offset]  = odb$fishyr[to.offset] - 1

      odb = odb[odb$fishyr >= 1996 ,] # years for which observer database are good
      #  odb$cw[odb$cw>175] = NA
      odb$tripset = paste( odb$trip, odb$set, sep="~")
      odb$cpue.kg.trap = ( odb$totmass*1000)/odb$num_hook_haul
      
      save(odb, file=fn, compress=T)

      return( "complete" )
    }


  }


