
  observer.db = function( DS, p=NULL, yrs=NULL, fn_root=project.datadirectory("bio.snowcrab"), 
    region="cfaall", yrs_show=7  ) {

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
     

    if (DS %in% c("bycatch_clean_data")) {
  
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
      
      obs[ , yr:=year(board_date) ]
 
      # cfa 4X has a fishing season that spans two years recode "yr" to "fishyr" to accomodate this
      cfa4x = polygon_inside(obs, aegis.polygons::polygon_internal_code("cfa4x"))
      to.offset = which( lubridate::month(obs$timestamp) >= 1 & lubridate::month(obs$timestamp) <= 7 )
      to.offset = sort(intersect(cfa4x, to.offset))
      obs$fishyr = obs$yr
      obs$fishyr[to.offset]  = obs$fishyr[to.offset] - 1

      obs = obs[obs$fishyr >= 1996 ,] # years for which observer database are good
 
      return(obs)
    }



    if (DS %in% c("bycatch_summary")) {
  
      obs = observer.db( DS="bycatch_clean_data", p=p,  yrs=yrs )  # Prepare at sea observed data
      # drop unid fish, "STONES AND ROCKS", "SEAWEED ALGAE KELP", "SNOW CRAB QUEEN", "CORAL", "SPONGES", "LEATHERBACK SEA TURTLE",  "BASKING SHARK", "SEALS", "WHALES"
      obs = obs[ !(speccd_id %in% c(90, 233, 900, 920, 8332, 8600, 9200, 9300, 9435 )), ]    # 711,412 
      obs[ , id:=paste(trip, set_no, sep="_") ]    # length(unique(obs$id))  32406
      llon = substring( paste(as.character(round(obs$lon, 2)), "00", sep=""), 2, 6)
      llat = substring( paste(as.character(round(obs$lat, 2)), "00", sep=""), 1, 5)
      obs[ , uid:=paste( cfv, llon, llat, fishyr, week(landing_date), sep="_") ]
     

      lgbk = logbook.db( DS="logbook", p=p, yrs=yrs ) # 71,194    
      names(lgbk) = tolower(names(lgbk))
      setDT(lgbk)
      lgbk[, id:=paste(licence,year,lat, lon, sep="_") ]
      llon = substring( paste(as.character(round(lgbk$lon, 2)), "00", sep=""), 2, 6)
      llat = substring( paste(as.character(round(lgbk$lat, 2)), "00", sep=""), 1, 5)
      lgbk[ , uid:=paste( cfv, llon, llat, year, week(date.landed), sep="_") ]
      lgbk$fishyr = lgbk$yr  # fishyr is coded as "yr" in logbook.db
      lgbk = lgbk[lgbk$fishyr >= 1996 ,] # years for which observer database are good


      i = polygon_inside( obs[,  c("lon", "lat")], region=region )
      j = polygon_inside( lgbk[, c("lon", "lat")], region=region )

      obs = obs[i,]
      lgbk = lgbk[i,]
 

      uid0 = unique( obs[, .(uid, fishyr, cfv, wk=week(board_date))] )

      obs = obs[!grep("NA", uid),]  # remove data with NA's in uid

      catch = data.table::dcast( obs, 
        formula =  uid ~ speccd_id, value.var="wgt", 
        fun.aggregate=sum, fill=0, drop=FALSE, na.rm=TRUE
      )  

      discards = data.table::dcast( obs, 
        formula =  uid ~ speccd_id, value.var="est_discard_wt", 
        fun.aggregate=sum, fill=0, drop=FALSE, na.rm=TRUE
      )  

      # kept weight is mostly snow crab as this is for the snow crab fishery 
      # non-zero (buy very low) kept species include: 2520, 2523, 2525, 2527, 10, 51  (probably for personal/direct use or sampling activity)
      # .. can otherwise be ignored as totals are mostly discarded
      kept = data.table::dcast( obs, 
        formula =  uid ~ speccd_id, value.var="est_kept_wt", 
        fun.aggregate=sum, fill=0, drop=FALSE, na.rm=TRUE
      )  

      # efficiency of snow crab capture
      eff = data.table(
        uid = kept$uid,
        eff = kept[["2526"]] / catch[["2526"]]
      ) 
      eff$loss = 1 - eff$eff
      eff = uid0[eff, on="uid"]

      eff_summ = eff[, .(loss=mean(loss, na.rm=TRUE), losssd=sd(loss, na.rm=TRUE)), by=.(fishyr) ]

      # catch rates kg per trap

      # number of sampling events
      effort = obs[ , .(effort=max(num_hook_haul, na.rm=TRUE)), by="uid" ] 
      effort[ !is.finite(effort), "effort"] = NA  #NA's   :1659
      
      cpue = effort[catch, on="uid"]
      cpue_uid = cpue$uid
      cpue = cpue[ , lapply(.SD, function(x) {x/effort}), .SDcols=patterns("[[:digit:]]+") ]

  #    cpue[,"2526"] = NULL
      bct  = colMeans(cpue, na.rm=TRUE)
      spec = colnames(cpue)
      
      tx = bio.taxonomy::taxonomy.recode( from="spec", to="taxa", tolookup=spec )
      species = tx$vern
      
      toshow = which(spec != "2526")  #drop as snow crab rates are very high
      spec = as.numeric(as.factor(spec[toshow]) )
    
      # fishery rates broken down by year
      lgyr = lgbk[, .(
          totaleffort = sum(effort, na.rm=TRUE),
          totallandings = sum(landings, na.rm=TRUE),
          cpue = sum(landings, na.rm=TRUE)/ sum(effort, na.rm=TRUE)
        ), 
        by=.(fishyr)
      ] 

      # return classifying variables to cpue kg/trap , standaradized to snow crab total catch rates 
      obs2 = cbind( uid=cpue_uid, cpue/cpue[["2526"]]  )  # cpu is total (kept+discard)
      obs2 = uid0[ obs2, on="uid"]

      # aggregate to year
      bycatch_table = obs2[order(fishyr), lapply(.SD, function(x) {mean(x[is.finite(x)])}), .SDcols=patterns("[[:digit:]]+"), by=.(fishyr) ]  

      bycatch_table = lgyr[ bycatch_table, on="fishyr"]
      bycatch_table = eff_summ[ bycatch_table, on="fishyr"]

      # drop no data years 
      bycatch_table = bycatch_table[ fishyr>2000, ]  

      # bycatch_table = zapsmall(bycatch_table)
      yrs = bycatch_table$fishyr
      bycatch_table$fishyr = NULL

      cpue_fraction = colMeans(bycatch_table[, .SD, .SDcols=patterns("[[:digit:]]+") ], na.rm=TRUE)

      # discards = totallandings * loss / (1-loss)
      # discard + kept = totallandings + totallandings * loss / (1-loss)
      #                = totallandings * ( (1 / (1-loss) )

      # compute by catch as a fraction of snow crab landings
      bycatch_table = bycatch_table[, lapply(.SD, function(x) {x*totallandings* ( 1 / (1-loss) ) }), .SDcols=patterns("[[:digit:]]+") ]

      bycatch_table[["totaleffort"]] = NULL
      bycatch_table[["totallandings"]] = NULL
  
      specs = names(bycatch_table)  
      tx = bio.taxonomy::taxonomy.recode( from="spec", to="taxa", tolookup=specs )

      bycatch_table = as.data.table( t(bycatch_table) )
      names( bycatch_table) = as.character( yrs )

      bycatch_table = zapsmall( bycatch_table, digits=9)
      # bycatch_table = zapsmall( bycatch_table )
      bycatch_table$spec = specs
      bycatch_table$cpue_fraction = cpue_fraction # mean fraction of landing per year (weight/weight)
      bycatch_table$species = str_to_title(tx$vern)
      bycatch_table$taxa = tx$tx

      # snow crab kept to this point as a sosuble check on computations
      yrss = year.assessment - (yrs_show-1):0

      to_show = c( "species", as.character( yrss ) )
      to_keep = setdiff( which(cpue_fraction >1e-9), which(names(cpue_fraction) =="2526" ))

      # check if years missing
      missing_years = setdiff( to_show, names(bycatch_table))
      j = length(missing_years)
      if ( j > 0) {
        for (i in 1:j) {
          vn = paste( "missing_years[", i, "]", sep="" )
          bycatch_table[, eval(parse(text=vn)) ] =NA 
        }
      }
      
      out = bycatch_table[to_keep, ..to_show] 
      i = 2:ncol(out)
      sums = data.table( species="Total", t(colSums( bycatch_table[to_keep, ..to_show][,..i] ) ) )
      out  = rbind(out, sums )

      landings = t(lgyr[ fishyr %in% yrss ,][["totallandings"]])
      landings = data.table( "Landings/Débarquements (kg)", landings )
      names(landings ) = c("species", yrss )
      out  = rbind(out, landings )

      effort2 = t(lgyr[ fishyr %in% yrss ,][["totaleffort"]])
      effort2 = data.table( "Effort (trap hauls/casiers levés)", effort2 ) 
      names(effort2 ) = c("species", yrss )
      out  = rbind(out, effort2 )


      coverage2 = uid0[ catch, on="uid"]
      coverage2$snowcrab = coverage2[["2526"]]
      coverage2 = coverage2[ fishyr %in% yrss, .(catch_sum=sum(snowcrab, na.rm=TRUE)), by=.(fishyr)]
      missing_years = setdiff( to_show, c("species", coverage2$fishyr) )
      j = length(missing_years)
      if ( j > 0) {
        for (i in 1:j) {
          coverage2 = rbind( coverage2, cbind( fishyr=as.numeric(missing_years[i]), catch_sum=NA ))
        }
      }
      coverage2 = coverage2[order(fishyr),][["catch_sum"]]
      coverage2 = data.table( "At sea observed catch / Débarquements observés en mer (kg)", t(coverage2) ) 
      names(coverage2 ) = c("species", yrss )
      out  = rbind(out, coverage2 )


      coverage = uid0[ effort, on="uid"]
      coverage = coverage[ fishyr %in% yrss, .(effort_sum=sum(effort, na.rm=TRUE)), by=.(fishyr)]
      missing_years = setdiff( to_show, c("species", coverage$fishyr) )
      j = length(missing_years)
      if ( j > 0) {
        for (i in 1:j) {
          coverage = rbind( coverage, cbind( fishyr=as.numeric(missing_years[i]), effort_sum=NA ))
        }
      }
      coverage = coverage[order(fishyr),][["effort_sum"]]
      coverage = data.table( "At sea observed effort / efforts observés en mer (trap hauls/casiers levés)", t(coverage) ) 
      names(coverage ) = c("species", yrss )
      out  = rbind(out, coverage )


      out$"Average/Moyen" = rowMeans(out[,.SD,.SDcols=patterns("[[:digit:]]+")], na.rm=TRUE)

      out = list( 
        eff_summ=eff_summ,
        bct = bct[toshow],
        bct_sc = round(bct[["2526"]], 1 ),
        spec = spec,
        species = str_to_title(species[toshow]),
        cpue_fraction= cpue_fraction, 
        bycatch_table = out
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
      odb$tripset = paste( odb$trip, odb$set_no, sep="~")
      odb$cpue.kg.trap = ( odb$totmass*1000)/odb$num_hook_haul
      
      save(odb, file=fn, compress=T)

      return( "complete" )
    }


  }


