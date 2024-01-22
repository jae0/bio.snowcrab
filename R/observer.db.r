
  observer.db = function( DS, p=NULL, yrs=NULL, fn_root=project.datadirectory("bio.snowcrab") ) {

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


