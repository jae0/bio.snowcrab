

  logbook.db = function( DS="logbook", p=NULL, yrs=NULL, fn_root=project.datadirectory("bio.snowcrab"), 
    region=NULL, sppoly=NULL, redo=FALSE ) {
 

		if (DS %in% c("rawdata.logbook", "rawdata.logbook.redo")) {

      fn.loc =  file.path( fn_root, "data", "logbook", "datadump" )

      dir.create( fn.loc, recursive = TRUE, showWarnings = FALSE )

			if (DS=="rawdata.logbook") {
				out = NULL
				for ( YR in yrs ) {
					fny = file.path( fn.loc, paste( YR, "rdz", sep="."))
					if (file.exists(fny)) {
						out = rbind( out, read_write_fast(fny) )
					}
				}
				return (out)
			}

			con=ROracle::dbConnect(DBI::dbDriver("Oracle"), dbname=oracle.snowcrab.server , 
        username=oracle.snowcrab.user, password=oracle.snowcrab.password, believeNRows=FALSE )

      logbook_null = read_write_fast( fn=file.path( fn.loc, paste( 1996,"rdz", sep=".")) )  # small empty file 

			for ( YR in yrs ) {
				fn = file.path( fn.loc, paste( YR,"rdz", sep="."))
				sqlquery = paste(
					"SELECT * from marfissci.marfis_crab ",
					"where target_spc=705",
					"AND EXTRACT(YEAR from DATE_LANDED) = ", YR )
				logbook = NULL
				#in following line replaced sqlQuery (Rrawdata) with  dbGetQuery (ROracle)
				logbook = try( ROracle::dbGetQuery(con, sqlquery ), silent=TRUE)
        if (inherits(logbook, "try-error")) {
          logbook = logbook_null
          message( YR, " ...No data? " )
        } else {
          print(YR)
        }
        read_write_fast( logbook, fn=fn )
				gc()  # garbage collection
			}
      ROracle::dbDisconnect(con)
      return (yrs)

		}


    # -------------------------


    if (DS %in% c("rawdata.licence.redo", "rawdata.licence" ) ) {

      fn.loc =  file.path( fn_root, "data", "logbook"  )

 			dir.create( fn.loc, recursive = TRUE, showWarnings = FALSE )

      filename.licence = file.path( fn.loc, "lic.datadump.rdz" )

      if (DS=="rawdata.licence") {
        lic = read_write_fast(filename.licence)
        return (lic)
      }

      con=ROracle::dbConnect(DBI::dbDriver("Oracle"),dbname=oracle.snowcrab.server , username=oracle.snowcrab.user, password=oracle.snowcrab.password, believeNRows=F)
      lic = ROracle::dbGetQuery(con, "select * from marfissci.licence_areas")
      read_write_fast(lic, fn=filename.licence )
      ROracle::dbDisconnect(con)
		}


    if (DS %in% c("rawdata.areas.redo", "rawdata.areas" ) ) {

      fn.loc =  file.path( fn_root, "data", "logbook"  )

			dir.create( fn.loc, recursive = TRUE, showWarnings = FALSE )

      filename.areas = file.path( fn.loc, "areas.datadump.rdz" )

      if (DS=="rawdata.areas") {
        areas = read_write_fast(filename.areas)
        return (areas)
      }

      con=ROracle::dbConnect(DBI::dbDriver("Oracle"),dbname=oracle.snowcrab.server , 
        username=oracle.snowcrab.user, password=oracle.snowcrab.password, believeNRows=F)
      areas = ROracle::dbGetQuery(con, "select * from marfissci.areas")
      read_write_fast(areas, fn=filename.areas)
      ROracle::dbDisconnect(con)
      return ("Complete")

    }


  # -------

    if (DS %in% c("logbook.filtered.positions", "logbook.filtered.positions.redo")) {

      # exclude data that have positions that are incorrect

      filename = file.path( project.datadirectory("bio.snowcrab"), "data",
         "logbook", "logbook.filtered.positions.rdz" )

      if (DS=="logbook.filtered.positions") {
        lgbk = read_write_fast( filename )
        return(lgbk)
      }

      lgbk = logbook.db( DS="logbook" )
      
  
      # alter position data that are strange in location .. land
  
      not_in_martimes = setdiff(1:nrow(lgbk), polygon_inside( lgbk[, c("lon", "lat")], region="cfaall") )

      too_deep = setdiff(1:nrow(lgbk), polygon_inside( lgbk[, c("lon", "lat")], region="isobath1000m") )

      no_position = which(is.na( lgbk$lon ) | is.na( lgbk$lat ) ) 


      # use the match/map syntax in bathymetry and filter out shallow sets .. < 10 m? TODO
      # keep those in the domain and deeper than depth=10 m
 
      pL = aegis.bathymetry::bathymetry_parameters( project_class="stmv"  )
      LUT= aegis_survey_lookuptable( aegis_project="bathymetry", 
        project_class="core", DS="aggregated_data", pL=pL )

      z = aegis_lookup( pL=pL,, LOCS=lgbk[, c("lon", "lat")], LUT=LUT,
        output_format="points" , 
        space_resolution=p$pres*2, variable_name="z.mean"  ) # core=="rawdata"

      too_shallow = which( z < 10 ) # negative = above land
       

      # check position to see if they are within the polygon designated
      icfa4x = polygon_inside( lgbk[, c("lon", "lat")], "cfa4x")
      icfanorth = polygon_inside( lgbk[, c("lon", "lat")], "cfanorth")
      icfasouth = polygon_inside( lgbk[, c("lon", "lat")], "cfasouth")

      lgbk$tag ="ok"


      bad_area = NULL
      bad_area = c(bad_area, icfa4x[which(lgbk$region[icfa4x] != "cfa4x" )] )
      bad_area = c(bad_area, icfanorth[which(lgbk$region[icfanorth] != "cfanorth" )] )
      bad_area = c(bad_area, icfasouth[which(lgbk$region[icfasouth] != "cfasouth" )] )
 
      bad = unique( c( bad_area, too_deep, too_shallow, not_in_martimes, no_position )) 


      if (0) {

        plot(plat~plon, lgbk[, ], pch=".", xlim=c(200,1100), ylim=c(4600, 5300))
        points(plat~plon, lgbk[too_deep, ], col="red" , pch=19, cex=0.8)
        points(plat~plon, lgbk[not_in_martimes, ], col="green" , pch=19, cex=0.8)
        points(plat~plon, lgbk[too_shallow, ], col="blue" , pch=19, cex=0.8)
        points(plat~plon, lgbk[bad_area, ], col="yellow" , pch=19, cex=0.8)
        points(plat~plon, lgbk[bad, ], col="purple" , pch=19, cex=0.8)
  
      }


      lgbk$tag[bad_area] = "bad_area"
      lgbk$tag[not_in_martimes] = "not_in_martimes"
      lgbk$tag[too_shallow] = "too_shallow"
      lgbk$tag[too_deep] = "too_deep"
      lgbk$tag[no_position] = "no_position"

      to_re_assign_location = intersect( bad, which(!is.na(lgbk$region)))
      message( "Re-assigning position as lon/lat is incorrect")



      for (i in to_re_assign_location) {

        n = 0

        j = which(
          lgbk$tag == "ok" &
          lgbk$region == lgbk$region[i] & 
          lgbk$year == lgbk$year[i] & 
          lgbk$licence == lgbk$licence[i]  
        )

        if (length(j)>2) {
            o = setdiff(j, i)
            n = o[which(!is.na(lgbk$effort[o]))]
        }

        if (length(n) < 2) {
          j = which(
            lgbk$tag == "ok" &
            lgbk$region == lgbk$region[i] & 
            lgbk$year == lgbk$year[i] 
          )
          if (length(j)>2) {
            o = setdiff(j, i)
            n = o[which(!is.na(lgbk$effort[o]))]
          }
        }

        if (length(n) < 2) {
          j = which(
            lgbk$tag == "ok" &
            lgbk$region == lgbk$region[i] & 
            lgbk$licence == lgbk$licence[i] 
          )
          if (length(j)>2) {
            o = setdiff(j, i)
            n = o[which(!is.na(lgbk$effort[o]))]
          }
        }
 
        if (length(n) < 2) {
            j = which(
              lgbk$tag == "ok" &
              lgbk$region == lgbk$region[i]  
            )
            if (length(j)>2) {            
              o = setdiff(j, i)
              n = o[which(!is.na(lgbk$effort[o]))]
            }
        }

        if(length(n) < 2){
           message("Problem unexpected in re-assigning position: ", lgbk[i,])
           next()
        }  

        k = NULL
        m = which(!is.na( lgbk$effort[n] ) ) 
        if (length(m) > 1) {  
          if (length(m) == 1) {
            k = n[m]
          } else {
            w = lgbk$effort[n[m]] / sum(lgbk$effort[n[m]], na.rm=TRUE)
            k = sample( n[m], 1, prob=w )
          }
        } else {
          m = which(!is.na( lgbk$effort[n] ) ) 
          w = lgbk$effort[n[m]] / sum(lgbk$effort[n[m]], na.rm=TRUE)
          k = sample( n[m], 1, prob=w )
        }

        if (is.null(k)) {
          j = which(
            lgbk$tag == "ok" &
            lgbk$region == lgbk$region[i]  
          )
          k = sample( j, 1 )
        }
        
        lgbk$lon[i]  = lgbk$lon[k] 
        lgbk$lat[i]  = lgbk$lat[k] 
        lgbk$plon[i] = lgbk$plon[k] 
        lgbk$plat[i] = lgbk$plat[k] 
      
      }
 
      read_write_fast( lgbk, fn=filename  )

      return(filename)

    }


  # -------


    if (DS %in% c("logbook", "logbook.redo")) {

      filename = file.path( project.datadirectory("bio.snowcrab"), "data", "logbook", "logbook.rdz" )

      if (DS=="logbook") {
        logbook = read_write_fast( filename )
        return(logbook)
      }
 
      # logbooks from marfissci tables (not _1996_) actually 2002 to present
      logbook = logbook.db( DS="rawdata.logbook", yrs=1996:(p$year.assessment + 1) )  # 2002 to present add one to catch 4X effort/landings in next year
      setDT(logbook)
      names( logbook ) = tolower( names( logbook ) )
      names( logbook ) = rename.bio.snowcrab.variables(names( logbook ))

      logbook$marfis_cfa = toupper( as.character( logbook$cfa ) )  # this designation from marfis
  
      logbook$soak.time = logbook$soak_days * 24  # make into hours
      logbook$trap.type = logbook$trap_type 
      logbook$status = ""  # "" "" ""

      # range(logbook$date.landed)
      ii = which(!is.finite( logbook$date.fished ) )
      if (length(ii) > 0) logbook$date.fished[ii] = logbook$date.landed[ii]  # to get approximate date

      logbook$landings = logbook$pro_rated_slip_wt_lbs * 0.45359  # convert to kg
      logbook$cpue = logbook$landings / logbook$effort
      logbook$depth = logbook$depth_fm*1.83

      logbook$lat =   round( as.numeric(substring(logbook$lat, 1,2)) + as.numeric(substring(logbook$lat, 3,6))/6000 ,6)
      logbook$lon = - round((as.numeric(substring(logbook$lon, 1,2)) + as.numeric(substring(logbook$lon, 3,6))/6000), 6)

      to.extract = c( "year","lat","lon","depth","landings","effort","soak.time",
                      "cpue","trap.type","cfv","status","licence", "vessel", "pot_size", 
                      "date.landed", "date.fished", "cfa" )

      logbook = logbook[, ..to.extract ]
      # range( logbook$date_fished, na.rm=TRUE )  # 2002 to present

      bring_in_obsolete_data = FALSE
      if ( bring_in_obsolete_data) {
        # though we will be dropping these older historical records (pre-2002), this show how to integrate it with the marfis-based data series
        logbooks_pre_marfis = logbook.db( DS="fisheries.historical" )  # very spotty .... 
        range( logbooks_pre_marfis$year) # 2002 to 2003 
        logbook = rbind( logbooks_pre_marfis[,..to.extract], logbook[,..to.extract] )  # no need to bind as they are overlapping and marfis is cleaner
      }

      # known errors:  manual fixes
      ix = which( round(logbook$lat)==46 & round(logbook$lon)==-5930 )

      if (length(ix) > 0)  logbook$lon[ ix ]  = logbook$lon[ ix ] / 100

      # determine fishing are using licence info
      lic = logbook.db( DS="rawdata.licence" )
      names(lic) = tolower( names(lic) )
      setDT(lic)
      setnames(lic, "area_id", "marfis_area_id")  # directly from MARFIS
 
      marfis_areas = logbook.db( DS="rawdata.areas" )
      names(marfis_areas) = tolower( names(marfis_areas) )
      setDT(marfis_areas)
      setnames(marfis_areas, "area_id", "marfis_area_id")  # directly from MARFIS
      setnames(marfis_areas, "area", "marfis_area")  # directly from MARFIS
      setnames(marfis_areas, "desc_eng", "marfis_area_desc")  # directly from MARFIS

      marfis_areas = marfis_areas[ , c("marfis_area_id", "marfis_area", "area_type_id", "marfis_area_desc" ) ]  # reduce size
      marfis_areas$marfis_area_desc = as.character(marfis_areas$marfis_area_desc ) 

      snowcrab.marfis_area = c(27, 304, 305, 306, 307, 308, 615, 616, 619, 620, 790, 791, 792, 793, 794, 795, 796, 921, 922, 923, 1058, 1059  )
      
      marfis_areas = marfis_areas[  marfis_area_id %in% snowcrab.marfis_area , ] 
      lic = lic[ marfis_area_id %in% snowcrab.marfis_area , ]  # reduce size to crab marfis_areas only
      
      lic = merge( lic, marfis_areas, by="marfis_area_id", all.x=T, all.y=F )
      lic$marfis_area = toupper( as.character( lic$marfis_area) )

      # data dump from marfissci.areas table (2011)
      # from 
      # marfis_areas = marfis_areas[ grep ("crab", marfis_areas$marfis_area_desc, ignore.case=T ) ,   ]  # reduce size
      # AREA_ID   AREA AREA_TYPE_ID                      DESC_ENG                                DESC_FRE PARENT_AREA_ID ACT_FLAG      CUSER
      #    304     20            9        CRAB FISHING AREA - 20  ZONE DE P?CHE DU CRABE DES NEIGES - 20             NA        Y CONVERSION
      #    305     21            9        CRAB FISHING AREA - 21  ZONE DE P?CHE DU CRABE DES NEIGES - 21             NA        Y CONVERSION
      #    306     22            9        CRAB FISHING AREA - 22  ZONE DE P?CHE DU CRABE DES NEIGES - 22             NA        Y CONVERSION
      #    615    22I            9 CRAB FISHING AREA - 22I INNER ZONE DE P?CHE DU CRABE - 22I INT?RIEUR>            306        Y     MARFIS
      #    619    22O            9 CRAB FISHING AREA - 22O OUTER ZONE DE P?CHE DU CRABE - 22O INT?RIEUR>            306        Y     MARFIS

      #    307     23            9        CRAB FISHING AREA - 23  ZONE DE P?CHE DU CRABE DES NEIGES - 23             NA        Y CONVERSION
      #    790    23B            9         CRAB FISHING AREA 23B ZONE DE P?CHE DU CRABE DES NEIGES - 23B             NA        Y   SCHLEITC
      #     791    23C            9         CRAB FISHING AREA 23C ZONE DE P?CHE DU CRABE DES NEIGES - 23C             NA        Y   SCHLEITC
      #    1058    23S            9       CRAB FISHING AREA - 23S            ZONE DE P?CHE DU CRABE - 23S             NA        Y     MARFIS
      #     792    23D            9         CRAB FISHING AREA 23D ZONE DE P?CHE DU CRABE DES NEIGES - 23D             NA        Y   SCHLEITC

      #     308     24            9        CRAB FISHING AREA - 24  ZONE DE P?CHE DU CRABE DES NEIGES - 24             NA        Y CONVERSION
      #     616    24E            9  CRAB FISHING AREA - 24E EAST        ZONE DE P?CHE DU CRABE - 24E EST            308        Y     MARFIS
      #     793    24B            9         CRAB FISHING AREA 24B ZONE DE P?CHE DU CRABE DES NEIGES - 24B             NA        Y   SCHLEITC
      #     794    24C            9       CRAB FISHING AREA - 24C ZONE DE P?CHE DU CRABE DES NEIGES - 24C             NA        Y   SCHLEITC
      #     795    24D            9         CRAB FISHING AREA 24D ZONE DE P?CHE DU CRABE DES NEIGES - 24D             NA        Y   SCHLEITC
      #    1059    24S            9       CRAB FISHING AREA - 24S            ZONE DE P?CHE DU CRABE - 24S             NA        Y     MARFIS

      #      27     4X            1            NAFO DIVISION - 4X             DIVISION DE L?OPANO - 4X             NA
      #     620    24W            9  CRAB FISHING AREA - 24W WEST      ZONE DE P?CHE DU CRABE - 24W OUEST            308        Y     MARFIS
      #     796    24H            9         CRAB FISHING AREA 24H              SNOW CRAB FISHING AREA 24H             NA        Y   SCHLEITC


      #     922  CFA24            0 OBSERVER AREA-SNOW CRAB CFA24           OBSERVER AREA-SNOW CRAB CFA24             NA        Y     MARFIS
      #     921  CFA23            0 OBSERVER AREA-SNOW CRAB cfa23           OBSERVER AREA-SNOW CRAB CFA23             NA        Y     MARFIS
      #     923 SURVEY            0   OBSERVER AREA-SNOW CRAB SUR             OBSERVER AREA-SNOW CRAB SUR             NA        Y     MARFIS
        

      # licence information
      
      # determine 4x lic id's
      
      lic_CFA4X = unique( lic$licence_id[ which( lic$marfis_area %in% c( "4X", "24W" ) ) ])  # known to be in CFA 4X
      lic_CFA24 = unique( lic$licence_id[ which( lic$marfis_area %in% c( "24A", "24B", "24C", "24D", "24E", "24S","24H","CFA24" ) ) ]) # known to be in CFA 24

      north = which( lic$marfis_area %in% c( "20", "21" ,"22", "22I", "22O" )   )  # use upper case
      cfa23 = which( lic$marfis_area %in% c( "23", "23A", "23B", "23C", "23D", "23S", "CFA23") )
      
      # "marfis_area=24" contains both CFA24 and CFA4X .. by default consider all to be part of CFA24
      #  and then recode after the fact those that belong to 4X---Feb 2015 add in the marfis_area 308 ie '24' to cfa24 only this was the big difference between ben's and my data
      #cfa24 = which( lic$marfis_area %in% c( "24A", "24B", "24C", "24D", "24E", "24S" ,'24') | ( (lic$marfis_area=="24") & ( lic$licence_id %in% lic_CFA24 ) ) ) 
      #cfa4x = which( lic$marfis_area %in% c( "4X", "24W" ) | ( (lic$marfis_area=="24") & ( lic$licence_id %in% lic_CFA4X ) 
      cfa24 = which( lic$marfis_area %in% c( "24A", "24B", "24C", "24D", "24E", "24S" ,'24') ) 
      cfa4x = which( lic$marfis_area %in% c( "4X", "24W" )) 

      lic$subarea = NA
      lic$subarea [north] = "cfanorth"
      lic$subarea [cfa23] = "cfa23"
      lic$subarea [cfa24] = "cfa24"
      lic$subarea [cfa4x] = "cfa4x"

      lic$region = NA
      lic$region [north] = "cfanorth"
      lic$region [cfa23] = "cfasouth"
      lic$region [cfa24] = "cfasouth"
      lic$region [cfa4x] = "cfa4x"

      lic = lic[ which( !is.na(lic$subarea)) , ]
      lic = lic[ - which(duplicated( lic$licence_id) ) ,]  # finalised licencing information
      lic$licence = as.character(lic$licence_id) 
      lic = lic[, c("licence", "marfis_area_id", "subarea", "region", "marfis_area", "marfis_area_desc") ]
   
      logbook = merge(logbook, lic, by="licence", all.x=T, all.y=F, sort=F)
 
      logbook = lonlat2planar( logbook,  proj.type=p$aegis_proj4string_planar_km )

      logbook$cfa_historical = gsub("[ABCDEF]", "", logbook$marfis_area)
      logbook$cfa_historical = tolower(paste("cfa", logbook$cfa_historical, sep="") )
      i = grep( "NA", logbook$cfa_historical )
      logbook$cfa_historical[i] = NA

      i = which(logbook$cfa_historical %in% "cfaslope")
      logbook$cfa_historical[i] = logbook$region[i]

      # prefer licence-based area designation as position (lon/lat) are often in error
      for (v in c( "subarea", "region", "cfa_historical" ) ) {

        i.missing = which( is.na( logbook[[v]] ) )
        if (length( i.missing) > 1) {
        # try to determine via geographics:
          G = rep( NA, length(i.missing) )
          FF = logbook[ i.missing , c("lon", "lat")]
          ids = unique( na.omit( logbook[[v]] ) )
          for (i in ids) {
            j = polygon_inside(FF, i)      
            if (length(j) > 0) {
              G[j] = i
              message( "Georeferencing to management unit: ", i,  length(j))
            }
          }
          logbook[[v]][ i.missing ] = G
        }
      
      }

      logbook$cfa_rawdata = logbook$cfa
      logbook$cfa = logbook$region # duplicate in case old functions require it ("cfa" is deprecated)

      # any remaining without proper area designation are likely errors 
      # they are kept in case other factors arise in future ...

      # cfa 4X has a fishing season that spans two years recode "yr" to accomodate this      i.cfa4x = which( logbook$region == "cfa4x" )
      to.offset = which( 
        logbook$region =="cfa4x" &
        lubridate::month(logbook$date.landed) >= 1 & 
        lubridate::month(logbook$date.landed) <= 6 
      )

      logbook$yr = logbook$year
      logbook$yr[to.offset] = logbook$yr[to.offset] - 1
      message( "Fishing 'yr' for CFA 4X has been set to starting year:: 2001-2002 -> 2001, etc.")
 
      # enforce bounds in effort and cpue
      oo = which( logbook$cpue > 650 * 0.454 )  # 600 - 650 lbs / trap is a real/reasonable upper limit
      if ( length(oo) > 0 ) logbook$cpue[oo] = NA

      pp = which ( logbook$effort > 240 ) # small traps were used at times with large 240 trap compliments
      if ( length(pp) > 0 ) logbook$effort[pp] = NA

      # merge historical data ... these are aggregated to subareas and so not completely useful 
      # except where aggregate values are desired
      fdb_hist = logbook.db( DS="aggregated_historical" )  # 1978 to 2005
      fdb_hist = fdb_hist[ year <= 2001, ]
      vns = names(logbook)
      logbook = rbind( fdb_hist[ , ..vns], logbook )
      logbook = logbook[ !is.na(region), ]  # missing are NFLD and Gulf, based on map positions
    
      read_write_fast(logbook, fn=filename  )  # this is for plotting maps, etc

      return( "Complete" )

    }

  # ---------------------

    if (DS %in% c( "fishing.grounds.annual", "fishing.grounds.global", "fishing.grounds.redo")) {

      fn1 = file.path(  project.datadirectory("bio.snowcrab"), "data", "logbook", "fishing.grounds.global.rdz")
      fn2 = file.path(  project.datadirectory("bio.snowcrab"), "data", "logbook", "fishing.grounds.annual.rdz")

      if (DS=="fishing.grounds.global" | DS=="fishing.grounds" ) {
        fg = read_write_fast( fn1 )
        return (fg)
      }

      if (DS=="fishing.grounds.annual") {
        fg = read_write_fast( fn2 )
        return (fg)
      }

      out = NULL
      fg.res = p$fisheries.grid.resolution

      x = logbook.db( DS="logbook.filtered.positions" )
      setDF(x)
      grid = spatial_grid(p=p, DS="planar.coords")

      x$plon = grid_internal( x$plon, grid$plon )
      x$plat = grid_internal( x$plat, grid$plat )
      yrs = c(T,F)
      message ("The following warnings are ok: JC ... just a few NA's created which will be removed..")
      for ( Y in yrs ) {
        if (Y) {
          fn = fn2
          x$gridid = paste(x$plon%/%fg.res*fg.res, x$plat%/%fg.res*fg.res, x$year, sep="." )
          ncols=3
        } else {
          fn = fn1
          x$gridid = paste(x$plon%/%fg.res*fg.res, x$plat%/%fg.res*fg.res, sep="." )
          ncols=2
        }
        v = "visits"
        x$visits=1
        w = x[is.finite(x[,v]),]
        tmp = as.data.frame( xtabs( as.integer(x[,v]) ~ as.factor(x[,"gridid"]), exclude="" ) )
        names(tmp) = c("gridid", "total.visits")
        out = tmp

        v = "landings"
        w = x[is.finite(x[,v]),]
        tmp = as.data.frame( xtabs( as.integer(x[,v]) ~ as.factor(x[,"gridid"]), exclude="" ) )
        names(tmp) = c("gridid", "total.landings")
        out = merge( out, tmp, by="gridid", all=T, sort=F)

        v = "effort"
        w = x[is.finite(x[,v]),]
        tmp = as.data.frame( xtabs( as.integer(x[,v]) ~ as.factor(x[,"gridid"]), exclude="" ) )
        names(tmp) = c("gridid", "total.effort")
        out = merge( out, tmp, by="gridid", all=T, sort=F)

        out$total.cpue = out$total.landings / out$total.effort

        tmp = matrix(unlist(strsplit(as.character(out$gridid), ".", fixed=T)), ncol=ncols, byrow=T)
        out$plat = as.numeric(tmp[,2])
        out$plon = as.numeric(tmp[,1])
        out$gridid = as.character( out$gridid )

        if (Y) out$yr = as.numeric(tmp[,3])

        fg = out[ which(is.finite(out$plat+out$plon)), ]
        read_write_fast( fg, fn=fn )
      }

      return( "Complete")
    }


  # -----------------------------


    if (DS %in% c("fisheries.historical", "fisheries.historical.redo" )) {
       
      fn = file.path( project.datadirectory("bio.snowcrab"), "data", "logbook", "logbook.historical.rdz" )

      if (DS=="fisheries.historical") {
        logbooks_pre_marfis = NULL
        if (file.exists(fn)) logbooks_pre_marfis = read_write_fast( fn)
        return( logbooks_pre_marfis )
      }

      historicaldataloc = file.path( project.datadirectory("bio.snowcrab"), "data", "logbook", "archive")
      files = c("logbooks1996.csv", "logbooks1997.csv", "logbooks1998.csv", "logbooks1999.csv", "logbooks2000.csv", "logbooks2001.csv")
      out = NULL
      for (g in files) {
        f = file.path(historicaldataloc, g )
        a = read.table(file=f, skip=1, as.is=T, strip.white=T, sep=";")
        a = a[,c(1:14)]
        names(a) = c("cfv", "areas", "status", "licence", "date.landed", "date.fished", "lat", "lon",
            "landings.kg", "landings.lbs", "effort", "soak.time", "cpue", "trap.type" )
        a$file = f
        a$yr = as.numeric(gsub("[[:alpha:]]", "", f))
        out = rbind(out, a)
      }

      out$date.landed = lubridate::mdy( out$date.landed, tz="America/Halifax" )
      out$date.fished = lubridate::mdy( out$date.fished, tz="America/Halifax" )
      ii = which(!is.finite( out$date.fished ) )
      if (length(ii) > 0) out$date.fished[ii] = out$date.landed[ii]
       
      # there are many issues with this data set ... 
      f4x = read.table(file=file.path(historicaldataloc, "logbooks4x.csv"), skip=1, as.is=T, strip.white=T, sep=";")
      names(f4x) = c("cfv", "areas",  "date.landed", "date.fished", "lat", "lon",
            "landings.kg", "landings.lbs", "effort", "soak.time", "cpue", "trap.type" )
        
      f4x$date.landed = lubridate::ymd( f4x$date.landed, tz="America/Halifax" )
      f4x$date.fished = lubridate::ymd( f4x$date.fished, tz="America/Halifax" )

      f4x$areas="4X"
      f4x$status = NA
      f4x$licence = NA
      f4x$file = "logbooks4x.csv"
      f4x$yr = lubridate::year(f4x$date.landed )
      f4x = f4x[, names(out) ]

      logbooks_pre_marfis = rbind (out, f4x)
      logbooks_pre_marfis$lon = -logbooks_pre_marfis$lon

      logbooks_pre_marfis$depth = NA
      logbooks_pre_marfis$year = logbooks_pre_marfis$yr
      logbooks_pre_marfis$landings = logbooks_pre_marfis$landings.kg
      logbooks_pre_marfis$cpue = logbooks_pre_marfis$landings / logbooks_pre_marfis$effort
  
      to.extract = c( "year","lat","lon","depth","landings","effort","soak.time",
                      "cpue","trap.type","cfv","status","licence",
                      "date.landed", "date.fished")
      
      logbooks_pre_marfis = logbooks_pre_marfis[, to.extract]
      logbooks_pre_marfis = logbooks_pre_marfis[ -which( !is.finite(logbooks_pre_marfis$year) ), ]
 
      logbooks_pre_marfis$pot_size = ""
      logbooks_pre_marfis$vessel = ""
      logbooks_pre_marfis$cfa = ""
  
      to.extract = c( "year","lat","lon","depth","landings","effort","soak.time",
                      "cpue","trap.type","cfv","status","licence", "vessel", "pot_size", 
                      "date.landed", "date.fished", "cfa" )
      
      setDT( logbooks_pre_marfis )

      logbooks_pre_marfis = logbooks_pre_marfis[, ..to.extract]
  
      read_write_fast(logbooks_pre_marfis, fn=fn )

      return( "completed")

    }   # end if historical



    # -------------------------


    if (DS %in% c("fisheries.complete", "fisheries.complete.redo" )) {

      fn = file.path( project.datadirectory("bio.snowcrab"), "data", "logbook", "logbook.complete.rdz" )

      if (DS=="fisheries.complete") {
        logbook = read_write_fast( fn)
        return( logbook )
      }

			logbook = logbook.db(DS="logbook.filtered.positions")

      nl0 = nrow( logbook )

      grid = spatial_grid(p=p, DS="planar.coords")

      logbook$plon = grid_internal( logbook$plon, grid$plon )
      logbook$plat = grid_internal( logbook$plat, grid$plat )

			logbook$timestamp = as.POSIXct( logbook$date.landed, tz="America/Halifax", origin=lubridate::origin  )
      logbook$timestamp = with_tz( logbook$timestamp, "UTC")

      logbook$dyear = lubridate::decimal_date( logbook$timestamp ) - lubridate::year(logbook$timestamp )

      logbook = logbook[which(logbook$yr<=p$year.assessment),]
			# bring in time invariant features:: depth
      logbook$z = logbook$depth
			logbook$depth = NULL
      oo =  which( logbook$z < 10 | logbook$z > 500 ) # screen out large z's
      if (length(oo) > 0 )  logbook$z[ oo ] = NA

      ii = which(!is.finite(logbook$z))
      if (length(ii)>0) {
        
        pL = aegis.bathymetry::bathymetry_parameters( project_class="stmv"  )
        LUT= aegis_survey_lookuptable( aegis_project="bathymetry", 
          project_class="core", DS="aggregated_data", pL=pL )
        logbook$z[ii] = aegis_lookup( 
          pL=pL, LOCS=logbook[ ii, c("lon", "lat")],   LUT=LUT,
          project_class="core", output_format="points", 
          space_resolution=p$pres*2, variable_name="z.mean"  ) # core=="rawdata"
      }
      logbook$z = log( logbook$z )

      ii = which( ! is.finite( logbook$z) )
      if (length(ii)>0) logbook = logbook[ -ii, ]

      # bring in time varing features:: temperature
      ii = which(!is.finite(logbook$t))
      if (length(ii)>0){
        
        pL = aegis.temperature::temperature_parameters( project_class="core", yrs=p$yrs )
        LUT= aegis_survey_lookuptable( aegis_project="temperature", 
          project_class="core", DS="aggregated_data", pL=pL )

        logbook$t[ii] = aegis_lookup( pL=pL, 
          LOCS=logbook[ ii, c("lon", "lat", "timestamp")],  LUT=LUT,
          project_class="core", output_format="points", 
          space_resolution=p$pres*2, time_resolution=p$tres*2, 
          year.assessment=p$year.assessment, 
          variable_name="t.mean", tz="wAmerica/Halifax"  )
      }

			read_write_fast( logbook, fn=fn )

      return  ("Complete")
    }

    # -------------------------

    if (DS %in% c("logbook.gridded",  "logbook.gridded.redo" ) ) {

      loc = file.path( project.datadirectory("bio.snowcrab"), "data", "logbook", "gridded.fishery.data" )
      dir.create( path=loc, recursive=T, showWarnings=F)

      if (DS == "logbook.gridded") {
        fn = file.path(loc, paste( "gridded.fishery", yrs, "rdz", sep=".") )
        gridded.fishery.data = read_write_fast(fn)
        return(gridded.fishery.data)
      }

      yy = logbook.db(DS="logbook")
      setDT(yy)
      yrs = sort( unique( yy$year))

      for ( y in yrs ) {

        fn = file.path(loc, paste( "gridded.fishery", y, "rdz", sep=".") )
        print (fn)

        # load logbook info: global summary
        fg = logbook.db( DS="fishing.grounds.global" )  # in dataframe fg
        fg0 = regrid.lonlat(old=fg, res=p$fisheries.grid.resolution, vr.to.sum=c("total.landings", "total.effort", "total.visits"))
        fg0$core.visits0 = ifelse( fg0$total.visits >= 5, 1, 0 )
        fg0$total.cpue = log( fg0$total.landings / fg0$total.effort + 1)
        fg0$total.landings = log( fg0$total.landings+1 )
        fg0$total.effort = log( fg0$total.effort+1 )
        fg0$total.visits = log( fg0$total.visits+1 )
        fg0 = fg0[ , c( "gridid", "core.visits0", "total.effort", "total.cpue", "total.landings", "total.visits" ) ]
        rm(fg)

        # load logbook info: annual
        fg = logbook.db( DS="fishing.grounds.annual"  )  # in dataframe fg
        ii = which(fg$yr==y)
        if (length(ii) <  5) next()
        fg = fg[ ii,]

        fg = regrid.lonlat(old=fg, res=p$fisheries.grid.resolution, vr.to.sum=c("total.landings", "total.effort", "total.visits") )
        fg$core.visits = ifelse(fg$total.visits >= 3, 1, 0)
        fg$core.visits = ifelse(fg$total.visits >= 3, 1, 0)
        fg$core.landings = ifelse(fg$total.landings >= 5, 1, 0)
        fg$core.effort = ifelse(fg$total.effort >= 100, 1, 0)
        fg$cpue = log( fg$total.landings / fg$total.effort + 1)
        fg$landings = log( fg$total.landings+1 )
        fg$effort = log( fg$total.effort+1 )
        fg$visits = log( fg$total.visits +1 )

        fg = fg[ , c( "gridid", "core.visits", "core.landings", "core.effort", "cpue", "landings", "effort", "visits" ) ]
        fg = fg[ is.finite(fg$cpue * fg$landings * fg$effort),]
        fg = fg[ which(fg$core.visits==1) , ]

        gridded.fishery.data = merge(fg0, fg, by="gridid", all.x=T, all.y=T, sort=F)
        read_write_fast( gridded.fishery.data, fn=fn )
      }
    } # end gridded fishieries data


  if ( DS=="aggregated_historical") {

    fn = file.path( project.datadirectory("bio.snowcrab"), "data", "fisheries", "aggregated_historical_data.rdz" )

    if (!redo) {
      out = NULL  
      out = read_write_fast(fn)
      return(out)
    }

    a = fread(file.path(  project.datadirectory("bio.snowcrab"), "data", "fisheries", "annual.landings.csv") )
    names(a) = c("yr","cfa_historical","nlicences","tac.tons","landings.tons","cpue.kg.trap","effort.100traps")
    numerics = c("yr","nlicences","tac.tons","landings.tons","cpue.kg.trap","effort.100traps")
    
    # create vars to match logbook format:
    a$year = a$yr
    a$date.fished = a$date.landed = lubridate::ymd( paste( a$yr, "07", "01", sep="-"), tz="America/Halifax" )  # force july 1
    a$lat = NA
    a$lon = NA
    a$plon = NA
    a$plat = NA
    a$depth = NA
    a$landings = a$landings.tons * 1000
    a$landings.tons = NULL
    a$effort = a$effort.100traps * 100
    a$effort.100traps = NULL
    a$soak.time = NA
    a$cpue = a$cpue.kg.trap
    a$cpue.kg.trap = NULL

    a$trap.type = NA
    a$cfv = NA
    a$status = NA
    a$licence = NA
    a$region = NA
    a$area_id = NA
    a$vessel = NA
    a$pot_size = NA
    a$marfis_area_id = NA
    a$marfis_area =NA
    a$marfis_area_desc =NA
 
    a$region = NA
    a$region[which(a$cfa_historical %in% c("cfa20", "cfa21", "cfa22", "cfanorth") )] = "cfanorth"
    a$region[which(a$cfa_historical %in% c("cfa23", "cfa24", "cfasouth", "cfaslope") )] = "cfasouth"
    a$region[which(a$cfa_historical %in% c("cfa4x") )] = "cfa4x"
    #cfaslope spans 23 and 24 ... cannot separate out unless we go back to marfis ... not sure if available in 2001 
    a$cfa = a$region # copy
    a$cfa_rawdata = a$cfa_historical # copy

    a$subarea = NA
    a$subarea[which(a$cfa_historical %in% c("cfa20", "cfa21", "cfa22", "cfanorth") )] = "cfanorth"
    a$subarea[which(a$cfa_historical %in% c("cfa23" ) )] = "cfa23"
    a$subarea[which(a$cfa_historical %in% c("cfa24" ) )] = "cfa24"
    a$subarea[which(a$cfa_historical %in% c("cfa4x") )] = "cfa4x"

    a$no.lic[ which(a$subarea=="cfaslope" & a$yr==2001) ] = 4
    a$no.lic[ which(a$subarea=="cfaslope" & a$yr==2002) ] = 4
    a$no.lic[ which(a$subarea=="cfaslope" & a$yr==2003) ] = 5

    read_write_fast( a, fn=fn )
    return(a)
  }



  if ( DS=="carstm_inputs") {

    # this cleans up catch and effort data for modelling, by dropping problematic records 
    # as such, the total catch and effort will not be consistent with fishery quotamanagement
    
  # maturity codes
  immature = 0
  mature = 1
  mat.unknown = 2

    # prediction surface
    crs_lonlat = st_crs(projection_proj4string("lonlat_wgs84"))
    if (is.null(sppoly)) sppoly = areal_units( p=p )  # will redo if not found
    sppoly = st_transform(sppoly, crs=crs_lonlat )
    # sppoly$data_offset = sppoly$sa
    
    areal_units_fn = attributes(sppoly)[["areal_units_fn"]]

    if (exists("carstm_directory", p)) {
      outputdir =  p$carstm_directory 
    } else {
      outputdir = file.path( p$modeldir, p$carstm_model_label )
    }
    outfn = paste( sep="_") # redundancy in case other files in same directory
    
    fn = file.path( outputdir, "carstm_inputs.rdz" )
 
    if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

    if (!redo)  {
      if (file.exists(fn)) {
        message( "Loading previously saved carstm_inputs ... ", fn)
        M = aegis::read_write_fast( fn )
        return( M )
      }
    }
    message( "Generating carstm_inputs ... ", fn)

    M = logbook.db( DS="logbook" )

    h = which( !is.na( M$lon + M$lat ) )
    M = M[h,]

    i = polygon_inside( M[, c("lon", "lat")], region="cfaall")
    M = M[i,]

    j = polygon_inside( M[, c("lon", "lat")], region="isobath1000m")
    M = M[j,]

    # additional constraint ..
    # remove data that are strange in location .. land
    crs_lonlat = st_crs( projection_proj4string("lonlat_wgs84") )

    mp = st_multipoint(cbind(p$corners$lon, p$corners$lat))
    bbox =  st_as_sfc(st_bbox( mp) )
    st_crs(bbox) =  crs_lonlat 

    coast = st_transform( st_as_sf(coastline_db( p=p )), crs=crs_lonlat )
    coast = (
        st_intersection( coast, bbox )
        %>% st_buffer(0.01)
        %>% st_union()
        %>% st_cast("POLYGON" )
        %>% st_union()
        %>% st_make_valid()
    )

    inside = st_points_in_polygons(
        pts =st_as_sf( M[,c("lon", "lat")], coords=c("lon","lat"), crs=crs_lonlat ),
        polys = coast, 
        method="sp::point.in.polygon"
    )
    onland = which(inside)
    if (length(onland)>0) M = M[-onland, ]

    # filter by depth ..
    # use the match/map syntax in bathymetry and filter out shallow sets .. < 10 m? TODO
    # keep those in the domain and deeper than depth=10 m

    pL = aegis.bathymetry::bathymetry_parameters( project_class="stmv"  )
    LUT= aegis_survey_lookuptable( aegis_project="bathymetry", 
    project_class="core", DS="aggregated_data", pL=pL )

    M$z = M$depth
    i = which( !is.finite( M$z ))
    if (length(i) > 0 ) {
      M$z[i] = aegis_lookup( pL=pL, 
          LOCS=M[i, c("lon", "lat")], LUT=LUT,
          output_format="points" , 
          space_resolution=p$pres*2, variable_name="z.mean"  ) # core=="rawdata"
    }

    aoi = which( M$z > 10 ) # negative = above land
    M = M[ aoi,]
 
    # only accept "correctly" positioned data within each subarea ... in case miscoded data have a large effect
    icfa4x = polygon_inside( M[, c("lon", "lat")], "cfa4x")
    icfanorth = polygon_inside( M[, c("lon", "lat")], "cfanorth")
    icfa23 = polygon_inside( M[, c("lon", "lat")], "cfa23")
    icfa24 = polygon_inside( M[, c("lon", "lat")], "cfa24")

    gooddata = sort( unique( c(icfa4x, icfanorth, icfa23, icfa24 ) ) )
    M = M[gooddata, ]
    
    M$timestamp = M$date.fished
    M$tiyr = lubridate::decimal_date(M$timestamp)
    M$dyear = M$tiyr - lubridate::year(M$timestamp) 
         # reduce size
    M = M[ which( M$lon > p$corners$lon[1] & M$lon < p$corners$lon[2]  & M$lat > p$corners$lat[1] & M$lat < p$corners$lat[2] ), ]
    # levelplot(z.mean~plon+plat, data=M, aspect="iso")
 
    M = M[ which( is.finite(M$dyear) ), ]

    # enforce bounds in effort and cpue
    oo = which( M$cpue > 650 * 0.454 )  # 600 - 650 lbs / trap is a real/reasonable upper limit
    if ( length(oo) > 0 ) M$cpue[oo] = NA

    pp = which ( M$effort > 240 ) # small traps were used at times with large 240 trap compliments
    if ( length(pp) > 0 ) M$effort[pp] = NA
 
    qn = quantile( M$cpue, p$quantile_bounds[2], na.rm=TRUE )
    ni = which( M$cpue > qn )

    # copy landings to "catch" .. which we modify as below
    M$catch = M$landings
    M$catch[ni] = floor( qn * M$effort[ni] )  # limit upper bound to one that is physically possible
    M = M[is.finite(M$catch),]
    
    M$data_offset = M$effort
    M = M[is.finite(M$data_offset),]

    lower_threshold =  1 # 1 kg per trap seems a reasonable cutoff ~ 1 large crab
    M$pa = 0
    M$pa[ which( M$catch > lower_threshold) ] = 1
   
    # data_offset is per unit trap haul 
    M = carstm_prepare_inputdata( 
      p=p, M=M, sppoly=sppoly,
      APS_data_offset=1, 
      retain_positions_outside_of_boundary = 25,  # centroid-point unit of p$aegis_proj4string_planar_km
      vars_to_retain=c(
        "licence", "cfv", "cfa", "region", "z", "catch", "data_offset", "pa",  # "soak.time", "trap.type", are motly missing data
      ) 
    )
  
    if ( exists("substrate.grainsize", M)) M$log.substrate.grainsize = log( M$substrate.grainsize )

    if (!exists("yr", M)) M$yr = M$year  
    
    # IMPERATIVE: 
    i =  which(!is.finite(M$z))
    j =  which(!is.finite(M$t)) 

    if (length(j)>0 | length(i)>0) {
      warning( "Some areal units that have no information on key covariates ... you will need to drop these and do a sppoly/nb reduction with areal_units_neighbourhood_reset() :")
          print( "Missing depths:")
      print(unique(M$AUID[i]) )
      print( "Missing temperatures:")
      print(unique(M$AUID[j] ) )
    }
 
    # imperative covariate(s)
    M = M[ which(is.finite(M$z)), ]  
    M = M[ which(is.finite(M$t)), ]  
 
    M$space = match( M$AUID, sppoly$AUID) # for bym/car .. must be numeric index matching neighbourhood graphs
    M$space_time = M$space  # copy for space_time component (INLA does not like to re-use the same variable in a model formula) 
    M$space_cyclic = M$space  # copy for space_time component (INLA does not like to re-use the same variable in a model formula) 

    M$time = match( M$year, p$yrs ) # copy for space_time component .. for groups, must be numeric index
    M$time_space = M$time    
 
    M$cyclic = match( M$dyri, discretize_data( span=c( 0, 1, p$nw) ) ) 
    M$cyclic_space = M$cyclic # copy cyclic for space - cyclic component .. for groups, must be numeric index
  
    read_write_fast( data=M, fn=fn )

    return( M )
  }

}

