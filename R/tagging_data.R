
tagging_data=function(
  toget="",
  fn.loc = "C:/home/jae/bio.data/snowcrab/data/movement", 
  tble=""
){

  # historical mark recapture data
  if (toget=="otn") {
 
    mark = tagging_data( "current", tble="sct_acoustic_releases" )

    ## Note there are data errors: 
    ## tag_code_space is xxx-yyy-zzz 
    ## there are cases where  tag_id_code is xxx-yyy;  then zzz needs to be derived from  tag_id_code
    ## there are cases where the zzz in tag_code_space is wrong (length(zzz)==5) .. again replace with  tag_id_code (which is correct)
    i = nchar(mark$tag_code_space)  # as of 2026, this should be exactly 14
    bad1 = which(i < 14)  
    if (length(bad1) > 0) { 
      mark$tag_code_space[bad1] = paste(mark$tag_code_space[bad1], "-", mark$tag_id_code[bad1], sep="")
    }
    
    bad2 = which(i > 14)  
    if (length(bad2) > 0) { 
      mark$tag_code_space[bad2] = gsub("-[^-]+$", "", mark$tag_code_space[bad2] )
      mark$tag_code_space[bad2] = paste(mark$tag_code_space[bad2], "-", mark$tag_id_code[bad2], sep="")
    }
    

    setnames(mark, "animal_id", "aid" )
    setnames(mark, "release_latitude", "lat")
    setnames(mark, "release_longitude", "lon")
    setnames(mark, "capture_latitude", "lat_at_initial_capture")
    setnames(mark, "capture_longitude", "lon_at_initial_capture")
    setnames(mark, "capture_depth_m", "z_at_initial_capture")
    setnames(mark, "utc_release_date_time", "timestamp")

    setnames(mark, "length_m", "cw")
    mark$cw = mark$cw * 1000  # m -> mm

    setnames(mark, "length2_m", "chela")
    mark$chela = mark$chela * 1000  # m -> mm

    setnames(mark, "weight_kg", "wgt")
    mark$wgt = mark$wgt  
    
    setnames(mark, "life_stage", "mat")
    setnames(mark, "age", "cc")
    setnames(mark, "tag_code_space", "tagid")
    
    mark = mark[,.(aid, tagid, timestamp, lat, lon, lat_at_initial_capture, lon_at_initial_capture, z_at_initial_capture, cw, chela, wgt, sex, mat, cc )]
    mark$tag = 0 # indicate intial mark

    mark$timestamp = lubridate::ymd_hms(mark$timestamp)

    # recaptures
    recapture = tagging_data( "current", tble="sct_acoustic_detections" )
    setnames(recapture, "catalognumber", "pid")
    setnames(recapture, "tagname", "tagid")

    recapture$aid = gsub( "ZSC-", "", recapture$pid )
    recapture = recapture[, .(pid, aid, bottom_depth, tagid, datecollected, longitude, latitude )]
    names(recapture) = c("pid", "aid", "z", "tagid", "timestamp", "lon", "lat")
    recapture$timestamp =  lubridate::mdy_hms(recapture$timestamp)
    i = which(recapture$lon < -65 )  # strange data  .. shallow too .. probably not snow crab
    recapture =recapture[-i,]
 
    recapture$z = as.numeric(recapture$z)
 

    # discretize to a daily location
    daily_dt = recapture[, .(
        pid = unique(pid),
        aid = unique(aid),
        lon = mean( lon, na.rm = TRUE),
        lat = mean( lat, na.rm = TRUE),
        z = mean( z, na.rm = TRUE)
      ), 
      by = .( tagid, timestamp=as.Date(timestamp) ) 
    ]
      
    recapture = daily_dt  # overwrite

    setorder(recapture, tagid, timestamp)     # Order the data chronologically by tagid and then timestamp
    recapture[, tag := 1:.N, by = tagid] # Assign a sequential ID ('tag') for each event within a 'tagid' group
 
    # not used ,,, modelling path should come much later
    # sct_accoustic_path = tagging_data( "current", tble="sct_accoustic_path" )
    # sct_accoustic_paths = tagging_data( "current", tble="sct_accoustic_paths" )
    # sct_accoustic_paths # n=6608
    # setnames(sct_accoustic_paths, "long", "lon" )
    # sct_accoustic_paths$cid = as.character( sct_accoustic_paths$cid )
    # acoustic_detections = sct_accoustic_path[ sct_accoustic_paths, on=.(pid, cid)]  # n=6608
    # acoustic_detections = sct_acoustic_detections[ acoustic_detections, on=.(pid)]

    # otn = list( mark=mark, recapture=recapture )

    # collapse into one table
    mark$timestamp = as.Date(mark$timestamp)

    otn = rbind( mark, recapture, fill=TRUE )
    otn = filter_dead(otn) # flags as is_dead (not moving)
    otn = otn[!is.na(tagid),]
    otn$z = as.numeric(otn$z)

    otn$datasource = "otn"
    otn[ nchar(tagid) > 10 ,]  # some tagid are too short ... investigate these eventually  .. dropping for now 
    
    return(otn)

  }

  if (toget=="bio") {

    tag_to_study_id = function(ids) {

      # determine study id from historical records
      out = rep(NA, length(ids))
      for (i in 1: length(ids)) {
        id = ids[i]
        study = NA
        if (id %in% c(0:600))      study=1
        if (id %in% c(2350:2399))  study=5
        if (id %in% c(1000:1600))  study=6
        if (id %in% c(2401:2403, 2411:2450)) study=7
        if (id %in% c(1601:2349))  study=8
        if (id %in% c(6000:6349))  study=9
        if (id %in% c(2716:2850))  study=10
        if (id %in% c(2456:2715))  study=11
        if (id %in% c(5050:5099))  study=12
        if (id %in% c(5100:5149))  study=13
        if (id %in% c(2851:2900, 3051:3262, 3446:3545))  study=14
        if (id %in% c(3263:3444, 3546:3600))  study=15
        if (id %in% c(4000:4246))  study=16
        if (id %in% c(7450:7542))  study=17
        if (id %in% c(7543:7699))  study=18
        if (id %in% c(4798:4999, 5250:5520, 6350:6999))  study=19

        if (id %in% c(5521:5808))  study=20
        if (id %in% c(4250:4480))  study=21
        if (id %in% c(4482:4797))  study=22
        if (id %in% c(7000:7449, 7700:7999))  study=23
        if (id %in% c(5809:5999, 8501:8570, 9229:9649))  study=24
        if (id %in% c(8571:8828))  study=25
        if (id %in% c(8829:9228))  study=26
        if (id %in% c(10482:10499, 11082:11334, 11350:11387))  study=27
        if (id %in% c(4301, 8040:8050, 8137:8150, 8201:8298, 9684:9741, 9748:9784, 10254:10406, 10468:10481, 11039:11081 ))  study=28
        if (id %in% c(4248:4249, 8000:8039, 8051:8068, 8070:8136, 8151:8200, 10407:10467, 11000:11036))  study=29
        if (id %in% c(8299:8500, 9650:9683))  study=30
        if (id %in% c(10032:10253))  study=31
        if (id %in% c(9742:9747, 9785:10031))  study=32
        if (id %in% c(10501:10530, 10609:10716))  study=33
        if (id %in% c(10531:10608, 10718:10798))  study=34
        if (id %in% c(10800:10949, 11400:11449))  study=35
        if (id %in% c(10950:10999, 11335:11349, 11388:11399, 11450:11499, 13000:13149))  study=36
        if (id %in% c(13150:13490))  study=37
        if (id %in% c(13491:13649, 15650:15673))  study=38
        if (id %in% c(13650:13834))  study=39
        if (id %in% c(13835:14109))  study=40
        if (id %in% c(14110:14226))  study=41
        if (id %in% c(14227:14306, 14308:14311, 14313:14314, 14316, 14318:14345, 14347:14430))  study=42
        if (id %in% c(14431:14880))  study=43
        if (id %in% c(15750:15779, 15800:15899, 15950:15999))  study=44
        if (id %in% c(15780:15799, 15900:15949))  study=45

        if (id %in% c(14881:14951))  study=46
        if (id %in% c(14952:15001))  study=47
        if (id %in% c(15002:15199))  study=48
        if (id %in% c(15400:15499, paste("t", 1605:1659, sep="")))  study=49
        if (id %in% c(15500:15649, 15674:15749, 16000:16711))  study=50
        if (id %in% c(16712:17199, 17300:17399, 17500:17607))  study=51

        if (id %in% c(15200:15399, paste("t", c(1217:1350, 1601:1603), sep="") ))  study=52
        if (id %in% c(17608:17849, 17953))  study=53
        if (id %in% c(17200:17299, 17850:17999))  study=54
        if (id %in% c(18043:18099))  study=55
        if (id %in% c(17401:17499, 18000:18042, 18100:18158,  paste("t", 1676:1694, sep="")))  study=56
        if (id %in% c(18226:18326, paste("t", 3026:3175, sep=""), paste("s", 99102:99196, sep="")))  study=57
        out[i] = study
      }
      return (out)

    }


    # marks
    sct_sample = tagging_data( "current", tble="sct_sample" )
    sct_sample_gulf = tagging_data( "current", tble="sct_sample_gulf" ) # empty
    sct_trip = tagging_data( "current", tble="sct_trip" )
    sct_trip_gulf = tagging_data( "current", tble="sct_trip_gulf" )  # empty

    setnames(sct_sample, "trip", "tripid")
    setnames(sct_trip, "trip_id", "tripid")

    mark = sct_trip[ sct_sample, on=.(tripid)]
    mark$z = as.numeric(mark$fathoms) * 1.8288
    setnames(mark, "release_date", "timestamp")
    setnames(mark, "lat_dd_dddd", "lat")
    setnames(mark, "long_dd_dddd", "lon")
    mark = mark[, .(tripid, sample_id, timestamp, lat, lon, z, comments)]

    mark$lon = as.numeric( mark$lon )
    mark$lat = as.numeric( mark$lat )
    mark[lon==0, lon:=NA]
    mark[lat==0, lat:=NA]
    mark[z==0, z:=NA]


    # add biologicals
    sct_bio = tagging_data( "current", tble="sct_bio" )
    names(sct_bio) = c("sample_id", "tagid", "cw", "chela", "cc", "durometer")
    mark = mark[sct_bio, on=.(sample_id)]

    mark[cc=="0", cc:=NA]


    # find historical marking data from tags:
    # use historical info to fill in the missing values:
    marked_historical = tagging_data( "marked_historical" )
    marked_historical$study_id = 1:nrow(marked_historical) # make study_id explicit
    recaptures_1996_2001 = tagging_data( "recaptures_1996_2001" )
    marked_1996_2001 = tagging_data( "marked_1996_2001" )

    i = which(is.na(mark$lat))
    studyids = tag_to_study_id( mark$tagid[i] )
    studyids = tag_to_study_id( gsub("G", "", mark$tagid[i]) )

    lookups = marked_historical[studyids,]
    mark$lat[i] =lookups$lat
    mark$lon[i] =lookups$lon
    mark$tripid[i] = paste("Historical", lookups$study_id, sep="_")
    mark$timestamp[i] = lubridate::ymd( paste(lookups$yr, "06", "01", sep="-") )

    # find historical marking data from tags:
    i = which(is.na(mark$timestamp))
    studyids = tag_to_study_id( gsub("G", "", mark$tagid[i]) )  # the G was added to the relational database causing issues
    lookups = marked_historical[studyids,]
    mark$lat[i] =lookups$lat
    mark$lon[i] =lookups$lon
    mark$tripid[i] = paste("Historical", lookups$study_id, sep="_")
    mark$timestamp[i] = lubridate::ymd( paste(lookups$yr, "06", "01", sep="-") )

    mark[ tripid=="Historical_NA", tripid:=NA ]  # these are likely tagid entry errors ... or Gulf tags
    mark$tag = 0 # intiial tagging event

    # recaptures

    # sct_bio_gulf = tagging_data( "current", tble="sct_bio_gulf" ) #nothing
    sct_capture = tagging_data( "current", tble="sct_capture" )
    sct_capture_gulf = tagging_data( "current", tble="sct_capture_gulf" )

    recapture = rbind( sct_capture, sct_capture_gulf )
    recapture$z = as.numeric(recapture$fathoms) * 1.8288
    setnames(recapture, "tag", "tagid")
    setnames(recapture, "lat_dd_dddd", "lat")
    setnames(recapture, "long_dd_dddd", "lon")
    setnames(recapture, "capture_date", "timestamp")
    setnames(recapture, "carapace_cond", "cc")
    recapture = recapture[, .(tagid, timestamp, lon, lat, z, cc)]

    setorder(recapture, tagid, timestamp)     # Order the data chronologically by tagid and then timestamp
    recapture[, tag := 1:.N, by = tagid] # Assign a sequential ID ('tag') for each event within a 'tagid' group


    recapture$lon = as.numeric( recapture$lon )
    recapture$lat = as.numeric( recapture$lat )
    recapture[lon==0, lon:=NA]
    recapture[lat==0, lat:=NA]
    recapture[z==0, z:=NA]
    recapture[cc=="0", cc:=NA]

    # bio = list( mark=mark, recapture=recapture )

    # collapse into one table
    bio = rbind(mark, recapture, fill=TRUE )  # failed to parse error -> no year found in initial data
    bio = bio[!is.na(tagid),]
    bio[z>350, z:=NA]  

    bio$chela = as.numeric(bio$chela)
    bio[chela == 0, chela:=NA]  
    bio[chela > 50, chela:=NA]  
    
    bio$cw = as.numeric(bio$cw)
    bio[cw == 0, cw:=NA]  
    bio[cw >161, cw:=NA]  

    bio$durometer = as.numeric(bio$durometer)
    bio[durometer == 0, durometer:=NA]  
    bio[durometer > 100, durometer:=NA]  

    too_long = c( "14552", "14557", "2074") # these were found in the analysis ...
    # bio_stats = summarize_tag_activity(bio)
    # likely data entry or submission errors
    # too_long = bio_stats[ duration_days > 2300, tagid] #  "14552" "14557" "2074"  
    
    bio = bio[ !(tagid %in% too_long), ] # drop them
 
    bio$timestamp = as.Date(bio$timestamp)
    bio$datasource = "bio" 

    return( bio )

  }


  if (toget=="data_dump_current") {
    # brent's database: Brent's library: https://github.com/brent0/SCtagging/blob/master/DESCRIPTION

    dir.create( fn.loc, recursive = TRUE, showWarnings = FALSE )

    require(DBI)
    require(ROracle)

    con = ROracle::dbConnect(
          DBI::dbDriver("Oracle"),
          dbname=oracle.snowcrab.server,
          username=oracle.snowcrab.user,
          password=oracle.snowcrab.password,
          believeNRows=F
        )

    # dump
    tbls = c(
        "SCT_ACCOUSTIC_PATH", "SCT_ACCOUSTIC_PATHS",
        "SCT_ACOUSTIC_DETECTIONS", "SCT_ACOUSTIC_RELEASES",
        "SCT_BIO", "SCT_BIO_GULF", "SCT_CAPTURE", "SCT_CAPTURE_GULF",
        "SCT_PATH", "SCT_PATHS", "SCT_PEOPLE", "SCT_PEOPLE_GULF",
        "SCT_SAMPLE", "SCT_SAMPLE_GULF", "SCT_TRIP", "SCT_TRIP_GULF",
        "SCT_ALL"
    )

    for (tb in tbls) {
     	sqlquery = paste( "select * from ", tb )
      o = NULL
     	o = ROracle::dbGetQuery(con, sqlquery)
      fn_out = file.path( fn.loc, paste(tolower(tb), ".rdz", sep="") )
      read_write_fast( data=o, fn=fn_out )
      gc()  # garbage collection
    }

    ROracle::dbDisconnect(con)
    return("data dump finished")
  }

  if (toget=="current") {

    if (tble =="") {
      # all spaghetti tags:
      res = read_write_fast( file.path( fn.loc, "sct_all.rdz" ) )
      names(res) = c(
          "tagid", "caplat", "caplon", "capdate", "caparea", "year",
          "relcode","csubarea", "rewarded", "sampdat", "area", "subarea",
          "sampyear","samplat", "samplon", "carapace", "chela", "cc"
      )

      library(data.table)
      res = as.data.table(res)

      # Convert numeric columns safely
      num_cols = c("caplat", "caplon", "samplat", "samplon", "carapace", "chela", "cc")
      res[, (num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols]

      # 1. Create Initial Tagging records (tag = 0)
      # We take the unique set of tagids and their initial sampling info
      marking_events = res[, .SD[1], by = tagid][, .(
        tagid = tagid,
        lon = samplon,
        lat = samplat,
        date = sampdat,
        tag = 0,
        carapace = carapace,
        chela = chela,
        cc = cc,
        relcode = relcode
      )]

      # 2. Create Recapture records
      # Each row in the original 'res' is a recapture event
      recapture_events = res[, .(
        tagid = tagid,
        lon = caplon,
        lat = caplat,
        date = capdate,
        carapace = carapace,
        chela = chela,
        cc = cc,
        relcode = relcode
      )]

      # Combine them
      tagdb = rbind(marking_events, recapture_events, fill = TRUE)

      # 3. Sort by tagid and Date, then re-calculate the tag index
      # This ensures that if there are multiple recaptures, they are indexed 1, 2, 3...
      setorder(tagdb, tagid, date)
      tagdb[, tag := 0:(.N - 1), by = tagid]

      # Final column ordering
      setcolorder(tagdb, c("tagid", "lon", "lat", "date", "tag", "carapace", "chela", "cc", "relcode"))

      # i = which( tagdb$lon == 0 | tagdb$lat==0)
      return(tagdb)
    } else {
      fn = file.path( fn.loc, paste(tolower(tble), ".rdz", sep="") )
      o = read_write_fast(fn)
      setDT(o)
      names(o) = tolower(names(o))
      return(o)
    }
  }

  if (toget=="marked_historical") {
    # marking event info as a series of ranges of tagid's
    marked2 =  read.table( file.path(fn.loc, "tags_summary1993_2005.csv" ), sep=";", header=T )
    return(marked2)
  }

  if (toget=="marked_1996_2001"){
    # marking events
    marked.file = "tags.1996_2001.csv"
    marked = fread( file.path( fn.loc, marked.file), sep=";", header=T )
    marked$Ncrabs = NULL
    f = which(marked$lon>-55)
    marked[f,] = NA
    marked$timestamp = lubridate::mdy( marked$date)
    v0 = c("tagID", "timestamp", "lon", "lat", "cw", "ch", "cc", "z.fm", "area" )
    marked = marked[, ..v0]
    names(marked) = tolower( names(marked ))
    return(marked)
  }

  if (toget=="recaptures_1996_2001") {

    # recaptures
    recaps.file="recaptures.csv"
    recaps = fread( file.path(fn.loc, recaps.file), sep=";", header=T )
    f = which(recaps$lon>-50 & recaps$lon < 100 )
    recaps[f,"lon"] = -recaps[f,"lon"]
    recaps = recaps[!is.na(recaps$yr) ,]

    f = which(recaps$month>12)
    months = recaps[f, "month"]
    recaps[f, "month"] = recaps[f, "date"]
    recaps[f, "date"] = months

    f = which(is.na(recaps$month))
    recaps[f, "month"] = 12 # assume a December capture

    f = which(recaps$date>31)
    recaps[f, "date"] = 1  # assume first day of the month

    f = which(is.na(recaps$date))
    recaps[f, "date"] = 1  # assume first day of the month

    recaps$timestamp = lubridate::ymd( dates.=paste(recaps$yr, recaps$month, recaps$date, sep="-"))
    recaps$julian = lubridate::yday(recaps$timestamp)

    cc3 = which(recaps$cc %in% c(
      "Intermediate", "Legal/hard(released)", "clean",
      "intermediate", "propre", " Legal (released)         ",
      "Clean / propre", "Intermediate/interm\351diaire",
      "Interm\351diate", "good condition",
      "very good condition", "Hard shell", "Legal (released)",
      "inter"))
    ccM = which(recaps$cc %in% c( "dirty", "Dirty / sale", "3M", "3m", "Mossy/Dirty", "mousseux",
                                  "intermediate (3M)", "sale", "mossy"))
    ccD = which(recaps$cc %in% c( "Undersize(released)      ", "Undersize(released)       ",
                                  " Undersize(released)      " ))
    cc5 = which(recaps$cc %in% c( "4 et 5" ))

    recaps[cc3,"cc"] = "3"
    recaps[ccM,"cc"] = "M"
    recaps[ccD,"cc"] = "Dw"
    recaps[cc5,"cc"] = "5"
    recaps$tagID = gsub("[ -]", "", tolower(recaps$tagID))

    v1 = c("tagID", "timestamp", "lon", "lat", "cw", "ch", "cc", "fisherman", "z.fm", "duro", "Comments")
    recaps = recaps[,..v1]
    names(recaps) = tolower( names(recaps ))

    return(recaps)
  }

}
