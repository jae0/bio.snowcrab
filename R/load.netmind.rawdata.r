
    load.netmind.rawdata = function(fn, f, set) {

      out=NULL
      netmind=NULL

      header = readLines(fn, n=20)
      outfmt = c(dates="year-m-d", times="h:m:s")

      if (length(header) < 20 ) return( out )
      localtime = netmindDate( header=header, outvalue="localtime" )

      tmpfile = "tmp.netmind"

      # remove "*" -- hard to do internally without fileswapping out and reloading
      tmp = readLines(fn)
      tmp = gsub("[*]", "", tmp )
      write( tmp, file=tmpfile )

      # skip 16 because first few records are sometimes incomplete
      netmind = read.table( file=tmpfile, sep="", as.is=T, colClasses="character", header=F, skip=16)
      file.remove(tmpfile)
      nr0 = nrow(netmind)
      if (nr0 < 30 ) return(out)

      nc0 = ncol(netmind)
      if ( nc0 < 12 ) return ( out ) # no net metrics stored ( only position and speed )  -- ignored
      if ( nc0 > 14 ) stop( fn)  # should not be the case unless new data streams start
      if ( nc0 == 12) {
        # older files do not have depths, add a dummy column
        # more modern data have 13 columns 2000 +
        netmind$depth= NA
      }
      #Temperature now available on Marport Sensors
 
      yr = as.numeric(unlist(strsplit(fn, "/"))[grep("archive", unlist(strsplit(fn, "/")))+1])
      if(nc0 == 14 & yr>2020){
        temperature = netmind[,13] 
      }
      else{
        if ( nc0 == 14 & as.numeric(netmind[1,1])<130000) {
          # a few files have wing spread as well
          # more modern data have 13 columns 2000 +
          netmind[,13] = NULL
        }
      }
      if(ncol(netmind) == 13){
        netmind[14] = NA
      }
      ### NOTE:: doorspread is actually wingspread ..
      ### the logging software calls it doorspread but it is actually wingspread.
     
       colnames(netmind) = c("ndate", "ntime", "lat.deg", "lat.min", "lat.orient",
                            "lon.deg", "lon.min", "lon.orient", "speed", "primary", "secondary", "doorspread", "depth", "temperature" )
      
      numbers = c("lat.deg", "lat.min", "lon.deg", "lon.min", "speed", "primary", "secondary", "doorspread", "depth", "temperature")
      netmind = factor2number(netmind, numbers)

      netmind$tinc = 1:nrow(netmind)

      # determine deepest point if possible, using a smoothed depth as variability due to incorrect pings are frequent
      deepest.point = NULL
      if (nc0 %in% c( 13, 14 ) ) {
        require(mgcv)
        z.gam = try ( gam( depth ~ s(tinc, k=5, bs="ts"), data=netmind, optimizer=c("outer", "nlm") ), silent=T )
        if ( ! "try-error" %in% class( z.gam )) {
          netmind$depth.smoothed = predict( z.gam, newdata=netmind, newdata.guaranteed=T  )
          deepest.point = which.max ( netmind$depth.smoothed)
        }
      }

      if ( nc0==12 | length(deepest.point) == 0 ) deepest.point = round( nrow(netmind) / 2 )

      netmind$lon = - (netmind$lon.deg + (netmind$lon.min / 60) )
      netmind$lat =    netmind$lat.deg + (netmind$lat.min / 60)
      netmind = netmind[, c("ndate", "ntime", "lat", "lon", "speed", "primary", "secondary", "doorspread", "depth", "temperature")]
      netmind$ndate = paste(substring(netmind$ndate,1,2), substring(netmind$ndate,3,4), substring(netmind$ndate,5,6), sep="-")
      netmind$ntime = paste(substring(netmind$ntime,1,2), substring(netmind$ntime,3,4), substring(netmind$ntime,5,6), sep=":")
      netmind$timestamp = lubridate::ymd_hms( paste( netmind$ndate, netmind$ntime), tz="UTC" )

      # netmind data stored in GMT/UTC from GPS; the offset varies depending upon season due to daylight savings time (3 or 4 hrs)
      # obtain time offset in hours
      time.offset = netmindDate( header=header, outvalue="timeoffset" )   #  rounded to hours of fractional days
      netmind$timestamp = netmind$timestamp + lubridate::hours(time.offset) # now in local time (America/Halifax)
      netmind.timestamprange = netmind$timestamp
      netmind.timestamp = netmind$timestamp[ deepest.point ]
      yr = lubridate::ymd_hms(netmind.timestamp)


      yr = lubridate::year(yr )
      if(months(netmind.timestamp) %in% c("January")){
        yr = lubridate::year( netmind.timestamp ) -1
      }

      line.localtime = grep("Local Time:", header, ignore.case=T  )
      line.ship = grep("Ship:", header, ignore.case=T  )
      line.comments = grep("Comments:", header, ignore.case=T )


      trip = gsub( "^.*Trip:", "", header[ line.ship ] )
      trip = gsub( "Tow:.*$", "", trip )
      trip = gsub( "[[:space:]]", "", trip )

      if ( ! grepl( "^S[[:digit:]]{8}$", trip, ignore.case=T ) )  { # not a standard code
        dy = paste( "00", as.character( days(netmind.timestamp) ), sep="")
        dy = substring( dy, nchar(dy)-1, nchar(dy) )
        mn = paste( "00", as.character( as.numeric(months(netmind.timestamp))), sep="")
        mn = substring( mn, nchar(mn)-1, nchar(mn) )
        yr = paste( "00", as.character( years(netmind.timestamp) ), sep="")
        yr = substring( yr, nchar(yr)-3, nchar(yr) )
        trip=paste("S", dy, mn, yr, sep="" )
      }

      setno = gsub( "^.*Tow:", "", header[ line.ship ] )
      setno = gsub( "[[:space:]]", "", setno )
      setno =  as.numeric( setno )

      station = unlist(strsplit( header[[1]], "\\", fixed=TRUE ))
      station = station[ length(station) ]
      station = gsub( "[[:alpha:]]", "", station)
      station = gsub( "[[:punct:]]", "", station)

      if(nchar(station)>4){
        station = unlist(strsplit( header[[1]], "/", fixed=TRUE ))
        station = station[ length(station) ]
        station = gsub( "[[:alpha:]]", "", station)
        station = gsub( "[[:punct:]]", "", station)

      }
      station = as.numeric( station )
      setxi = NULL

      if (is.null ( setxi ) ) {
        # check time first
        netmind.date.range = range( netmind.timestamprange )
        sets.in.date.range = which( set$timestamp >= netmind.date.range[1] & set$timestamp <= netmind.date.range[2] )
        if ( length( sets.in.date.range ) == 1 ) setxi= sets.in.date.range
      }

      if (is.null ( setxi ) ) {
        # check time first
        netmind.date.range = range( netmind.timestamprange )
        sets.in.date.range = which( set$timestamp >= netmind.date.range[1] & set$timestamp <= netmind.date.range[2]
          & set$set==setno )
        if ( length( sets.in.date.range ) == 1 ) setxi= sets.in.date.range
      }

      if (is.null ( setxi ) ) {
        # check time and station
        netmind.date.range = range( netmind.timestamprange )
        sets.in.date.range = which( set$timestamp >= netmind.date.range[1] & set$timestamp <= netmind.date.range[2]
          & set$trip==trip & set$station==station & set$set==setno )
        if ( length( sets.in.date.range ) == 1 ) setxi= sets.in.date.range
      }

      if (is.null ( setxi ) ) {
        # check time and station
        netmind.date.range = range( netmind.timestamprange )
        sets.in.date.range = which( set$timestamp >= netmind.date.range[1] & set$timestamp <= netmind.date.range[2]
          & set$station==station  )
        if ( length( sets.in.date.range ) == 1 ) setxi= sets.in.date.range
      }

      if (is.null ( setxi ) ) {
        # check distances
        ni = floor( nrow(netmind) / 2 )
        dx = abs( set$lon - netmind$lon[ni] )
        dy = abs( set$lat - netmind$lat[ni] )
        sets.in.spatial.range = which( dx < 5/60 & dy < 5/60
          & set$yr==yr )  # less than 5 minutes away
        if ( length( sets.in.spatial.range ) == 1 ) setxi= sets.in.spatial.range
      }

      if (is.null ( setxi ) ) {
        # check distances
        ni = floor( nrow(netmind) / 2 )
        dx = abs( set$lon1 - netmind$lon[ni] )
        dy = abs( set$lat1 - netmind$lat[ni] )
        sets.in.spatial.range = which( dx < 5/60 & dy < 5/60
          & set$yr==yr )  # less than 5 minutes away
        if ( length( sets.in.spatial.range ) == 1 ) setxi= sets.in.spatial.range
      }


      if (is.null ( setxi ) ) {
        # check distances
        ni = floor( nrow(netmind) / 2 )
        dx = abs( set$lon - netmind$lon[ni] )
        dy = abs( set$lat - netmind$lat[ni] )
        sets.in.spatial.range = which( dx < 5/60 & dy < 5/60
          & set$trip==trip )  # less than 5 minutes away
        if ( length( sets.in.spatial.range ) == 1 ) setxi= sets.in.spatial.range
      }

      if (is.null ( setxi ) ) {
        # check staion, distance, and time
        ni = floor( nrow(netmind) / 2 )
        dx = abs( set$lon - netmind$lon[ni] )
        dy = abs( set$lat - netmind$lat[ni] )
        sets.in.spatial.range = which( dx < 5/60 & dy < 5/60
          & set$station==station & set$trip==trip & set$set==setno )  # less than 10 minutes away
        if ( length( sets.in.spatial.range ) == 1 ) setxi= sets.in.spatial.range
      }

      if (is.null ( setxi ) ) {
        # check staion, distance, and time
        ni = floor( nrow(netmind) / 2 )
        dx = abs( set$lon - netmind$lon[ni] )
        dy = abs( set$lat - netmind$lat[ni] )
        sets.in.spatial.range = which( dx < 5/60 & dy < 5/60
          & set$station==station & set$trip==trip )  # less than 10 minutes away
        if ( length( sets.in.spatial.range ) == 1 ) setxi= sets.in.spatial.range
      }


      if  (is.null ( setxi ) ) return (NULL) # no matching data -- break

      setx = set[ setxi , ] # matching trip/set/station

      netmind_uid =  paste( "netmind", setx$trip, setx$set, setx$station, lubridate::hour(netmind.timestamp), lubridate::minute(netmind.timestamp), f, sep=".")

      filename = basename(fn)

      line.ship = grep("Ship:", header, ignore.case=T  )
      line.comments = grep("Comments:", header, ignore.case=T )

      comments = gsub("^Comments: ", "", header[ line.comments] )

      netmind$netmind_uid = netmind_uid

      metadata = data.frame( filename, netmind_uid, yr, netmind.timestamp, setx$trip, setx$set, setx$station, setx$Zx, setx$timestamp, comments, stringsAsFactors =FALSE )
      names( metadata ) = c("filename", "netmind_uid", "yr", "netmind_timestamp", "trip", "set", "station", "setZx", "set_timestamp",  "comments" )
      metadata$yr = as.numeric( as.character( metadata$yr ))

      basedata = netmind[ which( !is.na( netmind$netmind_uid) ) ,]

      out = list( metadata=metadata, basedata=basedata )

      return(out)
    }
