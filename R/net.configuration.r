
  net.configuration = function( N, t0=NULL, t1=NULL, set_timestamp=NULL, yr=NULL, plotdata=TRUE ) {

    # N is netmind data
    # t0 is current best estimate of start and end time
    # set_timestamp is timestamp from set-log .. alternate if there is no other useful time marker
    # note:: NETMIND data stream calls the wingspread as "DOORSPREAD".
    # It is not. It is wingspread.

    plotdir = project.datadirectory("bio.snowcrab", "data", "netmind", "figures" )

    # create default output should the following fail
    out = data.frame( slon=NA, slat=NA, distance=NA, spread=NA, spread_sd=NA,
      surfacearea=NA, vel=NA, vel_sd=NA, netmind_n=NA, t0=NA, t1=NA, dt=NA, yr=NA, 
      temperature.n=NA, temperature_sd.n=NA )

    n.req = 10

    if (!is.null(t0)) {
      if (length(t0)>1) {t0 = NULL}
      if (is.na(t0)) t0 = NULL
    }

    if (!is.null(t1)) {
      if (is.na(t1)) t1 = NULL
    }

    #changed this switch from depgth filed to lat as if there is not depth info can still run script
	  if ( length( which( is.finite( N$lat))) < n.req ) print(N[1,])

    problem = F

    # time checks
    if ( is.null(t0) ){
      if ( !is.null(set_timestamp) ) {
        t0 = set_timestamp  # no data in t0 ,, use set_timestamp as alternate
        set_timestamp =NULL # remove to cause no more effects
      }
    }

    #bad.list = c('netmind.S26092014.9.541.17.48.304',
    #             'netmind.S20092007.8.333.17.27.241' )
    bad.list = NULL
    bad.list = unique( c(bad.list, p$netmensuration.problems) )
    if ( N$netmind_uid[1] %in% bad.list) {
      return(out)
    }

    if ( any( is.null( t1 ) || is.null(t0) ) )  {
      # try to determine from netmind data if minilog/seadbird data methods have failed. .. not effective due to noise/and small data stream
      print(N$netmind_uid[1])
      M = N[, c("timestamp", "depth") ]

      oo = which( is.finite( M$depth ))
      if ( length(oo) < n.req ) return(out)

      if(!is.null(set_timestamp)) settimestamp=set_timestamp
      if(is.null(set_timestamp)) settimestamp=t0
      time.gate =  list( t0=settimestamp - dminutes(6), t1=settimestamp + dminutes(12) )

      bcp = list(id=N$netmind_uid[1], nr=nrow(M), YR=yr, trip = unlist(strsplit(N$netmind_uid[1], "\\."))[2], tdif.min=3, tdif.max=11, time.gate=time.gate,
                   depth.min=20, depth.range=c(-25,15), eps.depth=3,
                   smooth.windowsize=5, modal.windowsize=5, datasource = "snowcrab",
                   noisefilter.trim=0.05, noisefilter.target.r2=0.8, noisefilter.quants=c(0.05, 0.95) )

      if(yr<2007)bcp$from.manual.archive=FALSE # manual touchdown only done since 2007


      bcp = bottom.contact.parameters( bcp ) # add other default parameters
      #BC: Determine if this station was done yet, if not then we want user interaction.
      if(file.exists(file.path(bcp$from.manual.archive, "clicktouchdown_all.csv"))){
        manualclick = read.csv(file.path(bcp$from.manual.archive, "clicktouchdown_all.csv"), as.is=TRUE)
        station = unlist(strsplit(bcp$id, "\\."))[4]
        sta.ind = which(manualclick$station == station & manualclick$year == bcp$YR)
        if(length(sta.ind) == 1) bcp$user.interaction = FALSE
        else bcp$user.interaction = TRUE
        
      }
      
      dp.out.range = which(M$depth < .1 | M$depth > 360)
      if(length(dp.out.range)>0){
        M$depth[dp.out.range] = NA
      }
      M$wingspread = M$doorspread
      bc = NULL

      bc = bottom.contact( x=M, bcp=bcp )
      if ( is.null(bc) || (!is.null( bc$res) && ( ( !is.finite(bc$res$t0 ) || !is.finite(bc$res$t1 ) ) ) )) {
         bcp$noisefilter.target.r2=0.8
 
         bc = bottom.contact( x=M, bcp=bcp )
      }
      if ( is.null(bc) || (!is.null( bc$res) && ( ( !is.finite(bc$res$t0 ) || !is.finite(bc$res$t1 ) ) ) )) {
        # try once more with random settings
        bcp$noisefilter.inla.h =  0.05
        bcp$noisefilter.target.r2= 0.75
        bc = bottom.contact( x=M, bcp=bcp )
      }
      if ( is.null(bc) || (!is.null( bc$res) && ( ( !is.finite(bc$res$t0 ) || !is.finite(bc$res$t1 ) ) ) )) {
        # try once more with more random settings .. noise is high in netminds data
        bcp$eps.depth = 5
        M$depth = jitter( M$depth, amount = bcp$eps.depth/10 )
        bcp$noisefilter.inla.h =  0.1
        bc = bottom.contact( x=M, bcp=bcp )
      }
      if ( !is.null(bc) )  {
        if (plotdata) {
          bottom.contact.plot( bc )
          plotfn = file.path( plotdir, paste(bcp$id, "pdf", sep="." ) )
          print (plotfn)
          dev.flush()
          dev.copy2pdf( file=plotfn )
          graphics.off()
        }
      }
      if (is.null(t0) & !is.null(bc$bottom0) ) t0 = bc$bottom0
      if (is.null(t1) & !is.null(bc$bottom1) ) t1 = bc$bottom1
      N = N[ bc$bottom.contact , ]
    }

    if (all(is.na(t0))) t0=NA
    if (all(is.na(t1))) t1=NA

    if ( is.null(t1) || !is.finite(t1) ) {
      t1_tmp = NA  # rather than guessing, flag and then fill later
    } else {
      t1_tmp = t1
    }

    # if we are here, it is because a reasonable amount of data is present ..
    # do some more checks and get a first estimate of some parameters in case other errors/limits are found
    out$slon=N$lon[1]
    out$slat=N$lat[1]
    out$spread=mean( N$doorspread, na.rm=T ) / 1000  # though called doorspread it is actually wingspread
    out$spread_sd=sd( N$doorspread, na.rm=T ) /1000

    # first set the impossible door spreads to NA and then interpolate based upon a 5 and 95 % quantiles
    # Netmind has very noisy data
    hl = c( 3, 16 ) # m
    ihl = which( N$doorspread < hl[1] | N$doorspread > hl[2] )
    if (length (ihl)>0 ) N$doorspread[ihl] = NA
    quantiles.to.trim = c(0.05, 0.95)
    qnt =  quantile( N$doorspread[N$doorspread>0], quantiles.to.trim, na.rm=T)
    iqnt = which( N$doorspread < qnt[1] | N$doorspread > qnt[2] )
    N$doorspread[ iqnt ] = NA
    gooddoor =  which(   is.finite(N$doorspread) )
    baddoor =   which( ! is.finite(N$doorspread) )

    if (length(baddoor)>0) N$doorspread[ baddoor ] = NA
    if ( length( gooddoor) < n.req ) {
      #problem =T turned this off December 23, 2013 as if we have decent gps data now we can calc distance
      if ( length( gooddoor) >0 ) {
      	out$netmind_n <- length(gooddoor)
        out$spread = mean( N$doorspread[ gooddoor ], na.rm=T ) / 1000
        out$spread_sd = sd( N$doorspread[ gooddoor ], na.rm=T ) / 1000
      }
    }

    if ( length(t0)>1) t0 = NA
    out$t0 = t0
    out$t1 = t1

    out$dt = difftime( as.POSIXct(t1), as.POSIXct(t0) ) # minutes

    if(!is.na(out$t0))    out$yr = lubridate::year( out$t0)
    if(is.na(out$t0))    out$yr = lubridate::year(N$timestamp[1])
   
    itime =  which( N$timestamp >= t0  &  N$timestamp <= t1 )
    if ( length( itime) < n.req ) problem = T

    if (problem) {
      ## Need to check on these at some point
      out$t0 = t0
      out$t1 = NA
      out$dt = NA
      return(out)
    }

    # eOR checks
    
    
      # Higher GPS data accuracy in 2019 resulted in high approximation of distance travelled. This was due to frequent occasions of backwards
      # motion (assuming from boat roll and pitch) which was added to the total distance. When looking at previous years this also seem to have
      # but not to the same extent however removing the backwards motion would result in better data accuracy from year using e-sonar 2014+.
      # The solution is to trim the data by 3/4, or, with the current recording rate of every 4 seconds instead of every 1 second.
      
      # The above works when the metrics are repeated but not when we actually only record the values when we get a sensor hit. This does not 
      # work for the 2021 data. Trying to smooth latitude and longitude so that we dont have data loss.
      
      # This again did not work as paths were way out in some cases. Resolved by applying a windowed average smoothing
      # Keeping code for st_simplify in case
    itime.smoothing =  which( N$timestamp >= (t0 - lubridate::seconds(10))  &  N$timestamp <= (t1 + lubridate::seconds(10) ))
    
    rolling.ave = 6
    Ntemp = N[ itime.smoothing,]
    Ntemp$nlat = NA
    Ntemp$nlong = NA
    for(i in 1:nrow(Ntemp)){
      inds = (i-rolling.ave):(i+rolling.ave)
      inds = inds[which(inds>0)]
      inds = inds[which(inds<=nrow(Ntemp))]
      Ntemp$nlong[i] = mean(Ntemp$lon[inds])
      Ntemp$nlat[i] = mean(Ntemp$lat[inds])
    }
    
    
    #plot(Ntemp$lon,Ntemp$lat, type = "l", col = "#000000E6", lwd = 3)
    #lines(Ntemp$nlon,Ntemp$nlat, col = 'red', lwd = 2)
    
    itime =  which( Ntemp$timestamp >= t0  &  Ntemp$timestamp <= t1 )
    
    N = Ntemp[ itime, ]
    
      if(yr>=2014){
        
        
   #      require(sf)
   #      tolerence = 5
   #      N$new.lat = NA
   #      N$new.lon = NA
   #      xy = st_cast( st_multipoint(cbind(N$lon,N$lat), dim = "XY"), 'LINESTRING' ) 
   #      cxyll = st_sfc(xy, crs = "+init=epsg:4326" )
   #      c_xy = st_transform(cxyll, st_crs("+proj=utm +zone=20 +ellps=GRS80 +datum=NAD83 +units=m" ))
   #      N$utm.y = c_xy[[1]][,2]
   #      N$utm.x = c_xy[[1]][,1]
   #      xy_smoothed = st_simplify( c_xy, dTolerance=tolerence)
   #      tol.count = tolerence+1
   #      while(any(grepl("MULTILINESTRING", class(xy_smoothed)))){
   #        tol.count = tol.count+1
   #        print(paste("Moving tolerence to:", tol.count))
   #        xy_smoothed = st_simplify( c_xy, dTolerance=tol.count)
   #        if(tol.count > 100) break()
   #      }
   #      
   #   
   #      
   #      nodes =  which(paste(N$utm.y, N$utm.x) %in% paste(xy_smoothed[[1]][,2],xy_smoothed[[1]][,1]))
   #      nodes = nodes[!nodes %in% (nodes-1)]
   #      if(2 %in% nodes){
   #        nodes[which(nodes == 2)] = 1
   #      }
   #       important.nodes = nodes   
   #       
   #      n.xylon = xy_smoothed[[1]][,1]
   #      n.xylat = xy_smoothed[[1]][,2]
   #      if(!1%in%important.nodes){
   #        n.xylon = c(N$utm.x[1], n.xylon)
   #        n.xylat = c(N$utm.y[1], n.xylat)
   #        important.nodes = c(1,  important.nodes) 
   #      }
   #      if(!nrow(N)%in%important.nodes){
   #        n.xylon = c(n.xylon, N$utm.x[nrow(N)])
   #        n.xylat = c(n.xylat, N$utm.y[nrow(N)])
   #        important.nodes = c(important.nodes, nrow(N)) 
   #      }
   #      
   #      for(i in 1:(length(important.nodes)-1)){
   # 
   #      # from = which(paste(round(N$lon, 10), round(N$lat, 10)) == paste(round(n.xylon[i],10), round(n.xylat[i], 10)))
   #      #  to = which(paste(round(N$lon, 10), round(N$lat, 10)) == paste(round(n.xylon[i+1],10), round(n.xylat[i+1], 10)))
   #        from = important.nodes[i]
   #        to = important.nodes[i+1]
   #        fill.n = to-from
   #  
   #        tes = cbind(c(N$utm.y[important.nodes[i]],N$utm.y[important.nodes[i+1]]) , c(N$utm.x[important.nodes[i]], N$utm.x[important.nodes[i+1]]))
   #        tes.seg = st_cast( st_multipoint(tes, dim = "XY"), 'LINESTRING' ) 
   #        to.fill = st_sample(tes.seg, size = fill.n, type = "regular")
   #        N$new.lon[(from+1):to] = to.fill[[1]][,1]
   #        N$new.lat[(from+1):to] = to.fill[[1]][,2] 
   #      }
   #      N$new.lat[important.nodes] = n.xylon
   #      N$new.lon[important.nodes] = n.xylat
   #      
   #      
   #      xyt = st_cast( st_multipoint(cbind(N$new.lat, N$new.lon), dim = "XY"), 'LINESTRING' ) 
   #      cxyllt = st_sfc(xyt, crs = "+proj=utm +zone=20 +ellps=GRS80 +datum=NAD83 +units=m" )
   #      herex = st_transform(cxyllt, CRS("+init=epsg:4326"))
   #      
   #      
   # N$old.lat = N$lat
   # N$old.lon = N$lon
   # N$lat = herex[[1]][,2]
   # N$lon = herex[[1]][,1]

        N$lat = N$nlat
        N$lon = N$nlong
        
      }
      
      
      if(nrow(N)>1) {
   
        # t1 gives the approximate time of net lift-off from bottom
        # this is not good enough as there is a potential backdrift period before net lift off
        # this means distance trawled calculations must use the geo-positioning of the boat
        # at maximum tension and not the position at net movemnent up
        # (otherwise, this would introduce a net bias towards smaller swept areas).
        pos = c( "lon", "lat" )
        distance.from.start = geosphere::distCosine(N[1, pos], N[, pos])/1000  #  in km^2
        end = which.max( distance.from.start)
        if( !is.finite(end) ) end=nrow(N)
        n = N[ 1:end , ]

    
	      out$vel = mean(n$speed, na.rm=T, trim=0.1)
	      out$vel_sd = sd(n$speed, na.rm=T)
	      out$slon=n$lon[1]
	      out$slat=n$lat[1]
	      out$netmind_n=end

      delta.distance = NULL
	    n$distance = NA

	  
    if (nrow(n) > 10) {
  
      # integrate area:: piece-wise integration is used as there is curvature of the fishing track (just in case)
      delta.distance = geosphere::distGeo ( n[ 1:(end-1), pos ], n[ 2:end, pos ] ) / 1000 ## in meters convert to km
      n$distances = c( 0, cumsum( delta.distance  ) ) # cumsum used to do piecewise integration of distance

      
      # model/smooth/interpolate the spreads
      n$doorspread.predicted = NA
      ii = which( is.finite( n$doorspread ) )
      if ( length(ii) > n.req ) {
          # recall that doorspread is actually wingspread ...
         nn = as.data.table( n[ii,] )
         nn = nn[, median(doorspread, na.rm=TRUE), by=(distances)]
         n$doorspread.predicted = approx( x=nn$distances, y=nn$V1, xout=n$distances, method="linear", rule=2, na.rm=TRUE )$y
       #turned off gam model in December 20, 2013 giving unrealistic values for spread as the
       # new esnoar files have 0 and NA whereas older netmind are filled with previous value
			      	#gam.model = try( gam( doorspread ~ s(distances, k=5, bs="ts"), data=n[ii,], optimizer=c("outer", "nlm")), silent = T )
			        #if ( ! "try-error" %in% class( gam.model )) {
			        #  n$doorspread.predicted = predict( gam.model, newdata=n, newdata.guaranteed=T )
			        #} else {
			        #  n$doorspread.predicted = approx( x=n$distances, y=n$doorspread, xout=n$distances, method="linear", rule=2 )$y
			        #}
      }
      if ( length( which( is.finite( n$doorspread.predicted ) ) ) < 10 ) {
        n$doorspread.predicted = mean( n$doorspread , na.rm=T, trim=0.1 )
      }
      mean.doorspreads = ( n$doorspread.predicted[1:(end-1)] + n$doorspread.predicted[2:(end)] ) / 2 / 1000  # mean between two ends
      partial.area =  delta.distance * mean.doorspreads
      out$surfacearea = sum( partial.area )  # km^2
      out$surfacearea = abs(  out$surfacearea )

      out$spread = mean(n$doorspread.predicted, na.rm=T, trim=0.1)/1000  # in km
      spread_sd = sd(n$doorspread.predicted, na.rm=T )/1000
      if(!is.na(spread_sd) & spread_sd!=0) out$spread_sd = spread_sd #if just using the mean from above do not over write spread_sd
      out$distance=n$distances[end]

      # 2022 Temperature now available on Marport Sensors. Keep last 3/4 of bottom contact to help parse out latency  
      temp.to.keep = n$temperature[(length(n$temperature)/4):length(n$temperature)]
      if ( length(which(!is.na(temp.to.keep))) > 4){
        out$temperature.n = mean(temp.to.keep, na.rm=T, trim=0.1)
        out$temperature_sd.n = sd(temp.to.keep, na.rm=T)
      }
      else{
        out$temperature.n = NA
        out$temperature_sd.n = NA
      }
      
     
   }
      }
     print(out)
    return (out)
  }


