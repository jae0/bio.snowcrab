
- monthly basis
- 1/10 degree resolution
- no trap X soak days = trap days


  year.assessment = 2023

  p = bio.snowcrab::load.environment( year.assessment=year.assessment )
  
  yrs = 2019:2023

  # ------------------------------------------
  # Map:  Interpolated mean/geometric mean of various variables in the set data table
  #  p$do.parallel=F
  p$corners = data.frame(plon=c(220, 990), plat=c(4750, 5270) )

  #p$mapyears = year.assessment + c(-5:0 )

  outdir = "~/tmp"
  
  probs=c(0,0.975)
  plot_crs=st_crs( p$aegis_proj4string_planar_km )
   
  require(sf) 
  
  x = logbook.db( DS="logbook" )
  x = x [polygon_inside( x, region="isobath1000m"),]
  x = x[ which(x$effort <= 300) ,]
  x = x[ which(x$cpue < 500),]
  x$landings = x$landings/1000 

  x$year = year(x$date.fished)  
  x$month = month(x$date.fished)   
  
  x = x[ x$year %in% yrs , ]
  x = x[ x$month %in% 1:12 ,]

  x = st_as_sf( x, coords= c("lon", "lat"), crs=st_crs( projection_proj4string("lonlat_wgs84") ) )

  pres = 0.1
  y = as.data.table(st_coordinates(x)  )
  y$X = trunc(y$X / pres) * pres
  y$Y = trunc(y$Y / pres) * pres

  setDT(x)
  x = cbind(x, y)
  x$geometry = NULL
  x$AUID = paste(x$X, x$Y, sep="~")
  
  x = x[,.(AUID, X, Y, year, month, licence, cfv, year, landings, effort, soak.time, cpue)]
  x$uid = paste( x$AUID, x$year, x$month, sep="_" )
  x$soak_time_hrs = as.numeric(x$soak.time)
  x$effort_time = as.numeric(x$soak.time) * as.numeric(x$effort)

  X = x[ , 
      .(
          effort_mean = mean( effort, na.rm=TRUE),
          effort_sd = sd(  effort, na.rm=TRUE),
          soak_time_mean = mean( soak_time_hrs, na.rm=TRUE ),
          soak_time_sd = sd( soak_time_hrs, na.rm=TRUE ),
          no_vessels = length( unique(cfv)),
          effort_time_mean = mean( effort_time ),
          effort_time_sd = sd( effort_time )
      ), 
      by=.(uid)
  ]

  i = which(X$no_vessels <= 5 )
  X = X[i,]

  oo = data.table( matrix(unlist(strsplit(as.character(X$uid), "_")), ncol=3, byrow=T) )
  names(oo) = c("AUID", "year", "month")    
  oo = cbind(uid=X$uid, oo)
  out = oo[X, on=.(uid)]
  
  out$uid = NULL
  out$year = as.numeric(out$year)
  out$month = as.numeric(out$month)
  
  fwrite(out, file="~/tmp/fishery_data.csv")
    
