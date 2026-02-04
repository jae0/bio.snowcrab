-- Initial request in October 2024  

•	Data was aggregated to 10 km x 10 km grids and every 2 weeks for the years 2022 and 2023
•	The data format for "fishery_data.csv" was:
o	AUID
o	year
o	wk
o	effort_mean
o	soak_time_mean
o	Longitude
o	Latitude
 
•	And the labels description was:
•        AUID is the primary key to match to polygons in the shapefile.
•        year= year of fishing event
•        wk= week number ( discretized to every 2 weeks)
•        effort_mean = mean number of trap hauls ( x 10"3 trap hauls per 10x10 km grid)
•        soak_time_mean = mean soak time (hrs)
•        Lon, Lat -- centroid location


-- Additional information requested in Jan/Feb 2026:

    Please provide this data for the years 2024 and 2025 with the following additional requests:

o	Can the total number of trap hauls also be provided for each grid?  This could replace the ‘effort-mean’.

o	The total number of ‘trap-days’ for each grid.  This would be a calculation made as the product of 
traps hauled multiplied by soak time in days of each of the original landing records, 
which should be summed, not averaged, when producing grid cell summaries.

o	The October 2024 data did contain some records where soak days was not provided.  
In this case, can the total number of trap-hauls without soak-days be provided?  
This would permit an estimation of the number of trap-days by using the ‘mean-soak time’.
 

# - monthly basis
# - 1/10 degree resolution
# - no trap X soak days = trap days


  yrs = 2024:2025

  year.assessment = 2025

  p = bio.snowcrab::load.environment( year.assessment=year.assessment )

  p$corners = data.frame(plon=c(220, 990), plat=c(4750, 5270) )

  outdir = "~/tmp"
  
  probs=c(0,0.975)
  plot_crs=st_crs( p$aegis_proj4string_planar_km )
   
  require(sf) 
  
  x = logbook.db( DS="logbook" )
 
  x = x [polygon_inside( x, region="isobath1000m"),]
  x = x[ which(x$effort <= 300) ,]
  x = x[ which(x$cpue < 500),]
  
  x$landings = x$landings/1000  # 10^3 trap hauls 

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
  x$uid = paste( x$AUID, x$year, x$month, sep="_" )
  x$soak_time_days = as.numeric(x$soak.time/24)

  X = x[ , 
      .(
          effort_sum = sum( effort, na.rm=TRUE),  # 10^3 trap hauls per 10x10 km grid 
          trap_days_sum = sum( soak_time_days * effort , na.rm=TRUE), # 10^3 trap haul days per 10x10 km grid 
          soak_time_days_sum = sum( soak_time_days, na.rm=TRUE )   # total soak days per grid
      ), 
      by=.(uid)
  ]

  oo = data.table( matrix(unlist(strsplit(as.character(X$uid), "_")), ncol=3, byrow=T) )
  names(oo) = c("AUID", "year", "month")    
  oo = cbind(uid=X$uid, oo)
  out = oo[X, on=.(uid)]
  
  out$uid = NULL
  out$year = as.numeric(out$year)
  out$month = as.numeric(out$month)
  
  fwrite(out, file="~/tmp/fishery_data.csv")
    
          
