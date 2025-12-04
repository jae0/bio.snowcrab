

snowcrab.timeseries.db = function( DS="default", set=NULL, p=NULL, mau="region",
  trim=0, vn=NULL, sdci=FALSE, drop=NULL, save_results = TRUE ) {

  if (is.null(p)) p = bio.snowcrab::snowcrab_parameters()

  maus = management_areal_units(mau)  # labels etc
  auids = c( maus[["internal"]], "cfaall")
  
  tsoutdir = file.path( p$project.outputdir, "timeseries" )
  dir.create(tsoutdir, showWarnings=FALSE, recursive=TRUE)

  if (DS == "default") return( snowcrab.timeseries.db( DS="biologicals", p=p) )
 
  if (DS %in% c( "biologicals", "biologicals.redo" ) ) {

    fn = file.path( tsoutdir, paste("snowcrab_timeseries_", mau, ".rdz", sep="") )
    if (DS=="biologicals") {
      tsdata = NULL
      if (file.exists( fn) ) tsdata = read_write_fast(fn)
      return(tsdata)
    }

    dat = snowcrab.db( DS ="set.biologicals", p=p )

    dat$year = as.character(dat$yr)
    yrs = sort(unique(dat$yr))

    # only save if default .. if vn is supplied treat as a one-off
    if (is.null(vn)) {
      vn = setdiff( colnames(dat), c("trip", "set", "set_type", "station", "lon1", "lat1", 
        "towquality", "lon", "lat", "plon", "plat", "seabird_uid", "minilog_uid", "netmind_uid" ) )
    } else {
      save_results=FALSE    
    }
 
    lookup.table = snowcrab.db( p=p, DS="data.transforms" )

    
    tsdata = timeseries_simple( dat, auids, yrs, vn, lookup.table=lookup.table ) 

    if (save_results) {
      read_write_fast(data=tsdata, fn=fn )
      return( fn )
    } else {
      return(tsdata)
    }
  }


  # -------------------


  if (DS %in% c( "biologicals.2014", "biologicals.2014.redo" ) ) {
    #\\ "reduced" subset of stations found in 2014 ... to be comparable with smaller survey

    fn = file.path( tsoutdir, paste("snowcrab_timeseries_2014", mau, ".rdz", sep="") )
    if (DS=="biologicals.2014") {
      tsdata = NULL
      if (file.exists( fn) ) tsdata = read_write_fast(fn)
      return(tsdata)
    }


    dat = snowcrab.db( DS ="set.biologicals", p=p )

    dat$year = as.character(dat$yr)
    yrs = sort(unique(dat$yr))

    if (!is.null(drop)) { 
      dat = dat[ which(dat$station %in% unique( dat$station[ which(yr==drop) ] ) ),]
    }

    # only save if default .. if vn is supplied treat as a one-off
    if (is.null(vn)) {
      vn = setdiff( colnames(dat), c("trip", "set", "set_type", "station", "lon1", "lat1", 
        "towquality", "lon", "lat", "plon", "plat", "seabird_uid", "minilog_uid", "netmind_uid" ) )
    } else {
      save_results=FALSE    
    }

    lookup.table = snowcrab.db( p=p, DS="data.transforms" )

    tsdata = timeseries_simple( dat, auids, yrs, vn, lookup.table=lookup.table ) 
 
    return( tsdata )
  }


  # -------------------


  if (DS %in% c( "observer", "observer.redo" ) ) {
    #\\ "reduced" subset of stations found in 2014 ... to be comparable with smaller survey
    fn = file.path( tsoutdir, paste("snowcrab_observer_timeseries_", mau, ".rdz", sep="") )
    if (DS=="observer") {
      tsdata = NULL
      if (file.exists( fn) ) tsdata = read_write_fast(fn)
      return(tsdata)
    }
 
    dat = observer.db( DS="odb" )
    dat$yr = dat$fishyr #Bz March 2019- We want the fishing year (2018/19= 2018), not the calendar year of catch
    dat$year = as.character(dat$yr)
    #dat = dat[ which( dat$cw >= 95),] #BZ March 2019- We want to keep all observed animals, not just cw>95
    vn = c( "cw", "totmass", "abdomen", "chela", "shell", "durometer",  "cpue.kg.trap", "mass", "mat" )
    yrs = sort(unique(dat$yr))
  
    lookup.table = snowcrab.db( p=p, DS="data.transforms" )
 
    tsdata = timeseries_simple( dat, auids, yrs, vn, lookup.table=lookup.table ) 
 
    read_write_fast(data=tsdata, fn=fn )
    return( fn)
  }

  # -------------------


  if (DS %in% c( "groundfish.t", "groundfish.t.redo" ) ) {
    #\\ "reduced" subset of stations found in 2014 ... to be comparable with smaller survey
    fn = file.path( tsoutdir, "groundfish.t.rdz" )
    if (DS=="groundfish.t") {
      tsdata = NULL
      if (file.exists( fn) ) tsdata = read_write_fast(fn)
      return(tsdata)
    }
    tsdata = data.frame(r=NA,yrs=NA,V3='t',meanval=NA,se=NA,n=NA,ub=NA,lb=NA)
    h = groundfish_survey_db( DS='gshyd')
    g = groundfish_survey_db(DS='gsinf')
    g = g[,c('id','sdate','lon','lat','bottom_temperature')]
    h = h[,c('id','temp')]
    names(h)[2] <- 'bottom_temperature'
    f <- merge(g,h,by='id',all.x=T)
    i <- which(is.na(f$bottom_temperature.x) & !is.na(f$bottom_temperature.y))
    f[i,'bottom_temperature.x'] <- f[i,'bottom_temperature.y']
    f$yr <- as.numeric(format(f$sdate,'%Y'))
    f <- fishing.area.designations(f)
    ar <- unique(f$cfa)
    yy <- unique(f$yr)
    yy <- yy[order(yy)]
    for (r in ar) {
      for (yrs in yy) {
        y = f[which(f$yr == yrs & f$cfa ==r & !is.na(f$bottom_temperature.x)),]
        if(nrow(y)>3) {
          ym <- min(y$bottom_temperature.x[y$bottom_temperature.x>0])
          q = log(y$bottom_temperature.x+ym)
          m =  mean (q, na.rm=T)
          n = length(q)
          se = sd(q, na.rm=T)/ sqrt(n-1)
          meanval = exp(m)-ym
          ub = exp(m+se*1.96)-ym
          lb = exp(m-se*1.96)-ym
          j = as.data.frame(cbind(r, yrs, 't',meanval, se, n, ub, lb))
          tsdata <- rbind(tsdata,j)
        }
      }
    }
    #browser()
    colnames(tsdata) = c("region", "year", "variable","mean", "se","n", "ub", "lb")
    numbers = c("year", "mean", "se", "n", "ub", "lb")
    tsdata = factor2number(tsdata, numbers)
    tsdata <- tsdata[!is.na(tsdata$year),]

    read_write_fast(data=tsdata, fn=fn )
    return(fn)
  }

}
