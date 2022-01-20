## At the end of each survey year, a data folder for the recent survey needs to be placed in:
## bio.data/bio.snowcrab/data/seabird
## bio.data/bio.snowcrab/data/minilog
## bio.data/bio.snowcrab/data/esonar

## ideally, "click/touch" (manual touchdown / llift off) has already been completed externally with annual results appended to:
## C:/bio.data/bio.snowcrab/data/touchdown/results/clicktouchdown_all.csv
## Any stations without touchdown / liftoff will prompt for manual input when running code below
##


year.assessment = 2021


  if (0) {
    # if debugging:
    loadfunctions('bio.snowcrab')
    require(aegis)
    require(aegis.bathymetry)
    require(bio.snowcrab)

  }

p = bio.snowcrab::load.environment( year.assessment=year.assessment )


# get data tables from Oracle server and store local copies
if (obtain.database.snapshot) {

    if (0) {
      # should you want to store in a temporary location
      yrs=1996:2021
      alt_location = file.path( getwd(), "trawldata" )
      alt_location2 = file.path(getwd(), "logbook")
      alt_location3 = file.path(getwd(), "observer")
      snowcrab.db( DS="set.rawdata.redo", yrs=yrs, fn.root=alt_location ) 
      snowcrab.db( DS="det.rawdata.redo", yrs=yrs, fn.root=alt_location ) 
      snowcrab.db( DS="cat.rawdata.redo", yrs=yrs, fn.root=alt_location ) 
      logbook.db(  DS="rawdata.logbook.redo", yrs=yrs, fn.root=alt_location2 ) 
      observer.db( DS="rawdata.redo", yrs=yrs, fn.root=alt_location3 )
    }
  yrs=1996:2021
  yrs = p$year.assessment
  snowcrab.db( DS="set.rawdata.redo", yrs=yrs ) 
  snowcrab.db( DS="det.rawdata.redo", yrs=yrs ) 
  snowcrab.db( DS="cat.rawdata.redo", yrs=yrs ) 
  logbook.db(  DS="rawdata.logbook.redo", yrs=yrs )  
  logbook.db(  DS="rawdata.licence.redo" )  
  logbook.db(  DS="rawdata.areas.redo" )  
  observer.db( DS="rawdata.redo", yrs=yrs )
}


# -------------------------------------------------------------------------------------
# produce base data files from bio.snowcrab logbook database (marfis) and historical data
# and the observer databases:
# this needs to be done after the above datadumps to refresh locally stored databases

  if (make.fisheries.data) {
    observer.db( DS="odb.redo", p=p ) # 3 minutes
    logbook.db( DS="logbook.redo", p=p )
    logbook.db( DS="logbook.filtered.positions.redo", p=p )
    # fishing ground are used for determination of contraints for interpolation
    logbook.db( DS="fishing.grounds.redo",  p=p )
    logbook.db( DS="logbook.gridded.redo", p=p )
  }


# -------------------------------------------------------------------------------------
# create base set data and add all historical data fixes

  # sequence is important ... do not change
  # creates initial rdata and sqlite db
  snowcrab.db( DS="setInitial.redo", p=p ) # this is required by the seabird.db (but not minilog and netmind)
  snowcrab.db( DS="set.clean.redo", p=p )

    # few sanity checks on the initial data pulled from the raw tables
    problems = data.quality.check( type="stations", p=p)  # duplicates
    problems = data.quality.check( type="count.stations", p=p)
    problems = data.quality.check( type="position", p=p) #MG try checking the end position of the tow, if there is an error

# -------------------------------------------------------------------------------------
# process the net configuration and the temperatures from seabird, netmind, etc ..
# if just updating a single year, run the following, else all years will be run by default
  if (updating.year.assessment) {
    p$seabird.yToload = p$year.assessment
    p$minilog.yToload = p$year.assessment
    p$netmind.yToload = p$year.assessment
    p$esonar.yToload  = p$year.assessment
  }

  # Only do this line at the beginning to create the netmind format. Set 
  # redo.marp... to TRUE to convert the SDS csv file to a usable table 
  # which takes a long time. Only do once!!
  #convert.marport.sds_csv2netmind(yr=yrs, redo.marport_conversion = F )
      
    
  seabird.db( DS="load", Y=p$seabird.yToload ) # this begins 2012;duplicates are often due to seabird files not being restarted each morning

  minilog.db( DS="load", Y=p$minilog.yToload ) # minilog data series "begins" in 1999 -- 60 min?

  #netmind.db( DS='esonar2netmind.conversion',Y=p$esonar.yToload ) #may no longer need to be run BZ
  netmind.db( DS="load", Y=p$netmind.yToload) # netmind data series "begins" in 1998 -- 60 min?
  #JC note: 1998:2002 have about 60 files with no data, just a short header

  p$netmensuration.problems = c( "add trouble id's here") # add troublesome id's here .. eventually move into their respective functions
  seabird.db (DS="stats.redo", Y=p$seabird.yToload )
  minilog.db (DS="stats.redo", Y=p$minilog.yToload )  # note no depth in minilog any more (since 2014 ..  useful for temperature only)
  netmind.db (DS="stats.redo", Y=p$netmind.yToload )



# -------------------------------------------------------------------------------------
# merge in netmind, minilog, seabird, esonar data and do some sanity checks
# Can add any datachecks that might improve overall data quality
# BZ for 2017- If errors repeat and are not actually a problem, create an override

snowcrab.db( DS="set.clean.redo", p=p ) #Updated stats data, need to redo to update stats columns
    #problems = data.quality.check( type="minilog.mismatches", p=p )
    problems = data.quality.check( type="position.difference", p=p)
    problems = data.quality.check( type="minilog.load", p=p)
    problems = data.quality.check( type="minilog.dateproblems", p=p) #track down why ~all sets are giving mismatches
    problems = data.quality.check( type="minilog", p=p)   # Check for duplicate timestamps
    problems = data.quality.check( type="netmind.load", p=p)
    problems = data.quality.check( type="netmind.mismatches", p=p )
    problems = data.quality.check( type="tow.duration", p=p)
    problems = data.quality.check( type="tow.distance", p=p)
    problems = data.quality.check( type="seabird.mismatches", p=p )
    problems = data.quality.check( type="seabird.load", p=p)
    problems = data.quality.check( type="netmind.timestamp" , p=p)

# -------------------------------------------------------------------------------------


  if (0) {
    # local lookup tables are required for the following snowcrab.db steps
    # most do not need to be run, expect perhaps temperature

      pC = bio.snowcrab::snowcrab_parameters( project_class="carstm", yrs=1999:year.assessment )
      
      pB = aegis.bathymetry::bathymetry_parameters( p=parameters_reset(pC), project_class="carstm"  )
      M = aegis.bathymetry::bathymetry_db( p=pB, DS="aggregated_data" , redo=TRUE ) #this step can take ~20 minutes
      
      pS = substrate_parameters( p=parameters_reset(pC), project_class="carstm"  )
      M = aegis.substrate::substrate_db( p=pS, DS="aggregated_data" , redo=TRUE )
      
      pT = temperature_parameters( p=parameters_reset(pC), project_class="carstm"  )
      M = aegis.temperature::temperature_db( p=pT, DS="aggregated_data" , redo=TRUE )
   
  }


# -------------------------------------------------------------------------------------
# Finalize the data sets
# MG det.initial.redo updates and processes morphology. This code now identifies morphology errors, which must be
# checked with written logs, then sent to database and put in debugging here and re-run


  snowcrab.db( DS="det.initial.redo", p=p ) #if no morphology errors exist, you will get an error message. Confirm legitimacy.
  snowcrab.db( DS="det.georeferenced.redo", p=p )

  snowcrab.db( DS="cat.initial.redo", p=p )
  snowcrab.db( DS="cat.georeferenced.redo", p=p )

  snowcrab.db( DS="set.biologicals.redo", p=p )

  require(aegis.temperature)
  snowcrab.db( DS="set.complete.redo", p=p )


# -------------------------------------------------------------------------------------
# update a database of simple transformation ranges, etc.. for plotting range, etc.
  snowcrab.db( DS="data.transforms.redo", p=p)


# -------------------------------------------------------------------------------------
# create some simple/crude timeseries by each CFA
  snowcrab.timeseries.db( DS="observer.redo", p=p )
  snowcrab.timeseries.db( DS="biologicals.redo", p=p )
  snowcrab.timeseries.db( DS="groundfish.t.redo", p=p )
  # snowcrab.timeseries.db( DS="biologicals.2014.redo" )  # reduced subset that matches 2014 station id's ..

  # example: or to get a quick one for a few vars of interest, region of interest ... no saving to file
  snowcrab.timeseries.db( DS="biologicals.direct", p=p, regions='cfa4x', vn=c('R0.mass'), trim=0 )  # returns the data
