## At the end of each survey year, a data folder for the recent survey needs to be placed in:
## bio.data/bio.snowcrab/data/seabird
## bio.data/bio.snowcrab/data/minilog
## bio.data/bio.snowcrab/data/esonar

## ideally, "click/touch" (manual touchdown / llift off) has already been completed externally with annual results appended to:
## C:/bio.data/bio.snowcrab/data/touchdown/results/clicktouchdown_all.csv
## Any stations without touchdown / liftoff will prompt for manual input when running code below
##

  require(aegis)
  
  year.assessment = 2023

  p = bio.snowcrab::load.environment( year.assessment=year.assessment )

  if ( assimilate_rawdata_from_dfo_databases_to_local_rdata ){
    # choose years to do a data dump
    yrs = 1996:year.assessment # redo all years
    yrs = p$year.assessment  # redo just last year

    snowcrab.db( DS="set.rawdata.redo", yrs=yrs ) 
    snowcrab.db( DS="det.rawdata.redo", yrs=yrs ) 
    snowcrab.db( DS="cat.rawdata.redo", yrs=yrs ) 
    logbook.db(  DS="rawdata.logbook.redo", yrs=yrs )  
    logbook.db(  DS="rawdata.licence.redo" )  
    logbook.db(  DS="rawdata.areas.redo" )  
    observer.db( DS="rawdata.redo", yrs=yrs )
    observer.db( DS="bycatch.redo", yrs=yrs )  
  
  }


# -------------------------------------------------------------------------------------
# produce base data files from bio.snowcrab logbook database (marfis) and historical data
# and the observer databases:
# this needs to be done after the above datadumps to refresh locally stored databases

  if (make.fisheries.data) {
    observer.db( DS="odb.redo", p=p ) # 3 minutes
    observer.db( DS="bycatch_clean_data.redo", p=p, yrs=p$yrs ) # 3 minutes

    logbook.db( DS="logbook.redo", p=p )
    logbook.db( DS="logbook.filtered.positions.redo", p=p )
    # fishing ground are used for determination of contraints for interpolation
    logbook.db( DS="fishing.grounds.redo",  p=p )
    logbook.db( DS="logbook.gridded.redo", p=p )
  }


# -------------------------------------------------------------------------------------
# create base set data and add all historical data fixes
  # creates initial rdata and sqlite db
    snowcrab.db( DS="setInitial.redo", p=p ) # this is required by the seabird.db (but not minilog and netmind)

    # few sanity checks on the setInitial data pulled from the raw tables
    problems = data.quality.check( type="stations", p=p)  # duplicates
    problems = data.quality.check( type="count.stations", p=p)
    problems = data.quality.check( type="position", p=p) #MG try checking the end position of the tow, if there is an error
    problems = data.quality.check( type="position.difference", p=p)
  

# -------------------------------------------------------------------------------------
# process the net configuration and the temperatures from seabird, netmind, etc ..
  
  if (netmensuration) {

    # TODO:: this section needs to be moved to a static database system with incremental additions each year 

    if (updating.year.assessment) {
    # if just updating a single year, run the following, else all years will be run by default
      p$seabird.yToload = p$year.assessment
      p$minilog.yToload = p$year.assessment
      p$netmind.yToload = p$year.assessment
      p$esonar.yToload  = p$year.assessment
    }

    # Only do this line at the beginning to create the netmind format. Set 
    # redo.marp... to TRUE to convert the SDS csv file to a usable table 
    # which takes a long time. Only do once!!
    # convert.marport.sds_csv2netmind(yr=yrs, redo.marport_conversion = FALSE )
      
    seabird.db( DS="load", Y=p$seabird.yToload ) # this begins 2012;duplicates are often due to seabird files not being restarted each morning

    minilog.db( DS="load", Y=p$minilog.yToload ) # minilog data series "begins" in 1999 -- 60 min?

    #netmind.db( DS='esonar2netmind.conversion',Y=p$esonar.yToload ) #may no longer need to be run BZ
    netmind.db( DS="load", Y=p$netmind.yToload) # netmind data series "begins" in 1998 -- 60 min?
    #JC note: 1998:2002 have about 60 files with no data, just a short header
    
    p$netmensuration.problems = c( "add trouble id's here") # add troublesome id's here .. eventually move into their respective functions

    seabird.db (DS="stats.redo", Y=p$seabird.yToload )
    minilog.db (DS="stats.redo", Y=p$minilog.yToload )  # note no depth in minilog any more (since 2014 ..  useful for temperature only)
    netmind.db (DS="stats.redo", Y=p$netmind.yToload )
  
  }


# -------------------------------------------------------------------------------------
# merge in netmind, minilog, seabird, esonar data and do some sanity checks
  
  snowcrab.db( DS="set.clean.redo", p=p ) 
  
    # Identify any issues with set.clean
    # problems = data.quality.check( type="minilog.mismatches", p=p )
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
# QA/QC of morphology and catches

  # sanity check  morphology 
  # identifies morphology errors (they are written to logs), 
  # fix them and re-run .. if no morphology errors exist, you will get an error message. Confirm legitimacy.
  snowcrab.db( DS="det.initial.redo", p=p ) 

  # sanity check catches
  snowcrab.db( DS="cat.initial.redo", p=p )


# -------------------------------------------------------------------------------------

  if (0) {
    # this is a note to remind you: 
    # local empirical lookup tables are required for the following snowcrab.db steps
    # most do not need to be re-run, UNLESS you want reset the base data for each project
    # .. then run 01_* data steps for each respective project
    
      pC = bio.snowcrab::snowcrab_parameters( project_class="carstm", yrs=1999:year.assessment, areal_units_type="tesselation" )
      
      pB = aegis.bathymetry::bathymetry_parameters( p=parameters_reset(pC), project_class="carstm"  )
      M = aegis.bathymetry::bathymetry_db( p=pB, DS="aggregated_data" , redo=TRUE ) #this step can take ~20 minutes

    # ALSO: if this your first time around:
    # you will need to (re)-initialize polygons, etc:
    # these are found in:
    #  aegis.polygons -- 01_polygons.r 
    #  aegis.coastline -- 01_coastline.R

      pS = substrate_parameters( p=parameters_reset(pC), project_class="carstm"  )
      M = aegis.substrate::substrate_db( p=pS, DS="aggregated_data" , redo=TRUE )
      
    # temperature should however be re-run annually after assimilating new year's data
    # by running 01_temperature_data.R (relevant parts are copied below):
      pT = temperature_parameters( p=parameters_reset(pC), project_class="carstm"  )
      temperature_db( DS="bottom.annual.rawdata.redo", p=p, yr=1970:year.assessment )  # brent and amy's new db view
      o = temperature_db( DS="bottom.annual.redo", p=p,   yr=1970:year.assessment ) # use all years to improve spatial resolution 
      o = aegis.temperature::temperature_db( p=pT, DS="aggregated_data" , redo=TRUE )
   
    # also run:  aegis.survey  :: 01_survey_data.R to  02_survey_assimilate.r

    # this is for index building for 03.snowcrab*  :  aegis.speciescomposition:: 01_speciescomposition_carstm_default.R
      
  }


# -------------------------------------------------------------------------------------
# Finalize the data sets

  snowcrab.db( DS="det.georeferenced.redo", p=p )
  snowcrab.db( DS="cat.georeferenced.redo", p=p )
  snowcrab.db( DS="set.biologicals.redo", p=p )
  snowcrab.db( DS="set.complete.redo", p=p ) # note depth is log transformed here

  # this one is may be phased out as it is for quick generic plots... not used any more?
  snowcrab.db( DS="data.transforms.redo", p=p) # update a database of simple transformation ranges, etc.. for plotting range, etc.


# -------------------------------------------------------------------------------------
# create some simple/crude timeseries by each CFA using set.complete 
# -- TODO convert to data.table for speed /clarity
  snowcrab.timeseries.db( DS="observer.redo", p=p )
  snowcrab.timeseries.db( DS="biologicals.redo", p=p )


# -------------------------------------------------------------------------------------
# create areal_units (polygons) for biomass and size structure
   
  ps = snowcrab_parameters(
    project_class="carstm",
    yrs=1999:year.assessment,   
    areal_units_type="tesselation",
    carstm_model_label=  paste( "default", "fb", sep="_" )  # default for fb (fishable biomass)
  )
  
  xydata = snowcrab.db( p=ps, DS="areal_units_input", redo=TRUE )
  # xydata = snowcrab.db( p=ps, DS="areal_units_input" )


  additional_features = snowcrab_mapping_features(ps, redo=FALSE )  # for mapping (background) .. redo=TRUE if resetting colours etc
   

  # create constrained polygons with neighbourhood as an attribute
  sppoly = areal_units( p=ps, xydata=xydata, spbuffer=3, n_iter_drop=0, redo=TRUE, verbose=TRUE )  # this needs to match carstm related parameters in snowcrab_parameters

  # sppoly=areal_units( p=ps )

      plot(sppoly["AUID"])

      sppoly$dummyvar = ""
      xydata = st_as_sf( xydata, coords=c("lon","lat") )
      st_crs(xydata) = st_crs( projection_proj4string("lonlat_wgs84") )

      tmap_mode("plot")
      
      plt = 
        tm_shape(sppoly) +
          tm_borders(col = "slategray", alpha = 0.5, lwd = 0.5) + 
          tm_shape( xydata ) + tm_sf() +
          additional_features[["tmap"]] +
          tm_compass(position = c("right", "TOP"), size = 1.5) +
          tm_scale_bar(position = c("RIGHT", "BOTTOM"), width =0.1, text.size = 0.5) +
          tm_layout(frame = FALSE, scale = 2) +
          tm_shape( st_transform(polygons_rnaturalearth(), st_crs(sppoly) )) + 
          tm_borders(col = "slategray", alpha = 0.5, lwd = 0.5)

      dev.new(width=14, height=8, pointsize=20)
      plt


  


# Size-frequency distributions of snow crab cw from trawl data, broken down by maturity classes

  sppoly = areal_units(  p=ps )  # it uses sppoly
  xrange = c(10, 150)  # size range (CW)
  dx = 2 # width of carapace with discretization to produce "cwd"

  M = size_distributions(p=p, toget="base_data", pg=sppoly, xrange=xrange, dx=dx, redo=TRUE)

  # tabulate... non-zero counts ... must use add_zeros=TRUE to add them, on the fly
  M = size_distributions(p=p, toget="tabulated_data", xrange=xrange, dx=dx, redo=TRUE)

  years = as.character( c(-9:0) + year.assessment )
  M = size_distributions(p=p, toget="simple_direct", xrange=xrange, dx=dx, Y=years, redo=TRUE)



# -------------------------------------------------------------------------------------
# deprecated

  # snowcrab.timeseries.db( DS="groundfish.t.redo", p=p )  # deprecated to be removed shortly
  # snowcrab.timeseries.db( DS="biologicals.2014.redo" )  # reduced subset that matches 2014 station id's .. # deprecated

  # example: to get a quick view of a few vars of interest, region of interest ... no saving to file, but does return the data for viewing
  # snowcrab.timeseries.db( DS="biologicals.direct", p=p, regions='cfa4x', vn=c('R0.mass'), trim=0 )  


