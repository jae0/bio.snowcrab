
# Snow crab assessment data

This is a simple markdown document (without Quarto or Rmarkdown variations). It can be viewed in formatted form via:
  
  - a web browser open the webpage: [01_snowcrab_data.md](https://github.com/jae0/bio.snowcrab/tree/master/inst/markdown/01_snowcrab_data.md)   

  - a web browser open the local file directly: [01_snowcrab_data.md](../markdown/01_snowcrab_data.md) (you might need to install a browser add-in), or 
  
  - a text editor with markdown awareness (e.g. Rstudio, VSCode, etc. ).

This document is the roadmap of the full snow crab assessment. It is a nearly linear sequence as most of this sequence can be run automatically. But stop in the section *Prepare external data for lookup processing and spatiotemporal modelling* (~line 535) where external data need to be also created and merged back into this stream before completion of the full database as outlined below that section. 


## Prepare the data environment

  - Ensure the **year.assessment** is correct (to update for 4X in Feb, increment to capture new yaers logbook data)..

  - make sure your passwords are correct in the externally stored password file.

  - sequence is important due to multiple merges of data from various sources .. try not to skip steps in this script ... :)


```r

require(aegis)  # basic helper tools

year.assessment = 2025  # change this as appropriate

p = bio.snowcrab::load.environment( year.assessment=year.assessment )  # set up initial settings

# if something goes wrong run:  rlang::last_trace() # to show the trace

# only required if mapping
additional_features = snowcrab_mapping_features(p, plot_crs=projection_proj4string("lonlat_wgs84")) # ggplot background objects



```

## Core snow crab trawl survey data

First, ensure data have been added to the back-end storage relational databases (ISDB). This used to be double key-punched by the At-sea observer company. Currently this is loaded manually from electronic records by Brent Cameron. 

{**Brent**, please add a short description here, with links to your work, when you can}


After loading to back-end storage, we download the whole data system locally to manipulate and error check. There are several ways to proceed. 

  - Run in an R-session in **Linux** (on the server): This is challenging as ROracle does not compile easily. 

  - Run in an R-session in **MSWindows** with minimal dependencies (DBI, ROracle) to get all \*.rawdata.redo files below in one pass (after which you can skip the \*.rawdata.redo steps below). See:  
  
    - [99_data_extraction_jae.R](https://github.com/jae0/bio.snowcrab/tree/master/inst/scripts/99_data_extraction_jae.R)  

  - Run in an R-session in **MSWindows**: This requires you have bio.snowcrab and supporting libraries installed. Just run the following:

#### Snow crab trawl survey data

```r  

# choose appropriate years
yrs = p$year.assessment  # redo just last year
# yrs = 1996:p$year.assessment # redo all years
# yrs = p$year.assessment + -3:0 # redo last 4 years

snowcrab.db( DS="set.rawdata.redo", yrs=yrs ) # set data 
snowcrab.db( DS="det.rawdata.redo", yrs=yrs ) # determined data (individual sex, size, mat, etc)
snowcrab.db( DS="cat.rawdata.redo", yrs=yrs ) # aggregated weights and counts by set


```


#### At-sea-observed fishery data  

The storage is in the same data base system as the survey data: ISSDB

```r

# bring in raw data from ISDB backend as annual snapshots
observer.db( DS="rawdata.redo", yrs=yrs )
observer.db( DS="bycatch.redo", yrs=yrs )  # this is also an SQL call

# compute items of interest
observer.db( DS="odb.redo", p=p ) # 3 minutes
observer.db( DS="bycatch_clean_data.redo", p=p, yrs=p$yrs ) # 3 minutes

snowcrab.timeseries.db( DS="observer.redo", p=p, mau="region" )
snowcrab.timeseries.db( DS="observer.redo", p=p, mau="subarea" )

ts_years = 2004:p$year.assessment
ts_outdir = file.path( p$annual.results, "timeseries", "observer")
figure.timeseries.observer(p=p, outdir=ts_outdir, plotyears=ts_years, mau="region") # all variables
figure.timeseries.observer(p=p, outdir=ts_outdir, plotyears=ts_years, mau="subarea" ) # all variables



p$yrs_observer = c(p$year.assessment + c(-4:0))  # if you change this change yrs_observer 
# p$yrs_observer = p$yrs

# map them here as a quick check:
loc =  project.datadirectory("bio.snowcrab", "output", "maps", "observer.locations" )
map.observer.locations( p=p, basedir=loc, years=p$yrs_observer )  
 

# Map large charmismatic/megafauna in at sea obseved data:
# loc =  project.datadirectory("bio.snowcrab", "output", "maps", "observer.entanglements" ) 
# map.observer.entanglements(p=p, basedir=loc, years=p$yrs_observer, region = "cfaall" ) 



```


#### Fishery logbook data

Produce base data files from bio.snowcrab logbook database (marfis) and historical data

```r

# Update local local fisheries data table from quota report website 
# https://inter-j02.dfo-mpo.gc.ca/mqr/quotareports/snowcrab?rptyear=&year&rptnote=false&lang=en
web_fisheriesdata_update()

o = CA.getTable("TAC_SUBAREA")
o = CA.getTable("TAC_REGION")

fd = CA.getTable("fisheriesdata")  # check
str(fd)

# save locally to allow off-line running
o = snowcrab_tacs( "region", redo=TRUE )
o = snowcrab_tacs( "subarea", redo=TRUE )


# bring in raw data from back-end MARFIS databases as annual snapshots
logbook.db( DS="rawdata.logbook.redo", yrs=yrs ) 
logbook.db( DS="rawdata.licence.redo" ) 
logbook.db( DS="rawdata.areas.redo" )  

# clean up and format
logbook.db( DS="logbook.redo", p=p )
logbook.db( DS="logbook.filtered.positions.redo", p=p )

# if (0) {
  # deprecated ... ?
  # fishing ground are used for determination of constraints for interpolation (no longer used?)

  logbook.db( DS="fishing.grounds.redo",  p=p )
  logbook.db( DS="logbook.gridded.redo", p=p )
# }

# create summaries for differing area designations and time intervals (yearly, monthly, weekly, etc for reports)
# NOTE: this step requires the TAC database to be updated 

# -->>> note: labels and internal codes for Management Areak Units (maus)
# mau="region"  -> cfanorth, cfasouth, cfa4x
# mau="subarea" -> cfanorth, cfa23, cfa24, cfa4x

( maus = management_areal_units( mau="region" )   )  # to see contents

o = fishery_data( mau="region", redo=TRUE )  # region referes to "region" in logbook -- 3 areas
o = fishery_data( mau="subarea", redo=TRUE ) # subarea referes to "subarea" in logbook -- 4 areas


# map logbooks here as a quick check:
yrsplot = p$year.assessment + -3:0
# yrsplot = p$yrs

loc = project.datadirectory("bio.snowcrab", "output", "maps", "logbook.locations" )

# map all logbooks locations by year in whole domain
map.logbook.locations( p=p, basedir=loc, years=yrsplot )


# map last recent (two years) of logbook locations, by subarea and region 
lasttwoyears = p$year.assessment + (-1:0)
outdir = project.datadirectory("bio.snowcrab", "output", "maps", "logbook.locations", "recent" )
map.logbook.locations.by.area( p=p, basedir=outdir, years=lasttwoyears, mau="region"  ) 
map.logbook.locations.by.area( p=p, basedir=outdir, years=lasttwoyears, mau="subarea"  ) 


# timeseries:
fpts_loc = file.path( p$annual.results,  "timeseries", "fishery")  

figure.landings.timeseries( yearmax=p$year.assessment, outdir=fpts_loc, outfile="landings.ts",  plotmethod="withinset", mau="region" ) # with inset 
figure.effort.timeseries( yearmax=p$year.assessment, outdir=fpts_loc, outfile="effort.ts", mau="region"  )
figure.cpue.timeseries( yearmax=p$year.assessment, outdir=fpts_loc, outfile="cpue.ts", mau="region"  )


# save alternate with 23 and 24 split in another location: "subarea"
figure.landings.timeseries( yearmax=p$year.assessment, outdir=fpts_loc, outfile="landings.ts", mau="subarea" )
figure.effort.timeseries( yearmax=p$year.assessment, outdir=fpts_loc, outfile="effort.ts", mau="subarea"   )
figure.cpue.timeseries( yearmax=p$year.assessment, outdir=fpts_loc, outfile="cpue.ts", mau="subarea"  )


# seasonal figures
fishery_data_figures_seasonal(
  time_resolution="summary_weekly", 
  toget=c("cummulative_landings", "cummulative_effort", "cpue" ), 
  outdir=fpts_loc, 
  mau="region" 
)

fishery_data_figures_seasonal(
  time_resolution="summary_weekly", 
  toget=c("cummulative_landings", "cummulative_effort", "cpue" ), 
  outdir=fpts_loc, 
  mau="subarea" 
) 


# maps of fishery performance

fp_loc = file.path( p$project.outputdir, "maps", "logbook", "snowcrab", "annual" )

map.fisheries.data( 
  outdir=fp_loc, 
  probs=c(0,0.975), 
  additional_features=additional_features,
  outformat="png"
)

map.fisheries.data.alllocations(p=p, additional_features=additional_features)  # all locations

# singletons used for FSAR
figure.fisheries.timeseries( outdir=fp_loc, mau="region", region_id="cfanorth" ) 
figure.fisheries.timeseries( outdir=fp_loc, mau="region", region_id="cfasouth" ) 
figure.fisheries.timeseries( outdir=fp_loc, mau="region", region_id="cfa4x" ) 



```


#### Create base set-level information with fixes for year-specific issues, etc. 

This initial set-level information is required to sanity check net mensuration estimates:

```r

snowcrab.db( DS="setInitial.redo", p=p ) # this is required by the seabird.db (but not minilog and netmind)

# few sanity checks on the setInitial data pulled from the raw tables
problems = data.quality.check( type="stations", p=p)  # duplicates?
problems = data.quality.check( type="count.stations", p=p)
problems = data.quality.check( type="position", p=p) #MG try checking the end position of the tow, if there is an error
problems = data.quality.check( type="position.difference", p=p)


```


#### Net mensuration

Net configuration (and temperature, depth) has been recorded using several systems over the years. Currently, the following now tracks only net configuration. Temperature and depth are processed separately. 

At the end of each survey year, a data folder for the recent survey needs to be placed in:

  - bio.data/bio.snowcrab/data/seabird
  - bio.data/bio.snowcrab/data/minilog
  - bio.data/bio.snowcrab/data/esonar

Ideally, "click/touch" (manual touchdown / lift off) has already been completed externally with annual results appended to:

  - C:/bio.data/bio.snowcrab/data/touchdown/results/clicktouchdown_all.csv

TODO (anyone):  convert to an annual breakdown such that historical data is static and updates are simply additional files specific to the survey year

```
"station","start","end","depth","year","trip"
"125","2007-11-13 10:25:11","2007-11-13 10:33:53",144.7,2007,"S13112007"
```

This file contains all the manually determined start and stop times of each trip and set from 2007 to present. Any stations without touchdown / liftoff will use automatically determined start/stop times (see: bio.netmesuration).


{**Brent**, would you please add links to your system here and a short description or link to a Tech Doc}.

If updating a single year or a subset of years, specify the appropriate years here:

```r
yrs_to_load = p$year.assessment  # to operate on the year of assessment only
# yrs_to_load = 2007:p$year.assessment  # to redo a range of years 

p$seabird.yToload = yrs_to_load
p$minilog.yToload = yrs_to_load
p$netmind.yToload = yrs_to_load
p$esonar.yToload  = yrs_to_load

 
# seabird data begins in 2012; duplicates are often due to seabird files not being restarted each morning
seabird.db( DS="load", Y=p$seabird.yToload ) 

# minilog data series "begins" in 1999  .. no longer used .. 
minilog.db( DS="load", Y=p$minilog.yToload ) 

# netmind data series "begins" in 1998 
# JC note: 1998:2002 have about 60 files with no data, just a short header
netmind.db( DS="load", Y=p$netmind.yToload ) 


# add troublesome id's to the list here .. eventually move them into their respective functions
p$netmensuration.problems = c() 

# now compute statistics (means, sd, etc) where necessary
# this next section looks to be a bottleneck in terms of speed (esp seabird.db() ).. 
# any way to alter low-level data storage mechanism/merges to avoid redoing things? (Brent, Kate?)

fn_clicktouchdown_all = file.path( p$data_root, "data", "touchdown", "results", "clicktouchdown_all.csv" )

seabird.db( DS="stats.redo", Y=p$seabird.yToload, 
  fn_clicktouchdown_all=fn_clicktouchdown_all, skip_computations=TRUE )  # skip to only use manually determined solutions
minilog.db( DS="stats.redo", Y=p$minilog.yToload, 
  fn_clicktouchdown_all=fn_clicktouchdown_all, skip_computations=TRUE )  # note no depth in minilog any more (since 2014 ..  useful for temperature only)
netmind.db( DS="stats.redo", Y=p$netmind.yToload, 
  fn_clicktouchdown_all=fn_clicktouchdown_all, skip_computations=TRUE )


```


#### Merging data and compute aggregate variables

Merge in netmind, minilog, seabird, esonar data and store in set-level information ... and do some sanity checks.

NOTE: There are many historical data elements that require manual fixes to the raw data.


```r

snowcrab.db( DS="set.clean.redo", p=p ) 

# map them here as a quick check:
yrsplot = p$year.assessment + -1:0
loc = project.datadirectory("bio.snowcrab", "output", "maps", "survey.locations" )
map.survey.locations( p=p, basedir=loc, years=yrsplot ) # uses setClean

 
# Identify any issues with set.clean .. these need to be checked 

# minilogs no longer used: test not required
# problems = data.quality.check( type="minilog.mismatches", p=p )
# problems = data.quality.check( type="minilog.load", p=p)
# problems = data.quality.check( type="minilog.dateproblems", p=p) #track down why ~all sets are giving mismatches
# problems = data.quality.check( type="minilog", p=p)   # Check for duplicate timestamps

problems = data.quality.check( type="netmind.load", p=p)
problems = data.quality.check( type="netmind.mismatches", p=p )
problems = data.quality.check( type="netmind.timestamp" , p=p)


problems = data.quality.check( type="seabird.mismatches", p=p )
problems = data.quality.check( type="seabird.load", p=p)

problems = data.quality.check( type="tow.duration", p=p)
problems = data.quality.check( type="tow.distance", p=p)

snowcrab.db( DS="det.initial.redo", p=p ) 

# fishno or crabno not found ....
problems = data.quality.check( type="biologicals_fishno", p=p)
problems = data.quality.check( type="biologicals_duplicates", p=p)

# QA/QC of morphology (individuual level data)
# sanity check: identify morphology errors: fix them if any are found and re-run until satisfied.. 
problems = data.quality.check( type="biologicals_morphology", p=p)
names(problems)

snowcrab.db( DS="det.georeferenced.redo", p=p )  # merge set.clean


# sanity check aggregate catches
snowcrab.db( DS="cat.initial.redo", p=p )
snowcrab.db( DS="cat.georeferenced.redo", p=p )  # merge set.clean

# aggregate catches by snow crab categories via det.initial, cat.initial
snowcrab.db( DS="set.biologicals.redo", p=p )

# update a database of simple transformation ranges, etc.. for plotting range, etc.
snowcrab.db( DS="data.transforms.redo", p=p) 

# create some simple/crude timeseries by each CFA using set.biologicals and observer data 
snowcrab.timeseries.db( DS="biologicals.redo", p=p, mau="region" )
snowcrab.timeseries.db( DS="biologicals.redo", p=p, mau="subarea" )


```

End of core snow crab data assimilation.

#### Size-frequency 


Size frequency of carapace condition of at-sea-observed data

```r

outdir_obs_sizefreq = file.path( p$annual.results, "figures", "size.freq", "observer") 

figure.observed.size.freq( mau="region", years="all",  outdir=outdir_obs_sizefreq )

figure.observed.size.freq( mau="subarea", years="all",  outdir=outdir_obs_sizefreq )
 

```


Size frequency of **mature male** fraction from survey

```r
outdir_survey_sizefreq = file.path( p$annual.results, "figures", "size.freq", "carapacecondition" ) 

figure.sizefreq.carapacecondition( 
  X = snowcrab.db( p=p, DS="det.georeferenced" ), 
  cwbr=4, 
  vbar=95, 
  mau="region", 
  outdir= outdir_survey_sizefreq 
)  

figure.sizefreq.carapacecondition( 
  X = snowcrab.db( p=p, DS="det.georeferenced" ), 
  cwbr=4, 
  vbar=95, 
  mau="subarea", 
  outdir=outdir_survey_sizefreq  
)  


```
 
Size frequency distributions of snow crab carapace width from trawl data, broken down by maturity classes. This uses "set.clean" and "det.initial".

```r

# discretize size and compute arithmetic (den) and geometric mean (denl) areal densities by groups of years 

# first merge det and set with some QA/QC and area designations
M = size_distributions( p=p, toget="rawdata", mau="region", redo=TRUE )  

M = size_distributions( p=p, toget="rawdata", mau="subarea", redo=TRUE )  


 
# NENS, SENS, 4X
create_size_frequencies(
  p=p, 
  mau="region", 
  outdir=file.path( p$annual.results, "figures", "size.freq", "survey" ) 
)     

# NENS, CFA23, CFA24, 4X
create_size_frequencies(
  p=p, 
  mau="subarea" ,
  outdir = file.path( p$annual.results, "figures", "size.freq", "survey" )
)    


```
 

#### Survey-related timeseries

Generate some simple/generic timeseries before entry into the main Markdown reports. We do them here instead as they need to be created only once:
 
```r

ts_years = 2004:p$year.assessment
ts_outdir = file.path( p$annual.results, "timeseries", "survey")
  

figure.timeseries.survey(p=p, outdir=ts_outdir, plotyears=ts_years, mau="region") # all variables
figure.timeseries.survey(p=p, outdir=ts_outdir, plotyears=ts_years, mau="subarea" ) # all variables
 
# potential predators
species_predator = c(10, 11, 30, 40, 50, 64, 201, 202, 204, 610, 1203 )
bc_vars = c(paste("ms.mass", species_predator, sep='.'), paste("ms.no", species_predator, sep='.'))
figure.timeseries.survey(p=p, outdir=ts_outdir, plotyears=ts_years, variables=bc_vars, mau="region" )
figure.timeseries.survey(p=p, outdir=ts_outdir, plotyears=ts_years, variables=bc_vars, mau="subarea" ) # all variables

# potential competitors
species_competitors = c( 2521, 2511, 2211, 2523 ) #  2523=N stone crab
bc_vars = c(paste("ms.mass", species_competitors, sep='.'), paste("ms.no", species_competitors, sep='.'))
figure.timeseries.survey(p=p, outdir=ts_outdir, plotyears=ts_years, variables=bc_vars, mau="region" )
figure.timeseries.survey(p=p, outdir=ts_outdir, plotyears=ts_years, variables=bc_vars, mau="subarea" ) # all variables
  

over_ride_default_scaling = TRUE
if (over_ride_default_scaling) {

  figure.timeseries.survey(p=p, outdir=ts_outdir, plotyears=ts_years, variables="R0.mass", mau="region" ) 
  figure.timeseries.survey(p=p, outdir=ts_outdir, plotyears=ts_years, variables="R0.mass", mau="subarea" ) # all variables
  
  figure.timeseries.survey(p=p, outdir=ts_outdir, plotyears=ts_years, variables=c("sexratio.all","sexratio.mat","sexratio.imm"), mau="region" )
  figure.timeseries.survey(p=p, outdir=ts_outdir, plotyears=ts_years, variables=c("sexratio.all","sexratio.mat","sexratio.imm"), mau="subarea" ) # all variables
  
  figure.timeseries.survey(p=p, outdir=ts_outdir, plotyears=ts_years, variables="cw.male.mat.mean", backtransform=TRUE, mau="region" ) 
  figure.timeseries.survey(p=p, outdir=ts_outdir, plotyears=ts_years, variables="cw.male.mat.mean", backtransform=TRUE, mau="subarea" ) # all variables
  
}

create_deprecated_figures = FALSE
if (create_deprecated_figures) {
    
  # no longer relevant (incomplete) as temp is now created through temp db. and not in gshyd
  figure.timeseries.survey(p=p, outdir=ts_outdir, plotyears=ts_years,type='groundfish.t') # groundfish survey temperature

}




```


#### Survey-related maps

Need to generate some simple maps (mostly thin plate splines; previously, multilevel B-splines but latter seems to be poor extrapolations) before entry into the main Markdown reports. We do them here instead as they are very slow and need to be created only once:


```r
 
map_outdir = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual" )
map_years  = p$year.assessment + c(0:-3)

# pre-construct a prediction surface with depths to make mapping filter faster in map.set.information
predlocs = get_predlocs(p=p, redo=TRUE)  


create_maps_for_the_road_show = FALSE
if (create_maps_for_the_road_show) {
    # used? ... seems redundant ... probably delete ..
    road_show_vars = c("totmass.male.com", "totmass.female.mat", "R0.mass" )
    map.set.information( p=p, outdir=map_outdir, mapyears=map_years, variables=road_show_vars )

    map.set.information( p=p, outdir=map_outdir, mapyears=map_years, variables='t', log.variable=FALSE, theta=35)
}

# variables that should not be logged
set = snowcrab.db( p=p, DS="set.biologicals")
variables = bio.snowcrab::snowcrab.variablelist("all.data")
variables = intersect( variables, names(set) )

ratio_vars = c("sexratio.all", "sexratio.mat", "sexratio.imm")
nolog.variables = c("t", "z", "julian", variables[grep("cw", variables)])
log.variables = variables[ !variables %in% nolog.variables ]

mass_vars = log.variables[ grep('mass', log.variables)]
no_vars = log.variables[ grep('no', log.variables)]

map.set.information( p=p, outdir=map_outdir, mapyears=map_years, variables=nolog.variables, log.variable=FALSE, theta=35)

# logit transform for ratios
map.set.information( p=p, outdir=map_outdir, mapyears=map_years, variables=ratio_vars, log.variable=FALSE, theta=40)

map.set.information( p=p, outdir=map_outdir, mapyears=map_years, variables= mass_vars )

map.set.information( p=p, outdir=map_outdir, mapyears=map_years, variables= no_vars,  probs=c(0,0.975))


# **potential** predators
species_predator = c(10, 11, 30, 40, 50, 201, 202, 204 )
bc_vars = c(paste("ms.mass", species_predator, sep='.'), paste("ms.no", species_predator, sep='.'))
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )
map.set.information( p=p, outdir=outdir_bc, mapyears=map_years, variables=bc_vars, probs=c(0,0.975)) 


# **potential** competitors
species_competitors = c( 2521, 2511, 2211, 2523 ) # 2523=N stone crab
bc_vars = c(paste("ms.mass", species_competitors, sep='.'), paste("ms.no", species_competitors, sep='.'))
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )
map.set.information( p=p, outdir=outdir_bc, mapyears=map_years, variables=bc_vars, probs=c(0,0.975)) 


# all variables (geometric means)
# map.set.information( p, outdir=map_outdir) # takes a long time ... run only vars of interest

```

#### Create simple Markdown reports (Quarto)

Basic data QA/QC, assimilation is complete.

Now, generate a list of results for use in the 02_fishery_summary.md reports:

```r
# prepare data items and tables for quick load into reports:

o = snowcrab_load_key_results_to_memory( p=p, todo="fishery_results", mau="region",  redo=TRUE )  
o = snowcrab_load_key_results_to_memory( p=p, todo="fishery_results", mau="subarea", redo=TRUE ) 

o = snowcrab_load_key_results_to_memory( p=p, todo="survey", mau="region",  redo=TRUE )  
o = snowcrab_load_key_results_to_memory( p=p, todo="survey", mau="subarea", redo=TRUE ) 


```

The fishery related report can now be created in: 

- [02_fishery_summary.md](https://github.com/jae0/bio.snowcrab/tree/master/inst/markdown/02_fishery_summary.md) 

This can be done by running the following at a command prompt. It uses 'make' which follows the recipe in the associated Makefile found inside the markdown directory. It is repetitive. 

(Someone: please make the call shorter (it will require some minor shell programming).)


```bash
# ensure years and mau are correct (note they are both repeated twice): 

make quarto FN=02_fishery_summary.md YR=2025 MAU=region DATADIR=~/bio.data/bio.snowcrab DOCTYPE=html PARAMS="-P year_assessment:2025 -P mau:region -P todo:[fishery_results]" --directory=~/bio/bio.snowcrab/inst/markdown 


make quarto FN=02_fishery_summary.md YR=2025 MAU=subarea DATADIR=~/bio.data/bio.snowcrab DOCTYPE=html PARAMS="-P year_assessment:2025 -P mau:subarea -P todo:[fishery_results]" --directory=~/bio/bio.snowcrab/inst/markdown 

 
 
```


Basic trawl survey-related report can now be created in:

- [02_survey_summary.md](https://github.com/jae0/bio.snowcrab/tree/master/inst/markdown/02_survey_summary.md) 


```bash

# ensure years  and mau are correct (note they are both repeated twice): 
make quarto FN=02_survey_summary.md YR=2025 MAU=region DATADIR=~/bio.data/bio.snowcrab DOCTYPE=html PARAMS="-P year_assessment:2025 -P mau:region -P todo:[survey]" --directory=~/bio/bio.snowcrab/inst/markdown

# ensure years  and mau are correct (note they are both repeated twice): 
make quarto FN=02_survey_summary.md YR=2025 MAU=subarea DATADIR=~/bio.data/bio.snowcrab DOCTYPE=html PARAMS="-P year_assessment:2025 -P mau:subarea -P todo:[survey]" --directory=~/bio/bio.snowcrab/inst/markdown


```


## Prepare external data for lookup processing and spatiotemporal modelling

To this point all survey data have been assimilated. Next section prepares additional EXTERNAL DATA sources to act as lookup data for modelling purposes.


Local empirical lookup tables are required for the steps that follow. Most do not need to be re-run, UNLESS you want reset the base data for each project. 


#### Polygon data base

If this your first time using this system, you will need to see if the polygons defined are good enough for your use. If you need/want to (re)-initialize polygons, they can be found in:
  
  - [aegis.polygons/inst/scripts/01_polygons.md](https://github.com/jae0/aegis.polygons/tree/master/inst/scripts/01_polygons.md)  
  
  - [aegis.coastline/inst/scripts/01_coastline.md](https://github.com/jae0/aegis.coastline/tree/master/inst/scripts/01_coastline.md)  

NOTE: If you over-write them, then you will need to make sure that all the other polygon-based methods are consistent (projections,domain definition, etc). Beware.
 

NOTE: If there is a need to alter/create new polygons used for modelling size frequencies and abundance estimation, this can be done here:  
 

```r
# create areal_units (polygons) for biomass estimation and size structure
# this does not need to be run ... 
# ... if you do, then all analyses that use it (carstm-based results) would need to be re-run

you_are_sure_you_want_to_recreate_polygons = FALSE

if (you_are_sure_you_want_to_recreate_polygons) {

ps = snowcrab_parameters(
  project_class = "carstm",
  yrs = 1999:year.assessment,   
  areal_units_type = "tesselation",
  carstm_model_label = paste( "default", "fb", sep="_" )  # default for 'fb' (fishable biomass)
)

# this step uses "set.clean" 
xydata = snowcrab.db( p=ps, DS="areal_units_input", redo=TRUE )
# xydata = snowcrab.db( p=ps, DS="areal_units_input" )

# for mapping (background) .. redo=TRUE if resetting colours etc
additional_features = snowcrab_mapping_features(ps, redo=FALSE )  

# create constrained polygons with neighbourhood as an attribute
sppoly = areal_units( p=ps, xydata=xydata, n_iter_drop=0, redo=TRUE, verbose=TRUE )  # this needs to match carstm related parameters in snowcrab_parameters

# sppoly=areal_units( p=ps )  # to reload

plot(sppoly["AUID"])

sppoly = st_transform(sppoly, st_crs( projection_proj4string("lonlat_wgs84") ))
sppoly$dummyvar = ""

xydata = st_as_sf( xydata, coords=c("lon","lat") )
st_crs(xydata) = st_crs( projection_proj4string("lonlat_wgs84") )


# to map 
o = ggplot() + 
  geom_sf( data=sppoly, aes(), colour="slateblue", lwd=0.2, alpha=0.76)  +
  geom_sf( data=xydata, aes(), colour="darkblue", alpha=0.9, size=0.9 ) +
  coord_sf(xlim =ps$corners$lon, ylim =ps$corners$lat, expand = FALSE) +
  additional_features 

print(o)

} 

``` 

#### Bathymetry

The empirical lookup of **aggregated_data**, should be updated reasonably frequently. Every ~3-5 years? See the bathymetry project for more details:

  - [aegis.bathymetry/inst/scripts/01_bathymetry_data.r](https://github.com/jae0/aegis.bathymetry/tree/master/inst/scripts//01_bathymetry_data.R)

The model solutions for bathymetry comes from [**stmv**](https://github.com/jae0/stmv) smoothed solutions. There is no need to re-run the **carstm** solutions unless you want to replace the lookup system with them. 


#### Substrate grainsize

The empirical lookup of aggregated data however is useful and should be updated reasonably frequently. Every ~3-5 years?  See the substrate project for more details:

  - [aegis.substrate/inst/scripts/01_substrate_data.R](https://github.com/jae0/aegis.substrate/tree/master/inst/scripts/01_substrate_data.R)

 
The model solutions for substrate grain size comes from **stmv** smoothed solutions. There is no need to re-run the **carstm** solutions unless you want to replace the lookup system with them. 


#### Bottom temperature

Base data comes from many sources. Historical data was assimilated by OSD and copied here. Recent data are now coming from a separate project and stored in its' own relational database (lead by Amy Glass). 

{***Amy**, please add a short description and link to Tech Doc here}

After loading, we download annual time slices here for processing. See for more details:

  - [aegis.temperature/inst/scripts/01_temperature_data.md](https://github.com/jae0/aegis.temperature/tree/master/inst/scripts/01_temperature_data.md)
  

Bottom temperature modelling should however be re-run annually after assimilating new year's data. Currently, we use:
  
  - [aegis.temperature/inst/scripts/03_temperature_carstm.md](https://github.com/jae0/aegis.temperature/tree/master/inst/scripts/03_temperature_carstm.md)

NOTE: this is a long process and requires a lot of RAM (aim for 128 GB RAM and a large swap space) ... 1-2 days. Not sure if it will run on tethys as it only has 64GB RAM.
  

##### Area-specific timeseries of bottom temperatures:

```r
figure_area_based_extraction_from_carstm(DS="temperature", year.assessment )  # tis can only do done once we have an sppoly for snow crab and temperature/carstm has been completed
```



## Finalize data sets for modelling

We finalize the snow crab data after completing empirical data lookups of covariates (bathymetry and temperatures), in preparation for modelling.

```r

snowcrab.db( DS="set.complete.redo", p=p ) # note depth is log transformed here 


```




## Data from other surveys required for Index modelling with carstm

#### Groundfish

Namely groundfish surveys provide supplemental data for species composition analysis:

  - [aegis.survey/inst/scripts/01_survey_data.md](https://github.com/jae0/aegis.survey/tree/master/inst/scripts/01_survey_data.md).


#### Species composition 

Species composition is a covariate in the areal model for snow crab. To obtain a predictive surface, the following needs to be run: 
    
  - [aegis.speciescomposition/inst/scripts/03_speciescomposition_carstm.md](https://github.com/jae0/aegis.speciescomposition/tree/master/inst/scripts/03_speciescomposition_carstm.md)


## Continue with Biomass index modelling, Fishery model estimation, Quarto Reports and summaries used for CSAS documents

Once the above are complete, we can continue with biomass estimation etc.

  - [Index modelling: bio.snowcrab/inst/markdown/03_biomass_index_carstm.md](https://github.com/jae0/bio.snowcrab/tree/master/inst/markdown/03_biomass_index_carstm.md)

  - [Quarto report of ecosystem summary: bio.snowcrab/inst/markdown/04_ecosystem_summary.md](https://github.com/jae0/bio.snowcrab/tree/master/inst/markdown/04_ecosystem_summary.md)
 
  - [Biomass modelling via Julia: bio.snowcrab/inst/markdown/05_snowcrab_fishery_model_turing.md](https://github.com/jae0/bio.snowcrab/tree/master/inst/markdown/05_snowcrab_fishery_model_turing.md)
 
  - [Quarto report of fishery model results: bio.snowcrab/inst/markdown/06_fishery_model_summary.md](https://github.com/jae0/bio.snowcrab/tree/master/inst/markdown/06_fishery_model_summary.md)
 


## -- END --


## Deprecated

Unused but stored in case you need similar functionality ... 

```r

snowcrab.timeseries.db( DS="groundfish.t.redo", p=p )  # deprecated to be removed shortly
snowcrab.timeseries.db( DS="biologicals", drop=2014 )  # reduced subset that matches 2014 station id's  
snowcrab.timeseries.db( DS="biologicals", vn="R0.mass" )  # one-off
 
```
