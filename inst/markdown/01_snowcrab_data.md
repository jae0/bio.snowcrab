# Snow crab assessment data

This is a markdown document. It can be viewed in formatted form via:
  
  - a web browser open the webpage: [01_snowcrab_data.md](https://github.com/jae0/bio.snowcrab/tree/master/inst/markdown/01_snowcrab_data.md)   

  - a web browser open the local file directly: [01_snowcrab_data.md](../markdown/01_snowcrab_data.md) (you might need to install a browser add-in), or 
  
  - a text editor with markdown awareness (e.g. Rstudio, VSCode, etc. ).



## Prepare the data environment

  - Ensure the **year.assessment** is correct

  - make sure your passwords are correct in the externally stored password file.

  - sequence is important due to multiple merges of data from various sources .. try not to skip steps in this script ... :)


```r

require(aegis)  # basic helper tools

year.assessment = 2024  # change this as appropriate

p = bio.snowcrab::load.environment( year.assessment=year.assessment )  # set up initial settings

```

## Core snow crab trawl survey data

First, ensure data have been added to the back-end storage relational databases (ISSDB). This used to be double key-punched by the At-sea observer company. Currently this is loaded manually from electronic records by Brent Cameron. 

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
# yrs = 1996:year.assessment # redo all years

snowcrab.db( DS="set.rawdata.redo", yrs=yrs ) # set data 
snowcrab.db( DS="det.rawdata.redo", yrs=yrs ) # determined data (individual sex, size, mat, etc)
snowcrab.db( DS="cat.rawdata.redo", yrs=yrs ) # aggregated weights and counts by set

```


#### At-sea-observed fishery data  

The storage is in the same data base system as the survey data: ISSDB

```r
# bring in raw data from ISSDB backend as annual snapshots
observer.db( DS="rawdata.redo", yrs=yrs )
observer.db
# compute items of interest
observer.db( DS="bycatch.redo", yrs=yrs )  
observer.db( DS="odb.redo", p=p ) # 3 minutes
observer.db( DS="bycatch_clean_data.redo", p=p, yrs=p$yrs ) # 3 minutes

```


#### Fishery logbook data

Produce base data files from bio.snowcrab logbook database (marfis) and historical data

```r

# bring in raw data from back-end MARFIS databases as annual snapshots
logbook.db(  DS="rawdata.logbook.redo", yrs=yrs )  
logbook.db(  DS="rawdata.licence.redo" )  
logbook.db(  DS="rawdata.areas.redo" )  

# clean up and format
logbook.db( DS="logbook.redo", p=p )
logbook.db( DS="logbook.filtered.positions.redo", p=p )

# fishing ground are used for determination of contraints for interpolation (no longer used?)
logbook.db( DS="fishing.grounds.redo",  p=p )
logbook.db( DS="logbook.gridded.redo", p=p )

```


#### Create base set-level information with fixes for year-specific issues, etc. 

This initial set-level information is required to sanity check net mensuration estimates:

```r

snowcrab.db( DS="setInitial.redo", p=p ) # this is required by the seabird.db (but not minilog and netmind)

# few sanity checks on the setInitial data pulled from the raw tables
problems = data.quality.check( type="stations", p=p)  # duplicates
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


Any stations without touchdown / liftoff will prompt for manual input when running code below:


{**Brent**, would you please add links to your system here and a short description or link to a Tech Doc}.

If updating a single year or a subset of years, specify the appropriate years here:

```r
yrs_to_load = p$year.assessment  # to operate on the year of assessment only
# yrs_load = 2020:2024   # to redo a range of years 

p$seabird.yToload = yrs_to_load
p$minilog.yToload = yrs_to_load
p$netmind.yToload = yrs_to_load
p$esonar.yToload  = yrs_to_load

 
# seabird data begins in 2012; duplicates are often due to seabird files not being restarted each morning
seabird.db( DS="load", Y=p$seabird.yToload ) 

# minilog data series "begins" in 1999  
minilog.db( DS="load", Y=p$minilog.yToload ) 

# netmind data series "begins" in 1998 
# JC note: 1998:2002 have about 60 files with no data, just a short header
netmind.db( DS="load", Y=p$netmind.yToload ) 


# add troublesome id's to the list here .. eventually move into their respective functions
p$netmensuration.problems = c() 

# now compute statistics (means, sd, etc)
seabird.db (DS="stats.redo", Y=p$seabird.yToload )
minilog.db (DS="stats.redo", Y=p$minilog.yToload )  # note no depth in minilog any more (since 2014 ..  useful for temperature only)
netmind.db (DS="stats.redo", Y=p$netmind.yToload )

```


#### Merging data and compute aggregate variables

Merge in netmind, minilog, seabird, esonar data andstore in set-level information ... and do some sanity checks.

NOTE: There are many historical data elements that require manual fixes to the raw data.


```r

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


# QA/QC of morphology (individuual level data)
# sanity check: identify morphology errors (they are written to logs), fix them if any are found and re-run until satisfied.. 
# if no morphology errors exist, you will get an error message. .
snowcrab.db( DS="det.initial.redo", p=p ) 
snowcrab.db( DS="det.georeferenced.redo", p=p )  # merge set.clean

# sanity check aggregate catches
snowcrab.db( DS="cat.initial.redo", p=p )
snowcrab.db( DS="cat.georeferenced.redo", p=p )  # merge set.clean

# aggregate catches by snow crab categories via det.initial, cat.initial
snowcrab.db( DS="set.biologicals.redo", p=p )

# update a database of simple transformation ranges, etc.. for plotting range, etc.
snowcrab.db( DS="data.transforms.redo", p=p) 

# create some simple/crude timeseries by each CFA using set.biologicals and observer data 
snowcrab.timeseries.db( DS="biologicals.redo", p=p )
snowcrab.timeseries.db( DS="observer.redo", p=p )

```

End of core snow crab data assimilation.

#### Size-frequency 

Size frequency distributions of snow crab carapace width from trawl data, broken down by maturity classes. This uses "set.clean" and "det.initial".

```r
xrange = c(10, 150)  # size range (CW)
dx = 2 # width of carapace with discretization to produce "cwd"
years = as.character( c(-9:0) + year.assessment )

M = size_distributions(p=p, toget="crude", xrange=xrange, dx=dx, Y=years, redo=TRUE)

```
#### Generate timeseries

Generate some simple/generic timeseries before entry into the main Markdown reports. We do them here instead as they need to be created only once:

```r

ts_years = 2004:p$year.assessment
ts_outdir = file.path( p$annual.results, "timeseries", "survey")
 
figure.timeseries.survey(p=p, outdir=ts_outdir, plotyears=ts_years) # all variables
 
# potential predators
species_predator = c(10, 11, 30, 40, 50, 201, 202, 204 )
bc_vars = c(paste("ms.mass", species_predator, sep='.'), paste("ms.no", species_predator, sep='.'))
figure.timeseries.survey(p=p, outdir=ts_outdir, plotyears=ts_years, variables=bc_vars )

# potential competitors
species_competitors = c( 2521, 2511, 2211)
bc_vars = c(paste("ms.mass", species_competitors, sep='.'), paste("ms.no", species_competitors, sep='.'))
figure.timeseries.survey(p=p, outdir=ts_outdir, plotyears=ts_years, variables=bc_vars )


over_ride_default_scaling = FALSE
if (over_ride_default_scaling) {

  figure.timeseries.survey(p=p, outdir=ts_outdir, plotyears=ts_years, variables="R0.mass" ) 
  
  figure.timeseries.survey(p=p, outdir=ts_outdir, plotyears=ts_years, variables=c("sexratio.all","sexratio.mat","sexratio.imm"))

  figure.timeseries.survey(p=p, outdir=ts_outdir, plotyears=ts_years, variables="cw.male.mat.mean", backtransform=TRUE) 

}

create_deprecated_figures = FALSE
if (create_deprecated_figures) {
    
    # no longer relevant (incomplete) as temp is now created through temp db. and not in gshyd
    figure.timeseries.survey(p=p, outdir=ts_outdir, plotyears=ts_years,type='groundfish.t') # groundfish survey temperature

    # area-specific figures
    figure_area_based_extraction_from_carstm(DS="temperature", year.assessment )  # tis can only do done once we have an sppoly for snow crab and temperature/carstm ha been completed
}

```


#### Generate maps

Need to generate some simple maps before entry into the main Markdown reports. We do them here instead as they are very slow (up to a few hours) and need to be created only once:

```r
 
map_outdir = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual" )
map_years  = p$year.assessment + c(0:-3)
 
create_maps_for_the_road_show = FALSE
if (create_maps_for_the_road_show) {

    road_show_vars = c("totmass.male.com", "totmass.female.mat", "R0.mass" )
    map.set.information( p=p, outdir=map_outdir, mapyears=map_years, variables=road_show_vars )

    map.set.information( p=p, outdir=map_outdir, mapyears=map_years, variables='t', log.variable=FALSE, add.zeros=FALSE, theta=100)
}

# variables that shouldn't be logged
set = snowcrab.db( p=p, DS="set.biologicals")
variables = bio.snowcrab::snowcrab.variablelist("all.data")
variables = intersect( variables, names(set) )

nolog.variables = c("t","z", "julian", variables[grep("cw",variables)])
map.set.information( p=p, outdir=map_outdir, mapyears=map_years, variables=nolog.variables, log.variable=FALSE, add.zeros=FALSE, theta=35)

# logit transform for ratios
ratio_vars = c("sexratio.all","sexratio.mat","sexratio.imm")
map.set.information( p=p, outdir=map_outdir, mapyears=map_years, variables=ratio_vars, log.variable=FALSE, add.zeros=FALSE, theta=40)

mass_vars = variables[!variables%in%nolog.variables][grep('mass',variables[!variables%in%nolog.variables])]
map.set.information( p=p, outdir=map_outdir, mapyears=map_years, variables= mass_vars )

no_vars = variables[!variables%in%nolog.variables][grep('no',variables[!variables%in%nolog.variables])]
map.set.information( p=p, outdir=map_outdir, mapyears=map_years, variables= no_vars,  probs=c(0,0.975))


# potential predators
species_predator = c(10, 11, 30, 40, 50, 201, 202, 204 )
bc_vars = c(paste("ms.mass", species_predator, sep='.'), paste("ms.no", species_predator, sep='.'))
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )
map.set.information( p=p, outdir=outdir_bc, mapyears=map_years, variables=bc_vars, probs=c(0,0.975)) 


# potential competitors
species_competitors = c( 2521, 2511, 2211)
bc_vars = c(paste("ms.mass", species_competitors, sep='.'), paste("ms.no", species_competitors, sep='.'))
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )
map.set.information( p=p, outdir=outdir_bc, mapyears=map_years, variables=bc_vars, probs=c(0,0.975)) 


# all variables (geometric means)
# map.set.information( p, outdir=map_outdir) # takes a long time ... run only vars of interest

```

#### Create simple Markdown reports (Quarto)

Basic data QA/QC, assimilation is complete.    

Most fishery related figures and tables can now be created in: 

- [02_fishery_summary.md](https://github.com/jae0/bio.snowcrab/tree/master/inst/markdown/02_fishery_summary.md) 


Basic trawl survey-related figures and tables can also be created in:

- [02_survey_summary.md](https://github.com/jae0/bio.snowcrab/tree/master/inst/markdown/02_survey_summary.md) 




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
# this step uses "set.clean" 

ps = snowcrab_parameters(
  project_class="carstm",
  yrs=1999:year.assessment,   
  areal_units_type="tesselation",
  carstm_model_label=  paste( "default", "fb", sep="_" )  # default for fb (fishable biomass)
)

xydata = snowcrab.db( p=ps, DS="areal_units_input", redo=TRUE )
# xydata = snowcrab.db( p=ps, DS="areal_units_input" )


# for mapping (background) .. redo=TRUE if resetting colours etc
additional_features = snowcrab_mapping_features(ps, redo=FALSE )  


# create constrained polygons with neighbourhood as an attribute
sppoly = areal_units( p=ps, xydata=xydata, spbuffer=3, n_iter_drop=0, redo=TRUE, verbose=TRUE )  # this needs to match carstm related parameters in snowcrab_parameters

# sppoly=areal_units( p=ps )

plot(sppoly["AUID"])

sppoly$dummyvar = ""
xydata = st_as_sf( xydata, coords=c("lon","lat") )
st_crs(xydata) = st_crs( projection_proj4string("lonlat_wgs84") )

require(tmap)
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

NOTE: this is a long process and requires a lot of RAM (aim for 128 GB RAM anda large swap space) ... 1-2 days. Not sure if it will run on tethys as it only has 64GB RAM.
  

#### Data from other surveys (groundfish)

Namely groundfish surveys provide supplemental data for species composition analysis:

  - [aegis.survey/inst/scripts/01_survey_data.md](https://github.com/jae0/aegis.survey/tree/master/inst/scripts/01_survey_data.md).


#### Species composition 

Species composition is a covariate in the areal model for snow crab. To obtain a predictive surface, the following needs to be run: 
    
  - [aegis.speciescomposition/inst/scripts/03_speciescomposition_carstm.md](https://github.com/jae0/aegis.speciescomposition/tree/master/inst/scripts/03_speciescomposition_carstm.md)


## Finalize data sets for modelling

We finalize the snow crab data after completing empirical data lookups of covariates (bathymetry and temperatures), in preparation for modelling.

```r

snowcrab.db( DS="set.complete.redo", p=p ) # note depth is log transformed here 


```



## Deprecated

Unused but stored in case.

```r
snowcrab.timeseries.db( DS="groundfish.t.redo", p=p )  # deprecated to be removed shortly
snowcrab.timeseries.db( DS="biologicals.2014.redo" )  # reduced subset that matches 2014 station id's .. # deprecated

# example: to get a quick view of a few vars of interest, region of interest ... no saving to file, but does return the data for viewing
snowcrab.timeseries.db( DS="biologicals.direct", p=p, regions='cfa4x', vn=c('R0.mass'), trim=0 )  
```

