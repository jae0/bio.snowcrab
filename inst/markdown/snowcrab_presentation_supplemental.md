---
title: "Snow Crab, Scotian Shelf, Canada (NAFO Div. 4VWX)"
subtitle: "Supplemental information"
metadata-files:
  - _metadata.yml
params:
  year_assessment: 2024
  year_start: 1999
  data_loc:  "~/bio.data/bio.snowcrab"
  sens: 1
  debugging: FALSE
  model_variation: logistic_discrete_historical
  todo: [fishery_results,fishery_model,ecosystem]
--- 


<!-- 

make quarto FN=snowcrab_presentation_supplemental.md DOCTYPE=revealjs  PARAMS="-P year_assessment:2024 -P todo:[fishery_results,fishery_model,ecosystem,redo_data]"  --directory=~/bio/bio.snowcrab/inst/markdown

-->

 

```{r}
#| eval: true
#| output: false
#| echo: false
#| label: setup

  require(knitr)

  knitr::opts_chunk$set(
    root.dir = data_root,
    echo = FALSE,
    out.width="6.2in",
    fig.retina = 2,
    dpi=192
  )
  
```


<!-- 
# _load_results.qmd contains instructions in "todo" 
#  this is a shared R-script to boot strap and provide a consistent data interface
-->

{{< include _load_results.qmd >}}  

 
 
## Climate

## Connectivity: Oceanic currents {.c}

  - Warm, high salinity Gulf Stream from the S-SE along the shelf edge 

  - Cold, low salinity Labrador Current

  - Cold low salinity St. Lawrence outflow from the N-NE

  - Nearshore Nova Scotia current, running from the NE. 
  

```{r}
#| label: fig-ocean-currents
#| eval: true
#| echo: false 
#| output: true
#| fig-cap: "Ocean currents in the Maritimes.  Source: [DFO](https://www.dfo-mpo.gc.ca/oceans/publications/soto-rceo/2018/atlantic-ecosystems-ecosystemes-atlantiques/index-eng.html)"
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold 

fn = file.path( media_loc, c(
  "maritimes_currents.png" 
) )
include_graphics( fn ) 

```
 
## SST

```{r}
#| label: fig-rapid-climate-change
#| eval: true
#| echo: false 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold 
#| fig-cap: "Global surface (2 meter) air temperature.  Source: [The Crisis Report](https://richardcrim.substack.com/p/the-crisis-report-99) and [James E. Hansen](https://www.columbia.edu/~jeh1/mailings/2024/ICJ.PressBriefing.09December2024.pdf). Note 2023 was and El Nino year."
#| fig-subcap: 
#|   - "Anomalies relative to pre-industrial baseline."
#|   - "Seasonal variations by year."
#|   - "Annual temperatures."

fn = file.path( media_loc, c(
  "gst_anomaly.png",
  "gst_seasonal.png",
  "gst_ts.png"
) )

include_graphics( fn ) 

```
 
## Chl-a
  
```{r}
#| label: fig-ocean-productivity-chla
#| eval: true
#| echo: false 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold 
#| fig-cap: "Chlorphyll-a in the NW Atlantic.  Source: [Copernicus Marine Service](https://marine.copernicus.eu/access-data/ocean-monitoring-indicators/chlorophyll-and-primary-production)"
#| fig-subcap: 
#|   - "Full timeseries of surface Chl-a estimated from satellite imagery."
#|   - "Trends in surface Chl-a over time."

fn = file.path( media_loc, c(
  "copernicus_chla.png",
  "copernicus_chla_map.png"
) )
include_graphics( fn ) 

```
 

 
## Bathymetry {.c}
::: columns 

:::: column 

```{r bathymetry-map, out.width='100%', echo=FALSE, fig.align='center', fig.cap='Variations in log(depth; m) in the Scotian Shelf region.' }
bathydir = file.path( data_root, 'aegis', 'bathymetry', 'modelled', 'default', 'stmv', 'none_fft', 'z', 'maps', 'SSE' )
include_graphics( file.path( bathydir, 'bathymetry.z.SSE.png' ) )
```
::::
:::: column
```{r bathymetry-map2, out.width='100%', echo=FALSE, fig.align='center', fig.cap='Local SD log(depth; m) in the Scotian Shelf region.' }
bathydir = file.path( data_root, 'aegis', 'bathymetry', 'modelled', 'default', 'stmv', 'none_fft', 'z', 'maps', 'SSE' )
include_graphics( file.path( bathydir, 'bathymetry.b.sdSpatial.SSE.png' ) )
```

::::
:::


## Substrate {.c}
::: columns 

:::: column 

```{r substrate-map, out.width='100%', echo=FALSE, fig.align='center', fig.cap='Substrate grain size log(mm) variations in the Scotian Shelf region.' }
substrdir = file.path( data_root, 'aegis', 'substrate', 'maps', 'canada.east.highres' )
include_graphics( file.path(  substrdir, 'substrate.substrate.grainsize.canada.east.highres.png' ) )
```
::::
:::: column
```{r substrate-map2, out.width='100%', echo=FALSE, fig.align='center', fig.cap='Local SD substrate grain size log(mm) in the Scotian Shelf region.' }
substrdir = file.path( data_root, 'aegis', 'substrate', 'maps', 'canada.east.highres' )
include_graphics( file.path( substrdir, 'substrate.s.sdSpatial.canada.east.highres.png' ) )
```

::::
:::


## Bottom Temperature {.c} 


```{r}
#| label: fig-bottom-temperatures-timeseries
#| eval: true
#| echo: false 
#| output: true
#| fig-cap: "Bottom temperatures"
#| fig-subcap: 
#|   - "Annual variations in bottom temperature observed during the Snow Crab survey. The horizontal (black) line indicates the long-term, median temperature within each subarea. Error bars represent standard errors."
#|   - "Posterior densities of predicted average bottom temperatures from an analysis of historical temperature data using [carstm](https://github.com/jae0/carstm). Red horizontal line is at $7^\\circ$C."
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold 

tloc = file.path( data_loc, "assessments", year_assessment, "timeseries"  )

fns = c( 
  file.path("survey", "t.png"), 
  "temperature_bottom.png" 
)

include_graphics( file.path( tloc, fns) )

``` 
 

```{r}
#| label: fig-figures-temperature-bottom-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: "Bottom temperature ($^\\circ$C) observed during the Snow Crab survey."
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""
 
map_outdir = file.path( p$project.outputdir, "maps", "survey", "snowcrab", "annual" )
map_years  = year_assessment + c(0:-3)
  
fn = check_file_exists( file.path( 
  map_outdir, "t", paste( "t", map_years, "png", sep="." )  
) )

include_graphics( fn )
```
 

 
 
```{r}
#| label: fig-figures-temperature-bottom-map-predicted
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: "Predicted bottom temperature ($^\\circ$C) for 1 September fromhstorical re-analysis."
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - "" 
#|   - ""
#|   - ""
#|   - ""
#|   - "" 
#|   - ""
#|   - ""
#|   - ""
#|   - "" 
#|   - ""
#|   - ""
#|   - "" 

loc = file.path( data_root, 'aegis', 'temperature', 'modelled', 'default', 'maps' )
yrsplot =  year_assessment + c(0:-14)

fns = file.path( loc, paste( 'predictions.',  yrsplot, '.0.75',  '.png', sep='') )

include_graphics( fns )
```

## Movement

 
```{r}
#| label: fig-movement-tracks
#| eval: true
#| echo: false 
#| output: true
#| fig-cap: "Snow Crab movement"
#| fig-subcap: 
#|   - "Tracks from 1996-2004"
#|   - "Tracks from 2004 - present"
#|   - "Distance between mark and recapture (km)"
#|   - "Minimum speed for each mark-recapture event (km/month)"
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| layout: [[100], [100], [50,50] ]
  
fns = file.path( media_loc, c(
  "movement0.png", 
  "movement.png" ,
  "snowcrab_movement_distances.png", 
  "snowcrab_movement_rates.png" 
) )

include_graphics( fns ) 

``` 

## Habitat preferences
 

```{r}
#| label: fig-temp-depth
#| eval: true
#| echo: false 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "Habitat preferences associated with depth and temperature."
 
fn2=file.path( media_loc, "viable_habitat_depth_temp.png" )
include_graphics( c( fn2 ) ) 

```
 

```{r}
#| label: fig-viable-habitat-persistent
#| eval: true
#| echo: false 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "Persistent habitat, independent of temperature, time, etc."

fn1 = file.path( media_loc, "viable_habitat.png" ) 
include_graphics( c(fn1  ) ) 

```

## Co-occurring species


## Atlantic cod

```{r}
#| label: fig-atlcod-timeseries
#| eval: true
#| output: true
#| fig-cap: "Mean density of Atlantic cod log$_{10}$(no/km$^2$) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

species_predator = 10

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```



```{r}
#| label: fig-atlcod-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: Atlantic cod, density; log$_{10}$(no/km$^2$). 
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""
 
map_years  = year_assessment + c(0:-3)
  
species_predator = 10
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```

## Haddock

```{r}
#| label: fig-haddock-timeseries
#| eval: true
#| output: true
#| fig-cap: "Mean density of Haddock log$_{10}$(no/km$^2$) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}
species_predator = 11

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```

 

```{r}
#| label: fig-haddock-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: Haddock, density; log$_{10}$(no/km$^2$). 
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""
 
map_years  = year_assessment + c(0:-3)
  
species_predator = 11
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```



## Halibut

  
```{r}
#| label: fig-halibut-timeseries
#| eval: true
#| output: true
#| fig-cap: "Mean density of Halibut log$_{10}$(no/km$^2$) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

species_predator = 30

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```

 

```{r}
#| label: fig-halibut-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: Halibut, density; log$_{10}$(no/km$^2$). 
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""
 
map_years  = year_assessment + c(0:-3)
  
species_predator = 30
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```

## American plaice
  
```{r}
#| label: fig-amerplaice-timeseries
#| eval: true
#| output: true
#| fig-cap: "Mean density of American plaice log$_{10}$(no/km$^2$) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

species_predator = 40

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```
 

```{r}
#| label: fig-amerplaice-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: American plaice, density; log$_{10}$(no/km$^2$). 
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

map_years  = year_assessment + c(0:-3)
  
species_predator = 40
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```
 
## Striped Alantic Wolffish

```{r}
#| label: fig-stripatlwolffish-timeseries
#| eval: true
#| output: true
#| fig-cap: "Mean density of Striped Atlantic wolffish log$_{10}$(no/km$^2$) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

species_predator = 50

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```

 
```{r}
#| label: fig-stripatlwolffish-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: Striped Atlantic wolffish, density; log$_{10}$(no/km$^2$). 
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

map_years  = year_assessment + c(0:-3)
  
species_predator = 50
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```
 


## Thorny skate

  
```{r}
#| label: fig-thornyskate-timeseries
#| eval: true
#| output: true
#| fig-cap: "Mean density of Thorny skate log$_{10}$(no/km$^2$) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

species_predator = 201

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```

 

```{r}
#| label: fig-thornyskate-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: Thorny skate, density; log$_{10}$(no/km$^2$). 
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

map_years  = year_assessment + c(0:-3)
  
species_predator = 201
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```
 


## Northern shrimp


```{r}
#| label: fig-northernshrimp-timeseries
#| eval: true
#| output: true
#| fig-cap: "Mean density of Northern shrimp log$_{10}$(no/km$^2$) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

species_predator = 2211

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```

 

```{r}
#| label: fig-northernshrimp-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: Northern shrimp, density; log$_{10}$(no/km$^2$). 
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

map_years  = year_assessment + c(0:-3)
  
species_predator = 2211
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```
 

## Jonah crab


```{r}
#| label: fig-jonahcrab-timeseries
#| eval: true
#| output: true
#| fig-cap: "Mean density of Jonah crab log$_{10}$(no/km$^2$) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

species_predator = 2511

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```





```{r}
#| label: fig-jonahcrab-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: Jonah crab, density; log$_{10}$(no/km$^2$). 
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

map_years  = year_assessment + c(0:-3)
  
species_predator = 2511
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```
 


## Lesser toad crab

```{r}
#| label: fig-lyrecrab-timeseries
#| eval: true
#| output: true
#| fig-cap: "Mean density of Arctic Lyre crab (Lesser toad crab) log$_{10}$(no/km$^2$) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

species_predator = 2521

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```





```{r}
#| label: fig-lyrecrab-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: Arctic Lyre crab, density; log$_{10}$(no/km$^2$). 
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

map_years  = year_assessment + c(0:-3)
  
species_predator = 2521
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```
 

## Northern stone crab

```{r}
#| label: fig-nstonecrab-timeseries
#| eval: true
#| output: true
#| fig-cap: "Mean density of Northern stone crab log$_{10}$(no/km$^2$) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

species_predator = 2523

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```


```{r}
#| label: fig-nstonecrab-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: Northern stone crab, density; log$_{10}$(no/km$^2$). 
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

map_years  = year_assessment + c(0:-3)
  
species_predator = 2523
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```
 


## Potentially interacting species ("bycatch")
 
```{r}
#| label: fig-competitor-biplot-nens
#| eval: true
#| echo: false 
#| output: true
#| fig.show: hold 
#| fig-dpi: 144
#| fig-height: 8
#| fig-cap: "Bycatch as potentially interacting species in N-ENS (1999-2020). Relative location of interacting species (green) in the species composition ordination. Snow crab in orange."
 
o = BC[["cfanorth"]]    
lookup = bio.taxonomy::taxonomy.recode( from="spec", to="taxa", tolookup=o$specid )$vern
xlab = paste("PC1 (", pca$variance_percent[1], "%)", sep="" )
ylab = paste("PC2 (", pca$variance_percent[2], "%)", sep="" )
plot( PC2 ~ PC1, pcadata, type="n", xlab=xlab, ylab=ylab )
text( PC2 ~ PC1, labels=vern, data=pcadata, cex=0.75, col="slategrey"  )
i = grep("Snow crab", pcadata$vern, ignore.case=TRUE)
points( PC2 ~ PC1, pcadata[i,], pch=19, cex=3.0, col="darkorange" )
j = NULL
for (k in lookup) j = c(j, grep( k, pcadata$vern, ignore.case=TRUE))    
j = unique(j)
points( PC2 ~ PC1, pcadata[j,], pch=19, cex=2.0, col="lightgreen" )
text( PC2 ~ PC1, labels=vern, data=pcadata[j,], cex=0.75, col="darkgreen"  )
```

```{r}
#| label: fig-competitor-biplot-sens
#| eval: true
#| echo: false 
#| output: true
#| fig.show: hold 
#| fig-dpi: 144
#| fig-height: 8
#| fig-cap: "Bycatch as potentially interacting species in S-ENS (1999-2020). Relative location of interacting species (green) in the species composition ordination. Snow crab in orange."
 
o = BC[["cfasouth"]]
lookup = bio.taxonomy::taxonomy.recode( from="spec", to="taxa", tolookup=o$specid )$vern
xlab = paste("PC1 (", pca$variance_percent[1], "%)", sep="" )
ylab = paste("PC2 (", pca$variance_percent[2], "%)", sep="" )
plot( PC2 ~ PC1, pcadata, type="n", xlab=xlab, ylab=ylab )
text( PC2 ~ PC1, labels=vern, data=pcadata, cex=0.75, col="slategrey"  )
i = grep("Snow crab", pcadata$vern, ignore.case=TRUE)
points( PC2 ~ PC1, pcadata[i,], pch=19, cex=3.0, col="darkorange" )
j = NULL
for (k in lookup) j = c(j, grep( k, pcadata$vern, ignore.case=TRUE))    
j = unique(j)
points( PC2 ~ PC1, pcadata[j,], pch=19, cex=2.0, col="lightgreen" )
text( PC2 ~ PC1, labels=vern, data=pcadata[j,], cex=0.75, col="darkgreen"  )
```
  

```{r}
#| label: fig-competitor-biplot-4x
#| eval: true
#| echo: false 
#| output: true
#| fig.show: hold 
#| fig-dpi: 144
#| fig-height: 8
#| fig-cap: "Bycatch as potentially interacting species in 4X (1999-2020). Relative location of interacting species (green) in the species composition ordination. Snow crab in orange."
 
o = BC[["cfa4x"]]
lookup = bio.taxonomy::taxonomy.recode( from="spec", to="taxa", tolookup=o$specid )$vern
xlab = paste("PC1 (", pca$variance_percent[1], "%)", sep="" )
ylab = paste("PC2 (", pca$variance_percent[2], "%)", sep="" )
plot( PC2 ~ PC1, pcadata, type="n", xlab=xlab, ylab=ylab )
text( PC2 ~ PC1, labels=vern, data=pcadata, cex=0.75, col="slategrey"  )
i = grep("Snow crab", pcadata$vern, ignore.case=TRUE)
points( PC2 ~ PC1, pcadata[i,], pch=19, cex=3.0, col="darkorange" )
j = NULL
for (k in lookup) j = c(j, grep( k, pcadata$vern, ignore.case=TRUE))    
j = unique(j)
points( PC2 ~ PC1, pcadata[j,], pch=19, cex=2.0, col="lightgreen" )
text( PC2 ~ PC1, labels=vern, data=pcadata[j,], cex=0.75, col="darkgreen"  )
```
  
  
## Human interactions: Every known population is exploited worldwide {.c}

- Human consumption (males; > 95mm CW) due to higher meat yield and sexual dimorphism
- Bait in other fisheries (illegal in Canada, but still present)
- Fertilizer in agriculture
- Chitosan: glucosamine polysaccharide derived from chitin
    - agent to slow bleeding from wounds (Zhang et al. 2015, Moghadas et al. 2016)
    - agriculturally as natural fungicides and bactericides (Linden & Stoner 2007)
    - plastic replacement (Tampieri et al. 2003) 
    - battery electrolyte (Poosapati et al. 2021).

 
 

## Life history {.c}

```{r}
#| label: fig-photos-snowcrab-1
#| eval: true
#| echo: false 
#| output: true
#| fig-cap: "Snow Crab male, mature benthic form."
#| fig-subcap: 
#|   - ""
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| layout-ncol: 1

fns = file.path( media_loc, "snowcrab_male.png" )

include_graphics( fns ) 

```

 

## {#Photographs data-menu-title="Photographs"} 

```{r}
#| label: fig-photos-snowcrab-2
#| eval: true
#| echo: false 
#| output: true
#| fig-cap: "Snow Crab images. Note sexual dimorphism."
#| fig-subcap: 
#|   - "Pelagic zoea"
#|   - "Mating pair - note sexual dimorphism, with a smaller female. Note epibiont growth on female which suggests it is an older female (multiparous)."
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| layout-ncol: 2

fns = file.path( media_loc, c(
  "snowcrab_zoea.png", 
  "snowcrab_male_and_female.png"
) )

include_graphics( fns ) 

```

 
 


## Life history stages
 
```{r}
#| label: fig-lifehistory
#| eval: true
#| echo: false 
#| output: true
#| fig-cap: "Life history patterns of snow crab and approximate timing of the main life history stages of snow crab and size (carapace width; CW mm) and instar (Roman numerals). Size and timings are specific to the area of study and vary with environmental conditions, food availability and genetic variability."
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold

fn1=file.path( media_loc, "life_history.png" )
include_graphics( fn1 ) 

```

## Male growth stanzas 
 
```{r}
#| label: fig-growth-stanzas
#| eval: true
#| echo: false 
#| output: true
#| fig-cap: "The growth stanzas of the male component and decision paths to maturity and terminal moult. Black ellipses indicate terminally molted animals."
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
 
fn1=file.path( media_loc, "life_history_male.png" )
include_graphics( fn1 ) 

```

 
## Growth patterns inferred from size modes

```{r}
#| label: fig-growth-modes
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: "Modes (CW, ln mm) identified from survey data using *Kernel Density Estmation* of local moving data windows. Legend: sex|maturity|instar"
#| fig-subcap: 
#|   - "Female modes"
#|   - "Male modes" 


fns = file.path( media_loc, c(
  "density_f_imodes.png", 
  "density_m_imodes.png"
))

include_graphics( fns )
```
 
  
## Growth patterns inferred from modes
 
```{r}
#| label: fig-growth-modes-growth
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: "Inferred growth derived from *Kernel Mixture Models* (priors)."
#| fig-subcap: 
#|   - "Female growth trajectory"
#|   - "Male growth trajectory"

fns = file.path( media_loc, c(
  "plot_growth_female.png",
  "plot_growth_male.png"
))

include_graphics( fns )
```
 
 
 



## Life history: notable traits

  - Sexual dimorphism
  - Pelagic in larval stages; benthic in pre-adolescent and adult stages
  - Biannual, annual molts depending upon size/age and environmental conditions
  - Terminal molt to maturity and survives for up to 5 years
  - Total life span: up to 15 years. 
  - Stenothermic (narrow) temperature requirements < 6 C
  - Ontogenetic shifts in habitat preferences: warmer, complex substrates to mud.  
  - Cannibalism of immature crab by mature female snow crab are known 
  - Primiparous female (57.4 mm CW) produces between 35,000 to 46,000 eggs
  - Multiparous females more fecund (>100,000 eggs) 
  - Eggs extruded February and April and brooded for up to two years (> 80% in SSE follow an annual cycle) and hatched/released from April to June.
   

## Life history: Spatially and temporally complex
  - 14+ instars/life stages with different preferences; 
  - Multiply prey and predators (ontogenetic shift); many ecosystem changes 
  - Growth: biannual (up to instar 5), annual, biennial (skip moulters)
  - Reproduction: annual and biennial; strong sexual selection (dimorphism) and parental investment for up to 2 years; lags of up to 10 years
  - Mortality: pulsed due to moult vulnerability, cannibalism
  - Movement: seasonal, annual, ontogenetic 
  - Evolution: sexual dimorphism, parental investment of brood, stenothermic
  - Population structure integrates all these effects ("Cumulative Effects")




## Population structure integrates everything
  - Spatial structures: scales > 100's km; local crab holes ( 10's to 100's m ) 
    - Confluence of major oceanic currents (temperature, salinity, nutrients, primary and secondary production)
  - Temporal processes: scales > 10 yrs; daily, tidal currents,  
    - 10-15 y life cycle, seasonal, annual, multiyear, decadal, evolutionary, etc. 
  - Biological (individual, age classes, population, etc.): cross multiple space, time and organismal scales 
    - Ontogeny (developmental change and habitat requirements)
    - Life history (habitat preferences)
    - Inter-specific interactions (feeding relationships change with biological stage)
    - Interactions with humans

## Multiple processes at multiple scales 

  - Populations integrate these processes  
  - **Induces** auto-correlations at multiple spatial-temporal-organisational scales
  - We perceive this as "autocorrelation"
  - First law of geography: "everything is related to everything else, but near things are more related than distant things" (Tobler 1970) 
  - "Nonstationary" processes (mean and variance not stable in space and time)
  - Ignoring correlated errors causes inappropriate parameter estimates


## Observations error (precision, accuracy) 

  - Pragmatic compromise between costs vs information gain
  - SSE snow crab domain with ~400 stations and ~109,120  $\text{km}^{2}$ 
  - Each station represents 273 $\text{km}^{2}$ , but actually samples ~0.0039 $\text{km}^{2}$ 
    - A factor of 1:70,000 (spatially)
  - Each sample is about 5 minutes but represents a whole year (525,600 minutes)
    - A factor of 1:105,120
  - In space-time: 
    - A factor of ~ 1:10 billion (note: 1 Angstrom = 10$^{-10}$m) .. the same scale as an atom! 
    - A few samples will do if it is spatially homogenous ("well mixed" = constant forcings, no ecosystem variability, no spatiotemporal autocorrelation)
  - Spatial bias/aliasing
  - Temporal bias/aliasing ("Spring" before 2004)
  

## Solutions: Experimental design

  - Naive (Ignore everything) "Random" Sampling / Fixed Station Sampling (IID)
  - Stratified (Hope this works) "Random" Sampling -- depth is all that matters (IID within, IID between)
    - Ignore spatiotemporal structure 
    - Ignore Ecosystem variability
    - Assume **no autocorrelated errors**
  - Geostatistical design (gridded, pseudo-random)
    - Address spatial autocorrelation as a stationary process
    - Abundance is NOT first and second order stationary and non-Gaussian (distributional bias)
    - Ignores time (temporal discretization/aliasing)
    - Ignores ecosystem variability
    - Variogram solutions unstable when population size declines

## Solutions: Experimental design 2

  - Hierarchical process temporal and covariate as an external forcing and spatial as a local process (UKED)
    - Variogram instability across time/stage when abundance declined and distribution was spotty
    - NOT second-order stationary
  - Mosaic process:
    - Hierarchical first and second-order **non-stationary** model (modular)
    - Global Hurdle/INLA/GAM/GLM 
    - Local Random spatial effects (residuals) Matern (SPDE; FFT)
    - Independent space and time process
    - Slow
  - Model-based inference: Bayesian hierarchical mixed-effects model
    - Conditional Autoregressive Space-Time Models (CARSTM, separable: Kronecker product space BYM2, time AR1)
    - Hurdle/Poisson and Binomial number process, with Gaussian weight process


## Latent ecological process

- Real (latent, unobserved) ecological processes are generaly spatiotemporal processes  

- Observations from samples/surveys are used to infer the real (latent, unobserved) state

- Stock assessments, almost always, focus only upon a purely *temporal process* by integrating the latent spatiotemporal process  
 
- Experimental design assumes samples are random or random-stratified, almost always, as a purely *spatial process*

- This is a problem: 
  - temporal aliasing is likely (as it is usually ignored)
  - spatial aliasing is likely when sampling environments that cannot be randomly sampled 
    - many locations are not directly accessible to trawls (rocks and bedrock)
    - easier to sample locations (softer, gravel or mud substrates) coincide with preferred snow habitats


## Snow crab survey locations {.c}
   
```{r survey-locations-map, out.width='60%', fig.show='hold', fig.align='center', fig.cap= 'Snow Crab survey locations. No survey in 2020 (Covid-19) and incomplete 2022 (mechanical issues).' }
loc = file.path( data_loc, "output", "maps", "survey.locations" )
yrsplot = setdiff( year_assessment + c(0:-9), 2020)
fn6 = file.path( loc, paste( "survey.locations", yrsplot[6], "png", sep=".") )
fn5 = file.path( loc, paste( "survey.locations", yrsplot[5], "png", sep=".") )
fn4 = file.path( loc, paste( "survey.locations", yrsplot[4], "png", sep=".") )
fn3 = file.path( loc, paste( "survey.locations", yrsplot[3], "png", sep=".") )
fn2 = file.path( loc, paste( "survey.locations", yrsplot[2], "png", sep=".") )
fn1 = file.path( loc, paste( "survey.locations", yrsplot[1], "png", sep=".") )
include_graphics( c( fn2, fn1) )
# \@ref(fig:survey-locations-map)  
``` 

 
## {}

```{r surveydomain, out.width='100%', fig.show='hold', fig.align='center', fig.cap= 'Snow Crab survey prediction grid.' }
fn1 = file.path( media_loc, "carstm_prediction_domain.png" )
include_graphics( c( fn1) )
``` 

 
## Spatial clustering  
 
```{r}
#| label: fig-aggregation
#| eval: true
#| echo: false 
#| output: true
#| fig-cap: "Spider crab tend to cluster/aggregate in space."
#| fig-subcap: 
#|   - "Australian spider crab \\emph{Leptomithrax gaimardii} aggregation for moulting and migration."
#|   - "Alaska red king crab \\emph{Paralithodes camtschaticus} aggregation in Alaska for egg release, migrations."
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
 
fns = file.path( media_loc, c( 
  "australian_leptomithrax_gaimardii.png", 
  "kingcrab_aggregation.png" 
) )

include_graphics( fns ) 

```

- *Other* crab species show "Mounding" for protection from predation (larval, moulting and females).

- Narrow habitat preferences force them to move and cluster when environment is poor.

## Spatial clustering 2 
 

```{r}
#| label: fig-aggregation2
#| eval: true
#| echo: false 
#| output: true
#| fig-cap: "High density locations of Snow Crab, approximately 1 per square meter."
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
 
fn = file.path(p$project.outputdir, "maps", "map_highdensity_locations.png" )

include_graphics( fn ) 

```





## Sampling bias: depth
```{r bias-depth, out.width='40%', fig.show='hold', fig.align='center', fig.cap= 'Comparison of depths in snow crab domain ("predictions") and survey locations ("observations"). Bias: survey trawls preferentially sample deeper locations. Red line is overall average.' }

include_graphics( file.path( media_loc, "bias_depth.png" ) )
```
 
## Sampling bias: substrate grain size

```{r bias-substrate, out.width='40%', fig.show='hold', fig.align='center', fig.cap= 'Comparison of substrate grain size (log; mm) in snow crab domain ("predictions") and survey locations ("observations"). No bias observed. Red line is overall average.' }

include_graphics( file.path( media_loc, "bias_substrate.png" ) ) 
``` 

## Sampling bias: bottom temperature

```{r bias-temperature, out.width='40%', fig.show='hold', fig.align='center', fig.cap= 'Comparison of bottom temperature (degrees Celcius) in snow crab domain ("predictions") and survey locations ("observations"). Bias: surveys preferentially sample cold water bottoms. Red line is overall average.' }

include_graphics( file.path( media_loc, "bias_temp.png" ) ) 
``` 
 
## Sampling bias: species composition 1 (temperature related)

```{r bias-pc1, out.width='40%', fig.show='hold', fig.align='center', fig.cap= 'Comparison of species composition gradient (PC1) in snow crab domain ("predictions") and survey locations ("observations"). Bias: surveys preferentially sample cold water species. Red line is overall average.' } 

include_graphics( file.path( media_loc, "bias_pca1.png" ) ) 
``` 
  

## Sampling bias: species composition 2 (depth related)

```{r bias-pc2, out.width='40%', fig.show='hold', fig.align='center', fig.cap= 'Comparison of species composition gradient (PC2) in snow crab domain ("predictions") and survey locations ("observations"). Bias: surveys preferentially sample deeper species. Red line is overall average.' }

include_graphics( file.path( media_loc, "bias_pca2.png" ) ) 
``` 

## Sampling bias: summary

Indications of sampling bias relative to spatial domain of snow crab:

- deeper and more favourable locations are being sampled
- statistically overlapping but biologically important
 

## Generalized linear models (GLM)


$S={S_{1},...,S_{K}}$ is a set of $k=1,...,K$ non-overlapping (areal) units.

$\boldsymbol{Y}=(y_{1},...,y_{K})$ are observations on $S$, then:


$$\begin{aligned}
Y & \sim f(y|\Omega)
g(\mu) &=\boldsymbol{x}^{T}\boldsymbol{\beta}+\boldsymbol{O}+\boldsymbol{\varepsilon}.
\end{aligned}$$

- $\mu = \text{E}(Y)$  is the expected value of $\boldsymbol{Y}$

- $f(\cdot)$ indicates an exponential function  

- $\Omega$ is the set of the parameters of the function $f(\cdot)$ 

- $g(\cdot) = f^{-1}(\cdot)$ is the linearizing link function 
 
- $\boldsymbol{x=}(x_{kv})\in\Re^{K\times V}$ is the matrix of covariates

- $\boldsymbol{\beta}$ are the $V$ covariate parameters with MVN prior, with mean $\mu_{\beta}$ and diagonal variance matrix $\Sigma_{\beta}$
  
- $\boldsymbol{O=}(o_{1},...,o_{K}\boldsymbol{)}$ are offsets, if any
 
- $\boldsymbol{\varepsilon}=(\varepsilon_{1},...,\varepsilon_{K})$ are residual errors, if any

 
For each distributional family:
 
$Y\sim\text{Normal}(\mu,\sigma^{2})$ 

  - $\mu=\boldsymbol{x}^{T}\boldsymbol{\beta}+\boldsymbol{O}+\boldsymbol{\varepsilon}$

$Y\sim\text{Poisson}(\mu)$, and

  - $\text{ln}(\mu)=\boldsymbol{x}^{T}\boldsymbol{\beta}+\boldsymbol{O}+\boldsymbol{\varepsilon}$.

$Y\sim\text{Binomial}(\eta,\theta)$ 

  - $\text{ln}(\theta/(1-\theta))=\boldsymbol{x}^{T}\boldsymbol{\beta}+\boldsymbol{O}+\boldsymbol{\varepsilon}$,
  - $\eta$ is the vector of number of trials 
  - $\theta$ the vector of probabilities of success in each trial 
 

## **Stratified random sampling** assumes: 

$$\varepsilon_{s}  \sim\text{N}(0,~^{\varepsilon}\sigma_{s}^{2})$$

i.e., IID errors in space (time is usually ignored)

  
## Autocorrelated spatial errors 
 

**Kriging** extends $\varepsilon$ to spatial contraints ("variograms"). 

$$\begin{aligned}
Y_{t} & \sim\text{MVN}(\boldsymbol{\mu}_{t},\boldsymbol{\Sigma}_{t})\\
g(\mu_{t}) & =\boldsymbol{x_{t}}^{T}\boldsymbol{\beta_{t}}+\boldsymbol{\omega}_{t}+\boldsymbol{\varepsilon}_{t}\\
\varepsilon_{t} & \sim N(0,{}^{\varepsilon}\!\sigma_{t}^{2})\\
\omega_{t} & \sim\text{GP}(\boldsymbol{0},C(s_{t},s_{t}';^{\omega}\!\theta_{t}))\\
\boldsymbol{\Sigma_{t}} & =\left[C(\text{s}_{it},\text{s}_{jt};^{\omega}\!\theta_{t})\right]_{i,j=1}^{K}+{}^{\varepsilon}\!\sigma_{t}^{2}I_{S} \\
C(h)_{\text{Matérn}} & = ~^{\omega}\sigma^{2}\frac{1}{2^{\nu-1}\Gamma(\nu)}(\sqrt{2\nu}h/\phi)^{\nu}\ K_{\nu}(\sqrt{2\nu}h/\phi).
\end{aligned}$$
 
-   ignores temporal structure, focus upon spatial-only process

-   assumes first and second order stationarity (incorrect)

-   problems when abundance declines and spatial autocorrelation changes
    across time or becomes unstable or impossible to estimable


## Spatial autocorrelation
 
Spatial autocorrelation function describes the proportion
of the spatial variance $C(h=0)= ~^{\omega}\sigma^{2}$ decrease as distance
increases $h$.

Spatial covariance function $C(h)$ scaled by the total variance $C(0)$:

$\rho(h)=C(h)/C(0)$. 

The spatial covariance function:

$C(h)=C(s,s';\theta)$ 

expresses the tendency of observations closer together to be more similar to each other than those further away (Tobler's First law). 

Here, the distance, 

$h=\parallel s-s'\parallel$,

where $\parallel\cdot\parallel$ indicates a spatial norm which in $d=2$ spatial dimensions is simply the Euclidean distance:

$h=(\Delta\text{northing}^{2}+\Delta\text{easting}^{2})^{1/2}$. 



## Commonly used forms include:

$$\begin{matrix}C(h)_{\text{Spherical}} & = & \left\{ \begin{matrix}^{\omega}\sigma^{2}(1-\frac{3}{2}h/\phi+\frac{1}{2}(h/\phi)^{3}); & 0<h<=\phi\\
0; & h>\phi,
\end{matrix}\right.\\
C(h)_{\text{Exponential}} & = & ^{\omega}\sigma^{2}e^{-h/\phi},\\
C(h)_{\text{Gaussian}} & = & ^{\omega}\sigma^{2}e^{-(h/\phi)^{2}},\\
C(h)_{\text{Powered exponential}} & = & ^{\omega}\sigma^{2}e^{-|h/\phi|^{p}},\\
C(h)_{\text{Mat\'{e}rn}} & = & ^{\omega}\sigma^{2}\frac{1}{2^{\nu-1}\Gamma(\nu)}(\sqrt{2\nu}h/\phi)^{\nu}\ K_{\nu}(\sqrt{2\nu}h/\phi).
\end{matrix}$$

 
At zero distance,

$C(0)=\text{Cov}(Y_{s},Y_{s})=\text{Var}(Y_{s})={}^{\varepsilon}\sigma^{2}+^{\omega}\sigma^{2}$

(*i.e.*, global spatial variance, also called the sill), where

$^{\varepsilon}\sigma$ is the non-spatial, unstructured error (also
called the nugget), 

$^{\omega}\sigma$ is the spatially structured error (also called the partial sill), and

$^{\omega}\theta=\{\phi,\nu,p,\ldots\}$ are function-specific parameters

including $\phi$ the *range* parameter. 

$\Gamma(\cdot)$ is the Gamma function and 

$K_{\nu}(\cdot)$ is the Bessel function of the second kind with smoothness $\nu$. 

The Matérn covariance function is frequently used in the more recent literature as the shape of this function is more
flexible and known to be connected to a Gaussian spatial random process. The *stmv* approach  permits a local
estimation of these autocorrelation parameters and the spatial scale.
 

  

## STMV
 
```{r stmv-concept, out.width='40%', fig.show='hold', fig.align='center', fig.cap= 'Spatial distribution of data (blue dots) overlaid by a statistical grid. The $m$ nodes represent the centers of each local subdomain $S_{m}$ which extends to a distance (right-facing arrows; solid squares) that varies depending upon the underlying spatial variability of the data, defined as the distance at which the spatial autocorrelation drops to some small value (e.g., $\rho_{s}=0.1$). Data within this distance and parameters obtained from the local analysis are, under the assumption of second order stationarity' } 

include_graphics( file.path( media_loc, "stmv_concept.png" ) ) 
```
 


## Random field-based approximations using SPDE (INLA)
 
$$\begin{matrix}
Y & \sim & \text{N}(\mu_{st},^{\varphi}\!\sigma^{2})\\
g(\mu_{st}) & = & \boldsymbol{x}^{T}\boldsymbol{\beta}_{st}+\varphi_{st},\\
\varphi_{st} & \sim & \text{N}(0,^{\varphi}\!\sigma^{2}).
\end{matrix}$$

$$\begin{matrix}
\varphi_{mt}^{*} & = & \omega_{mt}+\varepsilon_{mt},\\
\omega_{mt} & \sim & \text{GP}(0,C(\mathbf{s},\mathbf{s}';\mathbf{^{\boldsymbol{\omega}}\theta}_{mt}=\{\nu_{mt},\phi_{mt},^{\omega}\sigma_{mt}\})),\\
\varepsilon_{mt} & \sim & \text{Normal}(0,{}^{\varepsilon}\sigma_{m}^{2}).
\end{matrix}$$
 
SPDE representation of a random field is efficiently computed from:

$$\frac{\partial}{\partial t} (\kappa(s)^{2}-\Delta)^{\alpha/2}(\tau(\mathbf{s})x(\mathbf{s},\mathbf{t}))=\mathcal{W}(\mathbf{s},\mathbf{t}),\;(\mathbf{s},\mathbf{t})\in\Omega\times\mathbb{R}$$
 

## autocorrelation

```{r autocorrelation, out.width='40%', fig.show='hold', fig.align='center', fig.cap= '' }

include_graphics( file.path( media_loc, "autocorrelation.png" ) ) 
``` 
\tiny

Matérn autocorrelation function, $\rho(h)=C(h)/C(0)$, the covariance function $C(h)$ scaled by the total variance $C(0)$, for two values of $\nu$ (dark lines). As $\nu$ increases $(\nu=100)$, it approaches the Gaussian curve (upper dark curve on the left side) while at smaller values $(\nu=0.5)$ the curve is exponential (lower dark curve on the left side). This flexibility has made it a popular choice in geostatistics. The associated semivariograms (scaled to unit variance) $\gamma(h)$ are shown in light stippled lines. Spatial scale is defined heuristically as the distance $h$ at which the autocorrelation falls to a low value (dashed horizontal line). The semivariance (also called a semivariogram), $\gamma(h)$, is more commonly used in the geostatistics literature, and is simply the covariance function $C(h)$ reflected on the horizontal axis of the global variance $C(0)$ such that $\gamma(h)=C(0)-C(h)=\frac{1}{2}\ \text{Var}[Y_{s}-Y_{s}']=^{\omega}\sigma^{2}[1-\rho(h)]$.

\normalsize

 
## CARSTM: Conditional Autoregressive Space-Time Models (via INLA)

\small

$$\begin{aligned}
Y &\sim f(\mu,{}^{\varepsilon}\!\sigma^{2})\\
g(\mu) &=\boldsymbol{x}^{T}\boldsymbol{\beta}+\boldsymbol{O}+\phi+\boldsymbol{\varepsilon}\\
\varepsilon &\sim\text{N}(0,^{\varepsilon}\!\sigma^{2})\\
\phi &\sim N(\boldsymbol{0},Q^{-1})\\
Q &=[\tau(D-\alpha W)]^{-1}\\
p(\phi_{i}|\phi_{i\sim j},\tau) &=N(\alpha\sum_{i\sim j}^{K}w_{ij}\phi_{j},\tau^{-1})
\end{aligned}$$

-   Bayesian Hierarchical models of time \| space (INLA)

-   Spatial autocorrelation modelled as a local autocorrelation
    (CAR/BYM)

-   Temporal autocorrelation modelled as a local AR1 (ARIMA)

-   Ecosystem variability enters as covariates (smooths)

\normalsize


## CARSTM: Temperature Posterior Predictive Check

```{r temp_ppc, out.width='50%', fig.show='hold', fig.align='center', fig.cap= '' }
loc = file.path( data_root, "aegis", "temperature", "modelled", "default" )
include_graphics( file.path( loc, "posterior_predictive_check.png" ) ) 
``` 

## CARSTM: Species composition Posterior Predictive Check

```{r pca_ppc, out.width='40%', fig.show='hold', fig.align='center', fig.cap= '' }
loc = file.path( data_root, "aegis", "speciescomposition", "modelled", "default" )
fn1 = file.path( loc, "pca1_posterior_predictive_check.png" )
fn2 = file.path( loc, "pca2_posterior_predictive_check.png" )
fn3 = file.path( loc, "pca3_posterior_predictive_check.png" )
include_graphics( c(fn1, fn2 ) ) 
``` 
 
## CARSTM: Snow crab Posterior Predictive Check

```{r snowcrab_ppc, out.width='32%', fig.show='hold', fig.align='center', fig.cap= '' }
loc = file.path( data_root, "bio.snowcrab", "modelled", "default_fb" )
fn1 = file.path( loc, "totno_posterior_predictive_check.png" )
fn2 = file.path( loc, "meansize_posterior_predictive_check.png" )
fn3 = file.path( loc, "pa_posterior_predictive_check.png" )
include_graphics( c(fn1, fn2, fn3) ) 
``` 
 


 

  

## Fishery model



$$\begin{aligned}
b_{t+1} & \sim N\left({{{b_{t}+r}b_{t}}{{({1-b_{t}})}-\mathit{Fishing}}_{t}K^{-1},\sigma_{p}}\right) \\
Y_{t} & \sim N\left({q^{-1}\mathit{Kb}_{t},\sigma_{o}}\right) \\
r & \sim N\left({1.0,0.1}\right) \\
K & \sim N\left({\kappa,0.1\cdot\kappa}\right) \\
q & \sim N\left({1.0,0.1}\right) \\
\end{aligned}$$

$\kappa={\lbrack{5.0,60,1.25}\rbrack}$

${b=B}K^{-1}$,  a real unobservable (latent) process


 
## Surplus (Shaeffer) production

```{r}
#| label: fig-logistic-surplus-production
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig-cap: "Surplus (Shaeffer) production. Each year is represented by a different colour."
#| fig-subcap:
#|   - "N-ENS"
#|   - "S-ENS"
#|   - "4X"

fns = paste("plot_surplus_production_", regions, ".png", sep="" )
fns = file.path( fm_loc, fns ) 

include_graphics( fns )
   
``` 


## Carrying capacity 

```{r}
#| label: fig-logistic-prior-posterior-K
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig-cap: "Prior-posterior comparisons of carrying capacity (K; kt)."
#| fig-subcap:
#|   - "N-ENS"
#|   - "S-ENS"
#|   - "4X"
   
fns = file.path( fm_loc, paste("plot_prior_K_", regions, ".png", sep="" ) ) 

include_graphics( fns )
   

``` 

## Intrinsic rate of increase

```{r}
#| label: fig-logistic-prior-posterior-r
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig-cap: "Prior-posterior comparisons of th iIntrinsic rate of biomass increase (r)."
#| fig-subcap:
#|   - "N-ENS"
#|   - "S-ENS"
#|   - "4X"

  
fns = file.path( fm_loc, paste("plot_prior_r_", regions, ".png", sep="" ) ) 

include_graphics( fns )

``` 

## Catchability coefficient
 
```{r}
#| label: fig-logistic-prior-posterior-q
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig-cap: "Prior-posterior comparisons of the catchability coefficient (q)."
#| fig-subcap:
#|   - "N-ENS"
#|   - "S-ENS"
#|   - "4X"

  
fns = file.path( fm_loc, paste("plot_prior_q1_", regions, ".png", sep="" ) ) 

include_graphics( fns )

``` 


## Observation error

```{r}
#| label: fig-logistic-prior-posterior-obserror
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig-cap: "Prior-posterior comparisons of Observation error (kt)."
#| fig-subcap:
#|   - "N-ENS"
#|   - "S-ENS"
#|   - "4X"

fns = file.path( fm_loc, paste("plot_prior_bosd_", regions, ".png", sep="" ) ) 

include_graphics( fns )

``` 
 
## Model process error

```{r}
#| label: fig-logistic-prior-posterior-processerror
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig-cap: "Prior-posterior comparisons of Model process error (kt)."
#| fig-subcap:
#|   - "N-ENS"
#|   - "S-ENS"
#|   - "4X"
 
fns = file.path( fm_loc, paste("plot_prior_bpsd_", regions, ".png", sep="" ) ) 

include_graphics( fns )

``` 

## State space 
  
```{r}
#| label: fig-logistic-state-space
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig-cap: "State space (kt): year vs year+1."
#| fig-subcap:
#|   - "N-ENS"
#|   - "S-ENS"
#|   - "4X"

fns = file.path( fm_loc, paste("plot_state_space_", regions, ".png", sep="" ) ) 

include_graphics( fns )

``` 
 
      
 
## References and further readings
 

Banerjee, S., Carlin, B. P., and Gelfand, A. E.. 2004. Hierarchical Modeling and Analysis for Spatial Data. Monographs on Statistics and Applied Probability. Chapman and Hall/CRC.

Besag, Julian. 1974. Spatial interaction and the statistical analysis of lattice systems. Journal of the Royal Statistical Society Series B (Methodological) 1974: 192-236.

Canada Gazette. 2022. [Regulations Amending the Fishery (General) Regulations. Part II, Volume 156, Number 8.](https://www.gazette.gc.ca/rp-pr/p2/2022/2022-04-13/html/sor-dors73-eng.html)

Canada Gazette. 2016. St. Anns Bank Marine Protected Area Regulations. Canada Gazette, Part I, Vol 150, Issue 51: 4143-4149.

Choi, J.S. 2020. A Framework for the assessment of Snow Crab (*Chioneocete opilio*) in Maritimes Region (NAFO Div 4VWX) . DFO Can. Sci. Advis. Sec. Res. Doc. 2020/nnn. v + xxx p.

Choi, J.S. 2022. Reconstructing the Decline of Atlantic Cod with the Help of Environmental Variability in the Scotian Shelf of Canada. bioRxiv. https://doi.org/10.1101/2022.05.05.490753.

Choi, J. S., and B. C. Patten. 2001. Sustainable Development: Lessons from the Paradox of Enrichment. Ecosystem Health 7: 163–77.

Choi, Jae S., B. Cameron, K. Christie, A. Glass, and E. MacEachern. 2022. Temperature and Depth Dependence of the Spatial Distribution of Snow Crab. bioRxiv. https://doi.org/10.1101/2022.12.20.520893.


Choi, Jae S. 2023. A Multi-Stage, Delay Differential Model of Snow Crab Population Dynamics in the Scotian Shelf of Atlantic Canada. bioRxiv. https://doi.org/10.1101/2023.02.13.528296.
 

DFO. 2013. [Integrated Fisheries Management Plan for Eastern Nova Scotia and 4X Snow Crab (*Chionoecetes Opillio*.)](http://www.dfo-mpo.gc.ca/fm-gp/peches-fisheries/ifmp-gmp/snow-crab-neige/snow-crab-neiges2013-eng.htm)


DFO. 2018. Stock Status Update of Atlantic Halibut (Hippoglossus hippoglossus) on the Scotian Shelf and Southern Grand Banks in NAFO Divisions 3NOPs4VWX5Zc. DFO Can. Sci. Advis. Sec. Sci. Resp. 2018/022.

Hebert M, Miron G, Moriyasu M, Vienneau R, and DeGrace P. Efficiency and ghost fishing of Snow Crab (Chionoecetes opilio) traps in the Gulf of St. Lawrence. Fish Res. 2001; 52(3): 143-153. 10.1016/S0165-7836(00)00259-9   
 
Riebler, A., Sørbye, S.H., Simpson D., and Rue, H. 2016. An intuitive Bayesian spatial model for disease mapping that accounts for scaling. Statistical methods in medical research 25: 1145-1165.

Simpson, D., Rue, H., Riebler, A., Martins, T.G., and Sørbye, SH. 2017. Penalising Model Component Complexity: A Principled, Practical Approach to Constructing Priors. Statist. Sci. 32: 1-28.
 
       sh Res. 2001; 52(3): 143-153. 10.1016/S0165-7836(00)00259-9   
 
Riebler, A., Sørbye, S.H., Simpson D., and Rue, H. 2016. An intuitive Bayesian spatial model for disease mapping that accounts for scaling. Statistical methods in medical research 25: 1145-1165.

Simpson, D., Rue, H., Riebler, A., Martins, T.G., and Sørbye, SH. 2017. Penalising Model Component Complexity: A Principled, Practical Approach to Constructing Priors. Statist. Sci. 32: 1-28.
 
       