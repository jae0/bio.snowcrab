---
title: "Snow Crab Advisory 4X"
subtitle: "Maritimes Region"

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


<!-- to render a presentation to a doctype:  DOCTYPE=revealjs, DOCTYPE=html, etc

make quarto FN=snowcrab_presentation_advisory_4x.md DOCTYPE=html PARAMS="-P year_assessment:2024 -P todo:[fishery_results,fishery_model,ecosystem,redo_data]"  --directory=~/bio/bio.snowcrab/inst/markdown

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
 
  load_bio_environment()

```

 

<!-- 
# NOTE:: _load_results.qmd contains instructions in "todo" 
# this is a shared R-script to boot strap and provide a consistent data interface
-->


{{< include _load_results.qmd >}}  



::: {.landscape}


## Terms of Reference

- Fishery performance

- Bycatch of non-target species

- Environmental and climate change considerations

- Stock status and trends

- Major sources of uncertainty, where applicable.
 

 
## Fishery performance


## Fishery statistics 
 

```{r}
#| label: tbl-fishery-performance-4X
#| echo: false 
#| eval: true
#| output: true
#| tbl-cap: "Fishery performance statistics: 4X -- years represent the starting year and currently ongoing."

r=3
  reg = regions[r]
  REG = reg_labels[r]
  ## cat(REG, "\n")
  oo = dt[ which(dt$Region==reg), c("Year", "Licenses", "TAC", "Landings", "Effort", "CPUE")] 

  names(oo) = c( "Year", "Licenses", "TAC (t)", "Landings (t)", "Effort (1000 th)", "CPUE (kg/th)" )

  gt::gt(oo) |> gt::tab_options(table.font.size = 20, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))



``` 

- Fishery ongoing at time of data capture
- TAC: decreased 64%
- Landings: decreased 78%
- Effort: decreased 90%
- Catch rates; increased 77%
 


## Effort

```{r}
#| label: fig-effort-timeseries
#| eval: true 
#| output: true
#| fig-dpi: 144
#| fig-height: 8
#| fig-cap: "Temporal variations in fishing effort."
 
if (params$sens==1) {
  ts_dir = file.path( data_loc, "assessments", year_assessment, "timeseries", "fishery" )
} else if (params$sens==2) {
  ts_dir = file.path( data_loc, "assessments", year_assessment, "timeseries", "fishery", "split" )
}

include_graphics( file.path( ts_dir, "effort.ts.png" ) )
 
```
 

```{r}
#| label: fig-effort-map
#| eval: true 
#| output: true
#| fig-dpi: 144
#| fig-height: 8
#| echo: false 
#| layout-ncol: 2
#| fig-cap: "Snow Crab fishing effort from fisheries logbook data for previous and current years (no X 10$^3$ per 10 km X 10 km grid)."
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""
 
loc0 = file.path( data_loc, "output", "maps", "logbook", "snowcrab", "annual", "effort" )
yrsplot = year_assessment + c(0:-3)
fn = file.path( loc0, paste( "effort", yrsplot, "png", sep=".") ) 

include_graphics( fn ) 

```
 


## Landings
 
```{r}
#| label: fig-landings-timeseries
#| eval: true 
#| output: true
#| fig-cap: "Landings (t) of Snow Crab on the SSE. For 4X, the year refers to the starting year of the season. Inset is a closeup view of the timeseries for N-ENS and 4X."
#| fig-dpi: 144
#| fig-height: 8

if (params$sens==1) {
  ts_dir = file.path( data_loc, "assessments", year_assessment, "timeseries", "fishery" )
} else if (params$sens==2) {
  ts_dir = file.path( data_loc, "assessments", year_assessment, "timeseries", "fishery", "split" )
}

include_graphics( file.path( ts_dir, "landings.ts.png" ) )
```
 
```{r}
#| label: fig-landings-map
#| eval: true 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| echo: false 
#| layout-ncol: 2
#| fig-cap: "Snow Crab landings from fisheries logbook data for previous and current years (tons per 10 km x 10 km grid)."
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

loc0 = file.path( data_loc, "output", "maps", "logbook", "snowcrab", "annual", "landings" )
yrsplot = year_assessment + c(0:-3)
fn = file.path( loc0, paste( "landings", yrsplot, "png", sep=".") ) 

include_graphics( fn ) 

``` 


## Catch rates
 
 
```{r}
#| label: fig-cpue-timeseries
#| eval: true 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig-cap: "Temporal variations in crude catch rates of Snow Crab (kg/trap haul)."
 
if (params$sens==1) {
  ts_dir = file.path( data_loc, "assessments", year_assessment, "timeseries", "fishery" )
} else if (params$sens==2) {
  ts_dir = file.path( data_loc, "assessments", year_assessment, "timeseries", "fishery", "split" )
}

include_graphics( file.path( ts_dir, "cpue.ts.png" ) ) 
```
 

```{r}
#| label: fig-cpue-map
#| eval: true 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| echo: false 
#| layout-ncol: 2
#| fig-cap: "Snow Crab crude catch rates on the Scotian Shelf for previous and current years. Units are kg/trap haul per 10 km x 10 km grid."
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

loc0= file.path( data_loc, "output", "maps", "logbook", "snowcrab", "annual", "cpue" )
yrsplot = year_assessment + c(0:-3)
fn = file.path( loc0, paste( "cpue", yrsplot, "png", sep=".") ) 

include_graphics( fn ) 
 
```  
 



## Bycatch of non-target species 


## At sea observed trips

```{r}
#| label: tbl-observed-summary
#| eval: true 
#| output: true
#| fig-dpi: 144
#| fig-height: 5
#| echo: false 
#| layout-ncol: 1
#| tbl-cap: "Table of observed data coverage"

fns = c( 
  "observersummary2.png"
)
 
include_graphics( file.path( media_supplementary, fns )  ) 

```
 

```{r}
#| label: fig-map-observer-locations
#| eval: true 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| echo: false 
#| layout-ncol: 2
#| fig-cap: "Snow Crab At-sea-observer locations."
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

loc = file.path( data_loc, "output", "maps", "observer.locations" )
yrsplot = year_assessment + c(0:-3) 
 
fns = paste( "observer.locations", yrsplot, "png", sep="." ) 
fn = file.path( loc, fns ) 

include_graphics( fn )

```
   

## Discard of (all) soft-shell crab map
 
 
  No data in 4X


```{r}
#| label: fig-observed-softshell-map
#| eval: true 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| echo: false 
#| layout-ncol: 2
#| fig-cap: "Map of observed soft shell locations."
# #| fig-subcap:  
# #|   - "N-ENS"
# #|   - "CFA23"
# #|   - "CFA24" 

fns = c( 
#  "nens_soft_crab_positions_68.png" #,
#  "cfa23_soft_crab_positions_68.png",
#  "cfa24_soft_crab_positions_68.png"
)
 
#include_graphics( file.path( media_supplementary, fns )  ) 

```
 


 

## Discard of non-target species ("Bycatch") 

 

```{r}
#| label: tbl-fishery-discard-effort-4x
#| eval: true
#| output: true
#| echo: false  
#| tbl-cap: "Bycatch (kg) estimated from fisheries effort. Dots indicate values less than 10 kg/year. Where species exist in a list but there is no data, this indicates some historical bycatch. The overall average is from 2004 to present." 

r = 3
  reg = regions[r]
  REG = reg_labels[r]
  # cat( REG, "\n")
  o = BC[[reg]]   
  oo = o$bycatch_table_effort
  oo[ oo==0 ] = NA
  oo[ is.na(oo) ] = "."
  gt::gt(oo) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))

 

```

- Historical average of 0.87% of landings; primarily other Crustacea (crab and lobster)  
 
## Environmental and climate change considerations


## Bottom temperatures


- Average bottom temperatures observed in the 2024 Snow Crab survey have returned to historical ranges, after worrisome highs in 2022. 

- Bottom temperatures are more stable in N-ENS than S-ENS; 4X exhibits the most erratic and highest annual mean bottom temperatures.


```{r}
#| label: fig-bottom-temperatures-timeseries
#| eval: true
#| echo: false 
#| output: true
#| fig-cap: "Annual variations in bottom temperature observed during the Snow Crab survey. The horizontal (black) line indicates the long-term, median temperature within each subarea. Error bars represent standard errors." 
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold 

tloc = file.path( data_loc, "assessments", year_assessment, "timeseries"  )

fns = c( 
  file.path("survey", "t.png")
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
#| label: fig-bottom-temperatures-timeseries-modelled
#| eval: true
#| echo: false 
#| output: true
#| fig-cap: "Temporal variations in bottom temperature from a historical reanalysis of temperature data. Red horizontal line is at $7^\\circ$C. Presented are 95% Credible Intervals of spatial variability in temperature at each time slice after adjustment for spatiotemporal autocorrelation."
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold 

tloc = file.path( data_loc, "assessments", year_assessment, "timeseries"  )

fns = c( 
  "temperature_bottom.png" 
)

include_graphics( file.path( tloc, fns) )

``` 
 

 

```{r}
#| label: fig-figures-temperature-bottom-map-modelled
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: "Bottom temperature ($^\\circ$C) estimated from a historical analysis of temperature data for 1 September."
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""
 
map_outdir = file.path( data_root, 'aegis', 'temperature', 'modelled', 'default', 'maps' )
map_years  = year_assessment + c(0:-3)
  
fn = check_file_exists( file.path( 
  map_outdir, paste( 'predictions.',  map_years, '.0.75',  '.png', sep='') 
) )

include_graphics( fn )
```

  
 
## Interspecific interactions

- Assumption: Elevated co-ocurrence in areal **densities** with snow crab trawl samples means greater potential encounter rates

```{r}
#| label: fig-speciescomposition-biplot
#| eval: true
#| echo: false 
#| output: true
#| fig.show: hold 
#| fig-dpi: 144
#| fig-height: 6
#| fig-cap: "Species ordination (PCA: eigenanalysis of correlation matrices). PC1 is associatd with bottom temperatures. PC2 is associated with depth. Snow crab is shown as an orange dot."
     
xlab = paste("PC1 (", pca$variance_percent[1], "%)", sep="" )
ylab = paste("PC2 (", pca$variance_percent[2], "%)", sep="" )
plot( PC2 ~ PC1, pcadata, type="n", xlab=xlab, ylab=ylab )
text( PC2 ~ PC1, labels=vern, data=pcadata, cex=0.7, col="slateblue"  )
i = grep("Snow crab", pcadata$vern, ignore.case=TRUE)
points( PC2 ~ PC1, pcadata[i,], pch=19, cex=3.0, col="darkorange" )

```


 

## Predators


```{r}
#| label: tbl-predators
#| echo: false
#| eval: true
#| output: true
#| tbl-cap: "Main predators based upon frequency of occuence of snow crab in finfish stomach samples, unadjusted for sampling effort."

gt::gt(counts[1:11,]) 

```


```{r}
#| label: fig-predator-biplot
#| eval: true
#| echo: false 
#| output: true
#| fig.show: hold 
#| fig-dpi: 144
#| fig-height: 6
#| fig-cap: "Main predators of snow crab on Scotian Shelf of Atlantic Canada (1999-2020). Relative location of snow crab predators (green) in the species composition ordination. Snow crab in orange. Of 58,287 finfish stomach samples, 159 had snow crab (0.28%). There is no information on snow crab diet in the database."
#| fig-subcap: 
#|   - "Total random effect (space and space-time) - PC1"
#|   - "Total random effect (space and space-time) - PC1"
#|   - "Mean annual PC1 score."
#|   - "Mean annual PC2 score."
#| layout-ncol: 1

# potential predators:
lookup= c( "cod", "halibut", "sculpin", "skate", "plaice", "hake", "wolffish", "atlantic cod", "atlantic halibut", "longhorn sculpin", "thorny skate", "striped atlantic wolffish", "haddock", "american plaice", "smooth skate", "winter skate", "white hake", "shorthorn sculpin", "eelpout newfoundland", "squirrel or red hake", "sea raven", "ocean pout", "barndoor skate", "Squid" )

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

- Overall, Atlantic Halibut densities have increased rapidly since 2010 (The Gully, Slope Edge and near Sable Island).

- Thorny skate densities have been increasing (especially in N-ENS and along the margins of Banquereau Bank)

- Striped Atlantic Wolffish densities have been high (declining in N-ENS since 2007, peaking in Laurentian Channel).

- Potential predators: in warmer and deeper (left and bottom of Snow Crab).


## Prey
 
- echinoderms
- polychaete worms (*Maldane*, *Nereis*), worm-like animals
- detritus (dead organic matter)
- large zooplankton, shrimp 
- juvenile crab (Rock Crab; Toad Crab; Lesser Toad Crab)
- Ocean Quahog (*Artica islandica*), bivalve molluscs (*Mytilus* sp, *Modiolus*, *Hiatella*)
- brittle stars (*Ophiura*, *Ophiopholis*)
- sea anemones (*Edwardsia*, *Metridium*). 


```{r}
#| label: fig-diet-biplot
#| eval: true
#| echo: false 
#| output: true
#| fig.show: hold 
#| fig-dpi: 144
#| fig-height: 6
#| fig-cap: "Relative location of snow crab prey (green) in the species composition ordination. Snow crab in orange. Most of the potental prey are found to the right of snow crab (i.e. colder-water species) at a variety of depths."
#| fig-subcap: 
#|   - "Total random effect (space and space-time) - PC1"
#|   - "Total random effect (space and space-time) - PC1"
#|   - "Mean annual PC1 score."
#|   - "Mean annual PC2 score."
#| layout-ncol: 1

# potential food items:
lookup= c( "echinoderm", "polychaete", "maldane", "nereis", "shrimp", "pandalus", "rock crab", "toad crab", "lesser toad crab", "quahog", "artica islandica", "mollusc", "mytilus", "modiolus", "hiatella", "starfish", "sea anemone", "brittle star", "sea star", "sea anemone", "ophiura", "ophiopholis", "edwardsia", "metridium", "euphasid" )

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

- Potential prey: in colder (right of Snow Crab). 
 

## Competitors
 

```{r}
#| label: fig-competitor-biplot
#| eval: true
#| echo: false 
#| output: true
#| fig.show: hold 
#| fig-dpi: 144
#| fig-height: 6
#| fig-cap: "Potential competitors of snow crab on Scotian Shelf of Atlantic Canada (1999-2020). Relative location of snow crab predators (green) in the species composition ordination. Snow crab in orange."
#| fig-subcap: 
#|   - "Total random effect (space and space-time) - PC1"
#|   - "Total random effect (space and space-time) - PC1"
#|   - "Mean annual PC1 score."
#|   - "Mean annual PC2 score."
#| layout-ncol: 1

# potential predators:

lookup= c( "pandalus", "Jonah Crab", "Atlantic Rock Crab" , "Toad Crab", "Hyas Coarctatus", "Northern Stone Crab" )   # add more here: most are not direct compeitors as they have slightly different depth/temp preferences


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

- Northern shrimp co-occur as they share similar habitat preferences. Numerical densities have declined after a peak in 2011, especially in S-ENS.

- Lesser toad crab is a co-occurring species and potential competitor. Their numbers have declined to low levels throughout, after a peak in densities in 2007 and 2016 in N-ENS.



## Viable habitat

```{r}
#| label: fig-fb-habitat-timeseries
#| eval: true
#| echo: false 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "Habitat viability (probability; fishable Snow Crab). Means and 95\\% Credible Intervals are presented."

loc = file.path( data_loc, "modelled", "default_fb", "aggregated_habitat_timeseries" )
include_graphics( file.path( loc, "habitat_M0.png") )

```


```{r}
#| label: fig-fb-habitat-map
#| eval: true
#| echo: false 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "Habitat viability (probability; fishable Snow Crab)."
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
#| layout: [[100], [100], [50,50], [50,50], [50,50], [50,50]]
  
loc = file.path( data_loc, "modelled", "default_fb", "predicted_habitat" )
vn = "habitat."
yrsplot =  year_assessment + c(0:-9)

fns = file.path( loc, paste( vn, yrsplot, ".png", sep="") )
include_graphics( fns )
 
```


- Nouanced: not just temperature, but joint distribution of temperature and other factors (substrate, depths, co-occuring species, ...)

- SSE is variable  Warm-Water incursions resulted in predators co-occupying Snow Crab habitats while simultaneously losing access to cold water preferring prey. 
   
 
- 4X showed a slight improvement in habitat, however, the overall trend has been downwards since 2010 . Even with amelioration of bottom temperatures in 2024, previous habitat space seems to have been overtaken by competitors and predators in 4X. Peaks in area south of Lunenburg.
  

## Stock status and trends

## Data sources: Trawl Survey

Marine Protected Areas (St Ann's, Gully) were not sampled
  
- CFA 4X: 24 stations completed 


```{r}
#| label: fig-survey-locations-map 
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| echo: false 
#| layout-ncol: 2
#| fig-cap: "Survey locations."
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

loc = file.path( data_loc, "output", "maps", "survey.locations" )
years = year_assessment + c(0:-3)
fn = check_file_exists( file.path( loc, paste( "survey.locations", years, "png", sep=".") ))
include_graphics( fn )
```

<!-- 

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
knitr::include_graphics( fn1 ) 

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
knitr::include_graphics( fn1 ) 

```


## Growth modes

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
-->

  
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
 
 
 

## Recruitment: males


```{r}
#| label: fig-sizefeq-male
#| eval: true
#| output: true
#| fig-cap: "Size-frequency (areal density; no/km$^2$) histograms by carapace width of male Snow Crab. The vertical line represents the legal size (95 mm). Immature animals are shown with light coloured bars, mature with dark."
#| fig-dpi: 144
#| fig-height: 10

if (params$sens==1) {
  sf_outdir = file.path( p$annual.results, "figures", "size.freq", "survey", "period1")
} else if (params$sens==2) {
  sf_outdir = file.path( p$annual.results, "figures", "size.freq", "survey_split", "period1")
}

include_graphics( file.path( sf_outdir,  "male.denl.png" ) )

```
  
- CFA 4X: erratic inter-annual patterns of low recruitment. The small mode near 68 mm CW (~instar 10 and so 1-3 years to entry to fishery)
   

## Reproduction (longer term)
 
```{r}
#| label: fig-sizefeq-female
#| eval: true
#| output: true
#| fig-cap: "Size-frequency (areal density; no/km$^2$) histograms by carapace width of female Snow Crab. Immature animals are shown with light coloured bars, mature with dark."
#| fig-dpi: 144
#| fig-height: 10


if (params$sens==1) {
  sf_outdir = file.path( p$annual.results, "figures", "size.freq", "survey", "period1")
} else if (params$sens==2) {
  sf_outdir = file.path( p$annual.results, "figures", "size.freq", "survey_split", "period1")
}

fn = file.path( sf_outdir, "female.denl.png" )

include_graphics( fn )
```

 
4X:
    - decline in numerical densities of both the mature and adolescent females since 2017
    - egg and larval production is expected to be moderate in the next year
    - A notable absence of immature females is evident in 4X: predation and habitat loss?
 
  

```{r}
#| label: fig-totno-female-mat-timeseries
#| eval: true
#| output: true
#| fig-cap: "The crude, unadjusted geometric mean of mature female density log$_{10}$(no/km$^2$) from the Snow Crab survey. Error bars represent 95\\% Confidence Intervals. Note the absence of data in 2020. Prior to 2004, surveys were conducted in the Spring."
#| fig-dpi: 144
#| fig-height: 4 
 
 
if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

fn = file.path( ts_outdir, paste("totno.female.mat", "png", sep=".") )
include_graphics( fn )

```
 

```{r}
#| label: fig-totno-female-mat-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: "Mature female density log$_{10}$(no/km$^2$) from the Snow Crab survey."
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""
 
map_outdir = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual" )
map_years  = year_assessment + c(0:-3)
  
fn = check_file_exists( file.path( 
  map_outdir, "totno.female.mat", paste( "totno.female.mat", map_years, "png", sep="." )  
) )

include_graphics( fn )
```

<!--

## Sex ratios

- Unbalanced sex ratios suggest that immature males are exposed to elevated natural mortality (predation, cannibalism, habitat loss?)

   
```{r}
#| label: fig-sexratio-mat-timeseries
#| eval: true
#| output: true
#| fig-cap: "The crude, unadjusted geometric mean of sex ratios (proportion female) of mature Snow Crab. Error bars represent 95\\% Confidence Intervals. Note the absence of data in 2020. Prior to 2004, surveys were conducted in the Spring."
#| fig-dpi: 144
#| fig-height: 4 
 
if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

fn = file.path( ts_outdir, paste("sexratio.mat", "png", sep=".") )
include_graphics( fn )

```
 

 

```{r}
#| label: fig-sexratio-mat-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: "Snow Crab survey sex ratios (proportion female) of mature Snow Crab."
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""
 

map_outdir = file.path( p$project.outputdir, "maps", "survey", "snowcrab", "annual" )
map_years  = year_assessment + c(0:-3)
  
fn = check_file_exists( file.path( 
  map_outdir, "sexratio.mat", paste( "sexratio.mat", map_years, "png", sep="." )  
) )

include_graphics( fn )
```
  
-->



## Fishable biomass density

- CFA 4X: densities have declined to historical lows; peak only in Mahonr Bay area (south of Lunenburg)
- Biomass density, however, does not equate to total biomass as the areas occupied by crab can contract, expand and shift with environmental conditions and ecosystem variability. 
 


```{r}
#| label: fig-R0-timeseries
#| eval: true
#| output: true
#| fig-cap: "The crude, unadjusted geometric mean fishable biomass density log$_{10}$(t/km$^2$) from the Snow Crab survey. Error bars represent 95\\% Confidence Intervals. Note the absence of data in 2020. Prior to 2004, surveys were conducted in the Spring."
#| fig-dpi: 144
#| fig-height: 4 
 
if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

fn = file.path( ts_outdir, paste("R0.mass", "png", sep=".") )
include_graphics( fn )

```
 


```{r}
#| label: fig-R0-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: "Snow Crab survey fishable component biomass density log$_{10}$(t/km$^2$)."
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""
 

map_outdir = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual" )
map_years  = year_assessment + c(0:-3)
  
fn = check_file_exists( file.path( 
  map_outdir, "R0.mass", paste( "R0.mass", map_years, "png", sep="." )  
) )

include_graphics( fn )
```

 
 

## Fishable biomass index

- CFA 4X: one peak (2009) increased marginally; area south of Sambro and Lunenburg


```{r}
#| label: fig-fbindex-timeseries
#| eval: true
#| echo: false 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "The fishable biomass index (t) predicted by CARSTM of Snow Crab survey densities. Error bars represent Bayesian 95\\% Credible Intervals. Note large errors in 2020 when there was no survey."

fn = file.path( data_loc, "modelled", "default_fb", "aggregated_biomass_timeseries" , "biomass_M0.png")
include_graphics( fn )

```
 

```{r}
#| label: fig-fbindex-map
#| eval: true
#| echo: false 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "Biomass index log~10(t/km$^2$) predicted from the Snow Crab survey."
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
#| layout: [[100], [100], [50,50], [50,50], [50,50], [50,50]]
  
loc = file.path( data_loc, "modelled", "default_fb", "predicted_biomass_densities" )
yrsplot =  year_assessment + c(0:-9)

fns = file.path( loc, paste( "biomass", yrsplot, "png", sep=".") )
include_graphics( fns )
```
  

## Modelled pre-fishery fishable biomass
 
- CFA 4X: 0.18t, relative to 0.14kt in the previous season


```{r}
#| label: fig-logisticPredictions
#| eval: true
#| echo: false
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "Fishable, posterior mean modelled biomass (pre-fishery; kt) are shown in dark orange. Light orange are posterior samples of modelled biomass (pre-fishery; kt) to illustrate the variability of the predictions. The biomass index (post-fishery, except prior to 2004) after model adjustment by the model catchability coefficient is in gray."
# #| fig-subcap:
# #|   - "N-ENS"
# #|   - "S-ENS"
# #|   - "4X"

loc = file.path( data_loc, "fishery_model", year_assessment, "logistic_discrete_historical" )
fns = file.path( loc, c(
  # "plot_predictions_cfanorth.png"
  #,
  #"plot_predictions_cfasouth.png",
  "plot_predictions_cfa4x.png"
) )

include_graphics( fns )

```

<!--

```{r}
#| label: fig-logisticPredictions2
#| eval: true
#| echo: false
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "Fishable, posterior mean modelled biomass (post-fishery; kt) are shown in dark orange. Light orange are posterior samples of modelled biomass (post-fishery; kt) to illustrate the variability of the predictions. The biomass index (post-fishery, except prior to 2004) after model adjustment by the model catchability coefficient is in gray."
#| fig-subcap:
#|   - "N-ENS"
#|   - "S-ENS"
#|   - "4X"

loc = file.path( data_loc, "fishery_model", year_assessment, "logistic_discrete_historical" )
fns = file.path( loc, c(
  "plot_predictions_postfishery_cfanorth.png",
  "plot_predictions_postfishery_cfasouth.png",
  "plot_predictions_postfishery_cfa4x.png"
) )

include_graphics( fns )

```

-->
 

## Posterior estimates of fishing mortality
 
- In 4X, the 2024–2025 season (ongoing) is estimated have F=0.052 (annual exploitation rate of 5%), while in the previous season it was F=0.27 (annual exploitation rate of 31%).
 
- Localized exploitation rates are likely higher, as not all areas for which biomass is estimated are fished (e.g., continental slope areas and western, inshore areas of CFA 24, 4X).



```{r}
#| label: fig-logisticFishingMortality
#| eval: true
#| echo: false
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "Time-series of modelled instantaneous fishing mortality from Model 1, for 4X. Samples of the posterior densities are presented, with the darkest line being the mean."
# #| fig-subcap:
# #|   - "N-ENS"
# #|   - "S-ENS"
# #|   - "4X"

odir = file.path( fishery_model_results, year_assessment, "logistic_discrete_historical" )
fns = file.path( odir, c(
  # "plot_fishing_mortality_cfanorth.png"
  #,
  #"plot_fishing_mortality_cfasouth.png",
  "plot_fishing_mortality_cfa4x.png"
))

include_graphics( fns )

```


## Reference Points

- Limit Reference Point (LRP) = K/4
- Upper Stock Reference (USR) = K/2
- Removal Reference (RR) = F_{MSY} = r/ 2

- 4X is very clearly in the “critical” zone


 
|   | N-ENS | S-ENS | 4X |
|----- | ----- | ----- | ----- |
| |  |  |  |
|q       | `r round(q_north, 3)` (`r round(q_north_sd, 3)`) | `r round(q_south, 3)` (`r round(q_south_sd, 3)`) | `r round(q_4x, 3)` (`r round(q_4x_sd, 3)`) |
|r       | `r round(r_north, 3)` (`r round(r_north_sd, 3)`) | `r round(r_south, 3)` (`r round(r_south_sd, 3)`) | `r round(r_4x, 3)` (`r round(r_4x_sd, 3)`) |
|K       | `r round(K_north, 2)` (`r round(K_north_sd, 2)`) | `r round(K_south, 2)` (`r round(K_south_sd, 2)`) | `r round(K_4x, 2)` (`r round(K_4x_sd, 2)`) |
|Prefishery Biomass   | `r round(B_north[t0], 2)` (`r round(B_north_sd[t0], 2)`) | `r round(B_south[t0], 2)`  (`r round(B_south_sd[t0], 2)`) | `r round(B_4x[t0], 2)`  (`r round(B_4x_sd[t0], 2)`)  |
|Fishing Mortality    | `r round(FM_north[t0], 3)` (`r round(FM_north_sd[t0], 3)`) | `r round(FM_south[t0], 3)` (`r round(FM_south_sd[t0], 3)`) | `r round(FM_4x[t0], 3)` (`r round(FM_4x_sd[t0], 3)`) |

: Reference points from the logistic biomass dynamics fishery model. K is Carrying capacity (kt); and r is Intrinsic rate of increase (non-dimensional). Note that FMSY (fishing mortality associated with 'Maximum Sustainable Yield') is r/2. Similarly, BMSY (biomass associated with 'Maximum Sustainable Yield') is K/2. SD is posterior Standard deviations.* {#tbl-reference-points}
 



```{r}
#| label: fig-ReferencePoints
#| eval: true
#| echo: false
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "Harvest control rules for the Scotian Shelf Snow Crab fisheries."

fn = file.path( media_loc, "harvest_control_rules.png")
include_graphics( fn )

```
 

```{r}
#| label: fig-logistic-hcr
#| eval: true
#| echo: false
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "Reference Points (fishing mortality and modelled biomass) from the Fishery Model for 4X. The large yellow dot indicates most recent year and the 95\\% CI. Not: the model does not account for illegal and unreported landings, and interspecific interactions. Prefishery."
# #| fig-subcap:
# #|   - "N-ENS"
# #|   - "S-ENS"
# #|   - "4X"

odir = file.path( fishery_model_results, year_assessment, "logistic_discrete_historical" )

fns = file.path( odir, c(
  # "plot_hcr_cfanorth.png" # ,
  # "plot_hcr_cfasouth.png",
  "plot_hcr_cfa4x.png"
) )

include_graphics( fns )

```
 




## Major sources of uncertainty

 
- Capture of soft-shell Snow Crab (handling mortality)

- Bycatch of Snow Crab in other fisheries (use as bait)

- Illegal, unreported, and unregulated fishing activities 

- Marine Protected Areas (MPAs) act as a refuge from fishing activities. However, positive effects upon other organisms (predators or prey) can have counter-balancing indirect effects. 

 
## Conclusions  

The SSE continues to experience rapid ecosystem and climatic variations. Under such conditions, it is prudent to be careful. 


4X: 
    - Recruitment is expected to be low for another three years. 
    - Viable habitat has been minimal for many years (since 2011). 
    - Total mortality is now in approximate balance with recruitment. 
    - In the “critical” zone.


## Acknowledgements
 
- Innovative Capital Investments Inc.

- Captain Coalie D’Eon and crew of the F/VRS Journey II for their expertise and provision of a safe and hospitable environment for the conduct of the survey.

- Snow Crab license holders and fishers of the SSE a stellar model of a collaborative, sustainable and precautionary co-management of a fishery


## END

:::