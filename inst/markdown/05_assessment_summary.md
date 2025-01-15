---
title: "Snow crab assessment summary"
author:
  - name: 
      given: Snow Crab Unit
      family: DFO Science
    # orcid: 0000-0003-3632-5723 
    # email: jae.choi@dfo-mpo.gc.ca
    # email: choi.jae.seok@gmail.com
    # corresponding: true
    affiliation: 
      - name: Bedford Institute of Oceanography, Fisheries and Oceans Canada
        city: Dartmouth
        state: NS
        # url: www.bio.gc.ca
date: last-modified
date-format: "YYYY-MM-D"
keywords: 
  - snow crab fishery model result and assessment of status
abstract: |
  Snow crab fishery model results and stock status assessment.
toc: true
toc-depth: 4
number-sections: true
highlight-style: pygments
# bibliography: media/references.bib  
# csl: media/canadian-journal-of-fisheries-and-aquatic-sciences.csl  # see https://www.zotero.org/styles for more
# license: "CC BY"
copyright: 
  holder: snow-crab-unit
  year: 2024
# citation: 
#  container-title: https://github.com/jae0/bio.snowcrab/
#  doi: NA
funding: "The snow crab scientific survey was funded by the snow crab fishers of Maritimes Region of Atlantic Canada."
editor:
  render-on-save: false
execute:
  echo: true
format:
  html: 
    code-fold: true
    code-overflow: wrap
    html-math-method: katex
    self-contained: true
    embed-resources: true
params:
  year_assessment: 2024
  year_start: 1999
  media_loc: "media"
  sens: 1
  debugging: FALSE
  model_variation: logistic_discrete_historical
---


 
<!-- 

## Preamble

Summary 3 of 3 -- This file is designed to be an HTML document that describes and summarizes the assessment of stock status. 


cd ~/bio/bio.snowcrab/inst/markdown

make quarto FN=05_assessment_summary YR=2024 SOURCE=~/bio/bio.snowcrab/inst/markdown WK=~/bio.data/bio.snowcrab/assessments DOCEXTENSION=html PARAMS="-P year_assessment:2024"


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
  # dev.args = list(type = "cairo"),
  fig.retina = 2,
  dpi=192
)

require(spsUtil)
quietly = spsUtil::quiet

require(ggplot2)
require(MBA)
require(gt)  # table formatting
require(aegis)  # basic helper tools
    
loadfunctions( "aegis")
loadfunctions( "bio.snowcrab")  # in case of local edits

year_assessment = params$year_assessment
model_variation = params$model_variation
media_loc = params$media_loc


p = load.environment( year.assessment=year_assessment )
p$corners = data.frame(plon=c(220, 990), plat=c(4750, 5270) )
p$mapyears = year_assessment + c(-5:0 )   # default in case not specified

year_previous = year_assessment - 1

yrs = 2000:year_assessment
years = as.character(yrs)

lregions = list(region=c("cfanorth", "cfasouth", "cfa4x"))
reg_labels = c("N-ENS", "S-ENS", "CFA 4X")  # formatted for label

if (params$sens==2) {
  lregions = list(subarea=c("cfanorth", "cfa23",  "cfa24", "cfa4x"))
  reg_labels = c("CFA 20-22", "CFA 23", "CFA 24", "CFA 4X")  # formatted for label
}

regions = unlist(lregions)
nregions = length(regions)

# directories
outtabledir = file.path( p$annual.results, "tables" )
SCD = project.datadirectory("bio.snowcrab")

# fishery_model_results = file.path( "/home", "jae", "projects", "dynamical_model", "snowcrab", "outputs" )
fishery_model_results = file.path( SCD, "fishery_model" )
  
```


&nbsp;  $~$  <br /> 


# Management areas

```{r}
#| label: fig-management-areas
#| eval: true
#| echo: false 
#| output: true
#| fig-cap: "The Scotian Shelf (NW Atlantic Ocean; NAFO Div. 4VWX). Shown are isobaths and major bathymetric features. Managed Crab Fishing Areas (CFAs; divided by dashed lines) include: NENS, SENS, 4X. SENS is further subdivided (dotted line) into 23 (NW) and 24 (SE)."
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold 

fn = file.path( media_loc, "snowcrab_cfas.png" )
knitr::include_graphics( fn ) 

```


&nbsp;  $~$  <br /> 


# Ecosystem

## Connectivity: Oceanic currents

```{r}
#| label: fig-ocean-currents
#| eval: true
#| echo: false 
#| output: true
#| fig-cap: "Ocean currents in the Martimes. <br />Source: https://www.dfo-mpo.gc.ca/oceans/publications/soto-rceo/2018/atlantic-ecosystems-ecosystemes-atlantiques/index-eng.html"
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold 

fn = file.path( media_loc, c(
  "maritimes_currents.png" 
) )
knitr::include_graphics( fn ) 

```

&nbsp;  $~$  <br /> 


## Bathymetry
  
```{r}
#| label: fig-bathymetry
#| eval: true
#| echo: false 
#| output: true
#| fig-cap: "Bathymetry (mean and standard deviation) in the Scotian Shelf region from an stmv analysis."
#| fig-subcap: 
#|   - "Predicted depth (log$_{10}$ m)"
#|   - "Predicted slope (log$_{10}$ m / m)"
#|   - "Predicted curvature (log$_{10}$ m / m$^2$)"
#|   - "Obervation standard deviation (log$_{10}$ m)"
#|   - "Spatial standard deviation (log$_{10}$ m)"
#|   - "Spatial range (log$_{10}$ km)"
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold 

bdir = file.path(data_root, 'aegis', 'bathymetry', 'modelled', 'default', 'stmv', 'none_fft', 'z', 'maps', 'canada.east.superhighres' )
fns = file.path( bdir, c(
  "bathymetry.z.canada.east.superhighres.png",
  "bathymetry.dZ.canada.east.superhighres.png",
  "bathymetry.ddZ.canada.east.superhighres.png",
  "bathymetry.b.sdObs.canada.east.superhighres.png",
  "bathymetry.b.sdSpatial.canada.east.superhighres.png",
  "bathymetry.b.localrange.canada.east.superhighres.png"
) )

knitr::include_graphics( fns) 

``` 

&nbsp;  $~$  <br /> 


## Substrate

```{r}
#| label: fig-substrate
#| eval: true
#| echo: false 
#| output: true
#| fig-cap: "Substrate grain size log(mm) variations (mean and standard deviation) in the Scotian Shelf region from a carstm analysis."
#| fig-subcap: 
#|   - "Mean prediction"
#|   - "Standard deviation"
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold 

substrdir = file.path( data_root, 'aegis', 'substrate', 'modelled', 'default', 'maps' )
fns = file.path( substrdir, c(
  "predictions.png",
  "space_re_total.png"
) )

knitr::include_graphics( fns) 

```


&nbsp;  $~$  <br /> 



## Bottom Temperature 

```{r}
#| label: fig-bottom-temperatures
#| eval: true
#| echo: false 
#| output: true
#| fig-cap: "Bottom temperatures"
#| fig-subcap: 
#|   - "Annual variations in bottom temperature observed during the Snow Crab survey. The horizontal (black) line indicates the long-term, median temperature within each subarea. Error bars represent standard errors."
#|   - "Posterior densities of predicted average bottom temperatures. Red horizontal line is at $7^\\circ$C."
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold 

tloc = file.path( SCD, 'assessments', year_assessment, 'timeseries'  )

fns = c( 
  file.path("survey", "t.png"), 
  "temperature_bottom.png" 
)

knitr::include_graphics( file.path( tloc, fns) )

```
   


&nbsp;  $~$  <br /> 


## Ecosystem change: Bottom Temperature ... {.c}

- Average bottom temperatures **observed** in the 2022 Snow Crab survey were near or above historical highs in all areas 

- Temperatures are more stable in N-ENS than S-ENS; 4X exhibits the most erratic and highest annual mean bottom temperatures. 

- Observed temperatures in the 2022 Snow Crab survey for S-ENS increased well above the average. 
 
- Average temperature increased well beyond the $7^\circ$C threshold in 4X. N-ENS and S-ENS also continued to experience historical highs in bottom temperature and elevated spatial variability of bottom temperatures.



## Ecosystem considerations: Bottom Temperature ... {.c}

```{r fig-bottom-temperatures-map, out.width='30%', echo=FALSE, fig.show='hold', fig.align='center', fig.cap = 'Spatial variations in bottom temperature estimated from a historical analysis of temperature data for 1 September.' }

loc = file.path( data_root, 'aegis', 'temperature', 'modelled', 'default', 'maps' )
yrsplot =  year_assessment + c(0:-10)
fn10 = file.path( loc, paste( 'predictions.',  yrsplot[10], '.0.75',  '.png', sep='') )
fn9  = file.path( loc, paste( 'predictions.',  yrsplot[9],  '.0.75',  '.png', sep='') )
fn8  = file.path( loc, paste( 'predictions.',  yrsplot[8],  '.0.75',  '.png', sep='') )
fn7  = file.path( loc, paste( 'predictions.',  yrsplot[7],  '.0.75',  '.png', sep='') )
fn6  = file.path( loc, paste( 'predictions.',  yrsplot[6],  '.0.75',  '.png', sep='') )
fn5  = file.path( loc, paste( 'predictions.',  yrsplot[5],  '.0.75',  '.png', sep='') )
fn4  = file.path( loc, paste( 'predictions.',  yrsplot[4],  '.0.75',  '.png', sep='') )
fn3  = file.path( loc, paste( 'predictions.',  yrsplot[3],  '.0.75',  '.png', sep='') )
fn2  = file.path( loc, paste( 'predictions.',  yrsplot[2],  '.0.75',  '.png', sep='') )
fn1  = file.path( loc, paste( 'predictions.',  yrsplot[1],  '.0.75',  '.png', sep='') )
knitr::include_graphics( c( fn3, fn2, fn1) )
# \@ref(fig:bottom-temperatures-map)
# *Spatial variations in bottom temperature estimated from a historical analysis of temperature data for 1 September.*
```

 * Groundfish surveys were not conducted in 2020 and 2022 in the snow crab domain.

 * Snow crab surveys were not conducted in 2020, and incomplete in 2022 for S-ENS.


 
## Ecosystem considerations: Bottom Temperature ... {.c}

Persistent spatial gradient of almost $15^\circ$C in bottom temperatures in the Maritimes Region. 

Variable due to confluence of the warm, high salinity Gulf Stream from the S-SE along the shelf edge; cold, low salinity Labrador Current; and cold low salinity St. Lawrence outflow from the N-NE, as well as a nearshore Nova Scotia current, running from the NE. 

```{r fig-bottom-temperatures-spatialeffect, out.width='35%', echo=FALSE, fig.align='center', fig.cap = 'Persistent spatial effect of bottom temperature, relative to the overall mean, after adjustment for spatiotemporal variability and autocorrelations. Time period from 1999 to present.' }
loc = file.path( data_root, 'aegis', 'temperature', 'modelled', 'default', 'maps' )
knitr::include_graphics( file.path( loc, 'space_re_total.png') )
# \@ref(fig:bottom-temperatures-spatialeffect)
```
 
## Ecosystem considerations: Bottom Temperature ... {.c}
  
```{r fig-bottom-temperatures, out.width='65%', echo=FALSE, fig.align='center', fig.cap = '' }
knitr::include_graphics( file.path( SCD, 'assessments', year_assessment, 'timeseries', 'temperature_bottom.png') )
# \@ref(fig:bottom-temperatures)
```
 
Temporal variations in bottom temperature estimated from a historical analysis of temperature data. Red horizontal line is at $7^\circ$C. Presented are 95\% Credible Intervals of spatial variability in temperature at each time slice, after adjustment for spatiotemporal autocorrelation.
 

## Ecosystem change: Predators {.c}

- SSE: many species changes, snow crab are long-lived so interact with many of them
- Stomach samples: Atlantic Halibut, Atlantic Wolffish, Thorny Skate as primary predators
- Skewed sex ratios (few mature females): possibly linked to differential predation mortality (mature females being much smaller and full of fat-rich gonads and eggs). 
- Halibut (DFO 2018) have significantly increased in abundance in the Region. 
- Elevated co-ocurrence in areal **densities** with snow crab trawl samples means greater encounter rates





## Species composition
```{r fig-speciesomposition0, echo=FALSE, out.width='75%', fig.align='center', fig.show='hold', fig.cap = 'Species ordination (PCA: eigenanalysis of correlation matrices). PC1 is associatd with bottom temperatures. PC2 is associated with depth. Snow crab is shown as an orange dot.' }
xlab = paste("PC1 (", pca$variance_percent[1], "%)", sep="" )
ylab = paste("PC2 (", pca$variance_percent[2], "%)", sep="" )
plot( PC2 ~ PC1, pcadata, type="n", xlab=xlab, ylab=ylab )
text( PC2 ~ PC1, labels=vern, data=pcadata, cex=0.7, col="slateblue"  )
i = grep("Snow crab", pcadata$vern, ignore.case=TRUE)
points( PC2 ~ PC1, pcadata[i,], pch=19, cex=3.0, col="darkorange" )
```

 

## Species composition PC1


```{r fig-speciesomposition1, echo=FALSE, out.width='45%', fig.align='center', fig.show='hold',  fig.cap = 'Species composition in space and time. Primary gradient is related to bottom temperatures. Groundfish surveys were not conducted in 2020 and 2022 in the snow crab domain. Snow crab surveys were not conducted in 2020, and incomplete in 2022 in S-ENS.' }
spc_loc = file.path( data_root, 'aegis', 'speciescomposition', 'modelled', 'default', 'maps' )
fn1 = file.path( spc_loc, 'pca1.space_re_total.png') 
fn2 = file.path( spc_loc, 'pca2.space_re_total.png')
ts_loc = file.path( data_root, 'aegis', 'speciescomposition', 'modelled', 'default', 'figures' )
fn3 = file.path( ts_loc, 'pca1_time.png') 
fn4 = file.path( ts_loc, 'pca2_time.png') 
knitr::include_graphics( c(fn1,  fn3  ) ) 
# \@ref(fig:habitat3)  
``` 

 

## Species composition PC2

```{r fig-speciesomposition2, echo=FALSE, out.width='45%', fig.align='center', fig.show='hold',  fig.cap = 'Species composition in space and time. Primary gradient is related to bottom temperatures. Groundfish surveys were not conducted in 2020 and 2022 in the snow crab domain. Snow crab surveys were not conducted in 2020, and incomplete in 2022 in S-ENS.' }
spc_loc = file.path( data_root, 'aegis', 'speciescomposition', 'modelled', 'default', 'maps' )
fn1 = file.path( spc_loc, 'pca1.space_re_total.png') 
fn2 = file.path( spc_loc, 'pca2.space_re_total.png')
ts_loc = file.path( data_root, 'aegis', 'speciescomposition', 'modelled', 'default', 'figures' )
fn3 = file.path( ts_loc, 'pca1_time.png') 
fn4 = file.path( ts_loc, 'pca2_time.png') 
knitr::include_graphics( c(fn2,  fn4  ) ) 
# \@ref(fig:habitat3)  
``` 
  
  

## Prey


```{r fig-diet, echo=FALSE, out.width='100%', fig.align='center', fig.show='hold',  fig.cap = 'Relative location of snow crab prey (green) in the species composition ordination. Snow crab in orange.' }
xlab = paste("PC1 (", pca$variance_percent[1], "%)", sep="" )
ylab = paste("PC2 (", pca$variance_percent[2], "%)", sep="" )
plot( PC2 ~ PC1, pcadata, type="n", xlab=xlab, ylab=ylab )
text( PC2 ~ PC1, labels=vern, data=pcadata, cex=0.75, col="slategrey"  )
i = grep("Snow crab", pcadata$vern, ignore.case=TRUE)
points( PC2 ~ PC1, pcadata[i,], pch=19, cex=3.0, col="darkorange" )
lookup= c( "echinoderm", "polychaete", "maldane", "nereis", "shrimp", "pandalus", "rock crab", "toad crab", "lesser toad crab", "quahog", "artica islandica", "mollusc", "mytilus", "modiolus", "hiatella", "starfish", "sea anemone", "brittle star", "sea star", "sea anemone", "ophiura", "ophiopholis", "edwardsia", "metridium", "euphasid" )
j = NULL
for (k in lookup) j = c(j, grep( k, pcadata$vern, ignore.case=TRUE))    
j = unique(j)
points( PC2 ~ PC1, pcadata[j,], pch=19, cex=2.0, col="lightgreen" )
text( PC2 ~ PC1, labels=vern, data=pcadata[j,], cex=0.75, col="darkgreen"  )
```


- echinoderms
- polychaete worms (*Maldane*, *Nereis*), worm-like animals
- detritus (dead organic matter)
- large zooplankton, shrimp 
- juvenile crab (Rock Crab; Toad Crab; Lesser Toad Crab)
- Ocean Quahog (*Artica islandica*), bivalve molluscs (*Mytilus* sp, *Modiolus*, *Hiatella*)
- brittle stars (*Ophiura*, *Ophiopholis*)
- sea anemones (*Edwardsia*, *Metridium*). 




## Predators {.c}


```{r fig-predator_ord, echo=FALSE, out.width='100%', fig.align='center', fig.show='hold',  fig.cap = 'Main predators of snow crab on Scotian Shelf of Atlantic Canada (1999-2020). Relative location of snow crab predators (green) in the species composition ordination. Snow crab in orange. Of 58,287 finfish stomach samples, 159 had snow crab (0.28%). There is no information on snow crab diet in the database.' }
xlab = paste("PC1 (", pca$variance_percent[1], "%)", sep="" )
ylab = paste("PC2 (", pca$variance_percent[2], "%)", sep="" )
plot( PC2 ~ PC1, pcadata, type="n", xlab=xlab, ylab=ylab )
text( PC2 ~ PC1, labels=vern, data=pcadata, cex=0.75, col="slategrey"  )
i = grep("Snow crab", pcadata$vern, ignore.case=TRUE)
points( PC2 ~ PC1, pcadata[i,], pch=19, cex=3.0, col="darkorange" )
lookup= c( "cod", "halibut", "sculpin", "skate", "plaice", "hake", "wolffish", "atlantic cod", "atlantic halibut", "longhorn sculpin", "thorny skate", "striped atlantic wolffish", "haddock", "american plaice", "smooth skate", "winter skate", "white hake", "shorthorn sculpin", "eelpout newfoundland", "squirrel or red hake", "sea raven", "ocean pout", "barndoor skate" )
j = NULL
for (k in lookup) j = c(j, grep( k, pcadata$vern, ignore.case=TRUE))    
j = unique(j)
points( PC2 ~ PC1, pcadata[j,], pch=19, cex=2.0, col="lightgreen" )
text( PC2 ~ PC1, labels=vern, data=pcadata[j,], cex=0.75, col="darkgreen"  )
```


```{r fig-predators2, echo=FALSE} 
kable( counts[1:11,], format="simple", row.names=FALSE)
```



 
##  Predators {.c}
  


```{r fig-predators3a, out.width='100%', echo=FALSE, fig.align='center', fig.cap = 'Location of predators with snow crab in stomach. 2000-2010', warning = FALSE, message = FALSE } 
ggplot() + 
    borders("world", fill = "lightgray", colour = "grey80") + 
    xlim(c( -65.5, -57.1)) + ylim(c(42.1, 47.1)) +
    geom_point(data=snowcrab_predators[year(timestamp) %in% c(2000:2010),], aes(x=slongdd, y=slatdd, colour=Predator ), size=2.5 ) +
    labs(x="", y="", caption="2000-2010") +
    theme(legend.position="inside", legend.position.inside =c(0.16, 0.65), legend.title=element_blank(), legend.text=element_text(size=7.0) 
) 
```

```{r fig-predators3b, out.width='100%', echo=FALSE, fig.align='center', fig.cap = 'Location of predators with snow crab in stomach. 2011-2020', warning = FALSE, message = FALSE } 
ggplot() + 
    borders("world", fill = "lightgray", colour = "grey80") + 
    xlim(c( -65.5, -57.1)) + ylim(c(42.1, 47.1)) +
    geom_point(data=snowcrab_predators[year(timestamp) %in% c(2011:2020),], aes(x=slongdd, y=slatdd, colour=Predator ), size=2.5 ) +
    labs(x="", y="", caption="2011-2020") +
    theme(legend.position="inside", legend.position.inside =c(0.16, 0.725), legend.title=element_blank(), legend.text=element_text(size=7.0) )
```




## Predators - Atlantic cod ... {.c}

```{r fig-cod-map, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Atlantic cod density log10(no/km$^2$) from the Snow Crab survey. Predation mortality more likely higher encounter rates (high densities co-occur).' }
loc = file.path( SCD, 'output', 'maps', 'survey', 'snowcrab', 'annual', 'bycatch', 'ms.no.10' )
yrsplot = setdiff( year_assessment + c(0:-9), 2020)
fn4 = file.path( loc, paste( 'ms.no.10', yrsplot[4], 'png', sep='.') )
fn3 = file.path( loc, paste( 'ms.no.10', yrsplot[3], 'png', sep='.') )
fn2 = file.path( loc, paste( 'ms.no.10', yrsplot[2], 'png', sep='.') )
fn1 = file.path( loc, paste( 'ms.no.10', yrsplot[1], 'png', sep='.') )
include_graphics( c(  fn3, fn2, fn1) )
# \@ref(fig:cod-map)  
```
 


## Predators - Atlantic cod {.c}

```{r fig-cod-timeseries, out.width='60%', echo=FALSE,  fig.align='center', fig.cap = 'Atlantic cod crude, unadjusted geometric mean numerical density (no/km$^2$) from annual Snow Crab survey. Error bars are 95\\%  Confidence Intervals.'}
include_graphics( file.path( SCD, 'assessments', year_assessment, 'timeseries', 'survey', 'ms.no.10.png') )
# \@ref(fig:cod-timeseries)
```


 
##  Predators - Atlantic Halibut  {.c}

```{r fig-halibut-map, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Halibut density log10(no/km$^2$) from the Snow Crab survey. Predation mortality more likely higher encounter rates (high densities co-occur).' }
loc = file.path( SCD, 'output', 'maps', 'survey', 'snowcrab', 'annual', 'bycatch', 'ms.no.30' )
yrsplot = setdiff( year_assessment + c(0:-9), 2020)
fn4 = file.path( loc, paste( 'ms.no.30', yrsplot[4], 'png', sep='.') )
fn3 = file.path( loc, paste( 'ms.no.30', yrsplot[3], 'png', sep='.') )
fn2 = file.path( loc, paste( 'ms.no.30', yrsplot[2], 'png', sep='.') )
fn1 = file.path( loc, paste( 'ms.no.30', yrsplot[1], 'png', sep='.') )
include_graphics( c(  fn3, fn2, fn1) )
# \@ref(fig:halibut-map)  
```
 

##  Predators - Atlantic Halibut ... {.c}

```{r fig-halibut-timeseries, out.width='50%', echo=FALSE,   fig.align='center', fig.cap = 'Atlantic Halibut crude, unadjusted geometric mean numerical density (no/km$^2$) from annual Snow Crab survey. Error bars are 95\\%  Confidence Intervals.' }
include_graphics( file.path( SCD, 'assessments', year_assessment, 'timeseries', 'survey', 'ms.no.30.png') )
# \@ref(fig:halibut-timeseries)
```


## Disease

```{r fig-bcd, out.width='80%', echo=FALSE,   fig.align='center', fig.cap = 'Bitter Crab Disease in Maritimes Region.' }
include_graphics( file.path( SCD, 'output', 'bcd.png') )
```
 

 


## Connectivity: Movement 
   
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

fns = file.path( media_loc, c(
  "movement0.png", 
  "movement.png" ,
  "snowcrab_movement_distances.png", 
  "snowcrab_movement_rates.png" 
) )

knitr::include_graphics( fns ) 

``` 

&nbsp;  $~$  <br /> 



# Life history {.c}

## Life history {.c}

```{r fig-photos, echo=FALSE, out.width='30%', fig.align='center', fig.show='hold',  fig.cap = 'Snow Crab pelagic Zoea, benthic male and mating pair. Note sexual dimorphism.' }
fn1=file.path( media_loc, "snowcrab_zoea.png" )
fn2=file.path( media_loc, "snowcrab_male.png" )
fn3=file.path( media_loc, "snowcrab_male_and_female.png" )
knitr::include_graphics( c(fn1, fn2, fn3) ) 
# \@ref(fig:photos)  
```

## Life history: stages{.c}
 
```{r fig-lifehistory, echo=FALSE, out.width='90%', fig.align='center', fig.cap = 'Life history patterns of snow crab and approximate timing of the main life history stages of snow crab and size (carapace width; CW mm) and instar (Roman numerals). Size and timings are specific to the area of study and vary with environmental conditions, food availability and genetic variability.' }
fn1=file.path( media_loc, "life_history.png" )
knitr::include_graphics( fn1 ) 

```

## Life history: male growth stanzas {.c}
 
```{r fig-lifehistory-male, echo=FALSE, out.width='50%', fig.align='center', fig.cap = 'The growth stanzas of the male component and decision paths to maturity and terminal moult. Black ellipses indicate terminally molted animals.' }
fn1=file.path( media_loc, "life_history_male.png" )
knitr::include_graphics( fn1 ) 

```

## Life history: growth modes{.c}

```{r fig-growth-modes, echo=FALSE, out.width='40%', fig.align='center', fig.cap = 'Modal analysis.' }
fn1=file.path( media_loc, "growth_summary.png" )
knitr::include_graphics( c(fn1) ) 

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
   

## Life history: movement (mark-recapture) {.c}
  
```{r fig-movementtracks33, echo=FALSE, out.width='60%', fig.align='center', fig.show='hold',  fig.cap = 'Snow Crab movement tracks from mark and recapture.' }
fn1=file.path( media_loc, "movement0.png" )
fn2=file.path( media_loc, "movement.png" )
knitr::include_graphics( c(fn1, fn2) ) 

```
 
```{r fig-movement33, echo=FALSE, out.width='55%', fig.align='center', fig.show='hold',  fig.cap = 'Snow Crab movement distances and rates.' }
fn1=file.path( media_loc, "snowcrab_movement_distances.png" )
fn2=file.path( media_loc, "snowcrab_movement_rates.png" )
knitr::include_graphics( c(fn1, fn2) ) 

``` 
 
 
## Life history: clustering  {.c}
 
```{r fig-aggregation, echo=FALSE, out.width='60%', fig.align='center', fig.show='hold',  fig.cap = 'Australian spider crab \\emph{Leptomithrax gaimardii} aggregation for moulting and migration and Alaska red king crab \\emph{Paralithodes camtschaticus} aggregation in Alaska for egg release, migrations.' }
fn1=file.path( media_loc, "australian_leptomithrax_gaimardii.png" )
fn2=file.path( media_loc, "kingcrab_aggregation.png" ) 
knitr::include_graphics( c(fn1, fn2 ) ) 

```

- *Other* crab species show "Mounding" for protection from predation (larval, moulting and females)

- Narrow habitat preferences force them to move and cluster when environment is poor


```{r fig-clustering, echo=FALSE, out.width='90%', fig.align='center', fig.show='hold',  fig.cap = 'High density locations of Snow Crab, approximately 1 per square meter.'}
fn = file.path(p$project.outputdir, "maps", "map_highdensity_locations.png" )
knitr::include_graphics( fn ) 

```


# Human interactions

## Human interactions: Every known population is exploited worldwide {.c}

- Human consumption (males; > 95mm CW) due to higher meat yield and sexual dimorphism
- Bait in other fisheries (illegal in Canada, but still present)
- Fertilizer in agriculture
- Chitosan: glucosamine polysaccharide derived from chitin
    - agent to slow bleeding from wounds (Zhang et al. 2015, Moghadas et al. 2016)
    - agriculturally as natural fungicides and bactericides (Linden & Stoner 2007)
    - plastic replacement (Tampieri et al. 2003) 
    - battery electrolyte (Poosapati et al. 2021).


  
## Human interactions: Management Approach
 
\begin{columns}
\begin{column}{.46\textwidth}
\begin{tiny}
```{r fig-area_map33, echo=FALSE, out.width='80%', fig.align='center', fig.cap = 'The Scotian Shelf (NW Atlantic Ocean; NAFO Div. 4VWX). Shown are isobaths and major bathymetric features. Managed Crab Fishing Areas (CFAs; divided by dashed lines) include: NENS, SENS, 4X. SENS is further subdivided (dotted line) into 23 (NW) and 24 (SE).' }
fn1=file.path( media_loc, "area_map.png" )
knitr::include_graphics( fn1 ) 
# \@ref(fig:area_map)  
```

\item Precautionary Approach, Fish Stock Provisions, 2022
\item Spatial refugia (slope edge, MPAs) 
\item Temporal refugia (fishing seasons) 
\item Biological refugia: most life stages protected 

  \item Conservative exploitation since mid-2000s 
  \item Spawning stock legally and completely protected 
  \item Market-driven protection for 10+ yrs  

\item Evidence-based decision making: trawl survey, assessment
\item Distributed knowledge network: traditional, historical, scientific  
\item Satellite VMS; biodegradeable mesh (ghost-fishing); weighted lines (entanglement), etc ... 

 

# Stock status  {.c}

## Stock status: survey  {.c}
 
- No survey in 2020 (Covid-19) and incomplete 2022 (mechanical issues). 

- Inshore areas of S-ENS were most affected. 

- N-ENS and CFA 4X were not affected.
  

## Stock status: Size structure 

Factors: early maturation, size-selective predation, fishing of largest individuals, start or end of a recruitment pulse, timing of survey. 


```{r fig-meansize-male-mat, out.width='50%', echo=FALSE, eval = FALSE, fig.align='center', fig.cap = 'Mean size of mature male Snow Crab (CW; mm) from surveys with 95\\% Confidence Intervals'}
include_graphics(  file.path( SCD, "assessments", year_assessment, "timeseries", "survey", "cw.mat.png" )  )

``` 

## Stock status: Recruitment
  
  \item Little to no recruitment is expected for the next 1-3 years in N-ENS.
  \item Moderate levels of recruitment are expected in S-ENS. 
  \item Low to moderate levels of recruitment are expected for 2 years in 4X.
  

## Stock status: Reproduction
   \item All areas had recruitment of female crab into the mature (egg-bearing) segment of the population from 2016-2022.
  \item In N-ENS for 2022, a decline in numerical densities, and low densities of adolescent females. 
  \item Egg and larval production is expected to be moderate to high in the next year in all areas except N-ENS.  



## Stock status: Reproduction ...

```{r fig-fmat-timeseries, out.width='50%', echo=FALSE, fig.align='center', fig.cap = 'Mature female density log$_{10}$(no/km$^2$) from the Snow Crab survey.'  }
include_graphics( file.path( SCD, "assessments", year_assessment, "timeseries", "survey", "totno.female.mat.png") )

```
 

## Stock status: Mature female 

Distributions are heterogeneous and often in shallower areas.  

```{r fig-fmat-map, echo=FALSE, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Mature female density log$_{10}$(no/km$^2$)  from the Snow Crab survey.' }
loc = file.path( SCD, "output", "maps", "survey", "snowcrab", "annual", "totno.female.mat" )
yrsplot = setdiff( year_assessment + c(0:-9), 2020)
fn4 = file.path( loc, paste( "totno.female.mat", yrsplot[4], "png", sep=".") )
fn3 = file.path( loc, paste( "totno.female.mat", yrsplot[3], "png", sep=".") )
fn2 = file.path( loc, paste( "totno.female.mat", yrsplot[2], "png", sep=".") )
fn1 = file.path( loc, paste( "totno.female.mat", yrsplot[1], "png", sep=".") )
include_graphics( c( fn3, fn2, fn1) )

```


## Stock status: Viable Habitat {.t}
 ```{r fig-habitat, echo=FALSE, out.width='52%', fig.align='center', fig.show='hold',  fig.cap = 'Persistent habitat.' }
fn1=file.path( media_loc, "viable_habitat.png" ) 
knitr::include_graphics( c(fn1  ) ) 

```

```{r fig-habitat2, echo=FALSE, out.width='40%', fig.align='center', fig.show='hold',  fig.cap = 'Habitat preferences.' }
fn2=file.path( media_loc, "viable_habitat_depth_temp.png" )
knitr::include_graphics( c( fn2 ) ) 

```
\end{column}
\end{columns}
\end{small}


<!-- 
## Stock status: Viable Habitat ... {.t}
 
\vspace{6mm}

```{r habitat3, echo=FALSE, eval = FALSE, out.width='30%', fig.align='center', fig.show='hold',  fig.cap = 'Persistent habitat.' }
fn1=file.path( media_loc, "vh_i.png" ) 
fn2=file.path( media_loc, "vh_mf.png" ) 
fn3=file.path( media_loc, "vh_mm.png" ) 
knitr::include_graphics( c(fn1, fn2, fn3  ) ) 
# \@ref(fig:habitat3)  
```  
-->  

## Stock status: Viable Habitat ... {.c}
  
```{r fig-fb-habitat-map, out.width='30%', fig.show='hold', fig.align='center', fig.cap= 'Habitat viability (probability; fishable Snow Crab).' }

loc = file.path( SCD, 'modelled', 'default_fb', 'predicted_habitat' )
vn = "habitat."
yrsplot =  year_assessment + c(0:-10)
fn10 = file.path( loc, paste( vn, yrsplot[10], '.png', sep='') )
fn9 = file.path( loc, paste( vn, yrsplot[9], '.png', sep='') )
fn8 = file.path( loc, paste( vn, yrsplot[8], '.png', sep='') )
fn7 = file.path( loc, paste( vn, yrsplot[7], '.png', sep='') )
fn6 = file.path( loc, paste( vn, yrsplot[6], '.png', sep='') )
fn5 = file.path( loc, paste( vn, yrsplot[5], '.png', sep='') )
fn4 = file.path( loc, paste( vn, yrsplot[4], '.png', sep='') )
fn3 = file.path( loc, paste( vn, yrsplot[3], '.png', sep='') )
fn2 = file.path( loc, paste( vn, yrsplot[2], '.png', sep='') )
fn1 = file.path( loc, paste( vn, yrsplot[1], '.png', sep='') )
include_graphics( c(  fn3, fn2, fn1) )

```


## Stock status: Viable Habitat ... {.c}

```{r fig-fb-habitat-timeseries, out.width='50%', echo=FALSE,   fig.align='center', fig.cap = 'Habitat viability (probability; fishable Snow Crab). Means and 95\\% Credible Intervals are presented.' }
loc = file.path( SCD, 'modelled', 'default_fb', 'aggregated_habitat_timeseries' )
include_graphics( file.path( loc, 'habitat_M0.png') )

```
 

  
 
## Stock status: Sex ratios (proportion female, mature)

\item Mostly male-dominated: larger size may be protective against predation?
\item Imbalance indicates differential mortality: predation, competition and fishing
\item In 4X, sex ratios are balanced.  

```{r fig-sexratio-mature, out.width='60%', echo=FALSE, fig.align='center', fig.cap = 'Timeseries of sex ratios.'  }
include_graphics( file.path( SCD, "assessments", year_assessment, "timeseries", "survey", "sexratio.mat.png") )
# \@ref(fig:sexratio-mature)
```


```{r fig-sexratio-map, echo=FALSE, out.width='25%', fig.show='hold', fig.align='center', fig.cap= 'Map of sex ratios.'}
yrsplot = setdiff( year_assessment + c(0:-4), 2020)
loc = file.path( SCD, "output", "maps", "survey", "snowcrab", "annual", "sexratio.mat" )
fn4 = file.path( loc, paste( "sexratio.mat", yrsplot[4], "png", sep=".") )
fn3 = file.path( loc, paste( "sexratio.mat", yrsplot[3], "png", sep=".") )
fn2 = file.path( loc, paste( "sexratio.mat", yrsplot[2], "png", sep=".") )
fn1 = file.path( loc, paste( "sexratio.mat", yrsplot[1], "png", sep=".") )
include_graphics( c( fn3, fn2, fn1 ) )
# \@ref(fig:sexratio-map) 
```
 
 

## Stock status: Biomass Density  {.c}
  

Note that high and low biomass density areas fluctuate with time 


```{r fig-fbgeomean-map, echo=FALSE, out.width='30%', fig.show='hold', fig.align='center', fig.cap= 'Snow Crab survey fishable component biomass density log~10(t/km$^2$). Note, there is no data in 2020.' }
loc = file.path( SCD, 'output', 'maps', 'survey', 'snowcrab', 'annual', 'R0.mass')
yrsplot =  setdiff(year_assessment + c(0:-9), 2020 ) 
fn6 = file.path( loc, paste( 'R0.mass', yrsplot[6], 'png', sep='.') )
fn5 = file.path( loc, paste( 'R0.mass', yrsplot[5], 'png', sep='.') )
fn4 = file.path( loc, paste( 'R0.mass', yrsplot[4], 'png', sep='.') )
fn3 = file.path( loc, paste( 'R0.mass', yrsplot[3], 'png', sep='.') )
fn2 = file.path( loc, paste( 'R0.mass', yrsplot[2], 'png', sep='.') )
fn1 = file.path( loc, paste( 'R0.mass', yrsplot[1], 'png', sep='.') )
include_graphics( c(  fn3, fn2, fn1) )
``` 



## Stock status: Biomass Density ... {.c}

- A peak in 2009 to 2014 and has since been declining in all areas. 
 



```{r fig-fbGMTS, out.width='50%', echo=FALSE, fig.align='center', fig.cap = 'The crude, unadjusted geometric mean fishable biomass density log~10(t/km$^2$) from the Snow Crab survey. Error bars represent 95\\% Confidence Intervals. Note the absence of data in 2020. Prior to 2004, surveys were conducted in the Spring.'}
fn = file.path(SCD,'assessments', year_assessment, 'timeseries','survey','R0.mass.png')
include_graphics( c(fn) )

```




## Stock status: Biomass Index (aggregate)
  
A contraction of spatial range in 4X and the western parts of S-ENS were also evident in 2021 to 2022. 
 
```{r fig-fbindex-map, echo=FALSE, out.width='30%', fig.show='hold', fig.align='center', fig.cap= 'Biomass index log~10(t/km$^2$) predicted from the Snow Crab survey.' }
loc = file.path( SCD, 'modelled', 'default_fb', 'predicted_biomass_densities' )
yrsplot =  year_assessment + c(0:-10)
fn10 = file.path( loc, paste( 'biomass', yrsplot[10], 'png', sep='.') )
fn9 = file.path( loc, paste( 'biomass', yrsplot[9], 'png', sep='.') )
fn8 = file.path( loc, paste( 'biomass', yrsplot[8], 'png', sep='.') )
fn7 = file.path( loc, paste( 'biomass', yrsplot[7], 'png', sep='.') )
fn6 = file.path( loc, paste( 'biomass', yrsplot[6], 'png', sep='.') )
fn5 = file.path( loc, paste( 'biomass', yrsplot[5], 'png', sep='.') )
fn4 = file.path( loc, paste( 'biomass', yrsplot[4], 'png', sep='.') )
fn3 = file.path( loc, paste( 'biomass', yrsplot[3], 'png', sep='.') )
fn2 = file.path( loc, paste( 'biomass', yrsplot[2], 'png', sep='.') )
fn1 = file.path( loc, paste( 'biomass', yrsplot[1], 'png', sep='.') )
include_graphics( c( fn3, fn2, fn1) )
```
 

## Stock status: Biomass Index (aggregate) ... {.c}
     
```{r fig-fbindex-timeseries, out.width='60%', echo=FALSE, fig.align='center', fig.cap = 'The fishable biomass index (t) predicted by CARSTM of Snow Crab survey densities. Error bars represent Bayesian 95\\% Credible Intervals. Note large errors in 2020 when there was no survey.' }
include_graphics( file.path( SCD, 'modelled', 'default_fb', 'aggregated_biomass_timeseries' , 'biomass_M0.png') )

```



<!-- 

  ## Stock status: Modelled Biomass (pre-fishery) {.c}
  
  N-ENS: #r round(B_north[t0], 2)` t in #r year_assessment`

    - #r round(B_north[t1], 2)` t in #r year_previous`. 

  S-ENS: #r round(B_south[t0], 2)` t in #r year_assessment`

    - #r round(B_south[t1], 2)` t in #r year_previous`. 

  4X:  #r round(B_4x[t0], 2)` t in #r year_assessment`-#r year_assessment+1`

    - #r round(B_4x[t1], 2)` t for the #r year_previous`-#r year_assessment` season. 
 
-->


## Stock status: Modelled Biomass (pre-fishery) ... {.c}
 
```{r fig-logisticPredictions, out.width='32%', echo=FALSE, fig.show='hold', fig.align='center', fig.cap = 'Model 1 fishable, posterior mean modelled biomass (pre-fishery; kt) are shown in dark orange for N-ENS, S-ENS and 4X (left, middle and right). Light orange are posterior samples of modelled biomass (pre-fishery; kt) to illustrate the variability of the predictions. The biomass index (post-fishery, except prior to 2004) after model adjustment by the model catchability coefficient is in gray.' } 
loc = file.path( SCD, 'fishery_model', year_assessment, 'logistic_discrete_historical' )
fn1 = file.path( loc, 'plot_predictions_cfanorth.png' ) 
fn2 = file.path( loc, 'plot_predictions_cfasouth.png' ) 
fn3 = file.path( loc, 'plot_predictions_cfa4x.png' ) 
include_graphics(c(fn1, fn2, fn3) )

``` 
 
## Stock status: Fishing Mortality  {.c}

N-ENS: #r round(FM_north[t0],3)` (annual exploitation rate of #r round(100*(exp(FM_north[t0])-1),2)`%) in #r year_assessment`

  - Up from the  #r year_previous` rate of #r round(FM_north[t1],3)` (annual exploitation rate of #r round(100*(exp(FM_north[t1])-1),1)`%)
 
S-ENS: #r round(FM_south[t0],3)` (annual exploitation rate of #r round(100*(exp(FM_south[t0])-1),1)`%) in #r year_assessment`

  - Decreasing marginally from the #r year_previous` rate of #r round(FM_south[t1],3)` (annual exploitation rate of #r round(100*(exp(FM_south[t1])-1),1)`%)

4X: #r round(FM_4x[t0],3)` (annual exploitation rate of #r round(100*(exp(FM_4x[t0])-1),1)`%) in #r year_assessment`-#r year_assessment+1` season 

  - Decreasing from the #r year_assessment-1`-#r year_assessment` season rate of #r round(FM_4x[t1],3)` (annual exploitation rate of #r round(100*(exp(FM_4x[t1])-1),1)`%)

Localized exploitation rates are likely higher, as not all areas for which biomass is estimated are fished. 
 
 

## Stock status: Fishing Mortality ... {.c}

```{r fig-logisticFishingMortality, out.width='32%', echo=FALSE,  fig.show='hold', fig.align='center', fig.cap = 'Time-series of modelled instantaneous fishing mortality from Model 1, for N-ENS (left), S-ENS (middle), and 4X (right). Samples of the posterior densities are presented, with the darkest line being the mean.' }
  odir = file.path( fishery_model_results, year_assessment, "logistic_discrete_historical" )
  fn1 = file.path( odir, "plot_fishing_mortality_cfanorth.png" ) 
  fn2 = file.path( odir, "plot_fishing_mortality_cfasouth.png" ) 
  fn3 = file.path( odir, "plot_fishing_mortality_cfa4x.png" ) 
include_graphics(c(fn1, fn2, fn3) ) 
```



## Stock status: Reference Points {.c}
 

```{r fig-ReferencePoints, out.width='40%', echo=FALSE, fig.align='center', fig.cap = 'Harvest control rules for the Scotian Shelf Snow Crab fisheries.' }
include_graphics( file.path( media_loc, 'harvest_control_rules.png') ) 
# \@ref(fig:ReferencePoints)
```


## Stock status: Reference Points ... {.c}

```{r fig-logistic-hcr, out.width='32%', echo=FALSE, fig.show='hold', fig.align='center', fig.cap = 'Reference Points (fishing mortality and modelled biomass) from Model 1, for N-ENS (left), S-ENS (middle), and 4X (right). The large yellow dot indicates most recent year and the 95\\% CI.' }
  odir = file.path( fishery_model_results, year_assessment, "logistic_discrete_historical" )
  fn1 = file.path( odir, 'plot_hcr_cfanorth.png' ) 
  fn2 = file.path( odir, 'plot_hcr_cfasouth.png' ) 
  fn3 = file.path( odir, 'plot_hcr_cfa4x.png' ) 
  include_graphics(c(fn1, fn2, fn3) )
#  \@ref(fig:logistic-hcr)
```

  - N-ENS is in the "cautious" zone
  - S-ENS is in the "healthy" zone
  - 4X is in the "critical" zone


 
## Fishery Model Summary  {.c}

|   | N-ENS | S-ENS | 4X |  
|----- | ----- | ----- | ----- |
| |  |  |  |
|q       | `r round(q_north, 3)` (`r round(q_north_sd, 3)`) | `r round(q_south, 3)` (`r round(q_south_sd, 3)`) | `r round(q_4x, 3)` (`r round(q_4x_sd, 3)`) |
|r       | `r round(r_north, 3)` (`r round(r_north_sd, 3)`) | `r round(r_south, 3)` (`r round(r_south_sd, 3)`) | `r round(r_4x, 3)` (`r round(r_4x_sd, 3)`) |
|K       | `r round(K_north, 2)` (`r round(K_north_sd, 2)`) | `r round(K_south, 2)` (`r round(K_south_sd, 2)`) | `r round(K_4x, 2)` (`r round(K_4x_sd, 2)`) |
|Prefishery Biomass   | `r round(B_north[t0], 2)` (`r round(B_north_sd[t0], 2)`) | `r round(B_south[t0], 2)`  (`r round(B_south_sd[t0], 2)`) | `r round(B_4x[t0], 2)`  (`r round(B_4x_sd[t0], 2)`)  |
|Fishing Mortality    | `r round(FM_north[t0], 3)` (`r round(FM_north_sd[t0], 3)`) | `r round(FM_south[t0], 3)` (`r round(FM_south_sd[t0], 3)`) | `r round(FM_4x[t0], 3)` (`r round(FM_4x_sd[t0], 3)`) |


\tiny
Note: Values in parentheses are Posterior standard deviations.
\normalsize



## Reference Points {.c}
 

```{r fig-ReferencePoints, out.width='60%', echo=FALSE, fig.align='center', fig.cap = 'Harvest control rules for the Scotian Shelf Snow Crab fisheries.' }
include_graphics( file.path( media_loc, 'harvest_control_rules.png') ) 

```


## Reference Points ... {.c}

```{r fig-logistic-hcr, out.width='29%', echo=FALSE, fig.show='hold', fig.align='center', fig.cap = 'Reference Points (fishing mortality and modelled biomass) from the Fishery Model, for N-ENS (left), S-ENS (middle), and 4X (right). The large yellow dot indicates most recent year and the 95\\% CI. Not: the model does not account for illegal and unreported landings, and interspecific interactions.' }
  odir = file.path( fishery_model_results, year_assessment, "logistic_discrete_historical" )
  fn1 = file.path( odir, 'plot_hcr_cfanorth.png' ) 
  fn2 = file.path( odir, 'plot_hcr_cfasouth.png' ) 
  fn3 = file.path( odir, 'plot_hcr_cfa4x.png' ) 
  include_graphics(c(fn1, fn2, fn3) ) 
```
  - N-ENS is in the "healthy" zone
  - S-ENS is in the "healthy" zone
  - 4X is in the "critical" zone
  
 

## Conclusions

The ESS ecosystem is still experiencing a lot of volatility and prudence is wise:

- Rapid ecosystem change (groundfish collapse in 1980s and 1990s)
- Climatic variability: very warm period from 2012-2023
- Snow crab: 
  - Can persist if extreme conditions are episodic by shifting spatial distributions    
  - 4X was possibly too long and did not have the habitat space
  - Climate variability can induce juxtaposition of warmer water species (prey and predators) with snow crab
  

## Conclusions: N-ENS

- Viable habitat stable near historical mean and marginally improved relative to 2022.  
- Female larval production peaked in 2017 and has since declined.
- Recruitment continues at low levels, a noticible gap in recruitment to the fishery persists. 
- Predation mortality (Cod, Halibut, Skate) is likely high and causing reduced recruitment to the fishery. 
- Fishing mortality has been increasing since 2019. 
- Fishable biomass continues a decreasing trend.
- PA template suggest "healthy zone". 
- A more careful harvest strategy would enable bridging this coming recruitment gap. 
- A reduced TAC is prudent. 


## Conclusions: S-ENS

- Viable habitat marginally improved relative to a historic low in 2022.  
- Female larval production peaked in 2023/2024.
- Recruitment to the fishery from a strong year class has begun, full entry by 2025. 
- Predation mortality (Halibut, Plaice, Sculpin, Skate) is likely high and causing reduced recruitment to the fishery. 
- Fishing mortality has been increasing since 2019. 
- Fishable biomass is at a historical low.
- PA template suggest "healthy zone". 
- Flexiblity in harvest strategy exists due to strong recuitment. 
- A status quo TAC is prudent until confirmation of recruitment. 



## Conclusions: 4X

- Viable habitat has been at historical lows since 2015.  
- Female larval production peaked in 2022.
- Recruitment to the fishery in the near-term is not evident. There is a pulse 4-5 years away. 
- Predation mortality (Halibut) and competition with other Crustacea is likely high and causing reduced recruitment to the fishery. 
- Fishing mortality has declined to negligible levels since a peak in 2019. 
- Fishable biomass is at a historical low.
- PA template suggest "critical zone". 
- Stock status is extremely poor. 
- Fishery closure until recovery is prudent. 

 
 
## References
 
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

 

 

## Conclusions

The ESS ecosystem is still experiencing a lot of volatility and prudence is wise:

- Rapid ecosystem change 
- Rapid climatic change
- Persist under extreme conditions if they are episodic
- Shifts in spatial distribution towards cooler and deeper waters  

Modelled solutions:

- A few of many possible views
- Over-emphasis of any one view and associated Reference Points is **not precautionary**. 


## Conclusions ...

- In N-ENS, though recruitment continues at low levels, a gap in future recruitment to the fishery is expected for the next 1-3 years in N-ENS. N-ENS still exists in the healthy zone and so still has flexibility in approach. However, fishing mortality may have been higher than desirable. A more conservative harvest strategy would enable bridging this coming recruitment gap. A reduced TAC is prudent. 

- In S-ENS, recruitment to the fishery is likely to continue at a moderate rate for the upcoming season. The S-ENS stock remains in the healthy zone. Exploitation rates derived from the fishery model have been declining in recent years. S-ENS has considerable flexibility in approach. All TAC directions are viable. 

- In 4X, low to moderate levels of recruitment are expected for 2 years. 4X exists in the "cautious zone". The area is also in the southern-most extent of Snow Crab distribution in the North Atlantic and viable habitat has been depressed for many years. A reduced TAC is prudent. 
  

## References
 
Banerjee, S., Carlin, B. P., and Gelfand, A. E.. 2004. Hierarchical Modeling and Analysis for Spatial Data. Monographs on Statistics and Applied Probability. Chapman and Hall/CRC.

Besag, Julian. 1974. Spatial interaction and the statistical analysis of lattice systems. Journal of the Royal Statistical Society Series B (Methodological) 1974: 192-236.

Canada Gazette. 2022. [Regulations Amending the Fishery (General) Regulations. Part II, Volume 156, Number 8.](https://www.gazette.gc.ca/rp-pr/p2/2022/2022-04-13/html/sor-dors73-eng.html)

Canada Gazette. 2016. St. Anns Bank Marine Protected Area Regulations. Canada Gazette, Part I, Vol 150, Issue 51: 4143-4149.

Choi, J.S. 2020. A Framework for the assessment of Snow Crab (*Chioneocete opilio*) in Maritimes Region (NAFO Div 4VWX) . DFO Can. Sci. Advis. Sec. Res. Doc. 2020/nnn. v + xxx p.

Choi, J.S. 2022. Reconstructing the Decline of Atlantic Cod with the Help of Environmental Variability in the Scotian Shelf of Canada. bioRxiv. https://doi.org/10.1101/2022.05.05.490753.

Choi, J. S., and B. C. Patten. 2001. Sustainable Development: Lessons from the Paradox of Enrichment. Ecosystem Health 7: 163–77.

Choi, Jae S., B. Cameron, K. Christie, A. Glass, and E. MacEachern. 2022. Temperature and Depth Dependence of the Spatial Distribution of Snow Crab. bioRxiv. https://doi.org/10.1101/2022.12.20.520893.


DFO. 2013. [Integrated Fisheries Management Plan for Eastern Nova Scotia and 4X Snow Crab (*Chionoecetes Opillio*.)](http://www.dfo-mpo.gc.ca/fm-gp/peches-fisheries/ifmp-gmp/snow-crab-neige/snow-crab-neiges2013-eng.htm)
 

DFO. 2018. Stock Status Update of Atlantic Halibut (Hippoglossus hippoglossus) on the Scotian Shelf and Southern Grand Banks in NAFO Divisions 3NOPs4VWX5Zc. DFO Can. Sci. Advis. Sec. Sci. Resp. 2018/022.

Hebert M, Miron G, Moriyasu M, Vienneau R, and DeGrace P. Efficiency and ghost fishing of Snow Crab (Chionoecetes opilio) traps in the Gulf of St. Lawrence. Fish Res. 2001; 52(3): 143-153. 10.1016/S0165-7836(00)00259-9   
 
Riebler, A., Sørbye, S.H., Simpson D., and Rue, H. 2016. An intuitive Bayesian spatial model for disease mapping that accounts for scaling. Statistical methods in medical research 25: 1145-1165.

Simpson, D., Rue, H., Riebler, A., Martins, T.G., and Sørbye, SH. 2017. Penalising Model Component Complexity: A Principled, Practical Approach to Constructing Priors. Statist. Sci. 32: 1-28.

 
# END  
 


<!--

To do:

-->
