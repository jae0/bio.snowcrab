---
title: "Snow crab ecosystem and life history"
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
  - snow crab ecosystem assessment 
abstract: |
  Snow crab ecosystem assessment summary.
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

make quarto FN=04_ecosystem_summary YR=2024 SOURCE=~/bio/bio.snowcrab/inst/markdown WK=~/bio.data/bio.snowcrab/assessments DOCEXTENSION=html PARAMS="-P year_assessment:2024"


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



sn_env = snowcrab_load_key_results_to_memory( year_assessment, debugging=params$debugging,  return_as_list=TRUE  ) 

attach(sn_env)

# predator diet data
diet_data_dir = file.path( SCD, "data", "diets" )
require(data.table) # for speed
require(lubridate)
require(stringr) 
require(gt)  # table formatting
library(janitor)
require(ggplot2)
require(aegis) # map-related 
require(bio.taxonomy)  # handle species codes

# assimilate the CSV data tables:
# diet = get_feeding_data( diet_data_dir, redo=TRUE )  # if there is a data update
diet = get_feeding_data( diet_data_dir, redo=FALSE )
tx = taxa_to_code("snow crab")  
# matching codes are 
#  spec    tsn                  tx                   vern tx_index
#1  528 172379        BENTHODESMUS           BENTHODESMUS     1659
#2 2522  98427        CHIONOECETES SPIDER QUEEN SNOW UNID      728
#3 2526  98428 CHIONOECETES OPILIO        SNOW CRAB QUEEN      729
# 2 and 3 are correct

snowcrab_predators = diet[ preyspeccd %in% c(2522, 2526), ]  # n=159 oservations out of a total of 58287 observations in db (=0.28% of all data)
snowcrab_predators$Species = code_to_taxa(snowcrab_predators$spec)$vern
snowcrab_predators$Predator = factor(snowcrab_predators$Species)

counts = snowcrab_predators[ , .(Frequency=.N), by=.(Species)]
setorderv(counts, "Frequency", order=-1)

# species composition
psp = speciescomposition_parameters( yrs=p$yrs, carstm_model_label="default" )
pca = speciescomposition_db( DS="pca", p=psp )  

pcadata = as.data.frame( pca$loadings )
pcadata$vern = stringr::str_to_title( taxonomy.recode( from="spec", to="taxa", tolookup=rownames( pcadata ) )$vern )

```


&nbsp;  $~$  <br /> 


# Area of Interest: Scotian Shelf Ecosystem

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

# Life history 

## Photographs

```{r}
#| label: fig-photos-snowcrab
#| eval: true
#| echo: false 
#| output: true
#| fig-cap: "Snow Crab images. Note sexual dimorphism."
#| fig-subcap: 
#|   - "Pelagic zoea"
#|   - "Male mature benthic form"
#|   - "Mating pair - note sexual dimorphism, with a smaller female. Note epibionet growth on female which suggests it is an older female (multiparous)."
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold

fns = file.path( media_loc, c(
  "snowcrab_zoea.png", 
  "snowcrab_male.png" ,
  "snowcrab_male_and_female.png" 
) )

knitr::include_graphics( fns ) 

```

&nbsp;  $~$  <br /> 


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


## Size modes (carapce width)


![Modes idenfied from survey data for females in all data.](media/density_f_imodes.png){#fig-density_f_imodes}

![Modes idenfied from survey data for males in all data.](media/density_m_imodes.png){#fig-density_m_imodes}


## Growth patterns inferred

![Growth in female snow crab.](media/plot_growth_female.png){#fig-plot_growth_female}

![Growth in male snow crab.](media/plot_growth_male.png){#fig-plot_growth_male}


 

 
 
## Spatial clustering  
 
```{r}
#| label: fig-aggregation
#| eval: true
#| echo: false 
#| output: true
#| fig-cap: "Australian spider crab \\emph{Leptomithrax gaimardii} aggregation for moulting and migration and Alaska red king crab \\emph{Paralithodes camtschaticus} aggregation in Alaska for egg release, migrations."
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
 
fns = file.path( media_loc, c( 
  "australian_leptomithrax_gaimardii.png", 
  "kingcrab_aggregation.png" 
) )

knitr::include_graphics( fns ) 

```

- *Other* crab species show "Mounding" for protection from predation (larval, moulting and females)

- Narrow habitat preferences force them to move and cluster when environment is poor

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

knitr::include_graphics( fn ) 

```

## Summary

- Sexual dimorphism
- Ontogenetic shifts (changes with age and stage) in habitat preferences: 
  - Pelagic in larval stages, with diurnal veritical migration  
  - Benthic in pre-adolescent and adult stages
  - Complex substrates for small snow crab to muddy 
  - Adult females and smaller immature benthic crab are heat tolerant (up to 12 C of intermittant exposure) 
  - Adult males are stenothermic (narrow in range; less than 6 C is preferred)  
- Biannual, annual and biennial molts depending upon size/age and environmental conditions
- Terminal molt to maturity and survives for up to 5 years thereafter
- Total life span: up to 13 years (females) and 15 years (males)
- Cannibalism of immature crab by mature female snow crab are known 
- Fecundity:
  - Primiparous - 57.4 mm CW can produce between 35,000 to 46,000 eggs
  - Multiparous - more than 100,000 eggs  
- Eggs extruded February and April and brooded for up to two years 
  - more than 80% in the Area of Interest follow an annual cycle and hatched/released from April to June.
- Movement and connectivity
  - pelagic and zoea stages -- no information: vertical migration related to temperature and light
  - adult females -- no information: possibly related to food availability, temperature and predator avoidance 
  - adult males is long-tailed in distribution and possibly related to mate finding and food availablity


# Ecosystem

## Bathymetry

Data derived from multiple sources including Canadian Hydrographic Service, Snow crab surveys, Ground fish surveys, AZMP surveys. Modelled with [stmv](https://github.com/jae0/stmv). Predicted surface and derived variables are used as covariates for habitat and abundance modelling.


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

bdir = file.path(data_root, "aegis", "bathymetry", "modelled", "default", "stmv", "none_fft", "z", "maps", "canada.east.superhighres" )
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
#| fig-cap: "Substrate grain size log$_{10}$ (mm) variations (mean and standard deviation) in the Scotian Shelf region from a [carstm](https://github.com/jae0/carstm) analysis."
#| fig-subcap: 
#|   - "Mean prediction"
#|   - "Standard deviation"
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold 

substrdir = file.path( data_root, "aegis", "substrate", "modelled", "default", "maps" )
fns = file.path( substrdir, c(
  "predictions.png",
  "space_re_total.png"
) )

knitr::include_graphics( fns) 

```

## Rapid climate change

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

knitr::include_graphics( fn ) 

```

&nbsp;  $~$  <br /> 

 

## Ocean productivity: Chlorophyll-a (satellite-based)
 
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
knitr::include_graphics( fn ) 

```

&nbsp;  $~$  <br /> 


## Ocean currents

A visualization of [The Thermohaline Circulation - The Great Ocean Conveyor Belt](https://cdn.jwplayer.com/previews/IKH7eFQc). 


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
knitr::include_graphics( fn ) 

```

&nbsp;  $~$  <br /> 

### Strength of the Atlantic Meridional Overturning Circulation (AMOC):

```{r}
#| label: fig-ocean-currents-amoc-strength
#| eval: true
#| echo: false 
#| output: true
#| fig-cap: "Strength of the Atlantic Meridional Overturning Circulation (AMOC). Top is full time-eries. Bottom is the annual average. Source: [Copernicus Marine Service](https://marine.copernicus.eu/access-data/ocean-monitoring-indicators/atlantic-meridional-overturning-circulation-amoc-timeseries)"
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold 

fn = file.path( media_loc, c(
  "GLOBAL_OMI_NATLANTIC_amoc_max26N_timeseries-hq.png" 
) )
knitr::include_graphics( fn ) 

```

&nbsp;  $~$  <br /> 


### Expectations if the Atlantic Meridional Overturning Circulation (AMOC) collapses:


```{r}
#| label: fig-ocean-currents-amoc-expectation
#| eval: true
#| echo: false 
#| output: true
#| fig-cap: "Expectations if the Atlantic Meridional Overturning Circulation (AMOC) collapses. Source: [Icelandic Met Office](https://en.vedur.is/media/ads_in_header/AMOC-letter_Final.pdf)"
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold 

fn = file.path( media_loc, c(
  "amoc_letter.png" 
) )
knitr::include_graphics( fn ) 

```

&nbsp;  $~$  <br /> 



## Bottom Temperature 

```{r}
#| label: fig-bottom-temperatures-ts
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

tloc = file.path( SCD, "assessments", year_assessment, "timeseries"  )

fns = c( 
  file.path("survey", "t.png"), 
  "temperature_bottom.png" 
)

knitr::include_graphics( file.path( tloc, fns) )

```
 
Temporal variations in bottom temperature estimated from an analysis of historical temperature data using  [carstm](https://github.com/jae0/carstm). Red horizontal line is at $7^\circ$C. Presented are 95\% Credible Intervals of spatial variability in temperature at each time slice, after adjustment for spatiotemporal autocorrelation. Temperatures are variable and warm relative to many other snow crab regions, due to confluence of the warm, high salinity Gulf Stream from the S-SE along the shelf edge; cold, low salinity Labrador Current; and cold low salinity St. Lawrence outflow from the N-NE, as well as a nearshore Nova Scotia current, running from the NE. 
  

```{r}
#| label: fig-bottom-temperatures-map
#| eval: true
#| echo: false 
#| output: true
#| fig.show: hold  
#| fig-dpi: 144
#| fig-height: 4
#| fig-cap: "Spatial variations in bottom temperature estimated from a historical analysis of temperature data for 1 September for specified years."
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
  
loc = file.path( data_root, "aegis", "temperature", "modelled", "default", "maps" )
yrsplot =  year_assessment + c(0:-9)

fns  = file.path( loc, paste( "predictions.",  yrsplot,  ".0.75",  ".png", sep="") )
knitr::include_graphics( fns )
 

``` 

- Groundfish surveys were not conducted in 2020 and 2022 in the snow crab domain.
- Snow crab surveys were not conducted in 2020, and incomplete in 2022 for S-ENS.

  
 
&nbsp;  $~$  <br /> 
 

## Species composition

This is a Principle Component Analysis (eigen decomposition of relative abundance on standard deviation scale or "Z-score") with additional constraintof being positive valued (right shifted) andincluding zero-values. 



### Biplot or ordination

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

 
### Time series
```{r}
#| label: fig-speciescomposition-ts
#| eval: true
#| echo: false 
#| output: true
#| fig.show: hold 
#| fig-dpi: 144
#| fig-height: 6
#| fig-cap: "Species composition in time. Primary gradient (PC1) is related to bottom temperatures; second (PC2) to depth. Groundfish surveys were not conducted in 2020 and 2022 in the snow crab domain. Snow crab surveys were not conducted in 2020, and incomplete in 2022 in S-ENS."
#| fig-subcap: 
#|   - "Mean annual PC1 score."
#|   - "Mean annual PC2 score."
 
pc = c(1, 2)

spc_loc = file.path( data_root, "aegis", "speciescomposition", "modelled", "default" )
fnpr = file.path( spc_loc, "figures", paste("pca", pc, "_time.png", sep="" ) )
knitr::include_graphics( fnpr ) 

``` 

### Maps: PC1

```{r}
#| label: fig-speciescomposition-map-pc1
#| eval: true
#| echo: false 
#| output: true
#| fig.show: hold 
#| fig-dpi: 144
#| fig-height: 6
#| fig-cap: "Species composition (PC1) in space. Groundfish surveys were not conducted in 2020 and 2022 in the snow crab domain. Snow crab surveys were not conducted in 2020, and incomplete in 2022 in S-ENS."
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""
   
yrsplot =  year_assessment + c(0:-3)
 
vn = "pca1"

spc_loc = file.path( data_root, "aegis", "speciescomposition", "modelled", "default", "maps" )
fns = file.path( spc_loc, paste( vn, "predictions", yrsplot, "png", sep=".") )

knitr::include_graphics( fns ) 

``` 

### Maps: PC2

```{r}
#| label: fig-speciescomposition-map-pc2
#| eval: true
#| echo: false 
#| output: true
#| fig.show: hold 
#| fig-dpi: 144
#| fig-height: 6
#| fig-cap: "Species composition (PC2) in space. Groundfish surveys were not conducted in 2020 and 2022 in the snow crab domain. Snow crab surveys were not conducted in 2020, and incomplete in 2022 in S-ENS."
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""
   
yrsplot =  year_assessment + c(0:-3)
 
vn = "pca2"

spc_loc = file.path( data_root, "aegis", "speciescomposition", "modelled", "default", "maps" )
fns = file.path( spc_loc, paste( vn, "predictions", yrsplot, "png", sep=".") )

knitr::include_graphics( fns ) 

``` 

## Predators 


Most predators of snow crab are species associated with warmer water conditions (PC1 left of snow crab;  exception: Skates, Scuplins and Ocean Pout)  and deeper in distribution. 


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
lookup= c( "cod", "halibut", "sculpin", "skate", "plaice", "hake", "wolffish", "atlantic cod", "atlantic halibut", "longhorn sculpin", "thorny skate", "striped atlantic wolffish", "haddock", "american plaice", "smooth skate", "winter skate", "white hake", "shorthorn sculpin", "eelpout newfoundland", "squirrel or red hake", "sea raven", "ocean pout", "barndoor skate" )

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
#| label: tab-predators
#| echo: false
#| eval: true
#| output: true
#| tab-cap: "Main predators based upon frequency of occuence of snow crab in finfish stomach samples, unadjusted for sampling effort."

gt::gt(counts[1:11,]) 

```


  
## Prey

Prey items in the literature include:

- echinoderms
- polychaete worms (*Maldane*, *Nereis*), worm-like animals
- detritus (dead organic matter)
- large zooplankton, shrimp 
- juvenile crab (Rock Crab; Toad Crab; Lesser Toad Crab)
- Ocean Quahog (*Artica islandica*), bivalve molluscs (*Mytilus* sp, *Modiolus*, *Hiatella*)
- brittle stars (*Ophiura*, *Ophiopholis*)
- sea anemones (*Edwardsia*, *Metridium*). 

Most of the potental prey are found to the right of snow crab (i.e. colder-water species) at a variety of depths.


```{r}
#| label: fig-diet-biplot
#| eval: true
#| echo: false 
#| output: true
#| fig.show: hold 
#| fig-dpi: 144
#| fig-height: 6
#| fig-cap: "Relative location of snow crab prey (green) in the species composition ordination. Snow crab in orange."
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



##  Competitors

Some potential competition from other detritivorous speciies such as shrimp are possible. However, they tend to have slightly colder-water preferences, with the exception of *Pandalus borealis* (Northern Shrimp).

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
lookup= c( "pandalus" )   # add more here:

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




## Disease

[Bitter crab disease](https://www.dfo-mpo.gc.ca/science/aah-saa/diseases-maladies/hematcb-eng.html) is a dinoflagellate (*Hematodinium*) that causes muscle degeneration. They are widespread (Alaska, NW Atlantic, Greenland) and usually found in warm-water, physiologically stressful conditions. In th Maritimes, it seems to be a low level background infection, found everywhere in the fishing grounds.


```{r}
#| label: tab-bcd
#| echo: false
#| eval: true
#| output: true
#| tab-cap: "Bitter Crab Disease in Maritimes Region."
#| fig.show: hold 
#| fig-dpi: 144
#| fig-height: 10

include_graphics( file.path( SCD, "output", "bcd.png") )

```
 

```{r}
#| label: fig-bcd-map
#| eval: true
#| echo: false 
#| output: true
#| fig-cap: "Bitter crab disease observations since 2008"
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
  
fns = file.path( media_loc, c(
  "BCD_map.png"  
) )

knitr::include_graphics( fns ) 

``` 

&nbsp;  $~$  <br /> 



 


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
#| layout: [[100], [100], [50,50] ]
  
fns = file.path( media_loc, c(
  "movement0.png", 
  "movement.png" ,
  "snowcrab_movement_distances.png", 
  "snowcrab_movement_rates.png" 
) )

knitr::include_graphics( fns ) 

``` 

&nbsp;  $~$  <br /> 




## Ecosystem changes: summary 

- SSE: many species changes, snow crab are long-lived so interact with many of them
- Stomach samples: Atlantic Halibut, Atlantic Wolffish, Thorny Skate as primary predators
- Skewed sex ratios (few mature females): possibly linked to differential predation mortality (mature females being much smaller and full of fat-rich gonads and eggs). 
- Halibut (DFO 2018) 

  - Increased in abundance in the Region
  - Elevated co-ocurrence in areal **densities** with snow crab trawl samples means greater encounter rates




# Human interactions

## Every known population is exploited worldwide 

- Human consumption (males; > 95mm CW) due to higher meat yield and sexual dimorphism
- Bait in other fisheries (illegal in Canada, but still present)
- Fertilizer in agriculture
- Chitosan: glucosamine polysaccharide derived from chitin
    - agent to slow bleeding from wounds (Zhang et al. 2015, Moghadas et al. 2016)
    - agriculturally as natural fungicides and bactericides (Linden & Stoner 2007)
    - plastic replacement (Tampieri et al. 2003) 
    - battery electrolyte (Poosapati et al. 2021).


  
## Management Approach
 
- Precautionary Approach, Fish Stock Provisions, 2022
- Spatial refugia (slope edge, MPAs) 
- Temporal refugia (fishing seasons) 
- Biological refugia: most life stages protected 

  - Conservative exploitation since mid-2000s 
  - Spawning stock legally and completely protected 
  - Market-driven protection for 10+ yrs  

- Evidence-based decision making: trawl survey, assessment
- Distributed knowledge network: traditional, historical, scientific  
- Satellite VMS; biodegradeable mesh (ghost-fishing); weighted lines (entanglement), etc. 


# Snow crab fishery summary
 
See the [Fishery performance and status in fishery summary](02_fishery_summary.html) for more detailed information.


# Science survey  

Stock status of snow crab is based largely upon a directed trawl survey that is conducted annually in the autumn. See the [Survey result summary](02_survey_summary.html) for more information of collected data and basic statisitics.


# Synthesis

## Viable Habitat

### Temperature and depth

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
knitr::include_graphics( c( fn2 ) ) 

```

### Spatial random effects

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
knitr::include_graphics( c(fn1  ) ) 

```

 
### Predicted habitat maps

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
  
loc = file.path( SCD, "modelled", "default_fb", "predicted_habitat" )
vn = "habitat."
yrsplot =  year_assessment + c(0:-9)

fns = file.path( loc, paste( vn, yrsplot, ".png", sep="") )
include_graphics( fns )
 
```
 
### Predicted habitat timeseries


```{r}
#| label: fig-fb-habitat-timeseries
#| eval: true
#| echo: false 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "Habitat viability (probability; fishable Snow Crab). Means and 95\\% Credible Intervals are presented."

loc = file.path( SCD, "modelled", "default_fb", "aggregated_habitat_timeseries" )
include_graphics( file.path( loc, "habitat_M0.png") )

```
 

## Biomass density  

### Biomass density maps

Note that high and low biomass density areas fluctuate with time 
 
```{r}
#| label: fig-fbgeomean-map
#| eval: true
#| echo: false 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "Snow Crab survey fishable component biomass density log~10(t/km$^2$). Note, there is no data in 2020."
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
  

loc = file.path( SCD, "output", "maps", "survey", "snowcrab", "annual", "R0.mass")
yrsplot =  setdiff(year_assessment + c(0:-9), 2020 ) 

fns = file.path( loc, paste( "R0.mass", yrsplot, "png", sep=".") )
include_graphics( fns )

``` 



## Biomass density timeseries

- A peak in 2009 to 2014 and has since been declining in all areas. 
  

```{r}
#| label: fig-fbGMTS
#| eval: true
#| echo: false 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "The crude, unadjusted geometric mean fishable biomass density log~10(t/km$^2$) from the Snow Crab survey. Error bars represent 95\\% Confidence Intervals. Note the absence of data in 2020. Prior to 2004, surveys were conducted in the Spring."

fn = file.path(SCD, "assessments", year_assessment, "timeseries","survey","R0.mass.png")

include_graphics( fn )

```



## Biomass Index maps
  
A contraction of spatial range in 4X and the western parts of S-ENS were also evident in 2021 to 2022. 
 

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
  
loc = file.path( SCD, "modelled", "default_fb", "predicted_biomass_densities" )
yrsplot =  year_assessment + c(0:-9)

fns = file.path( loc, paste( "biomass", yrsplot, "png", sep=".") )
include_graphics( fns )
```


## Biomass Index (aggregate)  
   
```{r}
#| label: fig-fbindex-timeeries
#| eval: true
#| echo: false 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "The fishable biomass index (t) predicted by CARSTM of Snow Crab survey densities. Error bars represent Bayesian 95\\% Credible Intervals. Note large errors in 2020 when there was no survey."

fn = file.path( SCD, "modelled", "default_fb", "aggregated_biomass_timeseries" , "biomass_M0.png")
include_graphics( fn )

```


# Fishery modelled results and advise

To be presented at CSAS (Feb 2025).


# References

Banerjee, S., Carlin, B. P., and Gelfand, A. E.. 2004. Hierarchical Modeling and Analysis for Spatial Data. Monographs on Statistics and Applied Probability. Chapman and Hall/CRC.

Besag, Julian. 1974. Spatial interaction and the statistical analysis of lattice systems. Journal of the Royal Statistical Society Series B (Methodological) 1974: 192-236.

Canada Gazette. 2022. [Regulations Amending the Fishery (General) Regulations. Part II, Volume 156, Number 8.](https://www.gazette.gc.ca/rp-pr/p2/2022/2022-04-13/html/sor-dors73-eng.html)

Canada Gazette. 2016. St. Anns Bank Marine Protected Area Regulations. Canada Gazette, Part I, Vol 150, Issue 51: 4143-4149.

Choi, J.S. 2020. A Framework for the assessment of Snow Crab (*Chioneocete opilio*) in Maritimes Region (NAFO Div 4VWX) . DFO Can. Sci. Advis. Sec. Res. Doc. 2020/nnn. v + xxx p.

Choi, J.S. 2022. Reconstructing the Decline of Atlantic Cod with the Help of Environmental Variability in the Scotian Shelf of Canada. bioRxiv. [https://doi.org/10.1101/2022.05.05.490753](https://doi.org/10.1101/2022.05.05.490753].

Choi, J. S., and B. C. Patten. 2001. Sustainable Development: Lessons from the Paradox of Enrichment. Ecosystem Health 7: 163–77.

Choi, Jae S., B. Cameron, K. Christie, A. Glass, and E. MacEachern. 2022. Temperature and Depth Dependence of the Spatial Distribution of Snow Crab. bioRxiv. [https://doi.org/10.1101/2022.12.20.520893](https://doi.org/10.1101/2022.12.20.520893).


Choi, Jae S. 2023. A Multi-Stage, Delay Differential Model of Snow Crab Population Dynamics in the Scotian Shelf of Atlantic Canada. bioRxiv. [https://doi.org/10.1101/2023.02.13.528296](https://doi.org/10.1101/2023.02.13.528296).

DFO. 2013. [Integrated Fisheries Management Plan for Eastern Nova Scotia and 4X Snow Crab (*Chionoecetes Opillio*.)](http://www.dfo-mpo.gc.ca/fm-gp/peches-fisheries/ifmp-gmp/snow-crab-neige/snow-crab-neiges2013-eng.htm)
 
DFO. 2018. Stock Status Update of Atlantic Halibut (Hippoglossus hippoglossus) on the Scotian Shelf and Southern Grand Banks in NAFO Divisions 3NOPs4VWX5Zc. DFO Can. Sci. Advis. Sec. Sci. Resp. 2018/022.

Hebert M, Miron G, Moriyasu M, Vienneau R, and DeGrace P. Efficiency and ghost fishing of Snow Crab (Chionoecetes opilio) traps in the Gulf of St. Lawrence. Fish Res. 2001; 52(3): 143-153. 10.1016/S0165-7836(00)00259-9   
 
Riebler, A., Sørbye, S.H., Simpson D., and Rue, H. 2016. An intuitive Bayesian spatial model for disease mapping that accounts for scaling. Statistical methods in medical research 25: 1145-1165.

Simpson, D., Rue, H., Riebler, A., Martins, T.G., and Sørbye, SH. 2017. Penalising Model Component Complexity: A Principled, Practical Approach to Constructing Priors. Statist. Sci. 32: 1-28.

 
# END  
 

<!--

# Stock status

- Survey timing changed from Spring to Autumn in 2004.

- No survey in 2020 (Covid-19) and incomplete 2022 (mechanical issues).

  - Inshore areas of S-ENS were most affected.
  - N-ENS and CFA 4X were not affected.


## Size structure

Factors: early maturation, size-selective predation, fishing of largest individuals, start or end of a recruitment pulse, timing of survey.


```{r}
#| label: fig-meansize-male-mat
#| eval: true
#| echo: false
#| output: true
#| fig-cap: "Mean size of mature male Snow Crab (CW; mm) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
fn =  file.path( SCD, "assessments", year_assessment, "timeseries", "survey", "cw.male.mat.mean.png" )
knitr::include_graphics( fn )

```

## Recruitment

- Little to no recruitment is expected for the next 1-3 years in N-ENS.
- Moderate levels of recruitment are expected in S-ENS.
- Low to moderate levels of recruitment are expected for 2 years in 4X.


## Reproduction

- All areas had recruitment of female crab into the mature (egg-bearing) segment of the population from 2016-2022.
- In N-ENS for 2022, a decline in numerical densities, and low densities of adolescent females.
- Egg and larval production is expected to be moderate to high in the next year in all areas except N-ENS.



## Reproduction (mature females)


```{r}
#| label: fig-fmat-timeseries
#| eval: true
#| echo: false
#| output: true
#| fig-cap: "Mature female density log$_{10}$(no/km$^2$) from the Snow Crab survey."
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold

include_graphics( file.path( SCD, "assessments", year_assessment, "timeseries", "survey", "totno.female.mat.png") )

```

Distributions are heterogeneous and often in shallower areas.


```{r}
#| label: fig-fmat-map
#| eval: true
#| echo: false
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "Mature female density log$_{10}$(no/km$^2$)  from the Snow Crab survey."
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


loc = file.path( SCD, "output", "maps", "survey", "snowcrab", "annual", "totno.female.mat" )
yrsplot = setdiff( year_assessment + c(0:-8), 2020)

fns = file.path( loc, paste( "totno.female.mat", yrsplot, "png", sep=".") )
include_graphics( fns )

```


## Sex ratios (proportion female, mature)

- Mostly male-dominated: larger size may be protective against predation?
- Imbalance indicates differential mortality: predation, competition and fishing
- In 4X, sex ratios are balanced.

```{r}
#| label: fig-sexratio-mature-timeeries
#| eval: true
#| echo: false
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "Timeseries of sex ratios."

fn = file.path( SCD, "assessments", year_assessment, "timeseries", "survey", "sexratio.mat.png")
include_graphics( fn )

```


```{r}
#| label: fig-sexratio-mature-map
#| eval: true
#| echo: false
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "Map of sex ratios."
#| fig-subcap:
#|   - ""
#|   - ""
#|   - ""
#|   - ""
#| layout: [[100], [100], [50,50] ]

yrsplot = setdiff( year_assessment + c(0:-3), 2020)

loc = file.path( SCD, "output", "maps", "survey", "snowcrab", "annual", "sexratio.mat" )

fns = file.path( loc, paste( "sexratio.mat", yrsplot, "png", sep=".") )

include_graphics( fn )

```




## Modelled Biomass (pre-fishery)

N-ENS: `{r} round(B_north[t0], 2)` t in `{r} year_assessment`

- `{r} round(B_north[t1], 2)` t in `{r} year_previous`.

S-ENS: `{r} round(B_south[t0], 2)` t in `{r} year_assessment`

- `{r} round(B_south[t1], 2)` t in `{r} year_previous`.

4X:  `{r} round(B_4x[t0], 2)` t in `{r} year_assessment`-`{r} year_assessment+1`

- `{r} round(B_4x[t1], 2)` t for the `{r} year_previous`-`{r} year_assessment` season.




## Modelled Biomass (pre-fishery)

```{r}
#| label: fig-logisticPredictions
#| eval: true
#| echo: false
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "Fishable, posterior mean modelled biomass (pre-fishery; kt) are shown in dark orange. Light orange are posterior samples of modelled biomass (pre-fishery; kt) to illustrate the variability of the predictions. The biomass index (post-fishery, except prior to 2004) after model adjustment by the model catchability coefficient is in gray."
#| fig-subcap:
#|   - "N-ENS"
#|   - "S-ENS"
#|   - "4X"

loc = file.path( SCD, "fishery_model", year_assessment, "logistic_discrete_historical" )
fns = file.path( loc, c(
  "plot_predictions_cfanorth.png",
  "plot_predictions_cfasouth.png",
  "plot_predictions_cfa4x.png"
) )

include_graphics( fns )

```

## Fishing Mortality

N-ENS: `{r} round(FM_north[t0],3)` (annual exploitation rate of `{r} round(100*(exp(FM_north[t0])-1),2)`%) in `{r} year_assessment`

  - Up from the  `{r} year_previous` rate of `{r} round(FM_north[t1],3)` (annual exploitation rate of `{r} round(100*(exp(FM_north[t1])-1),1)`%)

S-ENS: `{r} round(FM_south[t0],3)` (annual exploitation rate of `{r} round(100*(exp(FM_south[t0])-1),1)`%) in `{r} year_assessment`

  - Decreasing marginally from the `{r} year_previous` rate of `{r} round(FM_south[t1],3)` (annual exploitation rate of `{r} round(100*(exp(FM_south[t1])-1),1)`%)

4X: `{r} round(FM_4x[t0],3)` (annual exploitation rate of `{r} round(100*(exp(FM_4x[t0])-1),1)`%) in `{r} year_assessment`-`{r} year_assessment+1` season

  - Decreasing from the `{r} year_assessment-1`-`{r} year_assessment` season rate of `{r} round(FM_4x[t1],3)` (annual exploitation rate of `{r} round(100*(exp(FM_4x[t1])-1),1)`%)

Localized exploitation rates are likely higher, as not all areas for which biomass is estimated are fished.



## Fishing Mortality

```{r}
#| label: fig-logisticFishingMortality
#| eval: true
#| echo: false
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "Time-series of modelled instantaneous fishing mortality from Model 1, for N-ENS (left), S-ENS (middle), and 4X (right). Samples of the posterior densities are presented, with the darkest line being the mean."
#| fig-subcap:
#|   - "N-ENS"
#|   - "S-ENS"
#|   - "4X"

odir = file.path( fishery_model_results, year_assessment, "logistic_discrete_historical" )
fns = file.path( odir, c(
  "plot_fishing_mortality_cfanorth.png",
  "plot_fishing_mortality_cfasouth.png",
  "plot_fishing_mortality_cfa4x.png"
))

include_graphics( fns )

```





## Reference Points

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


## Fishery Model Summary

|   | N-ENS | S-ENS | 4X |
|----- | ----- | ----- | ----- |
| |  |  |  |
|q       | `r round(q_north, 3)` (`r round(q_north_sd, 3)`) | `r round(q_south, 3)` (`r round(q_south_sd, 3)`) | `r round(q_4x, 3)` (`r round(q_4x_sd, 3)`) |
|r       | `r round(r_north, 3)` (`r round(r_north_sd, 3)`) | `r round(r_south, 3)` (`r round(r_south_sd, 3)`) | `r round(r_4x, 3)` (`r round(r_4x_sd, 3)`) |
|K       | `r round(K_north, 2)` (`r round(K_north_sd, 2)`) | `r round(K_south, 2)` (`r round(K_south_sd, 2)`) | `r round(K_4x, 2)` (`r round(K_4x_sd, 2)`) |
|Prefishery Biomass   | `r round(B_north[t0], 2)` (`r round(B_north_sd[t0], 2)`) | `r round(B_south[t0], 2)`  (`r round(B_south_sd[t0], 2)`) | `r round(B_4x[t0], 2)`  (`r round(B_4x_sd[t0], 2)`)  |
|Fishing Mortality    | `r round(FM_north[t0], 3)` (`r round(FM_north_sd[t0], 3)`) | `r round(FM_south[t0], 3)` (`r round(FM_south_sd[t0], 3)`) | `r round(FM_4x[t0], 3)` (`r round(FM_4x_sd[t0], 3)`) |

Note: Values in parentheses are Posterior standard deviations.



## Reference Points

```{r}
#| label: fig-logistic-hcr
#| eval: true
#| echo: false
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "Reference Points (fishing mortality and modelled biomass) from the Fishery Model, for N-ENS (left), S-ENS (middle), and 4X (right). The large yellow dot indicates most recent year and the 95\\% CI. Not: the model does not account for illegal and unreported landings, and interspecific interactions."
#| fig-subcap:
#|   - "N-ENS"
#|   - "S-ENS"
#|   - "4X"

odir = file.path( fishery_model_results, year_assessment, "logistic_discrete_historical" )

fns = file.path( odir, c(
  "plot_hcr_cfanorth.png" ,
  "plot_hcr_cfasouth.png",
  "plot_hcr_cfa4x.png"
) )

include_graphics( fns )

```


- N-ENS is in the "healthy" zone
- S-ENS is in the "healthy" zone
- 4X is in the "critical" zone



# Conclusions

The ESS ecosystem is still experiencing a lot of volatility and prudence is wise:

- Rapid ecosystem change (groundfish collapse in 1980s and 1990s)
- Rapid climatic change: variability high, very warm period from 2012-2024
- Snow crab:
  - Can persist if extreme conditions are episodic by shifting spatial distributions
  - 4X was possibly too long and did not have the habitat space
  - Climate variability can induce juxtaposition of warmer water species (prey and predators) with snow crab
  - Some shifts in spatial distribution towards cooler and deeper waters

Modelled solutions:

- A few of many possible views
- Over-emphasis of any one view and associated Reference Points is **not precautionary**.


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
In N-ENS, though recruitment continues at low levels,
a gap in future recruitment to the fishery is expected for the next 1-3 years


 bridging this coming recruitment gap. A reduced TAC is prudent.


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

recruitment to the fishery is likely to continue at a moderate rate for the upcoming season


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

   low to moderate levels of recruitment are expected for 2 years.
   4X exists in the "cautious zone".

   habitat has been depressed for many years. A reduced TAC is prudent.





-->

