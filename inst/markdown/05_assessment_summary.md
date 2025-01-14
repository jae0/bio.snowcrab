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




# Ecosystem

## Connectivity: Oceanic currents

```{r}
#| eval: true
#| echo: false 
#| output: true
#| label: ocean-currents
#| fig-cap: "Ocean currents in the Martimes. Source: https://www.dfo-mpo.gc.ca/oceans/publications/soto-rceo/2018/atlantic-ecosystems-ecosystemes-atlantiques/index-eng.html"
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold 

fn = file.path( media_loc, "maritimes_currents.png" )
knitr::include_graphics( fn ) 

```

&nbsp;  $~$  <br /> 




## Connectivity: Movement 
   
```{r}
#| eval: true
#| echo: false 
#| output: true
#| label: movement-tracks
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

## Bathymetry


```{r}
#| eval: true
#| echo: false 
#| output: true
#| label: bathymetry
#| fig-cap: "Substrate grain size log(mm) variations (mean and standard deviation) in the Scotian Shelf region."
#| fig-subcap: 
#|   - "Mean"
#|   - "Standard deviation"
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold 


```


## Substrate

```{r}
#| eval: true
#| echo: false 
#| output: true
#| label: substrate
#| fig-cap: "Substrate grain size log(mm) variations (mean and standard deviation) in the Scotian Shelf region."
#| fig-subcap: 
#|   - "Mean"
#|   - "Standard deviation"
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold 

substrdir = file.path( data_root, 'aegis', 'substrate', 'maps', 'canada.east.highres' )
fns = file.path( substrdir, c(
  "substrate.substrate.grainsize.canada.east.highres.png",
  "substrate.s.sdSpatial.canada.east.highres.png"
) )

knitr::include_graphics( fns) 

```


&nbsp;  $~$  <br /> 



## Bottom Temperature 

```{r}
#| eval: true
#| echo: false 
#| output: true
#| label: bottom-temperatures
#| fig-cap: "Bottom temperatures"
#| fig-subcap: 
#|   - "Annual variati#| ons in bottom temperature observed during the Snow Crab survey. The horizontal (black) line indicates the long-term, median temperature within each subarea. Error bars represent standard errors."
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




 
<!--

To do:

-->
