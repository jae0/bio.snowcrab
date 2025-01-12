---
title: "Snow Crab, Scotian Shelf, Canada (NAFO Div. 4VWX) in 2023"
subtitle: "General Summary"
author: "Jae S. Choi"
# footnote: "jae.choi@dfo-mpo.gc.ca"
institute: "Bedford Institute of Oceanography, DFO Science"
# date: "`r format(Sys.time(), '%d %B, %Y')`"
output:
  beamer_presentation:
    theme: "metropolis"
    colortheme: "seagull"
    fonttheme: "professionalfonts"
    fig_caption: yes
    # latex_engine: pdflatex
    latex_engine: lualatex 
    keep_tex: true
classoption: 
  - aspectratio=169 #16:9 wide
  - t  # top align
header-includes: 
  - \usepackage{graphicx}
  - \usepackage[font={scriptsize}, labelfont={bf}]{caption}
  # - \usepackage{float}
  # - \usepackage{subfig}
  # - \newcommand{\btiny}{\begin{tiny}}
  # - \newcommand{\etiny}{\end{tiny}}
params:
  year.assessment: 2023
  media_loc: "media"
  debugging: FALSE
--- 


<!-- Preamble
 

This is a Markdown document ... To create HTML or PDF, etc, run: 

  make rmarkdown FN=snowcrab_presentation_general_summary YR=2023 SOURCE=~/projects/bio.snowcrab/inst/markdown WK=~/bio.data/bio.snowcrab/assessments  DOCTYPE=beamer_presentation DOCEXTENSION=pdf   # {via Rmarkdown}

  --- note: columns only works with beamer_document


  make quarto FN=snowcrab_presentation_general_summary YR=2023 SOURCE=~/projects/bio.snowcrab/inst/markdown WK=~/bio.data/bio.snowcrab/assessments DOCEXTENSION=html  # {via Quarto}


  make pdf FN=snowcrab_presentation_general_summary  # {via pandoc}

Alter year and directories to reflect setup or copy Makefile and alter defaults to your needs.
  
 
 
-->




<!-- Set up R-environment -->
 

```{r setup, include=FALSE}
  require(knitr)
  knitr::opts_chunk$set(
    root.dir = data_root,
    echo = FALSE,
    out.width="6.2in",
#     dev.args = list(type = "cairo"),
    fig.retina = 2,
    dpi=192
  )

  # inits and data loading (front load all required data)

  require(aegis)
  
  year.assessment = params$year.assessment
  year_previous = year.assessment - 1
  p = bio.snowcrab::load.environment( year.assessment=year.assessment )
  SCD = project.datadirectory("bio.snowcrab")
  media_loc = params$media_loc
  
  # fishery_model_results = file.path( "/home", "jae", "projects", "dynamical_model", "snowcrab", "outputs" )
  fishery_model_results = file.path( SCD, "fishery_model" )

  sn_env = snowcrab_load_key_results_to_memory( year.assessment, debugging=params$debugging,  return_as_list=TRUE  ) 

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

  # bycatch summaries
  o_cfaall = observer.db( DS="bycatch_summary", p=p,  yrs=p$yrs, region="cfaall" )
  o_cfanorth = observer.db( DS="bycatch_summary", p=p,  yrs=p$yrs, region="cfanorth" )   
  o_cfasouth = observer.db( DS="bycatch_summary", p=p,  yrs=p$yrs, region="cfasouth" )   
  o_cfa4x = observer.db( DS="bycatch_summary", p=p,  yrs=p$yrs, region="cfa4x" )   

```
  

## Components of Snow Crab status in Maritimes Region {.c}

- Life history 
- Ecosystem change
- Human interactions 
  - Fishery statistics
  - Other interactions
- Scientific survey results
  - Life history X Ecosystem change X Human interactions = Complexity 



# Life history {.c}

## Life history {.c}

```{r photos, echo=FALSE, out.width='30%', fig.align='center', fig.show='hold',  fig.cap = 'Snow Crab pelagic Zoea, benthic male and mating pair. Note sexual dimorphism.' }
fn1=file.path( media_loc, "snowcrab_zoea.png" )
fn2=file.path( media_loc, "snowcrab_male.png" )
fn3=file.path( media_loc, "snowcrab_male_and_female.png" )
knitr::include_graphics( c(fn1, fn2, fn3) ) 
# \@ref(fig:photos)  
```

## Life history: stages{.c}
 
```{r lifehistory, echo=FALSE, out.width='90%', fig.align='center', fig.cap = 'Life history patterns of snow crab and approximate timing of the main life history stages of snow crab and size (carapace width; CW mm) and instar (Roman numerals). Size and timings are specific to the area of study and vary with environmental conditions, food availability and genetic variability.' }
loc = file.path( Sys.getenv("HOME"), "projects", "dynamical_model", "snowcrab", "media" )
fn1=file.path( loc, "life_history.png" )
knitr::include_graphics( fn1 ) 
# \@ref(fig:lifehistory)  
```

## Life history: male growth stanzas {.c}
 
```{r lifehistory_male, echo=FALSE, out.width='50%', fig.align='center', fig.cap = 'The growth stanzas of the male component and decision paths to maturity and terminal moult. Black ellipses indicate terminally molted animals.' }
loc = file.path( Sys.getenv("HOME"), "projects", "dynamical_model", "snowcrab", "media" )
fn1=file.path( loc, "life_history_male.png" )
knitr::include_graphics( fn1 ) 
# \@ref(fig:lifehistory_male)  
```

## Life history: growth modes{.c}

```{r growth_modes, echo=FALSE, out.width='40%', fig.align='center', fig.cap = 'Modal analysis.' }
loc = file.path( Sys.getenv("HOME"), "bio.data", "bio.snowcrab", "output" )
fn1=file.path( loc, "size_structure", "growth_summary.png" )
knitr::include_graphics( c(fn1) ) 
# \@ref(fig:lifehistory_male)  
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
 
\begin{small}
\begin{columns}
\begin{column}{.48\textwidth}
```{r movementtracks, echo=FALSE, out.width='60%', fig.align='center', fig.show='hold',  fig.cap = 'Snow Crab movement tracks from mark and recapture.' }
fn1=file.path( media_loc, "movement0.png" )
fn2=file.path( media_loc, "movement.png" )
knitr::include_graphics( c(fn1, fn2) ) 
# \@ref(fig:movementtracks)  
```
\end{column}

\begin{column}{.48\textwidth}
```{r movement, echo=FALSE, out.width='55%', fig.align='center', fig.show='hold',  fig.cap = 'Snow Crab movement distances and rates.' }
fn1=file.path( media_loc, "snowcrab_movement_distances.png" )
fn2=file.path( media_loc, "snowcrab_movement_rates.png" )
knitr::include_graphics( c(fn1, fn2) ) 
# \@ref(fig:movement)  
```
\end{column}
\end{columns}
\end{small}
 
 
## Life history: clustering  {.c}


\begin{small}
\begin{columns}
\begin{column}{.48\textwidth}
```{r aggregation, echo=FALSE, out.width='60%', fig.align='center', fig.show='hold',  fig.cap = 'Australian spider crab \\emph{Leptomithrax gaimardii} aggregation for moulting and migration and Alaska red king crab \\emph{Paralithodes camtschaticus} aggregation in Alaska for egg release, migrations.' }
fn1=file.path( media_loc, "australian_leptomithrax_gaimardii.png" )
fn2=file.path( media_loc, "kingcrab_aggregation.png" ) 
knitr::include_graphics( c(fn1, fn2 ) ) 
# \@ref(fig:aggregation)  
```

- *Other* crab species show "Mounding" for protection from predation (larval, moulting and females)

- Narrow habitat preferences force them to move and cluster when environment is poor
\end{column}

\begin{column}{.48\textwidth}
```{r clustering, echo=FALSE, out.width='90%', fig.align='center', fig.show='hold',  fig.cap = 'High density locations of Snow Crab, approximately 1 per square meter.'}
fn = file.path(p$project.outputdir, "maps", "map_highdensity_locations.png" )
knitr::include_graphics( fn ) 
# \@ref(fig:aggregation)  
if (0) {
    # high density locations directly from databases
    M = snowcrab.db( DS="set.complete", p=p ) 
    setDT(M)
    i = which(M$totno.all > 2.5*10^5)
    H = M[i, .( plon, plat, towquality, dist, distance, surfacearea, vessel, yr, z, julian, no.male.all, no.female.all, cw.mean, totno.all, totno.male.imm, totno.male.mat, totno.female.imm, totno.female.mat, totno.female.primiparous, totno.female.multiparous, totno.female.berried)]
    H$log10density = log10(H$totno.all)
    library(ggplot2)
    cst = coastline_db( p=p, project_to=st_crs(pg) ) 
    isodepths = c(100, 200, 300)
    isob = isobath_db( DS="isobath", depths=isodepths, project_to=st_crs(pg))
    isob$level = as.factor( isob$level)
    plt = ggplot() +
      geom_sf( data=cst, show.legend=FALSE ) +
      geom_sf( data=isob, aes( alpha=0.1, fill=level), lwd=0.1, show.legend=FALSE) +
    geom_point(data=H, aes(x=plon, y=plat, colour=log10density), size=5) +
      coord_sf(xlim = c(270, 940 ), ylim = c(4780, 5200 )) +
    theme(legend.position="inside", legend.position.inside=c(0.08, 0.8)) 
    png(filename=fn, width=1000,height=600, res=144)
      (plt)
    dev.off()
}
```
\begin{block}{Uncertainty}
Historical snow crab high density locations  
\end{block}
\end{column}
\end{columns}
\end{small}

# Ecosystem change

## Ecosystem change: Predators {.c}

- SSE: many species changes, snow crab are long-lived so interact with many of them
- Stomach samples: Atlantic Halibut, Atlantic Wolffish, Thorny Skate as primary predators
- Skewed sex ratios (few mature females): possibly linked to differential predation mortality (mature females being much smaller and full of fat-rich gonads and eggs). 
- Halibut (DFO 2018) have significantly increased in abundance in the Region. 
- Elevated co-ocurrence in areal **densities** with snow crab trawl samples means greater encounter rates


##  Ecosystem change: Predators - Atlantic Halibut  {.c}

```{r halibut-timeseries, out.width='50%', echo=FALSE,   fig.align='center', fig.cap = 'Atlantic Halibut crude, unadjusted geometric mean numerical density (no/km$^2$) from annual Snow Crab survey. Error bars are 95\\%  Confidence Intervals.' }
include_graphics( file.path( SCD, 'assessments', year.assessment, 'timeseries', 'survey', 'ms.no.30.png') )
# \@ref(fig:halibut-timeseries)
```
 
##  Ecosystem change: Predators - Atlantic Halibut ... {.c}

```{r halibut-map, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Halibut density log10(no/km$^2$) from the Snow Crab survey.' }
loc = file.path( SCD, 'output', 'maps', 'survey', 'snowcrab', 'annual', 'bycatch', 'ms.no.30' )
yrsplot = setdiff( year.assessment + c(0:-9), 2020)
fn4 = file.path( loc, paste( 'ms.no.30', yrsplot[4], 'png', sep='.') )
fn3 = file.path( loc, paste( 'ms.no.30', yrsplot[3], 'png', sep='.') )
fn2 = file.path( loc, paste( 'ms.no.30', yrsplot[2], 'png', sep='.') )
fn1 = file.path( loc, paste( 'ms.no.30', yrsplot[1], 'png', sep='.') )
include_graphics( c(  fn3, fn2, fn1) )
# \@ref(fig:halibut-map)  
```

\begin{block}{Uncertainty}
Higher predation mortality seems likely (more encounters with warmer-water species)
\end{block}



## Ecosystem change: Predators - Thorny skate {.c}

```{r thornyskate-timeseries, out.width='60%', echo=FALSE,  fig.align='center', fig.cap = 'Thorny Skate crude, unadjusted geometric mean numerical density (no/km$^2$) from annual Snow Crab survey. Error bars are 95\\%  Confidence Intervals.'}
include_graphics( file.path( SCD, 'assessments', year.assessment, 'timeseries', 'survey', 'ms.no.201.png') )
# \@ref(fig:thornyskate-timeseries)
```

## Ecosystem change: Predators - Thorny skate ... {.c}

```{r thornyskate-map, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Thorny skate density log10(no/km$^2$) from the Snow Crab survey.' }
loc = file.path( SCD, 'output', 'maps', 'survey', 'snowcrab', 'annual', 'bycatch', 'ms.no.201' )
yrsplot = setdiff( year.assessment + c(0:-9), 2020)
fn4 = file.path( loc, paste( 'ms.no.201', yrsplot[4], 'png', sep='.') )
fn3 = file.path( loc, paste( 'ms.no.201', yrsplot[3], 'png', sep='.') )
fn2 = file.path( loc, paste( 'ms.no.201', yrsplot[2], 'png', sep='.') )
fn1 = file.path( loc, paste( 'ms.no.201', yrsplot[1], 'png', sep='.') )
include_graphics( c(  fn3, fn2, fn1) )
# \@ref(fig:thornyskate-map)  
```


\begin{block}{Uncertainty}
Higher predation mortality seems likely (more encounters with warmer-water species)
\end{block}



<!--
~~~{=comment}

##  Ecosystem change: Predators - Striped Atlantic Wolffish  {.c}



```{r Wolffish-timeseries, out.width='60%', echo=FALSE,   fig.align='center', fig.cap = 'Striped Atlantic Wolffish crude, unadjusted geometric mean numerical density (no/km$^2$) from annual Snow Crab survey. Error bars are 95\\%  Confidence Intervals.' }
include_graphics( file.path( SCD, 'assessments', year.assessment, 'timeseries', 'survey', 'ms.no.50.png') )
# \@ref(fig:Wolffish-timeseries)
```

##  Ecosystem change: Predators - Striped Atlantic Wolffish  ... {.c}


```{r Wolffish-map, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Striped Atlantic Wolffish density log10(no/km$^2$) from the Snow Crab survey.' }
loc = file.path( SCD, 'output', 'maps', 'survey', 'snowcrab', 'annual', 'bycatch', 'ms.no.50' )
yrsplot = setdiff( year.assessment + c(0:-9), 2020)
fn4 = file.path( loc, paste( 'ms.no.50', yrsplot[4], 'png', sep='.') )
fn3 = file.path( loc, paste( 'ms.no.50', yrsplot[3], 'png', sep='.') )
fn2 = file.path( loc, paste( 'ms.no.50', yrsplot[2], 'png', sep='.') )
fn1 = file.path( loc, paste( 'ms.no.50', yrsplot[1], 'png', sep='.') )
include_graphics( c( fn3, fn2, fn1) )
# \@ref(fig:Wolffish-map)  
```


## Ecosystem change: Competitor -  Lesser toad crab    {.c}


```{r lessertoadcrab-timeseries, out.width='60%', echo=FALSE,   fig.align='center', fig.cap = 'Lesser Toad Crab crude, unadjusted geometric mean numerical density (no/km$^2$) from annual Snow Crab survey. Error bars are 95\\%  Confidence Intervals.' }
include_graphics( file.path( SCD, 'assessments', year.assessment, 'timeseries', 'survey', 'ms.no.2521.png') )
# \@ref(fig:lessertoadcrab-timeseries)
```
 

## Ecosystem change: Competitor - Lesser toad crab ... {.c}


```{r lessertoadcrab-map, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Lesser Toad Crab density log10(no/km$^2$) from the Snow Crab survey.' }
loc = file.path( SCD, 'output', 'maps', 'survey', 'snowcrab', 'annual', 'bycatch', 'ms.no.2521' )
yrsplot = setdiff( year.assessment + c(0:-9), 2020)
fn4 = file.path( loc, paste( 'ms.no.2521', yrsplot[4], 'png', sep='.') )
fn3 = file.path( loc, paste( 'ms.no.2521', yrsplot[3], 'png', sep='.') )
fn2 = file.path( loc, paste( 'ms.no.2521', yrsplot[2], 'png', sep='.') )
fn1 = file.path( loc, paste( 'ms.no.2521', yrsplot[1], 'png', sep='.') )
include_graphics( c( fn3, fn2, fn1) )
# \@ref(fig:lessertoadcrab-map)  
```

~~~
-->  


## Ecosystem change: Co-occurring - Northern shrimp  {.c}


```{r Shrimp-timeseries, out.width='60%', echo=FALSE,   fig.align='center', fig.cap = 'Northern Shrimp crude, unadjusted geometric mean numerical density (n/$km^2$) from annual Snow Crab survey. Error bars are 95\\%  Confidence Intervals.' }
include_graphics( file.path( SCD, 'assessments', year.assessment, 'timeseries', 'survey', 'ms.no.2211.png') )
# \@ref(fig:Shrimp-timeseries)
```
 
## Ecosystem change: Co-occurring - Northern shrimp ... {.c}

```{r Shrimp-map, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Northern Shrimp density log10(no/km$^2$) from the Snow Crab survey.' }
loc = file.path( SCD, 'output', 'maps', 'survey', 'snowcrab', 'annual', 'bycatch', 'ms.no.2211' )
yrsplot = setdiff( year.assessment + c(0:-9), 2020)
fn4 = file.path( loc, paste( 'ms.no.2211', yrsplot[4], 'png', sep='.') )
fn3 = file.path( loc, paste( 'ms.no.2211', yrsplot[3], 'png', sep='.') )
fn2 = file.path( loc, paste( 'ms.no.2211', yrsplot[2], 'png', sep='.') )
fn1 = file.path( loc, paste( 'ms.no.2211', yrsplot[1], 'png', sep='.') )
include_graphics( c( fn3, fn2, fn1) )
# \@ref(fig:Shrimp-map)  
```

Shrimp with similar habitat preferences have declined, possibly due to large-scaled habitat variations and predation.

\begin{block}{Uncertainty}
Sampling was incomplete in 2020 and 2022 in S-ENS.
\end{block}

 


## Ecosystem change: species composition 

```{r speciesomposition, echo=FALSE, out.width='48%', fig.align='center', fig.show='hold',  fig.cap = 'Species composition in space and time. Primary gradient is related to bottom temperatures.' }
spc_loc = file.path( data_root, 'aegis', 'speciescomposition', 'modelled', 'default', 'maps' )
fn1 = file.path( spc_loc, 'pca1.space_re_total.png') 
fn2 = file.path( spc_loc, 'pca2.space_re_total.png')
ts_loc = file.path( data_root, 'aegis', 'speciescomposition', 'modelled', 'default', 'figures' )
fn3 = file.path( ts_loc, 'pca1_time.png') 
fn4 = file.path( ts_loc, 'pca2_time.png') 
knitr::include_graphics( c(fn1,  fn3  ) ) 
# \@ref(fig:habitat3)  
``` 
 

\begin{block}{Uncertainty}
Sampling was incomplete in 2020 and 2022 in S-ENS.
\end{block}


## Ecosystem change: Bottom Temperature {.c}

```{r bottom-temperatures-survey, out.width='50%', echo=FALSE, fig.align='center', fig.cap = 'Annual variations in bottom temperature observed during the Snow Crab survey. The horizontal (black) line indicates the long-term, median temperature within each subarea. Error bars represent standard errors.' }
knitr::include_graphics( file.path( SCD, 'assessments', year.assessment, 'timeseries', 'survey', 't.png') )
# \@ref(fig:bottom-temperatures-survey)
```

\begin{block}{Uncertainty}
Sampling was incomplete in 2020 and 2022 in S-ENS.
\end{block}


## Ecosystem change: Bottom Temperature ... {.c}

- Average bottom temperatures **observed** in the 2022 Snow Crab survey were near or above historical highs in all areas 

- Temperatures are more stable in N-ENS than S-ENS; 4X exhibits the most erratic and highest annual mean bottom temperatures. 

- Observed temperatures in the 2022 Snow Crab survey for S-ENS increased well above the average. 
 
- Average temperature increased well beyond the $7^\circ$C threshold in 4X. N-ENS and S-ENS also continued to experience historical highs in bottom temperature and elevated spatial variability of bottom temperatures.



## Ecosystem considerations: Bottom Temperature ... {.c}

```{r bottom-temperatures-map, out.width='30%', echo=FALSE, fig.show='hold', fig.align='center', fig.cap = 'Spatial variations in bottom temperature estimated from a historical analysis of temperature data for 1 September.' }

loc = file.path( data_root, 'aegis', 'temperature', 'modelled', 'default', 'maps' )
yrsplot =  year.assessment + c(0:-10)
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


\begin{block}{Uncertainty}
* Groundfish surveys were not conducted in 2020 and 2022 in the snow crab domain.

* Snow crab surveys were not conducted in 2020, and incomplete in 2022 for S-ENS.
\end{block}


 
## Ecosystem considerations: Bottom Temperature ... {.c}

Persistent spatial gradient of almost $15^\circ$C in bottom temperatures in the Maritimes Region. 

Variable due to confluence of the warm, high salinity Gulf Stream from the S-SE along the shelf edge; cold, low salinity Labrador Current; and cold low salinity St. Lawrence outflow from the N-NE, as well as a nearshore Nova Scotia current, running from the NE. 

```{r bottom-temperatures-spatialeffect, out.width='35%', echo=FALSE, fig.align='center', fig.cap = 'Persistent spatial effect of bottom temperature, relative to the overall mean, after adjustment for spatiotemporal variability and autocorrelations. Time period from 1999 to present.' }
loc = file.path( data_root, 'aegis', 'temperature', 'modelled', 'default', 'maps' )
knitr::include_graphics( file.path( loc, 'space_re_total.png') )
# \@ref(fig:bottom-temperatures-spatialeffect)
```
 
## Ecosystem considerations: Bottom Temperature ... {.c}
 
\begin{columns}
\begin{column}{.6\textwidth}
```{r bottom-temperatures, out.width='65%', echo=FALSE, fig.align='center', fig.cap = '' }
knitr::include_graphics( file.path( SCD, 'assessments', year.assessment, 'timeseries', 'temperature_bottom.png') )
# \@ref(fig:bottom-temperatures)
```
\end{column}
\begin{column}{.4\textwidth}

\vspace{12mm}

\begin{footnotesize}
\textbf{Figure}: Temporal variations in bottom temperature estimated from a historical analysis of temperature data. Red horizontal line is at $7^\circ$C. Presented are 95\% Credible Intervals of spatial variability in temperature at each time slice, after adjustment for spatiotemporal autocorrelation.
\end{footnotesize}
\end{column}
\end{columns}

 
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
```{r area_map, echo=FALSE, out.width='80%', fig.align='center', fig.cap = 'The Scotian Shelf (NW Atlantic Ocean; NAFO Div. 4VWX). Shown are isobaths and major bathymetric features. Managed Crab Fishing Areas (CFAs; divided by dashed lines) include: NENS, SENS, 4X. SENS is further subdivided (dotted line) into 23 (NW) and 24 (SE).' }
loc = file.path( Sys.getenv("HOME"), "projects", "dynamical_model", "snowcrab", "media" )
fn1=file.path( loc, "area_map.png" )
knitr::include_graphics( fn1 ) 
# \@ref(fig:area_map)  
```
\end{tiny}
\end{column}
\begin{column}{.52\textwidth}
\begin{footnotesize}
\begin{itemize}
\item Precautionary Approach, Fish Stock Provisions, 2022
\item Spatial refugia (slope edge, MPAs) 
\item Temporal refugia (fishing seasons) 
\item Biological refugia: most life stages protected 
\begin{itemize}
\begin{scriptsize}
  \item Conservative exploitation since mid-2000s 
  \item Spawning stock legally and completely protected 
  \item Market-driven protection for 10+ yrs  
\end{scriptsize}
\end{itemize}
\item Evidence-based decision making: trawl survey, assessment
\item Distributed knowledge network: traditional, historical, scientific  
\item Satellite VMS; biodegradeable mesh (ghost-fishing); weighted lines (entanglement), etc ... 
\end{itemize}
\end{footnotesize}
\end{column}
\end{columns}
 

## Human interactions: Fishing effort {.c}

Similar between `r year.assessment` and `r year_previous` in terms of spatial distribution. In S-ENS, there was, however, a minor spatial contraction to inshore areas and away from the area 23-24 boundary. 
 
\begin{tiny}
```{r effort-map, echo=FALSE, out.width='45%', fig.show='hold',  fig.align='center', fig.cap = 'Snow Crab fishing effort from fisheries logbook data for previous and current years. Units are No. $\\times 10^3$ per (10 km X 10 km) grid.' }
loc0= file.path( SCD, "output", "maps", "logbook", "snowcrab", "annual", "effort" )
fn1 = file.path( loc0, paste( "effort", year_previous,   "png", sep=".") ) 
fn2 = file.path( loc0, paste( "effort", year.assessment, "png", sep=".") ) 
include_graphics(  c(fn1, fn2) )
#  \@ref(fig:landings-map) 
```
\end{tiny}
 

## Human interactions: Fishing effort ...  {.c}
 
\begin{tiny}
```{r effort-timeseries, echo=FALSE, out.width='60%', fig.align='center', fig.cap = 'Temporal variations in fishing effort $\\times 10^3$ trap hauls.' }
fn1=file.path( SCD, "assessments", year.assessment, "timeseries", "fishery",   "effort.ts.png" )
knitr::include_graphics( fn1 ) 
# \@ref(fig:effort-timeseries)  
```
\end{tiny}

 
<!--
 
  ```{r table-fishery-nens, echo=FALSE, eval = FALSE }
  ii = which(dt$Region=="cfanorth")
  oo = dt[ii, c("Year", "Licenses", "TAC", "Landings", "Effort", "CPUE")] 
  kable( oo, format="simple", row.names=FALSE, align="cccccc",
  caption = "Fishery performance statistics in N-ENS. Units are: TACs and Landings (tons, t), Effort ($\\times 10^3$ trap hauls, th) and CPUE (kg/th).")
  # \@ref(tab:table-fishery-nens)
  ``` 

  ---
  
  ```{r table-fishery-sens, echo=FALSE, eval = FALSE }
  ii = which(dt$Region=="cfasouth")
  oo = dt[ii, c("Year", "Licenses", "TAC", "Landings", "Effort", "CPUE")] 
  kable( oo, format="simple", row.names=FALSE, align="cccccc",
  caption = "Fishery performance statistics in S-ENS. Units are: TACs and Landings (tons, t), Effort ($\\times 10^3$ trap hauls, th) and CPUE (kg/th).")
  # \@ref(tab:table-fishery-nens)
  ``` 

  ---
  
  ```{r table-fishery-4x, echo=FALSE, eval = FALSE }
  ii = which(dt$Region=="cfa4x")
  oo = dt[ii,c("Year", "Licenses", "TAC", "Landings", "Effort", "CPUE")]
  kable(oo, format="simple", row.names=FALSE, align="cccccc",
  caption = "Fishery performance statistics in 4X. Units are: TACs and Landings (tons, t), Effort ($\\times 10^3$ trap hauls, th) and CPUE (kg/th). There were no landings or TACs in 2018/2019 due to indications of low abundance.")
  # \@ref(tab:table-fishery-4x)
  ``` 


-->   

## Human interactions: Fishery landings and TACs {.c}

- The landings in N-ENS for 2022 and 2021 were similar in their spatial patterns.  

- The landings in 4X for 2022 were spatially more contracted than 2021.  

\begin{tiny}
```{r landings-map, echo=FALSE, out.width='45%', fig.show='hold', fig.align='center', fig.cap = 'Snow Crab landings from fisheries logbook data for previous and current years. Units are tons per 10 km x 10 km grid.' }
loc0= file.path( SCD, "output", "maps", "logbook", "snowcrab", "annual", "landings" )
fn1 = file.path( loc0, paste( "landings", year_previous,   "png", sep=".") ) 
fn2 = file.path( loc0, paste( "landings", year.assessment, "png", sep=".") ) 
knitr::include_graphics( c(fn1, fn2 ) )
#  \@ref(fig:landings-map)  
```
\end{tiny}


## Human interactions: Fishery landings and TACs ... {.c}
 
In 2022, landings in all areas were below respective TACs.
 
\begin{tiny}
```{r landings-timeseries, echo=FALSE, out.width='60%',  fig.align='center', fig.cap = 'Landings (t) of Snow Crab on the SSE. For 4X, the year refers to the starting year of the season.  Inset is a closeup view of the timeseries for N-ENS and 4X.'}
include_graphics( file.path( SCD, "assessments", year.assessment, "timeseries", "fishery",   "landings.ts.png" ) )
# \@ref(fig:landings-timeseries)
``` 
\end{tiny}



## Human interactions: Fishery catch rates
 
- Generally, the spatial extent of exploitation was smaller, many of the exploited area show elevated catch rates, 

- In 4X catch rates were lower in 2022.  

```{r cpue-map, echo=FALSE, out.width='45%', fig.show='hold', fig.align='center', fig.cap = 'Snow Crab crude catch rates on the Scotian Shelf for previous and current years. Units are kg/trap haul per 10 km x 10 km grid.' }
loc0= file.path( SCD, "output", "maps", "logbook", "snowcrab", "annual", "cpue" )
fn1 = file.path( loc0, paste( "cpue", year_previous,   "png", sep=".") ) 
fn2 = file.path( loc0, paste( "cpue", year.assessment, "png", sep=".") ) 
knitr::include_graphics( c(fn1, fn2 ) )
# \@ref(fig:cpue-map)  
```


## Human interactions: Fishery catch rates ...
 
\begin{tiny}
```{r cpue-timeseries, echo=FALSE, out.width='60%', fig.align='center', fig.cap = 'Temporal variations in crude catch rates of Snow Crab (kg per trap haul).'}
include_graphics( file.path( SCD, "assessments", year.assessment, "timeseries", "fishery",   "cpue.ts.png" ) ) 
# \@ref(fig:cpue-timeseries)  
```
\end{tiny}


## Human interactions: At-Sea-Observed information

- Target: 5% of landings

- In 2021, both N-ENS and 4X were not sampled by At-Sea-Observers. 

- In 2022, ~ 0.8 % of landings in 4X were sampled by At-Sea-Observers 
 
\begin{tiny}
```{r observer-locations-map, out.width='22%', fig.show='hold', fig.align='center', fig.cap= 'Snow Crab At-sea-observer locations.' }
loc = file.path( SCD, "output", "maps", "observer.locations" )
yrsplot = year.assessment + c(0:-4)
fn4 = file.path( loc, paste( "observer.locations", yrsplot[4], "png", sep=".") )
fn3 = file.path( loc, paste( "observer.locations", yrsplot[3], "png", sep=".") )
fn2 = file.path( loc, paste( "observer.locations", yrsplot[2], "png", sep=".") )
fn1 = file.path( loc, paste( "observer.locations", yrsplot[1], "png", sep=".") )
include_graphics( c( fn4, fn3, fn2, fn1) )
# \@ref(fig:observer-locations-map)   
```
\end{tiny}

\vspace{5mm}

Bycatch: last assessment was in 2017 and levels were << 1% by weight.

<!-- 

\begin{tiny}
```{r observer-CC, echo=FALSE, eval = FALSE, out.width='27%', fig.show='hold', fig.align='center', fig.cap = 'Size frequency distribution of Snow Crab sampled by At-sea-observers, broken down by Carapace Condition (CC). For 4X, the year refers to the starting year of the season. Vertical lines indicate 95 mm Carapace Width, the minimum legal commercial size.' }
  loc = file.path( SCD, "assessments", year.assessment, "figures", "size.freq", "observer")
  fn1 = file.path( loc, paste( "size.freqcfanorth", (year_previous), ".png", sep="" ) )
  fn2 = file.path( loc, paste( "size.freqcfanorth", (year.assessment  ), ".png", sep="" ) )
  fn3 = file.path( loc, paste( "size.freqcfasouth", (year_previous), ".png", sep="" ) )
  fn4 = file.path( loc, paste( "size.freqcfasouth", (year.assessment  ), ".png", sep="" ) )
  fn5 = file.path( loc, paste( "size.freqcfa4x", (year_previous), ".png", sep="" ) )
  fn6 = file.path( loc, paste( "size.freqcfa4x", (year.assessment  ), ".png", sep="" ) )
  include_graphics(  c(fn1, fn2, fn3, fn4, fn5, fn6) )
# \@ref(fig:observer-CC)  
```
\end{tiny}
 
-->


# Stock status  {.c}

## Stock status: survey  {.c}
 
- No survey in 2020 (Covid-19) and incomplete 2022 (mechanical issues). 

- Inshore areas of S-ENS were most affected. 

- N-ENS and CFA 4X were not affected.
 
\begin{tiny}
```{r survey-locations-map, out.width='26%', fig.show='hold', fig.align='center', fig.cap= 'Snow Crab survey locations.' }
loc = file.path( SCD, "output", "maps", "survey.locations" )
yrsplot = setdiff( year.assessment + c(0:-9), 2020)
fn6 = file.path( loc, paste( "survey.locations", yrsplot[6], "png", sep=".") )
fn5 = file.path( loc, paste( "survey.locations", yrsplot[5], "png", sep=".") )
fn4 = file.path( loc, paste( "survey.locations", yrsplot[4], "png", sep=".") )
fn3 = file.path( loc, paste( "survey.locations", yrsplot[3], "png", sep=".") )
fn2 = file.path( loc, paste( "survey.locations", yrsplot[2], "png", sep=".") )
fn1 = file.path( loc, paste( "survey.locations", yrsplot[1], "png", sep=".") )
include_graphics( c(fn3, fn2, fn1) )
# \@ref(fig:survey-locations-map)  
```
\end{tiny}


<!-- 
## Stock status: Size structure 

Factors: early maturation, size-selective predation, fishing of largest individuals, start or end of a recruitment pulse, timing of survey. 
\begin{scriptsize}
```{r meansize-male-mat, out.width='50%', echo=FALSE, eval = FALSE, fig.align='center', fig.cap = 'Mean size of mature male Snow Crab (CW; mm) from surveys with 95\\% Confidence Intervals'}
include_graphics(  file.path( SCD, "assessments", year.assessment, "timeseries", "survey", "cw.mat.png" )  )
# \@ref(fig:meansize-male-mat)
###
\end{scriptsize}
 
-->


## Stock status: Carapace condition of mature male crab 

\begin{small}
\begin{columns}
\begin{column}{.48\textwidth}

\vspace{4mm}

- CC5 crab of mature crab.

\vspace{12mm}

\begin{scriptsize}
\textbf{Figure}: Size-frequency of mature male Snow Crab by carapace width (mm) and carapace condition from surveys. Columns are years and rows are N-ENS (top), S-ENS(middle) and 4X(bottom).
\end{scriptsize}

\end{column}
\begin{column}{.48\textwidth}
\begin{tiny}
```{r sizefeq-male-survey-cc, out.width='45%', fig.show='hold', echo=FALSE, fig.align='center', fig.cap = ''}
  odir = file.path( SCD, "assessments", year.assessment, "figures", "size.freq", "carapacecondition" )
  fn1 = file.path( odir, "sizefreq.cfanorth.2019.png" ) 
  fn2 = file.path( odir, "sizefreq.cfasouth.2019.png" ) 
  fn3 = file.path( odir, "sizefreq.cfa4x.2019.png" ) 
  fn4 = file.path( odir, "sizefreq.cfanorth.2021.png" ) 
  fn5 = file.path( odir, "sizefreq.cfasouth.2021.png" ) 
  fn6 = file.path( odir, "sizefreq.cfa4x.2021.png" ) 
  fn7 = file.path( odir, "sizefreq.cfanorth.2022.png" ) 
  fn8 = file.path( odir, "sizefreq.cfasouth.2022.png" ) 
  fn9 = file.path( odir, "sizefreq.cfa4x.2022.png" ) 
  include_graphics(c( fn4, fn7, fn5, fn8, fn6, fn9) )
# \@ref(fig:sizefeq-male-survey-cc)
```
\end{tiny}
\end{column}
\end{columns}
\end{small}


## Stock status: Recruitment
 
\begin{small}
\begin{columns}
\begin{column}{.48\textwidth}

\vspace{4mm}
\begin{itemize}
  \item Little to no recruitment is expected for the next 1-3 years in N-ENS.
  \item Moderate levels of recruitment are expected in S-ENS. 
  \item Low to moderate levels of recruitment are expected for 2 years in 4X.
\end{itemize}

\vspace{4mm}
\begin{scriptsize}
\textbf{Figure}: Size-frequency (areal density; no/km$^2$) histograms by carapace width of male Snow Crab. The vertical line represents the legal size (95 mm). Immature animals are shown with light coloured bars, mature with dark.
\end{scriptsize}

\end{column}
\begin{column}{.48\textwidth}
```{r sizefeq-male, out.width='90%', echo=FALSE, fig.align='center', fig.cap = ''}
include_graphics(  file.path( SCD, "assessments", year.assessment, "figures", "size.freq", "survey", "male.denl.png" )  )
# \@ref(fig:sizefeq-male)
```
\end{column}
\end{columns}
\end{small}
 



## Stock status: Reproduction

\begin{small}
\begin{columns}
\begin{column}{.48\textwidth}


\begin{itemize}
  \item All areas had recruitment of female crab into the mature (egg-bearing) segment of the population from 2016-2022.
  \item In N-ENS for 2022, a decline in numerical densities, and low densities of adolescent females. 
  \item Egg and larval production is expected to be moderate to high in the next year in all areas except N-ENS. 
\end{itemize}

\vspace{2mm}
\begin{scriptsize}
\textbf{Figure}: Size-frequency (areal density; no/km$^2$) histograms by carapace width of female Snow Crab. Immature animals are shown with light coloured bars, mature with dark.
\end{scriptsize}

\end{column}
\begin{column}{.48\textwidth}
```{r sizefeq-female, out.width='90%', echo=FALSE, fig.align='center', fig.cap = ''}
include_graphics(  file.path( SCD, "assessments", year.assessment, "figures", "size.freq", "survey",  "female.denl.png" )  )
# \@ref(fig:sizefeq-female)
```
\end{column}
\end{columns}
\end{small}





## Stock status: Reproduction ...

```{r fmat-timeseries, out.width='50%', echo=FALSE, fig.align='center', fig.cap = 'Mature female density log$_{10}$(no/km$^2$) from the Snow Crab survey.'  }
include_graphics( file.path( SCD, "assessments", year.assessment, "timeseries", "survey", "totno.female.mat.png") )
# \@ref(fig:fmat-timeseries)
```
 

## Stock status: Mature female 

Distributions are heterogeneous and often in shallower areas.  

```{r fmat-map, echo=FALSE, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Mature female density log$_{10}$(no/km$^2$)  from the Snow Crab survey.' }
loc = file.path( SCD, "output", "maps", "survey", "snowcrab", "annual", "totno.female.mat" )
yrsplot = setdiff( year.assessment + c(0:-9), 2020)
fn4 = file.path( loc, paste( "totno.female.mat", yrsplot[4], "png", sep=".") )
fn3 = file.path( loc, paste( "totno.female.mat", yrsplot[3], "png", sep=".") )
fn2 = file.path( loc, paste( "totno.female.mat", yrsplot[2], "png", sep=".") )
fn1 = file.path( loc, paste( "totno.female.mat", yrsplot[1], "png", sep=".") )
include_graphics( c( fn3, fn2, fn1) )
# \@ref(fig:fmat-map)  
```


## Stock status: Viable Habitat {.t}

\begin{small}
\begin{columns}
\begin{column}{.48\textwidth}
```{r habitat, echo=FALSE, out.width='52%', fig.align='center', fig.show='hold',  fig.cap = 'Persistent habitat.' }
fn1=file.path( media_loc, "viable_habitat.png" ) 
knitr::include_graphics( c(fn1  ) ) 
# \@ref(fig:habitat)  
```
\end{column}
\begin{column}{.48\textwidth}
```{r habitat2, echo=FALSE, out.width='40%', fig.align='center', fig.show='hold',  fig.cap = 'Habitat preferences.' }
fn2=file.path( media_loc, "viable_habitat_depth_temp.png" )
knitr::include_graphics( c( fn2 ) ) 
# \@ref(fig:habitat2)  
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
  
```{r fb-habitat-map, out.width='30%', fig.show='hold', fig.align='center', fig.cap= 'Habitat viability (probability; fishable Snow Crab).' }

loc = file.path( SCD, 'modelled', 'default_fb', 'predicted_habitat' )
vn = "habitat."
yrsplot =  year.assessment + c(0:-10)
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
# \@ref(fig:fb-habitat-map)  
# *Figure XXX. Habitat viability (probability; fishable Snow Crab)* 
```


## Stock status: Viable Habitat ... {.c}

```{r fb-habitat-timeseries, out.width='50%', echo=FALSE,   fig.align='center', fig.cap = 'Habitat viability (probability; fishable Snow Crab). Means and 95\\% Credible Intervals are presented.' }
loc = file.path( SCD, 'modelled', 'default_fb', 'aggregated_habitat_timeseries' )
include_graphics( file.path( loc, 'habitat_M0.png') )
# \@ref(fig:fb-habitat-timeseries)
```
 

  
 
## Stock status: Sex ratios (proportion female, mature)

\begin{small}
\begin{columns}
\begin{column}{.48\textwidth}

\vspace{2mm}
\begin{itemize}
\item Mostly male-dominated: larger size may be protective against predation?
\item Imbalance indicates differential mortality: predation, competition and fishing
\item In 4X, sex ratios are balanced.  
\end{itemize}
\end{column}
\begin{column}{.48\textwidth}
```{r sexratio-mature, out.width='60%', echo=FALSE, fig.align='center', fig.cap = 'Timeseries of sex ratios.'  }
include_graphics( file.path( SCD, "assessments", year.assessment, "timeseries", "survey", "sexratio.mat.png") )
# \@ref(fig:sexratio-mature)
```
\end{column}
\end{columns}
\end{small}
 

```{r sexratio-map, echo=FALSE, out.width='25%', fig.show='hold', fig.align='center', fig.cap= 'Map of sex ratios.'}
yrsplot = setdiff( year.assessment + c(0:-4), 2020)
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


\begin{tiny}
```{r fbgeomean-map, echo=FALSE, out.width='30%', fig.show='hold', fig.align='center', fig.cap= 'Snow Crab survey fishable component biomass density log~10(t/km$^2$). Note, there is no data in 2020.' }
loc = file.path( SCD, 'output', 'maps', 'survey', 'snowcrab', 'annual', 'R0.mass')
yrsplot =  setdiff(year.assessment + c(0:-9), 2020 ) 
fn6 = file.path( loc, paste( 'R0.mass', yrsplot[6], 'png', sep='.') )
fn5 = file.path( loc, paste( 'R0.mass', yrsplot[5], 'png', sep='.') )
fn4 = file.path( loc, paste( 'R0.mass', yrsplot[4], 'png', sep='.') )
fn3 = file.path( loc, paste( 'R0.mass', yrsplot[3], 'png', sep='.') )
fn2 = file.path( loc, paste( 'R0.mass', yrsplot[2], 'png', sep='.') )
fn1 = file.path( loc, paste( 'R0.mass', yrsplot[1], 'png', sep='.') )
include_graphics( c(  fn3, fn2, fn1) )
``` 
\end{tiny}
 


## Stock status: Biomass Density ... {.c}

- A peak in 2009 to 2014 and has since been declining in all areas. 
 

\begin{tiny}

```{r fbGMTS, out.width='50%', echo=FALSE, fig.align='center', fig.cap = 'The crude, unadjusted geometric mean fishable biomass density log~10(t/km$^2$) from the Snow Crab survey. Error bars represent 95\\% Confidence Intervals. Note the absence of data in 2020. Prior to 2004, surveys were conducted in the Spring.'}
fn = file.path(SCD,'assessments', year.assessment, 'timeseries','survey','R0.mass.png')
include_graphics( c(fn) )
#\@ref(fig:fbGMTS)
```
\end{tiny}



## Stock status: Biomass Index (aggregate)
  
A contraction of spatial range in 4X and the western parts of S-ENS were also evident in 2021 to 2022. 
 
```{r fbindex-map, echo=FALSE, out.width='30%', fig.show='hold', fig.align='center', fig.cap= 'Biomass index log~10(t/km$^2$) predicted from the Snow Crab survey.' }
loc = file.path( SCD, 'modelled', 'default_fb', 'predicted_biomass_densities' )
yrsplot =  year.assessment + c(0:-10)
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
    
\begin{tiny}
```{r fbindex-timeseries, out.width='60%', echo=FALSE, fig.align='center', fig.cap = 'The fishable biomass index (t) predicted by CARSTM of Snow Crab survey densities. Error bars represent Bayesian 95\\% Credible Intervals. Note large errors in 2020 when there was no survey.' }
include_graphics( file.path( SCD, 'modelled', 'default_fb', 'aggregated_biomass_timeseries' , 'biomass_M0.png') )
# \@ref(fig:fbindex-timeseries)
```
\end{tiny}


<!-- 

  ## Stock status: Modelled Biomass (pre-fishery) {.c}
  
  N-ENS: #r round(B_north[t0], 2)` t in #r year.assessment`

    - #r round(B_north[t1], 2)` t in #r year_previous`. 

  S-ENS: #r round(B_south[t0], 2)` t in #r year.assessment`

    - #r round(B_south[t1], 2)` t in #r year_previous`. 

  4X:  #r round(B_4x[t0], 2)` t in #r year.assessment`-#r year.assessment+1`

    - #r round(B_4x[t1], 2)` t for the #r year_previous`-#r year.assessment` season. 
 
-->


## Stock status: Modelled Biomass (pre-fishery) ... {.c}

\begin{tiny}
```{r logisticPredictions, out.width='32%', echo=FALSE, fig.show='hold', fig.align='center', fig.cap = 'Model 1 fishable, posterior mean modelled biomass (pre-fishery; kt) are shown in dark orange for N-ENS, S-ENS and 4X (left, middle and right). Light orange are posterior samples of modelled biomass (pre-fishery; kt) to illustrate the variability of the predictions. The biomass index (post-fishery, except prior to 2004) after model adjustment by the model catchability coefficient is in gray.' } 
loc = file.path( SCD, 'fishery_model', year.assessment, 'logistic_discrete_historical' )
fn1 = file.path( loc, 'plot_predictions_cfanorth.png' ) 
fn2 = file.path( loc, 'plot_predictions_cfasouth.png' ) 
fn3 = file.path( loc, 'plot_predictions_cfa4x.png' ) 
include_graphics(c(fn1, fn2, fn3) )
# \@ref(fig:logisticPredictions)
```
\end{tiny}



<!-- 
## Stock status: Fishing Mortality  {.c}

N-ENS: #r round(FM_north[t0],3)` (annual exploitation rate of #r round(100*(exp(FM_north[t0])-1),2)`%) in #r year.assessment`

  - Up from the  #r year_previous` rate of #r round(FM_north[t1],3)` (annual exploitation rate of #r round(100*(exp(FM_north[t1])-1),1)`%)
 
S-ENS: #r round(FM_south[t0],3)` (annual exploitation rate of #r round(100*(exp(FM_south[t0])-1),1)`%) in #r year.assessment`

  - Decreasing marginally from the #r year_previous` rate of #r round(FM_south[t1],3)` (annual exploitation rate of #r round(100*(exp(FM_south[t1])-1),1)`%)

4X: #r round(FM_4x[t0],3)` (annual exploitation rate of #r round(100*(exp(FM_4x[t0])-1),1)`%) in #r year.assessment`-#r year.assessment+1` season 

  - Decreasing from the #r year.assessment-1`-#r year.assessment` season rate of #r round(FM_4x[t1],3)` (annual exploitation rate of #r round(100*(exp(FM_4x[t1])-1),1)`%)

Localized exploitation rates are likely higher, as not all areas for which biomass is estimated are fished. 
 

-->


## Stock status: Fishing Mortality ... {.c}

```{r logisticFishingMortality, out.width='32%', echo=FALSE,  fig.show='hold', fig.align='center', fig.cap = 'Time-series of modelled instantaneous fishing mortality from Model 1, for N-ENS (left), S-ENS (middle), and 4X (right). Samples of the posterior densities are presented, with the darkest line being the mean.' }
  odir = file.path( fishery_model_results, year.assessment, "logistic_discrete_historical" )
  fn1 = file.path( odir, "plot_fishing_mortality_cfanorth.png" ) 
  fn2 = file.path( odir, "plot_fishing_mortality_cfasouth.png" ) 
  fn3 = file.path( odir, "plot_fishing_mortality_cfa4x.png" ) 
include_graphics(c(fn1, fn2, fn3) )
# \@ref(fig:logisticFishingMortality)
```



## Stock status: Reference Points {.c}
 

```{r ReferencePoints, out.width='40%', echo=FALSE, fig.align='center', fig.cap = 'Harvest control rules for the Scotian Shelf Snow Crab fisheries.' }
include_graphics( file.path( params$media_loc, 'harvest_control_rules.png') ) 
# \@ref(fig:ReferencePoints)
```


## Stock status: Reference Points ... {.c}

```{r logistic-hcr, out.width='32%', echo=FALSE, fig.show='hold', fig.align='center', fig.cap = 'Reference Points (fishing mortality and modelled biomass) from Model 1, for N-ENS (left), S-ENS (middle), and 4X (right). The large yellow dot indicates most recent year and the 95\\% CI.' }
  odir = file.path( fishery_model_results, year.assessment, "logistic_discrete_historical" )
  fn1 = file.path( odir, 'plot_hcr_cfanorth.png' ) 
  fn2 = file.path( odir, 'plot_hcr_cfasouth.png' ) 
  fn3 = file.path( odir, 'plot_hcr_cfa4x.png' ) 
  include_graphics(c(fn1, fn2, fn3) )
#  \@ref(fig:logistic-hcr)
```

  - N-ENS is in the "healthy" zone
  - S-ENS is in the "healthy" zone
  - 4X is in the "cautious" zone


<!-- 
 
\begin{footnotesize}
*Table 6. Reference points from the logistic biomass dynamics fishery model (Model 1): K is Carrying capacity (kt); and r is Intrinsic rate of increase (non-dimensional). Note that FMSY (fishing mortality associated with 'Maximum Sustainable Yield') is r/2. Similarly, BMSY (biomass associated with 'Maximum Sustainable Yield') is K/2. SD is posterior Standard deviations.*
\end{footnotesize}


|       |  $K$ [SD] | $r$ [SD] |
| :---: |    :----:   | :---: |
| N-ENS | #r K_north` [#r K_north_sd`] | #r r_north` [#r r_north_sd`] |
| S-ENS | #r K_south` [#r K_south_sd`] | #r r_south` [#r r_south_sd`] |
| 4X    | #r K_4x`   [#r K_4x_sd`] | #r r_4x` [#r r_4x_sd`] |


--> 






<!--

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
 

-->



## References

\begin{tiny}


Banerjee, S., Carlin, B. P., and Gelfand, A. E.. 2004. Hierarchical Modeling and Analysis for Spatial Data. Monographs on Statistics and Applied Probability. Chapman and Hall/CRC.

Besag, Julian. 1974. Spatial interaction and the statistical analysis of lattice systems. Journal of the Royal Statistical Society Series B (Methodological) 1974: 192-236.

Canada Gazette. 2022. [Regulations Amending the Fishery (General) Regulations. Part II, Volume 156, Number 8.](https://www.gazette.gc.ca/rp-pr/p2/2022/2022-04-13/html/sor-dors73-eng.html)

Canada Gazette. 2016. St. Anns Bank Marine Protected Area Regulations. Canada Gazette, Part I, Vol 150, Issue 51: 4143-4149.

Choi, J.S. 2020. A Framework for the assessment of Snow Crab (*Chioneocete opilio*) in Maritimes Region (NAFO Div 4VWX) . DFO Can. Sci. Advis. Sec. Res. Doc. 2020/nnn. v + xxx p.

Choi, J.S. 2022. Reconstructing the Decline of Atlantic Cod with the Help of Environmental Variability in the Scotian Shelf of Canada. bioRxiv. https://doi.org/10.1101/2022.05.05.490753.

Choi, J. S., and B. C. Patten. 2001. Sustainable Development: Lessons from the Paradox of Enrichment. Ecosystem Health 7: 163–77.

Choi, Jae S., B. Cameron, K. Christie, A. Glass, and E. MacEachern. 2022. Temperature and Depth Dependence of the Spatial Distribution of Snow Crab. bioRxiv. https://doi.org/10.1101/2022.12.20.520893.


DFO. 2013. [Integrated Fisheries Management Plan for Eastern Nova Scotia and 4X Snow Crab (*Chionoecetes Opillio*.)](http://www.dfo-mpo.gc.ca/fm-gp/peches-fisheries/ifmp-gmp/snow-crab-neige/snow-crab-neiges2013-eng.htm)

\end{tiny}


## References ...

\begin{tiny}

DFO. 2018. Stock Status Update of Atlantic Halibut (Hippoglossus hippoglossus) on the Scotian Shelf and Southern Grand Banks in NAFO Divisions 3NOPs4VWX5Zc. DFO Can. Sci. Advis. Sec. Sci. Resp. 2018/022.

Hebert M, Miron G, Moriyasu M, Vienneau R, and DeGrace P. Efficiency and ghost fishing of Snow Crab (Chionoecetes opilio) traps in the Gulf of St. Lawrence. Fish Res. 2001; 52(3): 143-153. 10.1016/S0165-7836(00)00259-9   
 
Riebler, A., Sørbye, S.H., Simpson D., and Rue, H. 2016. An intuitive Bayesian spatial model for disease mapping that accounts for scaling. Statistical methods in medical research 25: 1145-1165.

Simpson, D., Rue, H., Riebler, A., Martins, T.G., and Sørbye, SH. 2017. Penalising Model Component Complexity: A Principled, Practical Approach to Constructing Priors. Statist. Sci. 32: 1-28.


\end{tiny}
 
# END  
 
