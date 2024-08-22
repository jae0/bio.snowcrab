---
title: "quick_test to try a page before adding to main document "
subtitle: "this makes debugging a page easier and faster"
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
  - \usepackage{float}
  - \usepackage{multicol}
  # - \usepackage{subfig}
  # - \newcommand{\btiny}{\begin{tiny}}
  # - \newcommand{\etiny}{\end{tiny}}
params:
  year.assessment: 2023
  media_loc: "media"
  debugging: FALSE
  loc_dde: ""
--- 


<!-- Preamble

This is a Markdown document ... To create HTML or PDF, etc, run: 


# for presentations to PDF (via beamer):
# note: section separation with '#' can confuse rmarkdown
  
  make rmarkdown FN=quick_test YR=2023 SOURCE=~/projects/bio.snowcrab/inst/markdown WK=~/bio.data/bio.snowcrab/assessments  DOCTYPE=beamer_presentation DOCEXTENSION=pdf # {via Rmarkdown}
 
  make rmarkdown FN=quick_test YR=2023 SOURCE=~/projects/bio.snowcrab/inst/markdown WK=~/bio.data/bio.snowcrab/assessments  DOCTYPE=html_document DOCEXTENSION=html # {via Rmarkdown}

# for html documents including presentations:
  make quarto FN=quick_test YR=2023 SOURCE=~/projects/bio.snowcrab/inst/markdown WK=~/bio.data/bio.snowcrab/assessments  DOCEXTENSION=pdf # {via Quarto}


  make pdf FN=quick_test  # {via pandoc}

Alter year and directories to reflect setup or copy Makefile and alter defaults to your needs.
   
Huge > huge > LARGE > Large > large > normalsize > small > footnotesize > scriptsize > tiny

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

  sn_env = snowcrab_load_key_results_to_memory( year.assessment, debugging=params$debugging, return_as_list=TRUE  ) 

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
  



# Ecosystem considerations


## Connectivity: Oceanic currents {.c}

```{r ocean_currents, echo=FALSE, out.width='55%', fig.align='center', fig.show='hold',  fig.cap = 'Ocean currents in the Martimes. Source: https://www.dfo-mpo.gc.ca/oceans/publications/soto-rceo/2018/atlantic-ecosystems-ecosystemes-atlantiques/index-eng.html' }
fn2=file.path( media_loc, "maritimes_currents.png" )
knitr::include_graphics( c(fn2) ) 
# \@ref(fig:movementtracks)  
``` 


## Connectivity: Movement {.c}

::: columns 
:::: column 
 
```{r movementtracks, echo=FALSE, out.width='71%', fig.align='center', fig.show='hold',  fig.cap = 'Snow Crab movement tracks from mark and recapture.' }
fn1=file.path( media_loc, "movement0.png" )
fn2=file.path( media_loc, "movement.png" )
knitr::include_graphics( c(fn1, fn2) ) 
# \@ref(fig:movementtracks)  
``` 
 
 
::::
:::: column

```{r movement, echo=FALSE, out.width='55%', fig.align='center', fig.show='hold',  fig.cap = 'Snow Crab movement distances and rates.' }
fn1=file.path( media_loc, "snowcrab_movement_distances.png" )
fn2=file.path( media_loc, "snowcrab_movement_rates.png" )
knitr::include_graphics( c(fn1, fn2) ) 
# \@ref(fig:movement)  
``` 
 
::::
:::


## Bathymetry {.c}
::: columns 

:::: column 

```{r bathymetry-map, out.width='100%', echo=FALSE, fig.align='center', fig.cap='Variations in log(depth; m) in the Scotian Shelf region.' }
bathydir = file.path( data_root, 'aegis', 'bathymetry', 'modelled', 'default', 'stmv', 'none_fft', 'z', 'maps', 'SSE' )
knitr::include_graphics( file.path( bathydir, 'bathymetry.z.SSE.png' ) )
```
::::
:::: column
```{r bathymetry-map2, out.width='100%', echo=FALSE, fig.align='center', fig.cap='Local SD log(depth; m) in the Scotian Shelf region.' }
bathydir = file.path( data_root, 'aegis', 'bathymetry', 'modelled', 'default', 'stmv', 'none_fft', 'z', 'maps', 'SSE' )
knitr::include_graphics( file.path( bathydir, 'bathymetry.b.sdSpatial.SSE.png' ) )
```

::::
:::


## Substrate {.c}
::: columns 

:::: column 

```{r substrate-map, out.width='100%', echo=FALSE, fig.align='center', fig.cap='Substrate grain size log(mm) variations in the Scotian Shelf region.' }
substrdir = file.path( data_root, 'aegis', 'substrate', 'maps', 'canada.east.highres' )
knitr::include_graphics( file.path(  substrdir, 'substrate.substrate.grainsize.canada.east.highres.png' ) )
```
::::
:::: column
```{r substrate-map2, out.width='100%', echo=FALSE, fig.align='center', fig.cap='Local SD substrate grain size log(mm) in the Scotian Shelf region.' }
substrdir = file.path( data_root, 'aegis', 'substrate', 'maps', 'canada.east.highres' )
knitr::include_graphics( file.path( substrdir, 'substrate.s.sdSpatial.canada.east.highres.png' ) )
```

::::
:::


## Bottom Temperature {.c}
::: columns 

:::: column 
\vspace{12mm}
```{r bottom-temperatures-survey, out.width='90%', echo=FALSE, fig.align='center', fig.cap = 'Annual variations in bottom temperature observed during the Snow Crab survey. The horizontal (black) line indicates the long-term, median temperature within each subarea. Error bars represent standard errors.' }
knitr::include_graphics( file.path( SCD, 'assessments', year.assessment, 'timeseries', 'survey', 't.pdf') )
# \@ref(fig:bottom-temperatures-survey)
```
::::
 
:::: column 
```{r bottom-temperatures, out.width='75%', echo=FALSE, fig.align='center', fig.cap = 'Posterior densities of predicted average bottom temperatures. Red horizontal line is at $7^\\circ$C.' }
knitr::include_graphics( file.path( SCD, 'assessments', year.assessment, 'timeseries', 'temperature_bottom.pdf') )
# \@ref(fig:bottom-temperatures)
```
::::
:::


## Bottom Temperature ... {.c}

```{r bottom-temperatures-map, out.width='30%', echo=FALSE, fig.show='hold', fig.align='center', fig.cap = 'Spatial variations in predicted (1 September) bottom temperature from 2021 (left) to 2023 (right). Groundfish surveys were not conducted in 2020 and 2022 in the snow crab domain. Snow crab surveys were not conducted in 2020, and incomplete in 2022 in S-ENS.' }

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

 
 
## Bottom Temperature ... {.c}
::: columns 
:::: column 

```{r bottom-temperatures-spatialeffect, out.width='90%', echo=FALSE, fig.align='center', fig.cap = 'Persistent spatial effect of bottom temperature, adjusting for spatiotemporal variability and autocorrelations. Time period from 1999 to present.' }
loc = file.path( data_root, 'aegis', 'temperature', 'modelled', 'default', 'maps' )
knitr::include_graphics( file.path( loc, 'space_re_total.png') )
# \@ref(fig:bottom-temperatures-spatialeffect)
```
 
::::
 
:::: column


\vspace{12mm}

- Persistent spatial gradient of  $>1^\circ$C in bottom temperatures in the Maritimes Region. 

- Variable due to confluence:

  - Warm, high salinity Gulf Stream from the S-SE along the shelf edge 

  - Cold, low salinity Labrador Current

  - Cold low salinity St. Lawrence outflow from the N-NE

  - Nearshore Nova Scotia current, running from the NE. 

::::
:::



## Species composition
```{r speciesomposition0, echo=FALSE, out.width='75%', fig.align='center', fig.show='hold', fig.cap = 'Species ordination (PCA: eigenanalysis of correlation matrices). PC1 is associatd with bottom temperatures. PC2 is associated with depth. Snow crab is shown as an orange dot.' }
xlab = paste("PC1 (", pca$variance_percent[1], "%)", sep="" )
ylab = paste("PC2 (", pca$variance_percent[2], "%)", sep="" )
plot( PC2 ~ PC1, pcadata, type="n", xlab=xlab, ylab=ylab )
text( PC2 ~ PC1, labels=vern, data=pcadata, cex=0.7, col="slateblue"  )
i = grep("Snow crab", pcadata$vern, ignore.case=TRUE)
points( PC2 ~ PC1, pcadata[i,], pch=19, cex=3.0, col="darkorange" )
```

 

## Species composition PC1


```{r speciesomposition1, echo=FALSE, out.width='45%', fig.align='center', fig.show='hold',  fig.cap = 'Species composition in space and time. Primary gradient is related to bottom temperatures. Groundfish surveys were not conducted in 2020 and 2022 in the snow crab domain. Snow crab surveys were not conducted in 2020, and incomplete in 2022 in S-ENS.' }
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

```{r speciesomposition2, echo=FALSE, out.width='45%', fig.align='center', fig.show='hold',  fig.cap = 'Species composition in space and time. Primary gradient is related to bottom temperatures. Groundfish surveys were not conducted in 2020 and 2022 in the snow crab domain. Snow crab surveys were not conducted in 2020, and incomplete in 2022 in S-ENS.' }
spc_loc = file.path( data_root, 'aegis', 'speciescomposition', 'modelled', 'default', 'maps' )
fn1 = file.path( spc_loc, 'pca1.space_re_total.png') 
fn2 = file.path( spc_loc, 'pca2.space_re_total.png')
ts_loc = file.path( data_root, 'aegis', 'speciescomposition', 'modelled', 'default', 'figures' )
fn3 = file.path( ts_loc, 'pca1_time.png') 
fn4 = file.path( ts_loc, 'pca2_time.png') 
knitr::include_graphics( c(fn2,  fn4  ) ) 
# \@ref(fig:habitat3)  
``` 
  
  


## Disease

```{r bcd, out.width='80%', echo=FALSE,   fig.align='center', fig.cap = 'Bitter Crab Disease in Maritimes Region.' }
include_graphics( file.path( SCD, 'output', 'bcd.png') )
```
 


## Entanglements of large megafauna 

```{r map-entanglements, out.width='70%', echo=FALSE, fig.align='center', fig.cap = 'Entanglements of large megafauna in the Maritimes Region. Key: whales (red), leatherback turtles (green), basking shark (blue).' }
region="cfaall"
o = observer.db( DS="bycatch_summary", p=p,  yrs=p$yrs, region=region )   
oss = o$oss  # subset for region of interest
# print("whale entaglements:")
whales = oss[ grep("whale", common, ignore.case=TRUE), ]
# print(whales[, .N, by=.(yr)] )
# print("leatherback entaglements:")
leatherback = oss[ grep("LEATHERBACK", common, ignore.case=TRUE), ]
# print(leatherback[, .N, by=.(yr)])
# print("basking sharks entaglements:")
basking_shark = oss[ grep("BASKING SHARK",  common, ignore.case=TRUE), ]
# print(basking_shark[, .N, by=.(yr)])
plot(lat~-lon, oss, pch=".", col="lightgray", xlim=c(-65.2, -57), ylim=c(42.9,47) )
points(lat~-lon, whales, pch=19, cex=1.5, col="darkred" )
points(lat~-lon, leatherback, pch=18, cex=1.5, col="darkgreen" )
points(lat~-lon, basking_shark, pch=17, cex=1.5, col="slateblue" )
```


## Bycatch Maritimes {.c}
 
```{r bycatch-cpue, out.width='70%', echo=FALSE, fig.align='center', fig.cap = 'Bycatch rates in Maritimes Region. Low levels attributable to trap design (top entry, conical, large mesh 5.25" knot-to-knot) permits escapement of non-target species.' }
o = o_cfaall    
o$bycatch_table[ o$bycatch_table==0 ] = NA
o$bycatch_table[ is.na(o$bycatch_table) ] = "."
o$bycatch_table_catch[ o$bycatch_table_catch==0 ] = NA
o$bycatch_table_catch[ is.na(o$bycatch_table_catch) ] = "."
plot( o$spec ~ o$bct, xlab = "At sea observed catch rate in snow crab fishery (kg/trap)", ylab="Species", type="p", cex=0.9, pch=19, col="darkorange", xlim=c(0, max(o$bct, na.rm=TRUE)*1.4), yaxt="n" )  
text( o$bct, o$spec,  labels=o$species, pos=4, srt=0 , cex=0.5, col="darkslateblue")
text( max(o$bct, na.rm=TRUE)*0.88, 2.5, labels=paste( "Snow crab CPUE (At sea obs., mean): ", o$bct_sc, " kg/trap"), col="darkred", cex=1.0 )
```


## Bycatch Maritimes ... {.c}
```{r bycatch-speciesordination_all, out.width='70%', echo=FALSE, fig.align='center', fig.cap = 'Bycatch as potentially interacting species in Maritimes.' }
o = o_cfaall    
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



    
    
## Bycatch N-ENS {.c}

\tiny
 
```{r bycatch-cpue-n-ens, out.width='12%', echo=FALSE, fig.align='center', fig.cap = 'Bycatch (kg) in N-ENS for those greater than 10 kg/year.' }
o = o_cfanorth    
o$bycatch_table = o$bycatch_table[ which(o$bycatch_table$"Average/Moyen" > 10 ),]
o$bycatch_table$"Average/Moyen" = round(o$bycatch_table$"Average/Moyen")
o$bycatch_table[ o$bycatch_table==0 ] = NA
o$bycatch_table[ is.na(o$bycatch_table) ] = "."
gt(o$bycatch_table)
``` 
\normalsize


## Bycatch N-ENS ... {.c}
```{r bycatch-speciesordination_n-ens, out.width='70%', echo=FALSE, fig.align='center', fig.cap = 'Bycatch as potentially interacting species in N-ENS.' }
o = o_cfanorth    
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

    
## Bycatch S-ENS {.c}

\tiny
 
```{r bycatch-cpue-s-ens, out.width='12%', echo=FALSE, fig.align='center', fig.cap = 'Bycatch (kg) in S-ENS for those greater than 10 kg/year' }
o = o_cfasouth    
o = observer.db( DS="bycatch_summary", p=p,  yrs=p$yrs, region=region )   
o$bycatch_table = o$bycatch_table[ which(o$bycatch_table$"Average/Moyen" > 10 ),]
o$bycatch_table$"Average/Moyen" = round(o$bycatch_table$"Average/Moyen")
o$bycatch_table[ o$bycatch_table==0 ] = NA
o$bycatch_table[ is.na(o$bycatch_table) ] = "."
gt(o$bycatch_table)
``` 
\normalsize



## Bycatch S-ENS ... {.c}
```{r bycatch-speciesordination_s-ens, out.width='70%', echo=FALSE, fig.align='center', fig.cap = 'Bycatch as potentially interacting species in S-ENS.' }
o = o_cfasouth
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

## Bycatch 4X {.c}

\tiny
 
```{r bycatch-cpue-4x, out.width='12%', echo=FALSE, fig.align='center', fig.cap = 'Bycatch (kg) in 4x for those greater than 10 kg/year.' }
o = o_cfa4x
o = observer.db( DS="bycatch_summary", p=p,  yrs=p$yrs, region=region )   
o$bycatch_table = o$bycatch_table[ which(as.numeric(o$bycatch_table$"Average/Moyen") > 10 ),]
o$bycatch_table$"Average/Moyen" = round(o$bycatch_table$"Average/Moyen")
o$bycatch_table[ o$bycatch_table==0 ] = NA
o$bycatch_table[ is.na(o$bycatch_table) ] = "."
gt(o$bycatch_table)
``` 
\normalsize



## Bycatch 4X ... {.c}
```{r bycatch-speciesordination_4x, out.width='70%', echo=FALSE, fig.align='center', fig.cap = 'Bycatch as potentially interacting species in 4X.' }
o = o_cfa4x
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



## Prey

::: columns 
:::: column 
```{r diet, echo=FALSE, out.width='100%', fig.align='center', fig.show='hold',  fig.cap = 'Relative location of snow crab prey (green) in the species composition ordination. Snow crab in orange.' }
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
::::
:::: column
- echinoderms
- polychaete worms (*Maldane*, *Nereis*), worm-like animals
- detritus (dead organic matter)
- large zooplankton, shrimp 
- juvenile crab (Rock Crab; Toad Crab; Lesser Toad Crab)
- Ocean Quahog (*Artica islandica*), bivalve molluscs (*Mytilus* sp, *Modiolus*, *Hiatella*)
- brittle stars (*Ophiura*, *Ophiopholis*)
- sea anemones (*Edwardsia*, *Metridium*). 
::::
:::



## Predators {.c}

::: columns 
:::: column 
```{r predator_ord, echo=FALSE, out.width='100%', fig.align='center', fig.show='hold',  fig.cap = 'Main predators of snow crab on Scotian Shelf of Atlantic Canada (1999-2020). Relative location of snow crab predators (green) in the species composition ordination. Snow crab in orange. Of 58,287 finfish stomach samples, 159 had snow crab (0.28%). There is no information on snow crab diet in the database.' }
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
::::
:::: column

```{r predators2, echo=FALSE} 
kable( counts[1:11,], format="simple", row.names=FALSE)
```
::::
:::


 
##  Predators {.c}
  
::: columns 
:::: column 

```{r predators3a, out.width='100%', echo=FALSE, fig.align='center', fig.cap = 'Location of predators with snow crab in stomach. 2000-2010', warning = FALSE, message = FALSE } 
ggplot() + 
    borders("world", fill = "lightgray", colour = "grey80") + 
    xlim(c( -65.5, -57.1)) + ylim(c(42.1, 47.1)) +
    geom_point(data=snowcrab_predators[year(timestamp) %in% c(2000:2010),], aes(x=slongdd, y=slatdd, colour=Predator ), size=2.5 ) +
    labs(x="", y="", caption="2000-2010") +
    theme(legend.position="inside", legend.position.inside =c(0.16, 0.65), legend.title=element_blank(), legend.text=element_text(size=7.0) 
) 
```
::::
 
:::: column
```{r predators3b, out.width='100%', echo=FALSE, fig.align='center', fig.cap = 'Location of predators with snow crab in stomach. 2011-2020', warning = FALSE, message = FALSE } 
ggplot() + 
    borders("world", fill = "lightgray", colour = "grey80") + 
    xlim(c( -65.5, -57.1)) + ylim(c(42.1, 47.1)) +
    geom_point(data=snowcrab_predators[year(timestamp) %in% c(2011:2020),], aes(x=slongdd, y=slatdd, colour=Predator ), size=2.5 ) +
    labs(x="", y="", caption="2011-2020") +
    theme(legend.position="inside", legend.position.inside =c(0.16, 0.725), legend.title=element_blank(), legend.text=element_text(size=7.0) )
```
::::
:::



## Predators - Atlantic cod ... {.c}

```{r cod-map, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Atlantic cod density log10(no/km$^2$) from the Snow Crab survey. Predation mortality more likely higher encounter rates (high densities co-occur).' }
loc = file.path( SCD, 'output', 'maps', 'survey', 'snowcrab', 'annual', 'bycatch', 'ms.no.10' )
yrsplot = setdiff( year.assessment + c(0:-9), 2020)
fn4 = file.path( loc, paste( 'ms.no.10', yrsplot[4], 'png', sep='.') )
fn3 = file.path( loc, paste( 'ms.no.10', yrsplot[3], 'png', sep='.') )
fn2 = file.path( loc, paste( 'ms.no.10', yrsplot[2], 'png', sep='.') )
fn1 = file.path( loc, paste( 'ms.no.10', yrsplot[1], 'png', sep='.') )
include_graphics( c(  fn3, fn2, fn1) )
# \@ref(fig:cod-map)  
```
 


## Predators - Atlantic cod {.c}

```{r cod-timeseries, out.width='60%', echo=FALSE,  fig.align='center', fig.cap = 'Atlantic cod crude, unadjusted geometric mean numerical density (no/km$^2$) from annual Snow Crab survey. Error bars are 95\\%  Confidence Intervals.'}
include_graphics( file.path( SCD, 'assessments', year.assessment, 'timeseries', 'survey', 'ms.no.10.pdf') )
# \@ref(fig:cod-timeseries)
```


 
##  Predators - Atlantic Halibut  {.c}

```{r halibut-map, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Halibut density log10(no/km$^2$) from the Snow Crab survey. Predation mortality more likely higher encounter rates (high densities co-occur).' }
loc = file.path( SCD, 'output', 'maps', 'survey', 'snowcrab', 'annual', 'bycatch', 'ms.no.30' )
yrsplot = setdiff( year.assessment + c(0:-9), 2020)
fn4 = file.path( loc, paste( 'ms.no.30', yrsplot[4], 'png', sep='.') )
fn3 = file.path( loc, paste( 'ms.no.30', yrsplot[3], 'png', sep='.') )
fn2 = file.path( loc, paste( 'ms.no.30', yrsplot[2], 'png', sep='.') )
fn1 = file.path( loc, paste( 'ms.no.30', yrsplot[1], 'png', sep='.') )
include_graphics( c(  fn3, fn2, fn1) )
# \@ref(fig:halibut-map)  
```
 

##  Predators - Atlantic Halibut ... {.c}

```{r halibut-timeseries, out.width='50%', echo=FALSE,   fig.align='center', fig.cap = 'Atlantic Halibut crude, unadjusted geometric mean numerical density (no/km$^2$) from annual Snow Crab survey. Error bars are 95\\%  Confidence Intervals.' }
include_graphics( file.path( SCD, 'assessments', year.assessment, 'timeseries', 'survey', 'ms.no.30.pdf') )
# \@ref(fig:halibut-timeseries)
```


<!--
~~~{=comment}


## Predators - Thorny skate {.c}

```{r thornyskate-timeseries, out.width='60%', echo=FALSE,  fig.align='center', fig.cap = 'Thorny Skate crude, unadjusted geometric mean numerical density (no/km$^2$) from annual Snow Crab survey. Error bars are 95\\%  Confidence Intervals.'}
include_graphics( file.path( SCD, 'assessments', year.assessment, 'timeseries', 'survey', 'ms.no.201.pdf') )
# \@ref(fig:thornyskate-timeseries)
```

## Predators - Thorny skate ... {.c}

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


##  Predators - Striped Atlantic Wolffish  {.c}
 
```{r Wolffish-timeseries, out.width='60%', echo=FALSE,   fig.align='center', fig.cap = 'Striped Atlantic Wolffish crude, unadjusted geometric mean numerical density (no/km$^2$) from annual Snow Crab survey. Error bars are 95\\%  Confidence Intervals.' }
include_graphics( file.path( SCD, 'assessments', year.assessment, 'timeseries', 'survey', 'ms.no.50.pdf') )
# \@ref(fig:Wolffish-timeseries)
```

##  Predators - Striped Atlantic Wolffish  ... {.c}


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


~~~
-->  


## Competitors/Prey - Lesser toad crab ... {.c}


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

## Competitor/Prey - Lesser toad crab    {.c}


```{r lessertoadcrab-timeseries, out.width='60%', echo=FALSE,   fig.align='center', fig.cap = 'Lesser Toad Crab crude, unadjusted geometric mean numerical density (no/km$^2$) from annual Snow Crab survey. Error bars are 95\\%  Confidence Intervals.' }
include_graphics( file.path( SCD, 'assessments', year.assessment, 'timeseries', 'survey', 'ms.no.2521.pdf') )
# \@ref(fig:lessertoadcrab-timeseries)
```
 

## Competitor/Prey - Northern shrimp {.c}

```{r Shrimp-map, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Northern Shrimp density log10(no/km$^2$) from the Snow Crab survey. Shrimp with similar habitat preferences have declined, possibly due to large-scaled habitat variations and predation. Sampling was incomplete in 2020 and 2022 in S-ENS.' }
loc = file.path( SCD, 'output', 'maps', 'survey', 'snowcrab', 'annual', 'bycatch', 'ms.no.2211' )
yrsplot = setdiff( year.assessment + c(0:-9), 2020)
fn4 = file.path( loc, paste( 'ms.no.2211', yrsplot[4], 'png', sep='.') )
fn3 = file.path( loc, paste( 'ms.no.2211', yrsplot[3], 'png', sep='.') )
fn2 = file.path( loc, paste( 'ms.no.2211', yrsplot[2], 'png', sep='.') )
fn1 = file.path( loc, paste( 'ms.no.2211', yrsplot[1], 'png', sep='.') )
include_graphics( c( fn3, fn2, fn1) )
# \@ref(fig:Shrimp-map)  
```


## Competitor/Prey - Northern shrimp ...  {.c}


```{r Shrimp-timeseries, out.width='60%', echo=FALSE,   fig.align='center', fig.cap = 'Northern Shrimp crude, unadjusted geometric mean numerical density (n/$km^2$) from annual Snow Crab survey. Error bars are 95\\%  Confidence Intervals.' }
include_graphics( file.path( SCD, 'assessments', year.assessment, 'timeseries', 'survey', 'ms.no.2211.pdf') )
# \@ref(fig:Shrimp-timeseries)
```
  


 
 
## Biomass Density  {.c}
  
 
\begin{tiny}
```{r fbgeomean-map, echo=FALSE, out.width='30%', fig.show='hold', fig.align='center', fig.cap= 'Snow Crab survey fishable component biomass density log~10(t/km$^2$). Note, there is no data in 2020. Note that high and low biomass density areas fluctuate with time.' }
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
 
    

## Biomass Density ... {.c}
 

\begin{tiny}

```{r fbGMTS, out.width='65%', echo=FALSE, fig.align='center', fig.cap = 'The crude, unadjusted geometric mean fishable biomass density log~10(t/km$^2$) from the Snow Crab survey. Error bars represent 95\\% Confidence Intervals. Note the absence of data in 2020. Prior to 2004, surveys were conducted in the Spring. A peak in 2009 to 2014 and has since been declining in all areas. '}
fn = file.path(SCD,'assessments', year.assessment, 'timeseries','survey','R0.mass.pdf')
include_graphics( c(fn) )
#\@ref(fig:fbGMTS)
```
\end{tiny}



## Viable Habitat {.t}

\begin{small}
\begin{columns}
\begin{column}{.48\textwidth}
```{r habitat, echo=FALSE, out.width='52%', fig.align='center', fig.show='hold',  fig.cap = 'Persistent habitat.' }
loc=file.path( SCD, 'modelled', 'default_fb', 'maps' )
fn1=file.path( loc, "Predicted_presence_absence.space_re_total.png" ) 
knitr::include_graphics( c(fn1  ) ) 
# \@ref(fig:habitat)  
```
\end{column}
\begin{column}{.48\textwidth}
```{r habitat2, echo=FALSE, out.width='40%', fig.align='center', fig.show='hold',  fig.cap = 'Habitat preferences.' }
loc=file.path( SCD, 'modelled', 'default_fb', 'figures' )
fn2=file.path( loc, "habitat_time.png" )
knitr::include_graphics( c( fn2 ) ) 
# \@ref(fig:habitat2)  
```
\end{column}
\end{columns}
\end{small}




## Viable Habitat ... {.c}
  
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


## Viable Habitat ... {.c}

```{r fb-habitat-timeseries, out.width='70%', echo=FALSE,   fig.align='center', fig.cap = 'Habitat viability (probability; fishable Snow Crab). Means and 95\\% Credible Intervals are presented.' }
loc = file.path( SCD, 'modelled', 'default_fb', 'aggregated_habitat_timeseries' )
include_graphics( file.path( loc, 'habitat_M0.png') )
# \@ref(fig:fb-habitat-timeseries)
```
 


    
## Biomass Index (aggregate)
 
```{r fbindex-map, echo=FALSE, out.width='30%', fig.show='hold', fig.align='center', fig.cap= 'Biomass index log$_{10}$(t/km$^2$) predicted from the Snow Crab survey.' }
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



 