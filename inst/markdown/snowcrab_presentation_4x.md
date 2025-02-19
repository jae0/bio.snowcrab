---
title: "Snow Crab, Scotian Shelf, Canada (NAFO Div. 4X)"
subtitle: "Science Assessment, Bedford Institute of Oceanography, Fisheries and Oceans Canada"
keywords: 
  - snow crab fishery assessment -- 4X focus
  - areal-unit modelling of numerical abundance, habitat, mean weight
  - convolution of posterior-simulations of biomass   
abstract: |
  Snow crab modelling of numerical abundance, mean weight and habitat using conditional autoregressive models. Biomass is estimated from these components via multiplication of the posterior draws. 4X subset update.
fontsize: 12pt
metadata-files:
  - _metadata.yml
params:
  year_assessment: 2024
  year_start: 1999
  data_loc:  "~/bio.data/bio.snowcrab"
  sens: 1
  debugging: FALSE
  model_variation: logistic_discrete_historical
--- 


<!-- 

make quarto FN=snowcrab_presentation_4x.md YR=2024 DATADIR=~/bio.data/bio.snowcrab DOCTYPE=revealjs  PARAMS="-P year_assessment:2024"  --directory=~/bio/bio.snowcrab/inst/markdown 

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

  # things to load into memory (in next step) via _load_results.qmd
  toget = c( "fishery_results", "fishery_model", "ecosystem" )  

```


<!-- 
# _load_results.qmd contains instructions to load data 
#  this is a shared R-script to boot strap and provide a consistent data interface
-->

{{< include _load_results.qmd >}}  

 

 
# Ecosystem considerations

## Connectivity: Oceanic currents {.c}

::: columns 
:::: column 

```{r ocean_currents, echo=FALSE, out.width='55%', fig.align='center', fig.show='hold',  fig.cap = 'Ocean currents in the Martimes. Source: https://www.dfo-mpo.gc.ca/oceans/publications/soto-rceo/2018/atlantic-ecosystems-ecosystemes-atlantiques/index-eng.html' }
fn2=file.path( media_loc, "maritimes_currents.png" )
knitr::include_graphics( c(fn2) ) 
# \@ref(fig:currents)  
``` 
::::
 
:::: column

- Variable due to confluence:

  - Warm, high salinity Gulf Stream from the S-SE along the shelf edge 

  - Cold, low salinity Labrador Current

  - Cold low salinity St. Lawrence outflow from the N-NE

  - Nearshore Nova Scotia current, running from the NE. 
 
::::
:::

 


## Bathymetry {.c}


```{r bathymetry-map, out.width='55%', echo=FALSE, fig.align='center', fig.cap='Variations in log(depth; m) in the Scotian Shelf region.' }
bathydir = file.path( data_root, 'aegis', 'bathymetry', 'modelled', 'default', 'stmv', 'none_fft', 'z', 'maps', 'SSE' )
knitr::include_graphics( file.path( bathydir, 'bathymetry.z.SSE.png' ) )
```




## Substrate {.c}


```{r substrate-map, out.width='55%', echo=FALSE, fig.align='center', fig.cap='Substrate grain size log(mm) variations in the Scotian Shelf region.' }
substrdir = file.path( data_root, 'aegis', 'substrate', 'maps', 'canada.east.highres' )
knitr::include_graphics( file.path(  substrdir, 'substrate.substrate.grainsize.canada.east.highres.png' ) )
```


## Bottom Temperature {.c}

\vspace{12mm}
```{r bottom-temperatures-survey, out.width='55%', echo=FALSE, fig.align='center', fig.cap = 'Annual variations in bottom temperature observed during the Snow Crab survey. The horizontal (black) line indicates the long-term, median temperature within each subarea. Error bars represent standard errors.' }
knitr::include_graphics( file.path( data_loc, 'assessments', year_assessment, 'timeseries', 'survey', 't.png') )
# \@ref(fig:bottom-temperatures-survey)
```


## Bottom Temperature ... {.c}

```{r bottom-temperatures, out.width='35%', echo=FALSE, fig.align='center', fig.cap = 'Posterior densities of predicted average bottom temperatures. Red horizontal line is at $7^\\circ$C.' }
knitr::include_graphics( file.path( data_loc, 'assessments', year_assessment, 'timeseries', 'temperature_bottom.png') )
# \@ref(fig:bottom-temperatures)
``` 

## Bottom Temperature ... {.c}

```{r bottom-temperatures-map, out.width='48%', echo=FALSE, fig.show='hold', fig.align='center', fig.cap = 'Spatial variations in predicted (1 September) bottom temperature from 2021 (left) to 2023 (right). Groundfish surveys were not conducted in 2020 and 2022 in the snow crab domain. Snow crab surveys were not conducted in 2020, and incomplete in 2022 in S-ENS.' }

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
knitr::include_graphics( c( fn2, fn1) )
# \@ref(fig:bottom-temperatures-map)
# *Spatial variations in bottom temperature estimated from a historical analysis of temperature data for 1 September.*
```

 





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
include_graphics( file.path( data_loc, 'output', 'bcd.png') )
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
Stomach sampling program

```{r predators2, echo=FALSE} 
gt( counts[1:11,] )  |> gt::tab_options(table.font.size = 10, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
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

```{r cod-map, out.width='48%', fig.show='hold', fig.align='center', fig.cap= 'Atlantic cod density log10(no/km$^2$) from the Snow Crab survey. Predation mortality more likely higher encounter rates (high densities co-occur).' }
loc = file.path( data_loc, 'output', 'maps', 'survey', 'snowcrab', 'annual', 'bycatch', 'ms.no.10' )
yrsplot = setdiff( year_assessment + c(0:-9), 2020)
fn4 = file.path( loc, paste( 'ms.no.10', yrsplot[4], 'png', sep='.') )
fn3 = file.path( loc, paste( 'ms.no.10', yrsplot[3], 'png', sep='.') )
fn2 = file.path( loc, paste( 'ms.no.10', yrsplot[2], 'png', sep='.') )
fn1 = file.path( loc, paste( 'ms.no.10', yrsplot[1], 'png', sep='.') )
include_graphics( c(  fn2, fn1) )
# \@ref(fig:cod-map)  
```
 


## Predators - Atlantic cod {.c}

```{r cod-timeseries, out.width='60%', echo=FALSE,  fig.align='center', fig.cap = 'Atlantic cod crude, unadjusted geometric mean numerical density (no/km$^2$) from annual Snow Crab survey. Error bars are 95\\%  Confidence Intervals.'}
include_graphics( file.path( data_loc, 'assessments', year_assessment, 'timeseries', 'survey', 'ms.no.10.png') )
# \@ref(fig:cod-timeseries)
```


 
##  Predators - Atlantic Halibut  {.c}

```{r halibut-map, out.width='48%', fig.show='hold', fig.align='center', fig.cap= 'Halibut density log10(no/km$^2$) from the Snow Crab survey. Predation mortality more likely higher encounter rates (high densities co-occur).' }
loc = file.path( data_loc, 'output', 'maps', 'survey', 'snowcrab', 'annual', 'bycatch', 'ms.no.30' )
yrsplot = setdiff( year_assessment + c(0:-9), 2020)
fn4 = file.path( loc, paste( 'ms.no.30', yrsplot[4], 'png', sep='.') )
fn3 = file.path( loc, paste( 'ms.no.30', yrsplot[3], 'png', sep='.') )
fn2 = file.path( loc, paste( 'ms.no.30', yrsplot[2], 'png', sep='.') )
fn1 = file.path( loc, paste( 'ms.no.30', yrsplot[1], 'png', sep='.') )
include_graphics( c( fn2, fn1) )
# \@ref(fig:halibut-map)  
```
 

##  Predators - Atlantic Halibut ... {.c}

```{r halibut-timeseries, out.width='50%', echo=FALSE,   fig.align='center', fig.cap = 'Atlantic Halibut crude, unadjusted geometric mean numerical density (no/km$^2$) from annual Snow Crab survey. Error bars are 95\\%  Confidence Intervals.' }
include_graphics( file.path( data_loc, 'assessments', year_assessment, 'timeseries', 'survey', 'ms.no.30.png') )
# \@ref(fig:halibut-timeseries)
```


<!--
~~~{=comment}


## Predators - Thorny skate {.c}

```{r thornyskate-timeseries, out.width='60%', echo=FALSE,  fig.align='center', fig.cap = 'Thorny Skate crude, unadjusted geometric mean numerical density (no/km$^2$) from annual Snow Crab survey. Error bars are 95\\%  Confidence Intervals.'}
include_graphics( file.path( data_loc, 'assessments', year_assessment, 'timeseries', 'survey', 'ms.no.201.png') )
# \@ref(fig:thornyskate-timeseries)
```

## Predators - Thorny skate ... {.c}

```{r thornyskate-map, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Thorny skate density log10(no/km$^2$) from the Snow Crab survey.' }
loc = file.path( data_loc, 'output', 'maps', 'survey', 'snowcrab', 'annual', 'bycatch', 'ms.no.201' )
yrsplot = setdiff( year_assessment + c(0:-9), 2020)
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
include_graphics( file.path( data_loc, 'assessments', year_assessment, 'timeseries', 'survey', 'ms.no.50.png') )
# \@ref(fig:Wolffish-timeseries)
```

##  Predators - Striped Atlantic Wolffish  ... {.c}


```{r Wolffish-map, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Striped Atlantic Wolffish density log10(no/km$^2$) from the Snow Crab survey.' }
loc = file.path( data_loc, 'output', 'maps', 'survey', 'snowcrab', 'annual', 'bycatch', 'ms.no.50' )
yrsplot = setdiff( year_assessment + c(0:-9), 2020)
fn4 = file.path( loc, paste( 'ms.no.50', yrsplot[4], 'png', sep='.') )
fn3 = file.path( loc, paste( 'ms.no.50', yrsplot[3], 'png', sep='.') )
fn2 = file.path( loc, paste( 'ms.no.50', yrsplot[2], 'png', sep='.') )
fn1 = file.path( loc, paste( 'ms.no.50', yrsplot[1], 'png', sep='.') )
include_graphics( c( fn3, fn2, fn1) )
# \@ref(fig:Wolffish-map)  
```


~~~
-->  



<!-- 
## Competitors/Prey - Lesser toad crab ... {.c}


```{r lessertoadcrab-map, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Lesser Toad Crab density log10(no/km$^2$) from the Snow Crab survey.' }
loc = file.path( data_loc, 'output', 'maps', 'survey', 'snowcrab', 'annual', 'bycatch', 'ms.no.2521' )
yrsplot = setdiff( year_assessment + c(0:-9), 2020)
fn4 = file.path( loc, paste( 'ms.no.2521', yrsplot[4], 'png', sep='.') )
fn3 = file.path( loc, paste( 'ms.no.2521', yrsplot[3], 'png', sep='.') )
fn2 = file.path( loc, paste( 'ms.no.2521', yrsplot[2], 'png', sep='.') )
fn1 = file.path( loc, paste( 'ms.no.2521', yrsplot[1], 'png', sep='.') )
include_graphics( c( fn3, fn2, fn1) )
# \@ref(fig:lessertoadcrab-map)  
```

## Competitor/Prey - Lesser toad crab    {.c}


```{r lessertoadcrab-timeseries, out.width='60%', echo=FALSE,   fig.align='center', fig.cap = 'Lesser Toad Crab crude, unadjusted geometric mean numerical density (no/km$^2$) from annual Snow Crab survey. Error bars are 95\\%  Confidence Intervals.' }
include_graphics( file.path( data_loc, 'assessments', year_assessment, 'timeseries', 'survey', 'ms.no.2521.png') )
# \@ref(fig:lessertoadcrab-timeseries)
```
 

## Competitor/Prey - Northern shrimp {.c}

```{r Shrimp-map, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Northern Shrimp density log10(no/km$^2$) from the Snow Crab survey. Shrimp with similar habitat preferences have declined, possibly due to large-scaled habitat variations and predation. Sampling was incomplete in 2020 and 2022 in S-ENS.' }
loc = file.path( data_loc, 'output', 'maps', 'survey', 'snowcrab', 'annual', 'bycatch', 'ms.no.2211' )
yrsplot = setdiff( year_assessment + c(0:-9), 2020)
fn4 = file.path( loc, paste( 'ms.no.2211', yrsplot[4], 'png', sep='.') )
fn3 = file.path( loc, paste( 'ms.no.2211', yrsplot[3], 'png', sep='.') )
fn2 = file.path( loc, paste( 'ms.no.2211', yrsplot[2], 'png', sep='.') )
fn1 = file.path( loc, paste( 'ms.no.2211', yrsplot[1], 'png', sep='.') )
include_graphics( c( fn3, fn2, fn1) )
# \@ref(fig:Shrimp-map)  
```


## Competitor/Prey - Northern shrimp ...  {.c}


```{r Shrimp-timeseries, out.width='60%', echo=FALSE,   fig.align='center', fig.cap = 'Northern Shrimp crude, unadjusted geometric mean numerical density (n/$km^2$) from annual Snow Crab survey. Error bars are 95\\%  Confidence Intervals.' }
include_graphics( file.path( data_loc, 'assessments', year_assessment, 'timeseries', 'survey', 'ms.no.2211.png') )
# \@ref(fig:Shrimp-timeseries)
```
  
-->

# Fishery performance


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


<!-- 
## Bycatch Maritimes {.c}
 
```{r bycatch-cpue, out.width='70%', echo=FALSE, fig.align='center', fig.cap = 'Bycatch rates in Maritimes Region. Low levels attributable to trap design (top entry, conical, large mesh 5.25" knot-to-knot) permits escapement of non-target species.' }
o = o_cfaall$bycatch_table_catch
o[ o==0 ] = NA
o[ is.na(o) ] = "."
gt(o)
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

-->

    
<!-- 
    
## Bycatch N-ENS {.c}
 
 
```{r bycatch-cpue-n-ens, out.width='12%', echo=FALSE, fig.align='center', fig.cap = 'Bycatch (kg) in N-ENS for those greater than 10 kg/year.' }
o = o_cfanorth$bycatch_table_catch

o = o[ which(o$"Average/Moyen" > 10 ),]
o$"Average/Moyen" = round(o$"Average/Moyen")
o[ o==0 ] = NA
o[ is.na(o) ] = "."
gt(o)
```  

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
 
 
```{r bycatch-cpue-s-ens, out.width='12%', echo=FALSE, fig.align='center', fig.cap = 'Bycatch (kg) in S-ENS for those greater than 10 kg/year' }
o = o_cfasouth$bycatch_table_catch    
o = o[ which(o$"Average/Moyen" > 10 ),]
o$"Average/Moyen" = round(o$"Average/Moyen")
o[ o==0 ] = NA
o[ is.na(o) ] = "."
gt(o)
```  


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
-->


## Bycatch 4X {.c}
Bycatch (kg) for those greater than 10 kg/year.
```{r bycatch-cpue-4x,  echo=FALSE, fig.align='center', fig.cap = 'Bycatch (kg) in 4X for those greater than 10 kg/year.' }

o = o_cfa4x$bycatch_table_catch 
o = o[ which(as.numeric(o$"Average/Moyen") > 10 ),]
o$"Average/Moyen" = round(o$"Average/Moyen")
o[ o==0 ] = NA
o[ is.na(o) ] = "."
gt(o)  |> gt::tab_options(table.font.size = 6, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
```  
 


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


## Fishery performance indices 

4X:

```{r}
#| label: table-fishery-4x-perf
#| tbl-cap: "Fishery performance statistics in 4X. Units are: TACs and Landings (tons, t), Effort ($\\times 10^3$ trap hauls, th) and CPUE (kg/th). There were no landings or TACs in 2018/2019 due to indications of low abundance. The 2022 season is ongoing."
#| eval: true
#| output: true

ii = which(dt$Region=="cfa4x")
oo = dt[ii, c("Year", "Licenses", "TAC", "Landings", "Effort", "CPUE")] 
gt::gt(oo) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
```
 

## At sea observed data: Carapace condition from observed data < 95mm CW

```{r}
#| label: setup-observer-data
#| eval: true
#| output: false
  odb = odb0[ cw < 95 & prodcd_id==0 & shell %in% c(1:5) & region %in% regions & sex==0, ]  # male
```
 
4X:

```{r}
#| eval: true
#| output: true
#| label: table-fishery-4x-sublegal
#| tbl-cap: "Fishery performance statistics in 4X. Distribution of at sea observations of males less than 95 mm CW by year and shell condition."
resX = dcast( odb0[ region=="cfa4x", .(N=.N), by=.(fishyr, shell) ], fishyr  ~ shell, value.var="N", fill=0, drop=FALSE, na.rm=TRUE )
if ("NA" %in% names(resX)) resX$"NA" = NULL
names(resX) = c("Year", "CC1", "CC2", "CC3", "CC4", "CC5" )
resX$Total = rowSums( resX[, 2:6 ], na.rm=TRUE)
resX[, 2:6 ] = round(resX[, 2:6 ] / resX$Total * 100, digits=2)
gt::gt(resX) |> gt::tab_options(table.font.size = 10, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
```


## At sea observed data: Carapace condition from observed data  >= 95mm CW

```{r}
#| eval: true
#| output: false
odb = odb0[ cw >= 95 & cw < 170  & prodcd_id==0 & shell %in% c(1:5) & region %in% regions & sex==0, ]  # male
```
 
4X:

```{r}
#| eval: true
#| output: true
#| label: table-fishery-4x-comm
#| tbl-cap: "Fishery performance statistics in 4X. Distribution of at sea observations of males greater than 95 mm CW by year and shell condition."
resX = dcast( odb[ region=="cfa4x", .(N=.N), by=.(fishyr, shell) ], fishyr  ~ shell, value.var="N", fill=0, drop=FALSE, na.rm=TRUE )
names(resX) = c("Year", "CC1", "CC2", "CC3", "CC4", "CC5" )
resX$Total = rowSums( resX[, 2:6 ], na.rm=TRUE)
resX[, 2:6 ] = round(resX[, 2:6 ] / resX$Total * 100, digits=2)
gt::gt(resX) |> gt::tab_options(table.font.size = 10, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
```

   
   
##  Percent soft from observed data

There are two possible definitions: 

- durometer < 68 (Soft, Total)
- carapace conditions 1 and 2 (SoftSC, TotalSC)


```{r}
#| eval: true
#| output: true
odb = odb0[ cw >= 95 & cw < 170  & prodcd_id==0 & shell %in% c(1:5) & region %in% regions & sex==0, ]  # male
shell_condition = odb[ !is.na(odb$region), .N, by=.(region, fishyr, shell) ]
shell_condition[, total:=sum(N, na.rm=TRUE), by=.(region, fishyr)]
shell_condition$percent = round(shell_condition$N / shell_condition$total, 3) * 100
shell_condition$Year = shell_condition$fishyr
```
 
4X:

```{r}
#| eval: true
#| output: true
#| label: table-fishery-4x-soft-durometer
#| tbl-cap: "Fishery performance statistics in 4X. Distribution of at sea observations of males soft-shelled  based on durometer (<68) and shell condition (1 and 2, SC)."
softX  = odb[ region=="cfa4x" & durometer <  68, .(Soft=.N), by=.(fishyr ) ] 
totalX = odb[ region=="cfa4x" & is.finite(durometer) , .(Total=.N), by=.(fishyr) ] 
resX = softX[totalX, on="fishyr"]
resX = resX[, .(Year=fishyr, Soft=round(Soft/Total*100,2), Total=Total) ]  
ssX = shell_condition[ region=="cfa4x" & shell %in% c(1,2), .(SoftSC=sum(percent), TotalSC=unique(total)[1]), by=.(Year)]
resX = resX[ssX, on="Year"]
gt::gt(resX) |> gt::tab_options(table.font.size = 8, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
```

   
    
  

<!--
# instars of interest: 11 and 12

# growth increment (assumming average weight in the midpoint of each increment)
growth.11.to.12 =  predict.mass.g.from.CW.mm( mean(CW.interval.male(12)) ) - predict.mass.g.from.CW.mm (mean(CW.interval.male(11)) )
(growth.11.to.12)
# = 419 g
#  12to13 = ~450
-->


# Survey results and stock status {.c}

  
##  Carapace condition from trawl data  >= 95mm CW  

```{r}
#| eval: true
#| output: false
det = snowcrab.db( p=p, DS="det.georeferenced" )
setDT(det)
det$fishyr = det$yr  ## the counting routine expectes this variable
det = det[ cw >= 95 ,]  # commerical sized crab only
years = sort( unique( det$yr ) )
det$region = NA
for ( reg in regions) {
  r = polygon_inside(x = det, region = aegis.polygons::polygon_internal_code(reg), planar=FALSE)
  det$region[r] = reg
}

```
 
4X:

```{r}
#| eval: true
#| output: true
#| label: table-survey-4X-comm
#| tbl-cap: "Distribution of 4X survey: males less than 95 mm CW by year and shell condition."
resX = dcast( det[ region=="cfa4x" & !is.na(shell), .(N=.N), by=.(fishyr, shell) ], fishyr  ~ shell, value.var="N", fill=0, drop=FALSE, na.rm=TRUE )
names(resX) = c("Year", "CC1", "CC2", "CC3", "CC4", "CC5" )
resX$Total = rowSums( resX[, 2:6 ], na.rm=TRUE)
resX[, 2:6 ] = round(resX[, 2:6 ] / resX$Total * 100, digits=2)
gt::gt(resX) |> gt::tab_options(table.font.size = 8, data_row.padding = gt::px(1), 
  summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
  footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
  row_group.padding = gt::px(1))
```


## Counts of stations in each area

```{r}
#| eval: true
#| output: true
#| label: table-survey-station-count
#| tbl-cap: "Survey station counts"
set = snowcrab.db(p=p, DS="set.clean")
setDT(set)
# check towquality .. this should always == 1
if (length( unique( set$towquality) ) != 1 ) print("error -- not good tows")
set$region = NA
for (reg in c( "cfanorth", "cfasouth", "cfa4x"  ) ) {
  d = polygon_inside(set[,c("lon","lat")], reg)
  set$region[d] = reg 
}
out = dcast( set[, .(N=.N), by=.(region, yr)], yr~region, value.var="N", fill=0, drop=FALSE, na.rm=TRUE )
out[,Total:=sum(cfanorth,cfasouth,cfa4x, na.rm=TRUE)]
out = out[, .(yr, cfanorth, cfasouth, cfa4x)]
names(out) = c("Year", "NENS", "SENS", "4X")
gt::gt(out) |> gt::tab_options(table.font.size = 10, data_row.padding = gt::px(1), 
  summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
  footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
  row_group.padding = gt::px(1))
```


## Connectivity: Movement {.c}

 
```{r movementtracks, echo=FALSE, out.width='65%', fig.align='center', fig.show='hold',  fig.cap = 'Snow Crab movement tracks from mark and recapture.' }
fn1=file.path( media_loc, "movement0.png" )
fn2=file.path( media_loc, "movement.png" )
knitr::include_graphics( c(fn1) ) 
# \@ref(fig:movementtracks)  
``` 

## Connectivity: Movement ... {.c}

 
```{r movementtracks2, echo=FALSE, out.width='71%', fig.align='center', fig.show='hold',  fig.cap = 'Snow Crab movement tracks from mark and recapture.' }
fn1=file.path( media_loc, "movement0.png" )
fn2=file.path( media_loc, "movement.png" )
knitr::include_graphics( c(fn2) ) 
# \@ref(fig:movementtracks2)  
``` 
 
## Connectivity: Movement {.c}
 

```{r movement_tr1, echo=FALSE, out.width='45%', fig.align='center', fig.show='hold',  fig.cap = 'Snow Crab movement distances and rates.' }
fn1=file.path( media_loc, "snowcrab_movement_distances.png" )
fn2=file.path( media_loc, "snowcrab_movement_rates.png" )
knitr::include_graphics( c(fn1) ) 
# \@ref(fig:movement_tr1)  
``` 
 

## Connectivity: Movement 2 {.c}
 

```{r movement_tr2, echo=FALSE, out.width='45%', fig.align='center', fig.show='hold',  fig.cap = 'Snow Crab movement distances and rates.' }
fn1=file.path( media_loc, "snowcrab_movement_distances.png" )
fn2=file.path( media_loc, "snowcrab_movement_rates.png" )
knitr::include_graphics( c(fn2) ) 
# \@ref(fig:movement_tr2)  
``` 
 





## Female Geometric mean density by size and maturity
    
```{r sizefeq-female, out.width='43%', echo=FALSE, fig.align='center', fig.cap = 'Geometric mean density distribution (no/km$^2$) by carapace width of immature (light) and mature (dark).'}
include_graphics(  file.path( data_loc, "assessments", year_assessment, "figures", "size.freq", "survey", "female.denl.png" )  )
# \@ref(fig:sizefeq-male)
```  

## Mature female maps

Distributions are heterogeneous and often in shallower areas.  

```{r fmat-map, echo=FALSE, out.width='48%', fig.show='hold', fig.align='center', fig.cap= 'Mature female density log$_{10}$(no/km$^2$)  from the Snow Crab survey.' }
loc = file.path( data_loc, "output", "maps", "survey", "snowcrab", "annual", "totno.female.mat" )
yrsplot = setdiff( year_assessment + c(0:-9), 2020)
fn4 = file.path( loc, paste( "totno.female.mat", yrsplot[4], "png", sep=".") )
fn3 = file.path( loc, paste( "totno.female.mat", yrsplot[3], "png", sep=".") )
fn2 = file.path( loc, paste( "totno.female.mat", yrsplot[2], "png", sep=".") )
fn1 = file.path( loc, paste( "totno.female.mat", yrsplot[1], "png", sep=".") )
include_graphics( c(  fn2, fn1) )
# \@ref(fig:fmat-map)  
```

## Mature female timeseries

```{r fmat-timeseries, out.width='75%', echo=FALSE, fig.align='center', fig.cap = 'Mature female density log$_{10}$(no/km$^2$) from the Snow Crab survey.'  }
include_graphics( file.path( data_loc, "assessments", year_assessment, "timeseries", "survey", "totno.female.mat.png") )
# \@ref(fig:fmat-timeseries)
```

## Male Geometric mean density by size and maturity
    
```{r sizefeq-male, out.width='43%', echo=FALSE, fig.align='center', fig.cap = 'Geometric mean density (no/km$^2$) by carapace width of immature (light) and mature (dark) male snow crab. The vertical line is legal size (95 mm).   '}
include_graphics(  file.path( data_loc, "assessments", year_assessment, "figures", "size.freq", "survey", "male.denl.png" )  )
# \@ref(fig:sizefeq-male)
```  

## Biomass Density  {.c}
 
\begin{tiny}
```{r fbgeomean-map, echo=FALSE, out.width='48%', fig.show='hold', fig.align='center', fig.cap= 'Snow Crab survey fishable component biomass density log~10(t/km$^2$). Note, there is no data in 2020. Note that high and low biomass density areas fluctuate with time.' }
loc = file.path( data_loc, 'output', 'maps', 'survey', 'snowcrab', 'annual', 'R0.mass')
yrsplot =  setdiff(year_assessment + c(0:-9), 2020 ) 
#fn6 = file.path( loc, paste( 'R0.mass', yrsplot[6], 'png', sep='.') )
fn5 = file.path( loc, paste( 'R0.mass', yrsplot[5], 'png', sep='.') )
fn4 = file.path( loc, paste( 'R0.mass', yrsplot[4], 'png', sep='.') )
fn3 = file.path( loc, paste( 'R0.mass', yrsplot[3], 'png', sep='.') )
fn2 = file.path( loc, paste( 'R0.mass', yrsplot[2], 'png', sep='.') )
fn1 = file.path( loc, paste( 'R0.mass', yrsplot[1], 'png', sep='.') )
include_graphics( c( fn2, fn1) )
``` 
\end{tiny}
 
    

## Biomass Density ... {.c}
 

\begin{tiny}

```{r fbGMTS, out.width='65%', echo=FALSE, fig.align='center', fig.cap = 'The crude, unadjusted geometric mean fishable biomass density log~10(t/km$^2$) from the Snow Crab survey. Error bars represent 95\\% Confidence Intervals. Note the absence of data in 2020. Prior to 2004, surveys were conducted in the Spring. A peak in 2009 to 2014 and has since been declining in all areas. '}
fn = file.path(data_loc, 'assessments', year_assessment, 'timeseries','survey','R0.mass.png')
include_graphics( c(fn) )
#\@ref(fig:fbGMTS)
```
\end{tiny}

 

## Viable Habitat ... {.c}
  
```{r fb-habitat-map, out.width='48%', fig.show='hold', fig.align='center', fig.cap= 'Habitat viability (probability; fishable Snow Crab).' }
loc = file.path( data_loc, 'modelled', 'default_fb', 'predicted_habitat' )
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
include_graphics( c(  fn2, fn1) )
# \@ref(fig:fb-habitat-map)  
# *Figure XXX. Habitat viability (probability; fishable Snow Crab)* 
```


## Viable Habitat ... {.c}

```{r fb-habitat-timeseries, out.width='70%', echo=FALSE,   fig.align='center', fig.cap = 'Habitat viability (probability; fishable Snow Crab). Means and 95\\% Credible Intervals are presented.' }
loc = file.path( data_loc, 'modelled', 'default_fb', 'aggregated_habitat_timeseries' )
include_graphics( file.path( loc, 'habitat_M0.png') )
# \@ref(fig:fb-habitat-timeseries)
```
 


    
## Biomass Index (aggregate)
 
```{r fbindex-map, echo=FALSE, out.width='48%', fig.show='hold', fig.align='center', fig.cap= 'Biomass index log$_{10}$(t/km$^2$) predicted from the Snow Crab survey.' }
loc = file.path( data_loc, 'modelled', 'default_fb', 'predicted_biomass_densities' )
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
include_graphics( c(   fn2, fn1) )
```




## Biomass Index (aggregate) ... {.c}
    
\begin{tiny}
```{r fbindex-timeseries, out.width='65%', echo=FALSE, fig.align='center', fig.cap = 'The fishable biomass index (t) predicted by CARSTM of Snow Crab survey densities. Error bars represent Bayesian 95\\% Credible Intervals. Note large errors in 2020 when there was no survey.' }
include_graphics( file.path( data_loc, 'modelled', 'default_fb', 'aggregated_biomass_timeseries' , 'biomass_M0.png') )
# \@ref(fig:fbindex-timeseries)
```
\end{tiny}




## Fishery Model Biomass (pre-fishery){.c}

\begin{tiny}
```{r logisticPredictions, out.width='45%', echo=FALSE, fig.show='hold', fig.align='center', fig.cap = 'Model 1 fishable, posterior mean modelled biomass (pre-fishery; kt) are shown in dark orange for N-ENS, S-ENS and 4X (left, middle and right). Light orange are posterior samples of modelled biomass (pre-fishery; kt) to illustrate the variability of the predictions. The biomass index (post-fishery, except prior to 2004) after model adjustment by the model catchability coefficient is in gray.' } 
loc = file.path( data_loc, 'fishery_model', year_assessment, 'logistic_discrete_historical' )
fn1 = file.path( loc, 'plot_predictions_cfanorth.png' ) 
fn2 = file.path( loc, 'plot_predictions_cfasouth.png' ) 
fn3 = file.path( loc, 'plot_predictions_cfa4x.png' ) 
include_graphics(c(fn3) )
# \@ref(fig:logisticPredictions)
```
\end{tiny}




## Fishing Mortality {.c}

```{r logisticFishingMortality, out.width='45%', echo=FALSE,  fig.show='hold', fig.align='center', fig.cap = 'Time-series of modelled instantaneous fishing mortality from Model 1, for N-ENS (left), S-ENS (middle), and 4X (right). Samples of the posterior densities are presented, with the darkest line being the mean.' }
  odir = file.path( fishery_model_results, year_assessment, "logistic_discrete_historical" )
  fn1 = file.path( odir, "plot_fishing_mortality_cfanorth.png" ) 
  fn2 = file.path( odir, "plot_fishing_mortality_cfasouth.png" ) 
  fn3 = file.path( odir, "plot_fishing_mortality_cfa4x.png" ) 
include_graphics(c( fn3) )
# \@ref(fig:logisticFishingMortality)
```


 
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
 

```{r ReferencePoints, out.width='60%', echo=FALSE, fig.align='center', fig.cap = 'Harvest control rules for the Scotian Shelf Snow Crab fisheries.' }
include_graphics( file.path( media_loc, 'harvest_control_rules.png') ) 
# \@ref(fig:ReferencePoints)
```


## Reference Points ... {.c}

```{r logistic-hcr, out.width='45%', echo=FALSE, fig.show='hold', fig.align='center', fig.cap = 'Reference Points (fishing mortality and modelled biomass) from the Fishery Model, for N-ENS (left), S-ENS (middle), and 4X (right). The large yellow dot indicates most recent year and the 95\\% CI. Not: the model does not account for illegal and unreported landings, and interspecific interactions.' }
  odir = file.path( fishery_model_results, year_assessment, "logistic_discrete_historical" )
  fn1 = file.path( odir, 'plot_hcr_cfanorth.png' ) 
  fn2 = file.path( odir, 'plot_hcr_cfasouth.png' ) 
  fn3 = file.path( odir, 'plot_hcr_cfa4x.png' ) 
  include_graphics(c( fn3) )
#  \@ref(fig:logistic-hcr)
```

4X is in the "critical" zone
  
 

## Conclusions: general

The ESS ecosystem is still experiencing a lot of volatility and prudence is wise:

- Rapid ecosystem change (groundfish collapse in 1980s and 1990s)
- Climatic variability: very warm period from 2012-2023
- Snow crab: 
  - Can persist if extreme conditions are episodic by shifting spatial distributions    
  - 4X was possibly too long and did not have the habitat space
  - Climate variability can induce juxtaposition of warmer water species (prey and predators) with snow crab
  

<!-- 
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

-->


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

\begin{tiny}


Banerjee, S., Carlin, B. P., and Gelfand, A. E.. 2004. Hierarchical Modeling and Analysis for Spatial Data. Monographs on Statistics and Applied Probability. Chapman and Hall/CRC.

Besag, Julian. 1974. Spatial interaction and the statistical analysis of lattice systems. Journal of the Royal Statistical Society Series B (Methodological) 1974: 192-236.

Canada Gazette. 2022. [Regulations Amending the Fishery (General) Regulations. Part II, Volume 156, Number 8.](https://www.gazette.gc.ca/rp-pr/p2/2022/2022-04-13/html/sor-dors73-eng.html)

Canada Gazette. 2016. St. Anns Bank Marine Protected Area Regulations. Canada Gazette, Part I, Vol 150, Issue 51: 4143-4149.

Choi, J.S. 2020. A Framework for the assessment of Snow Crab (*Chioneocete opilio*) in Maritimes Region (NAFO Div 4VWX) . DFO Can. Sci. Advis. Sec. Res. Doc. 2020/nnn. v + xxx p.

Choi, J.S. 2022. Reconstructing the Decline of Atlantic Cod with the Help of Environmental Variability in the Scotian Shelf of Canada. bioRxiv. https://doi.org/10.1101/2022.05.05.490753.

Choi, J. S., and B. C. Patten. 2001. Sustainable Development: Lessons from the Paradox of Enrichment. Ecosystem Health 7: 163–77.

Choi, Jae S., B. Cameron, K. Christie, A. Glass, and E. MacEachern. 2022. Temperature and Depth Dependence of the Spatial Distribution of Snow Crab. bioRxiv. https://doi.org/10.1101/2022.12.20.520893.


Choi, Jae S. 2023. A Multi-Stage, Delay Differential Model of Snow Crab Population Dynamics in the Scotian Shelf of Atlantic Canada. bioRxiv. https://doi.org/10.1101/2023.02.13.528296.

\end{tiny}


## References ...

\begin{tiny}

DFO. 2013. [Integrated Fisheries Management Plan for Eastern Nova Scotia and 4X Snow Crab (*Chionoecetes Opillio*.)](http://www.dfo-mpo.gc.ca/fm-gp/peches-fisheries/ifmp-gmp/snow-crab-neige/snow-crab-neiges2013-eng.htm)


DFO. 2018. Stock Status Update of Atlantic Halibut (Hippoglossus hippoglossus) on the Scotian Shelf and Southern Grand Banks in NAFO Divisions 3NOPs4VWX5Zc. DFO Can. Sci. Advis. Sec. Sci. Resp. 2018/022.

Hebert M, Miron G, Moriyasu M, Vienneau R, and DeGrace P. Efficiency and ghost fishing of Snow Crab (Chionoecetes opilio) traps in the Gulf of St. Lawrence. Fish Res. 2001; 52(3): 143-153. 10.1016/S0165-7836(00)00259-9   
 
Riebler, A., Sørbye, S.H., Simpson D., and Rue, H. 2016. An intuitive Bayesian spatial model for disease mapping that accounts for scaling. Statistical methods in medical research 25: 1145-1165.

Simpson, D., Rue, H., Riebler, A., Martins, T.G., and Sørbye, SH. 2017. Penalising Model Component Complexity: A Principled, Practical Approach to Constructing Priors. Statist. Sci. 32: 1-28.


\end{tiny}


# End

 

## Supplemental Information 
     