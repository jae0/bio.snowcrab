---
title: "Snow Crab fishery -- figures and tables"
author:
  - name: Snow crab group
    # orcid: 0000-0003-3632-5723 
    # email: jae.choi@dfo-mpo.gc.ca
    # email: choi.jae.seok@gmail.com
    # corresponding: true
    affiliation: 
      - name: Bedford Institute of Oceanography, Fisheries and Oceans Canada
        city: Dartmouth
        state: NS
        url: www.bio.gc.ca
date: 2024-08-17
keywords: 
  - snow crab trawl survey results
  - basic tables and figures
abstract: |
  Snow crab fishery performance. 
toc: true
number-sections: true
highlight-style: pygments
bibliography: media/references.bib  
# csl: media/canadian-journal-of-fisheries-and-aquatic-sciences.csl  # see https://www.zotero.org/styles for more
license: "CC BY"
copyright: 
  holder: Jae S. Choi
  year: 2024
citation: 
  container-title: https://github.com/jae0/bio.snowcrab/
  doi: NA
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
  pdf:
    pdf-engine: lualatex
  docx: default 
---

 

# Snow Crab fishery -- figures and tables
 
## Preamble

This is a markdown document. It can be viewed in formatted form via:
  
  - a web browser open the webpage: [02_fishery_summary.md](https://github.com/jae0/bio.snowcrab/tree/master/inst/markdown/02_fishery_summary.md)   

  - a web browser open the local file directly: [02_fishery_summary.md](../markdown/02_fishery_summary.md) (you might need to install a browser add-in), or 
  
  - a text editor with markdown awareness (e.g. Rstudio, VSCode, etc. ).


As this document uses the Quarto dialect of Markdown, you can easily create a report, by running one of the following in a shell command prompt or using Rstudio: 

```shell
# {via Quarto}
make quarto FN=02_fishery_summary YR=2024 SOURCE=~/bio/bio.snowcrab/inst/markdown WK=~/bio.data/bio.snowcrab/assessments  DOCEXTENSION=html 
 
```

Or, see the Makefile and alter defaults to your needs. As Quarto does not pass params easily. So you must adjust "params" in yaml at the top of this fle, or use to quarto command such as: 

```shell
quarto ... -P year.assessment:$(YR) -P media_loc:$(MEDIA) 
```

[See here for more YAML options.](https://quarto.org/docs/output-formats/all-formats.html)
 


## Set up environment

  - Ensure the **year.assessment** is correct

  - Ensure that the data pulls of survey data stored on the ISSDB is complete and assimilated to the end of [01_snowcrab_data.md](https://github.com/jae0/bio.snowcrab/tree/master/inst/markdown/01_snowcrab_data.md).
  



```{r}
#| label: setup
#| eval: true 
#| output: false

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
require(aegis)  # basic helper tools

year.assessment = 2024  # change this as appropriate
year_previous = year.assessment - 1

p = bio.snowcrab::load.environment( year.assessment=year.assessment )  

SCD = project.datadirectory("bio.snowcrab")
media_loc = project.codedirectory("bio.snowcrab", "inst", "markdown", "media")

require(gt)  # table formatting

outtabledir = file.path( p$annual.results, "tables" )

years = as.character(1996: year.assessment)

regions = c("cfanorth", "cfasouth", "cfa4x")
nregions = length(regions)

FD = fishery_data()  # mass in tonnes
fda = FD$summary_annual

dt = as.data.frame( fda[ which(fda$yr %in% c(year.assessment - c(0:10))),] )
dt =  dt[,c("region", "yr", "Licenses", "TAC", "landings", "effort", "cpue")] 
names(dt) = c("Region", "Year", "Licenses", "TAC", "Landings", "Effort", "CPUE") 
rownames(dt) = NULL

odb0 = setDT(observer.db("odb"))
odb0$region = NA
for ( reg in regions) {
  r = polygon_inside(x = odb0, region = aegis.polygons::polygon_internal_code(reg), planar=FALSE)
  odb0$region[r] = reg
}

# bycatch summaries
o_cfaall = observer.db( DS="bycatch_summary", p=p,  yrs=p$yrs, region="cfaall" )
o_cfanorth = observer.db( DS="bycatch_summary", p=p,  yrs=p$yrs, region="cfanorth" )   
o_cfasouth = observer.db( DS="bycatch_summary", p=p,  yrs=p$yrs, region="cfasouth" )   
o_cfa4x = observer.db( DS="bycatch_summary", p=p,  yrs=p$yrs, region="cfa4x" )   

```  

 
## Overview of sampling locations

### Observer locations
 
```{r}
#| label: map-observer-locations
#| eval: true 
#| output: true
#| fig-cap: "Snow Crab At-sea-observer locations."
#| fig-dpi: 144
#| fig-height: 4

loc = file.path( SCD, "output", "maps", "observer.locations" )
yrsplot = p$year.assessment + c(0:-3) 
 
quietly( map.observer.locations( p=p, basedir=loc, years=yrsplot ) )
 
fns = paste( "observer.locations", yrsplot, "png", sep="." ) 
fn = file.path( loc, fns ) 

include_graphics( fn )

``` 


### Logbook recorded locations
   
```{r}
#| label: map-logbook-locations
#| eval: true 
#| output: true
#| fig-cap: "Snow Crab logbook locations."
#| fig-dpi: 144
#| fig-height: 4

loc = file.path( SCD, "output", "maps", "logbook.locations" )
yrsplot = year.assessment + c(0:-3)

quietly( map.logbook.locations( p=p, basedir=loc, years=yrsplot ) )

fns = paste( "logbook.locations", yrsplot, "png", sep="." ) 
fn = file.path( loc, fns ) 

include_graphics( fn ) 

``` 


## Fishery performance

### NENS

```{r}
#| label: table-fishery-nens-perf
#| tbl-cap: "Fishery performance statistics in NENS. Units are: TACs and Landings (tons, t), Effort ($\\times 10^3$ trap hauls, th) and CPUE (kg/th)."
#| eval: true
#| output: true

ii = which(dt$Region=="cfanorth")
oo = dt[ii, c("Year", "Licenses", "TAC", "Landings", "Effort", "CPUE")] 
gt::gt(oo) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
```

### SENS

```{r}
#| label: table-fishery-sens-perf
#| tbl-cap: "Fishery performance statistics in SENS. Units are: TACs and Landings (tons, t), Effort ($\\times 10^3$ trap hauls, th) and CPUE (kg/th)."
#| eval: true
#| output: true

ii = which(dt$Region=="cfasouth")
oo = dt[ii, c("Year", "Licenses", "TAC", "Landings", "Effort", "CPUE")] 
gt::gt(oo) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
```

### 4X

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



### Summary figures 


Create the maps:

```{r}
#| label: create-maps-fishery-performance
#| eval: true
#| output: true

fp_loc = file.path( p$project.outputdir, "maps", "logbook","snowcrab","annual" )

quietly( 
  map.fisheries.data( 
    outdir=fp_loc, 
    probs=c(0,0.975),
    plot_crs=st_crs( p$aegis_proj4string_planar_km ),
    outformat="png"
  )
)
```

Create aggregate time-series:

```{r}
#| label: create-timeseries-fishery-performance
#| eval: true
#| output: true

fpts_loc = file.path( p$annual.results,  "timeseries", "fishery")


quietly( 
figure.landings.timeseries( yearmax=p$year.assessment, outdir=fpts_loc, outfile="landings.ts"  )
)


quietly( 
figure.effort.timeseries( yearmax=p$year.assessment, outdir=fpts_loc, outfile="effort.ts"  )
)


quietly( 
figure.cpue.timeseries( yearmax=p$year.assessment, outdir=fpts_loc, outfile="cpue.ts"  )
)

```


#### Landings

```{r}
#| label: landings-timeseries
#| eval: true 
#| output: true
#| fig-cap: "Landings (t) of Snow Crab on the SSE. For 4X, the year refers to the starting year of the season. Inset is a closeup view of the timeseries for N-ENS and 4X."
#| fig-dpi: 144
#| fig-height: 8

include_graphics( file.path( SCD, "assessments", year.assessment, "timeseries", "fishery", "landings.ts.png" ) )
``` 


```{r}
#| label: landings-map
#| eval: true 
#| output: true
#| fig-cap: "Snow Crab landings from fisheries logbook data for previous and current years (tons per 10 km x 10 km grid)."
#| fig-dpi: 144
#| fig-height: 4

loc0 = file.path( SCD, "output", "maps", "logbook", "snowcrab", "annual", "landings" )
fn1 = file.path( loc0, paste( "landings", year_previous,   "png", sep=".") ) 
fn2 = file.path( loc0, paste( "landings", year.assessment, "png", sep=".") ) 
knitr::include_graphics( c(fn1, fn2 ) )

```



#### Effort
 

```{r}
#| label: effort-timeseries
#| eval: true 
#| output: true
#| fig-cap: "Temporal variations in fishing effort $\\times 10^3$ trap hauls."
#| fig-dpi: 144
#| fig-height: 8

fn1=file.path( SCD, "assessments", year.assessment, "timeseries", "fishery",   "effort.ts.png" )
knitr::include_graphics( fn1 ) 
```


```{r}
#| label: effort-map
#| eval: true 
#| output: true
#| fig-cap: "Snow Crab fishing effort from fisheries logbook data for previous and current years (no $\\times 10^3$ per 10 km X 10 km grid)."
#| fig-dpi: 144
#| fig-height: 8
 
loc0= file.path( SCD, "output", "maps", "logbook", "snowcrab", "annual", "effort" )
fn1 = file.path( loc0, paste( "effort", year_previous,   "png", sep=".") ) 
fn2 = file.path( loc0, paste( "effort", year.assessment, "png", sep=".") ) 
include_graphics(  c(fn1, fn2) )

```
 


#### Catch rates


```{r}
#| label: cpue-timeseries
#| eval: true 
#| output: true
#| fig-cap: "Temporal variations in crude catch rates of Snow Crab (kg/trap haul)."
#| fig-dpi: 144
#| fig-height: 4
 
include_graphics( file.path( SCD, "assessments", year.assessment, "timeseries", "fishery",   "cpue.ts.png" ) ) 
```

 


```{r}
#| label: cpue-map
#| eval: true 
#| output: true
#| fig-cap: "Snow Crab crude catch rates on the Scotian Shelf for previous and current years. Units are kg/trap haul per 10 km x 10 km grid."
#| fig-dpi: 144
#| fig-height: 4
 
loc0= file.path( SCD, "output", "maps", "logbook", "snowcrab", "annual", "cpue" )
fn1 = file.path( loc0, paste( "cpue", year_previous,   "png", sep=".") ) 
fn2 = file.path( loc0, paste( "cpue", year.assessment, "png", sep=".") ) 
knitr::include_graphics( c(fn1, fn2 ) )
# \@ref(fig:cpue-map)  

```



## Carapace condition from observed data

#### \< 95mm CW

```{r}
#| label: setup-observer-data
#| eval: true
#| output: true

odb = odb0[ cw < 95 & prodcd_id==0 & shell %in% c(1:5) & region %in% regions & sex==0, ]  # male
```

NENS:

```{r}
#| label: table-fishery-nens-sublegal
#| tbl-cap: "Fishery performance statistics in NENS. Distribution of at sea observations of males less than 95 mm CW by year and shell condition."
#| eval: true
#| output: true

resN = dcast( odb0[ region=="cfanorth", .(N=.N), by=.(fishyr, shell) ], fishyr  ~ shell, value.var="N", fill=0, drop=FALSE, na.rm=TRUE )
if ( "NA" %in% names(resN) ) resN$"NA" = NULL
names(resN) = c("Year", "CC1", "CC2", "CC3", "CC4", "CC5" )
resN$Total = rowSums( resN[, 2:6 ], na.rm=TRUE)
resN[, 2:6 ] = round(resN[, 2:6 ] / resN$Total * 100, digits=2)
gt::gt(resN) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
```

SENS:

```{r}
#| eval: true
#| output: true
#| label: table-fishery-sens-sublegal
#| tbl-cap: "Fishery performance statistics in SENS. Distribution of at sea observations of males less than 95 mm CW by year and shell condition."

resS = dcast( odb0[ region=="cfasouth", .(N=.N), by=.(fishyr, shell) ], fishyr  ~ shell, value.var="N", fill=0, drop=FALSE, na.rm=TRUE )
if ( "NA" %in%  names(resS)) resS$"NA" = NULL
names(resS) = c("Year", "CC1", "CC2", "CC3", "CC4", "CC5" )
resS$Total = rowSums( resS[, 2:6 ], na.rm=TRUE)
resS[, 2:6 ] = round(resS[, 2:6 ] / resS$Total * 100, digits=2)
gt::gt(resS) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
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
gt::gt(resX) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
```

####  \>= 95mm CW

```{r}
#| eval: true
#| output: true
#| label: data-observer-95plus

odb = odb0[ cw >= 95 & cw < 170  & prodcd_id==0 & shell %in% c(1:5) & region %in% regions & sex==0, ]  # male
```

NENS:

```{r}
#| eval: true
#| output: true
#| label: table-fishery-nens-comm
#| tbl-cap: "Fishery performance statistics in NENS. Distribution of at sea observations of males greater than 95 mm CW by year and shell condition."

resN = dcast( odb[ region=="cfanorth", .(N=.N), by=.(fishyr, shell) ], fishyr  ~ shell, value.var="N", fill=0, drop=FALSE, na.rm=TRUE )
names(resN) = c("Year", "CC1", "CC2", "CC3", "CC4", "CC5" )
resN$Total = rowSums( resN[, 2:6 ], na.rm=TRUE)
resN[, 2:6 ] = round(resN[, 2:6 ] / resN$Total * 100, digits=2)
gt::gt(resN) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
```

SENS:

```{r}
#| eval: true
#| output: true
#| label: table-fishery-sens-comm
#| tbl-cap: "Fishery performance statistics in SENS. Distribution of at sea observations of males greater than 95 mm CW by year and shell condition."

resS = dcast( odb[ region=="cfasouth", .(N=.N), by=.(fishyr, shell) ], fishyr  ~ shell, value.var="N", fill=0, drop=FALSE, na.rm=TRUE )
names(resS) = c("Year", "CC1", "CC2", "CC3", "CC4", "CC5" )
resS$Total = rowSums( resS[, 2:6 ], na.rm=TRUE)
resS[, 2:6 ] = round(resS[, 2:6 ] / resS$Total * 100, digits=2)
gt::gt(resS) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
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
gt::gt(resX) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
```

#### Soft-shell

There are two possible definitions:

-   durometer \< 68 (Soft, Total)
-   carapace conditions 1 and 2 (SoftSC, TotalSC)

```{r}
#| eval: true
#| output: true
#| label: data-observer-sub95

odb = odb0[ cw >= 95 & cw < 170  & prodcd_id==0 & shell %in% c(1:5) & region %in% regions & sex==0, ]  # male
shell_condition = odb[ !is.na(odb$region), .N, by=.(region, fishyr, shell) ]
shell_condition[, total:=sum(N, na.rm=TRUE), by=.(region, fishyr)]
shell_condition$percent = round(shell_condition$N / shell_condition$total, 3) * 100
shell_condition$Year = shell_condition$fishyr
```

NENS:

```{r}
#| eval: true
#| output: true
#| label: table-fishery-nens-soft-durometer
#| tbl-cap: "Fishery performance statistics in NENS. Distribution of at sea observations of males soft-shelled based on durometer (<68) and shell condition (1 and 2, SC)."

softN  = odb[ region=="cfanorth" & durometer <  68, .(Soft=.N), by=.(fishyr ) ] 
totalN = odb[ region=="cfanorth" & is.finite(durometer) , .(Total=.N), by=.(fishyr) ] 
resN = softN[totalN, on="fishyr"]
resN = resN[, .(Year=fishyr, Soft=round(Soft/Total*100,2), Total=Total) ]  
ssN = shell_condition[ region=="cfanorth" & shell %in% c(1,2), .(SoftSC=sum(percent), TotalSC=unique(total)[1]), by=.(Year)]
resN = resN[ssN, on="Year"]
gt::gt(resN) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
```

SENS:

```{r}
#| eval: true
#| output: true
#| label: table-fishery-sens-soft-durometer
#| tbl-cap: "Fishery performance statistics in SENS. Distribution of at sea observations of males soft-shelled based on durometer (<68) and shell condition (1 and 2, SC)."

softS  = odb[ region=="cfasouth" & durometer <  68, .(Soft=.N), by=.(fishyr ) ] 
totalS = odb[ region=="cfasouth" & is.finite(durometer) , .(Total=.N), by=.(fishyr) ] 
resS = softS[totalS, on="fishyr"]
resS = resS[, .(Year=fishyr, Soft=round(Soft/Total*100,2), Total=Total) ]  
ssS = shell_condition[ region=="cfasouth" & shell %in% c(1,2), .(SoftSC=sum(percent), TotalSC=unique(total)[1]), by=.(Year)]
resS = resS[ssS, on="Year"]
gt::gt(resS) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
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
gt::gt(resX) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
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


 
#### Size frequency distributions by carapace condition

```{r}
#| label: create-size-frequency-carapace-condition
#| eval: true
#| output: true

quietly( 
  figure.observed.size.freq( regions = c("cfanorth", "cfasouth", "cfa4x"), years="all", outdir=file.path( p$annual.results, "figures", "size.freq", "observer")  )
)

```

NENS

```{r}
#| label: size-frequency-carapace-condition-observer-nens
#| eval: true 
#| output: true
#| fig-cap: "Size frequency distribution of Snow Crab sampled by At-sea-observers, broken down by Carapace Condition (CC). NENS. Vertical lines indicate 95 mm Carapace Width, the minimum legal commercial size."
#| fig-dpi: 144
#| fig-height: 4
#| layout-ncol: 2
  
loc = file.path( SCD, "assessments", year.assessment, "figures", "size.freq", "observer")

fn1 = file.path( loc, paste( "size.freqcfanorth", (year_previous), ".png", sep="" ) )
fn2 = file.path( loc, paste( "size.freqcfanorth", (year.assessment  ), ".png", sep="" ) ) 
 
include_graphics(  c(fn1, fn2 ) ) 

```

SENS

```{r}
#| label: size-frequency-carapace-condition-observer-sens
#| eval: true 
#| output: true
#| fig-cap: "Size frequency distribution of Snow Crab sampled by At-sea-observers, broken down by Carapace Condition (CC). SENS. Vertical lines indicate 95 mm Carapace Width, the minimum legal commercial size."
#| fig-dpi: 144
#| fig-height: 4
#| layout-ncol: 2

loc = file.path( SCD, "assessments", year.assessment, "figures", "size.freq", "observer")
 
fn1 = file.path( loc, paste( "size.freqcfasouth", (year_previous), ".png", sep="" ) )
fn2 = file.path( loc, paste( "size.freqcfasouth", (year.assessment  ), ".png", sep="" ) ) 
 
include_graphics(  c(fn1, fn2 ) ) 

```

4X

```{r}
#| label: size-frequency-carapace-condition-observer-4x
#| eval: true 
#| output: true
#| fig-cap: "Size frequency distribution of Snow Crab sampled by At-sea-observers, broken down by Carapace Condition (CC). 4X. The year refers to the starting year of the season; the current season is ongoing. Vertical lines indicate 95 mm Carapace Width, the minimum legal commercial size."
#| fig-dpi: 144
#| fig-height: 4
#| layout-ncol: 2
  
loc = file.path( SCD, "assessments", year.assessment, "figures", "size.freq", "observer")

fn1 = file.path( loc, paste( "size.freqcfa4x", (year_previous), ".png", sep="" ) )
fn2 = file.path( loc, paste( "size.freqcfa4x", (year.assessment  ), ".png", sep="" ) )
 
include_graphics(  c(fn1, fn2 ) ) 

```

  
  
## Discard rates (by-catch in fishery)

NENS:discard

```{r}
#| eval: true
#| output: true
#| warning: false
#| error: false 
#| label: table-fishery-nens-discard
#| tbl-cap: "Fishery performance statistics in NENS. Average by-catch discard rate by weight observed (kg/trap haul; and standard deviation, SD)."

region="cfanorth"
o = observer.db( DS="bycatch_summary", p=p,  yrs=p$yrs, region=region )   
resN = o$eff_summ[ order(fishyr), ]
names(resN) = c("Year", "Discards", "SD")
resN$Discards = round( resN$Discards*100, 2)
resN$SD = round( resN$SD*100, 2)
gt::gt(resN) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
  summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
  footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
  row_group.padding = gt::px(1))
```

SENS:

```{r}
#| eval: true
#| output: true
#| warning: false
#| error: false 
#| label: table-fishery-sens-discard
#| tbl-cap: "Fishery performance statistics in SENS. Average by-catch discard rate by weight observed (kg/trap haul; and standard deviation, SD)."

region="cfasouth"
o = observer.db( DS="bycatch_summary", p=p,  yrs=p$yrs, region=region )   
resS = o$eff_summ[ order(fishyr), ]
names(resS) = c("Year", "Discards", "SD")
resS$Discards = round( resS$Discards*100, 2)
resS$SD = round( resS$SD*100, 2)
gt::gt(resS) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
  summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
  footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
  row_group.padding = gt::px(1))
```

4X:

```{r}
#| eval: true
#| output: true
#| warning: false
#| error: false 
#| label: table-fishery-4x-discard
#| tbl-cap: "Fishery performance statistics in 4X. Average by-catch discard rate by weight observed (kg/trap haul; and standard deviation, SD)."

region="cfa4x"
o = observer.db( DS="bycatch_summary", p=p,  yrs=p$yrs, region=region )   
resX = o$eff_summ[ order(fishyr), ]
names(resX) = c("Year", "Discards", "SD")
resX$Discards = round( resX$Discards*100, 2)
resX$SD = round( resX$SD*100, 2)
gt::gt(resX) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
  summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
  footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
  row_group.padding = gt::px(1))
```


Approach: estimate bycatch from at seas observed data and project onto marfis data
 

```{r}
#| eval: true
#| output: false 
#| label: bycatch-setup

bycatch_dir = file.path( p$annual.results, "bycatch")
years = as.character(1996: year.assessment)

```

We estimate bycatch using at-sea-observed effort and catch and re-scaling these naively to total snow crab fishery effort associated with the observations by year: 

### All Maritimes

```{r}
#| eval: true
#| output: false
#| warning: false
#| error: false 

region="cfaall"
o = observer.db( DS="bycatch_summary", p=p,  yrs=p$yrs, region=region )   
o$bycatch_table[ o$bycatch_table==0 ] = NA
o$bycatch_table[ is.na(o$bycatch_table) ] = "."
o$bycatch_table_catch[ o$bycatch_table_catch==0 ] = NA
o$bycatch_table_catch[ is.na(o$bycatch_table_catch) ] = "."

```

#### Discard rates of snow crab, by weight

```{r}
#| label: figure-discard_maritimes
#| fig-cap: "At sea observed rate of snow crab discards relative to total catch (discard + kept), all Maritimes."
#| fig-dpi: 144
#| fig-height: 4
 
pl = ggplot( o$eff_summ, aes(x=fishyr, y=discard_rate, ymin=discard_rate-discard_rate_sd, ymax=discard_rate+discard_rate_sd) ) +
  geom_pointrange()  + # Vertical line with point in the middle
  geom_errorbar(width = 0.1, col="brown") + # Standard error bars
  geom_point(size = 1.5, col="darkred") +
  labs(x="Year", y="Discard rate of snow crab (At sea observed, by weight)" )
(pl)
``` 

#### Bycatch catch rates

```{r}
#| label: bycatch_rates_all
#| tbl-cap: "At sea observed bycatch rates in Maritimes"
#| fig-dpi: 144
#| fig-height: 15

plot( o$spec ~ o$bct, xlab = "At sea observed catch rate in snow crab fishery (kg/trap)", ylab="Species", type="p", cex=1.1, pch=19, col="darkorange", xlim=c(0, max(o$bct, na.rm=TRUE)*1.4), yaxt="n" )  
text( o$bct, o$spec,  labels=o$species, pos=4, srt=0 , cex=0.8, col="darkslateblue")
text( max(o$bct, na.rm=TRUE)*0.88, 2.5, labels=paste( "Snow crab CPUE (At sea obs., mean): ", o$bct_sc, " kg/trap"), col="darkred", cex=0.9 )
```

#### Bycatch table from effort

```{r}
#| label: bycatch_all
#| tbl-cap: "At sea observed bycatch estimates (kg) based upon effort rescaling of Maritimes Region snow crab fishery. Dots indicated low values. Where species exist in a list but there is no data, this indicates some historical bycatch. The average is only for the years shown."

gt::gt(o$bycatch_table) |> 
  gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
```


#### Bycatch table from catch

```{r}
#| label: bycatch_all_catch
#| tbl-cap: "At sea observed bycatch estimates (kg) based upon catch rescaling of Maritimes Region snow crab fishery. Dots indicated low values. Where species exist in a list but there is no data, this indicates some historical bycatch. The average is only for the years shown."

gt::gt(o$bycatch_table_catch) |> 
  gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
```

#### Entanglements of large megafauna 

```{r}
#| warning: false
#| error: false 
#| tbl-cap: "Entanglements Maritimes"
#| label: bycatch_entanglements_all

oss = o$oss  # subset for region of interest

print("whale entaglements:")
whales = oss[ grep("whale", common, ignore.case=TRUE), ]
print(whales[, .N, by=.(yr)] )

print("leatherback entaglements:")
leatherback = oss[ grep("LEATHERBACK", common, ignore.case=TRUE), ]
print(leatherback[, .N, by=.(yr)])

print("basking sharks entaglements:")
basking_shark = oss[ grep("BASKING SHARK",  common, ignore.case=TRUE), ]
print(basking_shark[, .N, by=.(yr)])
```

Map of locations of entanglements:

```{r}
#| warning: false
#| error: false 
#| fig-cap: "Entanglement locations in Maritimes. Whales (red), Leatherback turtles (green), Basking shark (blue)."
#| label: bycatch_entanglements_map_all
#| fig-dpi: 144
#| fig-height: 5

plot(lat~-lon, oss, pch=".", col="lightgray", xlim=c(-65.2, -57), ylim=c(42.9,47) )
points(lat~-lon, whales, pch=19, cex=1.5, col="darkred" )
points(lat~-lon, leatherback, pch=18, cex=1.5, col="darkgreen" )
points(lat~-lon, basking_shark, pch=17, cex=1.5, col="slateblue" )
```


### NENS

```{r}
#| eval: true
#| output: false
#| warning: false
#| error: false 

region="cfanorth"
o = observer.db( DS="bycatch_summary", p=p,  yrs=p$yrs, region=region )   
o$bycatch_table[ o$bycatch_table==0 ] = NA
o$bycatch_table[ is.na(o$bycatch_table) ] = "."
o$bycatch_table_catch[ o$bycatch_table_catch==0 ] = NA
o$bycatch_table_catch[ is.na(o$bycatch_table_catch) ] = "."

```

#### Discard rates of snow crab, by weight

```{r}
#| label: figure-discard_nens
#| fig-cap: "At sea observed rate of snow crab discards relative to total catch (discard + kept), NENS."
#| fig-dpi: 144
#| fig-height: 4

pl = ggplot( o$eff_summ, aes(x=fishyr, y=discard_rate, ymin=discard_rate-discard_rate_sd, ymax=discard_rate+discard_rate_sd) ) +
  geom_pointrange()  + # Vertical line with point in the middle
  geom_errorbar(width = 0.1, col="brown") + # Standard error bars
  geom_point(size = 1.5, col="darkred") +
  labs(x="Year", y="Discard rate of snow crab (At sea observed, by weight)" )
(pl)
``` 

#### Bycatch catch rates

```{r}
#| label: bycatch_rates_nens
#| tbl-cap: "At sea observed bycatch rates in NENS"
#| fig-dpi: 144
#| fig-height: 6

plot( o$spec ~ o$bct, xlab = "At sea observed catch rate in snow crab fishery (kg/trap)", ylab="Species", type="p", cex=1.1, pch=19, col="darkorange", xlim=c(0, max(o$bct, na.rm=TRUE)*1.4), yaxt="n" )  
text( o$bct, o$spec,  labels=o$species, pos=4, srt=0 , cex=0.8, col="darkslateblue")
text( max(o$bct, na.rm=TRUE)*0.88, 2.5, labels=paste( "Snow crab CPUE (At sea obs., mean): ", o$bct_sc, " kg/trap"), col="darkred", cex=0.9 )
```

#### Bycatch table from effort

```{r}
#| label: bycatch_NENS
#| tbl-cap: "At sea observed bycatch estimates (kg) based upon effort rescaling of NENS snow crab fishery. Dots indicated low values. Where species exist in a list but there is no data, this indicates some historical bycatch. The average is only for the years shown."

gt::gt(o$bycatch_table) |> 
  gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
```

#### Bycatch table from catch

```{r}
#| label: bycatch_NENS_catch
#| tbl-cap: "At sea observed bycatch estimates (kg) based upon catch rescaling of NENS snow crab fishery. Dots indicated low values. Where species exist in a list but there is no data, this indicates some historical bycatch. The average is only for the years shown."

gt::gt(o$bycatch_table_catch) |> 
  gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
```

#### Entanglements of large megafauna 

```{r}
#| warning: false
#| error: false 
#| label: bycatch_entanglements_ens
#| tbl-cap: "Entanglements NENS"

oss = o$oss  # subset for region of interst
whales = oss[ grep("whale", common, ignore.case=TRUE), ]
print(whales[, .(yr, lon, lat, uid, est_discard_wt)] )
leatherback = oss[ grep("LEATHERBACK", common, ignore.case=TRUE), ]
print(leatherback[, .(yr, lon, lat, uid, est_discard_wt)])
basking_shark = oss[ grep("BASKING SHARK",  common, ignore.case=TRUE), ]
print(basking_shark[, .(yr, lon, lat, uid, est_discard_wt)])
```


### SENS

```{r}
#| eval: true
#| output: false
#| warning: false
#| error: false 

region="cfasouth"
o = observer.db( DS="bycatch_summary", p=p,  yrs=p$yrs, region=region )   
o$bycatch_table[ o$bycatch_table==0 ] = NA
o$bycatch_table[ is.na(o$bycatch_table) ] = "."
o$bycatch_table_catch[ o$bycatch_table_catch==0 ] = NA
o$bycatch_table_catch[ is.na(o$bycatch_table_catch) ] = "."

```

#### Discard rates of snow crab, by weight

```{r}
#| label: figure-discard_sens
#| fig-cap: "At sea observed rate of snow crab discards relative to total catch (discard + kept), SENS."
#| fig-dpi: 144
#| fig-height: 4

pl = ggplot( o$eff_summ, aes(x=fishyr, y=discard_rate, ymin=discard_rate-discard_rate_sd, ymax=discard_rate+discard_rate_sd) ) +
  geom_pointrange()  + # Vertical line with point in the middle
  geom_errorbar(width = 0.1, col="brown") + # Standard error bars
  geom_point(size = 1.5, col="darkred") +
  labs(x="Year", y="Discard rate of snow crab (At sea observed, by weight)" )
(pl)
``` 

#### Bycatch catch rates

```{r}
#| label: bycatch_rates_sens
#| tbl-cap: "At sea observed bycatch rates in SENS"
#| fig-dpi: 144
#| fig-height: 14

plot( o$spec ~ o$bct, xlab = "At sea observed catch rate in snow crab fishery (kg/trap)", ylab="Species", type="p", cex=1.1, pch=19, col="darkorange", xlim=c(0, max(o$bct, na.rm=TRUE)*1.4), yaxt="n" )  
text( o$bct, o$spec,  labels=o$species, pos=4, srt=0 , cex=0.8, col="darkslateblue")
text( max(o$bct, na.rm=TRUE)*0.88, 2.5, labels=paste( "Snow crab CPUE (At sea obs., mean): ", o$bct_sc, " kg/trap"), col="darkred", cex=0.8 )
```

#### Bycatch table from effort

```{r}
#| label: bycatch_SENS
#| tbl-cap: "At sea observed bycatch estimates (kg) based upon effort rescaling of SENS snow crab fishery. Dots indicated low values. Where species exist in a list but there is no data, this indicates some historical bycatch. The average is only for the years shown."

gt::gt( o$bycatch_table ) |> 
  gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
```

#### Bycatch table from catch

```{r}
#| label: bycatch_SENS_catch
#| tbl-cap: "At sea observed bycatch estimates (kg) based upon catch rescaling of SENS snow crab fishery. Dots indicated low values. Where species exist in a list but there is no data, this indicates some historical bycatch. The average is only for the years shown."

gt::gt( o$bycatch_table_catch ) |> 
gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
  summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
  footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
  row_group.padding = gt::px(1))
```

```{r}
#| warning: false
#| error: false 
#| label: entanglements_SENS
#| tbl-cap: "Entanglements SENS"

# get large megafauna:
oss = o$oss  # subset for region of interst
whales = oss[ grep("whale", common, ignore.case=TRUE), ]
print(whales[, .(yr, lon, lat, uid, est_discard_wt)] )
leatherback = oss[ grep("LEATHERBACK", common, ignore.case=TRUE), ]
print(leatherback[, .(yr, lon, lat, uid, est_discard_wt)])
basking_shark = oss[ grep("BASKING SHARK",  common, ignore.case=TRUE), ]
print(basking_shark[, .(yr, lon, lat, uid, est_discard_wt)])

```



### 4X

```{r}
#| eval: true
#| output: false
#| warning: false
#| error: false 

region="cfa4x"
o = observer.db( DS="bycatch_summary", p=p,  yrs=p$yrs, region=region )   
o$bycatch_table[ o$bycatch_table==0 ] = NA
o$bycatch_table[ is.na(o$bycatch_table) ] = "."
o$bycatch_table_catch[ o$bycatch_table_catch==0 ] = NA
o$bycatch_table_catch[ is.na(o$bycatch_table_catch) ] = "."

```

#### Discard rates of snow crab, by weight

```{r}
#| label: figure-discard_4x
#| fig-cap: "At sea observed rate of snow crab discards relative to total catch (discard + kept), 4X."
#| fig-dpi: 144
#| fig-height: 4

pl = ggplot( o$eff_summ, aes(x=fishyr, y=discard_rate, ymin=discard_rate-discard_rate_sd, ymax=discard_rate+discard_rate_sd) ) +
  geom_pointrange()  + # Vertical line with point in the middle
  geom_errorbar(width = 0.1, col="brown") + # Standard error bars
  geom_point(size = 1.5, col="darkred") +
  labs(x="Year", y="Discard rate of snow crab (At sea observed, by weight)" )
(pl)
``` 


#### Bycatch catch rates

```{r}
#| label: bycatch_rates_4x
#| tbl-cap: "At sea observed bycatch rates in 4X"
#| fig-dpi: 144
#| fig-height: 6

plot( o$spec ~ o$bct, xlab = "At sea observed catch rate in snow crab fishery (kg/trap)", ylab="Species", type="p", cex=1.1, pch=19, col="darkorange", xlim=c(0, max(o$bct, na.rm=TRUE)*1.4), yaxt="n" )  
text( o$bct, o$spec,  labels=o$species, pos=4, srt=0 , cex=0.8, col="darkslateblue")
text( max(o$bct, na.rm=TRUE)*0.88, 2.5, labels=paste( "Snow crab CPUE (At sea obs., mean): ", o$bct_sc, " kg/trap"), col="darkred", cex=0.9 )
```

#### Bycatch table from effort

```{r}
#| label: bycatch_4X
#| tbl-cap: "At sea observed bycatch estimates (kg) based upon effort rescaling of 4X snow crab fishery. Dots indicated low values. Where species exist in a list but there is no data, this indicates some historical bycatch. The average is only for the years shown."

gt::gt(o$bycatch_table) |> 
gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
  summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
  footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
  row_group.padding = gt::px(1))
```


#### Bycatch table from catch

```{r}
#| label: bycatch_4X_catch
#| tbl-cap: "At sea observed bycatch estimates (kg) based upon catch rescaling of 4X snow crab fishery. Dots indicated low values. Where species exist in a list but there is no data, this indicates some historical bycatch. The average is only for the years shown."

gt::gt(o$bycatch_table_catch) |> 
gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
  summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
  footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
  row_group.padding = gt::px(1))
```

#### Entanglements of large megafauna 

```{r}
#| warning: false
#| error: false 
#| label: bycatch_entanglements_e4x
#| tbl-cap: "Entanglements 4X"

oss = o$oss  # subset for region of interst
whales = oss[ grep("whale", common, ignore.case=TRUE), ]
print(whales[, .(yr, lon, lat, uid, est_discard_wt)] )
leatherback = oss[ grep("LEATHERBACK", common, ignore.case=TRUE), ]
print(leatherback[, .(yr, lon, lat, uid, est_discard_wt)])
basking_shark = oss[ grep("BASKING SHARK",  common, ignore.case=TRUE), ]
print(basking_shark[, .(yr, lon, lat, uid, est_discard_wt)])
```



## To do next:

Estimate via carstm: requires a Poisson model of each species of catch (number) with offset of landings and covariates ... soon, ever?
 
