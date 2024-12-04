---
title: "Snow Crab fishery -- figures and tables"
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
  - snow crab trawl survey results
  - basic tables and figures
abstract: |
  Details of 2024 snow crab fishery performance. 
toc: true
toc-depth: 4
number-sections: true
highlight-style: pygments
# bibliography: media/references.bib  
# csl: media/canadian-journal-of-fisheries-and-aquatic-sciences.csl  # see https://www.zotero.org/styles for more
# license: "CC BY"
copyright: 
  holder: Jae S. Choi
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
  sens: 1
---

 

# Snow Crab fishery -- figures and tables

<!--

## Preamble

This is a markdown document. It can be viewed in formatted form via:
  
  - a web browser open the webpage: [02_fishery_summary.md](https://github.com/jae0/bio.snowcrab/tree/master/inst/markdown/02_fishery_summary.md)   

  - a web browser open the local file directly: [02_fishery_summary.md](../markdown/02_fishery_summary.md) (you might need to install a browser add-in), or 
  
  - a text editor with markdown awareness (e.g. Rstudio, VSCode, etc. ).


As this document uses the Quarto dialect of Markdown, you can easily create a report, by running one of the following in a shell command prompt or using Rstudio: 

```shell
# {via Quarto}

cd ~/bio/bio.snowcrab/inst/markdown

# sens as one group
make quarto FN=02_fishery_summary YR=2024 SOURCE=~/bio/bio.snowcrab/inst/markdown WK=~/bio.data/bio.snowcrab/assessments  DOCEXTENSION=html  PARAMS="sens:1"
 
# split sens into 23 and 24
make quarto FN=02_fishery_summary YR=2024 SOURCE=~/bio/bio.snowcrab/inst/markdown WK=~/bio.data/bio.snowcrab/assessments  DOCEXTENSION=html PARAMS="sens:2"


 
```

Or, see the Makefile and alter defaults to your needs. As Quarto does not pass params easily. So you must adjust "params" in yaml at the top of this fle, or use to quarto command such as: 

```shell
quarto ... -P year.assessment:$(YR) -P media_loc:$(MEDIA) 
```

[See here for more YAML options.](https://quarto.org/docs/output-formats/all-formats.html)
 


## Set up environment

  - Ensure the **year.assessment** is correct

  - Ensure that the data pulls of survey data stored on the ISSDB is complete and assimilated to the end of [01_snowcrab_data.md](https://github.com/jae0/bio.snowcrab/tree/master/inst/markdown/01_snowcrab_data.md).
  
-->



```{r}
#| label: setup
#| eval: true 
#| output: false
#| echo: false

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

# loadfunctions("bio.snowcrab")
source("~/bio/bio.snowcrab/R/observer.db.r")

SCD = project.datadirectory("bio.snowcrab")
media_loc = project.codedirectory("bio.snowcrab", "inst", "markdown", "media")

require(gt)  # table formatting

outtabledir = file.path( p$annual.results, "tables" )

years = as.character(1996: year.assessment)

lregions = list(region=c("cfanorth", "cfasouth", "cfa4x"))
reg_labels = c("N-ENS", "S-ENS", "CFA 4X")  # formatted for label

if (params$sens==2) {
  lregions = list(subarea=c("cfanorth", "cfa23",  "cfa24", "cfa4x"))
  reg_labels = c("CFA 20-22", "CFA 23", "CFA 24", "CFA 4X")  # formatted for label

}

vnr = names(lregions)
regions = unname( unlist(lregions) )
nregions = length(regions)

FD = fishery_data( regions=lregions)  # mass in tonnes

fda = FD$summary_annual
fdm = FD$summary_monthly
fdb = FD$summary_biweekly

dt = as.data.frame( fda[ which(fda$yr %in% c(year.assessment - c(0:10))),] )
dt =  dt[,c(vnr, "yr", "Licenses", "TAC", "landings", "effort", "cpue")] 
names(dt) = c("Region", "Year", "Licenses", "TAC", "Landings", "Effort", "CPUE") 
rownames(dt) = NULL

odb0 = setDT(observer.db("odb"))
odb0[[vnr]] = NA
for ( reg in regions) {
  r = polygon_inside(x = odb0, region = aegis.polygons::polygon_internal_code(reg), planar=FALSE)
  odb0[[vnr]][r] = reg
}

# bycatch summaries
BC = list()
for ( reg in regions) {

  BC[[reg]] = observer.db( DS="bycatch_summary", p=p,  yrs=p$yrs, region=reg  )  # using polygon test

  BC[[reg]]$bycatch_table_effort[ BC[[reg]]$bycatch_table_effort==0 ] = NA
  BC[[reg]]$bycatch_table_effort[ is.na(BC[[reg]]$bycatch_table_effort) ] = "."

  BC[[reg]]$bycatch_table_catch[ BC[[reg]]$bycatch_table_catch==0 ] = NA
  BC[[reg]]$bycatch_table_catch[ is.na(BC[[reg]]$bycatch_table_catch) ] = "."
 
}


 
```  

 
 

# Fishing locations recorded in logbooks
   
```{r}
#| label: map-logbook-locations
#| eval: true 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| echo: false 
#| layout-ncol: 2

# fig-cap: "Snow Crab logbook locations."

loc = file.path( SCD, "output", "maps", "logbook.locations" )
yrsplot = year.assessment + c(0:-3)
fns = paste( "logbook.locations", yrsplot, "png", sep="." ) 
fn = file.path( loc, fns ) 

include_graphics( fn ) 

``` 

Snow Crab logbook locations.

$~$


## Fishery performance

```{r}
#| echo: false
#| results: asis
#| tbl-cap: "Fishery performance statistics. Note: 4X years represent the starting year."
#| eval: true
#| output: true

for (r in 1:nregions) {
  reg = regions[r]
  REG = reg_labels[r]
  cat("#### ", REG, "\n")
  oo = dt[ which(dt$Region==reg), c("Year", "Licenses", "TAC", "Landings", "Effort", "CPUE")] 

  names(oo) = c( "Year", "Licenses", "TAC (kt)", "Landings (kt)", "Effort (1000 th)", "CPUE (kg/th)" )

  out = gt::gt(oo) |> gt::tab_options(table.font.size = 14, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
  print(out)
  cat("\n\n")
} 

``` 

Fishery performance statistics. Note: 4X years represent the starting year.

$~$


### Summary figures 


#### Landings

```{r}
#| label: landings-timeseries
#| eval: true 
#| output: true
#| fig-cap: "Landings (t) of Snow Crab on the SSE. For 4X, the year refers to the starting year of the season. Inset is a closeup view of the timeseries for N-ENS and 4X."
#| fig-dpi: 144
#| fig-height: 8

if (params$sens==1) {
  ts_dir = file.path( SCD, "assessments", year.assessment, "timeseries", "fishery" )
} else if (params$sens==2) {
  ts_dir = file.path( SCD, "assessments", year.assessment, "timeseries", "fishery", "split" )
}

include_graphics( file.path( ts_dir, "landings.ts.png" ) )
``` 

$~$


```{r}
#| label: landings-map
#| eval: true 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| echo: false 
#| layout-ncol: 2

# fig-cap: "Snow Crab landings from fisheries logbook data for previous and current years (tons per 10 km x 10 km grid)."

loc0 = file.path( SCD, "output", "maps", "logbook", "snowcrab", "annual", "landings" )
yrsplot = year.assessment + c(0:-3)
fn = file.path( loc0, paste( "landings", yrsplot, "png", sep=".") ) 

include_graphics( fn ) 

```

Snow Crab landings from fisheries logbook data for previous and current years (tons per 10 km x 10 km grid).

$~$


#### Effort
 

```{r}
#| label: effort-timeseries
#| eval: true 
#| output: true
#| fig-cap: "Temporal variations in fishing effort."
#| fig-dpi: 144
#| fig-height: 8



if (params$sens==1) {
  ts_dir = file.path( SCD, "assessments", year.assessment, "timeseries", "fishery" )
} else if (params$sens==2) {
  ts_dir = file.path( SCD, "assessments", year.assessment, "timeseries", "fishery", "split" )
}

include_graphics( file.path( ts_dir, "effort.ts.png" ) )
 

```


$~$


```{r}
#| label: effort-map
#| eval: true 
#| output: true
#| fig-dpi: 144
#| fig-height: 8
#| echo: false 
#| layout-ncol: 2

# fig-cap: "Snow Crab fishing effort from fisheries logbook data for previous and current years (no $\times 10^3$ per 10 km X 10 km grid)."
 
loc0 = file.path( SCD, "output", "maps", "logbook", "snowcrab", "annual", "effort" )
yrsplot = year.assessment + c(0:-3)
fn = file.path( loc0, paste( "effort", yrsplot, "png", sep=".") ) 

include_graphics( fn ) 

```
 
Snow Crab fishing effort from fisheries logbook data for previous and current years (x 10$~^3$ per 10 km x 10 km grid).

$~$


#### Catch rates


```{r}
#| label: cpue-timeseries
#| eval: true 
#| output: true
#| fig-cap: "Temporal variations in crude catch rates of Snow Crab (kg/trap haul)."
#| fig-dpi: 144
#| fig-height: 4
 
if (params$sens==1) {
  ts_dir = file.path( SCD, "assessments", year.assessment, "timeseries", "fishery" )
} else if (params$sens==2) {
  ts_dir = file.path( SCD, "assessments", year.assessment, "timeseries", "fishery", "split" )
}

include_graphics( file.path( ts_dir, "cpue.ts.png" ) ) 
```

 
$~$



```{r}
#| label: cpue-map
#| eval: true 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| echo: false 
#| layout-ncol: 2

# fig-cap: "Snow Crab crude catch rates on the Scotian Shelf for previous and current years. Units are kg/trap haul per 10 km x 10 km grid."

loc0= file.path( SCD, "output", "maps", "logbook", "snowcrab", "annual", "cpue" )
yrsplot = year.assessment + c(0:-3)
fn = file.path( loc0, paste( "cpue", yrsplot, "png", sep=".") ) 

include_graphics( fn ) 
 
```

Snow Crab crude catch rates on the Scotian Shelf for previous and current years. Units are kg/trap haul per 10 km x 10 km grid.

$~$


# At-sea Observed fishing locations
 
```{r}
#| label: map-observer-locations
#| eval: true 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| echo: false 
#| layout-ncol: 2

# fig-cap: "Snow Crab At-sea-observer locations."

loc = file.path( SCD, "output", "maps", "observer.locations" )
yrsplot = p$year.assessment + c(0:-3) 
 
fns = paste( "observer.locations", yrsplot, "png", sep="." ) 
fn = file.path( loc, fns ) 

include_graphics( fn )

``` 

Snow Crab At-sea-observer locations.

$~$


## At-sea-observed effort

```{r}
#| label: table-observer-number-trips
#| tbl-cap: "Number of at-sea observed events."
#| eval: true
#| output: true

oo = dcast( odb0[ fishyr>=2004,.(N=length(unique(tripset))), by=c(vnr, "fishyr")], 
  fishyr ~ get(vnr), value.var="N", fill=0, drop=FALSE, na.rm=TRUE )
if ( "NA" %in% names(oo) ) oo$"NA" = NULL
keep = c("fishyr", regions)
oo = oo[,..keep]
names(oo) = c("Year", reg_labels )
oo$Total = rowSums( oo[, 2:nregions ], na.rm=TRUE)

gt::gt(oo) |> gt::tab_options(table.font.size = 14, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
```


$~$


```{r}
#| label: table-observer-number-crab
#| tbl-cap: "Number of at-sea observed crab."
#| eval: true
#| output: true

oo = dcast( odb0[ fishyr>=2004,.(N=.N), by=c(vnr, "fishyr")], 
  fishyr ~ get(vnr), value.var="N", fill=0, drop=FALSE, na.rm=TRUE )
if ( "NA" %in% names(oo) ) oo$"NA" = NULL

keep = c("fishyr", regions)
oo = oo[,..keep]
names(oo) = c("Year", reg_labels )

oo[, 2:nregions] = round( oo[, 2:nregions] )
oo$Total = rowSums( oo[, 2:nregions ], na.rm=TRUE)

gt::gt(oo) |> gt::tab_options(table.font.size = 14, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
```

$~$


```{r}
#| label: table-observer-weight-crab
#| tbl-cap: "Total weight of at-sea observed crab (kg)."
#| eval: true
#| output: true

oo = dcast( odb0[ fishyr>=2004,.(N=sum( mass, na.rm=TRUE)/1000 ), by=c(vnr, "fishyr")], 
  fishyr ~ get(vnr), value.var="N", fill=0, drop=FALSE, na.rm=TRUE )
if ( "NA" %in% names(oo) ) oo$"NA" = NULL

keep = c("fishyr", regions)
oo = oo[,..keep]
names(oo) = c("Year", reg_labels )

oo[, 2:nregions] = round( oo[, 2:nregions] )
oo$Total = rowSums( oo[, 2:nregions ], na.rm=TRUE)

gt::gt(oo) |> gt::tab_options(table.font.size = 14, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
```



## Carapace condition from observed data

#### Carapace condition from observed data: \< 95mm CW

```{r}
#| echo: false
#| results: asis
#| label: table-observer-sublegal  
#| tbl-cap: "Fishery performance statistics: Distribution of at sea observations of males less than 95 mm CW by year and shell condition."
#| eval: true
#| output: true
#| layout-ncol: 2

odb = odb0[ cw < 95 & prodcd_id==0 & shell %in% c(1:5) & get(vnr) %in% regions & sex==0, ]  # male
  
for (r in 1:nregions){
  reg = regions[r]
  REG = reg_labels[r] 
  cat("#### ", REG, "\n")
  oo = dcast( odb0[ get(vnr)==reg, .(N=.N), by=.(fishyr, shell) ], fishyr  ~ shell, value.var="N", fill=0, drop=FALSE, na.rm=TRUE )
  if ( "NA" %in% names(oo) ) oo$"NA" = NULL
  
  keep = c("fishyr",  "1", "2", "3", "4", "5" )
  oo = oo[,..keep]
  names(oo) = c("Year", "CC1", "CC2", "CC3", "CC4", "CC5" )

  oo$Total = rowSums( oo[, 2:6 ], na.rm=TRUE)
  oo[, 2:6 ] = round(oo[, 2:6 ] / oo$Total * 100, digits=1)
  out = gt::gt(oo) |> gt::tab_options(table.font.size = 14, data_row.padding = gt::px(1), 
      summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
      footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
      row_group.padding = gt::px(1))
  print(out)
  cat("\n\n")
} 

``` 


Fishery performance statistics: Distribution of at sea observations of males less than 95 mm CW by year and shell condition.

$~$


#### Carapace condition from observed data: \>= 95mm CW

```{r}
#| echo: false
#| results: asis
#| label: table-observer-legal  
#| tbl-cap: "Fishery performance statistics: Distribution of at sea observations of males greater than 95 mm CW by year and shell condition."
#| eval: true
#| output: true
#| layout-ncol: 2

odb = odb0[ cw >= 95 & cw < 170  & prodcd_id==0 & shell %in% c(1:5) & get(vnr) %in% regions & sex==0, ]  # male
  
for (r in 1:nregions){
  reg = regions[r]
  REG = reg_labels[r]
  cat("#### ", REG, "\n")
  oo = dcast( odb0[ get(vnr)==reg, .(N=.N), by=.(fishyr, shell) ], fishyr  ~ shell, value.var="N", fill=0, drop=FALSE, na.rm=TRUE )
  if ( "NA" %in% names(oo) ) oo$"NA" = NULL
  keep = c("fishyr",  "1", "2", "3", "4", "5" )
  oo = oo[,..keep]
  names(oo) = c("Year", "CC1", "CC2", "CC3", "CC4", "CC5" )
  oo$Total = rowSums( oo[, 2:6 ], na.rm=TRUE)
  oo[, 2:6 ] = round(oo[, 2:6 ] / oo$Total * 100, digits=1)
  out = gt::gt(oo) |> gt::tab_options(table.font.size = 14, data_row.padding = gt::px(1), 
      summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
      footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
      row_group.padding = gt::px(1))
  print(out)
  cat("\n\n")
} 

``` 

Fishery performance statistics: Distribution of at sea observations of males greater than 95 mm CW by year and shell condition.

$~$


#### Soft-shell

There are two possible definitions:

-   durometer \< 68 (D)
-   carapace conditions 1 and 2 (CC)

```{r}
#| eval: true
#| output: true
#| echo: false
#| results: asis
#| label: data-observer-sub95
#| tbl-cap: "Fishery performance statistics: Distribution of at sea observations of males greater than 95 mm CW by year and shell condition." 
#| layout-ncol: 2

odb = odb0[ cw >= 95 & cw < 170  & prodcd_id==0 & shell %in% c(1:5) & get(vnr) %in% regions & sex==0, ]  # male
shell_condition = odb[ !is.na(get(vnr)), .N, by=c(vnr, "fishyr", "shell") ]
shell_condition[, total:=sum(N, na.rm=TRUE), by=c(vnr, "fishyr")]
shell_condition$percent = round(shell_condition$N / shell_condition$total, 3) * 100
shell_condition$Year = shell_condition$fishyr
 
for (r in 1:nregions){
  reg = regions[r]
  REG = reg_labels[r]
  cat("#### ", REG, "\n")

  soft  = odb[ get(vnr)==reg & durometer <  68, .(Soft=.N), by=.(fishyr ) ] 
  total = odb[ get(vnr)==reg & is.finite(durometer) , .(Total=.N), by=.(fishyr) ] 
  oo = soft[total, on="fishyr"]
  oo = oo[, .(Year=fishyr, Soft=round(Soft/Total*100, 1), Total=Total) ]  
  scond = shell_condition[ get(vnr)==reg & shell %in% c(1,2), .(SoftSC=sum(percent), TotalSC=unique(total)[1]), by=.(Year)]
  oo = oo[scond, on="Year"]

  names(oo) = c( "Year", "Soft (D)", "Total (D)", "Soft (CC)", "Total (CC)" )
  out = gt::gt(oo) |> gt::tab_options(table.font.size = 14, data_row.padding = gt::px(1), 
      summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
      footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
      row_group.padding = gt::px(1))
  print(out)
  cat("\n\n")
} 

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
#| echo: false
#| results: asis
#| eval: true
#| output: true 
#| label: size-frequency-carapace-condition-observer
#| fig-dpi: 144
#| fig-height: 6
#| layout-ncol: 2
#! fig-cap: "Size frequency distribution of Snow Crab sampled by At-sea-observers, broken down by Carapace Condition (CC). Vertical lines indicate 95 mm Carapace Width, the minimum legal commercial size."

odir = file.path( SCD, "assessments", year.assessment, "figures", "size.freq", "observer" )
years = p$year.assessment + c(0:-3) 
for (r in 1:nregions) {
  reg = regions[r]
  REG = reg_labels[r]
  cat("#### ", REG, "\n")
  fns = paste( "size.freq", reg,  years, "png", sep="." )  
  fn = check_file_exists(file.path( odir, fns ) )
  for (ff in fn) show_image(ff) 
  cat("\n\n")
}

``` 
  
$~$

  
## Total discards in fishery
 
```{r}
#| eval: true
#| output: true
#| echo: false
#| results: asis
#| label: table-fishery-discard-total 
#| tbl-cap: "Average by-catch discard rate by weight observed (kg/trap haul; and standard deviation, SD)." 
#| layout-ncol: 2
  

for (r in 1:nregions){
  reg = regions[r]
  REG = reg_labels[r]
  cat("#### ", REG, "\n")
 
  o = BC[[reg]]   
  oo = o$eff_summ[ order(fishyr), ]
  names(oo) = c("Year", "Discards", "SD")
  oo$Discards = round( oo$Discards*100, 1)
  oo$SD = round( oo$SD*100, 1)
  
  names(oo) = c("Year", "Discard mean", "Discard SD")
  out = gt::gt(oo) |> gt::tab_options(table.font.size = 14, data_row.padding = gt::px(1), 
      summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
      footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
      row_group.padding = gt::px(1))
  print(out)
  cat("\n\n")
} 

```
 



#### Discard rates of snow crab, by weight

```{r}
#| label: figure-discard_maritimes
#| fig-cap: "At sea observed rate of snow crab discards relative to total catch (discard + kept), all Maritimes."
#| fig-dpi: 144
#| fig-height: 6
 
oo = NULL
for (reg in regions) {
  o = BC[[reg]][["eff_summ"]]
  o$region = reg
  oo = rbind(oo, o)
}
oo$fishyr = jitter(oo$fishyr, amount=0.02)

color_map = c("#E69F00", "#56B4E9",  "#CC79A7" , "#D55E00", "#F0E442")[1:length(regions)]
shapes = c(15, 17, 19, 21, 23)[1:length(regions)]

pl = ggplot( oo, aes(x=fishyr, y=discard_rate, ymin=discard_rate-discard_rate_sd, ymax=discard_rate+discard_rate_sd, fill=region, colour=region) ) +
  geom_line( alpha=0.9, linewidth=1 ) +
  geom_point(aes(shape=region), size=5, alpha=0.7 )+
  geom_pointrange()  + # Vertical line with point in the middle
  geom_errorbar(width = 0.1, col="brown") + # Standard error bars
  geom_point(size = 1.5, col="darkred") +
  scale_colour_manual(values=color_map) +
  scale_fill_manual(values=color_map) +
  scale_shape_manual(values = shapes) +
  theme_light( base_size = 18) + 
  theme( legend.position="inside", legend.position.inside=c(0.2, 0.8), legend.title=element_blank()) 
  labs(x="Year", y="Discard rate of snow crab (At sea observed, by weight)" )
(pl)

``` 





```{r}
#| eval: true
#| output: false 
#| label: bycatch-setup

bycatch_dir = file.path( p$annual.results, "bycatch")
years = as.character(1996: year.assessment)

```

General approach: estimate bycatch from at seas observed data and project onto marfis data
 
### Estimates based on fisheries **effort**

Bycatch comes from at-sea-observed effort and catch. Rescale these naively to total snow crab fishery effort. 

```{r}
#| eval: true
#| output: true
#| echo: false
#| results: asis
#| label: table-fishery-discard-effort
#| tbl-cap: "Average by-catch discards (kg) using effort. Dots indicated low values. Where species exist in a list but there is no data, this indicates some historical bycatch. The average is only for the years shown." 
#| layout-ncol: 1
 
for (r in 1:nregions){
  reg = regions[r]
  REG = reg_labels[r]
  cat("#### ", REG, "\n")
  o = BC[[reg]]   
  oo = o$bycatch_table_effort
  oo[ oo==0 ] = NA
  oo[ is.na(oo) ] = "."
  out = gt::gt(oo) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
  print(out)
  cat("\n\n")
}


```

Average by-catch discards (kg) using effort. Dots indicated low values. Where species exist in a list but there is no data, this indicates some historical bycatch. The average is only for the years shown.


### Estimates based on fisheries **catch**

Bycatch comes from at-sea-observed effort and catch. Rescale these naively to total snow crab fishery catch. 

```{r}
#| eval: true
#| output: true
#| echo: false
#| results: asis
#| label: table-fishery-discard-catch 
#| tbl-cap: "Average by-catch discards (kg) using catch. Dots indicated low values. Where species exist in a list but there is no data, this indicates some historical bycatch. The average is only for the years shown." 
#| layout-ncol: 1
 
for (r in 1:nregions){
  reg = regions[r]
  REG = reg_labels[r]
  cat("#### ", REG, "\n")
 
  o = BC[[reg]]   
  oo = o$bycatch_table_catch 
  oo[ oo==0 ] = NA
  oo[ is.na(oo) ] = "."
  
  out = gt::gt(oo) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
  print(out)
  cat("\n\n")

}


```

Average by-catch discards (kg) using catch. Dots indicated low values. Where species exist in a list but there is no data, this indicates some historical bycatch. The average is only for the years shown.



```{r}
#| eval: false
#| label: bycatch_rates_all
#| tbl-cap: "At sea observed bycatch rates in Maritimes"
#| fig-dpi: 144
#| fig-height: 15

reg = "cfanorth"
o = BC[[reg]]

plot( o$spec ~ o$bct, xlab = "At sea observed catch rate in snow crab fishery (kg/trap)", ylab="Species", type="p", cex=1.1, pch=19, col="darkorange", xlim=c(0, max(o$bct, na.rm=TRUE)*1.4), yaxt="n" )  
text( o$bct, o$spec,  labels=o$species, pos=4, srt=0 , cex=0.8, col="darkslateblue")
text( max(o$bct, na.rm=TRUE)*0.88, 2.5, labels=paste( "Snow crab CPUE (At sea obs., mean): ", o$bct_sc, " kg/trap"), col="darkred", cex=0.9 )
```
 
 

#### Entanglements of large megafauna 

```{r}
#| warning: false
#| error: false 
#| tbl-cap: "Entanglements Maritimes"
#| label: bycatch_entanglements_all

reg = "cfanorth"
o = BC[[reg]]

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


<!--

## To do next:

Estimate via carstm: requires a Poisson model of each species of catch (number) with offset of landings and covariates ... soon, ever?
 
-->