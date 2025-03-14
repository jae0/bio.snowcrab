---
title: "Snow crab fishery summary" 
keywords: 
  - snow crab trawl survey results
  - basic tables and figures
abstract: |
  Details of 2024 snow crab fishery performance. 

metadata-files:
  - _metadata.yml

params:
  year_assessment: 2024
  year_start: 1999
  data_loc:  "~/bio.data/bio.snowcrab"
  sens: 1
  debugging: FALSE
  model_variation: logistic_discrete_historical
  todo: [fishery_results]
---

<!--

# Summary 1 of 4 -- This file is designed to be an HTML document that describes and summarizes the fishery performance. 


# sens as one group
make quarto FN=02_fishery_summary.md YR=2024 DATADIR=~/bio.data/bio.snowcrab DOCTYPE=html PARAMS="-P year_assessment:2024 -P sens:1 -P todo:[fishery_results,redo_data]"  --directory=~/bio/bio.snowcrab/inst/markdown 
 
# split sens into 23 and 24 (default behaviour)
make quarto FN=02_fishery_summary.md YR=2024 DATADIR=~/bio.data/bio.snowcrab DOCTYPE=html PARAMS="-P year_assessment:2024 -P sens:2 -P todo:[fishery_results,redo_data]" --directory=~/bio/bio.snowcrab/inst/markdown 
 
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

 

## Fishing locations

#### Fishing locations recorded in logbooks

```{r}
#| label: fig-map-logbook-locations
#| eval: true 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| echo: false 
#| layout-ncol: 2
#| fig-cap: "Snow Crab logbook locations."
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

loc = file.path( data_loc, "output", "maps", "logbook.locations" )
yrsplot = year_assessment + c(0:-3)
fns = paste( "logbook.locations", yrsplot, "png", sep="." ) 
fn = file.path( loc, fns ) 

include_graphics( fn ) 

```
 
$~$

#### Fishing locations recorded in logbooks: a closer look

```{r}
#| label: fig-map-logbook-locations-two-years
#| eval: true 
#| output: true
#| fig-dpi: 144
#| fig-height: 5
#| echo: false 
#| layout-ncol: 1
#| fig-cap: "Snow Crab logbook locations in the past two years."
#| fig-subcap: 
#|   - "N-ENS"
#|   - "S-ENS"
#|   - "4X"

fns = file.path( media_supplementary,  c(     
  "nens_past_two_years_fishing_positions.png",
  "sens_past_two_years_fishing_positions.png",
  "cfa4x_past_two_years_fishing_positions.png"
))
 
include_graphics(  fns  ) 

```
 
 
## Fishery performance

### Fishery performance: Summary

```{r}
#| label: tbl-fishery-performance
#| echo: false
#| results: asis
#| eval: true
#| output: true
#| tbl-cap: "Fishery performance statistics. Note: 4X years represent the starting year."
#| tbl-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

for (r in 1:nregions) {
  reg = regions[r]
  REG = reg_labels[r]
  cat(REG, "\n")
  oo = dt[ which(dt$Region==reg), c("Year", "Licenses", "TAC", "Landings", "Effort", "CPUE")] 

  names(oo) = c( "Year", "Licenses", "TAC (t)", "Landings (t)", "Effort (1000 th)", "CPUE (kg/th)" )

  out = gt::gt(oo) |> gt::tab_options(table.font.size = 14, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
  print(out)
  cat("\n\n")
} 

```
 
$~$

### Landings

#### Landings: annual

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

$~$

#### Landings: close-up

```{r}
#| label: fig-landings-decomposed
#| eval: true 
#| output: true
#| fig-dpi: 144
#| fig-height: 5
#| echo: false 
#| layout-ncol: 1
#| fig-cap: "Landings decomposed"
#| fig-subcap: 
#|   - "Fraction of landings in Spring"
#|   - "Fraction of landings in Summer"
#|   - "Fraction of landings in Winter"
#|   - "Cummulative landings"

fns = c( 
  "percent_spring_landings.png",
  "percent_summer_landings.png",
  "percent_winter_landings.png",
  "weekly_landing.png"
)
 
include_graphics( file.path( media_supplementary, fns )  ) 

```
 
$~$

#### Landings: spatial

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
 
$~$

### Effort

#### Effort: annual

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

$~$

#### Effort: spatial

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
 
$~$

### Catch rates

#### Catch rates: annual

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

$~$

#### Catch rates: close up

```{r}
#| label: fig-cpue-weekly
#| eval: true 
#| output: true
#| fig-dpi: 144
#| fig-height: 5
#| echo: false 
#| layout-ncol: 1
#| fig-cap: "CPUE (kg/trap haul) on a weekly basis"

fns = c( 
  "weekly_cpue_smoothed2.png"
)
 
include_graphics( file.path( media_supplementary, fns )  ) 

```


$~$

#### Catch rates: spatial

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

$~$

### At-sea Observed Data

#### At-sea Observed Data: Summary

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
 
#### At-sea Observed fishing locations

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
 
$~$

#### At-sea-observed effort: no. trips

```{r}
#| label: tbl-observer-number-trips
#| eval: true
#| output: true
#| tbl-cap: "Number of at-sea observed trips."

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

#### At-sea-observed effort: trap hauls 

```{r}
#| label: tbl-observer-number-traphauls
#| eval: true
#| output: true
#| tbl-cap: "Number of at-sea observed trap hauls."

odb0$th = paste(odb0$tripset, odb0$lat, odb0$lon)  ## <<< NOTE: needs a re-think 
oo = dcast( odb0[ fishyr>=2004,.(N=length(unique(th))), by=c(vnr, "fishyr")], 
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

#### At-sea-observed number of crab

```{r}
#| label: tbl-observer-number-crab
#| eval: true
#| output: true
#| tbl-cap: "Number of at-sea observed crab." 

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

#### At-sea-observed weight of crab

```{r}
#| label: tbl-observer-weight-crab
#| eval: true
#| output: true
#| tbl-cap: "Total weight of at-sea observed crab (kg)." 

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

#### Carapace condition from observed data: \< 95mm CW

```{r}
#| label: tbl-observer-sublegal  
#| echo: false
#| results: asis
#| eval: true
#| output: true
#| layout-ncol: 1
#| tbl-cap: "Fishery performance statistics: Distribution of at sea observations of males less than 95 mm CW by year and shell condition."
#| tbl-subcap: 
#|   - "(a)"
#|   - "(b)"
#|   - "(c)"
#|   - "(d)"

odb = odb0[ cw < 95 & prodcd_id==0 & shell %in% c(1:5) & get(vnr) %in% regions & sex==0, ]  # male
  
for (r in 1:nregions){
  reg = regions[r]
  REG = reg_labels[r] 
  cat( REG, "\n")
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
 
$~$

#### Carapace condition from observed data: \>= 95mm CW

```{r}
#| label: tbl-observer-legal  
#| echo: false
#| results: asis
#| eval: true
#| output: true
#| layout-ncol: 1
#| tbl-cap: "Fishery performance statistics: Distribution of at sea observations of males greater than 95 mm CW by year and shell condition."
#| tbl-subcap: 
#|   - "(a)"
#|   - "(b)"
#|   - "(c)"
#|   - "(d)"

odb = odb0[ cw >= 95 & cw < 170  & prodcd_id==0 & shell %in% c(1:5) & get(vnr) %in% regions & sex==0, ]  # male
  
for (r in 1:nregions){
  reg = regions[r]
  REG = reg_labels[r]
  cat( REG, "\n")
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
 
$~$

#### Soft-shell

There are two possible definitions:

-   durometer \< 68 (D)
-   carapace conditions 1 and 2 (CC)

```{r}
#| label: tbl-observer-softgt95
#| eval: true
#| output: true
#| echo: false
#| results: asis
#| layout-ncol: 1
#| tbl-cap: "Soft-shell incidence. There are two possible definitions of soft-shelled crab: (D) based on durometer measurements < 68 on the hardness scale; and (CC) based upon classification as carapace conditions 1 and 2." 
#| tbl-subcap: 
#|   - "(a)"
#|   - "(b)"
#|   - "(c)"
#|   - "(d)"

odb = odb0[ cw >= 95 & cw < 170  & prodcd_id==0 & shell %in% c(1:5) & get(vnr) %in% regions & sex==0, ]  # male
shell_condition = odb[ !is.na(get(vnr)), .N, by=c(vnr, "fishyr", "shell") ]
shell_condition[, total:=sum(N, na.rm=TRUE), by=c(vnr, "fishyr")]
shell_condition$percent = round(shell_condition$N / shell_condition$total, 3) * 100
shell_condition$Year = shell_condition$fishyr
 
for (r in 1:nregions){
  reg = regions[r]
  REG = reg_labels[r]
  cat( REG, "\n")

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
 

#### Soft shell locations

```{r}
#| label: fig-observed-softshell-map
#| eval: true 
#| output: true
#| fig-dpi: 144
#| fig-height: 5
#| echo: false 
#| layout-ncol: 1
#| fig-cap: "Map of observed soft shell locations."
#| fig-subcap:  
#|   - "N-ENS"
#|   - "CFA23"
#|   - "CFA24" 

fns = c( 
  "nens_soft_crab_positions_68.png",
  "cfa23_soft_crab_positions_68.png",
  "cfa24_soft_crab_positions_68.png"
)
 
include_graphics( file.path( media_supplementary, fns )  ) 

```
 
$~$

#### At-sea-observed size

```{r}
#| label: fig-observed-size
#| eval: true 
#| output: true
#| fig-dpi: 144
#| fig-height: 5
#| echo: false  
#| fig-cap: "At sea observed crab size"

fns = c( 
  "mean_cw_observed.png"
)
 
include_graphics( file.path( media_supplementary, fns )  ) 

```
 
$~$

#### Size frequency distributions by carapace condition

```{r}
#| label: fig-size-frequency-carapace-condition-observer
#| echo: false
#| results: asis
#| eval: true
#| output: true 
#| fig-dpi: 144
#| fig-height: 6
#| layout-ncol: 2
#! fig-cap: "Size frequency distribution of Snow Crab sampled by At-sea-observers, broken down by Carapace Condition (CC). Vertical lines indicate 95 mm Carapace Width, the minimum legal commercial size."
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

odir = file.path( data_loc, "assessments", year_assessment, "figures", "size.freq", "observer" )
years = year_assessment + c(0:-3) 
for (r in 1:nregions) {
  reg = regions[r]
  REG = reg_labels[r]
  cat( REG, "\n" )
  fns = paste( "size.freq", reg,  years, "png", sep="." )  
  fn = check_file_exists(file.path( odir, fns ) )
  for (ff in fn) show_image(ff) 
  cat("\n\n")
}

```

$~$

#### Total discards of snow crab, by weight

```{r}
#| label: tbl-fishery-discard-total 
#| eval: true
#| output: true
#| echo: false
#| results: asis
#| layout-ncol: 2
#| tbl-cap: "Average by-catch discard rate by weight observed (kg/trap haul; and standard deviation, SD)." 
#| tbl-subcap: 
#|   - "(a)"
#|   - "(b)"
#|   - "(c)"
#|   - "(d)"


for (r in 1:nregions){
  reg = regions[r]
  REG = reg_labels[r]
  cat( REG, "\n")
 
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

$~$


```{r}
#| label: fig-discard_maritimes
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
  labs(x="Year", y="Discard rate of snow crab (Observed, by weight)" ) + 
  theme( legend.position="inside", legend.position.inside=c(0.2, 0.8), legend.title=element_blank() )
  
(pl)

```

#### Bycatch of non-target species

General approach: estimate bycatch from at sea observed data and project
onto marfis data

##### Bycatch of non-target species: estimates based on fisheries **effort**

Bycatch comes from at-sea-observed effort and catch. Rescale these
naively to total snow crab fishery effort.

```{r}
#| label: tbl-fishery-discard-effort
#| eval: true
#| output: true
#| echo: false
#| results: asis
#| layout-ncol: 1
#| tbl-cap: "By-catch (kg) estimated from fisheries effort. Dots indicate low values. Where species exist in a list but there is no data, this indicates some historical bycatch. The overall average is only for the years shown." 
#| tbl-subcap: 
#|   - "(a)"
#|   - "(b)"
#|   - "(c)"
#|   - "(d)"

for (r in 1:nregions){
  reg = regions[r]
  REG = reg_labels[r]
  cat( REG, "\n")
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
 
##### Bycatch of non-target species: estimates based on fisheries **catch**

Bycatch comes from at-sea-observed effort and catch. Rescale these
naively to total snow crab fishery catch.

```{r}
#| label: tbl-fishery-discard-catch 
#| eval: true
#| output: true
#| echo: false
#| results: asis
#| layout-ncol: 1
#| tbl-cap: "By-catch (kg) estimated from landings. Dots indicate low values. Where species exist in a list but there is no data, this indicates some historical bycatch. The overall average is only for the years shown." 
#| tbl-subcap: 
#|   - "(a)"
#|   - "(b)"
#|   - "(c)"
#|   - "(d)"

for (r in 1:nregions){
  reg = regions[r]
  REG = reg_labels[r]
  cat( REG, "\n")
 
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

$~$


```{r}
#| label: fig-bycatch_rates_all
#| eval: false
#| tbl-cap: "At sea observed bycatch rates in Maritimes"
#| fig-dpi: 144
#| fig-height: 15

reg = "cfaall"
o = BC[[reg]]

plot( o$spec ~ o$bct, xlab = "At sea observed catch rate in snow crab fishery (kg/trap)", ylab="Species", type="p", cex=1.1, pch=19, col="darkorange", xlim=c(0, max(o$bct, na.rm=TRUE)*1.4), yaxt="n" )  
text( o$bct, o$spec,  labels=o$species, pos=4, srt=0 , cex=0.8, col="darkslateblue")
text( max(o$bct, na.rm=TRUE)*0.88, 2.5, labels=paste( "Snow crab CPUE (At sea obs., mean): ", o$bct_sc, " kg/trap"), col="darkred", cex=0.9 )
```

##### Entanglements of megafauna

```{r}
#| label: tbl-bycatch_entanglements_all
#| warning: false
#| error: false 
#| tbl-cap: "Entanglements of megafauna in Maritimes"
#| eval: true
#| output: true

region = "cfaall"

# observer data .. extra care needed as there are duplicated records, etc
oss = observer.db( DS="bycatch_clean_data", p=p,  yrs=yrs_observer )  # Prepare at sea observed data
i = polygon_inside( oss[,  c("lon", "lat")], region=region )
oss =  oss[i,]
 
whales = oss[grep("whale", common, ignore.case=TRUE), .(whales=.N), by=.(yr)] 
leatherback = oss[grep("leatherback",  common, ignore.case=TRUE), .(leatherback=.N), by=.(yr)] 
basking_shark = oss[ grep("BASKING SHARK",  common, ignore.case=TRUE), .(basking_shark=.N), by=.(yr)] 

out = data.table( yr=yrs_observer )
out = whales[out, on="yr"]
out = leatherback[out, on="yr"]
out = basking_shark[out, on="yr"]
out[ is.na(out) ] = 0

colnames(out) = c("Year", "Whale", "Leatherback turtle", "Basking shark")

gt::gt(out) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
  summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
  footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
  row_group.padding = gt::px(1))
  
```

##### Map of locations of megafauna entanglements:

```{r}
#| label: fig-bycatch_entanglements_map_all
#| warning: false
#| error: false 
#| fig-cap: "Entanglement locations in Maritimes of megafauna since 2000. Whales (red), Leatherback turtles (green), Basking shark (blue)."
#| fig-dpi: 144
#| fig-height: 5

loc = file.path( data_loc, "output", "maps", "observer.entanglements" )
fn = file.path( loc, "observed_bycatch_entanglements.png" )

include_graphics( fn ) 
   
```

```{=html}
<!--

## To do next:

Estimate via carstm: requires a Poisson model of each species of catch (number) with offset of landings and covariates ... soon, ever?
 
-->
```