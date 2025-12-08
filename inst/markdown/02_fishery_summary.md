---
title: "Snow crab fishery summary" 
keywords: 
  - snow crab trawl survey results
  - basic tables and figures
abstract: |
  Details of 2025 snow crab fishery performance. 

metadata-files:
  - _metadata.yml

params:
  year_assessment: 2025
  year_start: 1999
  data_loc:  "~/bio.data/bio.snowcrab"
  mau: "region"
  debugging: FALSE
  todo: [fishery_results]
---

<!--

# Summary 1 of 4 -- This file is designed to create an HTML document that describes and summarizes the fishery performance. 


# mau = "region" sens as one group
make quarto FN=02_fishery_summary.md YR=2025 MAU=region DATADIR=~/bio.data/bio.snowcrab DOCTYPE=html PARAMS="-P year_assessment:2025 -P mau:region -P todo:[fishery_results]"  --directory=~/bio/bio.snowcrab/inst/markdown 
 
# mau="subarea": sens into 23 and 24 (default behaviour)
make quarto FN=02_fishery_summary.md YR=2025 MAU=subarea DATADIR=~/bio.data/bio.snowcrab DOCTYPE=html PARAMS="-P year_assessment:2025 -P mau:subarea -P todo:[fishery_results]" --directory=~/bio/bio.snowcrab/inst/markdown 
 
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


loc = file.path( data_loc, "output", "maps", "logbook.locations", "recent" )

fns = file.path( loc,  c(     
  "logbook_locations_recent_cfanorth.png",
  "logbook_locations_recent_cfasouth.png",
  "logbook_locations_recent_cfa4x.png"
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


for (r in 1:maus[["n"]]) {
  reg = maus[["internal"]][r]
  REG = maus[["labels"]][r]
  cat(REG, "\n")
  oo = dt[ which( dt$Region==reg), c("Year", "Licenses", "TAC", "Landings", "Effort", "CPUE")]  # Region is hard-coded to mau

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
#| fig-dpi: 144
#| fig-height: 8
#| fig-cap: "Landings (t) of Snow Crab on the SSE. For 4X, the year refers to the starting year of the season. Inset is a closeup view of the timeseries for N-ENS and 4X."
#| fig-subcap: 
#|   - "Annual landings"
#|   - "Cummulative landings"

ts_dir = file.path( data_loc, "assessments", year_assessment, "timeseries", "fishery", params$mau )

fns = c("landings.ts.png", "fisheries_seasonal_cummulative_landings.png")
include_graphics( file.path( ts_dir, fns ) )

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
#| fig-subcap:
#|   - "Annual effort" 
#|   - "Cummulative effort"
 
ts_dir = file.path( data_loc, "assessments", year_assessment, "timeseries", "fishery", params$mau )

fns = c("effort.ts.png", "fisheries_seasonal_cummulative_effort.png")

include_graphics( file.path( ts_dir, fns ) )
 

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
#| fig-subcap: 
#|   - "Annual catch rates"
#|   - "Seasonal catch rates"
 

ts_dir = file.path( data_loc, "assessments", year_assessment, "timeseries", "fishery", params$mau )

fns = c( "cpue.ts.png", "fisheries_seasonal_cpue.png" )

include_graphics( file.path( ts_dir, fns ) ) 
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

 
oo = dcast( odb0[ fishyr>=(year_assessment - 5), .(N=length(unique(tripset))), by=c(mau, "fishyr")], 
  fishyr ~ get(mau), value.var="N", fill=0, drop=FALSE, na.rm=TRUE )
if ( "NA" %in% names(oo) ) oo$"NA" = NULL
keep = c("fishyr", maus[["internal"]])
oo = oo[,..keep]
names(oo) = c("Year", maus[["labels"]] )
oo$Total = rowSums( oo[, 2:maus[["n"]] ], na.rm=TRUE)

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
oo = dcast( odb0[ fishyr>=(year_assessment - 5),.(N=length(unique(th))), by=c(mau, "fishyr")], 
  fishyr ~ get(mau), value.var="N", fill=0, drop=FALSE, na.rm=TRUE )
if ( "NA" %in% names(oo) ) oo$"NA" = NULL
keep = c("fishyr", maus[["internal"]])
oo = oo[,..keep]
names(oo) = c("Year", maus[["labels"]] )
oo$Total = rowSums( oo[, 2:maus[["n"]] ], na.rm=TRUE)

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

oo = dcast( odb0[ fishyr>=(year_assessment - 5),.(N=.N), by=c(mau, "fishyr")], 
  fishyr ~ get(mau), value.var="N", fill=0, drop=FALSE, na.rm=TRUE )
if ( "NA" %in% names(oo) ) oo$"NA" = NULL

keep = c("fishyr", maus[["internal"]])
oo = oo[,..keep]
names(oo) = c("Year", maus[["labels"]] )

oo[, 2:maus[["n"]]] = round( oo[, 2:maus[["n"]]] )
oo$Total = rowSums( oo[, 2:maus[["n"]] ], na.rm=TRUE)

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

oo = dcast( odb0[ fishyr>=(year_assessment - 5),.(N=sum( mass, na.rm=TRUE)/1000 ), by=c(mau, "fishyr")], 
  fishyr ~ get(mau), value.var="N", fill=0, drop=FALSE, na.rm=TRUE )
if ( "NA" %in% names(oo) ) oo$"NA" = NULL

keep = c("fishyr", maus[["internal"]])
oo = oo[,..keep]
names(oo) = c("Year", maus[["labels"]] )

oo[, 2:maus[["n"]]] = round( oo[, 2:maus[["n"]]] )
oo$Total = rowSums( oo[, 2:maus[["n"]] ], na.rm=TRUE)

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

odb = odb0[ fishyr >= (year_assessment - 5) & cw < 95 & prodcd_id==0 & shell %in% c(1:5) & get(mau) %in% maus[["internal"]] & sex==0, ]  # male
  
for (r in 1:maus[["n"]]){
  reg = maus[["internal"]][r]
  REG = maus[["labels"]][r] 
  cat( REG, "\n")
  oo = dcast( odb0[ get(mau)==reg, .(N=.N), by=.(fishyr, shell) ], fishyr  ~ shell, value.var="N", fill=0, drop=FALSE, na.rm=TRUE )
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

odb = odb0[ fishyr >= (year_assessment - 5) & cw >= 95 & cw < 170  & prodcd_id==0 & shell %in% c(1:5) & get(mau) %in% maus[["internal"]] & sex==0, ]  # male
  
for (r in 1:maus[["n"]]){
  reg = maus[["internal"]][r]
  REG = maus[["labels"]][r]
  cat( REG, "\n")
  oo = dcast( odb0[ get(mau)==reg, .(N=.N), by=.(fishyr, shell) ], fishyr  ~ shell, value.var="N", fill=0, drop=FALSE, na.rm=TRUE )
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

odb = odb0[ fishyr >= (year_assessment - 5) & cw >= 95 & cw < 170  & prodcd_id==0 & shell %in% c(1:5) & get(mau) %in% maus[["internal"]] & sex==0, ]  # male
shell_condition = odb[ !is.na(get(mau)), .N, by=c(mau, "fishyr", "shell") ]
shell_condition[, total:=sum(N, na.rm=TRUE), by=c(mau, "fishyr")]
shell_condition$percent = round(shell_condition$N / shell_condition$total, 3) * 100
shell_condition$Year = shell_condition$fishyr
 
for (r in 1:maus[["n"]]){
  reg = maus[["internal"]][r]
  REG = maus[["labels"]][r]
  cat( REG, "\n")

  soft  = odb[ get(mau)==reg & durometer <  68, .(Soft=.N), by=.(fishyr ) ] 
  total = odb[ get(mau)==reg & is.finite(durometer) , .(Total=.N), by=.(fishyr) ] 
  oo = soft[total, on="fishyr"]
  oo = oo[, .(Year=fishyr, Soft=round(Soft/Total*100, 1), Total=Total) ]  
  scond = shell_condition[ get(mau)==reg & shell %in% c(1,2), .(SoftSC=sum(percent), TotalSC=unique(total)[1]), by=.(Year)]
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

ts_dir = file.path( data_loc, "assessments", year_assessment, "timeseries", "observer", params$mau )

fns = c( "cw.png" )
 
include_graphics( file.path( ts_dir, fns )  ) 

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

odir = file.path( data_loc, "assessments", year_assessment, "figures", "size.freq", "observer", params$mau )
years = year_assessment + c(0:-3) 
for (r in 1:maus[["n"]]) {
  reg = maus[["internal"]][r]
  REG = maus[["labels"]][r]
  cat( REG, "\n" )
  fns = paste( "size.freq", reg,  years, "png", sep="." )  
  fn = check_file_exists(file.path( odir, fns ) )
  for (ff in fn) show_image(ff) 
  cat("\n\n")
}

```

$~$

#### Total observed discards (by weight)

```{r}
#| label: tbl-fishery-discard-total 
#| eval: true
#| output: true
#| echo: false
#| results: asis
#| layout-ncol: 2
#| tbl-cap: "Average by-catch discard rate (all species, including snow crab) by weight observed (kg/trap haul; and standard deviation, SD)." 
#| tbl-subcap: 
#|   - "(a)"
#|   - "(b)"
#|   - "(c)"
#|   - "(d)"


for (r in 1:maus[["n"]]){
  reg = maus[["internal"]][r]
  REG = maus[["labels"]][r]
  cat( REG, "\n")
 
  o = BC[[reg]]   
  oo = o$eff_summ[ order(fishyr), ]
  oo = oo[ fishyr >= (year_assessment - 5) , ]

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
#| fig-cap: "At sea observed rate of snow crab discards relative to total catch (discard + kept), all Maritimes, with 1 SD error bars."
#| fig-dpi: 144
#| fig-height: 10
 
oo = NULL
for (reg in maus[["internal"]]) {
  o = BC[[reg]][["eff_summ"]]
  o$region = reg  # region is hard coded
  oo = rbind(oo, o)
}
oo$fishyr = jitter(oo$fishyr, amount=0.02)
oo$region = factor(oo$region, levels=maus[["internal"]], labels=maus[["labels"]], ordered=TRUE) # region hard-coded

color_map = maus[["color_map"]]
shapes = maus[["shapes"]]

pl = ggplot( oo, aes(x=fishyr, y=discard_rate, ymin=discard_rate-discard_rate_sd, ymax=discard_rate+ discard_rate_sd, group=region, fill=region, colour=region) ) +
  geom_line( alpha=0.9, linewidth=1 ) +
  geom_point(aes(shape=region), size=3, alpha=0.8 ) +  # position=position_dodge(width=0.4)
  geom_pointrange()  + # Vertical line with point in the middle
  geom_errorbar(width = 0.2 ) + # Standard error bars
  scale_colour_manual(values=color_map) +
  scale_fill_manual(values=color_map) +
  scale_shape_manual(values = shapes) +
  scale_x_continuous( breaks = seq(2000, 2030, by = 10)) +
  scale_y_continuous( breaks = seq(0, 1, by = 0.25)) +
  facet_wrap(~ region, ncol=1) + 
  labs(x="Year", y="Discard rate (fraction of landings)" ) + 
  theme_light( base_size = 16) + 
  theme( 
    axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.position="top", 
    # legend.position.inside=c(0.9, 0.9), 
    legend.title=element_blank() 
  )
  
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

for (r in 1:maus[["n"]]){
  reg = maus[["internal"]][r]
  REG = maus[["labels"]][r]
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

for (r in 1:maus[["n"]]){
  reg = maus[["internal"]][r]
  REG = maus[["labels"]][r]
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




