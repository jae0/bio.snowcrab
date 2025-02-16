---
title: "Snow Crab Assessment"
subtitle: "Scotian Shelf Ecosystem"
keywords: 
  - snow crab fishery assessment
  - areal-unit modelling of numerical abundance, habitat, mean weight
  - convolution of posterior-simulations of biomass   
abstract: |
  Snow crab modelling of numerical abundance, mean weight and habitat using conditional autoregressive models. Biomass is estimated from these components via multiplication of the posterior draws.   
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
 
make quarto FN=snowcrab_presentation_assessment.md DOCTYPE=revealjs  PARAMS="-P year_assessment:2024"  --directory=~/bio/bio.snowcrab/inst/markdown
 
-->

 

```{r}
#| eval: true
#| output: false
#| echo: false
#| label: setup

  require(knitr)

  opts_chunk$set(
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

## Terms of Reference

- Stock status and trends

- Environmental and climate change considerations

- Bycatch of non-target species

- Major sources of uncertainty, where applicable.

## Management
 
:::: {.columns}

::: {.column width="20%"}

- TAC  
- Season
- Effort
- IBQs
- 100% dockside monitoring
- mandatory logbooks
- At-sea monitoring by certified observers

:::

::: {.column width="5%"}
<!-- empty column to create gap -->
:::

::: {.column width="75%"}

```{r}
#| label: fig-management-areas
#| eval: true
#| echo: false
#| output: true
#| fig-cap: "The Scotian Shelf (NW Atlantic Ocean; NAFO Div. 4VWX). Shown are isobaths and major bathymetric features. Managed Crab Fishing Areas (CFAs; divided by dashed lines) include: NENS, SENS, 4X. SENS is further subdivided (dotted line) into 23 (NW) and 24 (SE)."
#| fig-dpi: 144
#| fig-height: 4


fn = file.path( media_loc, "snowcrab_cfas.png" )
include_graphics( fn ) 

```
:::

::::

## Conservation 

- Females (SSB) are not fished -- completely protected
- Immature crab are not fished and soft-shelled crab are avoided
- Generally conservative exploitation rates
- Spatial refugia in the Marine Protected Areas (MPAs) and continental shelf edge
- Bycatch reduction measures (season timing, biodegradable mesh,and area closures)
- Collaborative co-management -- fisheries-independent, evidence based decision making


## Data sources: Fishery logbooks


```{r}
#| label: fig-fishing-activity
#| eval: true 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig-cap: "All recorded fishing effort on the SSE (ln total effort)."

fn = file.path( media_loc, "snowcrab_fishing.png" )
include_graphics( fn ) 

```


## Effort

- N-ENS: increased 20% in 2024 relative to the previous year and spatially more dispersed. 
- S-ENS: increased 26.6% in 2024 relative to the previous year and spatially more dispersed.
- CFA 4X: decreased 90% in 2024-2025 relative to the previous year and spatially contracted.

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
 


## Landings

- N-ENS: declined marginally (3.9%) and spatially more dispersed 
- S-ENS: declined marginally (0.4%) and spatially more dispersed 
- CFA 4X: decreased (78.3%) but only a small fraction of TACs landed and fishery is ongoing
- TACs for 2024 were 938, 7,345 and 20 t, respectively.

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

## Catch rates

- N-ENS: decreased 20.4% to 2019 levels, declines inshore
- S-ENS: decreased 21.2% to 2017 levels, declines through, especially CFA 24
- CFA 4X: increased 77.3% to 2021 levels, but only a small fraction of TACs landed in Lunenburg area; fishery is ongoing
 
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
 

## Fishery performance: N-ENS
 

```{r}
#| label: tbl-fishery-performance-N-ENS
#| echo: false
#| results: asis
#| eval: true
#| output: true
#| tbl-cap: "Fishery performance statistics: N-ENS"

r=1
  reg = regions[r]
  REG = reg_labels[r]
  # cat(REG, "\n")
  oo = dt[ which(dt$Region==reg), c("Year", "Licenses", "TAC", "Landings", "Effort", "CPUE")] 

  names(oo) = c( "Year", "Licenses", "TAC (t)", "Landings (t)", "Effort (1000 th)", "CPUE (kg/th)" )

  out = gt::gt(oo) |> gt::tab_options(table.font.size = 20, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
  print(out)
 # cat("\n\n")

```
 



## Fishery performance: S-ENS
 

```{r}
#| label: tbl-fishery-performance-S-ENS
#| echo: false
#| results: asis
#| eval: true
#| output: true
#| tbl-cap: "Fishery performance statistics: S-ENS"

r=2
  reg = regions[r]
  REG = reg_labels[r]
  # cat(REG, "\n")
  oo = dt[ which(dt$Region==reg), c("Year", "Licenses", "TAC", "Landings", "Effort", "CPUE")] 

  names(oo) = c( "Year", "Licenses", "TAC (t)", "Landings (t)", "Effort (1000 th)", "CPUE (kg/th)" )

  out = gt::gt(oo) |> gt::tab_options(table.font.size = 20, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
  print(out)
 # cat("\n\n")

```
  
## Fishery performance: 4X
 
```{r}
#| label: tbl-fishery-performance-4X
#| echo: false
#| results: asis
#| eval: true
#| output: true
#| tbl-cap: "Fishery performance statistics: 4X -- years represent the starting year and currently ongoing."

r=3
  reg = regions[r]
  REG = reg_labels[r]
  # cat(REG, "\n")
  oo = dt[ which(dt$Region==reg), c("Year", "Licenses", "TAC", "Landings", "Effort", "CPUE")] 

  names(oo) = c( "Year", "Licenses", "TAC (t)", "Landings (t)", "Effort (1000 th)", "CPUE (kg/th)" )

  out = gt::gt(oo) |> gt::tab_options(table.font.size = 20, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
  print(out)
 # cat("\n\n")


``` 


## Data sources: At Sea Observers

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
   


## Discard of snow crab

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
  
 
## Discard of soft-shelf crab: N-ENS

- Unreliable (low coverage in 2024)

```{r}
#| label: tbl-observer-softgt95-nens
#| eval: true
#| output: true
#| echo: false
#| results: asis
#| layout-ncol: 1
#| tbl-cap: "Soft-shell incidence. There are two possible definitions of soft-shelled crab: (D) based on durometer measurements < 68 on the hardness scale; and (CC) based upon classification as carapace conditions 1 and 2." 


odb = odb0[ cw >= 95 & cw < 170  & prodcd_id==0 & shell %in% c(1:5) & get(vnr) %in% regions & sex==0, ]  # male
shell_condition = odb[ !is.na(get(vnr)), .N, by=c(vnr, "fishyr", "shell") ]
shell_condition[, total:=sum(N, na.rm=TRUE), by=c(vnr, "fishyr")]
shell_condition$percent = round(shell_condition$N / shell_condition$total, 3) * 100
shell_condition$Year = shell_condition$fishyr
 
 r=1
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

```


## Discard of soft-shelf crab: S-ENS


- Unreliable (low coverage in 2024)


```{r}
#| label: tbl-observer-softgt95-sens
#| eval: true
#| output: true
#| echo: false
#| results: asis
#| layout-ncol: 1
#| tbl-cap: "Soft-shell incidence. There are two possible definitions of soft-shelled crab: (D) based on durometer measurements < 68 on the hardness scale; and (CC) based upon classification as carapace conditions 1 and 2." 


odb = odb0[ cw >= 95 & cw < 170  & prodcd_id==0 & shell %in% c(1:5) & get(vnr) %in% regions & sex==0, ]  # male
shell_condition = odb[ !is.na(get(vnr)), .N, by=c(vnr, "fishyr", "shell") ]
shell_condition[, total:=sum(N, na.rm=TRUE), by=c(vnr, "fishyr")]
shell_condition$percent = round(shell_condition$N / shell_condition$total, 3) * 100
shell_condition$Year = shell_condition$fishyr
 
 r=2

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


```


## Discard of soft-shelf crab: 4X

- No coverage in 2024


```{r}
#| label: tbl-observer-softgt95-4x
#| eval: true
#| output: true
#| echo: false
#| results: asis
#| layout-ncol: 1
#| tbl-cap: "Soft-shell incidence. There are two possible definitions of soft-shelled crab: (D) based on durometer measurements < 68 on the hardness scale; and (CC) based upon classification as carapace conditions 1 and 2." 


odb = odb0[ cw >= 95 & cw < 170  & prodcd_id==0 & shell %in% c(1:5) & get(vnr) %in% regions & sex==0, ]  # male
shell_condition = odb[ !is.na(get(vnr)), .N, by=c(vnr, "fishyr", "shell") ]
shell_condition[, total:=sum(N, na.rm=TRUE), by=c(vnr, "fishyr")]
shell_condition$percent = round(shell_condition$N / shell_condition$total, 3) * 100
shell_condition$Year = shell_condition$fishyr
 
 r=3
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
 

```



```{r}
#| label: fig-observed-softshell-map
#| eval: true 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| echo: false 
#| layout-ncol: 2
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
 



 

## Discard of non-target species ("Bycatch"): N-ENS

- average of 0.02% of landings; primarily other Crustacea (crab and lobster)  

```{r}
#| label: tbl-fishery-discard-effort-nens
#| eval: true
#| output: true
#| echo: false
#| results: asis
#| layout-ncol: 1
#| tbl-cap: "Bycatch (kg) estimated from fisheries effort. Dots indicate low values. Where species exist in a list but there is no data, this indicates some historical bycatch. The overall average is only for the years shown." 

r = 1
  reg = regions[r]
  REG = reg_labels[r]
  # cat( REG, "\n")
  o = BC[[reg]]   
  oo = o$bycatch_table_effort
  oo[ oo==0 ] = NA
  oo[ is.na(oo) ] = "."
  out = gt::gt(oo) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
  print(out)
 # cat("\n\n")
 

```
 


## Discard of non-target species ("Bycatch"): S-ENS

- average of 0.03% of landings; primarily other Crustacea (crab and lobster)

```{r}
#| label: tbl-fishery-discard-effort-sens
#| eval: true
#| output: true
#| echo: false
#| results: asis
#| layout-ncol: 1
#| tbl-cap: "Bycatch (kg) estimated from fisheries effort. Dots indicate low values. Where species exist in a list but there is no data, this indicates some historical bycatch. The overall average is only for the years shown." 

r = 2
  reg = regions[r]
  REG = reg_labels[r]
  # cat( REG, "\n")
  o = BC[[reg]]   
  oo = o$bycatch_table_effort
  oo[ oo==0 ] = NA
  oo[ is.na(oo) ] = "."
  out = gt::gt(oo) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
  print(out)
 # cat("\n\n")
 

```
 


## Discard of non-target species ("Bycatch"): 4X

- average of 0.87% of landings; primarily other Crustacea (crab and lobster)  

```{r}
#| label: tbl-fishery-discard-effort-4x
#| eval: true
#| output: true
#| echo: false
#| results: asis
#| layout-ncol: 1
#| tbl-cap: "Bycatch (kg) estimated from fisheries effort. Dots indicate low values. Where species exist in a list but there is no data, this indicates some historical bycatch. The overall average is only for the years shown." 

r = 3
  reg = regions[r]
  REG = reg_labels[r]
  # cat( REG, "\n")
  o = BC[[reg]]   
  oo = o$bycatch_table_effort
  oo[ oo==0 ] = NA
  oo[ is.na(oo) ] = "."
  out = gt::gt(oo) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
  print(out)
 # cat("\n\n")
 

```
 

 

## Bitter crab disease

- Low level background infection rate of < 0.1% of the at-sea observed crab (-@tbl-bcd)
- Spatial distribution is widespread and usually in shallower locations (-@fig-bcd-map) 


```{r}
#| label: tbl-bcd
#| echo: false
#| eval: true
#| output: true
#| tbl-cap: "[Bitter Crab Disease in Maritimes Region](https://www.dfo-mpo.gc.ca/science/aah-saa/diseases-maladies/hematcb-eng.html) is an infection from a dinoflagellate (*Hematodinium*) that causes muscle degeneration. They are widespread (Alaska, NW Atlantic, Greenland) and usually found in warm-water, physiologically stressful conditions. In the Maritimes, it seems to be a low level background infection, found everywhere in the fishing grounds."
#| fig.show: hold 
#| fig-dpi: 144
#| fig-height: 10

include_graphics( file.path( data_loc, "output", "bcd.png") )

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
  
fns = file.path( media_loc, "BCD_map.png"  ) 
include_graphics( fns ) 

``` 

## Ecosystem


## Bottom Temperatures

- Average bottom temperatures observedin the 2024 Snow Crab survey have returned to historical ranges, after worrisome highs in 2022. 

- Bottom temperatures are more stable in N-ENS than S-ENS; 4X exhibits the most erratic and highest annual mean bottom temperatures.


```{r}
#| label: fig-bottom-temperatures-timeseries
#| eval: true
#| echo: false 
#| output: true
#| fig-cap: "Bottom temperatures"
#| fig-subcap: 
#|   - "Annual variations in bottom temperature observed during the Snow Crab survey. The horizontal (black) line indicates the long-term, median temperature within each subarea. Error bars represent standard errors."
#|   - "Posterior densities of predicted average bottom temperatures from an analysis of historical temperature data using [carstm](https://github.com/jae0/carstm). Red horizontal line is at $7^\\circ$C."
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold 

tloc = file.path( data_loc, "assessments", year_assessment, "timeseries"  )

fns = c( 
  file.path("survey", "t.png"), 
  "temperature_bottom.png" 
)

include_graphics( file.path( tloc, fns) )

``` 
 

```{r}
#| label: fig-figures-temperature-bottom-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: "Bottom temperature ($^\\circ$C) observed during the Snow Crab survey."
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""
 
map_outdir = file.path( p$project.outputdir, "maps", "survey", "snowcrab", "annual" )
map_years  = year_assessment + c(0:-3)
  
fn = check_file_exists( file.path( 
  map_outdir, "t", paste( "t", map_years, "png", sep="." )  
) )

include_graphics( fn )
```
 

 

## Interspecific interactions

- Potential prey: in colder (right of Snow Crab). 

- Potential predators: in warmer and deeper (left and bottom of Snow Crab). 

- Water incursions resulted in predators co-occupying Snow Crab habitats while simultaneously losing access to cold water preferring prey. 
 
- Overall, Atlantic Halibut densities have increased rapidly since 2010 (The Gully, Slope Edge and near Sable Island).

- Thorny skate densities have been increasing (especially in N-ENS and along the margins of Banquereau Bank)

- Striped Atlantic Wolffish densities have been high (declining in N-ENS since 2007, peaking in Laurentian Channel).

- Northern shrimp co-occur as they share similar habitat preferences. Numerical densities have declined after a peak in 2011, especially in S-ENS.

- Lesser toad crab is a co-occurring species and potential competitor. Their numbers have declined to low levels throughout, after a peak in densities in 2007 and 2016 in N-ENS.




```{r}
#| label: tbl-predators
#| echo: false
#| eval: true
#| output: true
#| tbl-cap: "Main predators based upon frequency of occuence of snow crab in finfish stomach samples, unadjusted for sampling effort."

gt::gt(counts[1:11,]) 

```




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


## PC1 -- temperature
 
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

include_graphics( fns ) 

``` 

 
```{r}
#| label: fig-speciescomposition-timeseries
#| eval: true
#| echo: false 
#| output: true
#| fig.show: hold 
#| fig-dpi: 144
#| fig-height: 6
#| fig-cap: "Species composition in time. Primary gradient (PC1) is related to bottom temperatures."
 
pc = c(1)

spc_loc = file.path( data_root, "aegis", "speciescomposition", "modelled", "default" )
fnpr = file.path( spc_loc, "figures", paste("pca", pc, "_time.png", sep="" ) )
include_graphics( fnpr ) 

``` 


## PC2 -- depth

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

include_graphics( fns ) 

``` 

 
```{r}
#| label: fig-speciescomposition-timeseries2
#| eval: true
#| echo: false 
#| output: true
#| fig.show: hold 
#| fig-dpi: 144
#| fig-height: 6
#| fig-cap: "Species composition in time. Second axis (PC2) is related to depth."
 
pc = c(2)

spc_loc = file.path( data_root, "aegis", "speciescomposition", "modelled", "default" )
fnpr = file.path( spc_loc, "figures", paste("pca", pc, "_time.png", sep="" ) )
include_graphics( fnpr ) 

``` 


## Predators

"cod", "halibut", "sculpin", "skate", "plaice", "hake", "wolffish", "atlantic cod", "atlantic halibut", "longhorn sculpin", "thorny skate", "striped atlantic wolffish", "haddock", "american plaice", "smooth skate", "winter skate", "white hake", "shorthorn sculpin", "eelpout newfoundland", "squirrel or red hake", "sea raven", "ocean pout", "barndoor skate" 

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
lookup= c( "cod", "halibut", "sculpin", "skate", "plaice", "hake", "wolffish", "atlantic cod", "atlantic halibut", "longhorn sculpin", "thorny skate", "striped atlantic wolffish", "haddock", "american plaice", "smooth skate", "winter skate", "white hake", "shorthorn sculpin", "eelpout newfoundland", "squirrel or red hake", "sea raven", "ocean pout", "barndoor skate", "Squid" )

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
 
- echinoderms
- polychaete worms (*Maldane*, *Nereis*), worm-like animals
- detritus (dead organic matter)
- large zooplankton, shrimp 
- juvenile crab (Rock Crab; Toad Crab; Lesser Toad Crab)
- Ocean Quahog (*Artica islandica*), bivalve molluscs (*Mytilus* sp, *Modiolus*, *Hiatella*)
- brittle stars (*Ophiura*, *Ophiopholis*)
- sea anemones (*Edwardsia*, *Metridium*). 


```{r}
#| label: fig-diet-biplot
#| eval: true
#| echo: false 
#| output: true
#| fig.show: hold 
#| fig-dpi: 144
#| fig-height: 6
#| fig-cap: "Relative location of snow crab prey (green) in the species composition ordination. Snow crab in orange. Most of the potental prey are found to the right of snow crab (i.e. colder-water species) at a variety of depths."
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


## Competitors
 

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

lookup= c( "pandalus", "Jonah Crab", "Atlantic Rock Crab" , "Toad Crab", "Hyas Coarctatus", "Northern Stone Crab" )   # add more here: most are not direct compeitors as they have slightly different depth/temp preferences


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


## Potentially interacting species ("bycatch")
 
```{r}
#| label: fig-competitor-biplot-nens
#| eval: true
#| echo: false 
#| output: true
#| fig.show: hold 
#| fig-dpi: 144
#| fig-height: 6
#| fig-cap: "Bycatch as potentially interacting species in N-ENS (1999-2020). Relative location of interacting species (green) in the species composition ordination. Snow crab in orange."
 
o = BC[["cfanorth"]]    
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

```{r}
#| label: fig-competitor-biplot-sens
#| eval: true
#| echo: false 
#| output: true
#| fig.show: hold 
#| fig-dpi: 144
#| fig-height: 6
#| fig-cap: "Bycatch as potentially interacting species in S-ENS (1999-2020). Relative location of interacting species (green) in the species composition ordination. Snow crab in orange."
 
o = BC[["cfasouth"]]
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
  

```{r}
#| label: fig-competitor-biplot-4x
#| eval: true
#| echo: false 
#| output: true
#| fig.show: hold 
#| fig-dpi: 144
#| fig-height: 6
#| fig-cap: "Bycatch as potentially interacting species in 4X (1999-2020). Relative location of interacting species (green) in the species composition ordination. Snow crab in orange."
 
o = BC[["cfa4x"]]
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
  
  
  



## Viable habitat

- Nouanced: not just temperature, but joint distribution of temperature and other factors (substrate, depths, co-occuring species, ...)

- SSE is variable

- N-ENS: Even with ameliorations in temperatures in 2024, overall habitat viability has declined to near historical lows (peaks in 2011 and 2020). Peaks are in the inner trench and GBH.

- S-ENS: Viable habitat is highest even though temperatures are more stable and cooler in N-ENS. Currently, declined to near historical lows (peaks in 2012). Peaks near Sable and Missaine.

- 4X showed a slight improvement in habitat, however, the overall trend has been downwards since 2010 . Even with amelioration of bottom temperatures in 2024, previous habitat space seems to have been overtaken by competitors and predators in 4X. Peaks in area south of Lunenburg.
 
 

```{r}
#| label: fig-fb-habitat-timeseries
#| eval: true
#| echo: false 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "Habitat viability (probability; fishable Snow Crab). Means and 95\\% Credible Intervals are presented."

loc = file.path( data_loc, "modelled", "default_fb", "aggregated_habitat_timeseries" )
include_graphics( file.path( loc, "habitat_M0.png") )

```


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
  
loc = file.path( data_loc, "modelled", "default_fb", "predicted_habitat" )
vn = "habitat."
yrsplot =  year_assessment + c(0:-9)

fns = file.path( loc, paste( vn, yrsplot, ".png", sep="") )
include_graphics( fns )
 
```

 

## Stock status

## Data sources: Trawl Survey
 
 
- Numerous vessel and captain changes have occured over the years
- 2004: transition from a Spring to Fall survey and mandate transfer from Moncton to BIO
- 2020: due to Covid-19 related to health and safety uncertainties, survey was not conducted 
- 2022: due to mechanical issues with the survey vessel, inshore areas of S-ENS were partially sampled  
- 2024: due to financial issues, areas associated with Marine Protected Areas (St Ann's, Gully) were not sampled
  - N-ENS: 58 stations completed 
  - S-ENS: 282 stations completed
  - CFA 4X: 24 stations completed


```{r}
#| label: fig-survey-locations-map 
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| echo: false 
#| layout-ncol: 2
#| fig-cap: "Survey locations."
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

loc = file.path( data_loc, "output", "maps", "survey.locations" )
years = year_assessment + c(0:-3)
fn = check_file_exists( file.path( loc, paste( "survey.locations", years, "png", sep=".") ))
include_graphics( fn )
```

## Growth

## Recruitment

## Recruitment: males


```{r}
#| label: fig-sizefeq-male
#| eval: true
#| output: true
#| fig-cap: "Size-frequency (areal density; no/km$^2$) histograms by carapace width of male Snow Crab. The vertical line represents the legal size (95 mm). Immature animals are shown with light coloured bars, mature with dark."
#| fig-dpi: 144
#| fig-height: 10

if (params$sens==1) {
  sf_outdir = file.path( p$annual.results, "figures", "size.freq", "survey")
} else if (params$sens==2) {
  sf_outdir = file.path( p$annual.results, "figures", "size.freq", "survey", "split")
}

include_graphics( file.path( sf_outdir,  "male.denl.png" ) )

```
 
- N-ENS: minor recruitment mode of immature male crab centered on 85 mm CW (~instar 10-11 so 1-3 years to entry into fishery) but subsequent year-classes are weak suggesting high natural mortality (predation?); little internal recruitment is expected for the next 3-4 years

- S-ENS: strong and stable size structure suggestive of a stable size/age structure; recruitment expected in 2025 and elevated soft-shell capture will require care

- CFA 4X: erratic inter-annual patterns of low recruitment. The small mode near 68 mm CW (~instar 10 and so 1-3 years toentry to fishery)
  


## Reproduction

- Continued recruitment of female crab into the mature (egg-bearing) stage of the population in all areas

- N-ENS and 4X: decline in numerical densities of both the mature and adolescent females since 2017. 

- S-ENS increase since 2021

- Egg and larval production is expected to be moderate to high in the next year in all areas except N-ENS

- A notable absence of immature females is evident in 4X: predation and habitat loss?

```{r}
#| label: fig-sizefeq-female
#| eval: true
#| output: true
#| fig-cap: "Size-frequency (areal density; no/km$^2$) histograms by carapace width of female Snow Crab. Immature animals are shown with light coloured bars, mature with dark."
#| fig-dpi: 144
#| fig-height: 10


if (params$sens==1) {
  sf_outdir = file.path( p$annual.results, "figures", "size.freq", "survey")
} else if (params$sens==2) {
  sf_outdir = file.path( p$annual.results, "figures", "size.freq", "survey", "split")
}

fn = file.path( sf_outdir, "female.denl.png" )

include_graphics( fn )
```


- N-ENS: densities (~potential egg production) has been low since 2021 and near historical lows; with high densities in inshore and Glace Bay hole
- S-ENS: densities (~potential egg production) declined marginally, especially offshore 
- CFA 4X: densities (~potential egg production) has continued to decline since 2021 and is near historial lows;highest densities close to Sambro and Lunenburg
- Reference: -@fig-totno-female-mat-timeseries, -@fig-totno-female-mat-map


```{r}
#| label: fig-totno-female-mat-timeseries
#| eval: true
#| output: true
#| fig-cap: "The crude, unadjusted geometric mean of mature female density log$_{10}$(no/km$^2$) from the Snow Crab survey. Error bars represent 95\\% Confidence Intervals. Note the absence of data in 2020. Prior to 2004, surveys were conducted in the Spring."
#| fig-dpi: 144
#| fig-height: 4 
 
 
if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

fn = file.path( ts_outdir, paste("totno.female.mat", "png", sep=".") )
include_graphics( fn )

```
 

```{r}
#| label: fig-totno-female-mat-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: "Mature female density log$_{10}$(no/km$^2$) from the Snow Crab survey."
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""
 
map_outdir = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual" )
map_years  = year_assessment + c(0:-3)
  
fn = check_file_exists( file.path( 
  map_outdir, "totno.female.mat", paste( "totno.female.mat", map_years, "png", sep="." )  
) )

include_graphics( fn )
```
  

## Sex ratios

- Unbalanced sex ratios suggest that immature males are exposed to elevated natural mortality (predation, cannibalism, habitat loss?)

   
```{r}
#| label: fig-sexratio-mat-timeseries
#| eval: true
#| output: true
#| fig-cap: "The crude, unadjusted geometric mean of sex ratios (proportion female) of mature Snow Crab. Error bars represent 95\\% Confidence Intervals. Note the absence of data in 2020. Prior to 2004, surveys were conducted in the Spring."
#| fig-dpi: 144
#| fig-height: 4 
 
if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

fn = file.path( ts_outdir, paste("sexratio.mat", "png", sep=".") )
include_graphics( fn )

```
 

 

```{r}
#| label: fig-sexratio-mat-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: "Snow Crab survey sex ratios (proportion female) of mature Snow Crab."
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""
 

map_outdir = file.path( p$project.outputdir, "maps", "survey", "snowcrab", "annual" )
map_years  = year_assessment + c(0:-3)
  
fn = check_file_exists( file.path( 
  map_outdir, "sexratio.mat", paste( "sexratio.mat", map_years, "png", sep="." )  
) )

include_graphics( fn )
```
  


## Fishable biomass density

- N-ENS: densities have declined; declines mostly inshore
- S-ENS: densities have increased; increases mostly offshore
- CFA 4X: densities have declined to historical lows; peak only in Mahonr Bay area (south of Lunenburg)

- Biomass density, however, does not equate to total biomass as the areas occupied by crab can contract, expand and shift with environmental conditions and ecosystem variability. 


```{r}
#| label: fig-R0-timeseries
#| eval: true
#| output: true
#| fig-cap: "The crude, unadjusted geometric mean fishable biomass density log$_{10}$(t/km$^2$) from the Snow Crab survey. Error bars represent 95\\% Confidence Intervals. Note the absence of data in 2020. Prior to 2004, surveys were conducted in the Spring."
#| fig-dpi: 144
#| fig-height: 4 
 
if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

fn = file.path( ts_outdir, paste("R0.mass", "png", sep=".") )
include_graphics( fn )

```
 


```{r}
#| label: fig-R0-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: "Snow Crab survey fishable component biomass density log$_{10}$(t/km$^2$)."
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""
 

map_outdir = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual" )
map_years  = year_assessment + c(0:-3)
  
fn = check_file_exists( file.path( 
  map_outdir, "R0.mass", paste( "R0.mass", map_years, "png", sep="." )  
) )

include_graphics( fn )
```


 



 

## Fishable biomass index

- N-ENS:  three peaks; decline to historical lows; decline throughout area
- S-ENS:  three peaks; increased marginally; increase offshore (marginally, north of Sable)
- CFA 4X: one peak (2009) increased marginally; area south of Sambro and Lunenburg
 

```{r}
#| label: fig-fbindex-timeseries
#| eval: true
#| echo: false 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "The fishable biomass index (t) predicted by CARSTM of Snow Crab survey densities. Error bars represent Bayesian 95\\% Credible Intervals. Note large errors in 2020 when there was no survey."

fn = file.path( data_loc, "modelled", "default_fb", "aggregated_biomass_timeseries" , "biomass_M0.png")
include_graphics( fn )

```
 

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
  
loc = file.path( data_loc, "modelled", "default_fb", "predicted_biomass_densities" )
yrsplot =  year_assessment + c(0:-9)

fns = file.path( loc, paste( "biomass", yrsplot, "png", sep=".") )
include_graphics( fns )
```
  

## Modelled fishable biomass

- N-ENS: prefishery biomass was 2.7t, relative to 3.4t in the previous year
- S-ENS: prefishery biomass was 41.5t, relative to 40.6t in the previous year
- CFA 4X: prefishery biomass was 0.18t, relative to 0.14t inthe previous season
 

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

loc = file.path( data_loc, "fishery_model", year_assessment, "logistic_discrete_historical" )
fns = file.path( loc, c(
  "plot_predictions_cfanorth.png",
  "plot_predictions_cfasouth.png",
  "plot_predictions_cfa4x.png"
) )

include_graphics( fns )

```
 
```{r}
#| label: fig-logisticPredictions2
#| eval: true
#| echo: false
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "Fishable, posterior mean modelled biomass (post-fishery; kt) are shown in dark orange. Light orange are posterior samples of modelled biomass (post-fishery; kt) to illustrate the variability of the predictions. The biomass index (post-fishery, except prior to 2004) after model adjustment by the model catchability coefficient is in gray."
#| fig-subcap:
#|   - "N-ENS"
#|   - "S-ENS"
#|   - "4X"

loc = file.path( data_loc, "fishery_model", year_assessment, "logistic_discrete_historical" )
fns = file.path( loc, c(
  "plot_predictions_postfishery_cfanorth.png",
  "plot_predictions_postfishery_cfasouth.png",
  "plot_predictions_postfishery_cfa4x.png"
) )

include_graphics( fns )

```
 
 

## Posterior estimates of fishing mortality

- N-ENS: F=0.32 (annual exploitation rate of 37%), up from F=0.26 (annual exploitation rate of 30%) in the previous year.

- S-ENS: F=0.17 (annual exploitation rate of 18.6%), F=0.18 (annual exploitation rate of 19%) in the previous year.

- In 4X, the 2024–2025 season (ongoing) is estimated have F=0.052 (annual exploitation rate of 5%), while in the previous season it was F=0.27 (annual exploitation rate of 31%).

- Localized exploitation rates are likely higher, as not all areas for which biomass is estimated are fished (e.g., continental slope areas and western, inshore areas of CFA 24, 4X).



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

- Limit Reference Point (LRP) = K/4
- Upper Stock Reference (USR) = K/2
- Removal Reference (RR) = F_{MSY} = r/ 2

- N-ENS is in the “cautious” zone, though with some important overlap with other zones

- S-ENS is in the “healthy” zone, though with someoverlap with other zones

- 4X is very clearly in the “critical” zone

 
|   | N-ENS | S-ENS | 4X |
|----- | ----- | ----- | ----- |
| |  |  |  |
|q       | `r round(q_north, 3)` (`r round(q_north_sd, 3)`) | `r round(q_south, 3)` (`r round(q_south_sd, 3)`) | `r round(q_4x, 3)` (`r round(q_4x_sd, 3)`) |
|r       | `r round(r_north, 3)` (`r round(r_north_sd, 3)`) | `r round(r_south, 3)` (`r round(r_south_sd, 3)`) | `r round(r_4x, 3)` (`r round(r_4x_sd, 3)`) |
|K       | `r round(K_north, 2)` (`r round(K_north_sd, 2)`) | `r round(K_south, 2)` (`r round(K_south_sd, 2)`) | `r round(K_4x, 2)` (`r round(K_4x_sd, 2)`) |
|Prefishery Biomass   | `r round(B_north[t0], 2)` (`r round(B_north_sd[t0], 2)`) | `r round(B_south[t0], 2)`  (`r round(B_south_sd[t0], 2)`) | `r round(B_4x[t0], 2)`  (`r round(B_4x_sd[t0], 2)`)  |
|Fishing Mortality    | `r round(FM_north[t0], 3)` (`r round(FM_north_sd[t0], 3)`) | `r round(FM_south[t0], 3)` (`r round(FM_south_sd[t0], 3)`) | `r round(FM_4x[t0], 3)` (`r round(FM_4x_sd[t0], 3)`) |

: Reference points from the logistic biomass dynamics fishery model. K is Carrying capacity (kt); and r is Intrinsic rate of increase (non-dimensional). Note that FMSY (fishing mortality associated with 'Maximum Sustainable Yield') is r/2. Similarly, BMSY (biomass associated with 'Maximum Sustainable Yield') is K/2. SD is posterior Standard deviations.* {#tbl-reference-points}
 



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
 

```{r}
#| label: fig-logistic-hcr
#| eval: true
#| echo: false
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "Reference Points (fishing mortality and modelled biomass) from the Fishery Model, for N-ENS (left), S-ENS (middle), and 4X (right). The large yellow dot indicates most recent year and the 95\\% CI. Not: the model does not account for illegal and unreported landings, and interspecific interactions. Prefishery."
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
 



## Other sources of uncertainty
 
- Capture of soft-shell Snow Crab (handling mortality)

- Bycatch of Snow Crab in other fisheries (use as bait)

- Illegal, unreported, and unregulated fishing activities 

- Marine Protected Areas (MPAs) act as a refuge from fishing activities. However, positive effects upon other organisms (predators or prey) can have counter-balancing indirect effects. 

 
## Conclusions  

- The SSE continues to experience rapid ecosystem and climatic variations. Under such conditions, it is prudent to be careful. 

- N-ENS: recruitment continues at low levels. Total mortality exceeded recruitment in 2024. N-ENS is likely in the “cautious” zone.

- S-ENS, recruitment to the fishery continues at a sustainable rate matching total mortality. The S-ENS stock remains in the “healthy” zone.

- 4X, recruitment is expected to be low for another three years. Viable habitat has been minimal for many years. Total mortality is now in approximate balance with recruitment. 4X is in the “critical” zone.


## Acknowledgements

- The Reviewers 

- Innovative Capital Investments Inc.

- Captain Coalie D’Eon and crew of the F/VRS Journey II for their expertise and provision of a safe and hospitable environment for the conduct of the survey.

- Snow Crab license holders and fishers of the SSE a stellar model of a collaborative, sustainable and precautionary co-management of a fishery


## END



## Supplemental Information 
  
 

## Fishery model diagnostics

## Surplus (Shaeffer) production

```{r}
#| label: fig-logistic-surplus-production
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig-cap: "Surplus (Shaeffer) production. Each year is represented by a different colour."
#| fig-subcap:
#|   - "N-ENS"
#|   - "S-ENS"
#|   - "4X"

fns = paste("plot_surplus_production_", regions, ".png", sep="" )
fns = file.path( fm_loc, fns ) 

include_graphics( fns )
   
``` 


## Carrying capacity 

```{r}
#| label: fig-logistic-prior-posterior-K
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig-cap: "Prior-posterior comparisons of carrying capacity (K; kt)."
#| fig-subcap:
#|   - "N-ENS"
#|   - "S-ENS"
#|   - "4X"
   
fns = file.path( fm_loc, paste("plot_prior_K_", regions, ".png", sep="" ) ) 

include_graphics( fns )
   

``` 

## Intrinsic rate of increase

```{r}
#| label: fig-logistic-prior-posterior-r
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig-cap: "Prior-posterior comparisons of th iIntrinsic rate of biomass increase (r)."
#| fig-subcap:
#|   - "N-ENS"
#|   - "S-ENS"
#|   - "4X"

  
fns = file.path( fm_loc, paste("plot_prior_r_", regions, ".png", sep="" ) ) 

include_graphics( fns )

``` 

## Catchability coefficient
 
```{r}
#| label: fig-logistic-prior-posterior-q
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig-cap: "Prior-posterior comparisons of the catchability coefficient (q)."
#| fig-subcap:
#|   - "N-ENS"
#|   - "S-ENS"
#|   - "4X"

  
fns = file.path( fm_loc, paste("plot_prior_q1_", regions, ".png", sep="" ) ) 

include_graphics( fns )

``` 


## Observation error

```{r}
#| label: fig-logistic-prior-posterior-obserror
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig-cap: "Prior-posterior comparisons of Observation error (kt)."
#| fig-subcap:
#|   - "N-ENS"
#|   - "S-ENS"
#|   - "4X"

fns = file.path( fm_loc, paste("plot_prior_bosd_", regions, ".png", sep="" ) ) 

include_graphics( fns )

``` 
 
## Model process error

```{r}
#| label: fig-logistic-prior-posterior-processerror
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig-cap: "Prior-posterior comparisons of Model process error (kt)."
#| fig-subcap:
#|   - "N-ENS"
#|   - "S-ENS"
#|   - "4X"
 
fns = file.path( fm_loc, paste("plot_prior_bpsd_", regions, ".png", sep="" ) ) 

include_graphics( fns )

``` 

## State space 
  
```{r}
#| label: fig-logistic-state-space
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig-cap: "State space (kt): year vs year+1."
#| fig-subcap:
#|   - "N-ENS"
#|   - "S-ENS"
#|   - "4X"

fns = file.path( fm_loc, paste("plot_state_space_", regions, ".png", sep="" ) ) 

include_graphics( fns )

``` 
 
## Climate

## Connectivity: Oceanic currents {.c}

  - Warm, high salinity Gulf Stream from the S-SE along the shelf edge 

  - Cold, low salinity Labrador Current

  - Cold low salinity St. Lawrence outflow from the N-NE

  - Nearshore Nova Scotia current, running from the NE. 
 
```{r ocean_currents, echo=FALSE, out.width='55%', fig.align='center', fig.show='hold',  fig.cap = 'Ocean currents in the Martimes. Source: https://www.dfo-mpo.gc.ca/oceans/publications/soto-rceo/2018/atlantic-ecosystems-ecosystemes-atlantiques/index-eng.html' }
fn2=file.path( media_loc, "maritimes_currents.png" )
include_graphics( c(fn2) ) 
# \@ref(fig:movementtracks)  
``` 


## Rapid climate change

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
include_graphics( fn ) 

```
 
## SST

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

include_graphics( fn ) 

```
 
## Chl-a
  
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
include_graphics( fn ) 

```
 

 
## Bathymetry {.c}
::: columns 

:::: column 

```{r bathymetry-map, out.width='100%', echo=FALSE, fig.align='center', fig.cap='Variations in log(depth; m) in the Scotian Shelf region.' }
bathydir = file.path( data_root, 'aegis', 'bathymetry', 'modelled', 'default', 'stmv', 'none_fft', 'z', 'maps', 'SSE' )
include_graphics( file.path( bathydir, 'bathymetry.z.SSE.png' ) )
```
::::
:::: column
```{r bathymetry-map2, out.width='100%', echo=FALSE, fig.align='center', fig.cap='Local SD log(depth; m) in the Scotian Shelf region.' }
bathydir = file.path( data_root, 'aegis', 'bathymetry', 'modelled', 'default', 'stmv', 'none_fft', 'z', 'maps', 'SSE' )
include_graphics( file.path( bathydir, 'bathymetry.b.sdSpatial.SSE.png' ) )
```

::::
:::


## Substrate {.c}
::: columns 

:::: column 

```{r substrate-map, out.width='100%', echo=FALSE, fig.align='center', fig.cap='Substrate grain size log(mm) variations in the Scotian Shelf region.' }
substrdir = file.path( data_root, 'aegis', 'substrate', 'maps', 'canada.east.highres' )
include_graphics( file.path(  substrdir, 'substrate.substrate.grainsize.canada.east.highres.png' ) )
```
::::
:::: column
```{r substrate-map2, out.width='100%', echo=FALSE, fig.align='center', fig.cap='Local SD substrate grain size log(mm) in the Scotian Shelf region.' }
substrdir = file.path( data_root, 'aegis', 'substrate', 'maps', 'canada.east.highres' )
include_graphics( file.path( substrdir, 'substrate.s.sdSpatial.canada.east.highres.png' ) )
```

::::
:::


## Bottom Temperature {.c} 

```{r bottom-temperatures-survey, out.width='90%', echo=FALSE, fig.align='center', fig.cap = 'Annual variations in bottom temperature observed during the Snow Crab survey. The horizontal (black) line indicates the long-term, median temperature within each subarea. Error bars represent standard errors.' }
include_graphics( file.path( data_loc, 'assessments', year_assessment, 'timeseries', 'survey', 't.png') )
# \@ref(fig:bottom-temperatures-survey)
``` 

```{r bottom-temperatures, out.width='75%', echo=FALSE, fig.align='center', fig.cap = 'Posterior densities of predicted average bottom temperatures. Red horizontal line is at $7^\\circ$C.' }
include_graphics( file.path( data_loc, 'assessments', year_assessment, 'timeseries', 'temperature_bottom.png') )
# \@ref(fig:bottom-temperatures)
``` 

 
```{r bottom-temperatures-map, out.width='30%', echo=FALSE, fig.show='hold', fig.align='center', fig.cap = 'Spatial variations in predicted (1 September) bottom temperature from 2021 (left) to 2023 (right). Groundfish surveys were not conducted in 2020 and 2022 in the snow crab domain. Snow crab surveys were not conducted in 2020, and incomplete in 2022 in S-ENS.' }

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
include_graphics( c( fn3, fn2, fn1) )
# \@ref(fig:bottom-temperatures-map)
# *Spatial variations in bottom temperature estimated from a historical analysis of temperature data for 1 September.*
```
 


## Movement

 
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

include_graphics( fns ) 

``` 

## Habitat preferences
 

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
include_graphics( c( fn2 ) ) 

```
 

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
include_graphics( c(fn1  ) ) 

```

## Co-occurring species


## Atlantic cod

```{r}
#| label: fig-atlcod-timeseries
#| eval: true
#| output: true
#| fig-cap: "Mean density of Atlantic cod log$_{10}$(no/km$^2$) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

species_predator = 10

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```



```{r}
#| label: fig-atlcod-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: Atlantic cod, density; log$_{10}$(no/km$^2$). 
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""
 
map_years  = year_assessment + c(0:-3)
  
species_predator = 10
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```

## Haddock

```{r}
#| label: fig-haddock-timeseries
#| eval: true
#| output: true
#| fig-cap: "Mean density of Haddock log$_{10}$(no/km$^2$) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}
species_predator = 11

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```

 

```{r}
#| label: fig-haddock-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: Haddock, density; log$_{10}$(no/km$^2$). 
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""
 
map_years  = year_assessment + c(0:-3)
  
species_predator = 11
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```



## Halibut

  
```{r}
#| label: fig-halibut-timeseries
#| eval: true
#| output: true
#| fig-cap: "Mean density of Halibut log$_{10}$(no/km$^2$) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

species_predator = 30

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```

 

```{r}
#| label: fig-halibut-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: Halibut, density; log$_{10}$(no/km$^2$). 
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""
 
map_years  = year_assessment + c(0:-3)
  
species_predator = 30
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```

## American plaice
  
```{r}
#| label: fig-amerplaice-timeseries
#| eval: true
#| output: true
#| fig-cap: "Mean density of American plaice log$_{10}$(no/km$^2$) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

species_predator = 40

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```
 

```{r}
#| label: fig-amerplaice-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: American plaice, density; log$_{10}$(no/km$^2$). 
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

map_years  = year_assessment + c(0:-3)
  
species_predator = 40
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```
 
## Striped Alantic Wolffish

```{r}
#| label: fig-stripatlwolffish-timeseries
#| eval: true
#| output: true
#| fig-cap: "Mean density of Striped Atlantic wolffish log$_{10}$(no/km$^2$) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

species_predator = 50

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```

 
```{r}
#| label: fig-stripatlwolffish-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: Striped Atlantic wolffish, density; log$_{10}$(no/km$^2$). 
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

map_years  = year_assessment + c(0:-3)
  
species_predator = 50
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```
 


## Thorny skate

  
```{r}
#| label: fig-thornyskate-timeseries
#| eval: true
#| output: true
#| fig-cap: "Mean density of Thorny skate log$_{10}$(no/km$^2$) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

species_predator = 201

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```

 

```{r}
#| label: fig-thornyskate-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: Thorny skate, density; log$_{10}$(no/km$^2$). 
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

map_years  = year_assessment + c(0:-3)
  
species_predator = 201
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```
 


# Northern shrimp


```{r}
#| label: fig-northernshrimp-timeseries
#| eval: true
#| output: true
#| fig-cap: "Mean density of Northern shrimp log$_{10}$(no/km$^2$) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

species_predator = 2211

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```

 

```{r}
#| label: fig-northernshrimp-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: Northern shrimp, density; log$_{10}$(no/km$^2$). 
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

map_years  = year_assessment + c(0:-3)
  
species_predator = 2211
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```
 

## Jonah crab


```{r}
#| label: fig-jonahcrab-timeseries
#| eval: true
#| output: true
#| fig-cap: "Mean density of Jonah crab log$_{10}$(no/km$^2$) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

species_predator = 2511

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```





```{r}
#| label: fig-jonahcrab-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: Jonah crab, density; log$_{10}$(no/km$^2$). 
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

map_years  = year_assessment + c(0:-3)
  
species_predator = 2511
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```
 


## Lesser toad crab

```{r}
#| label: fig-lyrecrab-timeseries
#| eval: true
#| output: true
#| fig-cap: "Mean density of Arctic Lyre crab (Lesser toad crab) log$_{10}$(no/km$^2$) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

species_predator = 2521

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```





```{r}
#| label: fig-lyrecrab-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: Arctic Lyre crab, density; log$_{10}$(no/km$^2$). 
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

map_years  = year_assessment + c(0:-3)
  
species_predator = 2521
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```
 

## Northern stone crab

```{r}
#| label: fig-nstonecrab-timeseries
#| eval: true
#| output: true
#| fig-cap: "Mean density of Northern stone crab log$_{10}$(no/km$^2$) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

species_predator = 2523

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```


```{r}
#| label: fig-nstonecrab-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: Northern stone crab, density; log$_{10}$(no/km$^2$). 
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

map_years  = year_assessment + c(0:-3)
  
species_predator = 2523
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
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


## Bycatch species in all Maritimes {.c}
 
```{r bycatch-cpue, out.width='70%', echo=FALSE, fig.align='center', fig.cap = 'Bycatch rates in Maritimes Region. Low levels attributable to trap design (top entry, conical, large mesh 5.25" knot-to-knot) permits escapement of non-target species.' }
o = BC[["cfaall"]]
o$bycatch_table[ o$bycatch_table==0 ] = NA
o$bycatch_table[ is.na(o$bycatch_table) ] = "."
o$bycatch_table_catch[ o$bycatch_table_catch==0 ] = NA
o$bycatch_table_catch[ is.na(o$bycatch_table_catch) ] = "."
plot( o$spec ~ o$bct, xlab = "At sea observed catch rate in snow crab fishery (kg/trap)", ylab="Species", type="p", cex=0.9, pch=19, col="darkorange", xlim=c(0, max(o$bct, na.rm=TRUE)*1.4), yaxt="n" )  
text( o$bct, o$spec,  labels=o$species, pos=4, srt=0 , cex=0.5, col="darkslateblue")
text( max(o$bct, na.rm=TRUE)*0.88, 2.5, labels=paste( "Snow crab CPUE (At sea obs., mean): ", o$bct_sc, " kg/trap"), col="darkred", cex=1.0 )
```

 
```{r bycatch-speciesordination_all, out.width='70%', echo=FALSE, fig.align='center', fig.cap = 'Bycatch as potentially interacting species in Maritimes.' }
o = BC[["cfaall"]]
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


 


       
 
## References and further readings
 

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
 
      
tian Shelf and Southern Grand Banks in NAFO Divisions 3NOPs4VWX5Zc. DFO Can. Sci. Advis. Sec. Sci. Resp. 2018/022.

Hebert M, Miron G, Moriyasu M, Vienneau R, and DeGrace P. Efficiency and ghost fishing of Snow Crab (Chionoecetes opilio) traps in the Gulf of St. Lawrence. Fish Res. 2001; 52(3): 143-153. 10.1016/S0165-7836(00)00259-9   
 
Riebler, A., Sørbye, S.H., Simpson D., and Rue, H. 2016. An intuitive Bayesian spatial model for disease mapping that accounts for scaling. Statistical methods in medical research 25: 1145-1165.

Simpson, D., Rue, H., Riebler, A., Martins, T.G., and Sørbye, SH. 2017. Penalising Model Component Complexity: A Principled, Practical Approach to Constructing Priors. Statist. Sci. 32: 1-28.
 
      
