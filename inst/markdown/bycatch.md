---
title: "By catch estimation"
author: "Jae S. Choi"
toc: true
number-sections: true
highlight-style: pygments
editor:
  render-on-save: true
format:
  html: 
    code-fold: true
    html-math-method: katex
    embed-resources: true
  pdf:
    pdf-engine: lualatex
  docx: default 
---



   
<!-- Preamble

This is a Markdown document ... To create HTML or PDF, etc, run: 
 
# for presentations to PDF (via beamer):
# note: section separation with '#' can confuse rmarkdown
  make rmarkdown FN=bycatch YR=2023 SOURCE=~/projects/bio.snowcrab/inst/markdown WK=~/bio.data/bio.snowcrab/assessments DOCTYPE=pdf_document DOCEXTENSION=pdf  # {via Rmarkdown}

  --- note: columns only works with beamer_document


# for html documents including presentations:
  make quarto FN=bycatch YR=2023 SOURCE=~/projects/bio.snowcrab/inst/markdown WK=~/bio.data/bio.snowcrab/assessments  DOCEXTENSION=html # {via Quarto}

# for generic pandoc docs
  make pdf FN=bycatch  # {via pandoc}

Alter year and directories to reflect setup or copy Makefile and alter defaults to your needs.
  
  
-->

 
 
<!--  NOTES: Make sure to have pulled observer data:

    year.assessment = 2023
    p = bio.snowcrab::load.environment( year.assessment=year.assessment )
    yrs = 1996:year.assessment # redo all years
    observer.db( DS="rawdata.redo", yrs=yrs )
    observer.db( DS="bycatch.redo", yrs=yrs )  
    observer.db( DS="odb.redo", p=p ) # 3 minutes
    observer.db( DS="bycatch_clean_data.redo", p=p, yrs=p$yrs ) # 3 minutes
 
-->
 
# Approach: estimate bycatch from at seas observed data and project onto marfis data

First set up environment. Can be made more efficient by doing intermediate saves (to do).
 

```{r}
#| eval: true
#| output: false
require(aegis)
year.assessment = 2023
p = bio.snowcrab::load.environment( year.assessment=year.assessment )
require(gt)  # table formatting
require(ggplot2)
bycatch_dir = file.path( p$annual.results, "bycatch")
years = as.character(1996: year.assessment)
if (0) {
  loadfunctions( "aegis")
  loadfunctions( "bio.snowcrab")  # in case of local edits
  require(data.table)
  obs = observer.db( DS="bycatch_clean_data", p=p,  yrs=p$yrs )  # At sea observed data is created on the fly, to access it directly 
}

```


## Naive estimation: directly from observations

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
 


## END 
