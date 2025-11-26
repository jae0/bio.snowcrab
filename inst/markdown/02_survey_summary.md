---
title: "Snow crab survey summary"
keywords: 
  - snow crab trawl survey results
  - basic tables and figures
abstract: |
  Snow crab demographic structure and environmental conditions. 

metadata-files:
  - _metadata.yml

params:
  year_assessment: 2025
  year_start: 1999
  data_loc:  "~/bio.data/bio.snowcrab"
  sens: 1
  debugging: FALSE
  model_variation: logistic_discrete_historical
  todo: [fishery_results,ecosystem]

---


<!--

# Summary 2 of 4 -- This file is designed to create an HTML document that describes and summarizes the results from the snow crab trawl survey. 


# sens as one group
make quarto FN=02_survey_summary.md YR=2025 DATADIR=~/bio.data/bio.snowcrab DOCTYPE=html PARAMS="-P year_assessment:2025 -P sens:1 -P todo:[fishery_results,ecosystem,redo_data]" --directory=~/bio/bio.snowcrab/inst/markdown
 
# split sens into 23 and 24 (default behaviour)
make quarto FN=02_survey_summary.md YR=2025 DATADIR=~/bio.data/bio.snowcrab DOCTYPE=html PARAMS="-P year_assessment:2025 -P sens:2 -P todo:[fishery_results,ecosystem,redo_data]" --directory=~/bio/bio.snowcrab/inst/markdown

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

 


## Overview of locations and counts


### Map: Survey locations

```{r}
#| label: fig-map-survey-locations 
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
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

 

### Survey counts

```{r}
#| label: tbl-survey-summary
#| echo: false
#| results: asis
#| eval: true
#| output: true
#| layout-ncol: 1
#| tbl-cap: "Summary of trawl survey"
#| tbl-subcap: 
#|   - "(a)"
#|   - "(b)"
#|   - "(c)"
#|   - "(d)"

for (r in 1:nregions) {
  reg = regions[r]
  REG = reg_labels[r]
  cat( REG, "\n")
  oo = set0[ region==reg, .(
    Nstations = .N, 
    Nmale = sum(no.male.all, na.rm=TRUE) ,
    Nfemale = sum(no.female.all, na.rm=TRUE),
    Mmalemature = sum(no.male.mat, na.rm=TRUE),
    Mfemalemature = sum(no.female.mat, na.rm=TRUE)
    ), by=.(yr)]
  oo$Total = oo$Nmale + oo$Nfemale 
  names(oo) = c("Year", "No. stations", "No. male", "No. female", "No. male mature", "No. female mature", "No. total"  )
  oo = oo[order(Year), ]

  out = gt::gt(oo) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
  print(out)
  cat("\n\n")
} 
 
```

### Carapace condition of males >= 95 mm CW

```{r}
#| label: tbl-survey-cc-mature
#| echo: false
#| results: asis
#| eval: true
#| output: true 
#| layout-ncol: 1
#| tbl-cap: "Carapace condition of males >= 95 mm CW captured in the survey"
#| tbl-subcap: 
#|   - "(a)"
#|   - "(b)"
#|   - "(c)"
#|   - "(d)"

det = det0[ cw >= 95 ,]  # commercial sized crab only
 
for (r in 1:nregions) {
  reg = regions[r]
  REG = reg_labels[r]
  cat( REG, "\n")
  oo = dcast( det[ region==reg & !is.na(shell), .(N=.N), by=.(fishyr, shell) ], fishyr  ~ shell, value.var="N", fill=0, drop=FALSE, na.rm=TRUE )

  keep = c("fishyr",  "1", "2", "3", "4", "5" )
  oo = oo[,..keep]
  names(oo) = c("Year", "CC1", "CC2", "CC3", "CC4", "CC5" )

  oo$Total = rowSums( oo[, 2:6 ], na.rm=TRUE)
  oo[, 2:6 ] = round(oo[, 2:6 ] / oo$Total * 100, digits=1)
  oo = oo[order(Year), ]
  out = gt::gt(oo) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
  print(out)
  cat("\n\n")
} 


```



### Size structure by carapace condition of mature males


```{r}
#| label: fig-sizefeq-male-mature-survey-cc
#| echo: false
#| results: asis
#| eval: true
#| output: true 
#| fig-dpi: 144
#| fig-height: 6
#| fig-cap: "Size-structure and carapace condition of mature males."

odir = file.path( data_loc, "assessments", year_assessment, "figures", "size.freq", "carapacecondition" )
years = year_assessment + c(0:-3) 
cat("$~$ \n\n")
for (r in 1:nregions) {
  reg = regions[r]
  REG = reg_labels[r]
  cat( REG, "\n\n" )
  fns = paste( "sizefreq", reg, years, "png", sep="." ) 
  fn = check_file_exists(file.path( odir, fns ) )
  for (ff in fn) show_image(ff)  
  cat("\n\n")
}

``` 

$~$

 

### Size structure and maturity of males


```{r}
#| label: fig-sizefeq-male
#| eval: true
#| output: true
#| fig-cap: "Size-frequency (areal density; no/km$^2$) histograms by carapace width of male Snow Crab. The vertical line represents the legal size (95 mm). Immature animals are shown with light coloured bars, mature with dark."
#| fig-dpi: 144
#| fig-height: 10

if (params$sens==1) {
  sf_outdir = file.path( p$annual.results, "figures", "size.freq", "survey", "period1")
} else if (params$sens==2) {
  sf_outdir = file.path( p$annual.results, "figures", "size.freq", "survey_split", "period1")
}

include_graphics( file.path( sf_outdir,  "male.denl.png" ) )

```



### Size structure maturity of females

 
```{r}
#| label: fig-sizefeq-female
#| eval: true
#| output: true
#| fig-cap: "Size-frequency (areal density; no/km$^2$) histograms by carapace width of female Snow Crab. Immature animals are shown with light coloured bars, mature with dark."
#| fig-dpi: 144
#| fig-height: 10

if (params$sens==1) {
  sf_outdir = file.path( p$annual.results, "figures", "size.freq", "survey", "period1")
} else if (params$sens==2) {
  sf_outdir = file.path( p$annual.results, "figures", "size.freq", "survey_split", "period1")
}


fn = file.path( sf_outdir, "female.denl.png" )

include_graphics( fn )
```

 


## Timeseries and maps of survey variables of interest


### Bottom temperature: trawl survey


```{r}
#| label: fig-temperature-bottom-ts
#| eval: true
#| output: true
#| fig-cap: "Annual variations in bottom temperature ($^\\circ$C) observed during the Snow Crab survey. The horizontal (black) line indicates the long-term, median temperature within each subarea. Error bars represent standard errors."
#| fig-dpi: 144
#| fig-height: 4 

if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

fn = file.path( ts_outdir, paste("t", "png", sep=".") )

include_graphics( fn )

```

$~$


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
 
$~$

### Sex ratios


```{r}
#| label: fig-sexratio-mat-ts
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

$~$


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
 
$~$


### Mature female
 

```{r}
#| label: fig-totno-female-mat-ts
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

$~$


### Fishable biomass 

```{r}
#| label: fig-R0-ts
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

$~$


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


$~$
    
### Fishable mean size 

```{r}
#| label: fig-cw-male-mat-ts
#| eval: true
#| output: true
#| fig-cap: "Mean size of mature male Snow Crab log$_{10}$(CW; mm) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 


if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

fn = file.path( ts_outdir, paste("cw.male.mat.mean", "png", sep=".") )
include_graphics( fn )

```

$~$
    

```{r}
#| label: fig-cw-male-mat-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: "Snow Crab survey fishable component mean carapace width; log$_{10}$(CW; mm)."
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""
 

map_outdir = file.path( p$project.outputdir, "maps", "survey", "snowcrab", "annual" )
map_years  = year_assessment + c(0:-3)
  
fn = check_file_exists( file.path( 
  map_outdir, "cw.male.mat.mean", paste( "cw.male.mat.mean", map_years, "png", sep="." )  
) )

include_graphics( fn )
```


$~$
   

 
## Potential Predators 

The main predators, based on literature and stomach content analysis, are: 

cod, haddock, halibut, plaice, wolfish, thornyskate, smoothskate, winterskate.

### Atlantic cod

  
```{r}
#| label: fig-atlcod-ts
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

$~$
   

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
 
$~$
   

### Haddock

  
```{r}
#| label: fig-haddock-ts
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

$~$
   

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


$~$
     
     
 

### Halibut

  
```{r}
#| label: fig-halibut-ts
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

$~$
   
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
 
$~$


 

### American plaice

  
```{r}
#| label: fig-amerplaice-ts
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

$~$
   

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
 
$~$

### Striped Atlantic wolffish

  
```{r}
#| label: fig-stripatlwolffish-ts
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

$~$
   

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
 
$~$
 
### Thorny skate

  
```{r}
#| label: fig-thornyskate-ts
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

$~$
   

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
 
$~$
 
### Smooth skate

  
```{r}
#| label: fig-smoothskate-ts
#| eval: true
#| output: true
#| fig-cap: "Mean density of Smooth skate log$_{10}$(no/km$^2$) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

species_predator = 202

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```

$~$
   

```{r}
#| label: fig-smoothskate-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: Smooth skate, density; log$_{10}$(no/km$^2$). 
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

map_years  = year_assessment + c(0:-3)
  
species_predator = 202
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```
 
$~$
 
### Winter skate

  
```{r}
#| label: fig-winterskate-ts
#| eval: true
#| output: true
#| fig-cap: "Mean density of Winter skate log$_{10}$(no/km$^2$) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

species_predator = 204

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```

$~$
   

```{r}
#| label: fig-winterskate-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: Winter skate, density; log$_{10}$(no/km$^2$). 
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

map_years  = year_assessment + c(0:-3)
  
species_predator = 204
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```
 

$~$
  


## Potential Competitors

The main potential predators, based on literature and overlpping distributions are: 

northernshrimp, jonahcrab, lessertoadcrab.
   

### Northern shrimp

  
```{r}
#| label: fig-northernshrimp-ts
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

$~$
   

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
 
$~$
  

### Jonah crab

Not exactly a competitor. Similar habitat except warmer areas so more an indicator of bottom temperatures.
  
```{r}
#| label: fig-jonahcrab-ts
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

$~$
   

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
 
$~$
  

### Arctic Lyre crab (Lesser toad crab)
 
Slightly more shallow environments than snow crab.


```{r}
#| label: fig-lyrecrab-ts
#| eval: true
#| output: true
#| fig-cap: "Mean density of Arctic Lyre crab log$_{10}$(no/km$^2$) from surveys with 95\\% Confidence Intervals."
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

$~$
   

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
 
$~$

   
### Northern stone crab
  
```{r}
#| label: fig-nstonecrab-ts
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

$~$
   

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
 
$~$

   
  
   
 <!--

  # ------------------------------------------
  # Timeseries: Larval brachyura from the SSIP data
  ##figure.timeseries.larvae( outdir=file.path(p$project.outputdir, "timeseries", "larvae") )

  # ------------------------------------------
  # Growth as a a function of instar for Scotian Shelf snow crab
  figure.growth.instar( outdir=file.path(p$project.outputdir, "growth") )


  # ------------------------------------------
  # Map: Larval distributions from the Scotian Shelf Ichtyoplankton Program data
  map.larvae( p=p, outdir=file.path(p$project.outputdir, "maps", "larvae"), conversions=conversions )


  # ------------------------------------------
  # Map: Spatial representation of maturity patterns of snow crab
  #MG Not sure we use these maps either, check with Adam and Jae
  # map.maturity( p, outdir=file.path(p$project.outputdir, "maps", "maturity"), newyear=T )

  res = maturity_region_year(p)  # timeseries of maturity

 


}
-->
