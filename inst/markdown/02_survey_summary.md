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
  mau: "region"
  debugging: FALSE
  todo: [survey]

---


<!--

# Summary 2 of 4 -- This file is designed to create an HTML document that describes and summarizes the results from the snow crab trawl survey. 


# sens as one group
make quarto FN=02_survey_summary.md YR=2025 MAU=region DATADIR=~/bio.data/bio.snowcrab DOCTYPE=html PARAMS="-P year_assessment:2025 -P mau:region -P todo:[survey]" --directory=~/bio/bio.snowcrab/inst/markdown
 
# split sens into 23 and 24 (default behaviour)
make quarto FN=02_survey_summary.md YR=2025 MAU=subarea DATADIR=~/bio.data/bio.snowcrab DOCTYPE=html PARAMS="-P year_assessment:2025 -P mau:subarea -P todo:[survey]" --directory=~/bio/bio.snowcrab/inst/markdown

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

for (r in 1:maus[["n"]]) {
  reg = maus[["internal"]][r]
  REG = maus[["labels"]][r]
  cat( REG, "\n")
  oo = set0[ set0[[mau]]==reg & yr >= (year_assessment - 5) , .(
    Nstations = .N, 
    Nmale = sum(no.male.all, na.rm=TRUE) ,
    Nfemale = sum(no.female.all, na.rm=TRUE),
    Mmalemature = sum(no.male.mat, na.rm=TRUE),
    Mfemalemature = sum(no.female.mat, na.rm=TRUE)
    ), by=.(yr)]
  oo$Total = oo$Nmale + oo$Nfemale 
  names(oo) = c("Year", "No. stations", "No. male", "No. female", "No. male mature", "No. female mature", "No. total"  )
  oo = oo[order(Year), ]

  print( table_format_simple(oo) )
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
 
for (r in 1:maus[["n"]]) {
  reg = maus[["internal"]][r]
  REG = maus[["labels"]][r]
  cat( REG, "\n")
  oo = dcast( det[ det[[mau]]==reg & !is.na(shell), 
    .(N=.N), by=.(fishyr, shell) ], fishyr  ~ shell, value.var="N", fill=0, drop=FALSE, na.rm=TRUE )

  keep = c("fishyr",  "1", "2", "3", "4", "5" )
  oo = oo[,..keep]
  names(oo) = c("Year", "CC1", "CC2", "CC3", "CC4", "CC5" )

  oo$Total = rowSums( oo[, 2:6 ], na.rm=TRUE)
  oo[, 2:6 ] = round(oo[, 2:6 ] / oo$Total * 100, digits=1)
  oo = oo[order(Year), ]
  oo = oo[ Year >= (year_assessment - 5) , ]

  print( table_format_simple(oo) )
  cat("\n\n")
} 


```



### Size structure by carapace condition of mature males


```{r}
#| label: fig-sizefeq-male-mature-survey-cc
#| echo: false
#| eval: true
#| output: true 
#| fig-dpi: 144
#| fig-height: 6
#| fig-cap: "Size-structure and carapace condition of mature males."
 

fn_grouped = file.path( data_loc, "assessments", year_assessment, paste( "size_freq_survey_", maus[["internal"]], ".png", sep="" ) )  # to create new files

odir = file.path( data_loc, "assessments", year_assessment, "figures", "size.freq", "carapacecondition"  )

years = year_assessment + c(0:-1) 

# save a modified version to disk 
for (r in 1:maus[["n"]]) {
  fns = paste( "sizefreq", maus[["internal"]][r], years, "png", sep="." ) 
  fn = check_file_exists(file.path( odir, fns ) )
  img = magick::image_append( c( 
      magick::image_read(fn[1]), 
      magick::image_read(fn[2]) 
    ))
  img = magick::image_annotate( img, maus[["labels"]][r], size=140, location="+0+0", font="Open Sans" )
  magick::image_write( img, path=fn_grouped[r], format="png" )
}

include_graphics(fn_grouped)
 

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


sf_outdir = file.path( p$annual.results, "figures", "size.freq", "survey", params$mau, "period1")

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



sf_outdir = file.path( p$annual.results, "figures", "size.freq", "survey", params$mau, "period1")

fn = file.path( sf_outdir, "female.denl.png" )

include_graphics( fn )
```

 


## Timeseries and maps of survey variables of interest


### Bottom temperature: trawl survey


```{r}
#| label: fig-temperature-bottom-ts
#| eval: true
#| output: true
#| fig-cap: "Annual variations in bottom temperature ($^\\circ$C) observed during the Snow Crab survey. The horizontal (black) line indicates the long-term, median temperature within each subarea. Error bars represent 1 standard error."
#| fig-dpi: 144
#| fig-height: 4 


ts_outdir = file.path( p$annual.results, "timeseries", "survey", params$mau )

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
#| fig-cap: "The crude, unadjusted geometric mean of sex ratios (proportion female) of mature Snow Crab. Error bars represent 1 SE. Note the absence of data in 2020. Prior to 2004, surveys were conducted in the Spring."
#| fig-dpi: 144
#| fig-height: 4 
 

ts_outdir = file.path( p$annual.results, "timeseries", "survey", params$mau)

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
#| fig-cap: "The crude, unadjusted geometric mean of mature female density log$_{10}$(no/km$^2$) from the Snow Crab survey. Error bars represent 1 SE. Note the absence of data in 2020. Prior to 2004, surveys were conducted in the Spring."
#| fig-dpi: 144
#| fig-height: 4 
 
 
ts_outdir = file.path( p$annual.results, "timeseries", "survey", params$mau)
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
#| fig-cap: "The crude, unadjusted geometric mean fishable biomass density log$_{10}$(t/km$^2$) from the Snow Crab survey. Error bars represent 1 SE. Note the absence of data in 2020. Prior to 2004, surveys were conducted in the Spring."
#| fig-dpi: 144
#| fig-height: 4 
 
 
ts_outdir = file.path( p$annual.results, "timeseries", "survey", params$mau)
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
    
### Mature male mean size (survey) 

```{r}
#| label: fig-cw-male-mat-ts
#| eval: true
#| output: true
#| fig-cap: "Mean size of mature male Snow Crab log$_{10}$(CW; mm) from surveys with 1 SE."
#| fig-dpi: 144
#| fig-height: 4 


ts_outdir = file.path( p$annual.results, "timeseries", "survey", params$mau)
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
#| fig-cap: "Mean density of Atlantic cod log$_{10}$(no/km$^2$) from surveys with 1 SE."
#| fig-dpi: 144
#| fig-height: 4 


ts_outdir = file.path( p$annual.results, "timeseries", "survey", params$mau)

species_code = 10

bc_vars = paste("ms.no", species_code, sep='.')
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
  
species_code = 10
bc_vars = paste("ms.no", species_code, sep='.')
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
#| fig-cap: "Mean density of Haddock log$_{10}$(no/km$^2$) from surveys with 1 SE."
#| fig-dpi: 144
#| fig-height: 4 


ts_outdir = file.path( p$annual.results, "timeseries", "survey", params$mau)
species_code = 11

bc_vars = paste("ms.no", species_code, sep='.')
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
  
species_code = 11
bc_vars = paste("ms.no", species_code, sep='.')
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
#| fig-cap: "Mean density of Halibut log$_{10}$(no/km$^2$) from surveys with 1SE."
#| fig-dpi: 144
#| fig-height: 4 


ts_outdir = file.path( p$annual.results, "timeseries", "survey", params$mau)
species_code = 30

bc_vars = paste("ms.no", species_code, sep='.')
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
  
species_code = 30
bc_vars = paste("ms.no", species_code, sep='.')
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
#| fig-cap: "Mean density of American plaice log$_{10}$(no/km$^2$) from surveys with 1 SE."
#| fig-dpi: 144
#| fig-height: 4 

ts_outdir = file.path( p$annual.results, "timeseries", "survey", params$mau)

species_code = 40

bc_vars = paste("ms.no", species_code, sep='.')
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
  
species_code = 40
bc_vars = paste("ms.no", species_code, sep='.')
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
#| fig-cap: "Mean density of Striped Atlantic wolffish log$_{10}$(no/km$^2$) from surveys with 1 SE."
#| fig-dpi: 144
#| fig-height: 4 

  ts_outdir = file.path( p$annual.results, "timeseries", "survey", params$mau)

species_code = 50

bc_vars = paste("ms.no", species_code, sep='.')
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
  
species_code = 50
bc_vars = paste("ms.no", species_code, sep='.')
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
#| fig-cap: "Mean density of Thorny skate log$_{10}$(no/km$^2$) from surveys with 1 SE."
#| fig-dpi: 144
#| fig-height: 4 

  ts_outdir = file.path( p$annual.results, "timeseries", "survey", params$mau)

species_code = 201

bc_vars = paste("ms.no", species_code, sep='.')
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
  
species_code = 201
bc_vars = paste("ms.no", species_code, sep='.')
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
#| fig-cap: "Mean density of Smooth skate log$_{10}$(no/km$^2$) from surveys with 1 SE."
#| fig-dpi: 144
#| fig-height: 4 

  ts_outdir = file.path( p$annual.results, "timeseries", "survey", params$mau)

species_code = 202

bc_vars = paste("ms.no", species_code, sep='.')
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
  
species_code = 202
bc_vars = paste("ms.no", species_code, sep='.')
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
#| fig-cap: "Mean density of Winter skate log$_{10}$(no/km$^2$) from surveys with 1 SE."
#| fig-dpi: 144
#| fig-height: 4 

  ts_outdir = file.path( p$annual.results, "timeseries", "survey", params$mau)

species_code = 204

bc_vars = paste("ms.no", species_code, sep='.')
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
  
species_code = 204
bc_vars = paste("ms.no", species_code, sep='.')
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
#| fig-cap: "Mean density of Northern shrimp log$_{10}$(no/km$^2$) from surveys with 1 SE."
#| fig-dpi: 144
#| fig-height: 4 

  ts_outdir = file.path( p$annual.results, "timeseries", "survey", params$mau)

species_code = 2211

bc_vars = paste("ms.no", species_code, sep='.')
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
  
species_code = 2211
bc_vars = paste("ms.no", species_code, sep='.')
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
#| fig-cap: "Mean density of Jonah crab log$_{10}$(no/km$^2$) from surveys with 1 SE."
#| fig-dpi: 144
#| fig-height: 4 

  ts_outdir = file.path( p$annual.results, "timeseries", "survey", params$mau)

species_code = 2511

bc_vars = paste("ms.no", species_code, sep='.')
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
  
species_code = 2511
bc_vars = paste("ms.no", species_code, sep='.')
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
#| fig-cap: "Mean density of Arctic Lyre crab log$_{10}$(no/km$^2$) from surveys with 1 SE."
#| fig-dpi: 144
#| fig-height: 4 

  ts_outdir = file.path( p$annual.results, "timeseries", "survey", params$mau)

species_code = 2521

bc_vars = paste("ms.no", species_code, sep='.')
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
  
species_code = 2521
bc_vars = paste("ms.no", species_code, sep='.')
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
#| fig-cap: "Mean density of Northern stone crab log$_{10}$(no/km$^2$) from surveys with 1 SE."
#| fig-dpi: 144
#| fig-height: 4 

  ts_outdir = file.path( p$annual.results, "timeseries", "survey", params$mau)

species_code = 2523

bc_vars = paste("ms.no", species_code, sep='.')
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
  
species_code = 2523
bc_vars = paste("ms.no", species_code, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```
 
$~$

   


   
### Northern sandlance
  
```{r}
#| label: fig-sandlance-ts
#| eval: true
#| output: true
#| fig-cap: "Mean density of Sandlance log$_{10}$(no/km$^2$) from surveys with 1 SE."
#| fig-dpi: 144
#| fig-height: 4 

  ts_outdir = file.path( p$annual.results, "timeseries", "survey", params$mau)

species_code = 610

bc_vars = paste("ms.no", species_code, sep='.')
# fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```

$~$
   

```{r}
#| label: fig-sandlance-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: Sandlance, density; log$_{10}$(no/km$^2$). 
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

map_years  = year_assessment + c(0:-3)
  
species_code = 610
bc_vars = paste("ms.no", species_code, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

# fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```
 
$~$

   
 

   
### Capelin
  
```{r}
#| label: fig-Capelin-ts
#| eval: true
#| output: true
#| fig-cap: "Mean density of Capelin log$_{10}$(no/km$^2$) from surveys with 1 SE."
#| fig-dpi: 144
#| fig-height: 4 

  ts_outdir = file.path( p$annual.results, "timeseries", "survey", params$mau)

 species_code = 64
 
  bc_vars = paste("ms.no", species_code, sep='.')
  fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
# include_graphics( fn )

```

$~$
   

```{r}
#| label: fig-Capelin-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: Capelin, density; log$_{10}$(no/km$^2$). 
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

map_years  = year_assessment + c(0:-3)
  
species_code = 64
bc_vars = paste("ms.no", species_code, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
# include_graphics( fn )
    
```
 
$~$

   
 
   
### Winter skate purse
  
```{r}
#| label: fig-winterskatepurse-ts
#| eval: true
#| output: true
#| fig-cap: "Mean density of winter skate purse log$_{10}$(no/km$^2$) from surveys with 1 SE."
#| fig-dpi: 144
#| fig-height: 4 

ts_outdir = file.path( p$annual.results, "timeseries", "survey", params$mau)

species_code = bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=1204 ) # 1192 (skate purses generic)

bc_vars = paste("ms.no", species_code, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```

$~$
   

```{r}
#| label: fig-winterskatepurse-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: Winter skate purse, density; log$_{10}$(no/km$^2$). 
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

map_years  = year_assessment + c(0:-3)
  
species_code = bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=1204 ) # 1192 (skate purses generic)

bc_vars = paste("ms.no", species_code, sep='.')
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
