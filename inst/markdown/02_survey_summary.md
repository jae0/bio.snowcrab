---
title: "Snow crab trawl survey -- figures and tables"
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
  Snow crab demographic structure and environmental conditions. 
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
    embed-resources: true
  pdf:
    pdf-engine: lualatex
  docx: default 
  beamer:
    pdf-engine: lualatex
---



# Snow crab trawl survey -- figures and tables

## Preamble

This is a markdown document. It can be viewed in formatted form via:
  
  - a web browser open the webpage: [02_survey_summary.md](https://github.com/jae0/bio.snowcrab/tree/master/inst/markdown/02_survey_summary.md)   

  - a web browser open the local file directly: [02_survey_summary.md](../markdown/02_survey_summary.md) (you might need to install a browser add-in), or 
  
  - a text editor with markdown awareness (e.g. Rstudio, VSCode, etc. ).


As this document uses the Quarto and Rmarkdown dialect of Markdown, you can  create reports, by running one of the following in a shell command prompt or using Rstudio. HTML is the most versatile at the moment: 


```shell
 
# {via Quarto}
cd ~/bio/bio.snowcrab/inst/markdown

make quarto FN=02_survey_summary YR=2024 SOURCE=~/bio/bio.snowcrab/inst/markdown WK=~/bio.data/bio.snowcrab/assessments DOCEXTENSION=html 
 

```

Or, see the Makefile and alter defaults to your needs. As Quarto does not pass params easily. So you must adjust "params" in yaml at the top of this fle, or use to quarto command such as: 

```shell
quarto ... -P year.assessment:$(YR) -P media_loc:$(MEDIA) 
```

[See here for more YAML options.](https://quarto.org/docs/output-formats/all-formats.html)
 
  

## Set up environment

  - Ensure the **year.assessment** is correct

  - Ensure that the data pulls of survey data stored on the ISSDB is complete and assimilated to the end of [01_snowcrab_data.md](https://github.com/jae0/bio.snowcrab/tree/master/inst/markdown/01_snowcrab-data.md).
  
 


```{r}
#| eval: true
#| output: false
#| label: setup

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
require(MBA)

require(aegis)  # basic helper tools
 
year.assessment = 2024  # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  change this as appropriate

year_previous = year.assessment - 1

p = bio.snowcrab::load.environment( year.assessment=year.assessment )  

SCD = project.datadirectory("bio.snowcrab")
media_loc = project.codedirectory("bio.snowcrab", "inst", "markdown", "media")

require(gt)  # table formatting

outtabledir = file.path( p$annual.results, "tables" )

years = as.character(1996: year.assessment)

regions = c("cfanorth", "cfasouth", "cfa4x")
nregions = length(regions)


yrs = 1996:year.assessment # redo all years

loadfunctions( "aegis")
loadfunctions( "bio.snowcrab")  # in case of local edits
 
p$corners = data.frame(plon=c(220, 990), plat=c(4750, 5270) )

p$mapyears = year.assessment + c(-5:0 )   # default in case not specified

  
```


## Overview of locations and counts


### Map: Survey locations

```{r}
#| eval: true
#| output: true
#| label: map-survey-locations 
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| echo: false 
#| layout-ncol: 2
  
loc = file.path( SCD, "output", "maps", "survey.locations" )
years = year.assessment + c(0:-3)
   
quietly(
  map.survey.locations( p=p, basedir=loc, years=years )
  # map.survey.locations( p=p, basedir=loc,  years=years, map.method="googleearth"  )
)

fn = check_file_exists( file.path( loc, paste( "survey.locations", years, "png", sep=".") ))
include_graphics( fn )
```
Survey locations.
 

### Counts of stations in each area

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
gt::gt(out) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
  summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
  footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
  row_group.padding = gt::px(1))
```

 
 

## Size frequency distributions by carapace condition


```{r}
#| label: create-size-frequency-moult-carapace-condition
#| eval: true
#| output: true

quietly( suppressWarnings(
  figure.sizefreq.carapacecondition( X=snowcrab.db( p=p, DS="det.georeferenced" ), cwbr=4, regions=c("cfanorth", "cfasouth", "cfa4x"), 
    outdir=file.path( p$annual.results, "figures", "size.freq", "carapacecondition" )  ) 
))

```

### Males \>= 95mm CW

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

NENS:

```{r}
#| eval: true
#| output: true
#| label: table-survey-nens-comm
#| tbl-cap: "Distribution of NENS survey: males less than 95 mm CW by year and shell condition."

resN = dcast( det[ region=="cfanorth" & !is.na(shell), .(N=.N), by=.(fishyr, shell) ], fishyr  ~ shell, value.var="N", fill=0, drop=FALSE, na.rm=TRUE )
names(resN) = c("Year", "CC1", "CC2", "CC3", "CC4", "CC5" )
resN$Total = rowSums( resN[, 2:6 ], na.rm=TRUE)
resN[, 2:6 ] = round(resN[, 2:6 ] / resN$Total * 100, digits=2)
gt::gt(resN) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
  summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
  footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
  row_group.padding = gt::px(1))
```


```{r}
#| eval: true
#| output: true
#| label: sizefeq-male-survey-cc-nens
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
  
# fig-cap: "Size-frequency of mature male Snow Crab by carapace width (mm) and carapace condition from surveys. NENS."

odir = file.path( SCD, "assessments", year.assessment, "figures", "size.freq", "carapacecondition" )

years = p$year.assessment + c(0:-3) 
cfa = "cfanorth"

fns = paste( "sizefreq", cfa, years, "png", sep="." ) 
fn = file.path( odir, fns ) 
fn = check_file_exists(fn)

include_graphics( fn )

```
Size-frequency of mature male Snow Crab by carapace width (mm) and carapace condition from surveys. NENS.

$~$

SENS:

```{r}
#| eval: true
#| output: true
#| label: table-survey-sens-comm
#| tbl-cap: "Distribution of SENS survey: males less than 95 mm CW by year and shell condition."
resS = dcast( det[ region=="cfasouth" & !is.na(shell), .(N=.N), by=.(fishyr, shell) ], fishyr  ~ shell, value.var="N", fill=0, drop=FALSE, na.rm=TRUE )
names(resS) = c("Year", "CC1", "CC2", "CC3", "CC4", "CC5" )
resS$Total = rowSums( resS[, 2:6 ], na.rm=TRUE)
resS[, 2:6 ] = round(resS[, 2:6 ] / resS$Total * 100, digits=2)
gt::gt(resS) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
  summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
  footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
  row_group.padding = gt::px(1))
```



```{r}
#| eval: true
#| output: true
#| label: sizefeq-male-survey-cc-sens
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2

# fig-cap: "Size-frequency of mature male Snow Crab by carapace width (mm) and carapace condition from surveys. SENS."

odir = file.path( SCD, "assessments", year.assessment, "figures", "size.freq", "carapacecondition" )

years = p$year.assessment + c(0:-3) 

cfa = "cfansouth"

fns = paste( "sizefreq", cfa, years, "png", sep="." ) 
fn = file.path( odir, fns ) 
fn = check_file_exists(fn)

include_graphics( fn )

```
Size-frequency of mature male Snow Crab by carapace width (mm) and carapace condition from surveys. SENS.

$~$


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
gt::gt(resX) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
  summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
  footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
  row_group.padding = gt::px(1))
```

 
```{r}
#| eval: true
#| output: true
#| label: sizefeq-male-survey-cc-4x
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2

# fig-cap: "Size-frequency of mature male Snow Crab by carapace width (mm) and carapace condition from surveys. 4X."

odir = file.path( SCD, "assessments", year.assessment, "figures", "size.freq", "carapacecondition" )

years = p$year.assessment + c(0:-3) 

cfa = "cfa4x"

fns = paste( "sizefreq", cfa, years, "png", sep="." ) 
fn = file.path( odir, fns ) 
fn = check_file_exists(fn)

include_graphics( fn )

```

Size-frequency of mature male Snow Crab by carapace width (mm) and carapace condition from surveys. 4X.

$~$

 

## Size-frequency distributions of snow crab cw, by sex and maturity


```{r}
#| label: create-size-frequency-sex-maturity
#| eval: true
#| output: true

# take subset in years
years = as.character( c(-9:0) + year.assessment )
# years = as.character(2004:2013)
# years = as.character(1996:2003)


regions=c("cfanorth", "cfasouth", "cfa4x")
# outdir=file.path( p$annual.results, "figures", "size.freq", "survey_1996_2003" )
# outdir=file.path( p$annual.results, "figures", "size.freq", "survey_2004_2013" )
outdir=file.path( p$annual.results, "figures", "size.freq", "survey" )

 
M = size_distributions(p=p, toget="crude", xrange=xrange, dx=dx, Y=years )

# NOTE :: these produce png files (instead of pdfs) change as required.
# den=arithmetic mean density, denl = geometric mean density  
quietly( 
  plot_histogram_carapace_width( M=M, years=years, regions=regions, plot_sex="female", yvar="denl", 
  outdir=outdir, cols = c("slategray", "gray95" ) )
)

quietly( 
  plot_histogram_carapace_width( M=M, years=years, regions=regions, plot_sex="male", yvar="denl", 
  outdir=outdir, cols = c("slategray", "gray95" ) )
)

  if (0) {
    # deprecated methods:
    histograms.size.maturity.update( outdir=file.path( p$annual.results, "figures", "size.freq", "survey"),  redo.data=T )
    histograms.size.maturity.single.area( outdir=file.path( p$annual.results, "figures", "size.freq", "survey"),  area='cfa4x',redo.data=T ) #area = cfanorth, cfasouth of cfa4x

    if (oneoff_2022){
      histograms.size.maturity_oneoff(p=p)
      figure.timeseries.survey_oneoff(p=p,
        outdir=file.path(p$annual.results, "timeseries", "survey", "oneoff"), 
        vlab="R0.mass", variables="totmass", plotyears=ts_years) # just R0 to see
    }
  }

```
 
 
```{r}
#| label: sizefeq-male
#| eval: true
#| output: true
#| fig-cap: "Size-frequency (areal density; no/km$^2$) histograms by carapace width of male Snow Crab. The vertical line represents the legal size (95 mm). Immature animals are shown with light coloured bars, mature with dark."
#| fig-dpi: 144
#| fig-height: 4 

fn = file.path( p$annual.results, "figures", "size.freq", "survey",  "male.denl.png" )
include_graphics( fn )
```
 
```{r}
#| label: sizefeq-female
#| eval: true
#| output: true
#| fig-cap: "Size-frequency (areal density; no/km$^2$) histograms by carapace width of female Snow Crab. Immature animals are shown with light coloured bars, mature with dark."
#| fig-dpi: 144
#| fig-height: 4 

fn = file.path( p$annual.results, "figures", "size.freq", "survey",  "female.denl.png" )
include_graphics( fn )
```

 


## Timeseries and maps of survey variables of interest


### Bottom temperature: trawl survey


```{r}
#| label: figures-temperature-bottom-ts
#| eval: true
#| output: true
#| fig-cap: "Annual variations in bottom temperature observed during the Snow Crab survey. The horizontal (black) line indicates the long-term, median temperature within each subarea. Error bars represent standard errors."
#| fig-dpi: 144
#| fig-height: 4 

ts_outdir = file.path( p$annual.results, "timeseries", "survey")

fn = file.path( ts_outdir, paste("t", "png", sep=".") )
include_graphics( fn )

```

```{r}
#| label: figures-temperature-bottom-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
 

map_outdir = file.path( p$project.outputdir, "maps", "survey", "snowcrab", "annual" )
map_years  = p$year.assessment + c(0:-3)
  
fn = check_file_exists( file.path( 
  map_outdir, "t", paste( "t", map_years, "png", sep="." )  
) )

include_graphics( fn )
```

Snow Crab survey bottom temperatures ($~^\circ$C). Note, there is no data in 2020.

$~$


### Sex ratios


```{r}
#| label: figures-sexratio-mat-ts
#| eval: true
#| output: true
#| fig-cap: "The crude, unadjusted geometric mean of sex ratios (proportion female) of mature Snow Crab. Error bars represent 95\\% Confidence Intervals. Note the absence of data in 2020. Prior to 2004, surveys were conducted in the Spring."
#| fig-dpi: 144
#| fig-height: 4 

ts_outdir = file.path( p$annual.results, "timeseries", "survey")

fn = file.path( ts_outdir, paste("sexratio.mat", "png", sep=".") )
include_graphics( fn )

```

```{r}
#| label: figures-sexratio-mat-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2

# fig-cap: "Snow Crab survey sex ratios (proportion female) of mature Snow Crab. Note, there is no data in 2020."

map_outdir = file.path( p$project.outputdir, "maps", "survey", "snowcrab", "annual" )
map_years  = p$year.assessment + c(0:-3)
  
fn = check_file_exists( file.path( 
  map_outdir, "sexratio.mat", paste( "sexratio.mat", map_years, "png", sep="." )  
) )

include_graphics( fn )
```

Snow Crab survey sex ratios (proportion female) of mature Snow Crab. Note, there is no data in 2020.

$~$


### Mature female
 

```{r}
#| label: figures-totno-female-mat-ts
#| eval: true
#| output: true
#| fig-cap: "The crude, unadjusted geometric mean of fature female density log$_{10}$(no/km$^2$) from the Snow Crab survey. Error bars represent 95\\% Confidence Intervals. Note the absence of data in 2020. Prior to 2004, surveys were conducted in the Spring."
#| fig-dpi: 144
#| fig-height: 4 

ts_outdir = file.path( p$annual.results, "timeseries", "survey")

fn = file.path( ts_outdir, paste("totno.female.mat", "png", sep=".") )
include_graphics( fn )

```

```{r}
#| label: figures-totno-female-mat-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2

# fig-cap: "Mature female density log$_{10}$(no/km$^2$) from the Snow Crab survey."

map_outdir = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual" )
map_years  = p$year.assessment + c(0:-3)
  
fn = check_file_exists( file.path( 
  map_outdir, "totno.female.mat", paste( "totno.female.mat", map_years, "png", sep="." )  
) )

include_graphics( fn )
```

Mature female density log$_{10}$(no/km$^2$) from the Snow Crab survey.

$~$


### Fishable biomass 

```{r}
#| label: figures-R0-ts
#| eval: true
#| output: true
#| fig-cap: "The crude, unadjusted geometric mean fishable biomass density log~10(t/km$^2$) from the Snow Crab survey. Error bars represent 95\\% Confidence Intervals. Note the absence of data in 2020. Prior to 2004, surveys were conducted in the Spring."
#| fig-dpi: 144
#| fig-height: 4 

ts_outdir = file.path( p$annual.results, "timeseries", "survey")

fn = file.path( ts_outdir, paste("R0.mass", "png", sep=".") )
include_graphics( fn )

```

```{r}
#| label: figures-R0-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2

# fig-cap: "Snow Crab survey fishable component biomass density log~10(t/km$^2$). Note, there is no data in 2020."
#| 
map_outdir = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual" )
map_years  = p$year.assessment + c(0:-3)
  
fn = check_file_exists( file.path( 
  map_outdir, "R0.mass", paste( "R0.mass", map_years, "png", sep="." )  
) )

include_graphics( fn )
```

Snow Crab survey fishable component biomass density log~10(t/km$^2$). Note, there is no data in 2020.

$~$
   


### Fishable mean size 

```{r}
#| label: figures-cw-male-mat-ts
#| eval: true
#| output: true
#| fig-cap: "Mean size of mature male Snow Crab log10(CW; mm) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

ts_outdir = file.path( p$annual.results, "timeseries", "survey")

fn = file.path( ts_outdir, paste("cw.male.mat.mean", "png", sep=".") )
include_graphics( fn )

```

```{r}
#| label: figures-cw-male-mat-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2

# fig-cap: "Snow Crab survey fishable component mean carapace width; log10(CW; mm). Note, there is no data in 2020."

map_outdir = file.path( p$project.outputdir, "maps", "survey", "snowcrab", "annual" )
map_years  = p$year.assessment + c(0:-3)
  
fn = check_file_exists( file.path( 
  map_outdir, "cw.male.mat.mean", paste( "cw.male.mat.mean", map_years, "png", sep="." )  
) )

include_graphics( fn )
```

Snow Crab survey fishable component mean carapace width; log10(CW; mm). Note, there is no data in 2020.

$~$
   

 
### Predators 

The main predators, based on literature and stomach content analysis, are: 

cod, haddock, halibut, plaice, wolfish, thornyskate, smoothskate, winterskate.

#### Atlantic cod

  
```{r}
#| label: figures-atlcod-ts
#| eval: true
#| output: true
#| fig-cap: "Mean density of Atlantic cod log10(no) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

ts_outdir = file.path( p$annual.results, "timeseries", "survey")
species_predator = 10

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```

```{r}
#| label: figures-atlcod-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2

map_years  = p$year.assessment + c(0:-3)
  
species_predator = 10
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```

Atlantic cod, mean density; log10(CW; mm). Note, there is no data in 2020.

$~$
   

#### Haddock

  
```{r}
#| label: figures-haddock-ts
#| eval: true
#| output: true
#| fig-cap: "Mean density of Haddock log10(no) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

ts_outdir = file.path( p$annual.results, "timeseries", "survey")
species_predator = 11

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```

```{r}
#| label: figures-haddock-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2

map_years  = p$year.assessment + c(0:-3)
  
species_predator = 11
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```

Haddock, mean density; log10(CW; mm). Note, there is no data in 2020.

$~$
     
     
 

#### Halibut

  
```{r}
#| label: figures-halibut-ts
#| eval: true
#| output: true
#| fig-cap: "Mean density of Halibut log10(no) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

ts_outdir = file.path( p$annual.results, "timeseries", "survey")
species_predator = 30

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```

```{r}
#| label: figures-halibut-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2

map_years  = p$year.assessment + c(0:-3)
  
species_predator = 30
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```

Halibut, mean density; log10(CW; mm). Note, there is no data in 2020.

$~$


 

#### American plaice

  
```{r}
#| label: figures-amerplaice-ts
#| eval: true
#| output: true
#| fig-cap: "Mean density of American plaice log10(no) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

ts_outdir = file.path( p$annual.results, "timeseries", "survey")
species_predator = 40

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```

```{r}
#| label: figures-amerplaice-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2

map_years  = p$year.assessment + c(0:-3)
  
species_predator = 40
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```

American plaice, mean density; log10(CW; mm). Note, there is no data in 2020.

$~$

#### Striped Atlantic wolffish

  
```{r}
#| label: figures-stripatlwolffish-ts
#| eval: true
#| output: true
#| fig-cap: "Mean density of Striped Atlantic wolffish log10(no) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

ts_outdir = file.path( p$annual.results, "timeseries", "survey")
species_predator = 50

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```

```{r}
#| label: figures-stripatlwolffish-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2

map_years  = p$year.assessment + c(0:-3)
  
species_predator = 50
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```

Striped Atlantic wolffish, mean density; log10(CW; mm). Note, there is no data in 2020.

$~$
 
#### Thorny skate

  
```{r}
#| label: figures-thornyskate-ts
#| eval: true
#| output: true
#| fig-cap: "Mean density of Thorny skate log10(no) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

ts_outdir = file.path( p$annual.results, "timeseries", "survey")
species_predator = 201

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```

```{r}
#| label: figures-thornyskate-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2

map_years  = p$year.assessment + c(0:-3)
  
species_predator = 201
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```

Thorny skate, mean density; log10(CW; mm). Note, there is no data in 2020.

$~$
 
#### Smooth skate

  
```{r}
#| label: figures-smoothskate-ts
#| eval: true
#| output: true
#| fig-cap: "Mean density of Smooth skate log10(no) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

ts_outdir = file.path( p$annual.results, "timeseries", "survey")
species_predator = 202

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```

```{r}
#| label: figures-smoothskate-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2

map_years  = p$year.assessment + c(0:-3)
  
species_predator = 202
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```

Smooth skate, mean density; log10(CW; mm). Note, there is no data in 2020.

$~$
 
#### Winter skate

  
```{r}
#| label: figures-winterskate-ts
#| eval: true
#| output: true
#| fig-cap: "Mean density of Winter skate log10(no) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

ts_outdir = file.path( p$annual.results, "timeseries", "survey")
species_predator = 204

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```

```{r}
#| label: figures-winterskate-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2

map_years  = p$year.assessment + c(0:-3)
  
species_predator = 204
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```

Winter skate, mean density; log10(CW; mm). Note, there is no data in 2020.

$~$
  


### Competitors

The main potential predators, based on literature and overlpping distributions are: 

northernshrimp, jonahcrab, lessertoadcrab.
   

#### Northern shrimp

  
```{r}
#| label: figures-northernshrimp-ts
#| eval: true
#| output: true
#| fig-cap: "Mean density of Northern shrimp log10(no) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

ts_outdir = file.path( p$annual.results, "timeseries", "survey")
species_predator = 2211

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```

```{r}
#| label: figures-northernshrimp-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2

map_years  = p$year.assessment + c(0:-3)
  
species_predator = 2211
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```

Northern shrimp, mean density; log10(CW; mm). Note, there is no data in 2020.

$~$
  

#### Jonah crab

Not exactly a competitor. Similar habitat except warmer areas so more an indicator of bottom temperatures.
  
```{r}
#| label: figures-jonahcrab-ts
#| eval: true
#| output: true
#| fig-cap: "Mean density of Jonah crab log10(no) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

ts_outdir = file.path( p$annual.results, "timeseries", "survey")
species_predator = 2511

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```

```{r}
#| label: figures-jonahcrab-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2

map_years  = p$year.assessment + c(0:-3)
  
species_predator = 2511
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```

Jonah crab, mean density; log10(CW; mm). Note, there is no data in 2020.

$~$
  

#### Arctic Lyre crab (Lesser toad crab)
 
Slightly more shallow environments than snow crab.


```{r}
#| label: figures-lyrecrab-ts
#| eval: true
#| output: true
#| fig-cap: "Mean density of Arctic Lyre crab log10(no) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

ts_outdir = file.path( p$annual.results, "timeseries", "survey")
species_predator = 2521

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```

```{r}
#| label: figures-lyrecrab-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2

map_years  = p$year.assessment + c(0:-3)
  
species_predator = 2521
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```

Arctic Lyre crab, mean density; log10(CW; mm). Note, there is no data in 2020.

$~$

   
 

<!--

#### Northern stone crab  2524
 
  
$~$

deprecated = TRUE
if (deprecated) {
  # the following are deprecated methods as of 2024 (JC), here only for reference .. most functionality is nowin the Quarto/Rmarkdown above. 
  

  # Tables obtained after completion of data assimilation and processing up to the end of "01.snowcrab.r"
 
  year.assessment = 2024

  p = bio.snowcrab::load.environment( year.assessment=year.assessment )
 

  require(gridExtra)
  library("xtable")
  library("R2HTML")

  odb0 = observer.db("odb")
  regions = c("cfanorth", "cfasouth", "cfa4x")
  nregions = length(regions)

  #------------------------------------------------
  #Fisheries statistics per region
  tabledir = file.path(project.datadirectory("bio.snowcrab"), "data", "fisheries")
  outtabledir= file.path(project.datadirectory("bio.snowcrab"), "assessments", p$year.assessment, "tables", "logbook")
  if(!dir.exists(tabledir)) dir.create(tabledir, recursive =T)
  if(!dir.exists(outtabledir)) dir.create(outtabledir, recursive =T)

  setwd(tabledir)

  NFS <- xtable(read.csv("NENS_FisherySummary.csv"))
  SFS <- xtable(read.csv("SENS_FisherySummary.csv"))
  Fx <- xtable(read.csv("4x_FisherySummary.csv"))

  setwd(outtabledir)
  print.xtable(NFS, type="latex", file="NENS_FisherySummary.tex")
  print.xtable(NFS, type="html", file="NENS_FisherySummary.html")

  print.xtable(SFS, type="latex", file="SENS_FisherySummary.tex")
  print.xtable(SFS, type="html", file="SENS_FisherySummary.html")

  print.xtable(Fx, type="latex", file="4x_FisherySummary.tex")
  print.xtable(Fx, type="html", file="4x_FisherySummary.html")

  #regions = c("cfaall")
  #regions = c("cfanorth", "cfasouth", "cfa4x")
  regions = c("cfanorth", "cfa23", "cfa24", "cfa4x")
  l = NULL
  for (r in regions) {
    res = get.fishery.stats.by.region( Reg=r) #need to add the TACs per year and number of licences
    #round the landings to ton
    #round CPUE to no decimal places
    #round the effort to per x1000 trap hauls
    print(r)
    print(res)
  }


  # ----------------------------------------
  #  Carapace condition from observed data  < 95mm CW

  outtabledir= file.path(project.datadirectory("bio.snowcrab"), "assessments", p$year.assessment, "tables", "observer")
  dir.create(outtabledir, recursive=TRUE)

  odb = odb0
  odb = odb[ which( odb$cw < 95 & odb$prodcd_id=="0" ) ,]
  regions = c("cfanorth", "cfasouth", "cfa4x")
  nregions = length(regions)
  years = sort( unique( odb$fishyr ) )

  res = NULL
  for (r in p$regions) {
    for (y in years) {
      out = proportion.cc (odb, region=r, year=y)
      res = rbind( res, cbind( r, y, t(out)) )
    }}

  cnames = c("region", "fishyr", c(1:5), "ntot")
  colnames(res) = cnames
  print(res)
  res = as.data.frame(res)
  res[is.na(res)] <- NA
  ct <- c("CC1", "CC2", "CC3", "CC4", "CC5", "Total")

  setwd(outtabledir)

  Rn= res[res$region=="cfanorth", 3:8]
  print(Rn)
  rownames(Rn) = years
  colnames(Rn) = ct
  print.xtable(Rn, type="latex", file="table.CC.small.north.obs.tex")
  HTML(Rn, file="table.CC.Small.north.obs.html")

  Rs= res[res$region=="cfasouth", 3:8]
  rownames(Rs) = years
  colnames(Rs) = ct
  print.xtable(Rs, type="latex", file="table.CC.small.south.obs.tex")
  HTML(Rs, file="table.CC.small.south.obs.html")

  Rx= res[res$region=="cfa4x", 3:8]
  rownames(Rx) = years
  colnames(Rx) = ct
  print.xtable(Rs, type="latex", file="table.CC.small.4x.obs.tex")
  HTML(Rx, file="table.CC.small.4x.obs.html")


  # ----------------------------------------
  #  Carapace condition from observed data >=95mm CW
  odb = odb0
  odb = odb[ which( odb$cw >= 95 & odb$cw < 170 & odb$prodcd_id=="0" ) ,]  # commerical sized crab only
  years = sort( unique( odb$fishyr ) )

  # get proportion by cc
  regions = c("cfanorth", "cfasouth", "cfa4x")
  years = sort( unique( odb$fishyr ) )

  res = NULL
  for (r in regions) {
    for (y in years) {
      out = proportion.cc (odb, region=r, year=y)
      res = rbind( res, cbind( r, y, t(out)) )
    }}

  cnames = c("region", "fishyr", c(1:5), "ntot")
  colnames(res) = cnames
  print(res)
  res = as.data.frame(res)
  res[is.na(res)] <- NA
  #  for (i in cnames[-1]) res[,i] = as.numeric(as.character((res[,i])))
  setwd(outtabledir)

  ct <- c("CC1", "CC2", "CC3", "CC4", "CC5")
  Rn = res[res$region=="cfanorth", 3:7]
  #Rn = as.matrix( res[ which(res$region=="cfanorth") , as.character(c(1:5)) ] )
  rownames(Rn) = years
  colnames(Rn) = ct
  print.xtable(Rn, type="latex", file="table.CC.large.north.obs.tex")
  HTML(Rn, file="table.CC.large.north.obs.html")

  #Rs = as.matrix( res[ which(res$region=="cfasouth") , as.character(c(1:5)) ] )
  Rs = res[res$region=="cfasouth", 3:7]
  rownames(Rs) = years
  colnames(Rs) = ct
  print.xtable(Rs, type="latex", file="table.CC.large.south.obs.tex")
  HTML(Rs, file="table.CC.large.south.obs.html")

  #Rx = as.matrix( res[ which(res$region=="cfa4x") , as.character(c(1:5)) ] )
  Rx = res[res$region=="cfa4x", 3:7]
  rownames(Rx) = years
  colnames(Rx) = ct
  print.xtable(Rx, type="latex", file="table.CC.large.4x.obs.tex")
  HTML(Rx, file="table.CC.large.4x.obs.html")

  # ----------------------------------------
  #  Percent soft from observed data

  odb = odb0
  odb = odb[ which( odb$cw > 95 & odb$cw < 170 & odb$prodcd_id=="0" ) ,]  # commercial crab
  years = sort( unique( odb$fishyr ) )


  res = NULL
  for (r in p$regions) {
    for (y in years) {
      out = proportion.soft (odb, region=r, year=y)
      res = rbind( res, cbind( r, y, t(out)) )
    }}

  cnames = c("region", "fishyr", "pr.soft", "nsoft", "ntot")
  colnames(res) = cnames
  print(res)
  res = as.data.frame(res)

  for (i in cnames[-1]) res[,i] = as.numeric(as.character((res[,i])))

  Rn = as.matrix( res[ which(res$region=="cfanorth") , ] )
  rownames(Rn) = years
  HTML(Rn, file="table.proportion.soft.north.obs.html")

  Rs = as.matrix( res[ which(res$region=="cfasouth" ),   ] )
  rownames(Rs) = years
  HTML(Rs, file="table.proportion.soft.south.obs.html")

  Rx = as.matrix( res[ which(res$region=="cfa4x") , ] )
  rownames(Rx) = years
  HTML(Rx, file="table.proportion.soft.4x.obs.html")



  # instars of interest: 11 and 12

  # growth increment (assumming average weight in the midpoint of each increment)
  growth.11.to.12 =  predict.mass.g.from.CW.mm( mean(CW.interval.male(12)) ) - predict.mass.g.from.CW.mm (mean(CW.interval.male(11)) )

  # = 419 g
  #  12to13 = ~450



  # Table of proportion discarded
  odb = observer.db("odb")
  regions = c("cfanorth", "cfasouth", "cfa4x")
  years = sort( unique( odb$fishyr ) )
  out = NULL
  for (r in regions) {
    for (y in years) {
      res = proportion.legal (odb, region=r, year=y)
      out = rbind(out, cbind( r, y, res[1], res[2], res[3] ) )
    } }
  out

  HTML(out, file="table.proportion.discarded.html")


  # ---------------------------------------- USED
  #  Carapace condition from trawl data  >= 95mm CW  ... not kriged .. simple proportions

  det0 = snowcrab.db( p=p, DS="det.georeferenced" )
  det0$fishyr = det0$yr  ## the counting routine expectes this variable

  det = det0[ which( det0$cw >= 95 ) ,]  # commerical sized crab only
  years = sort( unique( det$yr ) )

  res = NULL
  for (r in p$regions) {
    for (y in years) {
      out = proportion.cc (det, region=r, year=y)
      res = rbind( res, cbind( r, y, t(out)) )
    }}

  cnames = c("region", "fishyr", c(1:5), "ntot")
  colnames(res) = cnames
  print(res)
  res = as.data.frame(res)

  for (i in cnames[-1]) res[,i] = as.numeric(as.character((res[,i])))
  (res)

  HTML(res, file="table.CC.large.survey.html")

  # ------------------
  # counts of stations in each area

  # check towquality .. this should always == 1
  set = snowcrab.db(p=p, DS="set.clean")
  if (length( unique( set$towquality) ) != 1 ) print("error -- not good tows")

  out = data.frame(yr=sort( unique(set$yr )) )
  for (reg in c("cfaall", "cfanorth", "cfasouth","cfa4x"  ) ) {
    d = polygon_inside(set[,c("lon","lat")], reg)
    e = as.data.frame( xtabs(~yr, data=set[d,])  )
    names(e) = c("yr", reg)
    e$yr = as.numeric(as.character(e$yr) )
    out = merge(out, e, by="yr", all=T)
  }
  print(out)
 
  HTML(out, file="table.tow_counts_survey.html")



  # % mat calculations: deprecated  ... size analysis is now in 01_snowcrab.R

  # loc = file.path(sc.R, "size.data")
  # dir.create(path=loc, recursive=T, showWarnings=F)
  # outfilename = paste( c("mi", "mm", "fi", "fm"), "rdata", sep=".")
  # outfile = file.path(loc, paste(outfilename))
  # for (f in  outfile) load(f)


  # f.i = f.imm[which( rownames(f.imm)%in% sids ) ,]
  # f.i.means = apply(X=f.i, MARGIN=2, FUN=mean)
  # f.m = f.mat[which( rownames(f.mat)%in% sids ) ,]
  # f.m.means = apply(X=f.m, MARGIN=2, FUN=mean)

  # toplot = rbind(f.m.means, f.i.means)

  # ii = as.data.frame(t(toplot))
  # ii$cw = as.numeric(rownames(ii))
  # ii$pmat = ii[,1]/ (ii[,1]+ii[,2]) * 100

  # plot(ii$cw, ii$pmat)
  # abline(h=50)

  # str(ii)
  #`data.frame':   70 obs. of  4 variables:
  # $ f.m.means: num  0 0 0 0 0 ...
  # $ f.i.means: num   2.80  6.19 20.05 24.29 74.11 ...
  # $ cw       : num  12 14 16 18 20 22 24 26 28 30 ...
  # $ pmat     : num  0 0 0 0 0 ...
 

  # ----------------------------------------   NOT USED ____________
  #  Carapace condition from trawl data  < 95mm CW  ... not kriged .. simple proportions

  det0 = snowcrab.db( p=p, DS="det.georeferenced" )
  det0$fishyr = det0$yr  ## the counting routine expectes this variable

  det = det0[ which( det0$cw < 95 ) ,]  # commerical sized crab only
  years = sort( unique( det$yr ) )

  res = NULL
  for (r in p$regions) {
    for (y in years) {
      out = proportion.cc (det, region=r, year=y)
      res = rbind( res, cbind( r, y, t(out)) )
    }}

  cnames = c("region", "fishyr", c(1:5), "ntot")
  colnames(res) = cnames
  print(res)
  res = as.data.frame(res)

  for (i in cnames[-1]) res[,i] = as.numeric(as.character((res[,i])))
  (res)



# deprecated ?

  # % mat calculations:
  #as above but

  loc = file.path(sc.R, "size.data")
  dir.create(path=loc, recursive=T, showWarnings=F)
  outfilename = paste( c("mi", "mm", "fi", "fm"), "rdata", sep=".")
  outfile = file.path(loc, paste(outfilename))
  for (f in  outfile) load(f)


  f.i = f.imm[which( rownames(f.imm)%in% sids ) ,]
  f.i.means = apply(X=f.i, MARGIN=2, FUN=mean)
  f.m = f.mat[which( rownames(f.mat)%in% sids ) ,]
  f.m.means = apply(X=f.m, MARGIN=2, FUN=mean)

  toplot = rbind(f.m.means, f.i.means)

  ii = as.data.frame(t(toplot))
  ii$cw = as.numeric(rownames(ii))
  ii$pmat = ii[,1]/ (ii[,1]+ii[,2]) * 100

  plot(ii$cw, ii$pmat)
  abline(h=50)

  str(ii)
  #`data.frame':   70 obs. of  4 variables:
  # $ f.m.means: num  0 0 0 0 0 ...
  # $ f.i.means: num   2.80  6.19 20.05 24.29 74.11 ...
  # $ cw       : num  12 14 16 18 20 22 24 26 28 30 ...
  # $ pmat     : num  0 0 0 0 0 ...







  # ----------------------------------------   NOT USED ____________
  #  Carapace condition from trawl data  < 95mm CW  ... not kriged .. simple proportions

  det0 = snowcrab.db( p=p, DS="det.georeferenced" )
  det0$fishyr = det0$yr  ## the counting routine expectes this variable

  det = det0[ which( det0$cw < 95 ) ,]  # commerical sized crab only
  years = sort( unique( det$yr ) )

  res = NULL
  for (r in p$regions) {
    for (y in years) {
      out = proportion.cc (det, region=r, year=y)
      res = rbind( res, cbind( r, y, t(out)) )
    }}

  cnames = c("region", "fishyr", c(1:5), "ntot")
  colnames(res) = cnames
  print(res)
  res = as.data.frame(res)

  for (i in cnames[-1]) res[,i] = as.numeric(as.character((res[,i])))
  (res)




  ########################################################
  ########################  Retired figures ###################



  # ------------------------------------------
  # Map: Scotian Shelf with CFA lines and labels  .. using gmt
  # this is the basemap from map.r which is then post-labelled in sodipodi
  #  p$outdir = file.path(p$annual.results,"figures")
  #  p$outfile.basename = file.path(p$outdir, "map.CFAs")
  #  map.basemap.with.cfa.lines( p, conversions=c("ps2png")  )


  # ------------------------------------------
  # Habitat usage comparisons (bivariate) ... requires the full "set.rdata" database and "logbook.dZ.rdata" database
  # habitat.usage( usevar="totno.all", covariate="depth", outdir = file.path(p$annual.results, "habitat.templates") )

  # ------------------------------------------
  # Habitat usage comparisons (bivariate) ... requires the full "set.rdata" database and "logbook.dZ.rdata" database
  #habitat.usage( usevar="totno.all", covariate="temperature", outdir = file.path(p$annual.results, "habitat.templates") )

  # ------------------------------------------
  # Habitat usage comparisons (bivariate) ... requires the full "set.rdata" database and "logbook.dZ.rdata" database
  # habitat.usage( usevar="totno.all", covariate="bottom.slope", outdir = file.path(p$annual.results, "habitat.templates") )

  # ------------------------------------------
  # Habitat usage comparisons (bivariate) ... requires the full "set.rdata" database and "logbook.dZ.rdata" database
  #habitat.usage( usevar="totno.all", covariate="bottom.curvature", outdir = file.path(p$annual.results, "habitat.templates") )

  # ------------------------------------------
  # Habitat usage comparisons (bivariate) ... requires the full "set.rdata" database and "logbook.dZ.rdata" database
  #habitat.usage( usevar="totno.all", covariate="substrate", outdir = file.path(p$annual.results, "habitat.templates") )

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
