---
title: "Snow Crab Assessment"
subtitle: "Scotian Shelf Ecosystem"
author:
  - name: 
      given: Snow Crab Unit, DFO Science 
    #  family: Choi
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
toc: false
toc-depth: 4
number-sections: false
highlight-style: pygments
# bibliography: media/references.bib  
# csl: media/canadian-journal-of-fisheries-and-aquatic-sciences.csl  # see https://www.zotero.org/styles for more
# citation: 
#  container-title: https://github.com/jae0/bio.snowcrab/
#  doi: NA
keywords: 
  - snow crab fishery assessment
  - areal-unit modelling of numerical abundance, habitat, mean weight
  - convolution of posterior-simulations of biomass   
abstract: |
  Snow crab modelling of numerical abundance, mean weight and habitat using conditional autoregressive models. Biomass is estimated from these components via multiplication of the posterior draws.   
bibliography: references.bib  
# csl: media/canadian-journal-of-fisheries-and-aquatic-sciences.csl  # see https://www.zotero.org/styles for more
# license: "CC BY"
copyright: 
  holder: Jae S. Choi
  year: 2024
citation: 
  container-title: https://github.com/jae0/bio.snowcrab/
  doi: NA
funding: "The snow crab scientific survey was funded by the snow crab fishers of Maritimes Region of Atlantic Canada."
editor:
  render-on-save: false
format:
  docx:
    toc: true
    toc_depth: 2 
    highlight-style: github
    tbl-cap-location: top
    reference-doc: media/csas_template.docx
    # reference-doc: media/RES2024-eng.docx
  revealjs:
    theme: moon    
    self-contained: true
    scrollable: true
    smaller: true
    dpi: 144
  html: 
    code-fold: true
    code-overflow: wrap
    html-math-method: katex
    self-contained: true
    embed-resources: true
  pdf:
    pdf-engine: lualatex
  beamer:
    theme: "metropolis"
    colortheme: "seagull"
    fonttheme: "professionalfonts"
    fig_caption: yes
    latex_engine: lualatex 
    keep_tex: true 
# classoption: 
#  - aspectratio=169 #16:9 wide
#  - t  # top align
# header-includes: 
#   - \usepackage{graphicx}
 #  - \usepackage[font={scriptsize}, labelfont={bf}]{caption}
 #  - \usepackage{longtable}
 # - \usepackage{booktabs}
 # - \usepackage{caption}
 # - \usepackage{float}
 # - \usepackage{multicol}
  # - \usepackage{subfig}
  # - \newcommand{\btiny}{\begin{tiny}}
  # - \newcommand{\etiny}{\end{tiny}}
params:
  year_assessment: 2024
  year_start: 1999
  data_loc:  "~/bio.data/bio.snowcrab"
  sens: 1
  debugging: FALSE
  model_variation: logistic_discrete_historical
--- 


<!-- Preamble
 
  make quarto FN=snowcrab_presentation_assessment.md DOCTYPE=revealjs  PARAMS="-P year_assessment:2024"  --directory=~/bio/bio.snowcrab/inst/markdown

  make quarto FN=snowcrab_presentation_assessment.md DOCTYPE=html   PARAMS="-P year_assessment:2024"  --directory=~/bio/bio.snowcrab/inst/markdown
      
  make quarto FN=snowcrab_presentation_assessment.md DOCTYPE=beamer  PARAMS="-P year_assessment:2024"  --directory=~/bio/bio.snowcrab/inst/markdown 


-->




<!-- Set up R-environment -->

 
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

  require(flextable)
  require(gt)  # table formatting
  require(ggplot2)
  require(data.table) # for speed
  require(lubridate)
  require(stringr) 
  
  require(janitor)
  
  require(aegis)  # basic helper tools
  require(bio.taxonomy)  # handle species codes
  require(bio.snowcrab)

  require(spsUtil)
  quietly = spsUtil::quiet

  require(ggplot2)
  require(MBA)
  require(aegis)  # basic helper tools

  loadfunctions( "aegis")
  loadfunctions( "bio.snowcrab")  # in case of local edits

knit_print.flextable <- function (x, ...) {
    is_bookdown <- flextable:::is_in_bookdown()
    is_quarto <- flextable:::is_in_quarto()
    x <- flextable:::knitr_update_properties(x, bookdown = is_bookdown, quarto = is_quarto)
    if (is.null(knitr::pandoc_to())) {
        str <- flextable:::to_html(x, type = "table")
        str <- knitr:::asis_output(str)
    }
    else if (!is.null(getOption("xaringan.page_number.offset"))) {
        str <- knit_to_html:::knit_to_html(x, bookdown = FALSE, quarto = FALSE)
        str <- knitr:::asis_output(str, meta = html_dependencies_list(x))
    }
    else if (knitr:::is_html_output(excludes = "gfm")) {
        str <- knit_to_html:::knit_to_html(x, bookdown = is_bookdown, quarto = is_quarto)
        str <- knit_to_html:::raw_html(str, meta = html_dependencies_list(x))
    }
    else if (knitr:::is_latex_output()) {
        str <- flextable:::knit_to_latex(x, bookdown = is_bookdown, quarto = is_quarto)
        str <- flextable:::raw_latex(x = str, meta = unname(list_latex_dep(float = TRUE, 
            wrapfig = TRUE)))
    }
    else if (grepl("docx", knitr::opts_knit$get("rmarkdown.pandoc.to"))) {
        if (rmarkdown::pandoc_version() < numeric_version("2")) {
            stop("pandoc version >= 2 required for printing flextable in docx")
        }
        str <- flextable:::knit_to_wml(x, bookdown = is_bookdown, quarto = is_quarto)
        str <- gsub("</w:tbl>","</w:tbl><w:p></w:p>",str, fixed =  TRUE)
        str <- knitr:::asis_output(str)
    }
    else if (grepl("pptx", knitr::opts_knit$get("rmarkdown.pandoc.to"))) {
        if (rmarkdown::pandoc_version() < numeric_version("2.4")) {
            stop("pandoc version >= 2.4 required for printing flextable in pptx")
        }
        str <- flextable:::knit_to_pml(x)
        str <- knitr:::asis_output(str)
    }
    else {
        plot_counter <- getFromNamespace("plot_counter", "knitr")
        in_base_dir <- getFromNamespace("in_base_dir", "knitr")
        tmp <- fig_path("png", number = plot_counter())
        in_base_dir({
            dir.create(dirname(tmp), showWarnings = FALSE, recursive = TRUE)
            save_as_image(x, path = tmp, expand = 0)
        })
        str <- include_graphics(tmp)
    }
    str
}

registerS3method(
  "knit_print", 'flextable', knit_print.flextable, 
  envir = asNamespace("flextable") 
  # important to overwrite {flextable}s knit_print
)



knit_print.gt_tbl <- function (x, ..., inline = FALSE) {
    if (gt:::knitr_is_rtf_output()) {
        x <- gt:::as_rtf(x)
    }
    else if (knitr::is_latex_output()) {
        x <- gt:::as_latex(x)
    }
    else if (gt:::knitr_is_word_output()) {
       str <- gsub("</w:tbl>","</w:tbl><w:p></w:p>", as_word(x), fixed = TRUE)
        x <- knitr::asis_output(paste("\n\n``````{=openxml}", str, 
            "``````\n\n", sep = "\n"))
    }
    else {
        x <- htmltools:::as.tags(x, ...)
    }
    knitr::knit_print(x, ..., inline = FALSE)
}



registerS3method(
  "knit_print", 'gt_tbl', knit_print.gt_tbl, 
  envir = asNamespace("gt") 
  # important to overwrite {gt}s knit_print
)


  year_assessment = params$year_assessment
  
  year_start = params$year_start

  year_previous = year_assessment - 1

  model_variation = params$model_variation
  

data_loc= params$data_loc
media_loc = file.path( params$media_loc, "media" )


#### params and directories

  p = load.environment( year.assessment=year_assessment )  
  
  p$corners = data.frame(plon=c(220, 990), plat=c(4750, 5270) )

  p$mapyears = year_assessment + c(-5:0 )   # default in case not specified
  
  
  years = as.character(1996: year_assessment)
  yrs_observer = year_assessment + c(0:-4)
   
 
  # fishery_model_results = file.path( "/home", "jae", "projects", "dynamical_model", "snowcrab", "outputs" )
  fishery_model_results = file.path( data_loc, "fishery_model" )

 
  # note copied "roadshow figures" temporaily here ... figure creation should be be assimilated TODO
  media_supplementary = file.path( data_loc, "assessments",  year_assessment, "media_supplementary")

  outtabledir = file.path( p$annual.results, "tables" ) 

  lregions = list(region=c("cfanorth", "cfasouth", "cfa4x"))
  reg_labels = c("N-ENS", "S-ENS", "CFA 4X")  # formatted for label

  if (params$sens==2) {
    lregions = list(subarea=c("cfanorth", "cfa23",  "cfa24", "cfa4x"))
    reg_labels = c("CFA 20-22", "CFA 23", "CFA 24", "CFA 4X")  # formatted for label
  }

  vnr = names(lregions)
  regions = unname( unlist(lregions) )
  
  nregions = length(regions)

  
#### fishery 

  FD = fishery_data( regions=lregions)  # mass in tonnes

  fda = FD$summary_annual
  fdm = FD$summary_monthly
  fdb = FD$summary_biweekly

  dt = as.data.frame( fda[ which(fda$yr %in% c(year_assessment - c(0:10))),] )
  dt =  dt[,c(vnr, "yr", "Licenses", "TAC", "landings", "effort", "cpue")] 
  names(dt) = c("Region", "Year", "Licenses", "TAC", "Landings", "Effort", "CPUE") 
  rownames(dt) = NULL

  # observer data

  odb0 = setDT(observer.db("odb"))
  odb0[[vnr]] = NA
  for ( reg in regions) {
    r = polygon_inside(x = odb0, region = aegis.polygons::polygon_internal_code(reg), planar=FALSE)
    odb0[[vnr]][r] = reg
  }

  # bycatch summaries
  BC = list()
  for ( reg in c(regions, "cfaall")) {
    oo = observer.db( DS="bycatch_summary", p=p,  yrs=p$yrs, region=reg  )  # using polygon test
    if (is.null(oo)) { 
      message( "bycatch data is likely broken .. check" ) 
    }
      BC[[reg]] = oo
      BC[[reg]]$bycatch_table_effort[ BC[[reg]]$bycatch_table_effort==0 ] = NA
      BC[[reg]]$bycatch_table_effort[ is.na(BC[[reg]]$bycatch_table_effort) ] = "."
      BC[[reg]]$bycatch_table_catch[ BC[[reg]]$bycatch_table_catch==0 ] = NA
      BC[[reg]]$bycatch_table_catch[ is.na(BC[[reg]]$bycatch_table_catch) ] = "."
    
  }
  

  # predator diet data
  diet_data_dir = file.path( data_loc, "data", "diets" )
  
  # assimilate the CSV data tables:
  # diet = get_feeding_data( diet_data_dir, redo=TRUE )  # if there is a data update
  diet = get_feeding_data( diet_data_dir, redo=FALSE )
  tx = taxa_to_code("snow crab")  
  # matching codes are 
  #  spec    tsn                  tx                   vern tx_index
  #1  528 172379        BENTHODESMUS           BENTHODESMUS     1659
  #2 2522  98427        CHIONOECETES SPIDER QUEEN SNOW UNID      728
  #3 2526  98428 CHIONOECETES OPILIO        SNOW CRAB QUEEN      729
  # 2 and 3 are correct

  snowcrab_predators = diet[ preyspeccd %in% c(2522, 2526), ]  # n=159 oservations out of a total of 58287 observations in db (=0.28% of all data)
  snowcrab_predators$Species = code_to_taxa(snowcrab_predators$spec)$vern
  snowcrab_predators$Predator = factor(snowcrab_predators$Species)
  
  counts = snowcrab_predators[ , .(Frequency=.N), by=.(Species)]
  setorderv(counts, "Frequency", order=-1)
  
  # species composition
  psp = speciescomposition_parameters( yrs=p$yrs, carstm_model_label="default" )
  pca = speciescomposition_db( DS="pca", p=psp )  

  pcadata = as.data.frame( pca$loadings )
  pcadata$vern = stringr::str_to_title( taxonomy.recode( from="spec", to="taxa", tolookup=rownames( pcadata ) )$vern )



#### survey
 

  # recode region to selection above:

  set0 = snowcrab.db(p=p, DS="set.biologicals")
  setDT(set0)
  # check towquality .. this should always == 1
  if (length( unique( set0$towquality) ) != 1 ) print("error -- not good tows")
  set0$region = NA
  for (reg in regions ) {
    d = polygon_inside(set0[,c("lon","lat")], reg)
    set0$region[d] = reg 
  }


  # recode region to selection above:

  det0 = snowcrab.db( p=p, DS="det.georeferenced" )
  setDT(det0)
  det0$fishyr = det0$yr  ## the counting routine expects this variable
  det0$region = NA
  for ( reg in regions) {
    r = polygon_inside(x = det0, region = aegis.polygons::polygon_internal_code(reg), planar=FALSE)
    det0$region[r] = reg
  }
  

#### ecosystem


  # predator diet data
  diet_data_dir = file.path( data_loc, "data", "diets" )
  require(data.table) # for speed
  require(lubridate)
  require(stringr) 
  require(gt)  # table formatting
  library(janitor)
  require(ggplot2)
  require(aegis) # map-related 
  require(bio.taxonomy)  # handle species codes

  # assimilate the CSV data tables:
  # diet = get_feeding_data( diet_data_dir, redo=TRUE )  # if there is a data update
  diet = get_feeding_data( diet_data_dir, redo=FALSE )
  tx = taxa_to_code("snow crab")  
  # matching codes are 
  #  spec    tsn                  tx                   vern tx_index
  #1  528 172379        BENTHODESMUS           BENTHODESMUS     1659
  #2 2522  98427        CHIONOECETES SPIDER QUEEN SNOW UNID      728
  #3 2526  98428 CHIONOECETES OPILIO        SNOW CRAB QUEEN      729
  # 2 and 3 are correct

  snowcrab_predators = diet[ preyspeccd %in% c(2522, 2526), ]  # n=159 oservations out of a total of 58287 observations in db (=0.28% of all data)
  snowcrab_predators$Species = code_to_taxa(snowcrab_predators$spec)$vern
  snowcrab_predators$Predator = factor(snowcrab_predators$Species)

  counts = snowcrab_predators[ , .(Frequency=.N), by=.(Species)]
  setorderv(counts, "Frequency", order=-1)

  # species composition
  psp = speciescomposition_parameters( yrs=p$yrs, carstm_model_label="default" )
  pca = speciescomposition_db( DS="pca", p=psp )  

  pcadata = as.data.frame( pca$loadings )
  pcadata$vern = stringr::str_to_title( taxonomy.recode( from="spec", to="taxa", tolookup=rownames( pcadata ) )$vern )


#### fishery model

  fishery_model_results = file.path( data_loc, "fishery_model" )
    
  fm_loc = file.path( data_loc, 'fishery_model', year_assessment, model_variation )
    

  # as modelled years in fishery model can differ from iput data years, make sure  "years_model" is correct
  p$fishery_model_years = 2000:year_assessment
  sn_env = snowcrab_load_key_results_to_memory( year_assessment, years_model=p$fishery_model_years, return_as_list=TRUE  ) 

  attach(sn_env)
  
  
```
 


## Management

```{r}
#| label: fig-management-areas
#| eval: true
#| echo: false 
#| output: true
#| fig-cap: "The Scotian Shelf (NW Atlantic Ocean; NAFO Div. 4VWX). Shown are isobaths and major bathymetric features. Managed Crab Fishing Areas (CFAs; divided by dashed lines) include: NENS, SENS, 4X. SENS is further subdivided (dotted line) into 23 (NW) and 24 (SE)."
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold 

fn = file.path( media_loc, "snowcrab_cfas.png" )
knitr::include_graphics( fn ) 

```
 



## Data sources: Trawl Survey


```{r}
#| label: fig-survey-locations-map 
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
 



  
## Data sources: Fishery logbooks

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
   


## Fishery performance: N-ENS


```{r}
#| label: tbl-fishery-performance-N-ENS
#| echo: false
#| results: asis
#| eval: true
#| output: true
#| tbl-cap: "Fishery performance statistics. Note: 4X years represent the starting year and currently ongoing."

r=1
  reg = regions[r]
  REG = reg_labels[r]
  # cat(REG, "\n")
  oo = dt[ which(dt$Region==reg), c("Year", "Licenses", "TAC", "Landings", "Effort", "CPUE")] 

  names(oo) = c( "Year", "Licenses", "TAC (t)", "Landings (t)", "Effort (1000 th)", "CPUE (kg/th)" )

  out = gt::gt(oo) |> gt::tab_options(table.font.size = 14, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
  print(out)
  cat("\n\n")

```
 



## Fishery performance: S-ENS


```{r}
#| label: tbl-fishery-performance-S-ENS
#| echo: false
#| results: asis
#| eval: true
#| output: true
#| tbl-cap: "Fishery performance statistics. Note: 4X years represent the starting year and currently ongoing."

r=2
  reg = regions[r]
  REG = reg_labels[r]
  # cat(REG, "\n")
  oo = dt[ which(dt$Region==reg), c("Year", "Licenses", "TAC", "Landings", "Effort", "CPUE")] 

  names(oo) = c( "Year", "Licenses", "TAC (t)", "Landings (t)", "Effort (1000 th)", "CPUE (kg/th)" )

  out = gt::gt(oo) |> gt::tab_options(table.font.size = 14, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
  print(out)
  cat("\n\n")

```
  
## Fishery performance: 4X


```{r}
#| label: tbl-fishery-performance-4X
#| echo: false
#| results: asis
#| eval: true
#| output: true
#| tbl-cap: "Fishery performance statistics. Note: 4X years represent the starting year and currently ongoing."

r=3
  reg = regions[r]
  REG = reg_labels[r]
  # cat(REG, "\n")
  oo = dt[ which(dt$Region==reg), c("Year", "Licenses", "TAC", "Landings", "Effort", "CPUE")] 

  names(oo) = c( "Year", "Licenses", "TAC (t)", "Landings (t)", "Effort (1000 th)", "CPUE (kg/th)" )

  out = gt::gt(oo) |> gt::tab_options(table.font.size = 14, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
  print(out)
  cat("\n\n")


```
  

## Effort

- N-ENS: increased in 2024 relative to the previous year and spatially more dispersed. 
- S-ENS: increased in 2024 relative to the previous year and spatially more dispersed.
- CFA 4X: decreased in 2024-2025 relative to the previous year and spatially contracted.
- Reference: Figures -@fig-effort-timeseries,  -@fig-effort-map
 


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
 


## Landings

- N-ENS: declined marginally and spatially more dispersed 
- S-ENS: declined marginally and spatially more dispersed 
- CFA 4X: only a small fraction of TACs landed and fishery is ongoing
- Reference: Figure -@fig-landings-timeseries, -@fig-landings-map

 
 
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

## Landings 2

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

- N-ENS: decreased to 2019 levels, declines inshore
- S-ENS: decreased to 2017 levels, declines through, especially CFA 24
- CFA 4X: increased to 2021 levels, but only a small fraction of TACs landed in Sambro area; fishery is ongoing
- Reference: -@fig-cpue-timeseries, -@fig-cpue-map
    

## Catch rates 2


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

## Catch rates 3

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
 

## Fishery performance: Summary


- In `r year_assessment`, they were `r l_nens`, `r l_sens` and `r l_4x` t, in N-ENS, S-ENS and 4X (season ongoing), respectively. 

- Relative to `r year_previous`, they represent changes of `r dt_l_nens`%, `r dt_l_sens`% and  `r dt_l_4x`%, respectively. 

- Total Allowable Catches (TACs) for `r year_assessment` were `r tac_nens` t, `r tac_sens` t and `r tac_4x` t in N-ENS, S-ENS and 4X, respectively. 

Non-standardized fishery catch rates in `r year_assessment` were `r c_nens`, `r c_sens`  and `r c_4x` kg/trap haul in N-ENS, S-ENS and 4X, respectively; `r dt_c_nens` %, `r dt_c_sens` % and  `r dt_c_4x` % (season ongoing) relative to the previous year 




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
 

## Discard of snow crab

- N-ENS: increased since 2022, but within historical range; coverage was 1.3% of landings
- S-ENS: increased since 2022, but within historical range; coverage was 2.4% of landings
- CFA 4X: No data, fishery is ongoing 
- Reference: -@fig-discard_maritimes



## Discard of soft-shelf crab

- Generally, higher soft-shell indicates incoming recruitment to the fishery and their handling potential and unnecessary handling/discard mortality.

- Commercial catches of soft-shelled crab were `r cc_soft_nens`% (low sampling), `r cc_soft_sens`% (low sampling) and `r cc_soft_4x`% (no sampling; season ongoing) in N-ENS, S-ENS and 4X, respectively for `r year_assessment`. 

- In `r year_previous`, it was `r cc_soft_nens_p`% (no sampling), `r cc_soft_sens_p`% (low sampling) and `r cc_soft_4x_p`% (no sampling), respectively. 


- N-ENS: increased since 2022, but within historical range
- S-ENS: increased since 2022, but within historical range
- CFA 4X: No data, fishery is ongoing
- Reference: -@fig-observed-softshell-map, -@tbl-observer-softgt95
 
## Discard of soft-shelf crab: N-ENS


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
 




## Discard of non-target species ("Bycatch")

- N-ENS: average of 0.02% of landings; primarily other Crustacea (crab and lobster)  
- S-ENS: average of 0.03% of landings; primarily other Crustacea (crab and lobster)
- CFA 4X: average of 0.87% of landings; primarily other Crustacea (crab and lobster)  
- Reference: -@fig-map-observer-locations, -@tbl-fishery-discard-effort


## Discard of non-target species ("Bycatch"): N-ENS

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

- Low level background infection rate of < 0.1% of the surveyed crab (-@tbl-bcd)
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
 
 
## Bitter crab disease 2


```{r}
#| label: fig-bcd-map
#| eval: true
#| echo: false 
#| output: true
#| fig-cap: "Bitter crab disease observations since 2008"
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
  
fns = file.path( media_loc, c(
  "BCD_map.png"  
) )

knitr::include_graphics( fns ) 

``` 

 

## Survey indices

- Details of survey design and index modelling can be found in (Christie et al. 2024; Choi 2020;  -@fig-survey-locations-map )
- Numerous vessel and captain changes have occured over the years
- 2004: transition from a Spring to Fall survey and mandate transfer from Moncton to BIO
- 2020: due to Covid-19 related to health and safety uncertainties, survey was not conducted 
- 2022: due to mechanical issues with the survey vessel, inshore areas of S-ENS were partially sampled  
- 2024: due to financial issues, areas associated with Marine Protected Areas (St Ann's, Gully) were not sampled
  - N-ENS: 58 stations completed 
  - S-ENS: 282 stations completed
  - CFA 4X: 24 stations completed



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

  

## Recruitment: males

- N-ENS: recruitment gap noted in 2021 has been bridged (year-class centered at 68 mm CW in 2023 in now at 90 mm CW) but subsequent year-classes are weak suggesting high natural mortality (predation?); little internal recruitment is expected for the next 3-4 years
- S-ENS: strong and stable size structure suggestive of a stable size/age structure; recruitment expected in 2025 and elevated soft-shell capture will require care
- CFA 4X: strong years classes (30 mm CW in 2019) were progressing but in 2024 has been lost due to elevated natural morality (predation and high temperatures?); minimal internal recruitment expected 
- Reference: -@fig-sizefeq-male
  


## Recruitment: females

- N-ENS: mature fraction has senesced with minimal recruitment into mature phase; next large year class is at 25 mm CW
- S-ENS: strong and stable mature fraction; steady recruitment into mature phase expected in 2025  
- CFA 4X: mature fraction has been strong for much of the historical record; recruitment into mature phase is expected to be low and the female population to senesce
- Reference: -@fig-sizefeq-female
 
 
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


## Reproduction

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

## Reproduction


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

- N-ENS: sex ratio of mature crab are female have increased marginally; with increases in near shore and in shallower areas
- S-ENS: sex ratios have declined marginally; with declines occurring in Missaine and inshore 24
- CFA 4X: sex ratios have declned but still high; peaks in roseway and further downstream
- Reference: -@fig-sexratio-mat-timeseries, -@fig-sexratio-mat-map

The sex ratios (proportion female) of the mature component is particularly important as locally unbalanced ratios can impact upon encounter rates and ultimately reproductive success. In the ESS, there is generally a lack of females, in contrast to, for example, the Gulf of St-Lawrence where the reverse is often the case. Higher sex ratios are usually found in inshore and bottom slope areas. A decline in sex ratios has been observed since 2017 in N-ENS. In S-ENS the sex ratio increased from 20% in 2021 to just under 35% in 2022. In 4X, mature sex ratios are more balanced and currently near the 50% level.   


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
 


## Sex ratios


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
- Reference: -@fig-R0-timeseries, -@fig-R0-map


The fishable component is defined as Snow Crab that are male, mature, and larger than 95 mm CW. Fishable biomass density is the geometric mean biomass per unit swept area by the trawl. A peak in biomass density was observed in 2009 to 2014 and has since been declining in all areas. Note that high and low biomass density areas fluctuate with time. Biomass density, however, does not equate to total biomass as the areas occupied by crab can contract, expand and shift with environmental conditions and ecosystem variability. 


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


 

## Bottom Temperatures

- N-ENS: declined to historical lows since historical highs in 2022
- S-ENS: declined to historical lows since historical highs in 2022
- CFA 4X: declined to historical lows since historical highs in 2022
- Reference: -@fig-bottom-temperatures-timeseries, -@fig-figures-temperature-bottom-map

Average bottom temperatures **observed** in the 2022 Snow Crab survey were near or above historical highs in all areas. A general warming trend has been observed in the Snow Crab survey since the early 1990s on the Scotian Shelf (Choi et al. 2022). Temperatures are more stable in N-ENS than S-ENS; 4X exhibits the most erratic and highest annual mean bottom temperatures. Of particular note, the observed temperatures in the 2022 Snow Crab survey for S-ENS increased well above the average. This is expected to be a source of bias and certainty for S-ENS abundance predictions as it is outside the range normally encountered by the statistical models. Furthermore, the Groundfish surveys did not operate over Snow Crab grounds in 2020 and 2022 and so temperature data is also very sparse for the area of interest.

Upon aggregation and modelling of historical temperature data, the average temperature is found to have increased well beyond the $7^\circ$C threshold in 4X. N-ENS and S-ENS also continued to experience historical highs in bottom temperature and elevated spatial variability of bottom temperatures (). In particular, the spike observed in S-ENS (above) was less extreme and so likely due to a short term influx of warm water. 

Overall, since 1999, there has been a persistent spatial gradient of almost $15^\circ$C in bottom temperatures in the Maritimes Region. This large gradient in a spatially complex and temporally dynamic area makes assessment a particular challenge for a stenothermic organism such as Snow Crab (Choi et al. 2022). 



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

knitr::include_graphics( file.path( tloc, fns) )

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

- Predators, preys, competitors

  - N-ENS:
  - S-ENS:
  - CFA 4X:
  - Reference: 

Being long-lived, the influence of predators can be significant. Especially important are predators of the smaller immature and female Snow Crab. Increasing predation not only lowers the abundance and recruitment, it can also reduce the reproductive potential of Snow Crab and therefore long-term population cycles. N-ENS and S-ENS are well known to have skewed sex ratios with few mature females for extended periods of time, quite possibly linked to differential predation mortality (mature females being much smaller and full of fat-rich gonads and eggs). 

Based on stomach sampling, Atlantic Halibut, Atlantic Wolffish, Thorny Skate, and other skate species are known predators of Snow Crab. Large Atlantic Halibut with mature female Snow Crab in their stomachs have been reported. Anecdotal information of some seals having fed upon Snow Crab are also known. 

Some of these predators (e.g., Halibut; DFO 2018) have significantly increased in abundance in the Region. However, for all, the abundance and encounter rates in areas overlapping with snow crab habitat is more important, but this is not known. We do know from the bycatch in the Snow Crab survey that there are elevated areal **densities** with many snow crab trawl samples. This means that encounter rates will also likely increase and so too potentially predation mortality. However, high density does not equate to high abundance nor high total predation mortality; but this remains unknown and requires further analysis. The following presents information of areal density of co-occuring species; they are potential predators, competitors and prey. 

Atlantic Halibut densities have increased rapidly since 2010 (-@fig-halibut-timeseries; DFO 2018) on Snow Crab grounds. Most of these increases were towards The Gully, Slope Edge and near Sable Island (-@fig-halibut-map).

Thorny skate densities have been increasing as well (-@fig-thornyskate-timeseries), especially in N-ENS and along the margins of Banquereau Bank (-@fig-thornyskate-map). A minor decline in the densities have been seen in 4X.

Striped Atlantic Wolffish densities have been high, though declining in N-ENS since 2007 (-@fig-stripatlwolffish-timeseries). Highest densities were towards the Laurentian Channel (-@fig-stripatlwolffish-map).

Northern shrimp co-occur as they share similar habitat preferences and are also potential prey items of Snow Crab. Their numerical densities have declined after a peak in 2011, especially in S-ENS (-@fig-northernshrimp-timeseries, -@fig-northernshrimp-map). 

- Lesser toad crab is a co-occurring species and potential competitor. Their numbers have declined to low levels throughout, after a peak in densities in 2007 and 2016 in N-ENS (-@fig-lyrecrab-timeseries, -@fig-lyrecrab-map).

- Overall, higher predation mortality seems likely in N- and S-ENS and lower in 4X. Further shrimp with similar habitat preferences have declined, possibly due to large-scaled habitat variations and predation.



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

 
```{r}
#| label: fig-speciescomposition-timeseries
#| eval: true
#| echo: false 
#| output: true
#| fig.show: hold 
#| fig-dpi: 144
#| fig-height: 6
#| fig-cap: "Species composition in time. Primary gradient (PC1) is related to bottom temperatures; second (PC2) to depth. Groundfish surveys were not conducted in 2020 and 2022 in the snow crab domain. Snow crab surveys were not conducted in 2020, and incomplete in 2022 in S-ENS."
#| fig-subcap: 
#|   - "Mean annual PC1 score."
#|   - "Mean annual PC2 score."
 
pc = c(1, 2)

spc_loc = file.path( data_root, "aegis", "speciescomposition", "modelled", "default" )
fnpr = file.path( spc_loc, "figures", paste("pca", pc, "_time.png", sep="" ) )
knitr::include_graphics( fnpr ) 

``` 
 
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

knitr::include_graphics( fns ) 

``` 
 
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

knitr::include_graphics( fns ) 

``` 
 

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




## Viable habitat

- Snow Crab being cold water stenotherms, stability of environmental conditions is critical for their survival. However, it is not just temperature that is in control, but rather, where the temperature is good in combination with substrate, depths, co-occuring species, etc. The Maritimes Region being at the confluence of many oceanic currents renders the area highly variable. Rapid specis composition change and climate change and uncertainty exacerbates this situation. The viable habitat estimated for each area across time has shown some variations in the historical record. 

- 4X showed a significantly lower average viable habitat levels relative to the N-ENS and S-ENS.  

- A peak in average probability of observing fishable snow crab ("viable habitat") was observed in 2010 for 4X, 2011 for N-ENS and 2012 for S-ENS. Since 2015, the average viable habitat has declined to historical lows and remained so in 2022. 

- N-ENS: Pr has declined to near historical lows (peaks in 2011 and 2020); peaks in inner trench and GBH
- S-ENS: Pr has declined to near historical lows (peaks in 2012); peaks near Sable and missaine
- CFA 4X: declined to near historical lows (peaks in 2010); peak near Lunenburg
- Reference: -@fig-fb-habitat-timeseries, -@fig-fb-habitat-map

 

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

 
 

## Fishery model diagnostics


```{r}
#| label: fig-logistic-surplus-production
#| results: asis
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

fns = file.path( fm_loc, paste("plot_surplus_production_", regions, ".png", sep="" ) ) 

include_graphics( fns )
   
``` 




```{r}
#| label: fig-logistic-prior-posterior-K
#| results: asis
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



```{r}
#| label: fig-logistic-prior-posterior-r
#| results: asis
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

 
```{r}
#| label: fig-logistic-prior-posterior-q
#| results: asis
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



```{r}
#| label: fig-logistic-prior-posterior-obserror
#| results: asis
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
 


```{r}
#| label: fig-logistic-prior-posterior-processerror
#| results: asis
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
 
  
```{r}
#| label: fig-logistic-state-space
#| results: asis
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
 


 

## Fishable biomass index

- N-ENS: decline to historical lows; decline throughout area
- S-ENS: increased marginally; increase offshore (marginally, north of sable)
- CFA 4X: increased marginally; area south of sambro and lunenburg
- Reference:  -@fig-fbindex-timeseries, -@fig-fbindex-map

The fishable **biomass index** ( ) are optimistically high as the spatial expansion uses areal units with large surface areas, larger than the patchiness of Snow Crab distributions (spatial autocorrelation length is <20 km, on average; Choi 2020). As such, it should only be seen as a spatially and temporally comparable relative index of abundance.  

The spatial distribution of the biomass index has been consistent over the past six years, with a peak in overall biomass index in 2019 and 2020. Since then, a  reduction was observed throughout the region, with the exception of the core areas. A contraction of spatial range in 4X and the western parts of S-ENS were also evident in 2021 to 2022. Upon aggregation, the biomass index declined marginally in all areas  



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
  

## Posterior estimates of fishable biomass

- N-ENS: marginal increase in prefishery biomass
- S-ENS: marginal increase in prefishery biomass
- CFA 4X: marginal increase in prefishery biomass
- Reference: -@fig-logisticPredictions

The biomass index along with fishery removals are used to fit a *Logistic Biomass Dynamics Model* to determine fishable **modelled biomass** (biomass estimated from the fisheries model) and relevant biological reference points (i.e., carrying capacity and fishing mortality at maximum sustainable yield, or F~MSY~). In N-ENS, the modelled biomass (pre-fishery) of Snow Crab in `r year_assessment` was `r round(B_north[t0], 2)` t, relative to `r round(B_north[t1], 2)` t in `r year_previous`. In S-ENS, the `r year_assessment` modelled biomass (pre-fishery) was `r round(B_south[t0], 2)` t, relative to `r round(B_south[t1], 2)` t in `r year_previous`. In 4X, the modelled biomass (pre-fishery) for the `r year_assessment`-`r year_assessment+1` season was `r round(B_4x[t0], 2)` t, relative to `r round(B_4x[t1], 2)` t for the `r year_previous`-`r year_assessment` season. 
  


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

- N-ENS: increased
- S-ENS: decreased
- CFA 4X: decreased
- Reference: -@fig-logisticFishingMortality 

In N-ENS, the `r year_assessment` fishing mortality is estimated to have been `r round(FM_north[t0],3)` (annual exploitation rate of `r round(100*(exp(FM_north[t0])-1),2)`%), up from the  `r year_previous` rate of `r round(FM_north[t1],3)` (annual exploitation rate of `r round(100*(exp(FM_north[t1])-1),1)`%.

In S-ENS, the `r year_assessment` fishing mortality is estimated to have been `r round(FM_south[t0],3)` (annual exploitation rate of `r round(100*(exp(FM_south[t0])-1),1)`%), decreasing marginally from the `r year_previous` rate of `r round(FM_south[t1],3)` (annual exploitation rate of `r round(100*(exp(FM_south[t1])-1),1)`%. Localized exploitation rates are likely higher, as not all areas for which biomass is estimated are fished (e.g., continental slope areas and western, inshore areas of CFA 24).

In 4X, the `r year_assessment`-`r year_assessment+1` season (ongoing), fishing mortality is estimated to have been `r round(FM_4x[t0],3)` (annual exploitation rate of `r round(100*(exp(FM_4x[t0])-1),1)`%), decreasing from the `r year_assessment-1`-`r year_assessment` season rate of `r round(FM_4x[t1],3)` (annual exploitation rate of `r round(100*(exp(FM_4x[t1])-1),1)`%. Localized exploitation rates are likely higher, as not all areas for which biomass is estimated are fished. 


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

- N-ENS: In cautious zone; FM is high; 
- S-ENS: In healthy zone
- CFA 4X: In critical zone
- Reference: -@fig-logistic-hcr

Reference points are used to guide harvest strategies (Canada Gazette 2022; DFO 2013). More specifically, the state of the fishery relative to the following "Reference Points" are assessed. In terms of modelled biomass, Lower and Upper Stock Reference are 25% and 50% of carrying capacity which delineate "critical", "cautious" and "healthy" zones. In terms of exploitation rates, the Upper Removal Reference is the exploitation rate that we try not to go beyond; it is defined in term of the fishing mortality associated with Maximum Sustainable Yield (FMSY). In the biomass dynamics model, FMSY = $r /2$. As $r \approx 1$ for snow crab, FMSY $\approx 0.5$ is expected. 

The operational target exploitation changes depending upon the "zone" in which a population lands. When in the "healthy" zone, the rule of thumb has been to keep annual exploitation rates between 10% to 32% of the available biomass ($F = 0.11, 0.36$, respectively). In the "cautious" zone, the rule of thumb has been to keep annual exploitation rates between 0% to 20% ($F = 0, 0.22$, respectively). In the "critical" zone, fishery closure is considered until recovery is observed, where recovery indicates at a minimum, modelled biomass > LSR. Other biological and ecosystem considerations such as recruitment, spawning stock (female) biomass, size structure, sex ratios and environmental and ecosystem conditions, provide additional guidance within each range.

Model 1 estimates key Reference Points as shown in -@tbl-reference-points and -@fig-ReferencePoints. The related PA thresholds can be computed as:

  - Lower Stock Reference (LSR): $K/4$
  - Upper Stock Reference (USR): $K/2$
  - Upper Removal Reference (URR): keep fishing mortality below *FMSY* = $r/2$

Model 1 suggests the current state of the fishable components to be (-@fig-logistic-hcr):  

  - N-ENS is in the "healthy" zone
  - S-ENS is in the "healthy" zone
  - 4X is in the "cautious" zone

It should be emphasized that using these parameters assumes that the population dynamics are well described by the fishery model. This is, of course, **not true**. For example, the observation of fisheries landings is assumed to be known without error. This is not true as illegal and unreported exploitation occurs. These and other unaccounted factors (recruitment strength, environmental variability, predation intensity, disease) can easily bias parameter estimates. As such, **caution is required in using these reference points.** Other contextual reference points must be used in conjunction:

  - Strength of recruitment (short-term, long-term)
  - Strength of spawning stock (females)
  - Ecosystem variability (predator and prey trends and distributions) within norms
  - Habitat viability within norms
  - Availability of spatial and temporal refugia within norms 

We turn to some these additional factors in the next section.

 
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
 
- Capture of soft-shell Snow Crab is always a concern. Prompt and careful return of immature (small-clawed, non-terminally molted) crab to the water is an important conservation measure that will enhance the 2-3 year productivity of the fishable component.

- Bycatch of Snow Crab in other fisheries remains an area requiring attention. Anecdotal information from fishers suggest illegal retention of sublegal and female Snow Crab as bycatch and their use as bait. Illegal removals are also known to occur; however, the scale of such activities are not known.

- Illegal, unreported, and unregulated fishing activities are known to occur. Such activities hinder the application of a precautionary approach to the management of this resource and cause potential bias and uncertainty in Reference Point estimation.

- Marine Protected Areas (MPAs) continue to be developed (e.g., Canada Gazette 2016). The presence of a refuge from fishing activities is potentially a direct positive effect upon Snow Crab. However, positive effects upon other organisms (predators or prey) can have counter-balancing effects. The overall long-term effects of the MPAs upon Snow Crab are unknown.



# Conclusions  

- N-ENS:
- S-ENS:
- CFA 4X:
- Reference: 

The ESS ecosystem is still experiencing a lot of volatility driven by rapid ecosystem and climatic variations. Under such conditions, it is prudent to be careful. Further, the overall indications of population status suggest that Snow Crab are still able to persist under extreme conditions if they are episodic, albeit, with some shifts in spatial distribution towards cooler and deeper waters. 

The modelled solutions represent a few of many possible views of the dynamics of snow crab. Over-emphasis of any one of these modelled solutions and associated Reference Points in determining a strategy for fisheries management is **not prudent** and certainly not precautionary. 


- In N-ENS, though recruitment continues at low levels, a gap in future recruitment to the fishery is expected for the next 1-3 years in N-ENS. N-ENS still exists in the healthy zone and so still has flexibility in approach. However, fishing mortality may have been higher than desirable. A more conservative harvest strategy would enable bridging this coming recruitment gap. A reduced TAC is prudent. 

- In S-ENS, recruitment to the fishery is likely to continue at a moderate rate for the upcoming season. The S-ENS stock remains in the healthy zone. Exploitation rates derived from the fishery model have been declining in recent years. S-ENS has considerable flexibility in approach. All TAC directions are viable. 

- In 4X, low to moderate levels of recruitment are expected for 2 years. 4X exists in the "cautious zone". The area is also in the southern-most extent of Snow Crab distribution in the North Atlantic and viable habitat has been depressed for many years. A reduced TAC is prudent. 



  


## Supplemental Information 
  


## Connectivity: Oceanic currents {.c}

```{r ocean_currents, echo=FALSE, out.width='55%', fig.align='center', fig.show='hold',  fig.cap = 'Ocean currents in the Martimes. Source: https://www.dfo-mpo.gc.ca/oceans/publications/soto-rceo/2018/atlantic-ecosystems-ecosystemes-atlantiques/index-eng.html' }
fn2=file.path( media_loc, "maritimes_currents.png" )
knitr::include_graphics( c(fn2) ) 
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
knitr::include_graphics( fn ) 

```
 


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

knitr::include_graphics( fn ) 

```
 
  
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
knitr::include_graphics( fn ) 

```
 


## Connectivity: Movement {.c}

::: columns 
:::: column 
 
```{r movementtracks, echo=FALSE, out.width='71%', fig.align='center', fig.show='hold',  fig.cap = 'Snow Crab movement tracks from mark and recapture.' }
fn1=file.path( media_loc, "movement0.png" )
fn2=file.path( media_loc, "movement.png" )
knitr::include_graphics( c(fn1, fn2) ) 
# \@ref(fig:movementtracks)  
``` 
 
 
::::
:::: column

```{r movement, echo=FALSE, out.width='55%', fig.align='center', fig.show='hold',  fig.cap = 'Snow Crab movement distances and rates.' }
fn1=file.path( media_loc, "snowcrab_movement_distances.png" )
fn2=file.path( media_loc, "snowcrab_movement_rates.png" )
knitr::include_graphics( c(fn1, fn2) ) 
# \@ref(fig:movement)  
``` 
 
::::
:::


## Bathymetry {.c}
::: columns 

:::: column 

```{r bathymetry-map, out.width='100%', echo=FALSE, fig.align='center', fig.cap='Variations in log(depth; m) in the Scotian Shelf region.' }
bathydir = file.path( data_root, 'aegis', 'bathymetry', 'modelled', 'default', 'stmv', 'none_fft', 'z', 'maps', 'SSE' )
knitr::include_graphics( file.path( bathydir, 'bathymetry.z.SSE.png' ) )
```
::::
:::: column
```{r bathymetry-map2, out.width='100%', echo=FALSE, fig.align='center', fig.cap='Local SD log(depth; m) in the Scotian Shelf region.' }
bathydir = file.path( data_root, 'aegis', 'bathymetry', 'modelled', 'default', 'stmv', 'none_fft', 'z', 'maps', 'SSE' )
knitr::include_graphics( file.path( bathydir, 'bathymetry.b.sdSpatial.SSE.png' ) )
```

::::
:::


## Substrate {.c}
::: columns 

:::: column 

```{r substrate-map, out.width='100%', echo=FALSE, fig.align='center', fig.cap='Substrate grain size log(mm) variations in the Scotian Shelf region.' }
substrdir = file.path( data_root, 'aegis', 'substrate', 'maps', 'canada.east.highres' )
knitr::include_graphics( file.path(  substrdir, 'substrate.substrate.grainsize.canada.east.highres.png' ) )
```
::::
:::: column
```{r substrate-map2, out.width='100%', echo=FALSE, fig.align='center', fig.cap='Local SD substrate grain size log(mm) in the Scotian Shelf region.' }
substrdir = file.path( data_root, 'aegis', 'substrate', 'maps', 'canada.east.highres' )
knitr::include_graphics( file.path( substrdir, 'substrate.s.sdSpatial.canada.east.highres.png' ) )
```

::::
:::


## Bottom Temperature {.c}
::: columns 

:::: column 
\vspace{12mm}
```{r bottom-temperatures-survey, out.width='90%', echo=FALSE, fig.align='center', fig.cap = 'Annual variations in bottom temperature observed during the Snow Crab survey. The horizontal (black) line indicates the long-term, median temperature within each subarea. Error bars represent standard errors.' }
knitr::include_graphics( file.path( data_loc, 'assessments', year_assessment, 'timeseries', 'survey', 't.png') )
# \@ref(fig:bottom-temperatures-survey)
```
::::
 
:::: column 
```{r bottom-temperatures, out.width='75%', echo=FALSE, fig.align='center', fig.cap = 'Posterior densities of predicted average bottom temperatures. Red horizontal line is at $7^\\circ$C.' }
knitr::include_graphics( file.path( data_loc, 'assessments', year_assessment, 'timeseries', 'temperature_bottom.png') )
# \@ref(fig:bottom-temperatures)
```
::::
:::


## Bottom Temperature ... {.c}

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
knitr::include_graphics( c( fn3, fn2, fn1) )
# \@ref(fig:bottom-temperatures-map)
# *Spatial variations in bottom temperature estimated from a historical analysis of temperature data for 1 September.*
```

 
 
## Bottom Temperature ... {.c}
::: columns 
:::: column 

```{r bottom-temperatures-spatialeffect, out.width='90%', echo=FALSE, fig.align='center', fig.cap = 'Persistent spatial effect of bottom temperature, adjusting for spatiotemporal variability and autocorrelations. Time period from 1999 to present.' }
loc = file.path( data_root, 'aegis', 'temperature', 'modelled', 'default', 'maps' )
knitr::include_graphics( file.path( loc, 'space_re_total.png') )
# \@ref(fig:bottom-temperatures-spatialeffect)
```
 
::::
 
:::: column


\vspace{12mm}

- Persistent spatial gradient of  $>1^\circ$C in bottom temperatures in the Maritimes Region. 

- Variable due to confluence:

  - Warm, high salinity Gulf Stream from the S-SE along the shelf edge 

  - Cold, low salinity Labrador Current

  - Cold low salinity St. Lawrence outflow from the N-NE

  - Nearshore Nova Scotia current, running from the NE. 

::::
:::



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

knitr::include_graphics( fns ) 

``` 

## Habitat



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
knitr::include_graphics( c( fn2 ) ) 

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
knitr::include_graphics( c(fn1  ) ) 

```

## Co-occurring species

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


## Bycatch Maritimes {.c}
 
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


## Bycatch Maritimes ... {.c}
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


 


## Bycatch N-ENS ... {.c}
```{r bycatch-speciesordination_n-ens, out.width='70%', echo=FALSE, fig.align='center', fig.cap = 'Bycatch as potentially interacting species in N-ENS.' }
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

     

## Bycatch S-ENS ... {.c}
```{r bycatch-speciesordination_s-ens, out.width='70%', echo=FALSE, fig.align='center', fig.cap = 'Bycatch as potentially interacting species in S-ENS.' }
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
 

## Bycatch 4X ... {.c}
```{r bycatch-speciesordination_4x, out.width='70%', echo=FALSE, fig.align='center', fig.cap = 'Bycatch as potentially interacting species in 4X.' }
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

```{r predators2, echo=FALSE} 
kable( counts[1:11,], format="simple", row.names=FALSE)
```
::::
:::

       
 
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
 
      
