---
title: "Snow Crab, Scotian Shelf, Canada (NAFO Div. 4VWX) in 2023"
subtitle: "DDE model solutions"
author: "Snow Crab Group"
# author: "Jae S. Choi"
# footnote: "jae.choi@dfo-mpo.gc.ca"
institute: "Bedford Institute of Oceanography, DFO Science"
# date: "`r format(Sys.time(), '%d %B, %Y')`"
output:
  beamer_presentation:
    theme: "metropolis"
    colortheme: "seagull"
    fonttheme: "professionalfonts"
    fig_caption: yes
    # latex_engine: pdflatex
    latex_engine: lualatex 
    keep_tex: true
classoption: 
  - aspectratio=169 #16:9 wide
  - t  # top align
header-includes: 
  - \usepackage{graphicx}
  - \usepackage[font={scriptsize}, labelfont={bf}]{caption}
  # - \usepackage{float}
  # - \usepackage{subfig}
  # - \newcommand{\btiny}{\begin{tiny}}
  # - \newcommand{\etiny}{\end{tiny}}
params:
  year.assessment: 2023
  media_loc: "media"
  debugging: FALSE
  loc_dde: ""
--- 


<!-- Preamble


This is a Markdown document ... To create HTML or PDF, etc, run: 


  make quarto FN=dde YR=2023 SOURCE=~/projects/bio.snowcrab/inst/markdown WK=~/bio.data/bio.snowcrab/reports # {via Quarto}

  make rmarkdown FN=dde YR=2023 SOURCE=~/projects/bio.snowcrab/inst/markdown WK=~/bio.data/bio.snowcrab/reports {via Rmarkdown}

  make pdf FN=dde  # {via pandoc}

Alter year and directories to reflect setup or copy Makefile and alter defaults to your needs.
  
 

Huge > huge > LARGE > Large > large > normalsize > small > footnotesize > scriptsize > tiny


::: columns 
:::: column 

::::
 
:::: column

::::
:::


-->



<!-- Set up R-environment -->

```{r setup, include=FALSE}
  require(knitr)
  knitr::opts_chunk$set(
    root.dir = data_root,
    echo = FALSE,
    out.width="6.2in",
#     dev.args = list(type = "cairo"),
    fig.retina = 2,
    dpi=192
  )

  # inits and data loading (front load all required data)

  require(aegis)
  
  year.assessment = params$year.assessment
  year_previous = year.assessment - 1
  p = bio.snowcrab::load.environment( year.assessment=year.assessment )
  SCD = project.datadirectory("bio.snowcrab")
  media_loc = params$media_loc
  
  # fishery_model_results = file.path( "/home", "jae", "projects", "dynamical_model", "snowcrab", "outputs" )
  fishery_model_results = file.path( SCD, "fishery_model" )

  sn_env = snowcrab_load_key_results_to_memory( year.assessment, debugging=params$debugging, loc_dde=params$loc_dde, return_as_list=TRUE  ) 

  attach(sn_env)

  # predator diet data
  diet_data_dir = file.path( SCD, "data", "diets" )
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

  # bycatch summaries
  o_cfaall = observer.db( DS="bycatch_summary", p=p,  yrs=p$yrs, region="cfaall" )
  o_cfanorth = observer.db( DS="bycatch_summary", p=p,  yrs=p$yrs, region="cfanorth" )   
  o_cfasouth = observer.db( DS="bycatch_summary", p=p,  yrs=p$yrs, region="cfasouth" )   
  o_cfa4x = observer.db( DS="bycatch_summary", p=p,  yrs=p$yrs, region="cfa4x" )   

```
  
# DDE solutions 
   
##  (N-ENS, S-ENS, 4X)


```{r dde-predictions,  out.width='32%', fig.show='hold', fig.align='center', fig.cap='predictions' }
  loc = file.path( homedir, "projects", "dynamical_model", "snowcrab", "outputs",  year.assessment, "size_structured_dde_normalized" )
  fn1 = file.path( loc, "plot_predictions_cfanorth.pdf" ) 
  fn2 = file.path( loc, "plot_predictions_cfasouth.pdf" ) 
  fn3 = file.path( loc, "plot_predictions_cfa4x.pdf" ) 
  include_graphics(c(fn1, fn2, fn3) )
  # \@ref(fig:dde-predictions)
``` 
   
  
##  Fishery Footprint   {.c}
 
```{r dde-fisheryfootprint,  out.width='32%', fig.show='hold', fig.align='center', fig.cap='Fishery footprint (Model 2). N-ENS (left), S-ENS (middle), and 4X (right). Projections are based upon status quo TACs.'  }
  loc = file.path( homedir, "projects", "dynamical_model", "snowcrab", "outputs",  year.assessment, "size_structured_dde_normalized" )
  fn1 = file.path( loc, "plot_trace_footprint_projections_cfanorth__1.0__.pdf" ) 
  fn2 = file.path( loc, "plot_trace_footprint_projections_cfasouth__1.0__.pdf" ) 
  fn3 = file.path( loc, "plot_trace_footprint_projections_cfa4x__1.0__.pdf" ) 
  include_graphics(c(fn1, fn2, fn3) )
  # \@ref(fig:dde-fisheryfootprint)
``` 


##  (N-ENS, S-ENS, 4X)
```{r dde-hcr,  out.width='32%', fig.show='hold', fig.align='center', fig.cap='hcr' }
  loc = file.path( homedir, "projects", "dynamical_model", "snowcrab", "outputs",  year.assessment, "size_structured_dde_normalized" )
  fn1 = file.path( loc, "plot_hcr_cfanorth.pdf" ) 
  fn2 = file.path( loc, "plot_hcr_cfasouth.pdf" ) 
  fn3 = file.path( loc, "plot_hcr_cfa4x.pdf" ) 
  include_graphics(c(fn1, fn2, fn3) )
  # \@ref(fig:dde-hcr)
```

##  (N-ENS, S-ENS, 4X)
```{r dde-footprint-trace,   out.width='32%', fig.show='hold', fig.align='center', fig.cap='footprint trace' }
  loc = file.path( homedir, "projects", "dynamical_model", "snowcrab", "outputs",  year.assessment, "size_structured_dde_normalized" )
  fn1 = file.path( loc, "plot_footprint_trace_cfanorth.pdf" )
  fn2 = file.path( loc, "plot_footprint_trace_cfasouth.pdf" )
  fn3 = file.path( loc, "plot_footprint_trace_cfa4x.pdf" )
  include_graphics(c(fn1, fn2, fn3) )
  # \@ref(fig:dde-footprint-trace)
``` 


##  (N-ENS, S-ENS, 4X)
```{r dde-hcr-footprint,  out.width='32%', fig.show='hold', fig.align='center', fig.cap='hcr footprint' }
  loc = file.path( homedir, "projects", "dynamical_model", "snowcrab", "outputs",  year.assessment, "size_structured_dde_normalized" )
  fn1 = file.path( loc, "plot_hcr_footprint_cfanorth.pdf" )
  fn2 = file.path( loc, "plot_hcr_footprint_cfasouth.pdf" )
  fn3 = file.path( loc, "plot_hcr_footprint_cfa4x.pdf" )
  include_graphics(c(fn1, fn2, fn3) )
  # \@ref(fig:dde-hcr-footprint)
```

##  (N-ENS, S-ENS, 4X)
```{r dde-predictions-trace, out.width='32%', fig.show='hold', fig.align='center', fig.cap="predictions_trace" }
 loc = file.path( homedir, "projects", "dynamical_model", "snowcrab", "outputs", year.assessment, "size_structured_dde_normalized" )
 include_graphics( file.path( loc, "plot_predictions_trace_cfanorth.pdf" ) )
 include_graphics( file.path( loc, "plot_predictions_trace_cfasouth.pdf" ) ) 
 include_graphics( file.path( loc, "plot_predictions_trace_cfa4x.pdf" ) )
```

