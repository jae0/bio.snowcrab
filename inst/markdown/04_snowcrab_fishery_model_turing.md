---
title: "Snow Crab Fishery Model"
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
# date: "`r format(Sys.time(), '%d %B, %Y')`"
date: last-modified
date-format: "YYYY-MM-D"
keywords: 
  - snow crab fishery model
  - aggregate biomass and fishing mortality estimates
abstract: |
  Fishery model assimilation of fishery data and survey abundance using a biomass dynamics logistic surplus production model.   
toc: true
number-sections: true
highlight-style: pygments
highlight-style: pygments
number-sections: true
editor:
  render-on-save: false
  markdown:
    references:
      location: document
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
  year_assessment: 2024
  year_start: 1999
  media_loc: "media"
  sens: 1
  debugging: FALSE
  model_variation: logistic_discrete_historical
---
 
   
<!-- Preamble

Could be run as an automated process but probably better to run step wise in case of tweaks being needed. 
 
-->


## Purpose

Prepare aggregated time-series data for inference using a Bayesian fishery model. Currently, the default model is the discrete logistic formulation. It is labelled, "logistic_historical" in the following. There are other variations possible. See the calling code. 

---

## Prepare data for discrete logistic models

This step requires:the following to have been completed:

- [bio.snowcrab/inst/markdown/03_biomass_index_carstm.md](https://github.com/jae0/bio.snowcrab/inst/markdown/03_biomass_index_carstm.md) 
  
It creates the survey biomass index for each area.


## Model fitting and parameter inference

Fitting and inference from a latent "state-space" fishery dynamics model. This is implemented in (Julia)[https://julialang.org/]. Previous versions used STAN, JAGS and BUGS as platforms. The Julia/Turing implementation is the most stable and fast. 

First ensure you have Julia installed.  Then you can run the model inference in two ways: 

#### Run within Julia directly. 

This method is most flexible and [documented here](https://github.com/jae0/dynamical_model/blob/master/snowcrab/04_snowcrab_fishery_model.md), where more options are also shown.  But it will require learning the language Julia and Turing library, but is very simple.

```{julia}
#| eval: false 
#| output: false

  # ---- not clear if params can be passed to julia so make sure to update this
  model_variation = "logistic_discrete_historical"   
  year_assessment = 2024   

  yrs = 2000:year_assessment

  outformat = "png"  # for figures .. also, pdf, svg, etc...
  
  project_directory  = joinpath( homedir(), "bio", "bio.snowcrab", "inst", "julia" ) 
  bio_data_directory = joinpath( homedir(), "bio.data", "bio.snowcrab", "modelled", "default_fb" )
  outputs_dir = joinpath( homedir(), "bio.data", "bio.snowcrab", "fishery_model", string(year_assessment), model_variation )
 
  # load packages ... might need to re-run the following if this is your first time and also:   install_required_packages(pkgs)
  include( joinpath(project_directory, "startup.jl" ) )  # load libs
  
  # load core modelling functions
  include( joinpath(project_directory, "logistic_discrete_functions.jl" ) )  

  # load data; alter file path as required   
  o = load( joinpath( bio_data_directory, "biodyn_biomass.RData" ), convert=true)   
 
  # run the whole script (below) or run in parts inside the file one line at a time for more control
  aulab ="cfanorth"     
  include( joinpath(project_directory, "logistic_discrete.jl" ) )   

  aulab ="cfasouth"    
  include( joinpath(project_directory, "logistic_discrete.jl" ) )   

  aulab ="cfa4x"      
  include( joinpath(project_directory, "logistic_discrete.jl" ) )   


```



 **OR,**

#### Run within R and call Julia indirectly.

This requires the use of the JuliaCall R-library (install that too).  As this and the other Julia files exist in a Git repository and we do not want to add more files there, we copy the necessary files to a temporary work location such that new output files (figures, html documents, Julia logs, etc.) can be created there. 
   
The code follows but has been commented out from the report to keep it tidy ... you will have to open the file directly in a text editor and step through it.

<!--

```{r}
#| eval: false
#| output: false

  # Define location of data and outputs
  outputs_dir = file.path( data_root, "bio.snowcrab", "fishery_model", year_assessment, params$model_variation )
  work_dir = file.path( work_root, "fishery_model" )
  
  dir.create( outputs_dir, recursive=TRUE )
  dir.create( work_dir, recursive=TRUE)

  # copy julia scripts to a temp location where julia can write to (R-library is usually write protected)
  julia_scripts_location = system.file( "julia", package="bio.snowcrab" )
  file.copy( file.path(julia_scripts_location, "snowcrab_startup.jl"), work_dir, overwrite=TRUE  )
  file.copy( file.path(julia_scripts_location, "logistic_discrete_functions.jl"), work_dir, overwrite=TRUE )
  file.copy( file.path(julia_scripts_location, "logistic_discrete.jl"), work_dir, overwrite=TRUE )
 
      # to set up:
      # julia_setup(JULIA_HOME = "the folder that contains julia binary")
      # options(JULIA_HOME = "the folder that contains julia binary")
      # Set JULIA_HOME in command line environment.

 
  if ( !any( grepl("JuliaCall", o[,"Package"] ) ) ) install.packages("JuliaCall")

  # load JuliaCall interface
  library(JuliaCall)

  julia = try( julia_setup( install=FALSE, installJulia=FALSE ) )

  if ( inherits(julia, "try-error") ) {
    install_julia()
    julia = try( julia_setup( install=FALSE, installJulia=FALSE ) ) # make sure it works
    if ( inherits(julia, "try-error") )  stop( "Julia install failed, install manually?") 
  }

  ## transfer params to julia environment
  julia_assign( "year_assessment", year_assessment )  # copy data into julia session
  julia_assign( "bio_data_directory", data_root) 
  julia_assign( "outputs_dir", outputs_dir )  # copy data into julia session

  # julia_source cannot traverse directories .. temporarily switch directory

  currentwd = getwd() 
   
  setwd(work_dir) 

    # load/install libraries and setup directories 
    # .. if not all libraries are installed, you might need to re-run this a few times 
    # .. until there are no longer any pre-compilation messages
    julia_source( "snowcrab_startup.jl" )  
    
    # load core modelling functions
    julia_source( "logistic_discrete_functions.jl" )  

    # load data; alter file path as required
    julia_assign( "fndat", fndat )  # load data into julia session
    julia_command("o = load( fndat, convert=true) ")
 
    # Model and compute predictions for each area
    julia_command( "Random.seed!(1234)" )

    # modelling and outputs for each region 
    julia_assign( "yrs", p$yrs )   

    # cfanorth:    
    julia_assign( "aulab", "cfanorth" )   
    julia_source( "logistic_discrete.jl" )   

    # cfasouth:    
    julia_assign( "aulab", "cfasouth")  # copy data into julia session
    julia_source( "logistic_discrete.jl" )   # modelling and outputs

    # cfa4x:    
    julia_assign( "aulab", "cfa4x" )  # copy data into julia session
    julia_source( "logistic_discrete.jl" )   # modelling and outputs

  setwd(currentwd) # revert work directory

  # Example:
  # to move data into R
  # objects that might be useful:  m, num, bio, trace, Fkt, FR, FM 
  # bio = julia_eval("bio")  
```
 


-->


And that is it. Now we continue to look at some of the main results. This part gets processed as a Quarto document. See Premable (above for instructions on creating a report). 



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
  # dev.args = list(type = "cairo"),
  fig.retina = 2,
  dpi=192
)

require(spsUtil)

quietly = spsUtil::quiet


require(ggplot2)
require(MBA)

require(aegis)  # basic helper tools
 
year_assessment = params$year_assessment 

year_previous = year_assessment - 1

p = bio.snowcrab::load.environment( year.assessment=year_assessment )  

rootdir = file.path("/home", "jae" )
 
require(gt)  # table formatting

outtabledir = file.path( p$annual.results, "tables" )

years = as.character(1996: year_assessment)

lregions = list(region=c("cfanorth", "cfasouth", "cfa4x"))
reg_labels = c("N-ENS", "S-ENS", "CFA 4X")  # formatted for label
 
regions = unlist(lregions)
nregions = length(regions)
 
# in case of local changes
loadfunctions( "aegis")
loadfunctions( "bio.snowcrab")  # in case of local edits
 
p$corners = data.frame(plon=c(220, 990), plat=c(4750, 5270) )
p$mapyears = year_assessment + c(-5:0 )   # default in case not specified
 

loc = file.path( rootdir, "bio.data", "bio.snowcrab", 'fishery_model', year_assessment, params$model_variation )
 
```

  



## Stock status: Biomass Index (aggregate) {.c}
    
```{r}
#| results: asis
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| label: fbindex-timeseries

fn = file.path( rootdir, "bio.data", "bio.snowcrab", 'modelled', 'default_fb', 'aggregated_biomass_timeseries' , 'biomass_M0.png')

include_graphics( fn )

```

The fishable biomass index (t) predicted by CARSTM of Snow Crab survey densities. Error bars represent Bayesian 95\\% Credible Intervals. Note large errors in 2020 when there was no survey.
 

## Stock status: Modelled Biomass (pre-fishery) {.c}
```{r}
#| results: asis
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| label: logisticPredictions

 
for (r in 1:nregions){
  reg = regions[r]
  REG = reg_labels[r] 
  cat("#### ", REG, "\n")
  fn = file.path( loc, paste("plot_predictions_", reg, ".png", sep="" ) ) 
  show_image(fn) 
  cat("\n\n")
} 


```
Model 1 fishable, posterior mean modelled biomass (pre-fishery; kt) are shown in dark orange. Light orange are posterior samples of modelled biomass (pre-fishery; kt) to illustrate the variability of the predictions. The biomass index (post-fishery, except prior to 2004) after model adjustment by the model catchability coefficient is in gray.



## Stock status: Fishing Mortality ... {.c}

```{r}
#| results: asis
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| label: logisticFishingMortality



for (r in 1:nregions){
  reg = regions[r]
  REG = reg_labels[r] 
  cat("#### ", REG, "\n")
  fn = file.path( loc, paste("plot_fishing_mortality_", reg, ".png", sep="" ) ) 
  show_image(check_file_exists(fn)) 
  cat("\n\n")
} 

``` 
Time-series of modelled instantaneous fishing mortality from Model 1. Samples of the posterior densities are presented, with the darkest line being the mean.


## Stock status: Reference Points {.c}

```{r}
#| results: asis
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| label: logistic-hcr



for (r in 1:nregions){
  reg = regions[r]
  REG = reg_labels[r] 
  cat("#### ", REG, "\n")
  fn = file.path( loc, paste("plot_hcr_", reg, ".png", sep="" ) ) 
  show_image(check_file_exists(fn)) 
  cat("\n\n")
} 
  

``` 

Reference Points (fishing mortality and modelled biomass) from Model 1. The large yellow dot indicates most recent year and stars the 95\\% CI.



## Surplus production ("Shaeffer form") {.c}

```{r}
#| results: asis
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| label: logistic-surplus-production



for (r in 1:nregions){
  reg = regions[r]
  REG = reg_labels[r] 
  cat("#### ", REG, "\n")
  fn = file.path( loc, paste("plot_surplus_production_", reg, ".png", sep="" ) ) 
  show_image(check_file_exists(fn)) 
  cat("\n\n")
} 
  
``` 

Surplus production from Model 1. Each year is represented by a different colour.



 
## Prior-posterior comparsons: Carrying capacity (K; kt) {.c}

```{r}
#| results: asis
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| label: logistic-prior-posterior-K


for (r in 1:nregions){
  reg = regions[r]
  REG = reg_labels[r] 
  cat("#### ", REG, "\n")
  fn = file.path( loc, paste("plot_prior_K_", reg, ".png", sep="" ) ) 
  show_image(check_file_exists(fn)) 
  cat("\n\n")
} 
  

``` 

Prior-posterior comparisons from Model 1. Carrying capacity (K; kt).


## Prior-posterior comparsons: Intrinsic rate of biomass increase (r) {.c}

```{r}
#| results: asis
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| label: logistic-prior-posterior-r


for (r in 1:nregions){
  reg = regions[r]
  REG = reg_labels[r] 
  cat("#### ", REG, "\n")
  fn = file.path( loc, paste("plot_prior_r_", reg, ".png", sep="" ) ) 
  show_image(check_file_exists(fn)) 
  cat("\n\n")
} 
  

``` 

Prior-posterior comparisons from Model 1. Intrinsic rate of biomass increase (r).


 
## Prior-posterior comparsons: Catchability coefficient (q) {.c}

```{r}
#| results: asis
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| label: logistic-prior-posterior-q


for (r in 1:nregions){
  reg = regions[r]
  REG = reg_labels[r] 
  cat("#### ", REG, "\n")
  fn = file.path( loc, paste("plot_prior_q1_", reg, ".png", sep="" ) ) 
  show_image(check_file_exists(fn)) 
  cat("\n\n")
} 
  

``` 

Prior-posterior comparisons from Model 1. Catchability coefficient (q).



 
## Prior-posterior comparsons: Observation error  {.c}

```{r}
#| results: asis
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| label: logistic-prior-posterior-obserror


for (r in 1:nregions){
  reg = regions[r]
  REG = reg_labels[r] 
  cat("#### ", REG, "\n")
  fn = file.path( loc, paste("plot_prior_bosd_", reg, ".png", sep="" ) ) 
  show_image(check_file_exists(fn)) 
  cat("\n\n")
} 
  

``` 

Prior-posterior comparisons from Model 1. Observation error (kt).


 
## Prior-posterior comparsons: Model process error {.c}

```{r}
#| results: asis
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| label: logistic-prior-posterior-processerror


for (r in 1:nregions){
  reg = regions[r]
  REG = reg_labels[r] 
  cat("#### ", REG, "\n")
  fn = file.path( loc, paste("plot_prior_bpsd_", reg, ".png", sep="" ) ) 
  show_image(check_file_exists(fn)) 
  cat("\n\n")
} 
  

``` 

Prior-posterior comparisons from Model 1. Model process error.


 
## State space {.c}

```{r}
#| results: asis
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| label: logistic-state-space


for (r in 1:nregions){
  reg = regions[r]
  REG = reg_labels[r] 
  cat("#### ", REG, "\n")
  fn = file.path( loc, paste("plot_state_space_", reg, ".png", sep="" ) ) 
  show_image(check_file_exists(fn)) 
  cat("\n\n")
} 
  

``` 

State space.

