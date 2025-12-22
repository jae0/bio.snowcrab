# 05_snowcrab_fishery_model_turing.md

## Purpose

Prepare aggregated time-series data for inference using a Bayesian fishery model. Currently, the default model is the discrete logistic formulation. It is labeled, "logistic_historical" in the following. There are other variations possible. See the calling code.

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
  year_assessment = 2025

  yrs = 2000:year_assessment   ## This needs to be consistent with p$fishery_model_years in markdowns in R

  project_directory  = joinpath( homedir(), "bio", "bio.snowcrab", "inst", "julia" )
  bio_data_directory = joinpath( homedir(), "bio.data", "bio.snowcrab", "modelled", "default_fb" )

  model_variation = "logistic_discrete_historical"   
  outputs_dir = joinpath( homedir(), "bio.data", "bio.snowcrab", "fishery_model", string(year_assessment), model_variation )
  
  mkpath( outputs_dir )

  cd( outputs_dir )

  import Pkg  # or using Pkg
  Pkg.activate(outputs_dir)  # so now you activate the package
  Base.active_project()
  push!(LOAD_PATH, outputs_dir)  # add the directory to the load path, so it can be found

  pkgs = [
        "Pkg",  "Revise", "Logging", "StatsBase", "Statistics", "Distributions", "Random", "MCMCChains",
        "DataFrames", "JLD2", "CSV", "PlotThemes", "Colors", "ColorSchemes", "RData", "Setfield", "Memoization",
        "Plots",  "StatsPlots", "MultivariateStats", "ForwardDiff", "ADTypes", 
        "StaticArrays", "LazyArrays", "FillArrays", "LinearAlgebra", "MKL", "Turing"
  ]
        
  for pk in pkgs;
      @eval using $(Symbol(pk));
  end



  # load data; alter file path as required
 
  o = load( joinpath( bio_data_directory, "biodyn_biomass.rds" ), convert=true)
 
  do_model = true
  do_plots = true
    
  aulab ="cfanorth"
  include( joinpath(project_directory, "logistic_discrete.jl" ) )

  aulab ="cfasouth"
  include( joinpath(project_directory, "logistic_discrete.jl" ) )

  aulab ="cfa4x"
  include( joinpath(project_directory, "logistic_discrete.jl" ) )


```



[And that is it. Now we continue to look at some of the main results](06_assessment_summary.md).

 