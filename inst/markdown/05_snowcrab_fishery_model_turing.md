# 05_snowcrab_fishery_model_turing.md

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
  year_assessment = 2025   

  yrs = 2000:year_assessment   ## This needs to be consistent with p$fishery_model_years in markdowns in R

  outformat = "png"  # for figures .. also, pdf, svg, etc...
  
  project_directory  = joinpath( homedir(), "bio", "bio.snowcrab", "inst", "julia" ) 
  bio_data_directory = joinpath( homedir(), "bio.data", "bio.snowcrab", "modelled", "default_fb" )
  outputs_dir = joinpath( homedir(), "bio.data", "bio.snowcrab", "fishery_model", string(year_assessment), model_variation )
 
  # load packages ... might need to re-run the following if this is your first time and also:   install_required_packages(pkgs)
  include( joinpath(project_directory, "startup.jl" ) )  # load libs
  
  # load core modelling functions
  include( joinpath(project_directory, "logistic_discrete_functions.jl" ) )  

  # load data; alter file path as required   
  o = load( joinpath( bio_data_directory, "biodyn_biomass.rdz" ), convert=true)   
 
  # run the whole script (below) or run in parts inside the file one line at a time for more control
  aulab ="cfanorth"     
  include( joinpath(project_directory, "logistic_discrete.jl" ) )   

  aulab ="cfasouth"    
  include( joinpath(project_directory, "logistic_discrete.jl" ) )   

  aulab ="cfa4x"      
  include( joinpath(project_directory, "logistic_discrete.jl" ) )   


```



[And that is it. Now we continue to look at some of the main results](06_assessment_summary.md).



 **OR,**


#### Run within R and call Julia indirectly.

This requires the use of the JuliaCall R-library (install that too).  As this and the other Julia files exist in a Git repository and we do not want to add more files there, we copy the necessary files to a temporary work location such that new output files (figures, html documents, Julia logs, etc.) can be created there. 
   
The code follows but has been commented out from the report to keep it tidy ... you will have to open the file directly in a text editor and step through it.
 

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
    p$fishery_model_years = 2000:p$year.assessment   # this needs to be consistent with markdown R values !
    julia_assign( "yrs", p$fishery_model_years )   

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
 

 
