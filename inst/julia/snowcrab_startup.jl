
  
  # packages and directories
  
  print( "\n\nWARNING: if this is the initial run, it will take a while to precompile/download all libraries. \n\n" )
 

  import Pkg  # or using Pkg
  Pkg.activate(julia_scripts_location)  # so now you activate the package
  Base.active_project()  
  push!(LOAD_PATH, julia_scripts_location)  # add the directory to the load path, so it can be found
  cd( julia_scripts_location )  

  
  if ! @isdefined bio_data_directory 
    bio_data_directory = joinpath( homedir(), "bio.data" )  
  end

   
  if ! @isdefined outputs_directory 
    # tailor to your specific installation
    outputs_directory = joinpath( bio_data_directory, "bio.snowcrab",  "fishery_model" ) 
  end

  mkpath(outputs_directory)
  print( "outputs_directory: ", outputs_directory, "\n\n" )

  pkgs = [
    "Revise", "MKL", "Logging", "StatsBase", "Statistics", "Distributions", "Random", "Setfield", "Memoization",
    "ForwardDiff", "DataFrames", "JLD2", "CSV", "PlotThemes", "Colors", "ColorSchemes", "RData",  
    "Plots", "StatsPlots", "MultivariateStats", "StaticArrays", "LazyArrays", "FillArrays",
    "DynamicHMC", "Turing", "ModelingToolkit", "DifferentialEquations", "Interpolations", "LinearAlgebra"
  ]

 
  # for pk in pkgs; @eval using $(Symbol(pk)); end   # Pkg.add( pkgs ) # add required packages

  # automatic load of things related to all projects go here

  current_directory =  @__DIR__() 
  print( "Current directory is: ", current_directory, "\n\n" )

  pkgs_startup = [  
      "Revise", "Logging", "OhMyREPL",
      "Setfield", "Memoization",
      "DataFrames", "CSV", "RData",    
      "StatsBase", "Statistics", "Distributions", "Random", 
      "PlotThemes", "Colors", "ColorSchemes", 
      "Plots", "StatsPlots" 
  ]

  pkgs = unique!( [pkgs_startup; pkgs] )
   
  print( "Loading libraries:\n\n" )

  # load libs and check settings
  # pkgs are defined in snowcrab_startup.jl
  using Pkg
  for pk in pkgs; 
      if Base.find_package(pk) === nothing
          Pkg.add(pk)
      else
        @eval using $(Symbol(pk)); 
      end
  end   # Pkg.add( pkgs ) # add required packages

 