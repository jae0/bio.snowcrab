# packages 

print( "\n\nWARNING: If this is the initial run, it will take a while to precompile/download all libraries. \n\n" )


import Pkg  # or using Pkg
current_directory = @__DIR__()
Pkg.activate(current_directory)  # so now you activate the package
Base.active_project()  
push!(LOAD_PATH, current_directory)  # add the directory to the load path, so it can be found

pkgs = [
  "Revise", "OhMyREPL", "MKL", "Logging", "Setfield", "Memoization",
  "StatsBase", "Statistics", "Distributions", "Random", "MultivariateStats",
  "DataFrames", "RData", "JLD2", "CSV", 
  "PlotThemes", "Colors", "ColorSchemes", "Plots", "StatsPlots", 
  "StaticArrays", "LazyArrays", "FillArrays",
  "ForwardDiff", "DynamicHMC", "DifferentialEquations", "Interpolations", "LinearAlgebra", "Turing", "ModelingToolkit"
]
 
print( "\nWARNING: Consider updating libraries directly from within Julia as you have more control." )
print( "\nWARNING: Otherwise, you might need to re-run snowcrab_startup.jl several times to get libs to install." )

print( "\nLoading libraries:\n" )

# load libs and check settings
using Pkg
for pk in pkgs; 
    if Base.find_package(pk) === nothing
        Pkg.add(pk)
    else
        @eval using $(Symbol(pk)); 
    end
end   

# Pkg.update()


