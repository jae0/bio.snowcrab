# packages 

print( "\n\nWARNING: If this is the initial run, it will take a while to precompile/download all libraries. \n\n" )


import Pkg  # or using Pkg
current_directory = @__DIR__()
Pkg.activate(current_directory)  # so now you activate the package
Base.active_project()  
push!(LOAD_PATH, current_directory)  # add the directory to the load path, so it can be found

pkgs = [
    "Pkg",  "Revise", "Logging",
    "StatsBase", "Statistics", "Distributions", "Random", "Setfield", "Memoization",
    "MCMCChains",
    "DataFrames", "JLD2", "CSV", "PlotThemes", "Colors", "ColorSchemes", "RData",
    "Plots",  "StatsPlots", "MultivariateStats",
    "ForwardDiff", "ADTypes",  
    "StaticArrays", "LazyArrays", "FillArrays", "LinearAlgebra", "MKL", "Turing"
]
 

# load directly can cause conflicts due to same function names 
pkgtoskipload = ["CairoMakie", "CairoMakie", "PlotlyJS",  "PlotlyBase",  "PlotlyKaleido" ]


print( "\nWARNING: Consider updating libraries directly from within Julia as you have more control." )

print( "\nLoading libraries:\n" )
 
# load libs and check settings
# pkgs are defined in snowcrab_startup.jl
 
function install_required_packages(pkgs)    # to install packages
    for pk in pkgs; 
        if Base.find_package(pk) === nothing
            Pkg.add(pk)
        end
    end   # Pkg.add( pkgs ) # add required packages

    print( "Pkg.add( \"Bijectors\" , version => \"0.3.16\") # may be required \n" )

end

for pk in pkgs; 
    if !(pk in pkgtoskipload)
        @eval using $(Symbol(pk)); 
    end
end


println( """

Libraries loading:

If this is the initial run, it will take a while to precompile/download all libraries.
The variable, 'pkgs', contains the list of required libraries. To (re)-install required packages, run:

install_required_packages(pkgs)

""" )

 
# Pkg.update()


