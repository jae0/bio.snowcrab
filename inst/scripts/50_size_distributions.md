 

# estimate size structure
#   1. direct computation
#   2. modelled (factorial) == 1
#   3. modelled (space-time) using carstm


 ```R
# ---------------------
# define parameters:

    year.assessment = 2022
    p = bio.snowcrab::load.environment( year.assessment=year.assessment )
    loadfunctions( "bio.snowcrab")
    loadfunctions( "aegis")

    require(ggplot2)
    require(data.table)
 
    years = as.character(p$yrs)
    regions=c("cfanorth", "cfasouth", "cfa4x")

    plotdir=file.path( p$annual.results, "figures", "size.freq", "survey")



# ---------------------
# method 0: base data 

    # create base data of individual size records

    pg = areal_units( 
        p = snowcrab_parameters(
            project_class="carstm",
            yrs=1999:year.assessment,   
            areal_units_type="tesselation",
            carstm_model_label=  paste( "1999_present", "fb", sep="_" )
    ))
    dim(pg)

    Y = size_distributions(p=p, toget="base_data", pg=pg, redo=TRUE)

    Y = size_distributions(p=p, toget="tabulated_data", redo=TRUE)


   
# ---------------------
# method 1: direct computation  
    # aggregate across sid = sampling events
    # requires a lot of RAM: 128GB ... 
    # adding zeros

    M = size_distributions(p=p, toget="simple_direct" )

    plotdir=file.path( p$annual.results, "figures", "size.freq", "survey", "direct")

    plot_histogram_carapace_width( M=M, years=years, regions=regions, 
        plot_sex="female", 
        yvar="den",  # den=arithmetic mean density, denl = geometric mean density  
        outdir=plotdir 
    ) 
    
    plot_histogram_carapace_width( M=M, years=years, regions=regions, 
        plot_sex="female", 
        yvar="denl",  # den=arithmetic mean density, denl = geometric mean density  
        outdir=plotdir 
    )

    plot_histogram_carapace_width( M=M, years=years, regions=regions, 
        plot_sex="male", 
        yvar="den",   # den=arithmetic mean density, denl = geometric mean density  
        outdir=plotdir 
    )

    plot_histogram_carapace_width( M=M, years=years, regions=regions, 
        plot_sex="male", 
        yvar="denl",   # den=arithmetic mean density, denl = geometric mean density  
        outdir=plotdir 
    )





    if (0) {
        ggplot( M[ region==regions[a] & year==years[b] ,], aes(cwd, den, fill=mat, colour=mat) ) +
            #geom_ribbon(aes(ymin=density_lb, max=density_ub), alpha=0.2, colour=NA) +
            geom_bar(stat = "identity") +
            labs(x="cw", y="density", size = rel(1.5)) +
            # scale_y_continuous( limits=c(0, 300) )  
            theme_light( base_size = 22 ) 
    }

  

 
# ---------------------
# method 2: simple linear (gaussian) model via biglm .. too slow to use

    O = size_distributions(p=p, toget="linear_model" )
 
    ss = O[ region=="cfanorth" & year== 2017, which=TRUE]
 
    ggplot( O[ ss, ], aes(cwd, den, fill=mat, colour=mat) ) +
        # geom_ribbon(aes(ymin=density_lb, max=density_ub), alpha=0.2, colour=NA) +
        # geom_line() +
        geom_bar(stat = "identity") +
        labs(x="cw", y="density", size = rel(1.5)) +
        # scale_y_continuous( limits=c(0, 300) )  
        theme_light( base_size = 22 ) 



# ---------------------
# method 3: poisson model  via biglm .. problem is too large to compute
    
    # too slow to complete
    O = size_distributions(p=p, toget="poisson_glm" )
 
    plotdir=file.path( p$annual.results, "figures", "size.freq", "survey", "poisson_glm")
    
    regions = "cfanorth"
    plot_histogram_carapace_width( M=O$P, years=years, regions=regions, 
        plot_sex="male", 
        yvar="N",  # den=arithmetic mean density, denl = geometric mean density  
        outdir=plotdir 
    )



# ---------------------
# method 4: poisson via inla .. problem is too large to compute
     
    # adjust based upon RAM requirements and ncores
    require(INLA)
    inla.setOption(num.threads= floor( parallel::detectCores() / 2) )
    
    O = size_distributions(p=p, toget="poisson_inla" )
 
    plotdir=file.path( p$annual.results, "figures", "size.freq", "survey", "poisson_inla")
    
    regions = "cfanorth"
    plot_histogram_carapace_width( M=O$P, years=years, regions=regions, 
        plot_sex="male", 
        yvar="N",  # den=arithmetic mean density, denl = geometric mean density  
        outdir=plotdir 
    )


# ---------------------
# method 5 .. Basic Markov model 

    require(data.cube)  ## as a data.cube 

    # raw data in tabular form
    O = size_distributions(p=p, toget="data_cubes", Y="geometric"  )



# ---------------------
# Modal analysis

# iteratively go to every polygon and identify samples from it and surrounding areas and times (time window is 5 weeks :: +/- 2 weeks)
# and compute modes and troughs
 
    # time window is 11 weeks ~ 3 months 
    # a size bin must have 5 or more individuals and a total of 100 individuals in the whole sample
    
    # NOTE: for bw (bandwidth)
    # diff(range( log10(10:165) )) / 0.025 # = ~ 48 increments  
    # diff(range( log10(10:165) )) / 0.02  # = ~ 60 increments  
    # diff(range( log10(10:165) )) / 0.015  # = ~ 81 increments  
    # diff(range( log10(10:165) )) / 0.01 # 121 increments
    # there are < 13 modes expected so 48 increments or so should be enough to ID modes

    redo = FALSE
    # redo = TRUE
    
    Y = size_distributions(p=p, toget="size_modes", pg=pg, grad_method="Richardson",
        bw=0.01, sigdigits=3, ti_window = c(-6, 6), kexp=1,
        n_min=50, n_cutoff=3, lowpassfilter=0.0001, lowpassfilter2=0.0001, 
        n_neighbours=2, plot_solutions=TRUE, redo=TRUE)s
 
    # all = size_distributions(p=p, toget="size_modes" )
    f = size_distributions(p=p, toget="size_modes", group="female" )
    m = size_distributions(p=p, toget="size_modes", group="male" )
    mi = size_distributions(p=p, toget="size_modes", group="ifemale" )
    fi = size_distributions(p=p, toget="size_modes", group="imale" )

    
    hist( list_get(f,  "peaks"), breaks="fd")
    hist( list_get(m,  "peaks"), breaks="fd") 
    hist( list_get(mi, "peaks"), breaks="fd") 
    hist( list_get(fi, "peaks"), breaks="fd") 
    
    hist( list_get(f,  "troughs"), breaks="fd")
    hist( list_get(m,  "troughs"), breaks="fd") 
    hist( list_get(mi, "troughs"), breaks="fd") 
    hist( list_get(fi, "troughs"), breaks="fd") 
  
    hist( list_get(f,  "peak_values"), breaks="fd")
    hist( list_get(m,  "peak_values"), breaks="fd") 
    hist( list_get(mi, "peak_values"), breaks="fd") 
    hist( list_get(fi, "peak_values"), breaks="fd") 
    

    vf = identify_modes( 
        Z = list_get(f, "peaks"), 
        W = list_get(f, "peak_values"), 
        # T = list_get(f, "troughs"), 
        # V = list_get(f, "trough_values"), 
        lowpassfilter2=0.0001,
        bw=0.02, sigdigits=3, plot=TRUE) 

    vm = identify_modes( 
        Z = list_get(m, "peaks"), 
        W = list_get(m, "peak_values"), 
        # T = list_get(m, "troughs"), 
        # V = list_get(m, "trough_values"), 
        lowpassfilter2=0.0001,
        bw=0.02, sigdigits=3, plot=TRUE) 


    vmi = identify_modes( 
        Z = list_get(mi, "peaks"), 
        W = list_get(mi, "peak_values"), 
        # T = list_get(mi, "troughs"), 
        # V = list_get(mi, "trough_values"), 
        lowpassfilter2=0.0001,
        bw=0.02, sigdigits=3, plot=TRUE) 
 

    vfi = identify_modes( 
        Z = list_get(fi, "peaks"), 
        W = list_get(fi, "pkvalue"), 
        # T = list_get(fi, "troughs"), 
        #V = list_get(fi, "trough_values"),
        #override_range=c(1.5, 2.0), 
        lowpassfilter2=0.0001,
        bw=0.01, sigdigits=3, plot=TRUE) 
 


    out = data.table( rbind( 
        data.table( cw = 10^(vm$peaks),   cat="m"),
        data.table( cw = 10^(vf$peaks),   cat="f"),
        data.table( cw = 10^(vmi$peaks),  cat="mi"),
        data.table( cw = 10^(vfi$peaks),  cat="fi")
    ))
  
    fn = file.path( p$annual.results, "size_distributions_summary.RDS" )
    saveRDS( out, file=fn )



    growth = function(sex, instar) {
        if (sex=="f") cw = exp(2.198848 + 0.315026 * (instar - 4) )
        if (sex=="m") cw = exp(1.917564 + 0.298914 * (instar - 3) )
        return(cw)
    }


    cwbr = data.table( instar=1:14, ml=growth(sex="m", instar=1:14 ), fl=growth(sex="f", instar=1:14) )

    females = rbind(
        data.table( 
            cw = 10^(vfi$peaks),  
            cat="fimm",  
            instar=as.numeric( as.character( cut( 10^(vfi$peaks), breaks= cwbr$fl, labels=cwbr$instar[1:13] ) ) ) ),
        data.table( 
            cw = 10^(vf$peaks),  
            cat="fmat",  
            instar=as.numeric( as.character( cut( 10^(vf$peaks), breaks= cwbr$fl, labels=cwbr$instar[1:13] ) ) ) )
    )


o = females[cat=="fimm",] 
o = o[order(cw)]
n = nrow(o)
plot( o$cw[1:(n-1)] ~ o$cw[2:n] )
ol = loess( o$cw[1:(n-1)] ~ o$cw[2:n] )
op = predict(ol)
lines( op ~ o$cw[2:n]  )

od = diff(o$cw) / o$cw[-nrow(o)]

plot(cw ~ instar, females, type="n")
points(cw ~ instar, females[cat=="fimm"], col="brown", cex=2, pch=20)
points(cw ~ instar, females[cat=="fmat"], col="blue", cex=2, pch=20)


    males = rbind(
        data.table( 
            cw = 10^(vmi$peaks),  
            cat="mimm",  
            instar=as.numeric( as.character( cut( 10^(vmi$peaks), breaks= cwbr$ml, labels=cwbr$instar[1:13] ) ) ) ),
        data.table( 
            cw = 10^(vm$peaks),  
            cat="mmat",  
            instar=as.numeric( as.character( cut( 10^(vm$peaks), breaks= cwbr$ml, labels=cwbr$instar[1:13] ) ) ) )
    )


plot(cw ~ instar, males, type="n")
points(cw ~ instar, males[cat=="mimm"], col="brown", cex=2, pch=20)
points(cw ~ instar, males[cat=="mmat"], col="blue", cex=2, pch=20)

ggplot(df) +
  stat_density(aes(x=d, alpha = id), position = "stack", geom = "line", show.legend = F, color = "red") +
  stat_density(aes(x=d, linetype = id), position = "identity", geom = "line")+
  scale_alpha_manual(values = c(1,0,0,0))

 


```


# Model analysis with Kernel Mixture Modelling in Julia ...

See projects/model.abm/julia/MixtureModels.md for more info.

First prep the data in R:

```R

  Y = size_distributions( p=p, toget="base_data" )
  nb = attributes(pg)$nb$nbs
  au =pg$AUID

  out = list( Y=Y, nb=nb, au=au )
  fn = file.path( p$project.outputdir, "size_structure", "base_data_julia.rdata" )  

  save( out, file=fn )

```


```julia
       
 
pkgs = [
    "Revise", 
    "MKL", # "OpenBLAS32",
    "Logging", "Random", "Setfield", "Memoization", "ForwardDiff", 
    "RData", "DataFrames", "JLD2", "CSV", 
    "PlotThemes", "Colors", "ColorSchemes", "Plots", "StatsPlots", 
    "StatsBase", "Statistics", "Distributions", "KernelDensity",
    "MultivariateStats", "StaticArrays", "LazyArrays", "FillArrays",
    "Interpolations", "LinearAlgebra", "DynamicHMC", "Turing" 
]
     
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


# BLAS.get_config() # check which one is being used
# using OpenBLAS_jll  # should MKL fail
 


project_directory = joinpath( homedir(), "bio.data", "bio.snowcrab", "output" ) 

#   push!(LOAD_PATH, project_directory)  # add the directory to the load path, so it can be found
#   include( "startup.jl" )
#   # include( joinpath( project_directory, "startup.jl" ))    # alt
  
#   fn_env = joinpath( project_directory, "size_structured_dde_environment.jl" )  
#   include(  fn_env )
  

# colorscheme!("Monokai24bit") # for REPL
 
# theme(:default)  # defaults for graphics
#theme(:vibrant)
#theme(:bright)

# to start a graphics window
# gr(size=(1000,1000),legend=false,markerstrokewidth=0,markersize=4)
# gr()
 
    print( "\n\nThe following RData warnings Conversion of RData.RExtPtrcan be safely ignored. \n\n" )

    fndat  = joinpath( project_directory, "size_structure", "base_data_julia.rdata" )
    
    o = load( fndat, convert=true)

    # sex codes
    # male = 0
    # female = 1
    # sex.unknown = 2

    # # maturity codes
    # immature = 0
    # mature = 1
    # mat.unknown = 2

    Y  = o["out"]["Y"]
    Y.logcw = log.(Y.cw)
    Y.logsa = log.(Y.sa)

    # discretize time (weekly):
    sigdigits=3
    Y[!,:ti] = Y[!,:year] .+ round.(trunc.(Y[!,:julian] ./ 365. .* 52 ) ./ 52, digits=sigdigits) # weekly 

    # dimensionality of problem
    nb = o["out"]["nb"]
    aus = o["out"]["au"]
    yrs = 1999:2022
    weeks = 1:52
    sexes = ["0", "1"]
    mats = ["0", "1"]

    o = nothing

# Define a kernel mixture with n_nodes gaussian components and Estimate weights
 
    # maximum(Y.logcw) # 5.14166355650266
    # minimum(Y.logcw) # 2.138889000323256
    xr = (minimum(Y.logcw)*0.95, maximum(Y.logcw)*1.05 ) # (2.031944550307093, 5.398746734327793)
    bw = 0.05   # on log cw scale .. 0.2 seems ok .. lower reduces spurious echoes
    sd_modes = bw/5 # approx magnitude of variability of size classes
    nodes =  range(xr[1], stop=xr[2], step=bw)   # nodes
    n_nodes = length(nodes)  # 34
    ti_window = [-4, 4]  # weeks to include in moving window

#   sampler = Turing.NUTS(0.65)  # surprisingly slow s
    sampler = Turing.SMC()  # fast
    iterations = 1000
    n_chains =4 
     
    nmin = 30
  
    outdir = joinpath( project_directory, "size_structure", "posteriors_summaries" )
    mkpath(outdir) 

 

@model function kmm(x, nodes, n_nodes, sd_modes, N  )

    # Kernel Mixture model with nodes as pre-specified components
    # that cover the space from min_x to max_x
    # \alpha (alpha) = A concentration parameter of 1 results in all sets of probabilities being equally likely, i.e., in this case the Dirichlet distribution of dimension k is equivalent to a uniform distribution over a k-1-dimensional simplex.  
    
    sigmasq ~ filldist( InverseGamma(2, 3), n_nodes)  # variance prior 
    # inital guess of approx mean location of modal groups
     
    kernels = map( i -> truncated( Normal(nodes[i], sqrt(sigmasq[i])), lower=nodes[1], upper=nodes[n_nodes] ), 1:n_nodes )
    alpha ~ Dirichlet(n_nodes, 1.0)
    mixdist = MixtureModel(kernels, alpha)
    x ~ filldist(mixdist, N)
 end
 

showall(x) = show(stdout, "text/plain", x)

function sample_posteriors(Y; sexes="0", mats="0", aus="7", yrs=2009, weeks=32, 
    ti_window=[-1,1], sampler=Turing.NUTS(0.65), iterations=1000, n_chains=1, nmin=100, 
    outdir=joinpath( project_directory, "size_structure", "posteriors_summaries" ), overwrite=false )   

    if false
        # debug: 
        sex="0"; mat="0"; au="103"; yr=1999; week=14
    end

    out = nothing

    for yr in yrs
        fn = string( "posterior_summaries_", yr, ".csv" )
        fndat  = joinpath( outdir, fn )
        sample_the_data = true
        if isfile(fndat) 
            if !overwrite
                sample_the_data = false
            end
        end

        if sample_the_data

            for week in weeks
                mti = yr .+ round.( (week .+ ti_window) ./ 52, digits=3) # weekly   week + ti_window
                kt = findall( x -> x >= mti[1] && x <= mti[2], skipmissing(Y.ti) )
              if length(kt) >= nmin
            for au in aus
                ka = intersect( kt, findall( x -> x==au , skipmissing(Y.space_id) ) )
              if length(ka) >= nmin
            for sex in sexes
                ks = intersect( ka, findall( x -> x==sex, skipmissing(Y.sex) ) )
              if length(ks) >= nmin
            for mat in mats
                n =  intersect( ks, findall( x -> x==mat, skipmissing(Y.mat) ) )
                N =  length(n)
                tout = string("| sex: ", sex, "| mat: ", mat, "| au: ", au, "|year: ", yr, "| week: ", week, "| N: ", N ) 
                if N >= nmin 
                    @show tout
                    uu = kde( Y.logcw[n], bandwidth=bw, boundary=xr, npoints=1024, weights=1 ./ Y.sa[n] )
                    Neff = sum(1 ./ Y.sa[n])
                    vv = sample( uu.x, ProbabilityWeights(uu.density), 1024 )
                    model = kmm( vv, nodes, n_nodes, sd_modes, 1024)
                    if n_chains==1
                        chain = sample(model, sampler, iterations)
                    else
                        chain = sample(model, sampler, MCMCThreads(), iterations, n_chains)
                    end

                        if false
                            showall(chain)
                            v = fill(3.8387, n_nodes)
                            v = fill(4.7, n_nodes)
                            mp = predict( kmm(missing, v, n_nodes, sd_modes, 1024 ), chain)
                            k = mp.value.data[1:1000,:,:]
                            plot(legend=false, xlim=(nodes[[1,n_nodes]]) )
                            for i in 1:N 
                                density!(k[:,i,:], color="blue")
                            end
                            density!(Y.logcw[n])
                        end

                    res = DataFrame( summarystats(chain) )
                    res[!,:sex]  .= sex
                    res[!,:mat]  .= mat
                    res[!,:au]   .= au
                    res[!,:year] .= yr
                    res[!,:week] .= week
                    res[!,:N]    .= Neff
                    out = vcat(out, res)
                    res = nothing
                end
            end # mat
              end # sex if
            end # sex
              end # au if
            end  # au
              end # time if
            end   # time  
            
            CSV.write(fndat, res, writeheader=true)
            @show fndat 
     
        end
    end
end
    

sample_posteriors(Y, sexes=sexes, mats=mats, aus=aus, yrs=yrs, weeks=weeks, 
    ti_window=ti_window, sampler=sampler, iterations=iterations, n_chains=n_chains, nmin=nmin, 
    outdir=outdir, overwrite=true ) 



```


Now back in R for some additional analysis

```R

year.assessment = 2022
p = bio.snowcrab::load.environment( year.assessment=year.assessment )

require(hdf5r)

fname = file.path(p$project.outputdir, "size_structure", "posterior_summaries.hdf5")

fo = H5File$new(fname, mode = "r")
H5A(fo)
a = fo$attr_open("H5File")
fo$get_info()
a$attr_name()
a$get_space()
a$get_type()
a$get_storage_size()
a$read()
a

```
