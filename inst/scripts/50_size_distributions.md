 

# Estimate size structure

A number of methods: 

- Direct via counting (sums)
- Direct via areal density (arithmetic and geometric)
- Kernel density areal density (arithmetic and geometric)
- Modelled solutions: mostly too large of a problem to use


 ```R
    # Get data and format based upon parameters:

    year.assessment = 2022
    p = bio.snowcrab::load.environment( year.assessment=year.assessment )
    loadfunctions( "bio.snowcrab")
    loadfunctions( "aegis")

    require(ggplot2)
    require(data.table)
 
    plotdir=file.path( p$annual.results, "figures", "size.freq", "survey")
    
    years = as.character(1996: year.assessment)
    regions=c("cfanorth", "cfasouth", "cfa4x")

    # use generic fb polygons:
    pg = areal_units( 
        p = snowcrab_parameters(
            project_class="carstm",
            yrs=1999:year.assessment,   
            areal_units_type="tesselation",
            carstm_model_label=  paste( "1999_present", "fb", sep="_" )
    ))
    dim(pg)
 
    xrange = c(10, 160)  # size range (CW)

    dx = 2 # width of carapace with discretization to produce "cwd"

    # merge set and individual level data
    M = size_distributions(p=p, toget="base_data", pg=pg, xrange=xrange, dx=dx, redo=TRUE)


```

## Simple sums

Now that we have data ready, we can do some simple tabulations using data.tables (for speed):

```R
    # tabulate... non-zero counts ... must use add_zeros=TRUE to add them, on the fly
    M = size_distributions(p=p, toget="tabulated_data", xrange=xrange, dx=dx, redo=TRUE)
```


## Areal densities 

Areal densities are simply computed as well, but they need to make sure zero-valued results are included. Direct arithmetic and geometric means are simple. But to account for environmental covariates, a model-based approach is more flexible. Unfortunately the models below do not work due to large problem size and corresponding RAM/CPU bottlenecks. 

```R
    # Method 1 ..  equivalent to full factorial model without intercept.
    # directly compute areal densities (from above tabulations) 
    M = size_distributions(p=p, toget="simple_direct", xrange=xrange, dx=dx, Y=years, redo=TRUE)


    # take subset in years
    M = size_distributions(p=p, toget="simple_direct", xrange=xrange, dx=dx, Y=years )

    plotdir=file.path( p$annual.results, "figures", "size.freq", "survey", "direct")

    years_ss = as.character( c(-11:0) + year.assessment )
    plot_histogram_carapace_width( M=M, years=years_ss, regions=regions, 
        plot_sex="female", 
        yvar="den",  # den=arithmetic mean density, denl = geometric mean density  
        outdir=plotdir 
    ) 
    
    plot_histogram_carapace_width( M=M, years=years_ss, regions=regions, 
        plot_sex="female", 
        yvar="denl",  # den=arithmetic mean density, denl = geometric mean density  
        outdir=plotdir 
    )

    plot_histogram_carapace_width( M=M, years=years_ss, regions=regions, 
        plot_sex="male", 
        yvar="den",   # den=arithmetic mean density, denl = geometric mean density  
        outdir=plotdir 
    )

    plot_histogram_carapace_width( M=M, years=years_ss, regions=regions, 
        plot_sex="male", 
        yvar="denl",   # den=arithmetic mean density, denl = geometric mean density  
        outdir=plotdir 
    )
 

    if (0) {
        ss = M[ region=="cfanorth" & year== 2017, which=TRUE]

        ggplot( M[ ss ,], aes(cwd, den, fill=mat, colour=sex) ) +
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

```

So the above modelling attempts do not work. But the trick is to find an approach that will.


## Kernel-density based methods

Continuing with size analysis, we can use kernel density estimates of specific components: sex, maturity, year, time of year (season in quarters), region and set (sid).
First we construct kernel density estimates using a bandwidth of 0.025 units on a logarithmic scale. This corresponds to about 4 dx, where dx is the increment width of discretization (512 units in the xrange).  

These results are normalized by swept area to provide a density per unit area (1 km^2). 


```R
    # Normalization via weighted kernel density
   
    # loadfunctions( "bio.snowcrab")

    # key defaults that define kernal densities:
    np = 512  # # discretizations in fft
    xrange =c(10, 160)
   
    xr = round( log(xrange), digits=2 ) 
    dx = diff(xr)/(np-1)  # 0.00544
    xvals = seq( xr[1], xr[2], by=dx )

    # bw is on log scale ... approx (log) SD for each interval  ~ data_resolution is ~1 to 2 mm (observation error)
    # bw = 0.1 # ~ 20 dx ~ overly smooth
    bw = 0.025  # optimal for sparse data  <<<<<< DEFAULT for histograms >>>>>>
    # bw = 0.02 # 4 dx is too noisy for arithmetic but good for geometric means
    # bw = 0.015 # 2.8 dx is too noisy for arithmetic but good for geometric means
    # bw = 0.0165 # 3.0 dx is too noisy for arithmetic but good for geometric means
    # bw = 0.01 # 2 dx is too noisy for arithmetic but good for geometric means
    
    # years of interest:
    years = as.character( c(1996:year.assessment) )

    # sa-weighted kernel density by sid , sex, mat (with au and quarter)
    M = size_distributions(p=p, toget="kernel_density_weighted", bw=bw, np=np, xrange=xrange, Y=years,  redo=TRUE )   

     
``` 

### Size structure in continuous form:

Using these KD estimates, we can flexibly aggregate across any strata ... 

Aggregation by "agg_function" is arithmetic mean by default. Example below over-rides with a geometric mean after adding a small positive valued offset (smallest non-zero-value) that mimics the scale (magnitude) of observation errors. Here is amounts to 100 individuals /km^2, below which density is not defined. 

```R
    # get data
    M = size_distributions(p=p, toget="kernel_density_weighted", bw=bw, np=np, Y=years ) #subsets

    geometric_mean = function(x) {exp(mean( log(x), na.rm=TRUE) ) }

    K = aggregate_by( M, 
        agg_by = c("year", "sex", "mat", "region", "zi", "ti" ), 
        xvals= xvals,
        recale_density_to_numerical_density=TRUE, 
        agg_function = geometric_mean, 
        add_offset=TRUE 
    )

    # segregation by temp
    vn ="2022_0_1_cfanorth_100_-2"
    plot( K[[vn]] ~ logcw, K, type="l" )
    
    vn ="2022_0_1_cfanorth_100_6"
    lines( K[[vn]] ~ logcw, K  )
    
    vn ="2022_1_0_cfanorth_100_-2"
    lines( K[[vn]] ~ logcw, K  )
     
    vn ="2022_1_0_cfanorth_100_6"
    lines( K[[vn]] ~ logcw, K  )
     

    # segregation by depth
    vn ="2022_0_1_cfanorth_0_-2"
    plot( K[[vn]] ~ logcw, K, type="l" )
    
    vn ="2022_0_1_cfanorth_100_-2"
    lines( K[[vn]] ~ logcw, K  )
    
    vn ="2022_1_0_cfanorth_0_-2"
    lines( K[[vn]] ~ logcw, K  )
     
    vn ="2022_1_0_cfanorth_100_-2"
    lines( K[[vn]] ~ logcw, K  )
     


    K = aggregate_by( M, 
        agg_by = c( "year", "sex", "mat", "region" ), 
        xvals= xvals,
        recale_density_to_numerical_density=TRUE, 
        agg_function = geometric_mean, 
        add_offset=TRUE 
    )


    # annual male
    vn ="2022_0_0_cfa4x"
    plot( K[[vn]] ~ cw, K, type="l"     )
    
    vn ="2021_0_0_cfa4x"
    lines( K[[vn]] ~ cw, K , col="red" )
 
    vn ="2019_0_0_cfa4x"
    lines( K[[vn]] ~ cw, K , col="blue" )
  
    vn ="2018_0_0_cfa4x"
    lines( K[[vn]] ~ cw, K , col="green" )
 
 
    # annual female
    vn ="2022_1_0_cfa4x"
    plot( K[[vn]] ~ cw, K, type="l" )
    
    vn ="2021_1_0_cfa4x"
    lines( K[[vn]] ~ cw, K , col="red" )
 
    vn ="2019_1_0_cfa4x"
    lines( K[[vn]] ~ cw, K , col="blue" )
  
    vn ="2018_1_0_cfa4x"
    lines( K[[vn]] ~ cw, K , col="blue" )
 
```



# Modal analysis ... growth stanzas

Using the kernel density approach, we can compute on the normalized densities to identify size modes. 

The algorithm is simply to go to every polygon and identify samples from it and potentially surrounding areas and times (potentially time window) and compute modes and troughs from first and second order differentials of the smoothed kernel densities. For production: we use a small local spatial group with a time window...

However, we can use other factors to stratify, such as depth (zlevels are left bounds in meters) and temperature (tlevels are left bounds of temperature in Celcius). Here simple, low/med/high are implemented, but not used as the above small-area-based methods work well. 


```R 
    # moving average in space and time

    # band width needs to be a little bit larger than used for size-structure  .. distributions must be smoother for modes to be clear 
    bw = 0.05  

    # time window to span 9 weeks about center
    ti_window=c(-4,4)  
    
    # spatial window is nearest-neighbours in spatial graph

    # Choose:

    redo =TRUE
    redo =FALSE
    M = size_distributions(p=p, toget="kernel_density_weighted", moving_average=TRUE,
        pg=pg, ti_window=ti_window, 
        bw=bw, np=np, xrange =xrange, Y=years, redo=redo ) #subsets

    # aggregate based upon strata:
    K = aggregate_by( M, 
        agg_by = c( "year", "au", "sex", "mat" ),  # strata
        xvals= xvals,
        recale_density_to_numerical_density=TRUE,  ### keep normalized to reduce scale issues
        agg_function = function(x) {exp(mean( log(x), na.rm=TRUE) ) }, # geometric_mean 
        add_offset=TRUE 
    )

    # plot a K-slice: annual male, immature
    vn ="1996_140_1_0";  plot( K[[vn]] ~ cw, K, type="l" )
    
    # extract summary stats
    res = size_structure_extract( K, 
        aus=pg$AUID, 
        sexes=c("0", "1"), 
        mats=c("0", "1"), 
        years=years,
        X=xvals, # note bw is already determine by input data .. "kernel_density_weighted"
        lowpassfilter=0.0005, lowpassfilter2=0.0005, 
        dx=dx, sigdigits=3, plot=TRUE )
  
    fn = file.path( p$annual.results, "size_distributions_summary.RDS" )
    saveRDS( res, file=fn )
    res = readRDS(fn)

 
    # and some plots:

    vn = "peaks" 
    plot(density( unlist(res[[vn]][ s=="0" & m=="0", ..vn]), bw=0.05 ) )
 
    plot(density( unlist(res[[vn]][ s=="0" & m=="0"  , ..vn]), bw=0.05 ), col="blue" )
    lines(density( unlist(res[[vn]][ s=="0" & m=="1"  , ..vn]), bw=0.05 ), col="red" )
    
    plot(density( unlist(res[[vn]][ s=="0"  & y=="2022", ..vn]), bw=0.05 ), col="red" )
    lines(density( unlist(res[[vn]][ s=="0"  & y=="2021", ..vn]), bw=0.05 ), col="green" )
    lines(density( unlist(res[[vn]][ s=="0"  & y=="2019", ..vn]), bw=0.05 ), col="blue" )
    
    mi = identify_modes( 
        Z = unlist(res[["peaks"]][ s=="0" & m=="0" , peaks]),  
        lowpassfilter2=0.0001,
        dx=dx, bw=0.04, sigdigits=3, plot=TRUE) 
    
    # looks like variability is rather large and oversmoothed at bw=0.04 .. reduce bandwidth
    # looks like it misses oe at ~4.5
    mm = identify_modes( 
        Z = unlist(res[["peaks"]][ s=="0" & m=="1" , peaks]),  
        W = unlist(res[["peak_values"]][ s=="0" & m=="1" , peak_values]),  
        T = unlist(res[["troughs"]][ s=="0" & m=="1" , troughs]),  
        lowpassfilter2=0.001,
        dx=dx, bw=0.04, sigdigits=3, plot=TRUE) 

    
    fi = identify_modes( 
        Z = unlist(res[["peaks"]][ s=="1" & m=="0" , peaks]),  
        lowpassfilter2=0.0001,
        dx=dx, bw=0.04, sigdigits=3, plot=TRUE) 
    
    fm = identify_modes( 
        Z = unlist(res[["peaks"]][ s=="1" & m=="1" , peaks]),  
        lowpassfilter2=0.001,
        dx=dx, bw=0.04, sigdigits=3, plot=TRUE) 

    # collect them all 
    mds = rbind( 
        data.table( cw = exp(fm$peaks),   cat="fm"),
        data.table( cw = exp(fi$peaks),   cat="fi"),
        data.table( cw = exp(mm$peaks),   cat="mm"),
        data.table( cw = exp(mi$peaks),   cat="mi")
    )

    # growth tables from manual mode identification 
    growth = function(sex, instar) {
        if (sex=="f") cw = exp(2.198848 + 0.315026 * (instar - 4) )
        if (sex=="m") cw = exp(1.917564 + 0.298914 * (instar - 3) )
        return(cw)
    }


    cwbr = data.table( instar=1:14, ml=growth(sex="m", instar=1:14 ), fl=growth(sex="f", instar=1:14) )

    females = rbind(
        data.table( 
            logcw =  fi$peaks,  
            cat="fimm",  
            instar=as.numeric( as.character( cut( exp(fi$peaks), breaks= cwbr$fl, labels=cwbr$instar[1:13] ) ) ) ),
        data.table( 
            logcw = fm$peaks,  
            cat="fmat",  
            instar=as.numeric( as.character( cut( exp(fm$peaks), breaks= cwbr$fl, labels=cwbr$instar[1:13] ) ) ) )
    )
    


    females$cw = exp(females$logcw)
    females = females[order(cw)]
    nf = nrow(females)

    of = lm( females$logcw[2:nf] ~ females$logcw[1:(nf-1)] +  females$cat[2:nf] )

    summary(of)
    

    males = rbind(
        data.table( 
            logcw =  mi$peaks,  
            cat="mimm",  
            instar=as.numeric( as.character( cut( exp(mi$peaks), breaks= cwbr$ml, labels=cwbr$instar[1:13] ) ) ) ),
        data.table( 
            logcw = mm$peaks,  
            cat="mmat",  
            instar=as.numeric( as.character( cut( exp(mm$peaks), breaks= cwbr$ml, labels=cwbr$instar[1:13] ) ) ) )
    )

    males$cw = exp(males$logcw)
    males = males[order(cw)]
    nm = nrow( males )

    om = lm( males$logcw[2:nf]~ males$logcw[1:(nf-1)] + males$cat[2:nf]  )

    summary(om)

```

# Female growth of modes (t vs t-1)
 
Call:
lm(formula = females$logcw[2:nf] ~ females$logcw[1:(nf - 1)] + 
    females$cat[2:nf])

Residuals:
        1         2         3         4         5         6 
-5.42e-03  9.29e-03  8.06e-03 -2.35e-02  1.16e-02  2.82e-18 

Coefficients:
                          Estimate Std. Error t value Pr(>|t|)
(Intercept)                 0.1931     0.0545    3.54    0.038
females$logcw[1:(nf - 1)]   1.0372     0.0181   57.43  1.2e-05
females$cat[2:nf]fmat      -0.1310     0.0252   -5.21    0.014

Residual standard error: 0.017 on 3 degrees of freedom
Multiple R-squared:  0.999,	Adjusted R-squared:  0.999 
F-statistic: 2.66e+03 on 2 and 3 DF,  p-value: 1.34e-05


# Male growth of modes  (t vs t-1)
 
 
Call:
lm(formula = males$logcw[2:nf] ~ males$logcw[1:(nf - 1)] + males$cat[2:nf])

Residuals:
        1         2         3         4         5         6 
-2.39e-04  5.68e-03 -5.40e-03 -5.48e-03  5.44e-03 -1.14e-18 

Coefficients:
                        Estimate Std. Error t value Pr(>|t|)
(Intercept)              0.30153    0.02198    13.7  0.00084
males$logcw[1:(nf - 1)]  1.00027    0.00667   150.0  6.5e-07
males$cat[2:nf]mmat     -0.13564    0.00922   -14.7  0.00068

Residual standard error: 0.00635 on 3 degrees of freedom
Multiple R-squared:     1,	Adjusted R-squared:     1 
F-statistic: 1.74e+04 on 2 and 3 DF,  p-value: 8e-07


# And so more plots:

```R

    plot(0, 0, xlim=c(10,120), ylim=c(10,120), type="n", xlab="CW(t-1)", ylab="CW(t)" )

    lines( females$cw[2:nf] ~ females$cw[1:(nf-1)], col="orange" )
    points( females$cw[2:nf] ~ females$cw[1:(nf-1)], col="orange", pch=24 )

    lines( males$cw[2:nm] ~ males$cw[1:(nm-1)] , col="green" )
    points( males$cw[2:nm] ~ males$cw[1:(nm-1)] , col="green", pch=22 )
 
    odf = diff( females$cw) / females$cw[-nrow(females)]
    odm = diff( males$cw) /  males$cw[-nrow(males)] 

    mean(odf)  # 36 % increase each moult
    mean(odm)  # 27 % increase each moult

    plot(cw ~ instar, females, type="n")
    points(cw ~ instar, females[cat=="fimm"], col="brown", cex=2, pch=20)
    points(cw ~ instar, females[cat=="fmat"], col="blue", cex=2, pch=24)
    
    plot(cw ~ instar, males, type="n")
    points(cw ~ instar, males[cat=="mimm"], col="brown", cex=2, pch=20)
    points(cw ~ instar, males[cat=="mmat"], col="blue", cex=2, pch=24)
    abline(h=exp(4.5))  # location missed 

```


# END


```

--- NOT USED ... just a sketch of mode fit technique .. still a bit too slow to use

# Model analysis with Kernel Mixture Modelling in Julia ...

See projects/model.abm/julia/MixtureModels.md for more info.

First prep the data in R:

```R

  Y = size_distributions( p=p, toget="base_data" )
  nb = attributes(pg)$nb$nbs
  au =pg$AUID

  out = list( Y=years, nb=nb, au=au )
  fn = file.path( p$project.outputdir, "size_structure", "base_data_julia.rdata" )  

  save( out, file=fn )

```

tidyhouse

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
    # mature = 1s
    # mat.unknown = 2

    Y  = o["out"]["Y"]
    Y.logcw = log.(Y.cw)
    Y.logsa = log.(Y.sa)
    Y.wt = 1.0 ./ Y.sa

    sigdigits=3
    
    # data_resolution is ~1 to 2 mm (observations error)
    bw =  round(median(log.(22:142)-log.(20:140)), digits=2 )   # approx SD is interval

    # discretize time (weekly):s
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
    xr = round.( (minimum(Y.logcw), maximum(Y.logcw)), digits=2 ) # (2.031944550307093, 5.398746734327793)
     
    nodes =  range(xr[1], stop=xr[2], step=0.2 )   # we have about 14 classes max .. for size range about 10 *3
    n_nodes = length(nodes)  # 16
 
    xrez = round( maximum(diff(nodes)), digits=2 )
 
    ti_window = [-4, 4]  # weeks to include in moving window

#   sampler = Turing.NUTS(0.65)  # surprisingly slow s
    sampler = Turing.SMC()  # fast
    iterations = 1000
    n_chains =4 
     
    nmin = 30
  
    outdir = joinpath( project_directory, "size_structure", "posteriors_summaries" )
    mkpath(outdir) 

 
showall(x) = show(stdout, "text/plain", x)


@model function kmm(x, nodes, n_nodes, N )

    # Kernel Mixture model with nodes as pre-specified components
    # that cover the space from min_x to max_x
    # \alpha (alpha) = A concentration parameter of 1 results in all sets of probabilities being equally likely, i.e., in this case the Dirichlet distribution of dimension k is equivalent to a uniform distribution over a k-1-dimensional simplex.  
    
    sigmasq ~ filldist( truncated( InverseGamma(2, 3), lower=nodes[1], upper=nodes[n_nodes]), n_nodes)  # variance prior 
    # inital guess of approx mean location of modal groups
     
    kernels = map( i -> truncated( Normal(nodes[i], sqrt(sigmasq[i]) ), lower=nodes[1], upper=nodes[n_nodes] ), 1:n_nodes )
    # kernels = map( i -> Normal(nodes[i], sqrt(sigmasq[i]) ), 1:n_nodes )

    alpha ~ Dirichlet(n_nodes, 1.0)
    mixdist = MixtureModel(kernels, alpha)
    x ~ filldist(mixdist, N)
end
 


function normalize_size_structure(Y; sexes="0", mats="0", aus="7", yrs=2009, weeks=32, ti_window=[-1,1],
    bw=0.1,  nsmp=1024, np=256, xr=extrema(Y), nmin=5, 
    outdir=joinpath( project_directory, "size_structure", "posteriors_summaries" ), overwrite=false )   

    if false
        # debug: 
        sex="0"; mat="0"; au="103"; yr=1999; week=14
        sex="1"; mat="1"; au="304"; yr=1999; week=32
        sex="0"; mat="1"; au="318"; yr=1999; week=32
        
        nsmp=1024
    end

 
    for yr in yrs
        fn = string( "posterior_summaries_", yr, ".csv" )
        fnout  = joinpath( outdir, fn )
        sample_the_data = true
        if isfile(fnout) 
            if !overwrite
                sample_the_data = false
            end
        end

        if sample_the_data

            out1 = DataFrame( sex="-1", mat="-1", au="-1", year=-1, week=-1, Nsample=-1, Neffective=-1 )
            out2 = DataFrame( Tables.table(zeros(np)') )

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
                    uu = kde( Y.logcw[n], bandwidth=bw, boundary=xr, npoints=np, weights=Y.wt[n] )
                    # CairoMakie.lines!( uu , color="orange" )
                    push!( out1, [sex , mat, au, yr, week, N, round( sum( Y.wt[n]) ) ] )
                    push!( out2, uu.density' )
                    res = nothing
                end
            end # mat
              end # sex if
            end # sex
              end # au if
            end  # au
              end # time if
            end   # time  
                        
            CSV.write(fnout, hcat(out1, out2), writeheader=true)
            @show fnout 


        end
    end
end
    

normalize_size_structure(Y, sexes=sexes, mats=mats, aus=aus, yrs=yrs, weeks=weeks, 
    bw=bw, xr=xr, ti_window=ti_window, nmin=5, outdir=outdir, overwrite=true ) 





function sample_posteriors(Y; sexes="0", mats="0", aus="7", yrs=2009, weeks=32, 
    bw=0.1, xrez=bw*4.0, nsmp=1024, np=256, xr=extrema(Y),
    ti_window=[-1,1], sampler=Turing.NUTS(0.65), iterations=1000, n_chains=1, nmin=100, 
    outdir=joinpath( project_directory, "size_structure", "posteriors_summaries" ), overwrite=false )   

    if false
        # debug: 
        sex="0"; mat="0"; au="103"; yr=1999; week=14
        sex="1"; mat="1"; au="304"; yr=1999; week=32
        
        nsmp=1024
        xrez= round( maximum(diff(nodes)), digits=2 )
    end

    out = nothing

    for yr in yrs
        fn = string( "posterior_summaries_", yr, ".csv" )
        fnout  = joinpath( outdir, fn )
        sample_the_data = true
        if isfile(fnout) 
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
                    uu = kde( Y.logcw[n], bandwidth=bw, boundary=xr, npoints=np, weights=Y.wt[n] ) # base2 .. large enough to cover x .. fft
                    # CairoMakie.lines!( uu , color="orange" )
                         
                    Neff = round( sum( Y.wt[n]) )
                    vv = sample( uu.x, ProbabilityWeights(uu.density), nsmp )
                    model = kmm( vv, nodes, n_nodes, nsmp )
                    
                    chain = nothing
                    if n_chains==1
                        chain = sample(model, sampler, iterations)
                    else
                        chain = sample(model, sampler, MCMCThreads(), iterations, n_chains)
                    end

                        if false
                            showall(chain)
                            v = nodes
                            v = fill(4.14, n_nodes)
                            mp = predict( kmm(missing, v, n_nodes, nsmp  ), chain)
                            k = mp.value.data[1:1000,:,:]

                            # using CairoMakie

                            f = Figure()
                            ax = Axis(f[1, 1])

                            # xlim=(nodes[[1,n_nodes]]) 
                            for i in 1:N 
                                den = kde(vec(k[:,i,:]) )
                                CairoMakie.lines!( den, color="red" )
                            end
                            CairoMakie.lines!( kde(Y.logcw[n]), color="blue" )
                            f
 
                        end

                    if !isnothing(chain) 
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

                end
            end # mat
              end # sex if
            end # sex
              end # au if
            end  # au
              end # time if
            end   # time  
            
            CSV.write(fnout, out, writeheader=true)
            @show fnout 
            out = nothing

        end
    end
end
    

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
