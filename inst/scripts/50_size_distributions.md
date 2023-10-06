 

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

    pg = areal_units( 
        p = snowcrab_parameters(
            project_class="carstm",
            yrs=1999:year.assessment,   
            areal_units_type="tesselation",
            carstm_model_label=  paste( "1999_present", "fb", sep="_" )
    ))
    dim(pg)

   
if (0){
        data.cube.installed = try( require(data.cube), silent=TRUE )
        if (!data.cube.installed) {        
            # also stored on github
            install.packages("data.cube", repos = paste0("https://", c(
                "jangorecki.gitlab.io/data.cube",
                "cloud.r-project.org"
            )))
        }
          
        require(data.cube)

}


# ---------------------
# method 0: base data 

    # create base data of individual size records

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


# ---------------------
# method 5 .. Normalization via weighted kernel density

    # defaults:
    np = 512  # # discretizations in fft
    xrange =c(10, 160)
    xr = round( log(xrange), digits=2 ) 
    dx = diff(xr)/(np-1)  # 0.00544
    xvals = seq( xr[1], xr[2], by=dx )

    # bw is on log scale ... approx (log) SD for each interval  ~ data_resolution is ~1 to 2 mm (observation error)
    # bw = 0.1 # ~ 20 dx ~ overly smooth
    bw = 0.025  # optimal for sparse data
    # bw = 0.02 # 4 dx is too noisy for arithmetic but good for geometric means
    # bw = 0.015 # 2.8 dx is too noisy for arithmetic but good for geometric means
    # bw = 0.0165 # 3.0 dx is too noisy for arithmetic but good for geometric means
    # bw = 0.01 # 2 dx is too noisy for arithmetic but good for geometric means
    Y = c(1996:2023)
    xrange =c(10, 160)

    # sa-weighted kernel density by sid , sex, mat (with au and quarter)
    size_distributions(p=p, toget="kernel_density_weighted", pg=pg, bw=bw, np=np, xrange=xrange, redo=TRUE )  # 1996:present
    # size_distributions(p=p, toget="kernel_density_weighted", pg=pg, bw=bw, np=np, xrange =xrange, redo=TRUE, Y=Y )  # incremental updates
 
    # extract  
    O = size_distributions(p=p, toget="kernel_density_weighted", bw=bw, np=np, xrange =xrange, Y=Y ) #subsets
 


# ---------------------
# Modal analysis

# iteratively go to every polygon and identify samples from it and surrounding areas and times (time window is 5 weeks :: +/- 2 weeks)
# and compute modes and troughs
    
    Y = c(1996:2023) 
    # bw is on log scale ... approx (log) SD for each interval  ~ data_resolution is ~1 to 2 mm (observation error)
    bw = 0.025 # 4 dx is too noisy for arithmetic but good for geometric means
    np = 512  # # discretizations in fft
    xrange = c(10, 160)
    zlevels = c(0, 100)  # left bounds
    tlevels = c(-2, 6)  # left bounds
    sigdigits = 3

    xr = round( log(xrange), digits=2 ) 
    dx = diff(xr)/(np-1)  # 0.00544
    xvals = seq( xr[1], xr[2], by=dx )

    M = size_distributions(p=p, toget="kernel_density_weighted", bw=bw, np=np, xrange =xrange, Y=Y ) #subsets


    kdb = aggregate_by( M, 
        agg_by = c("year", "sex", "mat", "region", "zi", "ti" ), 
        xvals= xvals,
        recale_density_to_numerical_density=TRUE, 
        use_geometric=TRUE, 
        add_offset=TRUE 
    )
 
 
    res = extract_modes_by_factor( kdb,   
        sexes=c("0", "1"), 
        mats=c("0", "1"), 
        regions=c("cfanorth", "cfasouth", "cfa4x"),
        years=as.character(1996:2023),
        zlevels=zlevels, 
        tlevels=tlevels,
        X=xvals, 
        lowpassfilter=0.0001, lowpassfilter2=0.0001, 
        sigdigits=3, plot=TRUE )



    fn = file.path( p$annual.results, "size_distributions_summary.RDS" )
    saveRDS( res, file=fn )
    res = readRDS(fn)


    
    vn = "peaks"

    plot(density( unlist(res[[vn]][ s=="0" & m=="0", ..vn]), bw=0.05 ) )
 
    plot(density( unlist(res[[vn]][ s=="0" & m=="0" & t=="-2", ..vn]), bw=0.05 ), col="blue" )
    lines(density( unlist(res[[vn]][ s=="0" & m=="0" & t=="6", ..vn]), bw=0.05 ), col="red" )
 
    lines(density( unlist(res[[vn]][ s=="0" & m=="0" , ..vn]), bw=0.05 ), col="blue" )
    lines(density( unlist(res[[vn]][ s=="0" & m=="0" & z=="0", ..vn]), bw=0.05 ), col="red" )
   
    plot(density( unlist(res[[vn]][ s=="0"  & r=="cfa4x", ..vn]), bw=0.05 ), col="red" )
    lines(density( unlist(res[[vn]][ s=="0"  & r=="cfasouth", ..vn]), bw=0.05 ), col="green" )
    lines(density( unlist(res[[vn]][ s=="0"  & r=="cfanorth", ..vn]), bw=0.05 ), col="blue" )
    
    mi = identify_modes( 
        Z = unlist(res[["peaks"]][ s=="0" & m=="0" , peaks]),  
        lowpassfilter2=0.0001,
        bw=0.05, sigdigits=3, plot=TRUE) 
    
    mm = identify_modes( 
        Z = unlist(res[["peaks"]][ s=="0" & m=="1" , peaks]),  
        lowpassfilter2=0.0001,
        bw=0.05, sigdigits=3, plot=TRUE) 

    
    fi = identify_modes( 
        Z = unlist(res[["peaks"]][ s=="1" & m=="0" , peaks]),  
        lowpassfilter2=0.0001,
        bw=0.05, sigdigits=3, plot=TRUE) 
    
    fm = identify_modes( 
        Z = unlist(res[["peaks"]][ s=="1" & m=="1" , peaks]),  
        lowpassfilter2=0.0001,
        bw=0.05, sigdigits=3, plot=TRUE) 

    mds = rbind( 
        data.table( cw = exp(fm$peaks),   cat="fm"),
        data.table( cw = exp(fi$peaks),   cat="fi"),
        data.table( cw = exp(mm$peaks),   cat="mm"),
        data.table( cw = exp(mi$peaks),   cat="mi")
    )

    growth = function(sex, instar) {
        if (sex=="f") cw = exp(2.198848 + 0.315026 * (instar - 4) )
        if (sex=="m") cw = exp(1.917564 + 0.298914 * (instar - 3) )
        return(cw)
    }


    cwbr = data.table( instar=1:14, ml=growth(sex="m", instar=1:14 ), fl=growth(sex="f", instar=1:14) )

    females = rbind(
        data.table( 
            cw = exp(fi$peaks),  
            cat="fimm",  
            instar=as.numeric( as.character( cut( exp(fi$peaks), breaks= cwbr$fl, labels=cwbr$instar[1:13] ) ) ) ),
        data.table( 
            cw = exp(fm$peaks),  
            cat="fmat",  
            instar=as.numeric( as.character( cut( exp(fm$peaks), breaks= cwbr$fl, labels=cwbr$instar[1:13] ) ) ) )
    )
 

females$cw = log(females$cw)
females = females[order(cw)]
nf = nrow(females)

of = lm( females$cw[2:nf]~ females$cw[1:(nf-1)] )

summary(of)
 
Residuals:
       1        2        3        4        5        6 
-0.01972  0.00171 -0.00765  0.02150  0.06637 -0.06221 

Coefficients:
                       Estimate Std. Error t value Pr(>|t|)
(Intercept)              0.4660     0.1180    3.95    0.017
females$cw[1:(nf - 1)]   0.9447     0.0369   25.59  1.4e-05

Residual standard error: 0.0479 on 4 degrees of freedom
Multiple R-squared:  0.994,	Adjusted R-squared:  0.992 
F-statistic:  655 on 1 and 4 DF,  p-value: 1.39e-05
 
    males = rbind(
        data.table( 
            cw = exp(mi$peaks),  
            cat="mimm",  
            instar=as.numeric( as.character( cut( exp(mi$peaks), breaks= cwbr$ml, labels=cwbr$instar[1:13] ) ) ) ),
        data.table( 
            cw = exp(mm$peaks),  
            cat="mmat",  
            instar=as.numeric( as.character( cut( exp(mm$peaks), breaks= cwbr$ml, labels=cwbr$instar[1:13] ) ) ) )
    )

males$cw = log(males$cw)
males = males[order(cw)]
nm = nrow( males )

om = lm( males$cw[2:nf]~ males$cw[1:(nf-1)] )

summary(om)

Residuals:
        1         2         3         4         5         6 
-3.76e-05 -1.25e-02 -8.14e-03  2.34e-04  6.47e-02 -4.42e-02 

Coefficients:
                     Estimate Std. Error t value Pr(>|t|)
(Intercept)            0.3443     0.1068    3.22    0.032
males$cw[1:(nf - 1)]   0.9887     0.0306   32.33  5.5e-06

Residual standard error: 0.0399 on 4 degrees of freedom
Multiple R-squared:  0.996,	Adjusted R-squared:  0.995 
F-statistic: 1.05e+03 on 1 and 4 DF,  p-value: 5.46e-06
 
plot(0, 0, xlim=log(c(10,120)), ylim=log(c(10,120)), type="n" )

lines( females$cw[2:nf] ~ females$cw[1:(nf-1)], col="orange" )
points( females$cw[2:nf] ~ females$cw[1:(nf-1)], col="orange", pch=24 )

lines( males$cw[2:nm] ~ males$cw[1:(nm-1)] , col="green" )
points( males$cw[2:nm] ~ males$cw[1:(nm-1)] , col="green", pch=22 )

ol = loess( females$cw[2:nf] ~ females$cw[1:(nf-1)] )
op = predict(ol)
lines( op ~ females$cw[1:(nf-1)]   )

ol = loess(  males$cw[2:nm] ~ males$cw[1:(nm-1)] )
op = predict(ol)
lines( op ~ males$cw[1:(nm-1)]  )
 
odf = diff( exp(females$cw)) / exp(females$cw[-nrow(females)])
odm = diff( exp(males$cw)) / exp( males$cw[-nrow(males)] )

mean(odf)  # 34 % increase each moult
mean(odm)  # 33 % increase each moult

plot(cw ~ instar, females, type="n")
points(cw ~ instar, females[cat=="fimm"], col="brown", cex=2, pch=20)
points(cw ~ instar, females[cat=="fmat"], col="blue", cex=2, pch=20)
 
plot(cw ~ instar, males, type="n")
points(cw ~ instar, males[cat=="mimm"], col="brown", cex=2, pch=20)
points(cw ~ instar, males[cat=="mmat"], col="blue", cex=2, pch=20)
 


# ---------------------
# size structure continuous form

    Y = c(1996:2023) 
    # bw is on log scale ... approx (log) SD for each interval  ~ data_resolution is ~1 to 2 mm (observation error)
    bw = 0.025 # 4 dx is too noisy for arithmetic but good for geometric means
    np = 512  # # discretizations in fft
    xrange = c(10, 160)

    xr = round( log(xrange), digits=2 ) 
    dx = diff(xr)/(np-1)  # 0.00544
    xvals = seq( xr[1], xr[2], by=dx )

    zlevels = c(0, 100)  # left bounds
    tlevels = c(-2, 6)  # left bounds

    # get data
    M = size_distributions(p=p, toget="kernel_density_weighted", bw=bw, np=np, Y=Y ) #subsets

    kdb = aggregate_by( M, 
        agg_by = c("year", "sex", "mat", "region", "zi", "ti" ), 
        xvals= xvals,
        recale_density_to_numerical_density=TRUE, 
        use_geometric=TRUE, 
        add_offset=TRUE 
    )

    # segregation by temp
    vn ="2022_0_1_cfanorth_100_-2"
    plot( kdb[[vn]] ~ logcw, kdb, type="l" )
    
    vn ="2022_0_1_cfanorth_100_6"
    lines( kdb[[vn]] ~ logcw, kdb  )
    
    vn ="2022_1_0_cfanorth_100_-2"
    lines( kdb[[vn]] ~ logcw, kdb  )
     
    vn ="2022_1_0_cfanorth_100_6"
    lines( kdb[[vn]] ~ logcw, kdb  )
     

    # segregation by depth
    vn ="2022_0_1_cfanorth_0_-2"
    plot( kdb[[vn]] ~ logcw, kdb, type="l" )
    
    vn ="2022_0_1_cfanorth_100_-2"
    lines( kdb[[vn]] ~ logcw, kdb  )
    
    vn ="2022_1_0_cfanorth_0_-2"
    lines( kdb[[vn]] ~ logcw, kdb  )
     
    vn ="2022_1_0_cfanorth_100_-2"
    lines( kdb[[vn]] ~ logcw, kdb  )
     


    kdb = aggregate_by( M, 
        agg_by = c( "yr", "sex", "mat", "region" ), 
        xvals= xvals,
        recale_density_to_numerical_density=TRUE, 
        use_geometric=TRUE, 
        add_offset=TRUE 
    )


    # annual male
    vn ="2022_0_0_cfa4x"
    plot( kdb[[vn]] ~ cw, kdb, type="l"     )
    
    vn ="2021_0_0_cfa4x"
    lines( kdb[[vn]] ~ cw, kdb , col="red" )
 
    vn ="2019_0_0_cfa4x"
    lines( kdb[[vn]] ~ cw, kdb , col="blue" )
  
    vn ="2018_0_0_cfa4x"
    lines( kdb[[vn]] ~ cw, kdb , col="green" )
 
 
    # annual female
    vn ="2022_1_0"
    plot( kdb[[vn]] ~ cw, kdb, type="l" )
    
    vn ="2021_1_0"
    lines( kdb[[vn]] ~ cw, kdb , col="red" )
 
    vn ="2019_1_0"
    lines( kdb[[vn]] ~ cw, kdb , col="blue" )
  
    vn ="2018_1_0"
    lines( kdb[[vn]] ~ cw, kdb , col="blue" )
 

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

  out = list( Y=Y, nb=nb, au=au )
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
