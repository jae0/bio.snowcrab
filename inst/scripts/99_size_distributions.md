 

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

    require(ggplot2)
    require(data.table)
 
    years = as.character(p$yrs)
    regions=c("cfanorth", "cfasouth", "cfa4x")

    plotdir=file.path( p$annual.results, "figures", "size.freq", "survey")



# ---------------------
# method 0: base data 

    # create base data of individual size records
    Y = size_distributions(p=p, toget="base_data", redo=TRUE)

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
 
    pg = areal_units( 
        p = snowcrab_parameters(
            project_class="carstm",
            yrs=1999:year.assessment,   
            areal_units_type="tesselation",
            carstm_model_label=  paste( "1999_present", "fb", sep="_" )
    ))
    dim(pg)

    # time window is 11 weeks ~ 3 months 
    # a size bin must have 5 or more individuals and a total of 100 individuals in the whole sample
    
    # NOTE: for bw (bandwidth)
    # diff(range( log10(10:165) )) / 0.025 # = ~ 48 increments  
    # diff(range( log10(10:165) )) / 0.02  # = ~ 60 increments  
    # diff(range( log10(10:165) )) / 0.015  # = ~ 81 increments  
    # diff(range( log10(10:165) )) / 0.01 # 121 increments
    # there are < 13 modes expected so 48 increments or so should be enough to ID modes

    Y = size_distributions(p=p, toget="size_modes", pg=pg, grad_method="Richardson",
        bw=0.02, sigdigits=3, ti_window = c(-6, 6), kexp=2,
        n_min=30, n_cutoff=3, lowpassfilter=0.0001, lowpassfilter2=0.0001, 
        n_neighbours=1, plot_solutions=TRUE, redo=TRUE)
 
    Y = size_distributions(p=p, toget="size_modes", redo=FALSE )

    f  = size_structure_extract( Y$female  ) 
    m  = size_structure_extract( Y$male  ) ;
    mi = size_structure_extract( Y$imale  ) ;
    fi = size_structure_extract( Y$ifemale  ) ;
    
    # list_get = function(v, w) sapply( v, function(x) x[[w]] )
      
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
    

    u = f
    u = m
    u = mi
    u = fi

    v = identify_modes( list_get(u,  "peaks"), bw=0.02, sigdigits=3, 
        plot=TRUE, lowpassfilter=0.0001, lowpassfilter2=0.0001, grad_method="Richardson" )
    v = identify_modes( list_get(u,  "troughs"), bw=0.02, sigdigits=3, 
        plot=TRUE, lowpassfilter=0.0001, lowpassfilter2=0.0001, grad_method="Richardson" )
    v = identify_modes( list_get(u, "peak_values"), bw=0.02, sigdigits=3, 
        plot=TRUE, lowpassfilter=0.0001, lowpassfilter2=0.0001, grad_method="Richardson" )




    cwp = rbind( 
        data.table( cw = 10^(mm$peaks),   cat="m"),
        data.table( cw = 10^(mf$peaks),   cat="f"),
        data.table( cw = 10^(mmi$peaks),  cat="mi"),
        data.table( cw = 10^(mfi$peaks),  cat="fi")
    )

    cwt = rbind( 
        data.table( cw = 10^(mm$troughs), cat="m"), 
        data.table( cw = 10^(mf$troughs), cat="f"),
        data.table( cw = 10^(mmi$troughs), cat="mi"), 
        data.table( cw = 10^(mfi$troughs), cat="fi")
    )

    out = list( cwp=cwp, cwt=cwt ) 

    fn = file.path( p$annual.results, "size_distributions_inla_poisson.RDS" )
    saveRDS( out, file=fn )

    growth = function(sex, instar) {
        if (sex=="f") cw = exp(2.198848 + 0.315026 * (instar - 4) )
        if (sex=="m") cw = exp(1.917564 + 0.298914 * (instar - 3) )
        return(cw)
    }


    cwbr = data.table( instar=1:14, ml=growth(sex="m", instar=1:14 ), fl=growth(sex="f", instar=1:14) )

    females = rbind(
        data.table( 
            cw = 10^(mfi$peaks),  
            cat="fimm",  
            instar=as.numeric( as.character( cut( 10^(mfi$peaks), breaks= cwbr$fl, labels=cwbr$instar[1:13] ) ) ) ),
        data.table( 
            cw = 10^(mf$peaks),  
            cat="fmat",  
            instar=as.numeric( as.character( cut( 10^(mf$peaks), breaks= cwbr$fl, labels=cwbr$instar[1:13] ) ) ) )
    )


    males = rbind(
        data.table( 
            cw = 10^(mmi$peaks),  
            cat="mimm",  
            instar=as.numeric( as.character( cut( 10^(mmi$peaks), breaks= cwbr$ml, labels=cwbr$instar[1:13] ) ) ) ),
        data.table( 
            cw = 10^(mm$peaks),  
            cat="mmat",  
            instar=as.numeric( as.character( cut( 10^(mm$peaks), breaks= cwbr$ml, labels=cwbr$instar[1:13] ) ) ) )
    )

    

```


Now continue in Julia ...


```julia
     
    
```


