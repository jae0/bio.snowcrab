 

# estimate size structure
#   1. direct computation
#   2. modelled (factorial) == 1
#   3. modelled (space-time) using carstm
 
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

    pa = snowcrab_parameters(
        project_class="carstm",
        yrs=p$yrs,   
        areal_units_type="tesselation",
        carstm_model_label=  paste( "1999_present", "fb", sep="_" ),
        selection = list(
            type = "number",
            biologicals=list( spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 ) ),
            biologicals_using_snowcrab_filter_class="fb"
        )
    )

    pg = areal_units(p=pa)
    dim(pg)

    # time window is 11 weeks ~ 3 months 
    # a size bin must have 5 or more individuals and a total of 100 individuals in the whole sample
    Y = size_distributions(p=p, toget="size_modes", pg=pg, 
        bw=0.02, ti_window = c(-6, 6), n_min=30, lower_filter=3, num_levels=1, eps=0.0001,
        plot_solutions=TRUE, redo=TRUE)

    names(Y)

    i = extract_peaks( Y$immature, "peaks" ) ; hist(i, "fd")
    f = extract_peaks( Y$female, "peaks" ) ; hist(f, "fd")
    m = extract_peaks( Y$male, "peaks" ) ;hist(m, "fd") 

    mi = identify_modes( i, bw=0.01, plot=TRUE, eps=0 )
    mf = identify_modes( f, bw=0.01, plot=TRUE, eps=0 )
    mm = identify_modes( m, bw=0.01, plot=TRUE, eps=0 )

    mmm = 10^(mm$peaks)            
    fff = 10^(mf$peaks)
    iii = 10^(mi$peaks)

    cwmm = data.table( cw = 10^(mm$peaks),   cat="m")
    cwmt = data.table( cw = 10^(mm$troughs), cat="m")

    cwfm = data.table( cw = 10^(mf$peaks),   cat="f")
    cwft = data.table( cw = 10^(mf$troughs), cat="f")

    cwim = data.table( cw = 10^(mi$peaks),   cat="i")
    cwit = data.table( cw = 10^(mi$troughs), cat="i")

    cwp = rbind( cwmm, cwfm, cwim )
    cwt = rbind( cwmt, cwft, cwit )

require(ggplot2)

p1 = ggplot(cwp, aes(x = cw, y = cw, color = cat)) + geom_point()

