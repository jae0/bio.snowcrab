 

# estimate size structure
#   1. direct computation
#   2. modelled (factorial) == 1
#   3. modelled (space-time) using carstm
 
# ---------------------
# define parameters:

    year.assessment = 2022
    p = bio.snowcrab::load.environment( year.assessment=year.assessment )
        
    require(ggplot2)
    require(data.table)
 
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
  
    years = as.character(2012:2022)
    regions=c("cfanorth", "cfasouth", "cfa4x")

    outdir=file.path( p$annual.results, "figures", "size.freq", "survey")
    outdir = tempdir()

    ftype = "png"

    yvar ="den"
    yvar ="denl"
    
    plot_sex = "female"
    fn = file.path( outdir, paste( plot_sex, ftype, sep=".") )
    xlim= c(5, 85)
    
    plot_sex = "male"
    fn = file.path( outdir, paste( plot_sex, ftype, sep=".") )
    xlim= c(5, 140)

    plot_histogram_carapace_width( M=M, years=years, regions=regions, plot_sex=plot_sex, yvar=yvar, fn=fn, xlim=xlim )



    if (0) {
        ggplot( M[ region==regions[a] & year==years[b] ,], aes(cwd, den, fill=mat, colour=mat) ) +
            #geom_ribbon(aes(ymin=density_lb, max=density_ub), alpha=0.2, colour=NA) +
            geom_bar(stat = "identity") +
            labs(x="cw", y="density", size = rel(1.5)) +
            # scale_y_continuous( limits=c(0, 300) )  
            theme_light( base_size = 22 ) 
    }

 
    #  if (plot_sex =="male") mtext("MALE", side=3, outer=T, line=4, cex=1.4)
    #  if (plot_sex =="female") mtext("MALE", side=3, outer=T, line=4, cex=1.4)

    # require data.cube  ## as a data.cube 
    # M = size_distributions(p=p, toget="data_cubes", Y="simple_direct" )


 
# ---------------------
# method 2: simple linear (gaussian) model 

    M = size_distributions(p=p, toget="linear_model" )
 
    O = M[ region=="cfanorth" ,]
 
    ggplot( O[year== 2017  ,], aes(cwd, den, fill=mat, colour=mat) ) +
        # geom_ribbon(aes(ymin=density_lb, max=density_ub), alpha=0.2, colour=NA) +
        # geom_line() +
        geom_bar(stat = "identity") +
        labs(x="cw", y="density", size = rel(1.5)) +
        # scale_y_continuous( limits=c(0, 300) )  
        theme_light( base_size = 22 ) 



# ---------------------
# method 3: poisson model  via bigglm
    
    # too slow to complete
    M = size_distributions(p=p, toget="poisson_model" )
 


# ---------------------
# method 4: poisson via inla
 
 
    # adjust based upon RAM requirements and ncores
    require(INLA)
    inla.setOption(num.threads= floor( parallel::detectCores() / 2) )
    
    M = size_distributions(p=p, toget="poisson_inla" )
 
  
  
  
