 

# Estimate size structure and related parameters (growth, maturity, etc.)

A number of methods to explore: 

- Direct via counting (sums)
- Direct via areal density (arithmetic and geometric)
- Kernel density areal density (arithmetic and geometric)
- Modelled solutions: mostly too large of a problem to use
- Size at maturity
- Modal size groups

 ```R
    # Get data and format based upon parameters:

    year.assessment = 2022
    p = bio.snowcrab::load.environment( year.assessment=year.assessment )
    loadfunctions( "aegis")
    loadfunctions( "bio.snowcrab")

    require(ggplot2)
    require(data.table)
 
    survey_size_freq_dir = file.path( p$annual.results, "figures", "size.freq", "survey")
    
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
 
    # sex codes
    # male = 0
    # female = 1
    # sex.unknown = 2

    # # maturity codes
    # immature = 0
    # mature = 1
    # mat.unknown = 2
 
    xrange = c(10, 150)  # size range (CW)

    dx = 2 # width of carapace with discretization to produce "cwd"

  

```

## Base data

```R
    # merge set and individual level data
    M = size_distributions(p=p, toget="base_data", pg=pg, xrange=xrange, dx=dx, redo=TRUE)
```

## Simple sums

Now that we have data ready, we can do some simple tabulations using data.tables (for speed):

```R
    # tabulate... non-zero counts ... must use add_zeros=TRUE to add them, on the fly
    M = size_distributions(p=p, toget="tabulated_data", xrange=xrange, dx=dx, redo=TRUE)
```


## Hot spots of Areal densities 

Areal densities are simply computed as well, but they need to make sure zero-valued results are included. Direct arithmetic and geometric means are simple. But to account for environmental covariates, a model-based approach is more flexible. Unfortunately the models below do not work due to large problem size and corresponding RAM/CPU bottlenecks. 

```R
    # directly from databases
    M = snowcrab.db( DS="set.complete", p=p ) # note depth is log transformed here
    setDT(M)

    # high density locations
    i = which(M$totno.all > 2*10^5)
    H = M[i, .(uid, plon, plat, towquality, dist, distance, surfacearea, vessel, yr, z, julian, no.male.all, no.female.all, cw.mean, totno.all, totno.male.imm, totno.male.mat, totno.female.imm, totno.female.mat, totno.female.primiparous, totno.female.multiparous, totno.female.berried)]

    H$log10density = log10(H$totno.all)
 
    library(ggplot2)

    cst = coastline_db( p=p, project_to=st_crs(pg) ) 
        
    isodepths = c(100, 200, 300)
    isob = isobath_db( DS="isobath", depths=isodepths, project_to=st_crs(pg))
    isob$level = as.factor( isob$level)
        
    plt = ggplot() +
        geom_sf( data=cst, show.legend=FALSE ) +
        geom_sf( data=isob, aes( alpha=0.1, fill=level), lwd=0.1, show.legend=FALSE) +
        geom_point(data=H, aes(x=plon, y=plat, colour=log10density), size=5) +
        coord_sf(xlim = c(270, 940 ), ylim = c(4780, 5200 )) +
        theme(legend.position=c(0.08, 0.8), 
            axis.title.x=element_blank(),
            axis.title.y=element_blank())

    fn = file.path( p$project.outputdir, "maps",  "map_highdensity_locations.png" )
    png(filename=fn, width=1000,height=600, res=144)
        (plt)
    dev.off()

   
    # Method 1 ..  equivalent to full factorial model without intercept.
    # directly compute areal densities (from above tabulations) 
    M = size_distributions(p=p, toget="simple_direct", xrange=xrange, dx=dx, Y=years, redo=TRUE)

    # take subset in years
    M = size_distributions(p=p, toget="simple_direct", xrange=xrange, dx=dx, Y=years )

    outdir=file.path( survey_size_freq_dir, "direct")

    years_ss = as.character( c(-11:0) + year.assessment )
    plot_histogram_carapace_width( M=M, years=years_ss, regions=regions, 
        plot_sex="female", 
        yvar="den",  # den=arithmetic mean density, denl = geometric mean density  
        outdir=outdir 
    ) 
    
    plot_histogram_carapace_width( M=M, years=years_ss, regions=regions, 
        plot_sex="female", 
        yvar="denl",  # den=arithmetic mean density, denl = geometric mean density  
        outdir=outdir 
    )

    plot_histogram_carapace_width( M=M, years=years_ss, regions=regions, 
        plot_sex="male", 
        yvar="den",   # den=arithmetic mean density, denl = geometric mean density  
        outdir=outdir 
    )

    plot_histogram_carapace_width( M=M, years=years_ss, regions=regions, 
        plot_sex="male", 
        yvar="denl",   # den=arithmetic mean density, denl = geometric mean density  
        outdir=outdir 
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
 
    outdir=file.path( survey_size_freq_dir, "poisson_glm")
    
    regions = "cfanorth"
    plot_histogram_carapace_width( M=O$P, years=years, regions=regions, 
        plot_sex="male", 
        yvar="N",  # den=arithmetic mean density, denl = geometric mean density  
        outdir=outdir 
    )
 

# ---------------------
# method 4: poisson via inla .. problem is too large to compute
     
    # adjust based upon RAM requirements and ncores
    require(INLA)
    inla.setOption(num.threads= floor( parallel::detectCores() / 2) )
    
    O = size_distributions(p=p, toget="poisson_inla" )
 
    outdir=file.path( survey_size_freq_dir, "poisson_inla")
    
    regions = "cfanorth"
    plot_histogram_carapace_width( M=O$P, years=years, regions=regions, 
        plot_sex="male", 
        yvar="N",  # den=arithmetic mean density, denl = geometric mean density  
        outdir=outdir 
    )

```

So the above modelling attempts do not work. The trick is to find an approach that will.


## Size structure in continuous form: Kernel-density based methods

Continuing with size analysis, we can use kernel density estimates of specific components: sex, maturity, year, time of year (season in quarters), region and set (sid).
First we construct kernel density estimates using a bandwidth of 0.025 units on a logarithmic scale. This corresponds to about 4 dx, where dx is the increment width of discretization (512 units in the xrange).  

These results are normalized by swept area to provide a density per unit area (1 km^2). 


```R
    # Normalization via weighted kernel density
   
    # loadfunctions( "bio.snowcrab")

    # key defaults that define kernal densities:
    np = 512  # # discretizations in fft
    xrange =c(10, 150)
   
    xr = round( log(xrange), digits=2 ) 
    ldx = diff(xr)/(np-1)  # 0.00544
    xvals = seq( xr[1], xr[2], by=ldx )

    # bw is on log scale ... approx (log) SD for each interval  ~ data_resolution is ~1 to 2 mm (observation error)
    # bw = 0.1 # ~ 20 ldx ~ overly smooth
    bw = 0.05  # used for modal analysis
    # bw = 0.025  # optimal for sparse data  <<<<<< DEFAULT for histograms >>>>>>
    # bw = 0.02 # 4 ldx is too noisy for arithmetic but good for geometric means
    # bw = 0.015 # 2.8 ldx is too noisy for arithmetic but good for geometric means
    # bw = 0.0165 # 3.0 ldx is too noisy for arithmetic but good for geometric means
    # bw = 0.01 # 2 ldx is too noisy for arithmetic but good for geometric means
    
    # years of interest:
    years = as.character( c(1996:year.assessment) )
    
    ti_window=c(-4,4)
    # ti_window=c(-1,1)
    sigdigits = 3

    redo =TRUE
    redo = FALSE

   
    # sa-weighted kernel density by sid , sex, mat (with au and quarter)
    M = size_distributions(p=p, toget="kernel_density_weighted", bw=bw, np=np, 
    ldx=ldx, xrange=xrange, 
    Y=years, strata="yasm", pg=pg, sigdigits=sigdigits, ti_window=ti_window, redo=redo )   

    # sa-weighted kernel density by auid + ti, zi
    M = size_distributions(p=p, toget="kernel_density_weighted", bw=bw, np=np, 
    ldx=ldx, xrange=xrange, 
    Y=years, strata="smryzt", sigdigits=sigdigits, ti_window=ti_window, redo=redo )   
    
    # to reload:
    # M = size_distributions(p=p, toget="kernel_density_weighted", strata=strata, Y=years, pg=pg, bw=bw, np=np, sigdigits=sigdigits )
    
    

``` 

## Collect modes

Using these KD estimates, we can flexibly aggregate across any strata ... 

Aggregation by "agg_function" is arithmetic mean by default. Example below over-rides with a geometric mean after adding a small positive valued offset (smallest non-zero-value) that mimics the scale (magnitude) of observation errors. Here is amounts to 100 individuals /km^2, below which density is not defined. 
    
Using the kernel density approach, we can compute on the normalized densities to identify size modes. 

The algorithm is simply to go to every polygon and identify samples from it and potentially surrounding areas and times (potentially time window) and compute modes and troughs from first and second order differentials of the smoothed kernel densities. For production: we use a small local spatial group with a time window...

However, we can use other factors to stratify, such as depth (zlevels are left bounds in meters) and temperature (tlevels are left bounds of temperature in Celcius). Here simple, low/med/high are implemented, but not used as the above small-area-based methods work well. 


```R 

    bw = 0.05  # A bit larger than size structure ... smooth it
    np = 512
 
    lowpassfilter=0.0001
    lowpassfilter2=0.0001
    np = 512  # # discretizations in fft
    xrange =c(10, 150)
   
    xr = round( log(xrange), digits=2 ) 
    ldx = diff(xr)/(np-1)  # 0.00544
    xvals = seq( xr[1], xr[2], by=ldx )
 
    years = as.character( c(1996:year.assessment) )    # years of interest:
    ti_window=c(-4,4)     # time window to span 9 weeks about center
    # ti_window=c(-1,1)
    sigdigits = 3

    strata = "yasm"    # moving average in space and time
     # strata = "smryzt"  # qith temp and depth strata too  ...
  
    M = size_distributions(p=p, toget="kernel_density_modes", strata=strata, bw=bw, np=np, 
        Y=years, pg=pg, sigdigits=sigdigits, n_min=1,
        lowpassfilter=lowpassfilter, lowpassfilter2=lowpassfilter2, 
        redo=TRUE )
    
    # to reload:
    # M = size_distributions(p=p, toget="kernel_density_modes", strata=strata, bw=bw, np=np, sigdigits=sigdigits )
    
    plot(density~cw, MI$densities[sex=="0" & mat=="0" , ], pch="." )
    abline(v=MI$peaks[ sex=="0" & mat=="0", cw ], col="gray", lwd=0.5 )

    plot(density~cw, MI$densities[sex=="0" & mat=="1" , ], pch=".")
    abline(v=MI$peaks[ sex=="0" & mat=="1", cw ], col="gray" )

    plot(density~cw, MI$densities[sex=="1" & mat=="0" , ], pch=".")
    abline(v=MI$peaks[ sex=="1" & mat=="0", cw ], col="gray" )

    plot(density~cw, MI$densities[sex=="1" & mat=="1" , ], pch=".")
    abline(v=MI$peaks[ sex=="1" & mat=="1", cw ], col="gray" )
     
    
    O = MI[["peaks"]]
    
    # bad= O[ sex=="1" & exp(cw)>85, which=TRUE]
    # if (length(bad)>0) O= O[-bad,]

    O$instar = NA
    f = O[ sex=="1", which=TRUE ]
    O$instar[f] = cw_to_instar(O$cw[f], "f" ) 
    
    m = O[ sex=="0", which=TRUE ]
    O$instar[m] = cw_to_instar(O$cw[m], "m" ) 


    # and some plots:
   

    # fully aggregated view of modes:
    
    fn = file.path(survey_size_freq_dir, "modes_male_imm.png" )
    png(filename=fn, width=1000,height=600, res=144)
    mi = identify_modes( 
        Z = unlist(M[["peaks"]][ s=="0" & m=="0" , peaks]),  
        #W = unlist(M[["peak_values"]][ s=="0" & m=="0" , peak_values]),  
        lowpassfilter2=0.0001,
        xvals=xvals, dx=dx, bw=0.04, sigdigits=3, plot=TRUE) 
    dev.off()
    
    # looks like variability is rather large and oversmoothed at bw=0.04 .. reduce bandwidth
    # looks like it misses oe at ~4.5
    fn = file.path(survey_size_freq_dir, "modes_male_mat.png" )
    png(filename=fn, width=1000,height=600, res=144)
    mm = identify_modes( 
        Z = unlist(M[["peaks"]][ s=="0" & m=="1" , peaks]),  
        # W = unlist(M[["peak_values"]][ s=="0" & m=="1" , peak_values]),  
        # T = unlist(M[["troughs"]][ s=="0" & m=="1" , troughs]),  
        lowpassfilter2=0.001,
        xvals=xvals, dx=dx, bw=0.04, sigdigits=3, plot=TRUE) 
    dev.off()

    fn = file.path(survey_size_freq_dir, "modes_female_imm.png" )
    png(filename=fn, width=1000,height=600, res=144)
    fi = identify_modes( 
        Z = unlist(M[["peaks"]][ s=="1" & m=="0" , peaks]),  
        lowpassfilter2=0.0001,
        xvals=xvals, dx=dx, bw=0.04, sigdigits=3, plot=TRUE) 
    dev.off()
 
    fn = file.path(survey_size_freq_dir, "modes_female_mat.png" )
    png(filename=fn, width=1000,height=600, res=144)
    fm = identify_modes( 
        Z = unlist(M[["peaks"]][ s=="1" & m=="1" , peaks]),  
        lowpassfilter2=0.001,
        xvals=xvals, dx=dx, bw=0.04, sigdigits=3, plot=TRUE) 
    dev.off()
  
    # collect them all 
    mds = rbind( 
        data.table( cw = exp(fm$peaks),   cat="fm"),
        data.table( cw = exp(fi$peaks),   cat="fi"),
        data.table( cw = exp(mm$peaks),   cat="mm"),
        data.table( cw = exp(mi$peaks),   cat="mi")
    )

    # females
    plot(exp(peaks)~ jitter(instar+0.4), O[s=="1" & m=="1",], pch=".", col="red" ) # mature
    points(exp(peaks)~jitter(instar), O[s=="1" & m=="0",], pch=".", col="green") # immature
    
    abline(h=exp(fm$peaks), lwd=2, col="red")
    abline(h=exp(fi$peaks), col="green")
 
    # males
    plot(exp(peaks)~ jitter(instar+0.4), O[s=="0" & m=="1",], pch=".", col="blue" ) # mature
    points(exp(peaks)~jitter(instar), O[s=="0" & m=="0",], pch=".", col="orange") # immature
    
    abline(h=exp(mm$peaks), lwd=2, col="blue")
    abline(h=exp(mi$peaks), col="orange")
 
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

    # imm -> imm 
    i = which( males$cat =="mimm")
    imms = males[i,]
    ni = length(i) - 1
    oim = lm( imms$logcw[2:ni]~ imms$logcw[1:(ni-1)]  )

    summary(oim)

    # imm -> mat 
    i = which( males$cat =="mmat")
    im_ma = males[i,]
    j = which( males$cat =="mimm" & males$instar %in% (im_ma$instar-1) )
    im_mat_ma =  males[j,]
    im_mat_ma$instar =  im_mat_ma$instar + 1
    im = im_mat_ma[ im_ma, on="instar"]

    ni = length(i) - 1
    oim = lm( imms$logcw[2:ni]~ imms$logcw[1:(ni-1)]  )

    summary(oim)

    oimma = lm( im$i.logcw ~ im$logcw  )

    summary(oimma)


```

    Theses are the results:

    Female growth of modes (t vs t-1)
    
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



    Male growth of modes  (t vs t-1)

    immature:

    lm(formula = imms$logcw[2:ni] ~ imms$logcw[1:(ni - 1)])

    Residuals:
            1         2         3         4         5 
    -0.000239  0.005681 -0.005401 -0.005481  0.005440 

    Coefficients:
                        Estimate Std. Error t value Pr(>|t|)
    (Intercept)             0.30153    0.02198    13.7  0.00084
    imms$logcw[1:(ni - 1)]  1.00027    0.00667   150.0  6.5e-07

    Residual standard error: 0.00635 on 3 degrees of freedom
    Multiple R-squared:     1,	Adjusted R-squared:     1 
    F-statistic: 2.25e+04 on 1 and 3 DF,  p-value: 6.54e-07
 
    imm -> mature



    These results are entered directly into their respective Julia functions.

    And some more plots:

```R

    fn = file.path(survey_size_freq_dir, "modes_growth_increment.png" )
    png(filename=fn, width=1000,height=600, res=144)
 
    plot(0, 0, xlim=c(10,120), ylim=c(10,120), type="n", xlab="CW(t-1)", ylab="CW(t)" )

    lines( females$cw[2:nf] ~ females$cw[1:(nf-1)], col="orange" )
    points( females$cw[2:nf] ~ females$cw[1:(nf-1)], col="orange", pch=24 )

    lines( males$cw[2:nm] ~ males$cw[1:(nm-1)] , col="green" )
    points( males$cw[2:nm] ~ males$cw[1:(nm-1)] , col="green", pch=22 )
    dev.off()


    odf = diff( females$cw) / females$cw[-nrow(females)]
    odm = diff( males$cw) /  males$cw[-nrow(males)] 

    mean(odf)  # 33.5 % increase each moult
    mean(odm)  # 29.1 % increase each moult

    fn = file.path(survey_size_freq_dir, "modes_growth_female.png" )
    png(filename=fn, width=1000,height=600, res=144)
    
        plot(cw ~ instar, females, type="n")
        points(cw ~ instar, females[cat=="fimm"], col="brown", cex=2, pch=20)
        points(cw ~ instar, females[cat=="fmat"], col="blue", cex=2, pch=24)
    dev.off()
 
    fn = file.path(survey_size_freq_dir, "modes_growth_male.png" )
    png(filename=fn, width=1000,height=600, res=144)
        
        plot(cw ~ instar, males, type="n")
        points(cw ~ instar, males[cat=="mimm"], col="brown", cex=2, pch=20)
        points(cw ~ instar, males[cat=="mmat"], col="blue", cex=2, pch=24)
    dev.off()



```




# END
