 

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
    MI = M[["ysm"]][["densities"]]
    MO = M[["ysm"]][["peaks"]]

    plot(density~cw, MI[sex=="0" & mat=="0" , ], pch="." )
    abline(v=MO[ sex=="0" & mat=="0", cw ], col="gray", lwd=0.5 )

    plot(density~cw, MI[sex=="0" & mat=="1" , ], pch=".")  # NOTE misses the largest size group
    abline(v=MO[ sex=="0" & mat=="1", cw ], col="gray" )

    plot(density~cw, MI[sex=="1" & mat=="0" , ], pch=".")
    abline(v=MO[ sex=="1" & mat=="0", cw ], col="gray" )

    plot(density~cw, MI[sex=="1" & mat=="1" , ], pch=".")
    abline(v=MO[ sex=="1" & mat=="1", cw ], col="gray" )
    
  
    # collect point estimates 

    fn = file.path(survey_size_freq_dir, "modes_male_imm.png" )
    png(filename=fn, width=1000,height=600, res=144)
    mi = identify_modes( Z = unlist(MO[ sex=="0" & mat=="0" , cw]),  
        lowpassfilter2=0.0001, xvals=xvals, dx=dx, bw=0.05, sigdigits=3, plot=TRUE) 
    abline(v=4, col="orange", lwd=2, lty="dashed") # likely a nonmode
    dev.off()

    fn = file.path(survey_size_freq_dir, "modes_male_mat.png" )
    png(filename=fn, width=1000,height=600, res=144)
    mm = identify_modes( Z = unlist(MO[ sex=="0" & mat=="1" , cw]),  
        lowpassfilter2=0.0001, xvals=xvals, dx=dx, bw=0.05, sigdigits=3, plot=TRUE) 
    dev.off()

    fn = file.path(survey_size_freq_dir, "modes_female_imm.png" )
    png(filename=fn, width=1000,height=600, res=144)
    fi = identify_modes( Z = unlist(MO[ sex=="1" & mat=="0" , cw]),  
        lowpassfilter2=0.0001, xvals=xvals, dx=dx, bw=0.05, sigdigits=3, plot=TRUE) 
    dev.off()
 
    fn = file.path(survey_size_freq_dir, "modes_female_mat.png" )
    png(filename=fn, width=1000,height=600, res=144)
    fm = identify_modes( Z = unlist(MO[ sex=="1" & mat=="1" , cw]),  
        lowpassfilter2=0.0001, xvals=xvals, dx=dx, bw=0.05, sigdigits=3, plot=TRUE) 
    dev.off()
  
    mds = rbind( 
        data.table( logcw = fi$peaks,   sex="f", mat= "i" ),
        data.table( logcw = fm$peaks,   sex="f", mat= "m" ),
        data.table( logcw = mi$peaks,   sex="m", mat= "i" ),
        data.table( logcw = mm$peaks,   sex="m", mat= "m" )
    )
    mds$cw = exp(mds$logcw)

    plot( mds$logcw[ mds$sex=="f"])
    plot( mds$logcw[ mds$sex=="m"])  # one incorrect mode just under 4 ,.. remove
    
    bad = mds[  logcw > 3.95 & logcw < 4.05 & sex=="m" & mat=="i", which=TRUE ]
    mds = mds[-bad,]

    f = mds[ sex=="f", ][order(mat, cw),]
    f$seq = 1:nrow(f)
    plot( cw ~ seq, f)
    i = 5:7  # hyp: imm just under corresponding mature size
    arrows(f$seq[i], f$cw[i], f$seq[i+3], f$cw[i+3], length=0.2, col= 1:3)
    i = f[ mat=="i", which=TRUE]
    i = i[-length(i)]
    arrows(f$seq[i], f$cw[i], f$seq[i+1], f$cw[i+1], length=0.2 )

    m = mds[ sex=="m", ][order(mat, cw),]
    m$seq = 1:nrow(m)
    plot( cw ~ seq, m)
    i = 6:8 # hyp: imm just under corresponding mature size
    arrows(m$seq[i], m$cw[i], m$seq[i+4], m$cw[i+4], length=0.2, col= 1:3)
    i = m[ mat=="i", which=TRUE]
    i = i[-length(i)]
    arrows(m$seq[i], m$cw[i], m$seq[i+1], m$cw[i+1], length=0.2 )
    # last immature group -> maturity is missing .. add it below

    # assign instar: imm patterns seems simple
    mds$instar = NA
    
    # female
    ii = mds[sex=="f" & mat=="i", which=TRUE]
    mds$instar[ii] = cw_to_instar( mds$logcw[ii], "f" ) 

    jj = mds[sex=="f" & mat=="m", which=TRUE]
    for (j in jj) {
        k = mds[sex=="f" & mat=="i" & logcw < mds$logcw[j] , which=TRUE]
        mds$instar[j] = max( mds$instar[ k] ) + 1
    }

    # verify:
    f = mds[ sex=="f", ][order(mat, cw),]
    plot( cw ~ instar, f)
    j = f[ mat=="m", which=TRUE]
    for (i in j ){
        k = f[mat=="i" & instar==(f$instar[i] -1), which=TRUE ] 
        arrows(f$instar[i], f$cw[i], f$instar[k], f$cw[k], length=0.2, code=1, col="red")
    }
    i = f[ mat=="i", which=TRUE]
    i = i[-length(i)]
    arrows(f$instar[i], f$cw[i], f$instar[i+1], f$cw[i+1], length=0.2 )


    # male
    ii = mds[sex=="m" & mat=="i", which=TRUE]
    mds$instar[ii] = cw_to_instar( mds$logcw[ii], "m" ) 

    jj = mds[sex=="m" & mat=="m", which=TRUE]
    for (j in jj) {
        k = mds[sex=="m" & mat=="i" & logcw < mds$logcw[j] , which=TRUE]
        mds$instar[j] = max( mds$instar[ k] ) + 1
    }

    # verify:
    m = mds[ sex=="m", ][order(mat, cw),]
    plot( cw ~ instar, m)
    j = m[ mat=="m", which=TRUE]
    for (i in j ){
        k = m[mat=="i" & instar==(m$instar[i] -1), which=TRUE ] 
        arrows(m$instar[i], m$cw[i], m$instar[k], m$cw[k], length=0.2, code=1, col="blue")
    }
    i = m[ mat=="i", which=TRUE]
    i = i[-length(i)]
    arrows(m$instar[i], m$cw[i], m$instar[i+1], m$cw[i+1], length=0.2 )

 
    
    mds$pred = NA
    mds$logcw0 = NA

    i = mds[sex=="f" & mat=="i", which=TRUE ] 
    for (j in i) {
        k = mds[sex=="f" & mat=="i" & instar==(mds[j, instar]-1), which=TRUE ]
        if (length(k) >0) mds$logcw0[j] = mds$logcw[k]
     }
    
    i = mds[sex=="f" & mat=="m", which=TRUE ] 
    for (j in i) {
        k = mds[sex=="f" & mat=="i" & instar==(mds[j, instar]-1), which=TRUE ]
        if (length(k) >0) mds$logcw0[j] = mds$logcw[k]
    }

    i = mds[sex=="m" & mat=="i", which=TRUE ] 
    for (j in i) {
        k = mds[sex=="m" & mat=="i" & instar==(mds[j, instar]-1), which=TRUE ]
        if (length(k) >0) mds$logcw0[j] = mds$logcw[k]
     }
    
    i = mds[sex=="m" & mat=="m", which=TRUE ] 
    for (j in i) {
        k = mds[sex=="m" & mat=="i" & instar==(mds[j, instar]-1), which=TRUE ]
        if (length(k) >0) mds$logcw0[j] = mds$logcw[k]
    }


    of = lm( logcw ~ logcw0 * mat,  mds[ sex=="f",], na.action="na.omit")
    mds$pred[ which(mds$sex=="f")] = predict( of, mds[ sex=="f",] )
    summary(of)
    
    plot(logcw~logcw0,  mds[ sex=="f",])
    points(logcw~logcw0,  mds[ sex=="f" & mat=="i",], col="red")
    points(pred~logcw0,  mds[ sex=="f" & mat=="i",], col="green")
    points(pred~logcw0,  mds[ sex=="f" & mat=="m",], col="purple")
  
    om = lm( logcw ~ logcw0 * mat,  mds[ sex=="m",], na.action="na.omit")
      

    mds$pred[ which(mds$sex=="m")] = predict( om, mds[ sex=="m",] )
    summary(om)
    
    plot(pred~logcw0,  mds[ sex=="m",])
    points(logcw~logcw0,  mds[ sex=="m" & mat=="m",], col="red", pch="+")
    points(logcw~logcw0,  mds[ sex=="m" & mat=="i",], col="red")
    points(pred~logcw0,  mds[ sex=="m" & mat=="i",], col="green")
    points(pred~logcw0,  mds[ sex=="m" & mat=="m",], col="purple", pch="+")
   

    # add unobserved instars: 1:4 and 13 Male
    oif = lm( logcw~ instar, mds[sex=="f" & mat=="i", ], na.action="na.omit")
    omf = lm( logcw~ instar, mds[sex=="f" & mat=="m", ], na.action="na.omit")
    summary(oif) # Adjusted R-squared:  0.999
    summary(omf) # Adjusted R-squared:  0.977

    oim = lm( logcw~ instar, mds[sex=="m" & mat=="i", ], na.action="na.omit")
    omm = lm( logcw~ instar, mds[sex=="m" & mat=="m", ], na.action="na.omit")
    summary(oim) # Adjusted R-squared:  0.999
    summary(omm) # Adjusted R-squared:  0.999

    unobs= CJ(logcw=NA, sex=c("m", "f"), mat="i", cw=NA, instar=1:3, pred=NA, logcw0=NA)

    mds = rbind(mds, unobs)
    # mature male instar not represented (due to rarity vs instar 12)
    # add instar 13
    logcw12 = mds[sex=="m" & mat=="i" & instar==12, logcw]
    instar13 = data.table( NA, "m", "m", NA, 13, NA, logcw12)
    names(instar13) = names(mds)
    mds = rbind(mds, instar13 )
    
    mds$predicted = NA
    mds$predicted_se = NA

    i = mds[sex=="f" & mat=="i", which=TRUE]
    ip = predict(oif, mds[i,], se.fit=TRUE )
    mds$predicted[i] =ip$fit
    mds$predicted_se[i] =ip$se.fit

    i = mds[sex=="f" & mat=="m", which=TRUE]
    ip = predict(omf, mds[i,], se.fit=TRUE )
    mds$predicted[i] =ip$fit
    mds$predicted_se[i] =ip$se.fit
   
    i = mds[sex=="m" & mat=="i", which=TRUE]
    ip = predict(oim, mds[i,], se.fit=TRUE )
    mds$predicted[i] =ip$fit
    mds$predicted_se[i] =ip$se.fit

    i = mds[sex=="m" & mat=="m", which=TRUE]
    ip = predict(omm, mds[i,], se.fit=TRUE )
    mds$predicted[i] =ip$fit
    mds$predicted_se[i] =ip$se.fit

# full predicted pattern:
    mds$cwmean = exp(mds$predicted)
    mds$cwlb = exp(mds$predicted - 1.96*mds$predicted_se )
    mds$cwub = exp(mds$predicted + 1.96*mds$predicted_se )
  
    fn = file.path(survey_size_freq_dir, "modes_growth_female.png" )
    png(filename=fn, width=1000, height=600, res=144)
        f = mds[ sex=="f", ][order(mat, instar),]
        plot( cwmean ~ instar, f, type="p" )
        j = f[ mat=="m", which=TRUE]
        for (i in j ){
            k = f[mat=="i" & instar==(f$instar[i] -1), which=TRUE ] 
            arrows(f$instar[i], f$cwmean[i], f$instar[k], f$cwmean[k], length=0.2, code=1, col="red")
        }
        i = f[ mat=="i", which=TRUE]
        i = i[-length(i)]
        arrows(f$instar[i], f$cwmean[i], f$instar[i+1], f$cwmean[i+1], length=0.2 )
    dev.off()
 
    fn = file.path(survey_size_freq_dir, "modes_growth_male.png" )
    png(filename=fn, width=1000,height=600, res=144)
        m = mds[ sex=="m", ][order(mat, instar),]
        plot( cwmean ~ instar, m, type="p" )
        j = m[ mat=="m", which=TRUE]
        for (i in j ){
            k = m[mat=="i" & instar==(m$instar[i] -1), which=TRUE ] 
            arrows(m$instar[i], m$cwmean[i], m$instar[k], m$cwmean[k], length=0.2, code=1, col="blue")
        }
        i = m[ mat=="i", which=TRUE]
        i = i[-length(i)]
        arrows(m$instar[i], m$cwmean[i], m$instar[i+1], m$cwmean[i+1], length=0.2 )
 
    dev.off()


```

    Theses are the results:

    Female growth of modes (t vs t-1)
    
lm(formula = logcw ~ logcw0 * mat, data = mds[sex == "f", ], 
    na.action = "na.omit")

Residuals:
       2        3        4        5        6        7        8        9 
-0.01187  0.02915 -0.00847 -0.01420 -0.00507  0.01045 -0.03180  0.06122 
      10 
-0.02942 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept)   0.2101     0.0932    2.25    0.074
logcw0        1.0326     0.0293   35.28  3.4e-07
matm         -0.9082     0.3247   -2.80    0.038
logcw0:matm   0.1848     0.0843    2.19    0.080

Residual standard error: 0.0375 on 5 degrees of freedom
  (1 observation deleted due to missingness)
Multiple R-squared:  0.998,	Adjusted R-squared:  0.996 
F-statistic:  680 on 3 and 5 DF,  p-value: 6.02e-07


 
    Male growth of modes  (t vs t-1)

lm(formula = logcw ~ logcw0 * mat, data = mds[sex == "m", ], 
    na.action = "na.omit")

Residuals:
     Min       1Q   Median       3Q      Max 
-0.04790 -0.01198  0.00367  0.01192  0.03543 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept)   0.1954     0.0489    4.00   0.0052
logcw0        1.0302     0.0140   73.72  2.2e-11
matm         -0.2583     0.2521   -1.02   0.3398
logcw0:matm   0.0219     0.0606    0.36   0.7280

Residual standard error: 0.0272 on 7 degrees of freedom
  (1 observation deleted due to missingness)
Multiple R-squared:  0.999,	Adjusted R-squared:  0.999 
F-statistic: 2.28e+03 on 3 and 7 DF,  p-value: 7.9e-11


 summary(omf) # Adjusted R-squared:  0.977

Call:
lm(formula = logcw ~ instar, data = mds[sex == "f" & mat == "m", 
    ], na.action = "na.omit")

Residuals:
      1       2       3 
-0.0253  0.0507 -0.0253 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept)  -0.0117     0.4402   -0.03    0.983
instar        0.4090     0.0439    9.32    0.068

Residual standard error: 0.0621 on 1 degrees of freedom
Multiple R-squared:  0.989,	Adjusted R-squared:  0.977 
F-statistic: 86.9 on 1 and 1 DF,  p-value: 0.068

  summary(oif) # Adjusted R-squared:  0.999

Call:
lm(formula = logcw ~ instar, data = mds[sex == "f" & mat == "i", 
    ], na.action = "na.omit")

Residuals:
       1        2        3        4        5        6        7 
 0.02286 -0.01271  0.00171 -0.01086 -0.01943 -0.00900  0.02743 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept)  1.11986    0.02762    40.5  1.7e-07
instar       0.31157    0.00379    82.1  5.1e-09

Residual standard error: 0.0201 on 5 degrees of freedom
Multiple R-squared:  0.999,	Adjusted R-squared:  0.999 
F-statistic: 6.74e+03 on 1 and 5 DF,  p-value: 5.08e-09


summary(oim) # Adjusted R-squared:  0.999

Call:
lm(formula = logcw ~ instar, data = mds[sex == "m" & mat == "i", 
    ], na.action = "na.omit")

Residuals:
     Min       1Q   Median       3Q      Max 
-0.02999 -0.02159 -0.00889  0.01581  0.06051 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept)  1.16009    0.03428    33.9  5.1e-09
instar       0.30310    0.00408    74.3  2.1e-11

Residual standard error: 0.0316 on 7 degrees of freedom
Multiple R-squared:  0.999,	Adjusted R-squared:  0.999 
F-statistic: 5.53e+03 on 1 and 7 DF,  p-value: 2.1e-11


summary(omm) # Adjusted R-squared:  0.977

Call:
lm(formula = logcw ~ instar, data = mds[sex == "m" & mat == "m", 
    ], na.action = "na.omit")

Residuals:
       1        2        3 
 0.00367 -0.00733  0.00367 


 Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept)  0.56633    0.07005    8.08    0.078
instar       0.34300    0.00635   54.01    0.012

Residual standard error: 0.00898 on 1 degrees of freedom
  (1 observation deleted due to missingness)
Multiple R-squared:     1,	Adjusted R-squared:  0.999 
F-statistic: 2.92e+03 on 1 and 1 DF,  p-value: 0.0118

    These results are entered directly into their respective Julia functions.

    And some more plots:

```R




    mds$diff = exp(mds$logcw) - exp(mds$logcw0) 
    mds$diff_prop = mds$diff / exp(mds$logcw0)  # fractional increase

    fn = file.path(survey_size_freq_dir, "modes_growth_increment.png" )
    png(filename=fn, width=1000,height=600, res=144)
        plot(diff_prop ~ instar, mds, type="n")
        points(diff_prop ~ instar, mds[sex=="f" & mat=="i",], col="orange", pch=19 )
        points(diff_prop ~ instar, mds[sex=="f" & mat=="m",], col="darkred", pch=21 )
        points(diff_prop ~ instar, mds[sex=="m" & mat=="i",], col="green", pch=19 )
        points(diff_prop ~ instar, mds[sex=="m" & mat=="m",], col="darkblue", pch=21 )
    dev.off()

    mean(mds$diff_prop[ mds$sex=="f"], na.rm=TRUE)  # 30 % increase each moult
    mean(mds$diff_prop[ mds$sex=="m"], na.rm=TRUE)  # 30 % increase each moult

    mean(mds$diff_prop[ mds$sex=="f" & mds$mat=="i" ], na.rm=TRUE)  # 36.7 % increase each moult
    mean(mds$diff_prop[ mds$sex=="m" & mds$mat=="i" ], na.rm=TRUE)  # 34.9 % increase each moult

    mean(mds$diff_prop[ mds$sex=="f" & mds$mat=="m" ], na.rm=TRUE)  # 17 % increase each moult
    mean(mds$diff_prop[ mds$sex=="m" & mds$mat=="m" ], na.rm=TRUE)  # 16 % increase each moult



    plot(diff ~ instar, mds, type="n")
    points(diff ~ instar, mds[sex=="f" & mat=="i",], col="orange", pch=19 )
    points(diff ~ instar, mds[sex=="f" & mat=="m",], col="darkred", pch=21 )
    points(diff ~ instar, mds[sex=="m" & mat=="i",], col="green", pch=19 )
    points(diff ~ instar, mds[sex=="m" & mat=="m",], col="darkblue", pch=21 )
 
    mds$diffp = exp(mds$pred) - exp(mds$logcw0) 
    
    points(diffp ~ instar, mds[sex=="f" & mat=="i",], col="orange", pch=23 )
    points(diffp ~ instar, mds[sex=="f" & mat=="m",], col="darkred", pch=24 )
    points(diffp ~ instar, mds[sex=="m" & mat=="i",], col="green", pch=23 )
    points(diffp ~ instar, mds[sex=="m" & mat=="m",], col="darkblue", pch=24 )

    # save as rdata for use in julia
    fn = file.path(p$project.outputdir, "size_structure", "modes_summary.rdata")
    save(mds, file=fn)
    # load(fn)

```
 


## Julia functions for Modal analysis with Kernel Mixture Modelling  ...
 
Just a sketch of mode fit technique .. still a bit too slow to use

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

Bring data into Julia:

```julia
 
pkgs = [
    "Revise", 
    # "MKL", # "OpenBLAS32",
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

using OpenBLAS_jll  # should MKL fail
import LinearAlgebra, OpenBLAS32_jll
LinearAlgebra.BLAS.lbt_forward(OpenBLAS32_jll.libopenblas_path)




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
 
    print( "\n\nThe following RData warnings Conversion of RData can be safely ignored. \n\n" )

    fndat  = joinpath( project_directory, "size_structure", "base_data_julia.rdata" )
    o = load( fndat, convert=true)
 
    fnmds = joinpath( project_directory, "size_structure", "modes_summary.rdata")
    mds = load( fnmds, convert=true)
    mds = mds["mds"]

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
 
    outdir = joinpath( project_directory, "size_structure", "posteriors_summaries" )
    mkpath(outdir) 

 
## from projects/model.abm/julia/MixtureModels.md

using Distributions, Turing
using StatsPlots

# Unidimensional Kernel Mixture model with K pre-specified components
# that cover the space from min_x to max_x
# \alpha (α) = A concentration parameter of 1 (or k, the dimension of the Dirichlet distribution, by the definition used in the topic modelling literature) results in all sets of probabilities being equally likely, i.e., in this case the Dirichlet distribution of dimension k is equivalent to a uniform distribution over a k-1-dimensional simplex. This is not the same as what happens when the concentration parameter tends towards infinity. In the former case, all resulting distributions are equally likely (the distribution over distributions is uniform). In the latter case, only near-uniform distributions are likely (the distribution over distributions is highly peaked around the uniform distribution). Meanwhile, in the limit as the concentration parameter tends towards zero, only distributions with nearly all mass concentrated on one of their components are likely (the distribution over distributions is highly peaked around the k possible Dirac delta distributions centered on one of the components, or in terms of the k-dimensional simplex, is highly peaked at corners of the simplex).

 
showall(x) = show(stdout, "text/plain", x)


@model function kmm(x, nodes, n_nodes, N, sd_nodes )
    # Kernel Mixture model with nodes as pre-specified components
    # alpha = concentration parameter of results in all sets of probabilities being equally likely, i.e., in this case the Dirichlet distribution of dimension k is equivalent to a uniform distribution over a k-1-dimensional simplex.  
    sigmasq ~ filldist( truncated(InverseGamma(2, 3), 0.0, sd_nodes), n_nodes)  # variance prior 
    kernels = map( i -> Normal(nodes[i], sqrt(sigmasq[i]) ), 1:n_nodes ) 
    alpha ~ Dirichlet(n_nodes, 1.0)
    mixdist = MixtureModel(kernels, alpha)
    x ~ filldist(mixdist, N)
end
 

# Define a kernel mixture with n_nodes gaussian components and Estimate weights
 

fi = Y[ (Y.sex.=="0") .& (Y.mat.=="0") .& (round.(Y.year).==2022) , :logcw ]  # imm females
fm = Y[ (Y.sex.=="0") .& (Y.mat.=="1"), :logcw ]  # mat females
mi = Y[ (Y.sex.=="1") .& (Y.mat.=="0"), :logcw ]  # imm males
mm = Y[ (Y.sex.=="1") .& (Y.mat.=="1"), :logcw ]  # mat males

fin = mds[ (mds.sex.=="f") .& (mds.mat.=="i"), :predicted ] # imm females
fmn = mds[ (mds.sex.=="f") .& (mds.mat.=="m"), :predicted ] # mat females
min = mds[ (mds.sex.=="m") .& (mds.mat.=="i"), :predicted ] # imm males
mmn = mds[ (mds.sex.=="m") .& (mds.mat.=="m"), :predicted ] # mat males


# choose: 
data = fi
nodes = fin


# maximum(Y.logcw) # 5.14166355650266
# minimum(Y.logcw) # 2.138889000323256
xr = round.( (minimum(Y.logcw), maximum(Y.logcw)), digits=2 ) # (2.031944550307093, 5.398746734327793)
     
# nodes = range(xr[1], stop=xr[2], step=0.2 )   # we have about 14 classes max .. for size range about 10 *3
nodes = nodes[(nodes .>=xr[1]) .& (nodes .<=xr[2]) ]
 
xrez = round( maximum(diff(nodes)), digits=2 )  # assume ~ 3SD so 1sd = xrez/3

# ti_window = [-4, 4]  # weeks to include in moving window

#   sampler = Turing.NUTS(0.65)  # surprisingly slow s
    sampler = Turing.SMC()  # fast
    iterations = 1000
    n_chains =4 
     
    nmin = 30
  

N = size(data, 1)

n_nodes = length(nodes)  

sd_nodes = xrez / 3.0 # approx magnitude of variability of size classes see :predicted_se

model = kmm(data, nodes, n_nodes, N, sd_nodes )

iterations = 1000
chain = sample(model, sampler, iterations)

alpha = [ chain[1, Symbol("alpha[$k]"), l] for k in 1:n_nodes]

v = fill(nodes[4], n_nodes)

mp = predict( kmm(missing, v, n_nodes,  N, xr[1], xr[2]), chain)
 
k = mp.value.data[1:1000,:,:]

plot()
for i in 1:N 
  density!(k[:,i,:], color="red")
end
density!(data)
 

function kmm_samples( chain, modes, N ) 
    nch = size(chain)
    out = Any[];
    for a in 1:N
        g = rand(1:nch[1])  # sim
        l = rand(1:nch[3])  # chain
        sigmasq = [ chain[g, Symbol("sigmasq[$k]"), l] for k in 1:n_nodes]
        kernels = map( i -> Normal(nodes[i], sqrt(sigmasq[i]) ), 1:n_nodes ) 
        alpha = [ chain[g, Symbol("alpha[$k]"), l] for k in 1:n_nodes]
        mixdist = MixtureModel(kernels, alpha)
        push!(out, rand(mixdist ) )
    end  
    return out  
end

o = kmm_samples(chain, modes, 5000)
hist(o, bins=50)


function kmm_samples_separated( chain, modes, N ) 
    nch = size(chain)
    out = Any[];
    for a in 1:N
        g = rand(1:nch[1])  # sim
        l = rand(1:nch[3])  # chain
        sigmasq = [ chain[g, Symbol("sigmasq[$k]"), l] for k in 1:n_nodes]
        alpha = [ chain[g, Symbol("alpha[$k]"), l] for k in 1:n_nodes]
        sep = Any[]
        for b in 1:n_nodes
            kernels = map( i -> Normal(nodes[i], sqrt(sigmasq[i]) ), 1:n_nodes ) 
            alph = zeros(n_nodes) 
            alph[b] = 1
            mixdist = MixtureModel(kernels, alph)
            push!(sep, rand(mixdist ) )
        end
        push!(out, sep)
    end  
    return reduce(hcat, out)'  
end

o = kmm_samples_separated(chain, modes, 5000)


f = hist(o[:,1], xlim=(0,5) )
f = hist!(f, o[:,2], bins=50)
f = hist!(f, o[:,3], bins=50)



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
    

#normalize_size_structure(Y, sexes=sexes, mats=mats, aus=aus, yrs=yrs, weeks=weeks, 
#    bw=bw, xr=xr, ti_window=ti_window, nmin=5, outdir=outdir, overwrite=true ) 





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




# END
