 

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
    sppoly = areal_units( 
        p = snowcrab_parameters(
            project_class="carstm",
            yrs=1999:year.assessment,   
            areal_units_type="tesselation",
            carstm_model_label=  paste( "1999_present", "fb", sep="_" )
    ))
    dim(sppoly)
 
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
    M = size_distributions(p=p, toget="base_data", pg=sppoly, xrange=xrange, dx=dx, redo=TRUE)
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

    cst = coastline_db( p=p, project_to=st_crs(sppoly) ) 
        
    isodepths = c(100, 200, 300)
    isob = isobath_db( DS="isobath", depths=isodepths, project_to=st_crs(sppoly))
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

Giving up for now... :(



## Size structure in continuous form: Kernel-density based methods

Continuing with size analysis, we can use kernel density estimates of specific components: sex, maturity, year, time of year (season in quarters), region and set (sid).
First we construct kernel density estimates using a bandwidth of 0.025 units on a logarithmic scale. This corresponds to about 4 dx, where dx is the increment width of discretization (512 units in the xrange).  

These results are normalized by swept area to provide a density per unit area (1 km^-2). 


```R

# Normalization via weighted kernel density

# loadfunctions( "bio.snowcrab")

# key defaults that define kernal densities:
np = 512  # # discretizations in fft
# xrange =c(10, 150)

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
Y=years, strata="yasm", pg=sppoly, sigdigits=sigdigits, ti_window=ti_window, redo=redo )   

# sa-weighted kernel density by auid + ti, zi
M = size_distributions(p=p, toget="kernel_density_weighted", bw=bw, np=np, 
ldx=ldx, xrange=xrange, 
Y=years, strata="smryzt", sigdigits=sigdigits, ti_window=ti_window, redo=redo )   

# to reload:
# M = size_distributions(p=p, toget="kernel_density_weighted", strata=strata, Y=years, pg=sppoly, bw=bw, np=np, sigdigits=sigdigits )

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
    Y=years, pg=sppoly, sigdigits=sigdigits, n_min=1,
    lowpassfilter=lowpassfilter, lowpassfilter2=lowpassfilter2, 
    redo=TRUE )

# from these modes, find most frequent peaks with models as attributes
mds = size_distributions(p=p, toget="modal_groups", strata=strata, bw=bw, np=np, ldx=ldx, 
    sigdigits=sigdigits, redo=TRUE )

```
  

These are the results:

Female growth of modes (t vs t-1)

```output

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

```
 

Male growth of modes  (t vs t-1)

```output

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

```

Female mature simple model: cw vs instar

```output

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

```

Female immature simple model: cw vs instar 
```output
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
```




Male mature simple model: cw vs instar 
```output
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

```
Male immature simple model: cw vs instar 
```output
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
```





# Estimation of abundance at modes 

Density and variability estimation via Modal Kernel Mixture Models (KMM) is done in Julia: See projects/model.mixtures/kmm_snowcrab.md for more info. NOTE: this approach is still too slow to use operationally at each set level -- but is viable annually. But that would prevent further covariate modelling.  
 

Here instead, to estimate areal density, we use a knife-edged cut at midpoints between modal groups. This is imperfect as large groups can bleed into adjacent smaller groups. But it is consistent and simple.  We will use this for modelling: either via

- carstm (TODO -- that is to update data source as modelling approach is already complete )
- stmv (TODO -- that is to update data source as modelling approach is already complete )
- abm (see: projects/model.abm/julia/ )



```R
 
# density by stage (with zeros)
M = size_distributions(p=p, toget="tabulated_data_by_stage", pg=sppoly, xrange=xrange, dx=dx, redo=TRUE ) 

M$stage = as.character( M$stage )  # to make character matches simple below

r = "cfanorth"
r = "cfasouth"
r = "cfa4x"

plot(log(density)~ as.factor(region), M[  year==2022 & grepl("m|m", stage, fixed=TRUE ),])
    plot(log(density)~as.factor(stage) , M[region==r & grepl("m|i", stage, fixed=TRUE ),])
    plot(log(density)~as.factor(stage) , M[region==r & grepl("m|m", stage, fixed=TRUE ),])
    plot(log(density)~as.factor(stage) , M[region==r & grepl("f|i", stage, fixed=TRUE ),])
    plot(log(density)~as.factor(stage) , M[region==r & grepl("f|m", stage, fixed=TRUE ),])
 

plot(log(density)~ as.factor(year), M[region=="cfasouth" & grepl("m|i|12", stage, fixed=TRUE ),])
    plot(log(density)~ as.factor(year), M[region==r & grepl("f|m", stage, fixed=TRUE ),])
    plot(log(density)~ as.factor(year), M[region==r & grepl("f|i", stage, fixed=TRUE ),])
    plot(log(density)~ as.factor(year), M[region==r & grepl("m|m", stage, fixed=TRUE ),])
    plot(log(density)~ as.factor(year), M[region==r & grepl("m|i", stage, fixed=TRUE ),])
 




```




# Probability of observation: -- do not use (not really useful at this point)


- Create weight model for snow crab by smooth(cw), mat, sex, space, time, spacetime, depth, etc.. 

- then used to create observation weight for each observation in det tables.. 

- to be used for integration / simulation

- The model is based upon a CAR stucture

```R

  #require(ggplot2)
 
  require(bio.snowcrab)    
  loadfunctions("bio.snowcrab")

  year.assessment = 2022
  yrs = 1999:year.assessment
  years = as.character(yrs)

  spec_bio = bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 )

  snowcrab_filter_class = "all"    
  runlabel= "observation_weights"
 
  # params for probability of observation
  p = snowcrab_parameters( 
    project_class="carstm", 
    yrs=yrs,  
    areal_units_type="tesselation", 
    carstm_model_label= runlabel,  
    selection = list(
      type = "presence_absence",
      biologicals=list( spec_bio=spec_bio ),
      biologicals_using_snowcrab_filter_class=snowcrab_filter_class
    )
  )

  xrange = c(10, 150)  # size range (CW)
  dx = 10 # width of carapace with discretization to produce "cwd"


  survey_size_freq_dir = file.path( p$annual.results, "figures", "size.freq", "survey")
    
  regions=c("cfanorth", "cfasouth", "cfa4x")

  # use generic fb polygons:
  sppoly = areal_units( 
      p = snowcrab_parameters(
          project_class="carstm",
          yrs=1999:year.assessment,   
          areal_units_type="tesselation",
          carstm_model_label=  paste( "1999_present", "fb", sep="_" )
  )) 
  dim(sppoly)

  crs_lonlat = st_crs(projection_proj4string("lonlat_wgs84"))
  if (is.null(sppoly)) sppoly = areal_units( p=p )  # will redo if not found
  sppoly = st_transform(sppoly, crs=crs_lonlat )
  sppoly$data_offset = sppoly$sa

  plot(sppoly["AUID"])
 
  # prep params for carstm:

  p$space_name = sppoly$AUID 
  p$space_id = 1: nrow(sppoly)  # must match M$space

  p$time_name = as.character(p$yrs)
  p$time_id =  1: p$ny

  # p$cyclic_name = as.character(p$cyclic_levels)
  # p$cyclic_id = 1: p$nw

  
  #------------
    pH = p
    
    pH$label ="snowcrab_observation_probability"
    pH$carstm_model_label = model_label

    pH$variabletomodel = "pa"  

    pH$family = "binomial"   # "binomial",  # "nbinomial", "betabinomial", "zeroinflatedbinomial0" , "zeroinflatedbinomial0"
 
    pH$formula = formula( 
        pa ~ 1 
        + f(time, model="ar1", hyper=H$ar1)    
    #   + f( cyclic, model="seasonal", scale.model=TRUE, season.length=10, hyper=H$iid  )
        + f( inla.group( t, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) 
        + f( inla.group( z, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) 
    #   + f( inla.group( log.substrate.grainsize, method="quantile", n=11 ), model="rw2", scale.model=TRUE, hyper=H$rw2) 
        + f( cwd, model="rw2", scale.model=TRUE, hyper=H$rw2) 
        + f( space, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE,  hyper=H$bym2 ) 
        + f( space_time, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, group=time_space,  hyper=H$bym2, control.group=list(model="ar1", hyper=H$ar1_group)) 
    )

    mds = size_distributions(p=p, toget="modal_groups" )
    mds = mds[ instar>=3, ]
    pH$stage_name = mds$stage
    pH$stage_id = as.factor(mds$stage)
    pH$stage_value = mds$predicted  # logcw(mm)
    mds = NULL

    # -------------------
    
    pN = p
    pN$label ="snowcrab_abundance"
    pN$carstm_model_label = model_label

    pN$variabletomodel = "N"  

    pN$family = "poisson"   

    pN$formula = formula( 
        N ~ 1 + offset(sa )  # CARSTM does log transformation so do not log transform
        + f(time, model="ar1", hyper=H$ar1)    
    #   + f( cyclic, model="seasonal", scale.model=TRUE, season.length=10, hyper=H$iid  )
        + f( inla.group( t, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) 
        + f( inla.group( z, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) 
    #   + f( inla.group( log.substrate.grainsize, method="quantile", n=11 ), model="rw2", scale.model=TRUE, hyper=H$rw2) 
        + f( cwd, model="rw2", scale.model=TRUE, hyper=H$rw2) 
        + f( space, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE,  hyper=H$bym2 ) 
        + f( space_time, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, group=time_space,  hyper=H$bym2, control.group=list(model="ar1", hyper=H$ar1_group)) 
    )


#### contiinue with Habitat
 
  # prediction surface
  redo = FALSE
  # redo = TRUE
  M = size_distributions(p=p, toget="modal_groups_carstm_inputs", pg=sppoly, xrange=xrange, dx=dx , redo=redo )
  ntrials = round(1/M$sa)
  M=NULL


  # model pa using all data
  res = NULL; gc()
  res = carstm_model( p=p, 
    data='size_distributions(p=p, toget="modal_groups_carstm_inputs")', 
    sppoly=sppoly, 
    # theta = c( 0.926, 1.743, -0.401, 0.705, -2.574, 1.408, 2.390, 3.459, 3.321, -2.138, 3.083, -1.014, 3.558, 2.703 ),
    Ntrials=ntrials,
    nposteriors=5000,
    posterior_simulations_to_retain=c( "summary", "random_spatial", "predictions"), 
    family = p$family, 
    verbose=TRUE,
    #redo_fit=FALSE, 
    # debug = "fit",
    # control.family=list(control.link=list(model="logit")),  # default for binomial .. no need to specify
    # control.inla = list( strategy="laplace", int.strategy="eb" ),
    num.threads="4:3"
  )


  ( res$summary) 
   
    vn = "time"
    ts =  res$random[[vn]]
    IDx = as.numeric(ts$ID)
    plot( mean ~ IDx, ts, type="b", lwd=1.5, xlab=vn)
    lines( quant0.025 ~ IDx, ts, col="gray", lty="dashed")
    lines( quant0.975 ~ IDx, ts, col="gray", lty="dashed")

 
    vn = "cwd"
    ts =  res$random[[vn]]
    IDx = as.numeric(ts$ID)
    plot( mean ~ IDx, ts, type="b", lwd=1.5, xlab=vn)
    lines( quant0.025 ~ IDx, ts, col="gray", lty="dashed")
    lines( quant0.975 ~ IDx, ts, col="gray", lty="dashed")


    # inverse probability
    vn = "cwd"
    ts =  res$random[[vn]]
    IDx = as.numeric(ts$ID)
    plot( 1-mean ~ IDx, ts, type="b", lwd=1.5, xlab=vn)
    lines( 1-quant0.025 ~ IDx, ts, col="gray", lty="dashed")
    lines( 1-quant0.975 ~ IDx, ts, col="gray", lty="dashed")


    vn = "inla.group(t, method = \"quantile\", n = 9)"
    ts =  res$random[[vn]]
    IDx = as.numeric(ts$ID)
    plot( mean ~ IDx, ts, type="b", lwd=1.5, xlab=vn)
    lines( quant0.025 ~ IDx, ts, col="gray", lty="dashed")
    lines( quant0.975 ~ IDx, ts, col="gray", lty="dashed")

 
    vn = "inla.group(z, method = \"quantile\", n = 9)"
     ts =  res$random[[vn]]
    IDx = as.numeric(ts$ID)
    plot( mean ~ IDx, ts, type="b", lwd=1.5, xlab=vn)
    lines( quant0.025 ~ IDx, ts, col="gray", lty="dashed")
    lines( quant0.975 ~ IDx, ts, col="gray", lty="dashed")

 
    # maps of some of the results
    additional_features = features_to_add( 
        p=p, 
        area_lines="cfa.regions",
        isobaths=c( 100, 200, 300, 400, 500  ), 
        coastline =  c("canada", "united states of america"), 
        xlim=c(-80,-40), 
        ylim=c(38, 60) 
    )
  
    # map all bottom temps:
    outputdir = file.path(p$data_root, "maps", p$carstm_model_label )
    if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )
  
  
    graphics.off()
  
    # e.g. management lines, etc
    outputdir = file.path( p$modeldir, p$carstm_model_label, "maps" )
    if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )
     
    fn_root = paste("Predicted_habitat_probability_persistent_spatial_effect", sep="_")
    outfilename = file.path( outputdir, paste(fn_root, "png", sep=".") )
    

    vn = c( "random", "space", "combined" ) 
    
    # toplot = carstm_results_unpack( res, vn )
    # brks = pretty(  quantile(toplot[,"mean"], probs=c(0,0.975), na.rm=TRUE )  )
    
    brks = pretty(c(0, 1))
    plt = carstm_map(  res=res, vn=vn, 
      breaks = brks,
      colors=rev(RColorBrewer::brewer.pal(5, "RdYlBu")),
      additional_features=additional_features,
      title= "persistent spatial effect" ,
      outfilename=outfilename
    )  
    plt
 
    brks = pretty(c(1, 9))
    for (y in res$time_name ){
      for ( u in res$cyclic_name  ){
        fn_root = paste( "Bottom temperature",  as.character(y), as.character(u), sep="-" )
        outfilename = file.path( outputdir, paste( gsub(" ", "-", fn_root), "png", sep=".") )
        plt = carstm_map(  res=res, vn="predictions", tmatch=as.character(y), umatch=as.character(u),
          breaks=brks, 
          colors=rev(RColorBrewer::brewer.pal(5, "RdYlBu")),
          additional_features=additional_features,
          title=fn_root,  
          outfilename=outfilename
        )
      }
    }
      
 
```

# 1. Simplest model .. just CWD

```code
              mean    sd 0.025quant 0.5quant 0.975quant   mode kld
(Intercept) -7.629 0.014     -7.657   -7.629     -7.602 -7.629   0

Random effects:
  Name	  Model
    cwd RW2 model

Model hyperparameters:
                  mean    sd 0.025quant 0.5quant 0.975quant  mode
Precision for cwd 0.92 0.389      0.334    0.862       1.84 0.748

Deviance Information Criterion (DIC) ...............: 393406.04
Deviance Information Criterion (DIC, saturated) ....: 239126.71
Effective number of parameters .....................: 13.91

Watanabe-Akaike information criterion (WAIC) ...: 393403.56
Effective number of parameters .................: 11.43

Marginal log-Likelihood:  -196798.79 
 is computed 

$fixed_effects
                 mean        sd quant0.025  quant0.5 quant0.975          ID
(Intercept) 0.0004858 6.731e-06  0.0004724 0.0004858   0.000499 (Intercept)

$random_effects
        mean     sd quant0.025 quant0.5 quant0.975     ID
SD cwd 1.117 0.2528     0.7383    1.076      1.727 SD cwd
```


# 2.  CWD and time

```code
Fixed effects:
              mean   sd 0.025quant 0.5quant 0.975quant   mode   kld
(Intercept) -7.677 0.23     -8.245   -7.665     -7.192 -7.671 0.002

Random effects:
  Name	  Model
    time AR1 model
   cwd RW2 model

Model hyperparameters:
                     mean     sd 0.025quant 0.5quant 0.975quant  mode
Precision for time 16.931 12.229      2.195   13.847     47.433 6.655
Rho for time        0.876  0.093      0.645    0.900      0.987 0.961
Precision for cwd   0.935  0.401      0.360    0.866      1.909 0.737

Deviance Information Criterion (DIC) ...............: 390219.53
Deviance Information Criterion (DIC, saturated) ....: 235940.19
Effective number of parameters .....................: 35.17

Watanabe-Akaike information criterion (WAIC) ...: 390212.18
Effective number of parameters .................: 27.82

Marginal log-Likelihood:  -195241.48 
 is computed 

$fixed_effects
                 mean        sd quant0.025  quant0.5 quant0.975          ID
(Intercept) 0.0004753 0.0001134  0.0002629 0.0004686  0.0007496 (Intercept)

$random_effects
               mean      sd quant0.025 quant0.5 quant0.975           ID
SD time      0.3026 0.13807     0.1452   0.2655     0.6682      SD time
SD cwd       1.1047 0.23852     0.7255   1.0727     1.6579       SD cwd
Rho for time 0.8757 0.09279     0.6450   0.8974     0.9870 Rho for time

```

# CWD, time, space

```

Fixed effects:
              mean    sd 0.025quant 0.5quant 0.975quant   mode   kld
(Intercept) -8.006 0.268     -8.642   -7.993     -7.451 -7.999 0.004

Random effects:
  Name	  Model
    time AR1 model
   cwd RW2 model
   space BYM2 model

Model hyperparameters:
                      mean    sd 0.025quant 0.5quant 0.975quant  mode
Precision for time  14.287 8.864      2.488   12.389     35.552 7.468
Rho for time         0.876 0.080      0.684    0.894      0.982 0.943
Precision for cwd    0.929 0.391      0.350    0.867      1.859 0.747
Precision for space  2.047 0.207      1.661    2.040      2.475 2.032
Phi for space        0.962 0.031      0.880    0.971      0.996 0.987

Deviance Information Criterion (DIC) ...............: 377592.47
Deviance Information Criterion (DIC, saturated) ....: 223313.13
Effective number of parameters .....................: 411.58

Watanabe-Akaike information criterion (WAIC) ...: 377498.17
Effective number of parameters .................: 316.61

Marginal log-Likelihood:  -189196.38 
 is computed 

$fixed_effects
                 mean        sd quant0.025  quant0.5 quant0.975          ID
(Intercept) 0.0003455 9.994e-05  0.0001769 0.0003376   0.000579 (Intercept)

$random_effects
                mean      sd quant0.025 quant0.5 quant0.975            ID
SD time       0.3121 0.12086     0.1678   0.2811     0.6282       SD time
SD cwd        1.1080 0.24248     0.7353   1.0706     1.6812        SD cwd
SD space      0.7015 0.03539     0.6361   0.6999     0.7750      SD space
Rho for time  0.8759 0.08030     0.6845   0.8912     0.9818  Rho for time
Phi for space 0.9622 0.03103     0.8807   0.9705     0.9955 Phi for space

```

# CWD, time, space, t, z

```code

Fixed effects:
              mean    sd 0.025quant 0.5quant 0.975quant  mode   kld
(Intercept) -8.024 0.308     -8.734   -8.013     -7.382 -8.02 0.042

Random effects:
  Name	  Model
    time AR1 model
   inla.group(t, method = "quantile", n = 9) RW2 model
   inla.group(z, method = "quantile", n = 9) RW2 model
   cwd RW2 model
   space BYM2 model

Model hyperparameters:
                                                           mean      sd
Precision for time                                       12.417   9.608
Rho for time                                              0.893   0.085
Precision for inla.group(t, method = "quantile", n = 9) 201.720 316.983
Precision for inla.group(z, method = "quantile", n = 9)  13.276  17.753
Precision for cwd                                         0.938   0.396
Precision for space                                       2.215   0.233
Phi for space                                             0.963   0.035 

Deviance Information Criterion (DIC) ...............: 377264.52
Deviance Information Criterion (DIC, saturated) ....: 222985.18
Effective number of parameters .....................: 422.05

Watanabe-Akaike information criterion (WAIC) ...: 377167.21
Effective number of parameters .................: 324.06

Marginal log-Likelihood:  -189056.65 
 is computed 

$fixed_effects
                 mean        sd quant0.025  quant0.5 quant0.975          ID
(Intercept) 0.0003445 0.0001336  0.0001609 0.0003305  0.0006217 (Intercept)

$random_effects
                                               mean      sd quant0.025 quant0.5
SD time                                      0.3659 0.18391    0.16511  0.31451
SD inla.group(t, method = "quantile", n = 9) 0.1020 0.04642    0.03251  0.09558
SD inla.group(z, method = "quantile", n = 9) 0.3660 0.13999    0.13280  0.35545
SD cwd                                       1.1011 0.23576    0.72890  1.06839
SD space                                     0.6746 0.03573    0.61074  0.67212
Rho for time                                 0.8930 0.08510    0.67833  0.91436
Phi for space                                0.9630 0.03483    0.86922  0.97320 


```

```R

  M = size_distributions(p=p, toget="modal_groups_carstm_inputs", pg=sppoly, xrange=xrange, dx=dx )
  M = M[ -which(M$N == 0), ]
  M$sa[ !is.finite(M$sa) ] = 1
  M$logsa = log(M$sa)

  # model pa using all data
  res = NULL; gc()
  res = carstm_model( p=pN, 
    data=M, 
    sppoly=sppoly, 
    # theta = c( 0.926, 1.743, -0.401, 0.705, -2.574, 1.408, 2.390, 3.459, 3.321, -2.138, 3.083, -1.014, 3.558, 2.703 ),
    nposteriors=5000,
    posterior_simulations_to_retain=c( "summary", "random_spatial", "predictions"), 
    family = p$family, 
    verbose=TRUE,
    #redo_fit=FALSE, 
    # debug = "fit",
    # control.family=list(control.link=list(model="logit")),  # default for binomial .. no need to specify
    # control.inla = list( strategy="laplace", int.strategy="eb" ),
    num.threads="4:3"
  )


  ( res$summary) 
   
    vn = "time"
    ts =  res$random[[vn]]
    IDx = as.numeric(ts$ID)
    plot( mean ~ IDx, ts, type="b", lwd=1.5, xlab=vn)
    lines( quant0.025 ~ IDx, ts, col="gray", lty="dashed")
    lines( quant0.975 ~ IDx, ts, col="gray", lty="dashed")

 
    vn = "cwd"
    ts =  res$random[[vn]]
    IDx = as.numeric(ts$ID)
    plot( mean ~ IDx, ts, type="b", lwd=1.5, xlab=vn)
    lines( quant0.025 ~ IDx, ts, col="gray", lty="dashed")
    lines( quant0.975 ~ IDx, ts, col="gray", lty="dashed")

 
    vn = "inla.group(t, method = \"quantile\", n = 9)"
    ts =  res$random[[vn]]
    IDx = as.numeric(ts$ID)
    plot( mean ~ IDx, ts, type="b", lwd=1.5, xlab=vn)
    lines( quant0.025 ~ IDx, ts, col="gray", lty="dashed")
    lines( quant0.975 ~ IDx, ts, col="gray", lty="dashed")

 
    vn = "inla.group(z, method = \"quantile\", n = 9)"
    ts =  res$random[[vn]]
    IDx = as.numeric(ts$ID)
    plot( mean ~ IDx, ts, type="b", lwd=1.5, xlab=vn)
    lines( quant0.025 ~ IDx, ts, col="gray", lty="dashed")
    lines( quant0.975 ~ IDx, ts, col="gray", lty="dashed")

 
   H = carstm_model( p=pH, DS="carstm_modelled_summary", sppoly=sppoly  )
   N = carstm_model( p=pN, DS="carstm_modelled_summary", sppoly=sppoly  )
 

    vn = "cwd"
    ts =  H$random[[vn]] * N$random[[vn]]
    IDx = as.numeric(ts$ID)
    plot( mean ~ log(IDx), ts, type="b", lwd=1.5, xlab=vn)
    lines( quant0.025 ~ log(IDx), ts, col="gray", lty="dashed")
    lines( quant0.975 ~ log(IDx), ts, col="gray", lty="dashed")

```

# END
