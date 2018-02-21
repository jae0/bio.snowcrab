
  if (!exists("year.assessment")) {
    year.assessment=lubridate::year(Sys.Date()) -1
    year.assessment=lubridate::year(Sys.Date())
  }

# --------------------------------------------------------------
#  Ensure the following scripts complete without error:
#  these external dependencies permit lookup of data, for this script
#Require the following:
#snowcrab.db("complete.redo")

#Run aegis::(inst/scripts/11.surveys.r), need to update year within or run next line
  #system.file(package="aegis", "scripts", "11.surveys.r")

#Run aegis::(inst/scripts/temperature.r), need to update year within
  #system.file(package="aegis", "scripts", "07.temperature.R") or run next line

#Substrate and bathymetry can be run (as below) if suspect significant changes in one or both of these datasets

#BZ 2017 these lines below can directly run the indicators without goint to run individual scripts
#  system.file(package="aegis", "scripts", "07.temperature.R")
#  system.file(package="aegis", "scripts", "10.lookuptables.r")
#  system.file(package="aegis", "scripts", "11.surveys.r")
#  system.file(package="aegis", "scripts", "16.speciescomposition.R")

# 1. Define some additional starting parameters for debugging
#    choose various over-rides: these are initially defined in parameters.r


# --------------------------------------------------------------
# using environmental data ... estimate/lookup missing environmental data .. (t,z)
#BZ 2017 below lines shouldn't be required. datasets created in 01.snowcrab
#logbook.db( DS ="fisheries.complete.redo", p=p )
#snowcrab.db( DS ="set.complete.redo", p=p )

# -------------------------------------------------------------------------------------
# abundance .. positive valued data .. 
# takes about 5 hrs .. ~1 GB / process
# vn = "snowcrab.large.males_abundance"
# year.assessment = 2017
p = bio.snowcrab::load.environment( year.assessment=year.assessment ) 

# 11 hrs with these settings,
p = snowcrab_stmv( p=p, DS="parameters",
  selection=list(
    name = "snowcrab.large.males_abundance",
    type = "abundance",
    sex=0, # male
    mat=1, # do not use maturity status in groundfish data as it is suspect ..
    spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 ),
    len= c( 95, 200 )/10, #  mm -> cm ; aegis_db in cm
    drop.groundfish.data=TRUE # esp from 1970 to 1999 measurement of invertebrates was sporatic .. zero-values are dropped as they are unreliable
  ),
  DATA = 'snowcrab_stmv( p=p, DS="stmv_inputs" )',
  stmv_global_family = gaussian(link=log),
  stmv_local_modelengine = "twostep",
  stmv_twostep_space = "krige",
  stmv_gam_optimizer=c("outer", "bfgs") ,
  stmv_distance_statsgrid = 2, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
  stmv_distance_scale = 50,
  stmv_global_family = gaussian(link="log"),
  stmv_Y_transform =list(
        transf = function(x) {x/6675} ,
        invers = function(x) {x*6675}
  ) # transform data to unit interval to stabilize variance and speed up convergence
  # stmv_global_family = gaussian(link="identity"),
  # stmv_Y_transform =list(
  #       transf = function(x) {log(x/6675)} ,
  #       invers = function(x) {exp(x)*6675}
  # ) # transform data to unit interval to stabilize variance and speed up convergence
)

# range( INP$snowcrab.large.males_abundance )
# [1]   14.3 6675.0

# o = snowcrab_stmv(p=p, DS="stmv_inputs" )  # create fields for  

stmv( p=p, runmode=c("globalmodel", "interpolate" )  ) #  for a clean start
# stmv( p=p, runmode=c("globalmodel", "interpolate" ), use_saved_state=TRUE ) # for a restart  .. not working?

# if (really.finished) stmv_db( p=p, DS="cleanup.all" )

snowcrab_stmv( p=p, DS="predictions.redo" ) # warp predictions to other grids (if any)
snowcrab_stmv( p=p, DS="stmv.stats.redo" ) # warp stats to other grids (if any)
snowcrab_stmv( p=p, DS="complete.redo" )
snowcrab_stmv( p=p, DS="baseline.redo" )
snowcrab_stmv( p=p, DS="map.all" )

global_model = stmv_db( p=p, DS="global_model")
summary( global_model )
plot(global_model)

#Below lines are just model outputs for comparison sake
Family: gaussian 
Link function: identity 

Formula:
snowcrab.large.males_abundance ~ s(t, k = 3, bs = "ts") + s(tmean.climatology, 
    k = 3, bs = "ts") + s(tsd.climatology, k = 3, bs = "ts") + 
    s(log(z), k = 3, bs = "ts") + s(log(dZ), k = 3, bs = "ts") + 
    s(log(ddZ), k = 3, bs = "ts") + s(log.substrate.grainsize, 
    k = 3, bs = "ts") + s(pca1, k = 3, bs = "ts") + s(pca2, k = 3, 
    bs = "ts")

Parametric coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept)  -3.4506     0.0197    -175   <2e-16

Approximate significance of smooth terms:
                             edf Ref.df      F p-value
s(t)                       1.274      2  38.37 < 2e-16
s(tmean.climatology)       0.829      2   2.94  0.0044
s(tsd.climatology)         1.981      2 231.52 < 2e-16
s(log(z))                  1.634      2 141.57 < 2e-16
s(log(dZ))                 1.941      2  10.24 2.4e-05
s(log(ddZ))                1.696      2  17.52 2.3e-09
s(log.substrate.grainsize) 1.915      2  58.44 < 2e-16
s(pca1)                    1.990      2 220.02 < 2e-16
s(pca2)                    1.917      2 137.45 < 2e-16

R-sq.(adj) =  0.365   Deviance explained = 36.7%
GCV = 0.01104  Scale est. = 0.011016  n = 7255
 
# variation 2:
# Family: gaussian 
# Link function: log 
# 
# Formula:
# snowcrab.large.males_abundance ~ s(t, k = 3, bs = "ts") + s(tmean.climatology, 
#     k = 3, bs = "ts") + s(tsd.climatology, k = 3, bs = "ts") + 
#     s(log(z), k = 3, bs = "ts") + s(log(dZ), k = 3, bs = "ts") + 
#     s(log(ddZ), k = 3, bs = "ts") + s(log.substrate.grainsize, 
#     k = 3, bs = "ts") + s(pca1, k = 3, bs = "ts") + s(pca2, k = 3, 
#     bs = "ts")
# 
# Parametric coefficients:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)  -2.1317     0.0238   -89.4   <2e-16
# 
# Approximate significance of smooth terms:
#                                edf Ref.df      F p-value
# s(t)                       0.00044      2   0.00   0.196
# s(tmean.climatology)       1.99802      2  27.33 8.8e-13
# s(tsd.climatology)         1.73172      2  92.23 < 2e-16
# s(log(z))                  1.49345      2 162.71 < 2e-16
# s(log(dZ))                 1.91637      2  53.49 < 2e-16
# s(log(ddZ))                0.72659      2   1.67   0.028
# s(log.substrate.grainsize) 1.95224      2  44.19 < 2e-16
# s(pca1)                    1.99998      2 102.12 < 2e-16
# s(pca2)                    2.00000      2 125.48 < 2e-16
# 
# R-sq.(adj) =  0.231   Deviance explained = 23.3%
# GCV = 0.00014973  Scale est. = 0.00014942  n = 7255
# 


# -------------------------------------------------
# presence-absence
# this takes about 40 hrs ... and 5-6 GB /process 
# year.assessment = 2017
p = bio.snowcrab::load.environment( year.assessment=year.assessment )
p = snowcrab_stmv( p=p, DS="parameters",
  selection=list(
    name = "snowcrab.large.males_presence_absence",
    type = "presence_absence",
    sex=0, # male
    mat=1, # do not use maturity status in groundfish data as it is suspect ..
    spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 ),
    len= c( 95, 200 )/10, #  mm -> cm ; aegis_db in cm
    drop.groundfish.data=TRUE # esp from 1970 to 1999 measurement of invertebrates was sporatic .. zero-values are dropped as they are unreliable
  ),
  DATA = 'snowcrab_stmv( p=p, DS="stmv_inputs" )',
  stmv_global_family = binomial(),
  stmv_local_modelengine = "twostep",
  stmv_twostep_space = "krige",
  stmv_gam_optimizer=c("outer", "bfgs"),
  stmv_distance_statsgrid = 2, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** ),
  stmv_distance_scale = 50
)
# o = snowcrab_stmv(p=p, DS="stmv_inputs" )  # create fields for
stmv( p=p, runmode=c("globalmodel", "interpolate" ) ) # no global_model and force a clean restart

# stmv_db( p=p, DS="stmv.results" ) # save to disk for use outside stmv*, returning to user scale
# if (really.finished) stmv_db( p=p, DS="cleanup.all" )

snowcrab_stmv( p=p, DS="predictions.redo" ) # warp predictions to other grids
snowcrab_stmv( p=p, DS="stmv.stats.redo" ) # warp stats to other grids
snowcrab_stmv( p=p, DS="complete.redo" )
snowcrab_stmv( p=p, DS="baseline.redo" )
snowcrab_stmv( p=p, DS="map.all" )

global_model = stmv_db( p=p, DS="global_model")
summary( global_model )
plot(global_model)
Family: binomial 
Link function: logit 
Formula:
snowcrab.large.males_presence_absence ~ s(t, k = 3, bs = "ts") + 
    s(tmean.climatology, k = 3, bs = "ts") + s(tsd.climatology, 
    k = 3, bs = "ts") + s(log(z), k = 3, bs = "ts") + s(log(dZ), 
    k = 3, bs = "ts") + s(log(ddZ), k = 3, bs = "ts") + s(log.substrate.grainsize, 
    k = 3, bs = "ts") + s(pca1, k = 3, bs = "ts") + s(pca2, k = 3, 
    bs = "ts")

Parametric coefficients:
            Estimate Std. Error z value Pr(>|z|)
(Intercept)   2.2456     0.0433    51.9   <2e-16

Approximate significance of smooth terms:
                             edf Ref.df  Chi.sq p-value
s(t)                       1.993      2  359.98 < 2e-16
s(tmean.climatology)       1.999      2   22.05 1.5e-05
s(tsd.climatology)         1.999      2  537.87 < 2e-16
s(log(z))                  1.931      2 1560.45 < 2e-16
s(log(dZ))                 0.498      2    0.96    0.16
s(log(ddZ))                1.955      2   44.55 6.3e-11
s(log.substrate.grainsize) 1.831      2   80.98 < 2e-16
s(pca1)                    1.997      2  205.31 < 2e-16
s(pca2)                    1.997      2  574.22 < 2e-16

R-sq.(adj) =  0.598   Deviance explained = 55.4%
UBRE = -0.71562  Scale est. = 1         n = 25468
 


# collect all predictions into a single file and return:
p = bio.snowcrab::load.environment( year.assessment=year.assessment )

p$selection=list(
  name = "snowcrab.large.males_abundance",
  type = "abundance",
  sex=0, # male
  mat=1, # do not use maturity status in groundfish data as it is suspect ..
  spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 ),
  len= c( 95, 200 )/10, #  mm -> cm ; aegis_db in cm
  drop.groundfish.data=TRUE # esp from 1970 to 1999 measurement of invertebrates was sporatic .. zero-values are dropped as they are unreliable
)

p = snowcrab_stmv( p=p, DS="parameters" )

interpolation.db( DS="biomass.redo", p=p  )
interpolation.db( DS="biomass.map", p=p  )

K = interpolation.db( DS="timeseries", p=p  )
str(K)

# table.view( K )

figure.timeseries.snowcrab.habitat(p=p) # /bio.data/bio.snowcrab/assessments/2016/timeseries/interpolated/snowcrab.habitat.sa.png

figure.timeseries.snowcrab.habitat.temperatures(p=p) # /bio.data/bio.snowcrab/assessments/2016/timeseries/interpolated/mean.bottom.temp.snowcrab.habitat.png



# update data summaries of the above results
p$vars.tomodel="R0.mass"
biomass.summary.db("complete.redo", p=p) #Uses the model results to create a habitat area expanded survey index



### --------- prediction success:
set = snowcrab_stmv(p=p, DS="input_data", voi=p$selection$name )

S = set[ , c("plon", "plat") ]

ii = array_map( "xy->1", S, gridparams=p$gridparams )
bs = bathymetry.db(p=p, DS="baseline")
bb = array_map( "xy->1", bs, gridparams=p$gridparams )
im = match(  ii, bb )
it = match( set$yr, p$yrs )

bm = interpolation.db( DS="biomass", p=p  )
spred = bm$m[cbind(im, it)]  # approximate match (ignoring seasonality)

summary ( lm(spred~snowcrab.large.males_abundance, data=set, na.actio="na.omit" ) )
plot(spred~snowcrab.large.males_abundance, data=set )
cor(spred,set$snowcrab.large.males_abundance, use="complete.obs")


# determine presence absence(Y) and weighting(wt)
#      set$weekno = floor(set$julian / 365 * 52) + 1
#      set$dyear = floor(set$julian / 365 ) + 1
