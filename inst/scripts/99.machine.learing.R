
# machine learning on snow crab

  
Basic premise: 

- Identify data for inputs (features)
- Classify data by placement in healthy, poor, etc status

features to consider:

- number male > 95mm, 
- number male 60 - 95
- number female mature
- swept area
- SA of trawls (discretized)
- temperature
- time of year
- no. cod
- no. halibut
- no skate
- no shrimp
- fishery landings discretized area
- fishery effort discretized area
- fishery trap*days fishery
- SA of fishery (discretized)

ultimately estimate % of healthy in each zone with CI

  
  source( file.path( code_root, "bio_startup.R" )  )

  year.assessment = 2021

  p = bio.snowcrab::load.environment( year.assessment=year.assessment )
  set = snowcrab.db( DS="set.complete", p=p ) # note depth is log transformed here
  set$id = paste(set$trip, set$set, sep=".")

  # bring in params for carstm habitat analysis
  yrs = 1999:year.assessment
  spec_bio = bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 )

  snowcrab_filter_class = "fb"     # fishable biomass (including soft-shelled )
  runlabel = "1999_present_fb"

  # params of observation
  pH = snowcrab_parameters( 
    project_class="carstm", 
    yrs=yrs,  
    areal_units_type="tesselation", 
    family = "binomial",  # "binomial",  # "nbinomial", "betabinomial", "zeroinflatedbinomial0" , "zeroinflatednbinomial0"
    carstm_model_label= runlabel,  
    selection = list(
      type = "presence_absence",
      biologicals=list( spec_bio=spec_bio ),
      biologicals_using_snowcrab_filter_class=snowcrab_filter_class
    )
  )

  pN = snowcrab_parameters(
    project_class="carstm",
    yrs=yrs,   
    areal_units_type="tesselation",
    family="poisson",
    carstm_model_label= runlabel,  
    selection = list(
      type = "number",
      biologicals=list( spec_bio=spec_bio ),
      biologicals_using_snowcrab_filter_class=snowcrab_filter_class
    )
  )
 
  sppoly = areal_units( p=pN )

  areal_units_fn = attributes(sppoly)[["areal_units_fn"]]
  
  aufns = carstm_filenames( p=pN, returntype="carstm_modelled_fit", areal_units_fn=areal_units_fn )
  outfn = paste( gsub(".rdata", "", aufns), "aggregated_timeseries", "rdata", sep="." )
  load(outfn)
  
  ats = data.table(SM$RES)
  biom = melt( ats, id.vars="yrs", measure.vars=c( "cfanorth", "cfasouth", "cfa4x") )
  names(biom) = c("yr", "cfa", "biomass.kt")
  biom[ , hcr:=biomass.kt/max(biomass.kt, na.rm=TRUE) , by=cfa] # scale to max by region
  
  set = merge (set, biom, by=c("yr", "cfa") )
  set$biom_wt = set$biomass.kt / (set$totmass.male.com/ set$sa)   ## kt / (kg / km^2)
  set$biom_wt[ which(!is.finite(set$biom_wt) )] = 0
  set$region = as.numeric(as.factor(set$cfa))

  sset = data.table(set)
  sset = sset[ , mean(hcr, na.rm=TRUE), by=.(cfa, yr) ]
    
  dev.new(width=16, height=8, pointsize=20)
  ggplot( sset, aes(yr, V1, fill=cfa, colour=cfa) ) +
    geom_line() +
    labs(x="Year", y="HCR index", size = rel(1.5)) +
    # scale_y_continuous( limits=c(0, 300) )  
    theme_light( base_size = 22 ) 


  M = snowcrab.db( p=pN, DS="carstm_inputs", sppoly=sppoly )  # will redo if not found
  
  vn = "hcr"

  # cod, haddock, halibut, plaice, wolfish, thornyskate, smoothskate, winterskate, northernshrimp, jonahcrab, lessertoadcrab
  species = c(10, 11, 30, 40, 201, 50, 2521, 2511, 202, 204, 2211)
  sp_vns = paste( "ms.no.", species, sep="")

  features = c( "z", "t", "pca1", "pca2", "no.female.mat", "no.female.imm", "no.male.mat", "no.male.imm", 
    sp_vns,
    "landings", "effort", "cpue", "plon", "plat",
    "sa", "region" )
  
  tokeep = c( "id", vn, features)

  toscale = setdiff( features, c("id", "region", vn ) ) 


  dta = M[ which(M$tag=="observations"),  ]
  dta = merge( dta, set[which(!duplicated(set$id)) , ], by="id", suffixes=c("", ".set"), all.x=TRUE, all.y=FALSE ) 

  MO = dta[ , tokeep ]
  MO = na.omit(MO)
  dta = dta[ which(dta$id %in% MO$id), ]

  MO = MO[ , c(vn, features) ]
  nMO = nrow(MO)
  fraction_train = 0.6

  ntrain = trunc(nMO * fraction_train )
  ntest = nMO - ntrain

  jtrain = sample( nMO, ntrain )
  jtest = setdiff(1:nMO, jtrain)

  y_train = MO[jtrain, vn]
  y_test = MO[jtest, vn]
    

  MO = MO[, features] 
  # rr = to_categorical( MO$region )
  # rr = rr[, which(colSums(rr) > 0 )]
  # colnames(rr) = paste("region", c(1:3), sep="_" )
  # MO = cbind(MO, rr)
  # MO$region = NULL

  # MO[,toscale] = scale(MO[,toscale])
  
  MO$region = factor( MO$region, levels=c(1,2,3), labels=c("cfa4x", "cfanorth", "cfasouth") ) 
  MO$region = as.character( MO$region )

  require(tibble)
  MO = as_tibble(  MO  )

  train_df = MO[jtrain, ]
  test_df  = MO[jtest,  ]
 

  # merge other info to M


  # look up habitat score and quantize it: 0-5 or 0-9?

  # ML + predictions  
  
# input layer
if (0) {
  install.packages( "reticulate", deps=TRUE )
  install.packages( "keras", deps=TRUE )
  install.packages( "progress", deps=TRUE  ) 

}


library(ggplot2)

library(reticulate)
use_python("/usr/bin/python")

library(keras)
library(tfdatasets)


spec <- feature_spec( train_df, hcr ~ . )  
spec <- spec %>% 
  step_numeric_column( all_numeric(), normalizer_fn = scaler_standard() ) %>% 
  step_categorical_column_with_vocabulary_list(region, vocabulary_list=c("cfanorth", "cfasouth", "cfa4x") )  # (0..num_buckets)

spec = fit(spec)
# str( spec$dense_features() )
  
nc = ncol(train_df)

input <- layer_input_from_dataset( train_df  )
 
output =   
  layer_dense_features( input, feature_columns = spec$dense_features() ) %>% 
  layer_dense(units = nc*3, activation = "LeakyReLU") %>% 
  layer_dense(units = nc, activation = "LeakyReLU") %>%  
  layer_dense(units = trunc(nc/3), activation = "LeakyReLU") %>%  
  layer_dense(units = 3, activation = "LeakyReLU") %>%  
  layer_dense(units = 1, activation = "relu")  
 
model = keras_model(input, output)

# summary(model)
 
compile( model, 
  # optimizer = "adam",
  # optimizer = "nadam",
  optimizer = 'rmsprop',
 
  # loss = "kl_divergence",
  loss = "mse",
  # loss = "sparse_categorical_crossentropy",
  # loss = 'categorical_crossentropy',

  # metrics = c('accuracy')
  metrics = list("mean_absolute_error")
)

# train it;  
# Display training progress by printing a single dot for each completed epoch.
nepochs = 100

history <-  fit( model,
  x = train_df,
  y = y_train,
  epochs = nepochs,
  validation_data = list(x_val=test_df, yval=y_test),
  verbose = 0,
  callbacks = list( 
    print_dot_callback <- callback_lambda(
      on_epoch_end = function(epoch, logs) {
        nn = trunc(nepochs/10)
        if (epoch %% nn  == 0) cat(trunc(epoch/nn))
        # if (epoch %% 100 == 0) cat("\n")
      }
    )
  )
) 

plot(history)

# model performance on test data
out = model %>% evaluate(test_df, y_test, verbose = 0) 
out


# By default predict will return the output of the last Keras layer (probability for each class) 
pp = as.vector( model %>% predict(train_df ) )
plot( pp ~ y_train )
cor( pp, y_train, use="pairwise.complete.obs" )

test_predictions <- as.vector( model %>% predict(test_df  ) )
plot( test_predictions ~ y_test)
cor( test_predictions  , y_test, , use="pairwise.complete.obs" )   

MO$preds = NA
MO$preds[jtrain] = pp
MO$preds[jtest]  = test_predictions

MO$hcr = NA
MO$hcr[jtrain] = y_train
MO$hcr[jtest]  = y_test


Metrics::auc(MO$hcr, MO$preds)

out = as.data.table(MO)
out$year = dta$year

sout = out[ , mean(preds, na.rm=TRUE), by=.(region, year) ]
 
dev.new(width=16, height=8, pointsize=20)
ggplot( sout, aes(year, V1, fill=region, colour=region) ) +
  geom_line() +
  labs(x="Year", y="HCR index", size = rel(1.5)) +
  # scale_y_continuous( limits=c(0, 300) )  
  theme_light( base_size = 22 ) 



save_loc = file.path(work_root, "model")
save_model_tf(object = model, filepath = save_loc )

model <- load_model_tf(save_loc) 
