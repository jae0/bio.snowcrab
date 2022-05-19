
# machine learning on snow crab

  
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
 
  areal_units_fn = attributes(sppoly)[["areal_units_fn"]]
  
  aufns = carstm_filenames( p=pN, returntype="carstm_modelled_fit", areal_units_fn=areal_units_fn )
  outfn = paste( gsub(".rdata", "", aufns), "aggregated_timeseries", "rdata", sep="." )
  load(outfn)
  ats = data.table(SM$RES)
  biom = melt( ats, id.vars="yrs", measure.vars=c( "cfanorth", "cfasouth", "cfa4x") )
  names(biom) = c("yr", "cfa", "biomass.kt")
  
  set = merge (set, biom, by=c("yr", "cfa") )
  set$biom_wt = set$biomass.kt / (set$totmass.male.com/ set$sa)   ## kt / (kg / km^2)
  set$biom_wt[ which(!is.finite(set$biom_wt) )] = 0
  set$region = as.numeric(as.factor(set$cfa))

  p$corners = data.frame(plon=c(220, 990), plat=c(4750, 5270) )
  
  mapyears = p$yrs
  outdir = file.path( work_root )
  
  interpolate.method='mba'
  theta=p$pres*25
  ptheta=theta/2.3
  idp=2
  log.variable=TRUE
  add.zeros=TRUE
  minN=10
  probs=c(0.025, 0.975)

    mapyears = sort( unique(set$yr) )

    nplon = length( seq(min(p$corners$plon), max(p$corners$plon), by = p$pres) )
    nplat = length( seq(min(p$corners$plat), max(p$corners$plat), by = p$pres) )

    predlocs = spatial_grid(p) 
    predlocs = planar2lonlat( predlocs,  proj.type=p$aegis_proj4string_planar_km   )
    predlocs$z = aegis_lookup( parameters="bathymetry", LOCS=predlocs[, c("lon", "lat")], project_class="core", output_format="points" , DS="aggregated_data", variable_name="z.mean", space_resolution=p$pres ) # core=="rawdata"
    aoi = geo_subset( spatial_domain=p$spatial_domain, Z=predlocs )
    # predlocs = predlocs[ aoi, ]
  
    v = "biom_wt"

      for ( y in mapyears ) {
        # y = 2021

          outfn = paste( v,y, sep=".")
          outloc = file.path( outdir,v)
          ref=y

          set_xyz = set[ which(set$yr==y), c("plon","plat",v) ]
          names( set_xyz) = c("plon", "plat", "z")
          set_xyz = na.omit(subset(set_xyz,!duplicated(paste(plon,plat))))
          if(nrow(set_xyz)<minN)next() #skip to next variable if not enough data

          offset = exp(-8)   # offset fot log transformation
          er = exp(c(-8, 3))  # range of all years
 
          er=c(0,1)
          withdata = which(is.finite( set_xyz$z ))
 
          ler = er
          
          if (length(withdata) < 3) print(paste("skipped",v, y, "<3 data points to create map", sep="."))
          if (length(withdata) < 3) next()
          S = set_xyz[ withdata, c("plon", "plat") ]


          distances =  rdist( predlocs[ aoi, c("plon", "plat")], S)
          distances[ which(distances < ptheta) ] = NA
          shortrange = which( !is.finite( rowSums(distances) ) )
          ips = aoi[ shortrange ]
          
          if(log.variable){
            set_xyz$z = log(set_xyz$z+offset)
            ler=log(er+offset)
            #if(offset<1)if(shift) xyz$z = xyz$z + abs(log(offset))
          }

          datarange = seq( ler[1], ler[2], length.out=50)
          xyzi = na.omit(set_xyz)

          if(nrow(xyzi)<minN||is.na(er[1]))next() #skip to next variable if not enough data

          #!# because 0 in log space is actually 1 in real space, the next line adds the log of a small number (offset)
          #!# surrounding the data to mimic the effect of 0 beyond the range of the data
          if(add.zeros)  xyzi =na.omit( zeroInflate(set_xyz,corners=p$corners,type=2,type.scaler=0.5,eff=log(offset),blank.dist=20) )

          if(interpolate.method=='mba'){
            u= MBA::mba.surf(x=xyzi[,c("plon","plat", "z")], nplon, nplat, sp=TRUE   )
            res = cbind( predlocs[ips, 1:2], u$xyz.est@data$z[ips] )
          }
        
          if(interpolate.method=='tps'){
            # broken ?
            u= fastTps(x=xyzi[,c("plon","plat")] , Y=xyzi[,'z'], theta=theta )
            res = cbind( predlocs[ips, 1:2], predict(u, xnew=predlocs[ips, 1:2]))
          }
          if(interpolate.method=='idw'){
            # broken?
            require(gstat)
            u = gstat(id = "z", formula = z ~ 1, locations = ~ plon + plat, data = xyzi, set = list(idp = idp))
            res = predict(u, predlocs[ips, 1:2])[,1:3]
          }
          #print(summary(set_xyz))
          #print(summary(res))

          xyz = res
          names( xyz) = c("plon", "plat", "z")
          #if(shift)xyz$z = xyz$z - abs(log(offset))

          cols = colorRampPalette(c("darkblue","cyan","green", "yellow", "orange","darkred", "black"), space = "Lab")

          xyz$z[xyz$z>ler[2]] = ler[2]
          if(ratio)xyz$z[xyz$z<ler[1]] = ler[1]

          ckey=NULL
          if(log.variable){
            # create labels for legend on the real scale
            labs=as.vector(c(1,2,5)%o%10^(-4:5))
            labs=labs[which(labs>er[1]&labs<er[2])]
            ckey=list(labels=list(at=log(labs+offset),labels=labs,cex=2))
          }

          dir.create (outloc, showWarnings=FALSE, recursive =TRUE)
          annot=ref
          filename=file.path(outloc, paste(outfn, "png", sep="."))
          print(filename)
          png( filename=filename, width=3072, height=2304, pointsize=40, res=300 )
          lp = aegis_map( xyz, xyz.coords="planar", depthcontours=TRUE, pts=set_xyz[,c("plon","plat")],
            annot=annot, annot.cex=4, at=datarange , col.regions=cols(length(datarange)+1),
            colpts=F, corners=p$corners, display=F, colorkey=ckey, plotlines="cfa.regions" )
          print(lp)
          dev.off()

      }
   
   



  sppoly = areal_units( p=pN )
  M = snowcrab.db( p=pN, DS="carstm_inputs", sppoly=sppoly )  # will redo if not found


  vn = "biom_wt"

  features = c( "z", "t", "pca1", "pca2", "no.female.mat", "no.female.imm", "no.male.mat", "no.male.imm", "ms.no.10", "ms.no.30", "sa", "region" )
  
  tokeep = c( "id", vn, features)

  MO = M[M$tag=="observations", which(names(M) %in% tokeep )]
  SO = set[ , which( names(set)  %in% tokeep ) ]

  MO = merge( MO, SO, by="id", suffixes=c("", ".set") ) 
  
  MO = MO[,c("biom_wt", features) ]
  MO = na.omit(MO)

  nMO = nrow(MO)
  ntrain = trunc(nMO * 0.75)
  ntest = nMO - ntrain

  jtrain = sample( nMO, trunc(nMO * 0.75) )
  jtest = setdiff(1:nMO, jtrain)

    
  require(tibble)

  y_train = MO[jtrain, vn]
  y_test = MO[jtest, vn]

  MO = as_tibble(  MO[, features] )

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


spec <- feature_spec( train_df, features ~ . ) %>% 
  step_numeric_column(all_numeric(), normalizer_fn = scaler_standard()) %>% 
  fit()
spec = dense_features(spec)

# feature_spec interface implemented in the tfdatasets package for normalization
# feature_columns also has options

layer <- layer_dense_features(
  feature_columns = dense_features(spec), 
  dtype = tf$float32
)
layer(train_df)

input <- layer_input_from_dataset( train_df  )
 
output =   
  layer_dense_features(input, spec) %>% 
  layer_dense(units = 128, activation = "relu") %>% 
  layer_dense(units = 64, activation = "relu") %>% 
  layer_dense(units = 1, activation = "relu")  
 
model = keras_model(input, output)
plot_model(model)

summary(model)
 
compile( model, 
  # optimizer = optimizer_rmsprop(),
  # optimizer = "adam",
    optimizer = 'rmsprop',
  #  optimizer = "adam",

  loss = "kl_divergence",
  # loss = "mse",
  # loss = "sparse_categorical_crossentropy",
  # loss = 'categorical_crossentropy',

  metrics = c('accuracy')
  # metrics = list("mean_absolute_error")
)


# train it;  
# Display training progress by printing a single dot for each completed epoch.

history <-  fit( model,
  x = train_df,
  y = y_train,
  epochs = 200,
  validation_split = 0.2,
  verbose = 0,
  callbacks = list( 
    print_dot_callback <- callback_lambda(
      on_epoch_end = function(epoch, logs) {
        if (epoch %% 50 == 0) cat("\n")
        cat(".")
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
plot( test_predictions[ , 1] ~ y_test)
cor( test_predictions[ , 1] , y_test, , use="pairwise.complete.obs" )   

save_model_tf(object = model, filepath = "model")

model <- load_model_tf("model")




   