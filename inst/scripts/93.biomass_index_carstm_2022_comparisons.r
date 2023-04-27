
  source( file.path( code_root, "bio_startup.R" )  )
   
  require(bio.snowcrab)   # loadfunctions("bio.snowcrab") 

  year.assessment = 2022
  
  yrs = 1999:year.assessment
  spec_bio = bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 )
  snowcrab_filter_class = "fb"     # fishable biomass (including soft-shelled )  "m.mat" "f.mat" "imm"

  loadfunctions("bio.snowcrab")
  
  # key name

  runlabel= paste( "1999_present", snowcrab_filter_class, sep="_" )
    carstm_modelengine = "inla" 

  # variation 1 is about removing spe comp and temp
  runlabel= paste( "1999_2022_variation1", snowcrab_filter_class, sep="_" )
    carstm_modelengine = "inla.reduced1"
  
  # removing stations in 2015
  runlabel= paste( "1999_2022_variation2", snowcrab_filter_class, sep="_" )
    carstm_modelengine = "inla"
  
  # variation 3 is about removing spatial and spatiotemporal effects
  runlabel= paste( "1999_2022_variation3", snowcrab_filter_class, sep="_" )
    carstm_modelengine = "inla.reduced3"


  # params for number
  pN = snowcrab_parameters(
    project_class="carstm",
    yrs=yrs,   
    areal_units_type="tesselation",
    family = "nbinomial",
    carstm_model_label= runlabel,  
    carstm_modelengine = carstm_modelengine,
    selection = list(
      type = "number",
      biologicals=list( spec_bio=spec_bio ),
      biologicals_using_snowcrab_filter_class=snowcrab_filter_class
    )
  )

  # params for mean size .. mostly the same as pN
  pW = snowcrab_parameters(
    project_class="carstm",
    yrs=yrs,   
    areal_units_type="tesselation",
    family =  "gaussian",
    carstm_model_label= runlabel,  
    carstm_modelengine = carstm_modelengine,
    selection = list(
      type = "meansize",
      biologicals=list( spec_bio=spec_bio ),
      biologicals_using_snowcrab_filter_class=snowcrab_filter_class
    )
  )

  # params for probability of observation
  pH = snowcrab_parameters( 
    project_class="carstm", 
    yrs=yrs,  
    areal_units_type="tesselation", 
    family = "binomial",  # "binomial",  # "nbinomial", "betabinomial", "zeroinflatedbinomial0" , "zeroinflatednbinomial0"
    carstm_model_label= runlabel,  
    carstm_modelengine = carstm_modelengine,
    selection = list(
      type = "presence_absence",
      biologicals=list( spec_bio=spec_bio ),
      biologicals_using_snowcrab_filter_class=snowcrab_filter_class
    )
  )
  
  # use what was defined in the main script
  sppoly=areal_units( p=pN )

  additional_features = snowcrab_features_tmap(pN)  # for mapping below
 
  tmap_mode("plot")
   


# -------------------------------------------------
# Part 1 -- construct basic parameter list defining the main characteristics of the study
# require(aegis)


  M = snowcrab.db( p=pN, DS="carstm_inputs", sppoly=sppoly  )  # will redo if not found
  # can just copy datafile: carstm_inputs_snowcrab~tesselation~1~snowcrab~24~1~none~snowcrab_managementareas.rdata
  # into working directory to speed things up

# extract results
  p = pN  #  bio.snowcrab::snowcrab_parameters( project_class="carstm", yrs=1999:year.assessment )
  p = pW 
  p = pH 
 
  io = which(M$tag=="observations")
  ip = which(M$tag=="predictions")
  iq = unique( c( which( M$totno > 0), ip ) )  # for totno
  iw = unique( c( which( M$totno > 5), ip ) )  # need a good sample to estimate mean size
  iqo = intersect(iq, io)
  

  M$density = M$totno / M$data_offset 
  
  if (p==pN) Mfit = M[iq, ]
  if (p==pW) Mfit = M[iw, ]
  if (p==pH) Mfit = M

  iobs = which(Mfit$tag=="observations")

  
  fit = carstm_model( p=p, DS="carstm_modelled_fit",  sppoly = sppoly )  # extract currently saved model fit

  obs = Mfit$density[iobs]
  obs = Mfit$meansize[iobs]
  obs = Mfit$pa[iobs]

  preds = fit$summary.fitted.values[["mean"]][iobs]
  
  plot( log(preds) ~ log(obs) )
  # plot( (preds) ~ jitter(obs) )


  outp =c( 
    pearson=cor( preds,  obs ), 
    spearman=cor(preds, obs, method = "spearman" ),
    dic=fit$dic$dic, 
    p.eff=fit$dic$p.eff,
    waic=fit$waic$waic,
    waic.p.eff=fit$waic$p.eff,
    mlk=fit$mlik[["log marginal-likelihood (Gaussian)",1]],
    n = round(length(iobs),0)
  )
 

  # # plot(fit)
  # plot( fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
  # names(fit$marginals.hyperpar)
  # plot( fit$marginals.hyperpar$"Phi for space_time", type="l")  # posterior distribution of phi nonspatial dominates
  # plot( fit$marginals.hyperpar$"Precision for space_time", type="l")
  # plot( fit$marginals.hyperpar$"size for the nbinomial observations (1/overdispersion)", type="l")
  # # fit = NULL

 
  res = carstm_model( p=p, DS="carstm_modelled_summary",  sppoly = sppoly ) # to load currently saved results

  options(width =140)
  options(scipen =4, digits=4)

  res$summary$fixed_effects$cv = round( res$summary$fixed_effects$sd / res$summary$fixed_effects$mean * 100, 2 )
  
  res$summary$random_effects$cv = round( res$summary$random_effects$sd / res$summary$random_effects$mean * 100, 2 )
  outeff = rbind( 
    round( res$summary$fixed_effects[, c("mean", "sd", "cv") ], 5),
    round( res$summary$random_effects[, c("mean", "sd", "cv")] , 5 )
  )
  rownames(outeff) = gsub( ', method = \"quantile\", n = 11)', '', rownames(outeff), fixed=TRUE)
  rownames(outeff) = gsub( 'inla.group(', '', rownames(outeff), fixed=TRUE)
  rownames(outeff) = gsub( 'size for the nbinomial observations (1/overdispersion)', 'nbin 1/overdispersion', rownames(outeff), fixed=TRUE)

  out = list( params=outeff, gof=round(outp,3)  )

  out


 
  

  runlabel= paste( "1999_present", snowcrab_filter_class, sep="_" )
    carstm_modelengine = "inla" 
 
  # params for number
  pN0 = snowcrab_parameters(
    project_class="carstm",
    yrs=yrs,   
    areal_units_type="tesselation",
    family = "nbinomial",
    carstm_model_label= runlabel,  
    carstm_modelengine = carstm_modelengine,
    selection = list(
      type = "number",
      biologicals=list( spec_bio=spec_bio ),
      biologicals_using_snowcrab_filter_class=snowcrab_filter_class
    )
  )

  # params for mean size .. mostly the same as pN
  pW0 = snowcrab_parameters(
    project_class="carstm",
    yrs=yrs,   
    areal_units_type="tesselation",
    family =  "gaussian",
    carstm_model_label= runlabel,  
    carstm_modelengine = carstm_modelengine,
    selection = list(
      type = "meansize",
      biologicals=list( spec_bio=spec_bio ),
      biologicals_using_snowcrab_filter_class=snowcrab_filter_class
    )
  )

  # params for probability of observation
  pH0 = snowcrab_parameters( 
    project_class="carstm", 
    yrs=yrs,  
    areal_units_type="tesselation", 
    family = "binomial",  # "binomial",  # "nbinomial", "betabinomial", "zeroinflatedbinomial0" , "zeroinflatednbinomial0"
    carstm_model_label= runlabel,  
    carstm_modelengine = carstm_modelengine,
    selection = list(
      type = "presence_absence",
      biologicals=list( spec_bio=spec_bio ),
      biologicals_using_snowcrab_filter_class=snowcrab_filter_class
    )
  )
  


  # removing stations in 2015
  runlabel= paste( "1999_2022_variation2", snowcrab_filter_class, sep="_" )
    carstm_modelengine = "inla"
  
  # params for number
  pN2 = snowcrab_parameters(
    project_class="carstm",
    yrs=yrs,   
    areal_units_type="tesselation",
    family = "nbinomial",
    carstm_model_label= runlabel,  
    carstm_modelengine = carstm_modelengine,
    selection = list(
      type = "number",
      biologicals=list( spec_bio=spec_bio ),
      biologicals_using_snowcrab_filter_class=snowcrab_filter_class
    )
  )

  # params for mean size .. mostly the same as pN
  pW2 = snowcrab_parameters(
    project_class="carstm",
    yrs=yrs,   
    areal_units_type="tesselation",
    family =  "gaussian",
    carstm_model_label= runlabel,  
    carstm_modelengine = carstm_modelengine,
    selection = list(
      type = "meansize",
      biologicals=list( spec_bio=spec_bio ),
      biologicals_using_snowcrab_filter_class=snowcrab_filter_class
    )
  )

  # params for probability of observation
  pH2 = snowcrab_parameters( 
    project_class="carstm", 
    yrs=yrs,  
    areal_units_type="tesselation", 
    family = "binomial",  # "binomial",  # "nbinomial", "betabinomial", "zeroinflatedbinomial0" , "zeroinflatednbinomial0"
    carstm_model_label= runlabel,  
    carstm_modelengine = carstm_modelengine,
    selection = list(
      type = "presence_absence",
      biologicals=list( spec_bio=spec_bio ),
      biologicals_using_snowcrab_filter_class=snowcrab_filter_class
    )
  )
  


  SM0 = aggregate_simulations( 
    fn=carstm_filenames( pN0, returnvalue="filename", fn="aggregated_timeseries" ), 
    redo=FALSE
  )$cfasouth
  SM0 = SM0[ which(rownames(SM0)=="2015"),]
  SM0 = data.frame( B=SM0, Scenario="0")

  SM2 = aggregate_simulations( 
    fn=carstm_filenames( pN2, returnvalue="filename", fn="aggregated_timeseries" ), 
    redo=FALSE
  )$cfasouth
  SM2 = SM2[ which(rownames(SM2)=="2015"),]
  SM2= data.frame( B=SM2, Scenario="1")   # called scenario 1 in doc

  out =  rbind(SM0, SM2)
  
  ggplot(out, aes(x = B, fill = Scenario)) +                       # Draw overlaying histogram
  geom_histogram(position = "identity", alpha = 0.52, bins = 50)

  
  

## specific to are in question:

  fnout = file.path( pN2$modeldir, pN2$carstm_model_label, "unsampled_polygon.dat" )
  X = read.csv  (fnout) # outline of unsampled region
  X$X= NULL
  X = as.matrix( rbind(X, X[1,]))
  
  X = st_sfc( st_multipoint( X ) )
  X = st_cast(X, "POLYGON" )
  X = st_simplify(X)
  X = st_make_valid(X)
  st_crs(X) = crs_lonlat
  X = st_transform( X, st_crs( sppoly ) )
  X = st_cast( st_as_sf(X), "MULTIPOLYGON" )

  ooo = st_intersection( X, sppoly )
  ooo$surfacearea = st_area( ooo )
  sppoly[[ "aoi_sa_km2" ]] = 0
  j = match( ooo$AUID, sppoly$AUID )      
  if (length(j) > 0)  sppoly[[ "aoi_sa_km2" ]][j] = ooo$surfacearea
   
  sims0 = carstm_posterior_simulations( pN=pN0, pW=pW0, pH=pH0, sppoly=sppoly, pa_threshold=0.05, qmax=0.99 )
  sims0 = sims0  / 10^6 # 10^6 kg -> kt;; kt/km^2
  
  SM0_aoi = colSums( sims0 * sppoly$aoi_sa_km2, na.rm=TRUE )
  SM0_aoi = SM0_aoi[ which(rownames(SM0_aoi)=="2015"),]

 
  sims2 = carstm_posterior_simulations( pN=pN2, pW=pW2, pH=pH2, sppoly=sppoly, pa_threshold=0.05, qmax=0.99 )
  sims2 = sims2  / 10^6 # 10^6 kg -> kt;; kt/km^2
  
  SM2_aoi = colSums( sims2 * sppoly$aoi_sa_km2, na.rm=TRUE )
  SM2_aoi = SM2_aoi[ which(rownames(SM2_aoi)=="2015"),]


  SM0= data.frame( B=SM0_aoi, Scenario="0")   # called scenario 1 in doc
  SM2= data.frame( B=SM2_aoi, Scenario="1")   # called scenario 1 in doc

  out =  rbind(SM0, SM2)
  
  ggplot(out, aes(x = B, fill = Scenario)) +                       # Draw overlaying histogram
  geom_histogram(position = "identity", alpha = 0.52, bins = 50)






  outputdir = file.path( p$modeldir, p$carstm_model_label, "residuals" )

  # extract results and examine
  # prediction surface
  crs_lonlat = st_crs(projection_proj4string("lonlat_wgs84"))

  sppoly = areal_units( p=p )  # will redo if not found
  sppoly = st_transform(sppoly, crs=crs_lonlat )

  # do this immediately to reduce storage for sppoly (before adding other variables)
  S = snowcrab.db( p=p, DS="biological_data" )  # will redo if not found .. not used here but used for data matching/lookup in other aegis projects that use bathymetry
  S$tiyr=lubridate::decimal_date(S$timestamp)

  # S$totno = S$totno_adjusted / S$cf_set_no   # convert density to counts
  # S$totwgt = S$totwgt_adjusted / S$cf_set_mass # convert density to total wgt
  # S$data_offset = 1 / S$cf_set_no  ## offset only used in poisson model

  # reduce size
  S = S[ which( S$lon > p$corners$lon[1] & S$lon < p$corners$lon[2]  & S$lat > p$corners$lat[1] & S$lat < p$corners$lat[2] ), ]
  # levelplot(z.mean~plon+plat, data=S, aspect="iso")

  S$AUID = st_points_in_polygons(
    pts = st_as_sf( S, coords=c("lon","lat"), crs=crs_lonlat ),
    polys = sppoly[, "AUID"],
    varname="AUID"
  )

  S = S[!is.na(S$AUID),]

  names(S)[which(names(S)=="yr") ] = "year"
  # S = S[ which(S$year %in% p$yrs), ]
  # S$tiyr = lubridate::decimal_date ( S$timestamp )
  # S$dyear = S$tiyr - S$year

  obsMM = M[M$tag=="observations",]
  plocs = M[M$tag=="predictions",]

  obs = S
  obs$density = obs$totwgt / obs$data_offset
 
  rr = as.data.table(sims)
 
  rr$AUID = as.character( rr$space)
  # rr$AUID = match( rr$AUID, sppoly$AUID)

  rr$year = as.numeric( as.character( rr$time) )
  st = rr[,.(mean=mean(value)), by=.(AUID, year) ]

  obs = merge( obs, st, by=c("year", "AUID"), all.x=TRUE, all.y=FALSE )
  obs$mean[ !is.finite(obs$mean) ] = 0
  obs$resid =  obs$mean - obs$density
  obs$resid_per_set = obs$resid * obs$data_offset
  obs$yr = obs$year

  plot(mean~density, obs)


  vn = "resid"
  vn = "resid_per_set"
  #er = range( obs[,vn], na.rm=T) * c(0.95, 1.05)
  er = c(-100, 100)

  resol = p$pres

  B = bathymetry_db(p=p, DS="baseline")  # 1 km (p$pres )

  for ( y in  2000:2018 ) {
      ii = which( obs$yr==y & is.finite(obs[,vn] ))
      if ( length(ii) > 3 ) {
      dir.create( file.path( outputdir, "residuals", p$carstm_model_label), recursive=TRUE, showWarnings =FALSE)
      fn = file.path( outputdir, "residuals", p$carstm_model_label, paste( "residuals", y, "png", sep=".") )
      png( filename=fn, width=3072, height=2304, pointsize=40, res=300 )
       lp = map_simple( toplot=obs[ ii, c("plon","plat", vn) ], plotarea=B, resol=1, theta=15, filterdistances=7.5, vn=vn, annot=paste("Residuals", y), er=er )
       print(lp)
      dev.off()
      print(fn)
    }
  }

