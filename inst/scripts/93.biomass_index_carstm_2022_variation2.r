


# -------------------------------------------------
# Snow crab --- Areal unit modelling Hurdle / Delta model  
# combination of three models via posterior simulation
# 1. Poisson on positive valued numbers offset by swept are
# 2. Meansize in space and time 
# 3  Presence-absence
# the convolution of all three after simulation is called a Hurdle or Delta model
# -------------------------------------------------
 

# this is copied from 03.biomass_index_carstm.r 
# but stripped down with modifications for comparisons demanded by 2023 CSAS review
# actual comparisons begin toward the bottom of this file (lines 500+). 
# This next section is generating the fitted results of the deletion


# -------------------------------------------------
# Part 1 -- construct basic parameter list defining the main characteristics of the study

  source( file.path( code_root, "bio_startup.R" )  )
   
  require(bio.snowcrab)   # loadfunctions("bio.snowcrab") 

  year.assessment = 2022
  
  yrs = 1999:year.assessment
  spec_bio = bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 )
  snowcrab_filter_class = "fb"     # fishable biomass (including soft-shelled )  "m.mat" "f.mat" "imm"

  
  # key name 
  carstm_model_label= paste( "1999_2022_variation2", snowcrab_filter_class, sep="_" )

  # params for number
  pN = snowcrab_parameters(
    project_class="carstm",
    yrs=yrs,   
    areal_units_type="tesselation",
    family = "nbinomial",
    carstm_model_label= carstm_model_label,  
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
    carstm_model_label= carstm_model_label,  
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
    carstm_model_label= carstm_model_label,  
    selection = list(
      type = "presence_absence",
      biologicals=list( spec_bio=spec_bio ),
      biologicals_using_snowcrab_filter_class=snowcrab_filter_class
    )
  )
  
  # use what was defined in the main script
  sppoly=areal_units( p=pN )

  additional_features = snowcrab_mapping_features(pN)  # for mapping below
 
  tmap_mode("plot")
   



# -- identify area to drop

# filter out candidate time period
# this scenario tries to mimic the effect of the incomplete sampling in 2022
# using stations that were within 5 km of the stations completed in 2022 

  fnout = file.path( pN$modeldir, pN$carstm_model_label, "unsampled_polygon.dat" )
  
  if (!file.exists(fnout)) {
    project.library( "aegis", "aegis.polygons" )
    RLibrary( "raster"  )
    message("FIXE ME::: deprecated libs, use sf/stars")
    # plot background map and using mouse interaction define region:
    # left mouse click to register, right mouse click to finish (or [Esc] is using Rstudio)
    maps::map( database="worldHires", regions=c("Canada", "USA") ,
       xlim=c(-66, -55 ), ylim=c(42, 47), fill=FALSE, plot=TRUE )
    set = snowcrab.db( DS="set.clean")
    points(lat~lon, set[ set$yr==2021,], pch=19, col="red" ) 
    points(lat~lon, set[ set$yr==2022,], pch=22, cex=2 ) 
    X = locator(type="o" )
    X = as.data.frame( X)
    colnames(X) = c("lon", "lat")
    X = rbind(X, X[1,])
    lines (X)
    write.csv(X, fnout)
  }
  
  X = read.csv(fnout) # outline of unsampled region
  
  set = snowcrab.db( DS="set.clean")
  set$id = paste(set$trip, set$set, sep=".")
  u = which(set$yr == 2022)

  Z = set[set$yr==2015, ]  # year to test removal
  a = which(point.in.polygon(Z$lon, Z$lat, X$lon, X$lat) != 0 )
  todrop = Z$id[a]
  tokeep = Z$id[-a] 
   
   
  if (0) {
      # plot these locations
      require(ggplot2)
      crs_domain = st_crs( "+proj=utm +ellps=WGS84 +zone=20 +units=km"  )
      domain = st_union( st_as_sf(Z, coords=c("lon", "lat")) )
      st_crs(domain) = st_crs(projection_proj4string("lonlat_wgs84"))
      domain = st_transform(domain, crs_domain )
        bb = st_bbox(domain)

      isobaths = c( 50, 100, 150, 200, 250, 300, 350, 400 )
      isobs = aegis.bathymetry::isobath_db( depths=isobaths, project_to=crs_domain )
      isobs = st_intersection(isobs, domain)
      coastline = st_transform( polygons_rnaturalearth(countries=c("United States of America", "Canada"),
        xlim=c(-80,-40), ylim=c(38, 60)), st_crs(crs_domain) )   
      
      plt = ggplot() +
        geom_point(data=Z, aes(x=plon, y=plat, colour="gray", alpha=1.0) ) +
        geom_point(data=Z[a,], aes(x=plon, y=plat, colour="darkorange", alpha=1.0) ) +
        coord_sf(xlim = c(bb[c("xmin", "xmax")]), ylim = c(bb[c("ymin", "ymax")]) ) +
        coord_fixed()
    
    dev.new(width=14, height=8, pointsize=20)
    print(plt)  
    
  }


# ------------------------------------------------
# Part 2 -- spatiotemporal statistical model

  if ( spatiotemporal_model ) {

    # total numbers
    sppoly = areal_units( p=pN )
    M = snowcrab.db( p=pN, DS="carstm_inputs", sppoly=sppoly  )  # will redo if not found
    # can also just copy datafile: carstm_inputs_snowcrab~tesselation~1~snowcrab~24~1~none~snowcrab_managementareas.rdata
    # into working directory to speed things up

    # variation2: drop number, size info:: this forces prediction of these points 
    Z_M = match(todrop, M$id) # in yr selected (n=83)
    Z_N = match(tokeep, M$id) # in yr selected (n=330)

    sd(M$totwgt[Z_M] / M$data_offset[Z_M], na.rm=T) #    [1] 5.4
    sd(M$totwgt[Z_N] / M$data_offset[Z_N], na.rm=T)  #    [1] 3.8

    sd(M$totwgt[Z_M] / M$data_offset[Z_M], na.rm=T) #    [1] 9.65 .. se= 9.65 /sqrt(83) = 1.06
    sd(M$totwgt[Z_N] / M$data_offset[Z_N], na.rm=T)  #    [1] 8.46 .. se=8.46/sqrt(330) = 0.81
 
    sd(M$totwgt[Z_M] / M$data_offset[Z_M], na.rm=T) / sqrt(83) #    [1] 9.65 .. se= 9.65 /sqrt(83) = 1.06
    sd(M$totwgt[Z_N] / M$data_offset[Z_N], na.rm=T) /sqrt(330) #    [1] 8.46 .. se=8.46/sqrt(330) = 0.81
 
    
    io = which(M$tag=="observations")
    ip = which(M$tag=="predictions")
    iq = unique( c( which( M$totno > 0), ip ) )
    iw = unique( c( which( M$totno > 5), ip ) )  # need a good sample to estimate mean size
    
    Mdropped = M[ Z_M, ] 
    
    M[ Z_M, "totno"] = NA
    M[ Z_M, "pa"] = NA
    M[ Z_M, "meansize"] = NA
 

    # number 
    fit = NULL; gc()
    fit = carstm_model( p=pN, data=M[ iq, ], sppoly=sppoly, 
      posterior_simulations_to_retain="predictions", improve.hyperparam.estimates=TRUE
    )

    MO = M[ iq, ]  
    ii = match( todrop, MO$id  )
    jj = match( todrop, Mdropped$id  )
    observed = Mdropped[jj, "totno"] / Mdropped[jj, "data_offset"] 
    fitted = fit$summary.fitted.values[["mean"]] [ii]
    plot( fitted ~ observed )
    abline(0,1)

    cor( fitted, observed, use="pairwise.complete.obs" )
    cor( fitted, observed, use="pairwise.complete.obs", "spearman" )

    # mean size
    fit = NULL; gc()
    fit = carstm_model( p=pW, data=M[ iw, ], sppoly = sppoly, 
      posterior_simulations_to_retain="predictions", improve.hyperparam.estimates=TRUE,
      control.inla = list( strategy="laplace", int.strategy="eb" )
    ) 

    MO = M[ iw, ]  
    ii = match( todrop, MO$id  )
    jj = match( todrop, Mdropped$id  )
    observed = Mdropped[jj, "meansize"]  
    fitted = fit$summary.fitted.values[["mean"]] [ii]
    plot( fitted ~ observed )
    abline(0,1)
    cor( fitted, observed, use="pairwise.complete.obs" )
    cor( fitted, observed, use="pairwise.complete.obs", "spearman" )


    # model pa using all data
    fit = NULL; gc()
    fit = carstm_model( p=pH, data=M, sppoly=sppoly, 
      posterior_simulations_to_retain="predictions", improve.hyperparam.estimates=TRUE,
      # control.family=list(control.link=list(model="logit")),  # default
      control.inla = list( strategy="laplace", int.strategy="eb" )
    )

    MO = M 
    ii = match( todrop, MO$id  )
    jj = match( todrop, Mdropped$id  )
    observed = Mdropped[jj, "pa"]  
    fitted = fit$summary.fitted.values[["mean"]] [ii]
    plot( fitted ~ jitter(observed) )
     abline(0,1)
    cor( fitted, observed, use="pairwise.complete.obs" )
    cor( fitted, observed, use="pairwise.complete.obs", "spearman" )


    ### continue as in 03_biomass_index_carstm.R .. alter save locations

 
  }  # end spatiotemporal model




# ----------------------
# Part 3: assimilation of models


  assimilate_numbers_and_size = TRUE

  if (assimilate_numbers_and_size ) {

    ### continue as in 03_biomass_index_carstm.R  .. alter save locations

  }  # end assimilate size and numbers


  
# end
