


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
  runlabel= paste( "1999_2022_variation2", snowcrab_filter_class, sep="_" )

  # params for number
  pN = snowcrab_parameters(
    project_class="carstm",
    yrs=yrs,   
    areal_units_type="tesselation",
    family = "nbinomial",
    carstm_model_label= runlabel,  
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



    for ( selection in c("number", "meansize", "presence_absence") ) {
      
      # map all :
      if ( selection=="number" ) {
        p=pN
        res = carstm_model( p=p, DS="carstm_modelled_summary",  sppoly = sppoly ) # to load currently saved results
        outputdir = file.path( p$modeldir, p$carstm_model_label, "predicted.numerical.densities" )
        ylab = "Number"
        if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )
        fn_root_prefix = "Predicted_numerical_abundance"
        fn_root =  "Predicted_numerical_abundance_persistent_spatial_effect" 
        outfilename = file.path( outputdir, paste(fn_root, "png", sep=".") )
        title= paste( snowcrab_filter_class, "Number; no./m^2"  )
      }
      if ( selection=="meansize") {
        p=pW
        res = carstm_model( p=p, DS="carstm_modelled_summary",  sppoly = sppoly ) # to load currently saved results
        ylab = "Mean weight"
        outputdir = file.path( p$modeldir, p$carstm_model_label, "predicted.meansize" )
        if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )
        fn_root_prefix = "Predicted_meansize"
        fn_root =  "Predicted_meansize_persistent_spatial_effect" 
        outfilename = file.path( outputdir, paste(fn_root, "png", sep=".") )
        title= paste( snowcrab_filter_class, "Mean weight; kg" ) 
      }
      if ( selection=="presence_absence" ) {
        p=pH
        res = carstm_model( p=p, DS="carstm_modelled_summary",  sppoly = sppoly ) # to load currently saved results
        ylab = "Probability"
        outputdir = file.path( p$modeldir, p$carstm_model_label, "predicted.presence_absence" )
        if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )
        fn_root_prefix = "Predicted_presence_absence"
        fn_root =  "Predicted_presence_absence_persistent_spatial_effect" 
        outfilename = file.path( outputdir, paste(fn_root, "png", sep=".") )
        title= paste( snowcrab_filter_class, "Probability")  
      }


      if (0) {
           # choose:  
          #p = pN
          #p = pW
          #p = pH

        # extract results
        fit = carstm_model( p=p, DS="carstm_modelled_fit",  sppoly = sppoly )  # extract currently saved model fit
        fit$summary$dic$dic
        fit$summary$dic$p.eff
        plot(fit)
        plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
        plot( fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
        plot( fit$marginals.hyperpar$"Phi for space_time", type="l")  # posterior distribution of phi nonspatial dominates
        plot( fit$marginals.hyperpar$"Precision for space_time", type="l")
        plot( fit$marginals.hyperpar$"Precision for setno", type="l")
        fit = NULL

        res = carstm_model( p=p, DS="carstm_modelled_summary",  sppoly = sppoly ) # to load currently saved results

        # quick plots
        vn=c( "random", "space", "re" )
        vn=c( "random", "spacetime", "re" )
        vn="predictions"  # numerical density (km^-2)

        tmatch= as.character(year.assessment)

        carstm_map(  res=res, vn=vn, tmatch=tmatch, 
            sppoly = sppoly, 
            palette="-RdYlBu",
            plot_elements=c(  "compass", "scale_bar", "legend" ),
            additional_features=additional_features,
            title =paste( vn, paste0(tmatch, collapse="-"), "no/m^2"  )
        )

      }

      vn = c( "random", "space", "re" ) 
      toplot = carstm_results_unpack( res, vn )
      brks = pretty(  quantile(toplot[,"mean"], probs=c(0,0.975), na.rm=TRUE )  )


      plt = carstm_map(  res=res, vn=vn, 
        sppoly = sppoly, 
        breaks = brks,
        palette="-RdYlBu",
        plot_elements="",
        # c(  "compass", "scale_bar", "legend" ),
        #        title= title
        additional_features=additional_features,
        outfilename=outfilename
      )  
      plt
    

      vn="predictions"
      toplot = carstm_results_unpack( res, vn )
      brks = pretty(  quantile(toplot[,,"mean"], probs=c(0,0.975), na.rm=TRUE )  )

      for (y in res$time ){
        tmatch = as.character(y)
        fn_root = paste(fn_root_prefix, paste0(tmatch, collapse="-"), sep="_")
        outfilename = file.path( outputdir, paste(fn_root, "png", sep=".") )

        plt = carstm_map(  res=res, vn=vn, tmatch=tmatch,
          sppoly = sppoly, 
          breaks =brks,
          palette="-RdYlBu",
          plot_elements="",
          # plot_elements=c(   "compass", "scale_bar", "legend" ),
          additional_features=additional_features,
          title = tmatch, 
#          title=paste(fn_root_prefix, snowcrab_filter_class,  paste0(tmatch, collapse="-") )
          outfilename=outfilename
        )
        plt
        print(outfilename)
      
      }

      # plots with 95% PI
      oeffdir = file.path( outputdir, fn_root_prefix, "effects" )
      if ( !file.exists(oeffdir)) dir.create( oeffdir, recursive=TRUE, showWarnings=FALSE )

      (fn = file.path( oeffdir, "time.png"))
      png( filename=fn, width=1024, height=1024, pointsize=12, res=196 )
        carstm_plotxy( res, vn=c( "res", "random", "time" ), 
          type="b",  xlab="Year", ylab=ylab, h=0, cex=1.25, cex.axis=1.25, cex.lab=1.25   )
      dev.off()

      (fn = file.path( oeffdir, "cyclic.png"))
      png( filename=fn, width=1024, height=1024, pointsize=12, res=196 )
        carstm_plotxy( res, vn=c( "res", "random", "cyclic" ), 
          type="b", col="slategray", pch=19, lty=1, lwd=2.5,  
          xlab="Season", ylab=ylab, cex=1.25, cex.axis=1.25, cex.lab=1.25   )
      dev.off()


      (fn = file.path( oeffdir, "temperature.png"))
      png( filename=fn, width=1024, height=1024, pointsize=12, res=196 )
        carstm_plotxy( res, vn=c( "res", "random", "inla.group(t, method = \"quantile\", n = 11)" ), 
          type="b", col="slategray", pch=19, lty=1, lwd=2.5 ,
          xlab="Bottom temperature (degrees Celsius)", ylab=ylab, cex=1.25, cex.axis=1.25, cex.lab=1.25 )
      dev.off()


      (fn = file.path( oeffdir, "pca1.png"))
      png( filename=fn, width=1024, height=1024, pointsize=12, res=196 )
        carstm_plotxy( res, vn=c( "res", "random", "inla.group(pca1, method = \"quantile\", n = 11)" ), 
          type="b", col="slategray", pch=19, lty=1, lwd=2.5 ,
          xlab="PCA1", ylab=ylab, cex=1.25, cex.axis=1.25, cex.lab=1.25 )
      dev.off()

      (fn = file.path( oeffdir, "pca2.png"))
      png( filename=fn, width=1024, height=1024, pointsize=12, res=196 )
        carstm_plotxy( res, vn=c( "res", "random", "inla.group(pca2, method = \"quantile\", n = 11)" ), 
          type="b", col="slategray", pch=19, lty=1, lwd=2.5 ,
          xlab="PCA2", ylab=ylab, cex=1.25, cex.axis=1.25, cex.lab=1.25 )
      dev.off()

      (fn = file.path( oeffdir, "depth.png"))
      png( filename=fn, width=1024, height=1024, pointsize=12, res=196 )
        carstm_plotxy( res, vn=c( "res", "random", "inla.group(z, method = \"quantile\", n = 11)" ), 
          type="b", col="slategray", pch=19, lty=1, lwd=2.5  ,
          xlab="Depth (m)", ylab=ylab, cex=1.25, cex.axis=1.25, cex.lab=1.25 )
      dev.off()


      # (fn = file.path( outputdir, "substrate.png"))
      # png( filename=fn, width=1024, height=1024, pointsize=12, res=196 )
      #   carstm_plotxy( res, vn=c( "res", "random", "inla.group(substrate.grainsize, method = \"quantile\", n = 11)" ), 
      #     type="b", col="slategray", pch=19, lty=1, lwd=2.5  ,
      #     xlab="Substrate grain size (mm)", ylab=ylab, cex=1.25, cex.axis=1.25, cex.lab=1.25 )
      # dev.off()

      if (0) {
        fit = carstm_model( p=pW, DS="carstm_modelled_fit",  sppoly = sppoly ) # to load currently saved results
    
        plot( fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
        plot( fit$marginals.hyperpar$"Phi for space_time", type="l")  # posterior distribution of phi nonspatial dominates
        plot( fit$marginals.hyperpar$"Precision for space_time", type="l")
        plot( fit$marginals.hyperpar$"Precision for setno", type="l")
      }
    }

    } #end each seletion

  }  # end spatiotemporal model




# ----------------------
# Part 3: assimilation of models


  assimilate_numbers_and_size = TRUE

  if (assimilate_numbers_and_size ) {

    # wgts_max = 1.1 # kg, hard upper limit
    # N_max = NULL
    # #  quantile( M$totno[ipositive]/M$data_offset[ipositive], probs=0.95, na.rm=TRUE )  
    
    # posterior sims 
  
    sims = carstm_posterior_simulations( pN=pN, pW=pW, pH=pH,  pa_threshold=0.05, qmax=0.99 )
    sims = sims  / 10^6 # 10^6 kg -> kt;; kt/km^2
   
    SM = aggregate_simulations( 
      sims=sims, 
      sppoly=sppoly, 
      fn=carstm_filenames( pN, returnvalue="filename", fn="aggregated_timeseries" ), 
      yrs=pN$yrs, 
      method="sum", 
      redo=TRUE 
    ) 

    if (0) {
      # to compute habitat prob
      sims = carstm_posterior_simulations( pH=pH,  pa_threshold=0.05, qmax=0.99 )
      SM = aggregate_simulations( 
        sims=sims, 
        sppoly=sppoly, 
        fn=carstm_filenames( pN, returnvalue="filename", fn="aggregated_timeseries" ), 
        yrs=pN$yrs, 
        method="mean", 
        redo=TRUE 
      ) 
      outputdir = file.path( carstm_filenames( pN, returnvalue="output_directory"), "aggregated_habitat_timeseries" )
      RES= SM$RES
 
    }      
    
    RES= SM$RES
    # RES = aggregate_simulations( fn=carstm_filenames( pN, returnvalue="filename", fn="aggregated_timeseries" ) )$RES

    # note: using pN, even though this is biomass 
    
    outputdir = file.path( carstm_filenames( pN, returnvalue="output_directory"), "aggregated_biomass_timeseries" )

    if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )


    ( fn = file.path( outputdir, "cfa_all.png") )
    png( filename=fn, width=3072, height=2304, pointsize=12, res=300 )
      plot( cfaall ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Biomass index (kt)", xlab="")
      lines( cfaall_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
      lines( cfaall_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    dev.off()

    ( fn = file.path( outputdir, "cfa_south.png") )
    png( filename=fn, width=3072, height=2304, pointsize=12, res=300 )
      plot( cfasouth ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Biomass index (kt)", xlab="")
      lines( cfasouth_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
      lines( cfasouth_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    dev.off()

    ( fn = file.path( outputdir, "cfa_north.png") )
    png( filename=fn, width=3072, height=2304, pointsize=12, res=300 )
      plot( cfanorth ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Biomass index (kt)", xlab="")
      lines( cfanorth_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
      lines( cfanorth_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    dev.off()

    ( fn = file.path( outputdir, "cfa_4x.png") )
    png( filename=fn, width=3072, height=2304, pointsize=12, res=300 )
      plot( cfa4x ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Biomass index (kt)", xlab="")
      lines( cfa4x_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
      lines( cfa4x_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    dev.off()


    regions = c("cfanorth", "cfasouth",  "cfa4x" )
    region_label = c("N-ENS", "S-ENS", "4X")
 
    a= cbind( "cfanorth", RES[,c("yrs", "cfanorth", "cfanorth_lb", "cfanorth_ub")] )
    b= cbind( "cfasouth", RES[,c("yrs", "cfasouth", "cfasouth_lb", "cfasouth_ub")] )
    c= cbind( "cfa4x", RES[,c("yrs", "cfa4x", "cfa4x_lb", "cfa4x_ub")] )
    names(a) = names(b) = names(c) = c("region", "year", "mean", "lb", "ub")
    tdb = rbind(a, b, c)

    tdb$region = factor(tdb$region, levels=regions, labels =region_label)
    tdb = tdb[(which(!is.na(tdb$region))), ]
   
    fn = file.path( outputdir, "biomass_M0.png" )
    
    require(ggplot2)
    library(ggbreak) 

    color_map = c("#E69F00", "#56B4E9",  "#CC79A7" )  
    
    out = ggplot(tdb, aes(x=year, y=mean, fill=region, colour=region)) +
      geom_line( alpha=0.9, linewidth=1.2 ) +
      geom_point(aes(shape=region), size=3, alpha=0.7 ) +
      geom_errorbar(aes(ymin=lb,ymax=ub), linewidth=0.8, alpha=0.8, width=0.3)  +
      labs(x=NULL, y=NULL) +
      # labs(x="Year", y="Biomass index (kt)", size = rel(1.5)) +
      scale_colour_manual(values=color_map) +
      scale_fill_manual(values=color_map) +
      scale_shape_manual(values = c(15, 17, 19)) +
      theme_light( base_size = 22) + 
      theme( legend.position=c(0.75, 0.9), legend.title=element_blank()) +
      scale_y_break(c(14, 28), scales = 1)
      
      # scale_y_continuous( limits=c(0, 300) )  
      ggsave(filename=fn, plot=out, device="png", width=12, height = 8)
 


    # map it ..mean density
    sppoly = areal_units( p=pN )  # to reload

    vn = paste("biomass", "predicted", sep=".")

    outputdir = file.path( carstm_filenames( pN, returnvalue="output_directory"), "predicted_biomass_densitites" )

    if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

    B = apply( sims, c(1,2), mean ) 
    B[ which(!is.finite(B)) ] = NA

    brks = pretty( log10( quantile( B[], probs=c(0.05, 0.95), na.rm=TRUE )* 10^6)  )
  
    additional_features = snowcrab_mapping_features(pN)  # for mapping below

    for (i in 1:length(pN$yrs) ) {
      y = as.character( pN$yrs[i] )
      sppoly[,vn] = log10( B[,y]* 10^6 )
      outfilename = file.path( outputdir , paste( "biomass", y, "png", sep=".") )
      plt =  carstm_map(  sppoly=sppoly, vn=vn,
          breaks=brks,
          additional_features=additional_features,
          title=y,
          # title=paste( "log_10( Predicted biomass density; kg/km^2 )", y ),
          palette="-RdYlBu",
          plot_elements=c( "compass", "scale_bar", "legend" ), 
          outfilename=outfilename
      )
      plt
      
    }
  
  }  # end assimilate size and numbers


 
  
# prep data for discrete version
# Rdata files are ready load them through julia and model
# fishery_model_data_inputs( year.assessment=year.assessment,  type="biomass_dynamics", for_julia=TRUE  )


# end
