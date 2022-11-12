


# -------------------------------------------------
# Snow crab --- Areal unit modelling Hurdle / Delta model  
# combination of three models via posterior simulation
# 1. Poisson on positive valued numbers offset by swept area
# 2. Meansize in space and time 
# 3  Presence-absence
# the convolution of all three after simulation is called a Hurdle or Delta model
# -------------------------------------------------
 

# TODO::: move plotting calls to self-contained functions:


# -------------------------------------------------
# Part 1 -- construct basic parameter list defining the main characteristics of the study

source( file.path( code_root, "bio_startup.R" )  )
 
require(bio.snowcrab)   # loadfunctions("bio.snowcrab") 
require(tmap)

year.assessment = 2021
yrs = 1999:year.assessment
spec_bio = bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 )

for (snowcrab_filter_class in c( "M0", "M1", "M2", "M3", "M4", "f.mat") ) {

      # snowcrab_filter_class = "M0"     # "fishable biomass" (excluding soft-shelled )
      # snowcrab_filter_class = "M1"     # some fraction expected to enter M0 next year
      # snowcrab_filter_class = "M2"     # some fraction expected to enter M1 next year / M0 in 2 years
      # snowcrab_filter_class = "M3"     # some fraction expected to enter M2 next year / M0 in 3 years
      # snowcrab_filter_class = "M4"     # some fraction expected to enter M3 next year / M0 in 4 years
      # snowcrab_filter_class = "f.mat"
      
      # snowcrab_filter_class = "imm"  # note poisson will not work due to var inflation .. nbinomial is a better choice 
      # snowcrab_filter_class = "m.mat"
    
      
      runlabel= paste( "1999_present", snowcrab_filter_class, sep="_" )
    
      # poisson works too but variance is not exactly poisson (higher than mean)
      Nfamily = switch( snowcrab_filter_class, 
        M0 = "poisson",  # to match index for historical approach, nbinomial is less dampened
        M1 = "nbinomial",
        M2 = "nbinomial",
        M3 = "nbinomial",
        M4 = "nbinomial",
        f.mat = "nbinomial"
      )

      # params for number
      pN = snowcrab_parameters(
        project_class="carstm",
        yrs=yrs,   
        areal_units_type="tesselation",
        family = Nfamily,  
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
    
      
      if (areal_units) {
        # polygon structure:: create if not yet made
        # for (au in c("cfanorth", "cfasouth", "cfa4x", "cfaall" )) plot(polygon_managementareas( species="snowcrab", au))
        xydata = snowcrab.db( p=pN, DS="areal_units_input", redo=TRUE )
        xydata = snowcrab.db( p=pN, DS="areal_units_input" )
        
        sppoly = areal_units( p=pN, xydata=xydata[ which(xydata$yr %in% pN$yrs), ], redo=TRUE, verbose=TRUE )  # create constrained polygons with neighbourhood as an attribute
        sppoly=areal_units( p=pN )
      
        plot(sppoly["npts"])
    
        sppoly$dummyvar = ""
        xydata = st_as_sf( xydata, coords=c("lon","lat") )
        st_crs(xydata) = st_crs( projection_proj4string("lonlat_wgs84") )
      
        tmap_mode("plot")
    
        require("tmap")
        
        tmout = 
          tm_shape(sppoly) +
            tm_borders(col = "slategray", alpha = 0.5, lwd = 0.5) + 
            tm_shape( xydata ) + tm_sf() +
            additional_features +
            tm_compass(position = c("right", "TOP"), size = 1.5) +
            tm_scale_bar(position = c("RIGHT", "BOTTOM"), width =0.1, text.size = 0.5) +
            tm_layout(frame = FALSE, scale = 2) +
            tm_shape( st_transform(polygons_rnaturalearth(), st_crs(sppoly) )) + 
            tm_borders(col = "slategray", alpha = 0.5, lwd = 0.5)
    
        dev.new(width=14, height=8, pointsize=20)
        tmout
    
      }
    
    
      if (temperature_figures) {
      
        # area-specific figures
        figure_area_based_extraction_from_carstm(DS="temperature" )  # can only do done once we have an sppoly for snow crab
      
        # full domain:
        # default paramerters (copied from 03_temperature_carstm.R )
        require(aegis.temperature)
        params = list( 
          temperature = temperature_parameters( 
            project_class="carstm", 
            yrs=1970:year.assessment, 
            carstm_model_label="1970_present"
          ) 
        )
    
        tss = aegis_lookup(  
          parameters=params["temperature"], 
          LOCS=expand.grid( AUID=sppoly$AUID, timestamp= yrs + 0.75 ), LOCS_AU=sppoly, 
          project_class="carstm", output_format="areal_units", 
          variable_name=list( "predictions" ), statvars=c("mean", "sd"), space_resolution=pN$pres,
          returntype = "data.table"
        ) 
      
      }
    
    
      sppoly=areal_units( p=pN )
      additional_features = snowcrab_features_tmap(pN)  # for mapping below
      
      M = snowcrab.db( p=pN, DS="carstm_inputs", sppoly=sppoly, redo=TRUE )  # will redo if not found
    
      
      
      
    
    # ------------------------------------------------
    # Part 2 -- spatiotemporal statistical model
      spatiotemporal_model = TRUE
      if ( spatiotemporal_model ) {
    
        # total numbers
        sppoly = areal_units( p=pN )
        M = snowcrab.db( p=pN, DS="carstm_inputs", sppoly=sppoly  )  # will redo if not found
      
        io = which(M$tag=="observations")
        ip = which(M$tag=="predictions")
        iq = unique( c( which( M$totno > 0), ip ) )
        iw = unique( c( which( M$totno > 5), ip ) )  # need a good sample to estimate mean size
    
        
        # number 
        fit = NULL; gc()
        fit = carstm_model( p=pN, data=M[ iq, ], sppoly=sppoly, 
          posterior_simulations_to_retain="predictions", improve.hyperparam.estimates=TRUE
        )
    
        # mean size
        fit = NULL; gc()
        fit = carstm_model( p=pW, data=M[ iw, ], sppoly = sppoly, 
          posterior_simulations_to_retain="predictions", improve.hyperparam.estimates=TRUE,
          control.inla = list( strategy="laplace", int.strategy="eb" )
        ) 
    
        # model pa using all data
        fit = NULL; gc()
        fit = carstm_model( p=pH, data=M, sppoly=sppoly, 
          posterior_simulations_to_retain="predictions", improve.hyperparam.estimates=TRUE,
          # control.family=list(control.link=list(model="logit")),  # default
          control.inla = list( strategy="laplace", int.strategy="eb" )
        )
    
    
        if (0) {
          # examine fits: 
          # choose:  
          p = pN
          p = pW
          p = pH
    
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
    
          # quick plots
          res = carstm_model( p=p, DS="carstm_modelled_summary",  sppoly = sppoly ) # to load currently saved results
          vn=c( "random", "space", "combined" )
          vn=c( "random", "spacetime", "combined" )
          vn="predictions"  # numerical density (km^-2)
    
          tmatch= as.character(year.assessment)
    
          carstm_map(  res=res, vn=vn, tmatch=tmatch, 
              sppoly = sppoly, 
              palette="-RdYlBu",
              plot_elements=c(  "compass", "scale_bar", "legend" ),
              additional_features=additional_features,
              title =paste( vn, paste0(tmatch, collapse="-"), "no/m^2"  )
          )
    
          # express habitat relationships 
          p = pH
          res = carstm_model( p=p, DS="carstm_modelled_summary",  sppoly = sppoly ) # to load currently saved results
    
          o = carstm_2D_effects_probability( 
            res,
            xvar = "inla.group(t, method = \"quantile\", n = 11)",  
            yvar = "inla.group(z, method = \"quantile\", n = 11)", 
            xgrid = seq( -1, 10.5, by=0.5),
            ygrid = seq( 25, 350, by=25),
            xslice = 4,
            yslice = -200,
            nx=200, ny=200,
            theta = 140,
            phi = 15
          )
        
          # use a larger domain than sppoly for the following estimate:
          # sppoly is constrained to sampled locations, and so missing a lot of the inshore areas
      
          x11()
          crs_plot = st_crs( sppoly )
          domain = polygon_managementareas( species="maritimes" )
          domain = st_transform( domain, crs_plot )
          data_mask = st_union( sppoly[which(sppoly$filter==1),1] ) 
          # all = st_union( domain, data_mask )
          nearshore = st_cast( st_difference( domain, data_mask ), "POLYGON")[1]
          domain_new = st_union( data_mask, nearshore )
            
    
          o = carstm_optimal_habitat( 
            res = res,
            xvar = "inla.group(t, method = \"quantile\", n = 11)",  
            yvar = "inla.group(z, method = \"quantile\", n = 11)",
            depths=switch( snowcrab_filter_class, 
              fb = c(100, 350),
              imm = c( 160, 350),
              f.mat = c(100, 160),
              m.mat = c(160, 300)
            ),
            probability_limit = 0.25,
            nsims = 100,
            domain=domain_new 
          ) 
          
          dev.new();
          print( o["depth_plot"] )
      
            if (0) {
              u = readRDS('/home/jae/tmp/temp_depth_habitat.RDS')
              dev.new()
              plot( habitat~yr, u, type="b", ylim=c(0.1, 0.33))
              lines( habitat_lb~yr, u)
              lines( habitat_ub~yr, u)
              abline(v=1993)
              abline(v=2012)
            
              dev.new()
              plot( habitat_sa~yr, u, type="b" )
              lines( habitat_sa_lb~yr, u)
              lines( habitat_sa_ub~yr, u)
              abline(v=1993)
              abline(v=2012)
      
              ll = loess(habitat~yr, u, span=0.25 )
              pp = predict( ll, u )
              lines(pp ~ u$yr)
      
            }
      
    
          outputdir = file.path( p$modeldir, p$carstm_model_label )
          fn_optimal = file.path( outputdir, "optimal_habitat_temperature_depth_effect.RDS" )
          saveRDS( o, file=fn_optimal, compress=FALSE )
          o = readRDS(fn_optimal)
    
          library(ggplot2)
    
          dev.new(width=14, height=8, pointsize=20)
          ggplot( o[["temperature_depth"]], aes(yr, habitat ) ) +
            geom_ribbon(aes(ymin=habitat_lb, max=habitat_ub), alpha=0.2, colour=NA) +
            geom_line() +
            labs(x="Year", y="Habitat probabtility", size = rel(1.5)) +
            # scale_y_continuous( limits=c(0, 300) )  
            theme_light( base_size = 22 ) 
          
    
          dev.new(width=14, height=8, pointsize=20)
          ggplot( o[["temperature_depth"]], aes(yr, habitat_sa ) ) +
            geom_ribbon(aes(ymin=habitat_sa_lb, max=habitat_sa_ub), alpha=0.2, colour=NA) +
            geom_line() +
            labs(x="Year", y=bquote("Habitat surface area;" ~ km^2), size = rel(1.5)) +
            # scale_y_continuous( limits=c(0, 300) )  
            theme_light( base_size = 22 ) 
    
            
        }
    
        # generic maps and figures:
    
        for ( data_class in c( "number", "weight", "habitat") ) {
        
          if ( data_class == "number" ) {
            p=pN
            outputdir = file.path( p$modeldir, p$carstm_model_label, "predicted.numerical.densitites" )
            ylab = "Number"
            fn_root_prefix = "Predicted_numerical_abundance"
            fn_root =  "Predicted_numerical_abundance_persistent_spatial_effect" 
            title= paste( snowcrab_filter_class, "Number; no./m^2"  )
          }
          if ( data_class == "meansize" ) {
            p=pW
            ylab = "Mean weight"
            outputdir = file.path( p$modeldir, p$carstm_model_label, "predicted.meansize" )
            fn_root_prefix = "Predicted_meansize"
            fn_root =  "Predicted_meansize_persistent_spatial_effect" 
            title= paste( snowcrab_filter_class, "Mean weight; kg" ) 
          }
          if ( data_class == "presence_absence" ) {
            p=pH
            ylab = "Probability"
            outputdir = file.path( p$modeldir, p$carstm_model_label, "predicted.presence_absence" )
            fn_root_prefix = "Predicted_presence_absence"
            fn_root =  "Predicted_presence_absence_persistent_spatial_effect" 
            title= paste( snowcrab_filter_class, "Probability")  
          }
                
          res = carstm_model( p=p, DS="carstm_modelled_summary",  sppoly = sppoly ) # to load currently saved results
          
          if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )
          outfilename = file.path( outputdir, paste(fn_root, "png", sep=".") )
    
    
          # spatial effects
          vn = c( "random", "space", "combined" ) 
          toplot = carstm_results_unpack( res, vn )
          brks = pretty(  quantile(toplot[,"mean"], probs=c(0,0.975), na.rm=TRUE )  )
    
          tmout = carstm_map(  res=res, vn=vn, 
            sppoly = sppoly, 
            breaks = brks,
            palette="-RdYlBu",
            plot_elements="",
            # c(  "compass", "scale_bar", "legend" ),
            #        title= title
            additional_features=additional_features,
            outfilename=outfilename
          )  
          tmout
        
          # predictions
          vn="predictions"
          toplot = carstm_results_unpack( res, vn )
          brks = pretty(  quantile(toplot[,,"mean"], probs=c(0,0.975), na.rm=TRUE )  )
    
          for (y in res$time ){
            tmatch = as.character(y)
            fn_root = paste(fn_root_prefix, paste0(tmatch, collapse="-"), sep="_")
            outfilename = file.path( outputdir, paste(fn_root, "png", sep=".") )
    
            tmout = carstm_map(  res=res, vn=vn, tmatch=tmatch,
              sppoly = sppoly, 
              breaks =brks,
              palette="-RdYlBu",
              plot_elements="",
              # plot_elements=c(   "compass", "scale_bar", "legend" ),
              additional_features=additional_features,
    #          title=paste(fn_root_prefix, snowcrab_filter_class,  paste0(tmatch, collapse="-") )
              outfilename=outfilename
            )
            tmout
            print(outfilename)
          
          }
    
          # plots with 95% PI
          outputdir = file.path( p$modeldir, p$carstm_model_label, "effects" )
          if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )
    
          (fn = file.path( outputdir, "time.png"))
          png( filename=fn, width=1024, height=1024, pointsize=12, res=196 )
            carstm_plotxy( res, vn=c( "res", "random", "time" ), 
              type="b",  xlab="Year", ylab=ylab, h=0, cex=1.25, cex.axis=1.25, cex.lab=1.25   )
          dev.off()
    
          (fn = file.path( outputdir, "cyclic.png"))
          png( filename=fn, width=1024, height=1024, pointsize=12, res=196 )
            carstm_plotxy( res, vn=c( "res", "random", "cyclic" ), 
              type="b", col="slategray", pch=19, lty=1, lwd=2.5,  
              xlab="Season", ylab=ylab, cex=1.25, cex.axis=1.25, cex.lab=1.25   )
          dev.off()
    
    
          (fn = file.path( outputdir, "temperature.png"))
          png( filename=fn, width=1024, height=1024, pointsize=12, res=196 )
            carstm_plotxy( res, vn=c( "res", "random", "inla.group(t, method = \"quantile\", n = 11)" ), 
              type="b", col="slategray", pch=19, lty=1, lwd=2.5 ,
              xlab="Bottom temperature (degrees Celsius)", ylab=ylab, cex=1.25, cex.axis=1.25, cex.lab=1.25 )
          dev.off()
    
    
          (fn = file.path( outputdir, "pca1.png"))
          png( filename=fn, width=1024, height=1024, pointsize=12, res=196 )
            carstm_plotxy( res, vn=c( "res", "random", "inla.group(pca1, method = \"quantile\", n = 11)" ), 
              type="b", col="slategray", pch=19, lty=1, lwd=2.5 ,
              xlab="PCA1", ylab=ylab, cex=1.25, cex.axis=1.25, cex.lab=1.25 )
          dev.off()
    
          (fn = file.path( outputdir, "pca2.png"))
          png( filename=fn, width=1024, height=1024, pointsize=12, res=196 )
            carstm_plotxy( res, vn=c( "res", "random", "inla.group(pca2, method = \"quantile\", n = 11)" ), 
              type="b", col="slategray", pch=19, lty=1, lwd=2.5 ,
              xlab="PCA2", ylab=ylab, cex=1.25, cex.axis=1.25, cex.lab=1.25 )
          dev.off()
    
          (fn = file.path( outputdir, "depth.png"))
          png( filename=fn, width=1024, height=1024, pointsize=12, res=196 )
            carstm_plotxy( res, vn=c( "res", "random", "inla.group(z, method = \"quantile\", n = 11)" ), 
              type="b", col="slategray", pch=19, lty=1, lwd=2.5  ,
              xlab="Depth (m)", ylab=ylab, cex=1.25, cex.axis=1.25, cex.lab=1.25 )
          dev.off()
          
          if (0) {
            # (fn = file.path( outputdir, "substrate.png"))
            # png( filename=fn, width=1024, height=1024, pointsize=12, res=196 )
            #   carstm_plotxy( res, vn=c( "res", "random", "inla.group(substrate.grainsize, method = \"quantile\", n = 11)" ), 
            #     type="b", col="slategray", pch=19, lty=1, lwd=2.5  ,
            #     xlab="Substrate grain size (mm)", ylab=ylab, cex=1.25, cex.axis=1.25, cex.lab=1.25 )
            # dev.off()
      
            fit = carstm_model( p=pW, DS="carstm_modelled_fit",  sppoly = sppoly ) # to load currently saved results
          
            plot( fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
            plot( fit$marginals.hyperpar$"Phi for space_time", type="l")  # posterior distribution of phi nonspatial dominates
            plot( fit$marginals.hyperpar$"Precision for space_time", type="l")
            plot( fit$marginals.hyperpar$"Precision for setno", type="l")
          }
        }
    
    
    
      }  # end spatiotemporal model
    
}


# ----------------------
# Part 3: assimilation of models


  assimilate_numbers_and_size = TRUE

  if (assimilate_numbers_and_size ) {

    # wgts_max = 1.1 # kg, hard upper limit
    # N_max = NULL
    # #  quantile( M$totno[ipositive]/M$data_offset[ipositive], probs=0.95, na.rm=TRUE )  
    
    # posterior sims 
  
    sims = carstm_posterior_simulations( pN=pN, pW=pW, pH=pH, sppoly=sppoly, pa_threshold=0.05, qmax=0.99 )
    sims = sims  / 10^6 # 10^6 kg -> kt;; kt/km^2

    
    SM = aggregate_biomass_from_simulations( 
      sims=sims, 
      sppoly=sppoly, 
      fn=carstm_filenames( pN, returnvalue="filename", fn="aggregated_timeseries" ), 
      yrs=pN$yrs, 
      method="sum", 
      redo=TRUE 
    ) 
    
    RES= SM$RES
    # RES = aggregate_biomass_from_simulations( fn=carstm_filenames( pN, returnvalue="filename", fn="aggregated_timeseries" ) )$RES

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


    # map it ..mean density
    sppoly = areal_units( p=pN )  # to reload

    vn = paste("biomass", "predicted", sep=".")

    outputdir = file.path( carstm_filenames( pN, returnvalue="output_directory"), "predicted_biomass_densitites" )

    if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

    B = apply( sims, c(1,2), mean ) 
    B[ which(!is.finite(B)) ] = NA

    brks = pretty( log10( quantile( B[], probs=c(0.05, 0.95), na.rm=TRUE )* 10^6)  )
  
    additional_features = snowcrab_features_tmap(pN)  # for mapping below

    for (i in 1:length(pN$yrs) ) {
      y = as.character( pN$yrs[i] )
      sppoly[,vn] = log10( B[,y]* 10^6 )
      outfilename = file.path( outputdir , paste( "biomass", y, "png", sep=".") )
      tmout =  carstm_map(  sppoly=sppoly, vn=vn,
          breaks=brks,
          additional_features=additional_features,
          title=paste( "log_10( Predicted biomass density; kg/km^2 )", y ),
          palette="-RdYlBu",
          plot_elements=c( "compass", "scale_bar", "legend" ), 
          outfilename=outfilename
      )
      tmout
      
    }
  
  }  # end assimilate size and numbers


 
  


# end
