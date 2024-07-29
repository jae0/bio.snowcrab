
### OPTIONAL -- exploratory and not essential to assessment

# -------------------------------------------------
# Snow crab --- Areal unit modelling Hurdle / Delta model  
# combination of three models via posterior simulation
# 1. Poisson on positive valued numbers offset by swept area
# 2. Meansize in space and time 
# 3  Presence-absence
# the convolution of all three after simulation is called a Hurdle or Delta model
# -------------------------------------------------
  

# -------------------------------------------------
# Part 1 -- construct basic parameter list defining the main characteristics of the study

source( file.path( code_root, "bio_startup.R" )  )
 
require(bio.snowcrab)   # loadfunctions("bio.snowcrab") 
require(ggplot2)
require(Matrix)
require(spam)

# save results to a location outside of bio.data as this is not operational (yet) 
carstm_results_directory = file.path( homedir, "projects", "dynamical_model", "snowcrab", "data" )

year.assessment = 2023
yrs = 1999:year.assessment

p = bio.snowcrab::load.environment( year.assessment=year.assessment )
  
spec_bio = bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 )

temperature_figures_redo = FALSE
areal_units_redo = FALSE
    
assimilate_numbers_and_size = TRUE

additional_features = snowcrab_mapping_features(p, redo=FALSE )  # for mapping below

sppoly_tweaks = list(
  # vary params by variable as data densities vary for these size/age/sex groups 
  areal_units_constraint_ntarget= list( M0=8, M1=10, M2=14, M3=14, M4=14, f.mat=8 ),
  n_iter_drop=list( M0=0, M1=1, M2=1, M3=1, M4=1, f.mat=0  )
)


theta_init = list(
  notes = "These are solutions from 2023",
  M0=list(
    N = c(1.665, 2.838, 2.051, 0.526, 3.657, 0.486, 5.049, 5.303, 5.405, 6.716, 0.426, 3.264, 0.965, 1.921, 1.886 ),
    W = c(5.891, 8.441, 0.859, 2.837, 9.942, 7.441, 11.249, 11.576, 12.614, 11.109, 6.548, 3.713, 5.805, 3.408, 1.509),
    H = c(1.030, 1.657, 2.818, 1.320, -3.475, 3.014, 3.470, -1.314,  -1.903, -0.495, -1.821, 2.891)
  ),
  M1=list(     
    N = c(2.792, 3.282, 1.661, 4.182, 0.001, 3.627, 0.163, 5.787, 4.922, 6.469, 1.619, 2.147, 1.550, 0.761, 1.469 ),
    W = c(6.917, 8.171, 2.694, 10.034, 0.001, 9.370, 6.347, 11.923, 12.879, 11.565, 9.176, -0.398, 7.224, 2.632, 2.508),
    H = c(0.567, 1.432, 2.229, -0.001, 2.236, 2.329, 3.437, 3.780, 3.664, -1.571, 4.327, -0.768, 3.812, 2.014)
  ), 
  M2=list(    
    N = c(1.077, 2.255, 1.459, 2.106, 2.320, 2.582, -2.764, 4.167, 4.841, 5.640, 0.626, 2.340, 0.923, 1.219, 1.549),
    W = c(7.611, 9.986, 1.079, 9.742, 0.038, 9.362, 7.743, 9.320, 12.403, 11.791, 26.247, -2.238, 26.324, -0.897, -0.001 ),
    H = c(0.516, 1.407, 3.556, -0.005, 0.734, 2.281, 2.091, 2.842, 4.989, -1.260, 2.401, -1.098, 3.209, 2.494 )
  ),   
  M3=list(    
    N = c(1.213,  2.078, 1.044,   3.881,  3.972,  3.600,  4.320,  4.287,  1.130, -3.137,  0.288, -2.936, 1.006),
    W = c(9.209, 10.730, 1.203, 11.130, 0.194, 9.036, 8.687, 11.692, 14.667, 15.657, 10.529, 0.263, 9.239, 1.766, 0.932 ), 
    H = c(1.121,  1.153, 1.278,   -0.082, -1.959, 4.228,   4.422, -0.916, -1.424, 0.336, -1.371, 1.940)
  ),
  M4=list(    
    N = c(1.081, 2.657, 1.034, 3.567, 3.515, 3.742, 3.888, 4.504, 1.315, -2.853, 0.253, -3.166, 0.857),
    W = c(10.515, 12.125, 1.303, 13.383, 13.479, 11.835, 16.365, 16.582, 12.323, -2.478, 11.765, -2.205, 1.229),
    H = c(1.195, 0.740, 1.266, 1.021, 1.516, -3.517, 0.302, 4.378, 5.130, -1.010, 2.072, -0.574, 2.290, 1.935)
  ), 
  f.mat=list(
    N = c(0.556, 1.846, 1.482, 4.167, 1.589, 4.093, 3.084, 4.054, 4.285, 3.997, 0.775, -2.991, -0.007, -0.378, 1.859), 
    W = c(9.498, 12.436, 0.148, 11.175, 12.500, 10.804, 15.664, 15.261,  8.366, -1.689,  9.333, -1.629, 1.941),
    H = c(0.839, 1.329, -0.324, 1.926, -4.342, 4.939, 4.525, -0.768, -3.873, -0.515, -2.004, 2.498)
  )
)

    

for (snowcrab_filter_class in c(  "M0", "M1", "M2", "M3", "M4", "f.mat" ) ) {

  # snowcrab_filter_class = "M0"     # "fishable biomass" (excluding soft-shelled )
  # snowcrab_filter_class = "M1"     # some fraction expected to enter M0 next year
  # snowcrab_filter_class = "M2"     # some fraction expected to enter M1 next year / M0 in 2 years
  # snowcrab_filter_class = "M3"     # some fraction expected to enter M2 next year / M0 in 3 years
  # snowcrab_filter_class = "M4"     # some fraction expected to enter M3 next year / M0 in 4 years
  # snowcrab_filter_class = "f.mat"
  
  # snowcrab_filter_class = "imm"  # note poisson will not work due to var inflation .. nbinomial is a better choice 
  # snowcrab_filter_class = "m.mat"

  
  carstm_model_label= paste( "default", snowcrab_filter_class, sep="_" )

  # poisson works too but variance is not exactly poisson (higher than mean)
  Nfamily = switch( snowcrab_filter_class, 
    M0 = "nbinomial",   
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
    carstm_model_label= carstm_model_label,  
    carstm_directory = file.path(carstm_results_directory, carstm_model_label ),
    theta = theta_init[[snowcrab_filter_class]][["N"]],
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
    carstm_directory = file.path( carstm_results_directory, carstm_model_label  ),
    theta = theta_init[[snowcrab_filter_class]][["W"]],
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
    carstm_directory = file.path(carstm_results_directory, carstm_model_label ),
    theta = theta_init[[snowcrab_filter_class]][["H"]],
    selection = list(
      type = "presence_absence",
      biologicals=list( spec_bio=spec_bio ),
      biologicals_using_snowcrab_filter_class=snowcrab_filter_class
    )
  )
    
  pN$areal_units_constraint_ntarget = sppoly_tweaks[["areal_units_constraint_ntarget"]][[snowcrab_filter_class]]
  pN$n_iter_drop = sppoly_tweaks[["n_iter_drop"]][[snowcrab_filter_class]]


  if (areal_units_redo) {
    # polygon structure:: create if not yet made
    # for (au in c("cfanorth", "cfasouth", "cfa4x", "cfaall" )) plot(polygon_managementareas( species="snowcrab", au))
            
    xydata = snowcrab.db( p=pN, DS="areal_units_input", redo=TRUE )
    xydata = snowcrab.db( p=pN, DS="areal_units_input" )

    sppoly = areal_units( p=pN, xydata=xydata[ which(xydata$yr %in% pN$yrs), ], redo=TRUE, verbose=TRUE )  # create constrained polygons with neighbourhood as an attribute

    plot(sppoly["npts"])

    sppoly$dummyvar = ""
    xydata = st_as_sf( xydata, coords=c("lon","lat") )
    st_crs(xydata) = st_crs( projection_proj4string("lonlat_wgs84") )
  
    tmap_mode("plot")

    require("tmap")
    
    plt = 
      tm_shape(sppoly) +
        tm_borders(col = "slategray", alpha = 0.5, lwd = 0.5) + 
        tm_shape( xydata ) + tm_sf() +
        additional_features[["tmap"]] +
        tm_compass(position = c("right", "TOP"), size = 1.5) +
        tm_scale_bar(position = c("RIGHT", "BOTTOM"), width =0.1, text.size = 0.5) +
        tm_layout(frame = FALSE, scale = 2) +
        tm_shape( st_transform(polygons_rnaturalearth(), st_crs(sppoly) )) + 
        tm_borders(col = "slategray", alpha = 0.5, lwd = 0.5)

    dev.new(width=14, height=8, pointsize=20)
    plt

  }


  sppoly=areal_units( p=pN )

  pN$space_name = sppoly$AUID 
  pN$space_id = 1:nrow(sppoly)  # must match M$space

  pN$time_name = as.character(pN$yrs)
  pN$time_id =  1:pN$ny

  pN$cyclic_name = as.character(pN$cyclic_levels)
  pN$cyclic_id = 1:pN$nw

  pW$space_name = sppoly$AUID 
  pW$space_id = 1:nrow(sppoly)  # must match M$space

  pW$time_name = as.character(pW$yrs)
  pW$time_id =  1:pW$ny

  pW$cyclic_name = as.character(pW$cyclic_levels)
  pW$cyclic_id = 1:pW$nw

  pH$space_name = sppoly$AUID 
  pH$space_id = 1:nrow(sppoly)  # must match M$space

  pH$time_name = as.character(pH$yrs)
  pH$time_id =  1:pH$ny

  pH$cyclic_name = as.character(pH$cyclic_levels)
  pH$cyclic_id = 1:pH$nw

  if (temperature_figures_redo) {
  
    # area-specific figures
    figure_area_based_extraction_from_carstm(DS="temperature", year.assessment )  # can only do done once we have an sppoly for snow crab
  
    # full domain:
    # default paramerters (copied from 03_temperature_carstm.R )  
    require(aegis.temperature)
    params = list( 
      temperature = temperature_parameters( 
        project_class="carstm", 
        yrs=1970:year.assessment, 
        carstm_model_label="default"
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


    
  M = snowcrab.db( p=pN, DS="carstm_inputs", sppoly=sppoly, redo=TRUE )  # will redo if not found

  
  
  

  # ------------------------------------------------
  # Part 2 -- spatiotemporal statistical model
  spatiotemporal_model = TRUE
  if ( spatiotemporal_model ) {

    # total numbers
    
    sppoly=areal_units( p=pN )
    
    M = snowcrab.db( p=pN, DS="carstm_inputs", sppoly=sppoly  )  # will redo if not found
  
    io = which(M$tag=="observations")
    ip = which(M$tag=="predictions")
    iq = unique( c( which( M$totno > 0), ip ) )
    iw = unique( c( which( M$totno > 3), ip ) )  # need a good sample to estimate mean size

    # number 
    res = NULL; gc()
    res = carstm_model( p=pN, data=M[ iq, ], sppoly=sppoly, 
      nposteriors=5000,
      posterior_simulations_to_retain=c( "summary", "random_spatial", "predictions"), 
      control.inla = list( int.strategy="eb", strategy="adaptive", h=0.01 ),  # simplified.laplace
      # redo_fit=FALSE, 
      verbose=TRUE,
      num.threads="4:3"  
    )

    # mean size
    res = NULL; gc()
    res = carstm_model( p=pW, data=M[ iw, ], sppoly = sppoly, 
      nposteriors=5000,
      posterior_simulations_to_retain=c( "summary", "random_spatial", "predictions"), 
      control.inla = list( int.strategy="eb", strategy="adaptive", h=0.01 ),  # simplified.laplace
      # redo_fit=FALSE, 
      verbose=TRUE,
      num.threads="4:3"   
    ) 

    # model pa using all data
    res = NULL; gc()
    res = carstm_model( p=pH, data=M, sppoly=sppoly, 
      nposteriors=5000,
      posterior_simulations_to_retain=c( "summary", "random_spatial", "predictions"), 
      # control.family=list(control.link=list(model="logit")),  # default
      control.inla = list( int.strategy="eb", strategy="adaptive", h=0.01 ),  # simplified.laplace
      # redo_fit=FALSE, 
      verbose=TRUE,
      num.threads="4:3"  
    )
    res = NULL; gc()
 
  }
  
  # end spatiotemporal model

  # some maps and plots

    for (vns in c( "number", "meansize", "habitat") ) {
 
      if ( vns=="number" ) {
        p=pN
        ylab = "Number"
        outputdir = file.path( p$modeldir, p$carstm_model_label, "predicted.numerical.densities" )
        if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )
        fn_root_prefix = "Predicted_numerical_abundance"
        fn_root =  "Predicted_numerical_abundance_persistent_spatial_effect" 
        outfilename = file.path( outputdir, paste(fn_root, "png", sep=".") )
        title= paste( snowcrab_filter_class, "Number; no./m^2"  )
      }
      if ( vns=="meansize") {
        p=pW
        ylab = "Mean weight"
        outputdir = file.path( p$modeldir, p$carstm_model_label, "predicted.meansize" )
        if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )
        fn_root_prefix = "Predicted_meansize"
        fn_root =  "Predicted_meansize_persistent_spatial_effect" 
        outfilename = file.path( outputdir, paste(fn_root, "png", sep=".") )
        title= paste( snowcrab_filter_class, "Mean weight; kg" ) 
      }
      if ( vns=="habitat" ) {
        p=pH
        ylab = "Probability"
        outputdir = file.path( p$modeldir, p$carstm_model_label, "predicted.presence_absence" )
        if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )
        fn_root_prefix = "Predicted_presence_absence"
        fn_root =  "Predicted_presence_absence_persistent_spatial_effect" 
        outfilename = file.path( outputdir, paste(fn_root, "png", sep=".") )
        title= paste( snowcrab_filter_class, "Probability")  

          # if (vns=="habitat") {
          #   # to compute habitat prob
          #   sims = carstm_posterior_simulations( pH=pH, pa_threshold=0.05, qmax=0.95 )
          #   SM = aggregate_simulations( 
          #     sims=sims, 
          #     sppoly=sppoly, 
          #     fn=carstm_filenames( pN, returnvalue="filename", fn="aggregated_timeseries" ), 
          #     yrs=pN$yrs, 
          #     method="mean", 
          #     redo=TRUE 
          #   ) 
          #   outputdir = file.path( carstm_filenames( pN, returnvalue="output_directory"), "aggregated_habitat_timeseries" )
          #   ylabel = "Habitat probability"
          #   fn_ts = "habitat_M0.png"
          #   vn = paste("habitat", "predicted", sep=".")
          #   outputdir2 = file.path( carstm_filenames( pN, returnvalue="output_directory"), "predicted_habitat" )

          # }         
      
      }

      res = carstm_model( p=p, DS="carstm_modelled_summary",  sppoly = sppoly ) # to load currently saved results

      plots_from_fit = FALSE
      if (plots_from_fit) {
        # extract results
        fit = carstm_model( p=p, DS="modelled_fit",  sppoly = sppoly )  # extract currently saved model fit
        fit$summary$dic$dic
        fit$summary$dic$p.eff
        plot(fit)
        plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
        plot( fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
        plot( fit$marginals.hyperpar$"Phi for space_time", type="l")  # posterior distribution of phi nonspatial dominates
        plot( fit$marginals.hyperpar$"Precision for space_time", type="l")
        plot( fit$marginals.hyperpar$"Precision for setno", type="l")
        fit = NULL
      }

      
      if ( vns=="habitat" ) {
        habitat_2D = FALSE
        if (habitat_2D) {
        
          o = carstm_2D_effects_probability( 
            res,
            xvar = "inla.group(t, method = \"quantile\", n = 13)",  
            yvar = "inla.group(z, method = \"quantile\", n = 13)", 
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
            xvar = "inla.group(t, method = \"quantile\", n = 13)",  
            yvar = "inla.group(z, method = \"quantile\", n = 13)",
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
            u = aegis::read_write_fast('/home/jae/tmp/temp_depth_habitat.RDS')
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
          read_write_fast( data=o, file=fn_optimal )
          o = aegis::read_write_fast(fn_optimal)

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
      }
 
      # maps
      vn = c( "random", "space", "re_total" ) 
      toplot = carstm_results_unpack( res, vn )
      brks = pretty(  quantile(toplot[,"mean"], probs=c(0,0.975), na.rm=TRUE )  )

      carstm_map(  res=res, vn=vn, 
        sppoly = sppoly, 
        breaks = brks,
        palette="-RdYlBu",
        plot_elements="",
        additional_features=additional_features,
        outfilename=outfilename
      )  
        

      vn="predictions"
      toplot = carstm_results_unpack( res, vn )
      brks = pretty(  quantile(toplot[,,"mean"], probs=c(0,0.975), na.rm=TRUE )  )

      for (y in res$time_id ){
        tmatch = as.character(y)
        fn_root = paste(fn_root_prefix, paste0(tmatch, collapse="-"), sep="_")
        outfilename = file.path( outputdir, paste(fn_root, "png", sep=".") )

        carstm_map(  res=res, vn=vn, tmatch=tmatch,
          sppoly = sppoly, 
          breaks =brks,
          palette="-RdYlBu",
          plot_elements="",
          additional_features=additional_features,
          title = y, 
#          title=paste(fn_root_prefix, snowcrab_filter_class,  paste0(tmatch, collapse="-") )
          outfilename=outfilename
        )
        # print(outfilename)
      
      }

      # plots with 95% PI
      res = carstm_model( p=p,  DS="carstm_summary" )  # parameters in p and direct summary
      
      oeffdir = file.path( outputdir, fn_root_prefix, "effects" )
      if ( !file.exists(oeffdir)) dir.create( oeffdir, recursive=TRUE, showWarnings=FALSE )

      # annual component 
      res$time$yr = as.numeric(p$time_name[res$time$ID])
      plt = ggplot( res$time, aes(x=yr, y=mean)) +  geom_line(color="gray", linewidth=1.75) + geom_point( size=3, color="slategray") +
        geom_errorbar(aes(ymin=quant0.025, ymax=quant0.975), color="slategray",  linewidth=1.0, width=0.5) +
        labs(title="Annual component", x="Year", y =ylab) +
        theme_light( base_size=22) 
      print(plt)    
      (fn_plt = file.path( oeffdir, "time.png"))
      ggsave(filename=fn_plt, plot=plt, dpi=144, width=12, height = 8)

      # seasonal component 
      res$cyclic$seas = as.numeric( p$cyclic_name[res$cyclic$ID] )
      plt = ggplot( res$cyclic, aes(x=seas, y=mean) ) + geom_line(color="gray", linewidth=1.75) + geom_point( size=3, color="slategray") +
        geom_errorbar(aes(ymin=quant0.025, ymax=quant0.975), color="slategray",  linewidth=1.0, width=0.05) +
        labs(title="Seasonal component", x="Year (fraction)", y =ylab) +
        theme_light( base_size=22) 
      print(plt)          
      (fn_plt = file.path( oeffdir, "cyclic.png"))
      ggsave(filename=fn_plt, plot=plt, dpi=144, width=12, height = 8)


      # relationship with depth
      vn = "inla.group(z, method = \"quantile\", n = 13)"
      plt = ggplot( res[[vn]], aes(x=ID, y=mean ))+ geom_line(color="gray", linewidth=1.75) + geom_point( size=3, color="slategray") +
        geom_errorbar(aes(ymin=quant0.025, ymax=quant0.975), color="slategray",  linewidth=1.0, width=5) +
        labs(title="Depth component", x="Depth (m)", y =ylab) +
        theme_light( base_size=22) 
      print(plt)    
      (fn_plt = file.path( oeffdir, "depth.png" )) 
      ggsave(filename=fn_plt, plot=plt, dpi=144, width=12, height = 8)

  
      vn = "inla.group(t, method = \"quantile\", n = 13)"
      plt = ggplot( res[[vn]], aes(x=ID, y=mean ))+ geom_line(color="gray", linewidth=1.75) + geom_point( size=3, color="slategray") +
        geom_errorbar(aes(ymin=quant0.025, ymax=quant0.975), color="slategray",  linewidth=1.0, width=5) +
        labs(title="Temperature component", x="Bottom temperature (deg C)", y =ylab) +
        theme_light( base_size=22) 
      print(plt)    
      (fn_plt = file.path( oeffdir, "temperature.png")) 
      ggsave(filename=fn_plt, plot=plt, dpi=144, width=12, height = 8)


      vn= "inla.group(pca1, method = \"quantile\", n = 13)" 
      plt = ggplot( res[[vn]], aes(x=ID, y=mean ))+ geom_line(color="gray", linewidth=1.75) + geom_point( size=3, color="slategray") +
        geom_errorbar(aes(ymin=quant0.025, ymax=quant0.975), color="slategray",  linewidth=1.0, width=5) +
        labs(title="PCA1 component", x="PCA1", y =ylab) +
        theme_light( base_size=22) 
      print(plt)    
      (fn_plt = file.path( oeffdir, "pca1.png")) 
      ggsave(filename=fn_plt, plot=plt, dpi=144, width=12, height = 8)


      vn = "inla.group(pca2, method = \"quantile\", n = 13)" 
      plt = ggplot( res[[vn]], aes(x=ID, y=mean ))+ geom_line(color="gray", linewidth=1.75) + geom_point( size=3, color="slategray") +
        geom_errorbar(aes(ymin=quant0.025, ymax=quant0.975), color="slategray",  linewidth=1.0, width=5) +
        labs(title="PCA2 component", x="PCA2", y =ylab) +
        theme_light( base_size=22) 
      print(plt)    
      (fn_plt = file.path( oeffdir, "pca2.png")) 
      ggsave(filename=fn_plt, plot=plt, dpi=144, width=12, height = 8)

      if (0) {
        vn = "inla.group(substrate.grainsize, method = \"quantile\", n = 13)"  
        plt = ggplot( res[[vn]], aes(x=ID, y=mean ))+ geom_line(color="gray", linewidth=1.75) + geom_point( size=3, color="slategray") +
          geom_errorbar(aes(ymin=quant0.025, ymax=quant0.975), color="slategray",  linewidth=1.0, width=5) +
          labs(title="Depth component", x="Depth (m)", y =ylab) +
          theme_light( base_size=22) 
        print(plt)    
        (fn_plt = file.path( oeffdir, "substrate.png")) 
        ggsave(filename=fn_plt, plot=plt, dpi=144, width=12, height = 8)
      }




 
    }





    # ----------------------
    # Part 3: assimilation of models


    if (assimilate_numbers_and_size ) {

      if (snowcrab_filter_class == "M0" ) {

        # wgts_max = 1.1 # kg, hard upper limit
        # N_max = NULL
        # #  quantile( M$totno[ipositive]/M$data_offset[ipositive], probs=0.95, na.rm=TRUE )  
        
        # posterior sims 
        
          sims = carstm_posterior_simulations( pN=pN, pW=pW, pH=pH, pa_threshold=0.05, qmax=0.99 )
          sims = sims  / 10^6 # 10^6 kg -> kt;; kt/km^2

          
          SM = aggregate_simulations( 
            sims=sims, 
            sppoly=sppoly, 
            fn=carstm_filenames( pN, returnvalue="filename", fn="aggregated_timeseries" ), 
            yrs=pN$yrs, 
            method="sum", 
            redo=TRUE 
          ) 
          
          RES= SM$RES
          # RES = aggregate_simulations( fn=carstm_filenames( pN, returnvalue="filename", fn="aggregated_timeseries" ) )$RES

          outputdir = file.path( carstm_results_directory, "aggregated_biomass_timeseries" )

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



          regions = c( "cfanorth", "cfasouth", "cfa4x" )
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
      #     legend.position=c( 0.1, 0.9 ),
            # labs(x="Year", y="Biomass index (kt)", size = rel(1.5)) +
            scale_colour_manual(values=color_map) +
            scale_fill_manual(values=color_map) +
            scale_shape_manual(values = c(15, 17, 19)) +
            theme_light( base_size = 22) + 
            theme( legend.position=c(0.75, 0.9), legend.title=element_blank()) +
            scale_y_break(c(14, 28), scales = 1)
            
            # scale_y_continuous( limits=c(0, 300) )  
            ggsave(filename=fn, plot=out, width=12, height = 8)
      


          pN$areal_units_constraint_ntarget = sppoly_tweaks[["areal_units_constraint_ntarget"]][[snowcrab_filter_class]]
          pN$n_iter_drop = sppoly_tweaks[["n_iter_drop"]][[snowcrab_filter_class]]
              
          # map it ..mean density
          sppoly=areal_units( p=pN )

          vn = paste("biomass", "predicted", sep=".")

          outputdir = file.path( carstm_filenames( pN, returnvalue="output_directory"), "predicted_biomass_densitites" )

          if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

          B = apply( sims, c(1,2), mean ) 
          B[ which(!is.finite(B)) ] = NA

          brks = pretty( log10( quantile( B[], probs=c(0.05, 0.95), na.rm=TRUE )* 10^6)  )
        
          for (i in 1:length(pN$yrs) ) {
            y = as.character( pN$yrs[i] )
            sppoly[,vn] = log10( B[,y]* 10^6 )
            outfilename = file.path( outputdir , paste( "biomass", y, "png", sep=".") )
            carstm_map(  
                sppoly=sppoly, 
                vn=vn,
                breaks=brks,
                additional_features=additional_features,
                legend.position=c( 0.1, 0.9 ),
                annotation=y,
                # annotation=paste( "log_10( Predicted biomass density; kg/km^2 )", y ),
                colors=rev(RColorBrewer::brewer.pal(5, "RdYlBu")),
                outfilename=outfilename
            ) 
            
          }
      }

    }  # end assimilate size and numbers

  }  # end for loop categories



# if prepping data for continuous version (julia): 
#   if (grepl("size_structured", model_variation)) {
      # fishery landings has a weekly time step = 2/52 ~ 0.0385 ~ 0.04  X dt=0.01 seems to work best
        # save directory is defined at the start:
        # carstm_results_directory = file.path( homedir, "projects", "dynamical_model", "snowcrab", "data" )

        fishery_model_data_inputs( year.assessment=year.assessment, type="size_structured_numerical_dynamics",
          for_julia=TRUE, sppoly_tweaks=sppoly_tweaks, save_location=carstm_results_directory)
#  }
 
 
  


# end
