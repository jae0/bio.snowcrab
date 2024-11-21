


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
# this variation 3 is about removing spatial and spatiotemporal effects


# -------------------------------------------------
# Part 1 -- construct basic parameter list defining the main characteristics of the study

  source( file.path( code_root, "bio_startup.R" )  )
   
  require(bio.snowcrab)   # loadfunctions("bio.snowcrab") 

  year.assessment = 2022
  
  yrs = 1999:year.assessment
  spec_bio = bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 )
  snowcrab_filter_class = "fb"     # fishable biomass (including soft-shelled )  "m.mat" "f.mat" "imm"

  
  # key name 
  carstm_model_label= paste( "1999_2022_variation3", snowcrab_filter_class, sep="_" )

  # params for number
  pN = snowcrab_parameters(
    project_class="carstm",
    yrs=yrs,   
    areal_units_type="tesselation",
    family = "nbinomial",
    carstm_modelengine = "inla.reduced3",
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
    carstm_modelengine = "inla.reduced3",
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
    carstm_modelengine = "inla.reduced3",
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
   

# ------------------------------------------------
# Part 2 -- spatiotemporal statistical model

  if ( spatiotemporal_model ) {

    # total numbers
    sppoly = areal_units( p=pN )
    M = snowcrab.db( p=pN, DS="carstm_inputs", sppoly=sppoly  )  # will redo if not found
    # can just copy datafile: carstm_inputs_snowcrab~tesselation~1~snowcrab~24~1~none~snowcrab_managementareas.rdata
    # into working directory to speed things up

  
    io = which(M$tag=="observations")
    ip = which(M$tag=="predictions")
    iq = unique( c( which( M$totno > 0), ip ) )
    iw = unique( c( which( M$totno > 5), ip ) )  # need a good sample to estimate mean size
    
    ### continue as in 03_biomass_index_carstm.R


  }  # end spatiotemporal model




# ----------------------
# Part 3: assimilation of models


  assimilate_numbers_and_size = TRUE

  if (assimilate_numbers_and_size ) {

    ### continue as in 03_biomass_index_carstm.R  .. alter save locations

  }  # end assimilate size and numbers


  

## Plot maps of residuals of numbers per set obs vs pred


#To add a title to any carstm_map, please see below example
#carstm_map( res=res, vn=vn, main=list(label="my plot title", cex=2) )


# -------------------------------------------------
# Part 1 -- construct basic parameter list defining the main characteristics of the study
# require(aegis)

  p = pN  #  bio.snowcrab::snowcrab_parameters( project_class="carstm", yrs=1999:year.assessment )

  outputdir = file.path( p$modeldir, p$carstm_model_label, "residuals" )

  # extract results and examine

  fit =  carstm_model( p=p, DS="modelled_fit" )  # extract currently saved model fit
  summary(fit)

  res = carstm_model( p=p, DS="carstm_modelled_summary"  )

  # prediction surface
  crs_lonlat = st_crs(projection_proj4string("lonlat_wgs84"))

  sppoly = areal_units( p=p )  # will redo if not found
  sppoly = st_transform(sppoly, crs=crs_lonlat )

  # do this immediately to reduce storage for sppoly (before adding other variables)
  M = snowcrab.db( p=p, DS="biological_data" )  # will redo if not found .. not used here but used for data matching/lookup in other aegis projects that use bathymetry
  M$tiyr=lubridate::decimal_date(M$timestamp)

  # M$totno = M$totno_adjusted / M$cf_set_no   # convert density to counts
  # M$totwgt = M$totwgt_adjusted / M$cf_set_mass # convert density to total wgt
  # M$data_offset = 1 / M$cf_set_no  ## offset only used in poisson model

  # reduce size
  M = M[ which( M$lon > p$corners$lon[1] & M$lon < p$corners$lon[2]  & M$lat > p$corners$lat[1] & M$lat < p$corners$lat[2] ), ]
  # levelplot(z.mean~plon+plat, data=M, aspect="iso")

  M$AUID = st_points_in_polygons(
    pts = st_as_sf( M, coords=c("lon","lat"), crs=crs_lonlat ),
    polys = sppoly[, "AUID"],
    varname="AUID"
  )

  M = M[!is.na(M$AUID),]

  names(M)[which(names(M)=="yr") ] = "year"
  # M = M[ which(M$year %in% p$yrs), ]
  # M$tiyr = lubridate::decimal_date ( M$timestamp )
  # M$dyear = M$tiyr - M$year

  MM = res$M

  obsMM = MM[MM$tag=="observations",]
  plocs = MM[MM$tag=="predictions",]

  obs = M
  if (p$variabletomodel=="totno") {
    obs$density = obs$totno / obs$data_offset
    rr = as.data.frame.table(res$totno.predicted)
  }
  if (p$variabletomodel=="totwgt") {
    obs$density = obs$totwgt / obs$data_offset
    rr = as.data.frame.table(res$totwgt.predicted)
  }

  rr$AUID = as.character( rr$AUID)

  region.id = slot(sppoly, "region.id" )

  rr$space = match( rr$AUID, region.id )
  rr$AUID = rr$space

  rr$year = as.numeric( as.character( rr$year) )

  obs = merge( obs, rr, by=c("year", "AUID"), all.x=TRUE, all.y=FALSE )
  obs$Freq[ !is.finite(obs$Freq) ] = 0
  obs$resid =  obs$Freq - obs$density
  obs$resid_per_set = obs$resid * obs$data_offset
  obs$yr = obs$year


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


# end
