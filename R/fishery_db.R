
fishery_db = function(p=NULL, DS=NULL, sppoly=NULL) {
  
  # this function cleans up catch and effort data for modelling, by dropping problematic records 
  # as such, the total catch and effort will not be consistent with fishery quotamanagement
   
  # maturity codes
  immature = 0
  mature = 1
  mat.unknown = 2

  if ( DS=="carstm_inputs") {

    # prediction surface
    crs_lonlat = st_crs(projection_proj4string("lonlat_wgs84"))
    if (is.null(sppoly)) sppoly = areal_units( p=p )  # will redo if not found
    sppoly = st_transform(sppoly, crs=crs_lonlat )
    # sppoly$data_offset = sppoly$sa
    
    areal_units_fn = attributes(sppoly)[["areal_units_fn"]]

    if (exists("carstm_directory", p)) {
      outputdir =  p$carstm_directory 
    } else {
      outputdir = file.path( p$modeldir, p$carstm_model_label )
    }
    outfn = paste( sep="_") # redundancy in case other files in same directory
    
    fn = file.path( outputdir, "carstm_inputs.rdz" )
 
    if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

    if (!redo)  {
      if (file.exists(fn)) {
        message( "Loading previously saved carstm_inputs ... ", fn)
        M = aegis::read_write_fast( fn )
        return( M )
      }
    }
    message( "Generating carstm_inputs ... ", fn)

    M = logbook.db( DS="logbook" )

    h = which( !is.na( M$lon + M$lat ) )
    M = M[h,]

    i = polygon_inside( M[, c("lon", "lat")], region="cfaall")
    M = M[i,]

    j = polygon_inside( M[, c("lon", "lat")], region="isobath1000m")
    M = M[j,]

    # additional constraint ..
    # remove data that are strange in location .. land
    crs_lonlat = st_crs( projection_proj4string("lonlat_wgs84") )

    mp = st_multipoint(cbind(p$corners$lon, p$corners$lat))
    bbox =  st_as_sfc(st_bbox( mp) )
    st_crs(bbox) =  crs_lonlat 

    coast = st_transform( st_as_sf(coastline_db( p=p )), crs=crs_lonlat )
    coast = (
        st_intersection( coast, bbox )
        %>% st_buffer(0.01)
        %>% st_union()
        %>% st_cast("POLYGON" )
        %>% st_union()
        %>% st_make_valid()
    )

    inside = st_points_in_polygons(
        pts =st_as_sf( M[,c("lon", "lat")], coords=c("lon","lat"), crs=crs_lonlat ),
        polys = coast,  
        method="sp::point.in.polygon"
    )
    onland = which (inside)
    if (length(onland)>0) M = M[-onland, ]

    # filter by depth ..
    # use the match/map syntax in bathymetry and filter out shallow sets .. < 10 m? TODO
    # keep those in the domain and deeper than depth=10 m

    pL = aegis.bathymetry::bathymetry_parameters( project_class="stmv"  )
    LUT= aegis_survey_lookuptable( aegis_project="bathymetry", 
    project_class="core", DS="aggregated_data", pL=pL )

    M$z = M$depth
    i = which( !is.finite( M$z ))
    if (length(i) > 0 ) {
      M$z[i] = aegis_lookup( pL=pL, 
          LOCS=M[i, c("lon", "lat")], LUT=LUT,
          output_format="points" , 
          space_resolution=p$pres*2, variable_name="z.mean"  ) # core=="rawdata"
    }

    aoi = which( M$z > 10 ) # negative = above land
    M = M[ aoi,]
 
    # only accept "correctly" positioned data within each subarea ... in case miscoded data have a large effect
    icfa4x = polygon_inside( M[, c("lon", "lat")], "cfa4x")
    icfanorth = polygon_inside( M[, c("lon", "lat")], "cfanorth")
    icfa23 = polygon_inside( M[, c("lon", "lat")], "cfa23")
    icfa24 = polygon_inside( M[, c("lon", "lat")], "cfa24")

    gooddata = sort( unique( c(icfa4x, icfanorth, icfa23, icfa24 ) ) )
    M = M[gooddata, ]
    
    M$timestamp = M$date.fished
    M$tiyr = lubridate::decimal_date(M$timestamp)
    M$dyear = M$tiyr - lubridate::year(M$timestamp) 
         # reduce size
    M = M[ which( M$lon > p$corners$lon[1] & M$lon < p$corners$lon[2]  & M$lat > p$corners$lat[1] & M$lat < p$corners$lat[2] ), ]
    # levelplot(z.mean~plon+plat, data=M, aspect="iso")
 
    M = M[ which( is.finite(M$dyear) ), ]

    # enforce bounds in effort and cpue
    oo = which( M$cpue > 650 * 0.454 )  # 600 - 650 lbs / trap is a real/reasonable upper limit
    if ( length(oo) > 0 ) M$cpue[oo] = NA

    pp = which ( M$effort > 240 ) # small traps were used at times with large 240 trap compliments
    if ( length(pp) > 0 ) M$effort[pp] = NA
 
    qn = quantile( M$cpue, p$quantile_bounds[2], na.rm=TRUE )
    ni = which( M$cpue > qn )

    # copy landings to "catch" .. which we modify as below
    M$catch = M$landings
    M$catch[ni] = floor( qn * M$effort[ni] )  # limit upper bound to one that is physically possible
    M = M[is.finite(M$catch),]
    
    M$data_offset = M$effort
    M = M[is.finite(M$data_offset),]

    lower_threshold =  1 # 1 kg per trap seems a reasonable cutoff ~ 1 large crab
    M$pa = 0
    M$pa[ which( M$catch > lower_threshold) ] = 1
   
    # data_offset is per unit trap haul 
    M = carstm_prepare_inputdata( 
      p=p, M=M, sppoly=sppoly,
      APS_data_offset=1, 
      retain_positions_outside_of_boundary = 25,  # centroid-point unit of p$aegis_proj4string_planar_km
      vars_to_retain=c(
        "licence", "cfv", "cfa", "cfa0", "z", "catch", "data_offset", "pa",  # "soak.time", "trap.type", are motly missing data
      ) 
    )
  
    if ( exists("substrate.grainsize", M)) M$log.substrate.grainsize = log( M$substrate.grainsize )

    if (!exists("yr", M)) M$yr = M$year  
    
    # IMPERATIVE: 
    i =  which(!is.finite(M$z))
    j =  which(!is.finite(M$t)) 

    if (length(j)>0 | length(i)>0) {
      warning( "Some areal units that have no information on key covariates ... you will need to drop these and do a sppoly/nb reduction with areal_units_neighbourhood_reset() :")
          print( "Missing depths:")
      print(unique(M$AUID[i]) )
      print( "Missing temperatures:")
      print(unique(M$AUID[j] ) )
    }
 
    # imperative covariate(s)
    M = M[ which(is.finite(M$z)), ]  
    M = M[ which(is.finite(M$t)), ]  
 
    M$space = match( M$AUID, sppoly$AUID) # for bym/car .. must be numeric index matching neighbourhood graphs
    M$space_time = M$space  # copy for space_time component (INLA does not like to re-use the same variable in a model formula) 
    M$space_cyclic = M$space  # copy for space_time component (INLA does not like to re-use the same variable in a model formula) 

    M$time = match( M$year, p$yrs ) # copy for space_time component .. for groups, must be numeric index
    M$time_space = M$time    
 
    M$cyclic = match( M$dyri, discretize_data( span=c( 0, 1, p$nw)  ) )  # as integer
    M$cyclic_space = M$cyclic # copy cyclic for space - cyclic component .. for groups, must be numeric index
  
    read_write_fast( data=M, fn=fn )

    return( M )
  }

}