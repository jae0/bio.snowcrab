code_to_taxa = function(tx){
	taxonomy.recode( from="spec", tolookup=tx ) 
}


get_feeding_data = function( data_dir, redo=FALSE ) {

  # Hi Jae,

  # Here is the requested data. There is some variability in the resolution of records, hence you’ll find that some predators may have full weight information while others do not, and the same will apply to prey species.  This will limit how informative this data will be and how far you can use these analytically.

  # The new database schema includes unique identifiers to ensure records are as intended (prior to this there was a lot of duplication of records).  These are as described:
  # SET_SEQ – relates to the set information details
  # PRED_SEQ – relates to the predator fish the stomach was taken from
  # PREY_SEQ – relates to the identified prey within a stomach

  # Some column headers may need further description:
  # DEPTH – average set depth in meters
  # GEAR – set fishing gear code (list attached)
  # SPEC – predator id code (species list attached)
  # FWT – predator weight in grams
  # FLEN – predator length in cm
  # STOWGT/EMPTYWGT – stomach weights before and after processing in grams
  # FULLNESS – level of fullness associated with stomachs, where codes are described from 0 through 6:
  # 0             empty - no food contents
  # 1             less than ¼ full
  # 2             ¼ to ½ full
  # 3             ½ to ¾ full
  # 4             ¾ full to full
  # 5             everted
  # 6             regurgitated
  # FGEN – Gender code of predator:
  # 0             Unknown
  # 1             Male
  # 2             Female
  # PREYSPECCD – prey id code (species list attached)
  # PWT – prey weight in grams
  # PLEN – prey length in cm
  # DIGESTION – level of digestion of prey, where:
  # 1             Good Condition
  # 2             Partly Digested
  # 3             Well Digested
  # 4             Unidentifiable
  # 9             Unidentified digestion state

  fn = file.path( data_dir, "feeding_network.rdz" ) 
  diet = NULL
  if (!redo) {
    if ( file.exists(fn)) {
      diet = aegis::read_write_fast(fn)
      return( diet )
    }
  }

  diet = fread( file.path( data_dir, "FH.DataRequest.Dec.2022.csv" ) )
  names(diet) = tolower(names(diet))
  diet = diet[ which(!is.na(preyspeccd )), ]

  species = fread( file.path( data_dir, "FH.Species.List.OCT2022.csv" ), fill=TRUE ) # comments contain multiple separators, header added to compensate 
  names(species) = tolower(names(species))
 
  gear = fread( file.path( data_dir, "FH.GEAR.Dec.2022.csv" ) )
  names(gear) = tolower(names(gear))
 

  # trim to living organisms   
  # partial matches
  todrop_labels = c("UNID ", "BAIT ")
  todrop = NULL
  for (i in todrop_labels) {
    j = grep( i, species$common )
    todrop = c( todrop, species$speccd[ j ] )
  }
  todrop =   unique(todrop)  
  species = species[ -which(speccd %in% todrop), ]

  # exact matches
  todrop_labels = c( "MUCUS", "BAIT", "MUD", "FLUID", "SAND", "WATER", "OIL (CRUDE)", "", "FISH REMAINS (NS)",
     "EGGS (NS)", "PURSE (SKATE NS)", "FISH EGGS", "EGGS (SCULPIN NS)", "FISH LARVAE (NS)", "EGGS (CRAB NS)", "EGGS (TUBE WORM NS)",
     "EGGS (SNAIL/SLUG NS)", "EGGS (SQUID NS)", "CRUSTACEA", "CRAB", "WORM CASTS", "BONE FISH",
     "WORM CAST", "CRUSTACEANS (NS)", "EGGS (FISH NS)"  
    )
  todrop = NULL
  for (i in todrop_labels) {
    j = which( i == species$common )
    todrop = c( todrop, species$speccd[ j ] )
  }
  todrop = unique(as.numeric(todrop)  )
  ii = which(species$speccd %in% todrop)
  if (length(ii) > 0) species = species[ -ii, ]
 
  todrop_manual = c( "DECAPODA", "CRAB", "CANCER" )
  todrop_cd = c(2241, 629, 630 ) 
  ii = which(species$speccd %in% todrop_cd)  # 93
  if (length(ii) > 0) species = species[ -ii, ]

  # strange or uninformative data:
  ii = which( species$speccd == 7208 ) # walrus eaten by small hake and ocean pout .. likely an error
  if (length(ii) > 0) diet = diet[-ii,]

  #recode to parsimony
  require(bio.taxonomy)

  diet$spec = bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=diet$spec  )
  diet$preyspeccd = bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=diet$preyspeccd  )

  species$speccd = bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=species$speccd  )
  ii = which(duplicated(species$speccd))
  if (length(ii) > 0) species = species[-ii, ]

  # reduce species list to what is observed in diets
  species$speccd = as.numeric(species$speccd)
  species = species[ is.finite(speccd)]
  species = species[ !is.na(phylum) ]

  sps_uniq = unique( c(unique( diet$preyspeccd), unique(diet$spec ) ) )
  species = species[ speccd %in% sps_uniq ]
  
  ii = which( diet$preyspeccd %in% species$speccd) 
  if (length(ii) > 0) diet = diet[ ii, ]
  
  # strange or uninformative data:
  # ii = which( diet$preyspeccd == 7208 ) # walrus eaten by small hake and ocean pout .. likely an error
  # if (length(ii) > 0) diet = diet[-ii,]

  tocollapse = c(
    "TAPEWORMS", "SEA URCHIN", "BRITTLE STAR", "SEA URCHIN",
    "SEA SPIDER", "OCTOPUS", "NUTCLAMS", "SEA MOUSE", "FLABELLIGERID WORMS", "POLYCHAETE",
    "ISOPOD", "SKELETON SHRIMP"
  )
  for (i in tocollapse) {
    k = grep(i, species$common, fixed=TRUE )
    if (length(k)> 1) {
      id0 = k[1]
      for (j in k) { 
        m = grep( species$speccd[j], diet$spec )
        if ( length(m)>0 ) {
          diet$spec[m] = species$speccd[id0] 
        }
        species[j,] = species[id0,] 
      }
    }
  }

  species = species[ which(!duplicated(species)) , ]
   

  diet$stime =  paste("0000", diet$stime, sep="") 
  nchars = nchar( diet$stime )
  diet$stime = paste( substr( diet$stime, (nchars-3), nchars ), "00", sep="" )
  diet$timestamp = lubridate::dmy_hms( paste( diet$sdate, diet$stime), tz="UTC" )

  diet = diet[ is.finite(slatdd + slongdd),]
  diet = diet[, .(datasource, mission, setno, timestamp, bottom_temperature, depth, slatdd, slongdd, nafo_zone, nafo_subunit, spec, fshno, fwt, flen, stowgt, emptywgt, fullness, 
   preyspeccd, pwt, plen, pnum, digestion )]
 
  read_write_fast( data=diet, fn=fn )
  return( diet )
}


survey_data = function( data_dir, redo=FALSE, cthreshold = 0.005 ) {
    # fixed time snapshot to work with stomach data
    
  fn = file.path( data_dir, "survey_data.rdz" ) 
  res = NULL
  if (!redo) {
    if ( file.exists(fn)) {
      res = aegis::read_write_fast(fn)
      return( res )
    }
  }

  
  year.assessment = 2022  # last year of stomach data above
  yrs = 1995:year.assessment
  
  carstm_model_label="1995_2022"

  require(aegis.speciescomposition)

  p = speciescomposition_parameters( yrs=yrs, carstm_model_label=carstm_model_label )

  res = survey_data_prepare(p=p, cthreshold = cthreshold) # bring in aegis_survey data first

  p0 = speciescomposition_parameters(
    project_class="carstm",
    data_root = project.datadirectory( "aegis", "speciescomposition" ),
    variabletomodel = "",  # will b eover-ridden .. this brings in all pca's and ca's
    carstm_model_label = carstm_model_label,
    carstm_model_label = carstm_model_label,
    inputdata_spatial_discretization_planar_km = 0.5,  # km controls resolution of data prior to modelling to reduce data set and speed up modelling
    inputdata_temporal_discretization_yr = 1/52,  # ie., every 1 weeks .. controls resolution of data prior to modelling to reduce data set and speed up modelling
    year.assessment = max(yrs),
    yrs = yrs, 
    spatial_domain = "SSE",  # defines spatial area, currenty: "snowcrab" or "SSE"
    areal_units_resolution_km = 1, # km dim of lattice ~ 1 hr
    areal_units_proj4string_planar_km = aegis::projection_proj4string("utm20"),  # coord system to use for areal estimation and gridding for carstm
    areal_units_type = "tesselation",       
    areal_units_overlay = "none" 
  )
  
  crs_lonlat = st_crs(projection_proj4string("lonlat_wgs84"))
  xydata = res$set[, .(lon, lat, timestamp, yr)]
  sppoly = areal_units( p=p0, xydata=xydata, spbuffer=5, n_iter_drop=3, verbose=TRUE )  # to force create

  # sppoly = areal_units( p=p0)
  sppoly = st_transform(sppoly, st_crs(crs_lonlat))
  
  diet = get_feeding_data( data_dir )
  Dpts = st_as_sf( diet[,c("slongdd","slatdd")], coords=c("slongdd","slatdd"), crs=crs_lonlat )

  # observations
  # length(unique(diet$id))  # 4558

  diet$AUID = st_points_in_polygons(
    pts = Dpts,
    polys = sppoly[, "AUID"],
    varname = "AUID"
  )

  diet$pnum[ !is.finite(diet$pnum)] = 1  # where missing assume minimal number
 
  diet$id = paste( diet$mission, diet$setno, sep=".")
  diet$id2 = paste(diet$id, diet$spec, sep="." )
  diet$id3 = paste(diet$id2, diet$fshno, sep="." )
  diet$id4 = paste(diet$id3, diet$preyspeccd, sep="." )
 
  diet_aggregated = diet[, .(Nprey=sum(pnum)), by=.(id4) ]
   
  diet_aggregated = diet_aggregated[ 
    diet[ which(!duplicated(id4)),.(spec, fshno, fwt, flen, preyspeccd, pwt, id, id2, id3, id4)], 
    on="id4" ]

  # match to clean groundfish database
  res$sc$spec = bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=res$sc$spec_bio  )
  res$sc$id2 = paste(res$sc$id, res$sc$spec, sep="." )
      
  diet_aggregated = res$sc[diet_aggregated, on=.(id2) ]
  # species ID mismatches (predator catch) .. likely data entry errors ... 125 aggregated records .. dropping for now
  ii = which( is.na(diet_aggregated$spec_bio) ) 
  if (length(ii)> 0) {
    possible.error = unique(diet_aggregated$id2[ii])
    diet_aggregated = diet_aggregated[-ii,]
    jj = which(diet_aggregated$id2 %in% possible.error) ## 205 entries
    if (length(jj)>0) {
      # diet_aggregated[jj, ]
      diet_aggregated = diet_aggregated[-jj, ]
    }
  }

  todrop = which( grepl("i\\.", names(diet_aggregated)) )
  diet_aggregated[,(todrop) := NULL]
  
  redundant = which( names(diet_aggregated) %in% c("datasource", "mission", "setno", "bottom_temperature", "depth", 
    "slatdd", "slongdd", "stowgt", "emptywgt", "fullness", "digestion", "plen", "id2", "id3", "id4" ) )
  diet_aggregated[,(redundant) := NULL]

  diet_aggregated$bottom_temperature = diet_aggregated$t
  diet_aggregated$depth = diet_aggregated$z
  
  res$sppoly = sppoly
  res$p = p  #basic survey data
  res$p0 = p0  # carstm config
  res$diet = diet_aggregated
  
  read_write_fast( data=res, fn=fn )

  return(res) 

}


diet_matrix = function(rdt, cthreshold=0) {
  # create (id = auid,time) X (sp= prey species) matrix
  
  m = data.table::dcast( setDT(rdt), 
    formula =  id + spec ~ preyspeccd, value.var="zscore", 
    fun.aggregate=mean, fill=NA, drop=FALSE, na.rm=TRUE
  )  # mean is just to keep dcast happy

  id = paste(m$id, m$spec, sep=".")
  m$id = NULL
  m$spec = NULL
  sps = colnames(m)

  m = as.matrix(m[]) 
  dimnames(m)[[1]] = id
 
  # remove low counts (absence) in the timeseries  .. species (cols) only
  # cthreshold = 0.005  # 0.01%  -> 97;  0.05 => 44 species; 0.001 => 228 species; 0.005 -> 146 
  m [ m < cthreshold ] = 0
  m [ !is.finite(m) ] = 0

  # reduce:: remove taxa/locations with very low counts (~cthreshold)
  finished = FALSE
  while( !(finished) ) {
    oo = colSums( ifelse( is.finite(m), 1, 0 ), na.rm=TRUE)
    i = unique( which(rowSums(m, na.rm=TRUE) == 0 ) )
    uu = colSums(m, na.rm=TRUE)/oo
    j = unique( which(!is.finite( uu ) | uu < cthreshold ) )
    if ( ( length(i) == 0 ) & (length(j) == 0 ) ) finished=TRUE
    if (length(i) > 0 ) m = m[ -i , ]
    if (length(j) > 0 ) m = m[ , -j ]
  }
 
  return(m)
}



predator_prey_matrix = function(rdt, cthreshold=0) {
  # create (id = auid,time) X (sp= prey species) matrix
    
  m = data.table::dcast( setDT(rdt), 
    formula = preyspeccd ~ spec, value.var="zscore", 
    fun.aggregate=mean, fill=NA, drop=FALSE, na.rm=TRUE
  )  # mean is just to keep dcast happy

  id = paste( m$preyspeccd, sep=".")

  m$preyspeccd = NULL
  sps = colnames(m)

  m = as.matrix(m[]) 
  dimnames(m)[[1]] = id
 
  # remove low counts (absence) in the timeseries  .. species (cols) only
  # cthreshold = 0.005  # 0.01%  -> 97;  0.05 => 44 species; 0.001 => 228 species; 0.005 -> 146 
  m [ m < cthreshold ] = 0
  m [ !is.finite(m) ] = 0

  # reduce:: remove taxa/locations with very low counts (~cthreshold)
  finished = FALSE
  while( !(finished) ) {
    oo = colSums( ifelse( is.finite(m), 1, 0 ), na.rm=TRUE)
    i = unique( which(rowSums(m, na.rm=TRUE) == 0 ) )
    uu = colSums(m, na.rm=TRUE)/oo
    j = unique( which(!is.finite( uu ) | uu < cthreshold ) )
    if ( ( length(i) == 0 ) & (length(j) == 0 ) ) finished=TRUE
    if (length(i) > 0 ) m = m[ -i , ]
    if (length(j) > 0 ) m = m[ , -j ]
  }
 
  return(m)
}

 


taxa_to_code = function(tx){

	idval = taxonomy.recode( from="taxa.fast", tolookup=tx ) # look up species id from both itis and local taxonomy.db
	taxalist = taxonomy.recode( from="spec", tolookup=idval )
	if (length(idval) > 1){
		message("There is more than one possible match. Review the list and use the appropriate 'spec' numeric code instead." )
	} else if (length(idval)==1) {
		message("A unique match. If incorrect, send a more generic name." )
		print( taxalist )
	} else{
		message("No matches in local list, trying ITIS:")
		idval = taxonomy.recode( from="taxa", tolookup=tx ) # look up species id from both itis and local taxonomy.db
		taxalist = taxonomy.recode( from="tsn", to="taxa", tolookup= idval[[1]]$tsn )
	}
	print( taxalist )
	return(idval)
}


my_predator = function(  diet, tx ) {

	if (is.numeric(tx)) {
		idval = tx
	} else {
		idval = taxa_to_code(tx)
		if (length(idval) != 1) stop( "Review options and re-enter another name or enter 'spec' code directly" )
	}

	predator = diet[ preyspeccd==idval, ]

  if (nrow(predator) == 0) {
    message("No information of predator items found.")
    return(NULL)
  }

  predator_summary = predator[, .(nobs=.N, 
      meanpreywgt=mean(pwt, na.rm=TRUE), 
      meanpreynum=mean(pnum, na.rm=TRUE),
      meanpredatorwgt=mean(fwt, na.rm=TRUE),
      meantemp=mean(bottom_temperature, na.rm=TRUE),
      meandepth=mean(depth, na.rm=TRUE)
    ), by="spec" ]
	 
	predator_names = taxonomy.recode( from="spec", tolookup=predator_summary$spec )
	setDT(predator_names)
  predator_summary = predator_names[predator_summary, on="spec"]
  setorder(predator_summary, -nobs)
  print(predator_summary)
	attr(predator, "predator_summary") = predator_summary
 
	return(predator)
}



my_prey = function(  diet, tx, preyvector=NULL ) {

	if (is.numeric(tx)) {
		idval = tx
	} else {
		idval = taxa_to_code(tx)
		if (length(idval) != 1) stop( "Review options and re-enter another name or enter 'spec' code directly" )
	}

	prey = diet[ spec==idval, ]
  
  if (nrow(prey) == 0) {
    message("No information of prey items found.")
    return(NULL)
  }

  prey_summary = prey[, .(nobs=.N, 
      meanpreywgt=mean(pwt, na.rm=TRUE), 
      meanpredatorwgt=mean(fwt, na.rm=TRUE),
      meantemp=mean(bottom_temperature, na.rm=TRUE),
      meandepth=mean(depth, na.rm=TRUE)
    ), by="preyspeccd" ]
	 
	prey_names = taxonomy.recode( from="spec", tolookup=prey_summary$preyspeccd )
	setDT(prey_names)
  setnames(prey_names, "spec", "preyspeccd")
  prey_summary = prey_names[prey_summary, on="preyspeccd"]
  setorder(prey_summary, -nobs)
  
  print(prey_summary)
	attr(prey, "prey_summary") = prey_summary
 
  if (!is.null(preyvector)) {
    
    pv = as.numeric(preyvector)

    pr = prey[ which(is.finite(preyspeccd)),.(Nprey=sum(Nprey, na.rm=TRUE)), by=.(preyspeccd, id)]   
    pr = pr[CJ(preyspeccd=pv, id=unique(pr$id), unique=TRUE  ), on=.(preyspeccd, id) ]

    m = data.table::dcast( setDT(pr), 
      formula =  id ~ preyspeccd, value.var="Nprey", 
      fun.aggregate=mean, fill=NA, drop=FALSE, na.rm=TRUE
    )  # mean is just to keep dcast happy

    id = m$id 
    m$id = NULL 
    setcolorder(m, preyvector )

    m = as.matrix(m)
    m[!is.finite(m)] = 0

    return( m )
  }

	return(prey)
}


