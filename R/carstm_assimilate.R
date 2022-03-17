
carstm_assimilate = function( pN=NULL, pW=NULL, pH=NULL, pB=NULL, sppoly=NULL, redo=FALSE ) {

  operation = NULL
  
  if (!is.null(pN)) operation = c( operation, "number" )
  if (!is.null(pW)) operation = c( operation, "meansize" )
  if (!is.null(pH)) operation = c( operation, "presence_absence" )
  if (!is.null(pB)) operation = c( operation, "biomass" )
 
  if (is.null(operation)) {
    if (!redo) stop( "The operation is ambiguous, check parameters." )
  }

  if (is.null(sppoly)) sppoly = areal_units( p=p )
  areal_units_fn = attributes(sppoly)[["areal_units_fn"]]

  aufns = carstm_filenames( p=p, returntype="carstm_modelled_fit", areal_units_fn=areal_units_fn )
  # same file naming as in carstm ..
  outputdir = dirname( aufns )

  if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

  fn_no  =  paste( gsub(".rdata", "", aufns), "sims_number", "rdata", sep="." )
  fn_bio =  paste( gsub(".rdata", "", aufns), "sims_biomass", "rdata", sep="." )
  fn_pa  =  paste( gsub(".rdata", "", aufns), "sims_pa", "rdata", sep="." )
  fn_wgts  =  paste( gsub(".rdata", "", aufns), "sims_meanweight", "rdata", sep="." )

  nums = NULL
  biom = NULL
  pa = NULL
  wgts = NULL

  if (!redo) {
    if ( DS=="number" ) {
      if (file.exists(fn_no)) load( fn_no )
      return( nums )
    }

    if ( DS=="biomass" ) {
      if (file.exists(fn_bio)) load( fn_bio )
      return( biom )
    }

    if ( DS=="presence_absence")  {
      if (file.exists(fn_pa)) load( fn_pa )
      return( pa )
    }

    if ( DS=="meansize" ) {
      if (file.exists(fn_wgts)) load( fn_wgts )
      return( wgts )
    }
  }



  # construct meansizes matrix used to convert number to weight
  if ( "presence_absence" %in% operation ) {
    pa = carstm_model( p=pH, DS="carstm_modelled_summary", sppoly=sppoly  )
    pa = pa[[ "predictions_posterior_simulations" ]]
    pa[!is.finite(pa)] = NA
    # pa = inverse.logit(pa)
    # pa[!is.finite(pa)] = NA
    attr( pa, "unit") = "probability"
    save( pa, file=fn_pa, compress=TRUE )
    if (length(operation)==1) return(pa)
  }


  if ( "meansize" %in% operation) {
    wgts = carstm_model( p=pW, DS="carstm_modelled_summary", sppoly=sppoly  )
    wgts = wgts[[ "predictions_posterior_simulations"  ]]
    wgts[!is.finite(wgts)] = NA
    if (exists("wgts_max", pW)) {
      i = which( wgts > pW$wgts_max )
      if (length(i) > 0 ) wgts[ i ] = pW$wgts_max
    }
    attr( wgts, "unit") = "kg"
    save( wgts, file=fn_wgts, compress=TRUE )
    if (length(operation)==1) return(wgts)
  }


  if ( "biomass" %in% operation) {
    biom = carstm_model( p=pB, DS="carstm_modelled_summary", sppoly=sppoly  )
    biom = biom[[ "predictions_posterior_simulations" ]] 
    biom[!is.finite(biom)] = NA
    if (exists("B_max", pB)) {
      i = which( biom > pB$B_max )
      if (length(i) > 0 ) biom[ i ] = pB$B_max
    }
    # biom = biom / 10^6  # kg / km^2 -> kt / km^2
    attr(biom, "unit") = "kg/km^2"
    save( biom, file=fn_bio, compress=TRUE )
    if (length(operation)==1) return(biom)
  }


  if ( "number" %in% operation ) {
    nums = carstm_model( p=pN, DS="carstm_modelled_summary", sppoly=sppoly  )
    nums = nums[[ "predictions_posterior_simulations" ]]    
    nums[!is.finite(nums)] = NA
    if (exists("N_max", pN)) {
      i = which( nums > pN$N_max )
      if (length(i) > 0 ) nums[ i ] = pN$N_max
    }
    # nums = nums / 10^6  # n/km2 ->  M n  / km^2
    attr(nums, "unit") = "n/km^2"
    save( nums, file=fn_no, compress=TRUE )
    if (length(operation)==1) return(nums)
  }

 
  if (length(operation) == 2 ) {
    if ("meansize" %in% operation ) {
      if ( "biomass" %in% operation ) {
        nums = biom / wgts   # (n * 10^6) / km^2
        attr(nums, "unit") = "n/km^2"
        return(nums)
      }
      if ( "number" %in% operation ) {
        # numerical: density no / m^2  -->>  (no. km^2)
        biom = nums * wgts # * 10^6 / 10^6  # cancels out .. kg / km^2 -> kt / km^2
        attr(biom, "unit") = "kg/km^2"
        return(biom)
      }
    } 
    
    if ("pa" %in% operation) {
      if ( "biomass" %in% operation ) {
        biom = biom * pa
        attr(biom, "unit") = "kg/km^2"
        return(biom)
      }
      if ("number" %in% operation) {
        nums = nums * pa
        attr(nums, "unit") = "n/km^2"
        return(nums)
      }
    }
  }

  if (length(operation) == 3 ) {
    if ("meansize" %in% operation ) {
      if ( "biomass" %in% operation ) {
        if ( "presence_absence" %in% operation ) {
          nums = biom*pa / wgts
          return(nums)
        }
      }
      if ( "number" %in% operation ) {
        if ( "presence_absence" %in% operation ) {
          biom = nums*pa*wgts
          return(biom)
        }
      }
    }
  }

  return(operation)


}
