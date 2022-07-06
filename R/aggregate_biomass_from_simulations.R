
aggregate_biomass_from_simulations = function( sims, sppoly,  yrs, fn=NULL, method="sum", redo=FALSE ) {

  if ( ! redo ) {
    SM = NULL
    if (!is.null( fn) ) {
      if (file.exists(fn)) load( fn)
      return( SM )
    }
  }

  SM =  list()

  if (method =="sum") {
    SM$cfaall    = colSums( sims * sppoly$au_sa_km2, na.rm=TRUE )
    SM$cfanorth  = colSums( sims * sppoly$cfanorth_surfacearea, na.rm=TRUE )
    SM$cfasouth  = colSums( sims * sppoly$cfasouth_surfacearea, na.rm=TRUE )
    SM$cfa23     = colSums( sims * sppoly$cfa23_surfacearea, na.rm=TRUE )
    SM$cfa24     = colSums( sims * sppoly$cfa24_surfacearea, na.rm=TRUE )
    SM$cfa4x     = colSums( sims * sppoly$cfa4x_surfacearea, na.rm=TRUE )
  }

  if (method == "mean") {
    SM$cfaall    = colSums( sims * sppoly$au_sa_km2/ sum(sppoly$au_sa_km2), na.rm=TRUE )
    SM$cfanorth  = colSums( sims * sppoly$cfanorth_surfacearea/ sum(sppoly$cfanorth_surfacearea), na.rm=TRUE )
    SM$cfasouth  = colSums( sims * sppoly$cfasouth_surfacearea/ sum(sppoly$cfasouth_surfacearea), na.rm=TRUE )
    SM$cfa23     = colSums( sims * sppoly$cfa23_surfacearea/ sum(sppoly$cfa23_surfacearea), na.rm=TRUE )
    SM$cfa24     = colSums( sims * sppoly$cfa24_surfacearea/ sum(sppoly$cfa24_surfacearea), na.rm=TRUE )
    SM$cfa4x     = colSums( sims * sppoly$cfa4x_surfacearea/ sum(sppoly$cfa4x_surfacearea), na.rm=TRUE )
  }

  RES = data.frame( yrs = yrs )
  RES$cfaall = apply( simplify2array(SM[["cfaall"]]), 1, mean )
  RES$cfaall_sd = apply( simplify2array(SM[["cfaall"]]), 1, sd )
  RES$cfaall_median = apply( simplify2array(SM[["cfaall"]]), 1, median )
  RES$cfaall_lb = apply( simplify2array(SM[["cfaall"]]), 1, quantile, probs=0.025 )
  RES$cfaall_ub = apply( simplify2array(SM[["cfaall"]]), 1, quantile, probs=0.975 )

  RES$cfanorth = apply( simplify2array(SM[["cfanorth"]]), 1, mean )
  RES$cfanorth_sd = apply( simplify2array(SM[["cfanorth"]]), 1, sd )
  RES$cfanorth_median = apply( simplify2array(SM[["cfanorth"]]), 1, median )
  RES$cfanorth_lb = apply( simplify2array(SM[["cfanorth"]]), 1, quantile, probs=0.025 )
  RES$cfanorth_ub = apply( simplify2array(SM[["cfanorth"]]), 1, quantile, probs=0.975 )

  RES$cfasouth = apply( simplify2array(SM[["cfasouth"]]), 1, mean )
  RES$cfasouth_sd = apply( simplify2array(SM[["cfasouth"]]), 1, sd )
  RES$cfasouth_median = apply( simplify2array(SM[["cfasouth"]]), 1, median )
  RES$cfasouth_lb = apply( simplify2array(SM[["cfasouth"]]), 1, quantile, probs=0.025 )
  RES$cfasouth_ub = apply( simplify2array(SM[["cfasouth"]]), 1, quantile, probs=0.975 )

  RES$cfa23 = apply( simplify2array(SM[["cfa23"]]), 1, mean )
  RES$cfa23_sd = apply( simplify2array(SM[["cfa23"]]), 1, sd )
  RES$cfa23_median = apply( simplify2array(SM[["cfa23"]]), 1, median )
  RES$cfa23_lb = apply( simplify2array(SM[["cfa23"]]), 1, quantile, probs=0.025 )
  RES$cfa23_ub = apply( simplify2array(SM[["cfa23"]]), 1, quantile, probs=0.975 )

  RES$cfa24 = apply( simplify2array(SM[["cfa24"]]), 1, mean )
  RES$cfa24_sd = apply( simplify2array(SM[["cfa24"]]), 1, sd )
  RES$cfa24_median = apply( simplify2array(SM[["cfa24"]]), 1, median )
  RES$cfa24_lb = apply( simplify2array(SM[["cfa24"]]), 1, quantile, probs=0.025 )
  RES$cfa24_ub = apply( simplify2array(SM[["cfa24"]]), 1, quantile, probs=0.975 )

  RES$cfa4x = apply( simplify2array(SM[["cfa4x"]]), 1, mean )
  RES$cfa4x_sd = apply( simplify2array(SM[["cfa4x"]]), 1, sd )
  RES$cfa4x_median = apply( simplify2array(SM[["cfa4x"]]), 1, median )
  RES$cfa4x_lb = apply( simplify2array(SM[["cfa4x"]]), 1, quantile, probs=0.025 )
  RES$cfa4x_ub = apply( simplify2array(SM[["cfa4x"]]), 1, quantile, probs=0.975 )

  SM$RES = RES 

  if (!is.null(fn)) save( SM, file=fn, compress=TRUE )

  return(SM)
}
