
  observer_summarize = function( toget=c("shell_condition", "fraction_observed"), regions=c("cfanorth", "cfasouth", "cfa4x") ) {
    # carapace condition .. could have used the stuff in 02_fisures and table but it was convoluted
    # should replace that mess with this eventually 
    out = list()

    male = 0
    odb = observer.db("odb")
    setDT(odb)
    odb = odb[ which( odb$sex==male & odb$cw >= 95 & odb$cw < 170 & odb$prodcd_id=="0" & is.finite(odb$shell) ) ,]  # commerical sized crab only
    odb$region = NA
    for ( reg in regions ) {
      r = polygon_inside(x=odb, region=aegis.polygons::polygon_internal_code(reg), planar=FALSE)
      if (length(r)> 0)  odb$region[r] = reg
    }

    if (any( grepl("shell_condition", toget)) ) {
      shell_condition = odb[ !is.na(odb$region), .N, by=.(region, fishyr, shell) ]
      shell_condition[, total:=sum(N, na.rm=TRUE), by=.(region, fishyr)]
      shell_condition$percent = round(shell_condition$N / shell_condition$total, 3) * 100
      out$shell_condition=shell_condition
    }

    if (any( grepl("fraction_observed", toget)) ) {

      FD = fisherydata_summary()
      setDT(FD)

      if (!exists("totmass_kept", odb)) odb$totmass_kept = odb$totmass  # temporary fix until observer db raw data gets refreshed (JC 2022)
      message("Make sure all observer local raw tables have been refreshed")
      
      fraction_observed = odb[ !is.na(odb$region), .(no_traps=.N, kept=sum(totmass_kept, na.rm=TRUE), caught=sum(totmass, na.rm=TRUE) ), by=.(region, fishyr) ]
      setnames(fraction_observed, "fishyr", "yr")

      fraction_observed = fraction_observed[ FD, on=.(region, yr)]
      # in percentages:

      fraction_observed$observed_landings_pct = round( fraction_observed$kept / fraction_observed$landings *100, 2)  # values in metric tonnes
      fraction_observed$observed_effort_pct = round( (fraction_observed$no_traps/1000) / fraction_observed$effort * 100,2)  #  effort is in 1000 th

      fraction_observed$observed_landings_pct[ which(!is.finite(fraction_observed$observed_landings_pct)) ] = 0
      fraction_observed$observed_effort_pct[ which(!is.finite(fraction_observed$observed_effort_pct)) ] = 0

      out$fraction_observed = fraction_observed 

    }

    return(out)
  }

