 
  get.fishery.stats.by.region = function( Reg="cfaall", y=NULL ) {
    
    landings = snowcrab_landings_db()
    
    if (is.null(y)) y = sort(unique(landings$yr) )
    out = data.frame( yr=y )

    if (Reg=="cfaall")  region = sort( unique(landings$cfa) ) # all data
    if (Reg=="cfanorth") region = c("cfa20", "cfa21", "cfa22", "cfanorth", "north")
    if (Reg=="cfasouth") region = c("cfa23", "cfa24", "cfasouth", "cfaslope")
    if (Reg=="cfa4x") region = "cfa4x"
    if (Reg=="cfa23") region = "cfa23"
    if (Reg=="cfa24") region = "cfa24"
    

    lnd = landings[ which(landings$cfa %in% region) ,]
    #lnd = landings[ which(landings$cfa0 %in% region) ,]
    l = aggregate( lnd$landings, list(yr=lnd$yr), function(x) sum(x, na.rm=T))
    names(l) = c("yr", "landings")
    
    out = merge(out, l, by="yr", all.x=T, all.y=F, sort=T)

    lnd$cpue_direct = lnd$landings / lnd$effort
    lnd$cpue_direct[ which( lnd$cpue_direct > (650*0.454)) ] = NA  # same rule as in snowcrab_landings_db -- 650lbs/trap is a reasonable upper limit

    cpue = aggregate( lnd$cpue, list(yr=lnd$yr), function(x) mean(x, na.rm=T))
    names(cpue) = c("yr", "cpue")
    
    # force known errors to be overwritten: TODO need to check why it deviates (think landings occurred in Gulf Region for a few years and they were lost from Marfis)
    if (Reg=='cfanorth') {
      message( "Forcing CPUE in 2013 and 2014. Likely due to missing landings in MARFIS --\n landings occurred in Gulf REGION -- this means landings are likely off too for those years. \n TODO: fix this")
      ii = which(res$yr==2014)
      if (length(ii)==1) res[ii,'cpue'] = 104.5
      
      ii=which(res$yr==2013)
      if (length(ii)==1) res[ii,'cpue'] <- 106
    }

    out = merge (out, cpue, by="yr", all.x=T, all.y=F, sort=T)
    
    out$effort = out$landings / out$cpue  ## estimate effort level as direct estimates are underestimates (due to improper logbook records)
    rownames(out) = out$yr
    
    return(out)
    
  }


