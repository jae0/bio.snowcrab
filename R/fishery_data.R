fishery_data = function(  
    toget="",
    regions = list( region=c("cfanorth", "cfasouth", "cfa4x") ), 
    # region_label = c("N-ENS", "S-ENS", "4X"),
    yrs=NULL, 
    redo=FALSE
    ) {

    # aggregated fishery data
 
    out = NULL
    
    if (!is.list(regions) ) {
        message("Expecting a named list where the name identifies the variable to filter upon, eg:\n" )
        message("  'region=c(...)' for cfanorth, cfasouth, cfa4x, in region (default)")
        message("  'subarea=c(...)' for cfa23, cfa24, in subarea")
        message("  'marfis_area=c(...)' in marfis designations")
        regions = list(region = regions)
    }

    vn = names(regions)

    fn = file.path( project.datadirectory("bio.snowcrab"), "output", paste("fishery_data_summary_by_", vn, ".RDS", sep="") )

    if (!redo) {
        if (file.exists(fn)) out = readRDS(fn)
        if (toget != "") out = out[[toget]]
        return(out)
    }

    # 1978 to present
    fstat = setDT( logbook.db(DS="logbook")  )

    fstat = fstat[ fstat[[vn]] %in% regions[[vn]] , ]

    # note: "yr" is fishing year, in 4x: 1999-2000 is yr=1999
    if (is.null(yrs)) yrs = min(fstat$yr, na.rm=TRUE) : max(fstat$yr, na.rm=TRUE)
    fstat = fstat[ yr %in% yrs, ]  
       

    # Y  & fraction_observed 
    Y = fstat[ , .(
        landings = sum(landings, na.rm=TRUE),
        cpue = mean(cpue, na.rm=TRUE),
        cpue_sd=sd(cpue, na.rm=TRUE)
    ), by=c("yr", vn) ]

    Y$effort = Y$landings / Y$cpue  ## estimate effort level as direct estimates are underestimates (due to improper logbook records)
    Y = Y[order(get(vn), yr),]

    Y$landings = Y$landings / 1000  # t
    Y$effort = Y$effort / 1000  # 1000 th

    # please update this (it is embeded in a function right now but needs to be moved to a CSV or DB lookup)
    tacs = snowcrab_tacs( vn ) 
    setDT( tacs )

    Y = merge(Y, tacs, by=c(vn, "yr"), all.x=TRUE, all.y=FALSE)

    Y$landings = round(Y$landings )
    Y$effort = round(Y$effort, 1 )
    Y$cpue = round(Y$cpue )
    Y$TAC = round(Y$TAC)
    Y$Licenses = as.numeric(Y$Licenses)
    Y = Y[order(Y[[vn]], Y[["yr"]] ),]

    out[["summary_annual"]] = Y
    
    
    # summary_monthly (monthly stats)
    M = NA
    
    fstat$month = month( as.POSIXct(fstat$date.landed) )

    M = fstat[, .(
        landings=sum(landings, na.rm=TRUE), 
        cpue=mean(cpue, na.rm=TRUE),
        cpue_sd=sd(cpue, na.rm=TRUE)
    ), by=c(vn, "yr", "month") ]
        
    M = M[ which(!is.na(M[[vn]] )) , ]

    M$effort =  M$landings / M$cpue  ## estimate effort level as direct estimates are underestimates (due to improper logbook records)
#       setDF(M)

    M$landings = M$landings / 1000  # t
    M$effort = M$effort / 1000  # 1000 th

    M$landings = round(M$landings )
    M$effort = round(M$effort, 1 )
    M$cpue = round(M$cpue )

    # add breaks in aug for 4x to plot correctly
    M1 = CJ( R="cfa4x", month=8, yr=yrs, cpue=NA, landings=NA, effort=NA, cpue_sd=NA)
    setnames(M1, "R", vn )  # set to correct variable name
    M = rbind( M, M1 )

    # add breaks  for north to plot correctly: month=6 (separate spring-fall)
    M2 = CJ( R="cfanorth", month=6.5, yr=yrs, cpue=NA, landings=NA, effort=NA, cpue_sd=NA)
    setnames(M2, "R", vn )  # set to correct variable name
    M = rbind( M, M2 )
    M = M[ order(get(vn), yr, month),]
    out[["summary_monthly"]] = M
 

    # summary_weekly stats

    W = NA
    
    fstat$week = week( as.POSIXct(fstat$date.landed) )

    W = fstat[, .(
        landings=sum(landings, na.rm=TRUE), 
        cpue=mean(cpue, na.rm=TRUE),
        cpue_sd=sd(cpue, na.rm=TRUE)
    ), by=c(vn, "yr", "week") ]
        
    W = W[ which(!is.na(W[[vn]])) , ]

    W$effort =  W$landings / W$cpue  ## estimate effort level as direct estimates are underestimates (due to improper logbook records)
#       setDF(W)

    W$landings = W$landings / 1000  # t
    W$effort = W$effort / 1000  # 1000 th

    W$landings = round(W$landings )
    W$effort = round(W$effort, 1 )
    W$cpue = round(W$cpue )

    # add breaks in aug for 4x to plot correctly
    W1 = CJ( R="cfa4x", week=floor(8/12*52)+0.5, yr=yrs, cpue=NA, landings=NA, effort=NA, cpue_sd=NA)
    setnames(W1, "R", vn )  # set to correct variable name
    W = rbind( W, W1)

# add breaks in aug for north to plot correctly
    W2 = CJ( R="cfanorth", week=26.5, yr=yrs, cpue=NA, landings=NA, effort=NA, cpue_sd=NA)
    setnames(W2, "R", vn )  # set to correct variable name
    W = rbind( W, W2)

    W = W[ order(get(vn), yr, week),]
    out[["summary_weekly"]] = W
 

    # summary_biweekly stats

    WW = NA
    
    fstat$week = floor( week( as.POSIXct(fstat$date.landed) ) / 2 )  * 2

    WW = fstat[, .(
        landings=sum(landings, na.rm=TRUE), 
        cpue=mean(cpue, na.rm=TRUE),
        cpue_sd=sd(cpue, na.rm=TRUE)
    ), by=c(vn, "yr", "week") ]
        
    WW = WW[ which(!is.na(WW[[vn]])) , ]

    WW$effort =  WW$landings / WW$cpue  ## estimate effort level as direct estimates are underestimates (due to improper logbook records)
#       setDF(WW)

    WW$landings = WW$landings / 1000  # t
    WW$effort = WW$effort / 1000  # 1000 th

    WW$landings = round(WW$landings )
    WW$effort = round(WW$effort, 1 )
    WW$cpue = round(WW$cpue )

    # add breaks in aug for 4x to plot correctly
    WW1 = CJ( R="cfa4x", week=floor(8/12*52), yr=yrs, cpue=NA, landings=NA, effort=NA, cpue_sd=NA)
    setnames(WW1, "R", vn )  # set to correct variable name
    WW = rbind( WW, WW1)

# add breaks in aug for north to plot correctly
    WW2 = CJ( R="cfanorth", week=26, yr=yrs, cpue=NA, landings=NA, effort=NA, cpue_sd=NA)
    setnames(WW2, "R", vn )  # set to correct variable name
    WW = rbind( WW, WW2)

    WW = WW[ order(get(vn), yr, week),]
    out[["summary_biweekly"]] = W
 

     # shell condition
    male = 0
    odb = observer.db("odb")
    setDT(odb)
    odb[[vn]] = NA
    for ( reg in regions[[vn]] ) {
        r = polygon_inside(x=odb, region=aegis.polygons::polygon_internal_code(reg), planar=FALSE)
        if (length(r)> 0)  odb[[vn]][r] = reg
    } 
    odb = odb[ !is.na(get(vn)) , ]
 
    odb = odb[ sex==male & cw >= 95 & cw < 170 & prodcd_id=="0" & is.finite(shell)  ,]  # commerical sized crab only

    shell_condition = odb[ , .N, by=c(vn, "fishyr", "shell") ]
    shell_condition[, total:=sum(N, na.rm=TRUE), by=c(vn, "fishyr")]
    shell_condition$percent = round(shell_condition$N / shell_condition$total, 3) * 100
    out[["shell_condition"]] = shell_condition


    # fraction_observed 

    if (!exists("totmass_kept", odb)) odb$totmass_kept = odb$totmass  # temporary fix until observer db raw data gets refreshed (JC 2022)
    message("Make sure all observer local raw tables have been refreshed")
    
    fraction_observed = odb[ !is.na(odb[[vn]]), .(no_traps=.N, kept=sum(totmass_kept, na.rm=TRUE), caught=sum(totmass, na.rm=TRUE) ), by=c(vn, "fishyr") ]
    setnames(fraction_observed, "fishyr", "yr")

    fraction_observed = fraction_observed[ out$summary_annual, on=c(vn, "yr")]
    # in percentages:

    fraction_observed$observed_landings_pct = round( fraction_observed$kept / fraction_observed$landings *100, 2)  # values in metric tonnes
    fraction_observed$observed_effort_pct = round( (fraction_observed$no_traps/1000) / fraction_observed$effort * 100,2)  #  effort is in 1000 th

    fraction_observed$observed_landings_pct[ which(!is.finite(fraction_observed$observed_landings_pct)) ] = 0
    fraction_observed$observed_effort_pct[ which(!is.finite(fraction_observed$observed_effort_pct)) ] = 0

    out[["fraction_observed"]] = fraction_observed 

    saveRDS( out, file=fn, compress=TRUE ) 

    return(out)

}

