fishery_data = function(  
    mau = "region",  # management areal unit variable name to filter upon
    toget="",
    yrs=NULL, 
    redo=FALSE
    ) {

 
    # aggregated fishery data
 
    out = NULL

    fn = file.path( project.datadirectory("bio.snowcrab"), "output", paste("fishery_data_summary_by_", mau, ".rdz", sep="") )

    if (!redo) {
        if (file.exists(fn)) out = read_write_fast(fn)
        if (toget != "") out = out[[toget]]
        return(out)
    }

    maus = management_areal_units(mau)  # labels etc


    # 1978 to present
    fstat = setDT( logbook.db(DS="logbook")  )

    fstat = fstat[ fstat[[mau]] %in% maus[["internal"]], ]

    # note: "yr" is fishing year, in 4x: 1999-2000 is yr=1999
    if (is.null(yrs)) yrs = min(fstat$yr, na.rm=TRUE) : max(fstat$yr, na.rm=TRUE)
    fstat = fstat[ yr %in% yrs, ]  
       

    # Y  & fraction_observed 
    Y = fstat[ , .(
        landings = sum(landings, na.rm=TRUE),
        cpue = mean(cpue, na.rm=TRUE),
        cpue_sd=sd(cpue, na.rm=TRUE)
    ), by=c("yr", mau) ]

    Y$effort = Y$landings / Y$cpue  ## estimate effort level as direct estimates are underestimates (due to improper logbook records)
    Y = Y[order(get(mau), yr),]

    Y$landings = Y$landings / 1000  # t
    Y$effort = Y$effort / 1000  # 1000 th

    # please update this (it is embeded in a function right now but needs to be moved to a CSV or DB lookup)
    tacs = snowcrab_tacs( mau ) 
    setDT( tacs )

    Y = merge(Y, tacs, by=c(mau, "yr"), all.x=TRUE, all.y=FALSE)

    Y$landings = round(Y$landings )
    Y$effort = round(Y$effort, 1 )
    Y$cpue = round(Y$cpue )
    Y$TAC = round(Y$TAC)
    Y$Licenses = as.numeric(Y$Licenses)
    Y = Y[order(Y[[mau]], Y[["yr"]] ),]

    out[["summary_annual"]] = Y
    
    
    # summary_monthly (monthly stats)
    M = NA
    
    fstat$month = month( as.POSIXct(fstat$date.landed) )

    M = fstat[, .(
        landings=sum(landings, na.rm=TRUE), 
        cpue=mean(cpue, na.rm=TRUE),
        cpue_sd=sd(cpue, na.rm=TRUE)
    ), by=c(mau, "yr", "month") ]
        
    M = M[ which(!is.na(M[[mau]] )) , ]

    M$effort =  M$landings / M$cpue  ## estimate effort level as direct estimates are underestimates (due to improper logbook records)
#       setDF(M)

    M$landings = M$landings / 1000  # t
    M$effort = M$effort / 1000  # 1000 th

    M$landings = round(M$landings )
    M$effort = round(M$effort, 1 )
    M$cpue = round(M$cpue )

    # add breaks in aug for 4x to plot correctly
    M1 = CJ( R="cfa4x", month=8, yr=yrs, cpue=NA, landings=NA, effort=NA, cpue_sd=NA)
    setnames(M1, "R", mau )  # set to correct variable name
    M = rbind( M, M1 )

    # add breaks  for north to plot correctly: month=6 (separate spring-fall)
    M2 = CJ( R="cfanorth", month=6.5, yr=yrs, cpue=NA, landings=NA, effort=NA, cpue_sd=NA)
    setnames(M2, "R", mau )  # set to correct variable name
    M = rbind( M, M2 )
    M = M[ order(get(mau), yr, month),]
    out[["summary_monthly"]] = M
 

    # summary_weekly stats

    W = NA
    
    fstat$week = week( as.POSIXct(fstat$date.landed) )

    W = fstat[, .(
        landings=sum(landings, na.rm=TRUE), 
        cpue=mean(cpue, na.rm=TRUE),
        cpue_sd=sd(cpue, na.rm=TRUE)
    ), by=c(mau, "yr", "week") ]
        
    W = W[ which(!is.na(W[[mau]])) , ]

    W$effort =  W$landings / W$cpue  ## estimate effort level as direct estimates are underestimates (due to improper logbook records)
#       setDF(W)

    W$landings = W$landings / 1000  # t
    W$effort = W$effort / 1000  # 1000 th

    W$landings = round(W$landings )
    W$effort = round(W$effort, 1 )
    W$cpue = round(W$cpue )

    # add breaks in aug for 4x to plot correctly
    W1 = CJ( R="cfa4x", week=floor(8/12*52)+0.5, yr=yrs, cpue=NA, landings=NA, effort=NA, cpue_sd=NA)
    setnames(W1, "R", mau )  # set to correct variable name
    W = rbind( W, W1)

# add breaks in aug for north to plot correctly
    W2 = CJ( R="cfanorth", week=26.5, yr=yrs, cpue=NA, landings=NA, effort=NA, cpue_sd=NA)
    setnames(W2, "R", mau )  # set to correct variable name
    W = rbind( W, W2)

    W = W[ order(get(mau), yr, week),]
    out[["summary_weekly"]] = W
  
    # summary_biweekly stats

    WW = NA
    
    fstat$week = floor( week( as.POSIXct(fstat$date.landed) ) / 2 )  * 2

    WW = fstat[, .(
        landings=sum(landings, na.rm=TRUE), 
        cpue=mean(cpue, na.rm=TRUE),
        cpue_sd=sd(cpue, na.rm=TRUE)
    ), by=c(mau, "yr", "week") ]
        
    WW = WW[ which(!is.na(WW[[mau]])) , ]

    WW$effort =  WW$landings / WW$cpue  ## estimate effort level as direct estimates are underestimates (due to improper logbook records)
#       setDF(WW)

    WW$landings = WW$landings / 1000  # t
    WW$effort = WW$effort / 1000  # 1000 th

    WW$landings = round(WW$landings )
    WW$effort = round(WW$effort, 1 )
    WW$cpue = round(WW$cpue )

    # add breaks in aug for 4x to plot correctly
    WW1 = CJ( R="cfa4x", week=floor(8/12*52), yr=yrs, cpue=NA, landings=NA, effort=NA, cpue_sd=NA)
    setnames(WW1, "R", mau )  # set to correct variable name
    WW = rbind( WW, WW1)

# add breaks in aug for north to plot correctly
    WW2 = CJ( R="cfanorth", week=26, yr=yrs, cpue=NA, landings=NA, effort=NA, cpue_sd=NA)
    setnames(WW2, "R", mau )  # set to correct variable name
    WW = rbind( WW, WW2)

    WW = WW[ order(get(mau), yr, week),]
    out[["summary_biweekly"]] = W
 

     # shell condition
    male = 0
    odb = observer.db("odb")
    setDT(odb)

    if (!exists(mau, odb)) {
        odb[[mau]] = NA
        for ( reg in maus[["internal"]] ) {
            r = polygon_inside(x=odb, region=reg, planar=FALSE)
            if (length(r)> 0)  odb[[mau]][r] = reg
        }
    }

    odb = odb[ !is.na(get(mau)) , ]
 
    odb = odb[ sex==male & cw >= 95 & cw < 170 & prodcd_id=="0" & is.finite(shell)  ,]  # commerical sized crab only

    shell_condition = odb[ , .N, by=c(mau, "fishyr", "shell") ]
    shell_condition[, total:=sum(N, na.rm=TRUE), by=c(mau, "fishyr")]
    shell_condition$percent = round(shell_condition$N / shell_condition$total, 3) * 100
    out[["shell_condition"]] = shell_condition


    # fraction_observed 

    if (!exists("totmass_kept", odb)) odb$totmass_kept = odb$totmass  # temporary fix until observer db raw data gets refreshed (JC 2022)
    message("Make sure all observer local raw tables have been refreshed")

    stop( "not finished" )

    fraction_observed = odb[ 
        !is.na(odb[[mau]]), 
        .(  no_traps= sum(num_hook_haul, na.rm=TRUE) ,
            kept=sum(totmass_kept, na.rm=TRUE), 
            caught=sum(totmass, na.rm=TRUE) ,
            trips = length(unique(trip))
        ), 
        by=c(mau, "fishyr") 
    ]
    
    setnames(fraction_observed, "fishyr", "yr")

    fraction_observed = fraction_observed[ out$summary_annual, on=c(mau, "yr")]
    
    # in percentages:
message("NOTE: units are not right for observed fraction by weight ... check all fields")


    fraction_observed$observed_landings_pct = round( fraction_observed$kept / fraction_observed$landings *100, 2)  # values in metric tonnes
    fraction_observed$observed_effort_pct = round( (fraction_observed$no_traps/1000) / fraction_observed$effort * 100,2)  #  effort is in 1000 th

    fraction_observed$observed_landings_pct[ which(!is.finite(fraction_observed$observed_landings_pct)) ] = 0
    fraction_observed$observed_effort_pct[ which(!is.finite(fraction_observed$observed_effort_pct)) ] = 0

    out[["fraction_observed"]] = fraction_observed 

    read_write_fast( out, fn=fn ) 

    return(out)

}

