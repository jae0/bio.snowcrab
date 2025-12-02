fishery_data_figures_subannual = function( 
    yrsplot = NULL,   
    FD=NULL, 
    toget="cummulative_landings",
    mau = "region",
    time_resolution = "summary_weekly"
  ) {
     

  
  maus = management_areal_units( mau=mau )  
  
  if (is.null(FD))  FD = fishery_data( mau=mau )  # mass in tonnes
  
  if (time_resolution=="summary_annual") {
    message( "These are seasonal figures. An annual resolution makes no sense. ")
    message( "Returning just the data. ")
    return(FD[[time_resolution]])     
  }

  M = FD[[time_resolution]]

  if (is.null(yrsplot))  yrsplot = max(M$yr) + c(-2:0)
  M = M[ yr %in% yrsplot,]

  # year offset to correct the sequence for 4X
  if (exists("week", M)) {
    i4x = which(M$region=="cfa4x" & M$week < 34 )
    M$week[i4x] = 52 + M$week[i4x]
    M$yr[i4x] = M$yr[i4x] + 1  # revert fishing year to calendar year

    M$julian = M$week / 52
    M$timestamp = M$yr + M$julian
    M$timestamp[i4x] = M$timestamp[i4x] - 1
  }

  reglkup = data.table(region = maus[["internal"]], reg = maus[["labels"]]) 
  M = reglkup[M, on="region"] 

  M$reg = factor(M$reg, levels=maus[["labels"]]  )
  M$yr = factor(M$yr)
  M = M[ order( reg,  yr, week, deMeasing=FALSE), ]

  out = list(FD=FD, toget=toget, CALL=match.call())

  if (toget %in% "cummulative_landings") {
    M[ is.na(landings), "landings"] = 0 

    M[, landings_cummulative := cumsum(landings), by=.(yr, reg)]
    M[, landings_cummulative := landings_cummulative/max(landings_cummulative, na.rm=TRUE), by=.(yr, reg)]

    out[["cummulative_landings"]] = ggplot( M, aes(julian, landings_cummulative, shape=reg, colour=yr, group=interaction( yr, reg)) ) + 
        geom_point( size=2.5 ) + 
        geom_line() + 
        xlim( 0 , 1.3) +
        xlab("Year fraction") + 
        ylab("Cummulative catch") + 
        theme_bw() + 
        theme(legend.title = element_blank(), legend.position = "bottom", legend.key = element_rect(colour = NA, fill = NA),)
  }

  if (toget %in% "cpue") {

    out[["cpue"]] = ggplot( M, aes(timestamp, cpue, colour=reg, , group=interaction(yr, reg)) ) + 
        geom_point( shape="circle", size=2.5 ) + 
        geom_line() + 
        xlim( min(yrsplot)-0.025 , max(yrsplot)+ 1.025) +
        xlab("Year") + 
        ylab("Catch Rate (kg / trap)") + 
        theme_bw() + 
        theme(legend.title = element_blank(), legend.position = "bottom", legend.key = element_rect(colour = NA, fill = NA),)
  }
 
  return(out)
}      
