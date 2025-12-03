fishery_data_figures_seasonal = function( 
    yrsplot = NULL,   
    FD=NULL, 
    toget=c("cummulative_landings", "cpue" ),
    mau = "region",
    outdir=NULL,
    time_resolution = "summary_weekly"
  ) {
     
  dir.create( outdir, recursive=TRUE, showWarnings=FALSE )

  fnbase = "fisheries_seasonal"

  require(ggplot2)

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
    i4x = which(M[[mau]]=="cfa4x" & M$week < 34 )
    M$week[i4x] = 52 + M$week[i4x]
    M$yr[i4x] = M$yr[i4x] + 1  # revert fishing year to calendar year

    M$julian = M$week / 52
    M$timestamp = M$yr + M$julian
    M$timestamp[i4x] = M$timestamp[i4x] - 1
  }

  setnames(M, mau, "internal")

  reglkup = data.table(internal = maus[["internal"]], area = maus[["labels"]]) 
  M = reglkup[M, on="internal"] 

  M$area = factor(M$area, levels=maus[["labels"]]  )
  M$yr = factor(M$yr)
  M = M[ order( area,  yr, week, decreasing=FALSE), ]

  out = list(FD=FD, toget=toget, CALL=match.call())

  if ("cummulative_landings" %in% toget) {
    M[ is.na(landings), "landings"] = 0 

    M[, landings_cummulative := cumsum(landings), by=.(yr, area)]
    M[, landings_cummulative := landings_cummulative/max(landings_cummulative, na.rm=TRUE), by=.(yr, area)]

    out[["cummulative_landings"]] = ggplot( M, aes(julian, landings_cummulative, shape=area, colour=yr, group=interaction( yr, area)) ) + 
        geom_point( size=2.5 ) + 
        geom_line() + 
        xlim( 0 , 1.3) +
        xlab("Year fraction") + 
        ylab("Cummulative landings") + 
        theme_bw() + 
        theme(legend.title = element_blank(), legend.position = "bottom", legend.key = element_rect(colour = NA, fill = NA))

        fn = file.path( outdir, paste( fnbase, "_", "cummulative_landings", ".png", sep="" ) )

        ggsave(filename=fn, plot=out[["cummulative_landings"]], device="png", width=12, height = 8)
  }

  if ("cpue" %in% toget) {

    out[["cpue"]] = ggplot( M, aes(timestamp, cpue, colour=area, group=interaction(yr, area)) ) + 
        geom_point( shape="circle", size=2.5 ) + 
        geom_line() + 
        xlim( min(yrsplot)-0.025 , max(yrsplot)+ 1.025) +
        xlab("Year") + 
        ylab("Catch Rate (kg / trap)") + 
        theme_bw() + 
        theme(legend.title = element_blank(), legend.position = "bottom", legend.key = element_rect(colour = NA, fill = NA))
    
    fn = file.path( outdir, paste( fnbase, "_", "cpue", ".png", sep="" ) )

    ggsave(filename=fn, plot=out[["cpue"]], device="png", width=12, height = 8)

  }
 
  if ("cummulative_effort" %in% toget) {
    M[ is.na(effort), "effort"] = 0 

    M[, effort_cummulative := cumsum(effort), by=.(yr, area)]
    M[, effort_cummulative := effort_cummulative/max(effort_cummulative, na.rm=TRUE), by=.(yr, area)]

    out[["cummulative_effort"]] = ggplot( M, aes(julian, effort_cummulative, shape=area, colour=yr, group=interaction( yr, area)) ) + 
        geom_point( size=2.5 ) + 
        geom_line() + 
        xlim( 0 , 1.3) +
        xlab("Year fraction") + 
        ylab("Cummulative effort") + 
        theme_bw() + 
        theme(legend.title = element_blank(), legend.position = "bottom", legend.key = element_rect(colour = NA, fill = NA))

        fn = file.path( outdir, paste( fnbase, "_", "cummulative_effort", ".png", sep="" ) )

        ggsave(filename=fn, plot=out[["cummulative_effort"]], device="png", width=12, height = 8)
  }

  return(outdir)
}      
