

figure.fisheries.timeseries = function( 
    outdir=NULL,  
    mau = "region",
    region_id = NULL
  )  {

  dir.create( outdir, recursive=T, showWarnings=F  )

  fnbase = paste( "fisheries_timeseries_single", region_id, sep="_" )

  types = c("effort", "landings", "cpue" )
  fn = file.path( outdir, paste( fnbase, "_", types, ".png", sep="" ) ) 

  require(ggplot2)

  maus = management_areal_units( mau=mau )  
  j = which(maus[["internal"]] == region_id) 

  color_map = maus[["color_map"]][1]
  shapes = maus[["shapes"]][1]
  
  FD = NULL
  FD = fishery_data( mau=mau )
  AN = FD[["summary_annual"]]

  # force use of "region" as variable name ...
  AN$region = factor(AN[[mau]], levels=maus[["labels"]] )

  AN$effort = AN$effort / 1000
  AN$landings = AN$landings / 1000
 
  i = which( AN$region==region_id )
  AN = AN[i,]

   
  oE = ggplot(AN, aes(x=yr, y=effort, fill=region, colour=region)) +
    geom_line( alpha=0.9, linewidth=1 ) +
    geom_point(aes(shape=region), size=5, alpha=0.7 )+
    labs(x="Year / Année", y="Effort (1000 trap haul) /\n Débarquements (1000 casiers levé)" ) +
    scale_colour_manual(values=color_map) +
    scale_fill_manual(values=color_map) +
    scale_shape_manual(values = shapes) +
    theme_light( base_size = 22) + 
    theme( legend.position="inside", legend.position.inside=c(0.25, 0.9), legend.title=element_blank()) 

  ggsave(filename=fn[1], plot=oE, device="png", width=12, height = 8)


  oL = ggplot(AN, aes(x=yr, y=landings, fill=region, colour=region)) +
    geom_line( alpha=0.9, linewidth=1.2 ) +
    geom_point(aes(shape=region), size=5, alpha=0.7 )+
    labs(x="Year / Année", y="Landings (t) / Débarquements (t)" ) +
    theme_light( base_size = 22) + 
    # color_map = c("#E69F00", "#56B4E9",  "#CC79A7" )
    scale_colour_manual(values=color_map) +
    scale_fill_manual(values=color_map) +
    scale_shape_manual(values = shapes) +
    theme( legend.position="inside", legend.position.inside=c(0.85, 0.9), legend.title=element_blank()) 

  ggsave(filename=fn[2], plot=oL, device="png", width=12, height = 8)


  oC = ggplot(AN, aes(x=yr, y=cpue, fill=region, colour=region)) +
    geom_line( alpha=0.9, linewidth=1 ) +
    geom_point(aes(shape=region), size=5, alpha=0.7 )+
    labs(x="Year / Année", y="Catch rate (kg/trap) /\n Taux de prise (kg/casier levé)" ) +
    # color_map = c("#E69F00", "#56B4E9",  "#CC79A7" )
    scale_colour_manual(values=color_map) +
    scale_fill_manual(values=color_map) +
    scale_shape_manual(values =shapes) +
    theme_light( base_size = 22) + 
    theme( legend.position="inside", legend.position.inside=c(0.1, 0.85), legend.title=element_blank()) 

  ggsave(filename=fn[3], plot=oC, device="png", width=12, height = 8)


  return( fn )

}

