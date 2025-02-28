

figure.fisheries.timeseries = function( 
  outdir=NULL,  
  region =  "cfanorth" , 
  region_label = "N-ENS" )  {

  dir.create( outdir, recursive=T, showWarnings=F  )

  fnbase = paste( "fisheries_timeseries_single", region, sep="_" )

  types = c("effort", "landings", "cpue" )
  fn = file.path( outdir, paste( fnbase, "_", types, ".png", sep="" ) ) 

  require(ggplot2)

  color_map = c("#E69F00", "#56B4E9",  "#CC79A7" , "#D55E00", "#F0E442")[1]
  shapes = c(15, 17, 19, 21, 23)[1]
  
  k = NULL
  k = fishery_data( toget="summary_annual"  )
  i = which( k$region==region )
  k = k[i,]

  k$region = region_label
  k$region = factor(k$region, levels=region_label)

  k$effort = k$effort / 1000
  k$landings = k$landings / 1000
    
  oE = ggplot(k, aes(x=yr, y=effort, fill=region, colour=region)) +
    geom_line( alpha=0.9, linewidth=1 ) +
    geom_point(aes(shape=region), size=5, alpha=0.7 )+
    labs(x="Year / Année", y="Effort (1000 trap haul) /\n Débarquements (1000 casiers levé)", size = rel(1.5)) +
    scale_colour_manual(values=color_map) +
    scale_fill_manual(values=color_map) +
    scale_shape_manual(values = shapes) +
    theme_light( base_size = 22) + 
    theme( legend.position="inside", legend.position.inside=c(0.25, 0.9), legend.title=element_blank()) 

  ggsave(filename=fn[1], plot=oE, device="png", width=12, height = 8)


  oL = ggplot(k, aes(x=yr, y=landings, fill=region, colour=region)) +
    geom_line( alpha=0.9, linewidth=1.2 ) +
    geom_point(aes(shape=region), size=5, alpha=0.7 )+
    labs(x="Year / Année", y="Landings (t) / Débarquements (t)", size = rel(1.5)) +
    theme_light( base_size = 22) + 
    # color_map = c("#E69F00", "#56B4E9",  "#CC79A7" )
    scale_colour_manual(values=color_map) +
    scale_fill_manual(values=color_map) +
    scale_shape_manual(values = shapes) +
    theme( legend.position="inside", legend.position.inside=c(0.85, 0.9), legend.title=element_blank()) 

  ggsave(filename=fn[2], plot=oL, device="png", width=12, height = 8)


  oC = ggplot(k, aes(x=yr, y=cpue, fill=region, colour=region)) +
    geom_line( alpha=0.9, linewidth=1 ) +
    geom_point(aes(shape=region), size=5, alpha=0.7 )+
    labs(x="Year / Année", y="Catch rate (kg/trap) /\n Taux de prise (kg/casier levé)", size = rel(1.5)) +
    # color_map = c("#E69F00", "#56B4E9",  "#CC79A7" )
    scale_colour_manual(values=color_map) +
    scale_fill_manual(values=color_map) +
    scale_shape_manual(values =shapes) +
    theme_light( base_size = 22) + 
    theme( legend.position="inside", legend.position.inside=c(0.1, 0.85), legend.title=element_blank()) 

  ggsave(filename=fn[3], plot=oC, device="png", width=12, height = 8)


  return( fn )

}

