
 
map.observer.entanglements = function(p, basedir=tempdir(), years=p$yrs_observer, region="cfaall", 
  plot_crs=st_crs("EPSG:32620"),  # UTM20N ) {
  
  require(ggplot2)

  oss = observer.db( DS="bycatch_clean_data", p=p,  yrs=years )  # Prepare at sea observed data
  i = polygon_inside( oss[,  c("lon", "lat")], region=region )
  oss =  oss[i,]

  whales = oss[ grep("whale", common, ignore.case=TRUE), which=TRUE ]
  leatherback = oss[ grep("leatherback", common, ignore.case=TRUE), which=TRUE ]
  basking_shark = oss[ grep("basking shark",  common, ignore.case=TRUE), which=TRUE ]
  
  W = oss[ whales, .(whales=.N), by=.(yr)] 
  L = oss[ leatherback, .(leatherback=.N), by=.(yr)] 
  B = oss[ basking_shark, .(basking_shark=.N), by=.(yr)] 

  out = data.table( yr=years )
  out = W[out, on="yr"]
  out = L[out, on="yr"]
  out = B[out, on="yr"]
  out[ is.na(out) ] = 0

  colnames(out) = c("Year", "Whale", "Leatherback turtle", "Basking shark")

  print("Entaglements:")
  print( out )

  additional_features = snowcrab_mapping_features(p, redo=FALSE ) [["ggplot"]] [["layers"]]
  
  title = "Entanglements"

  vn_label ="fish_yr"

  toplot = st_as_sf( oss[, c("lon", "lat", "trip", "fishyr") ], coords= c("lon", "lat"), crs=st_crs( projection_proj4string("lonlat_wgs84") ) )
  toplot = st_transform(toplot, plot_crs )

  toplot$species = NA
  toplot$species[whales] = "Whale"
  toplot$species[leatherback] = "Leatherback turtle"
  toplot$species[basking_shark] = "Basking shark"
 
  xr = p$corners$plon
  yr = p$corners$plat

  plt = ggplot() +
    geom_sf( data = toplot[which(is.na(toplot$species)),], alpha=0.6, lwd=0, cex=0.75, col="grey" ) +  
    geom_sf( data=toplot[which(!is.na(toplot$species)),], aes(colour=species ), alpha=0.9, size=5 ) +
    additional_features  +
    coord_sf(xlim =xr, ylim =yr, expand = FALSE, crs=plot_crs) +
    theme(
      axis.line=element_blank(),
      axis.ticks=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(), 
      legend.position = "bottom",
      #legend.title = element_blank(),
      panel.background = element_rect(fill =NA),
      panel.border=element_blank(),
      panel.grid.major = element_line(color = "grey"),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(), 
      plot.caption = element_text(hjust = 0, size = 14)
  )
  
  print(plt)
 
  if (!file.exists(basedir)) dir.create( basedir,  showWarnings = FALSE, recursive = TRUE )

  fn = file.path( basedir, "observed_bycatch_entanglements.png" )

  ggsave(filename=fn, plot=plt,  width=12, height = 8)

  return(fn)
}