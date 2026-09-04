plot_tag_trajectories = function(df, plot_crs=st_crs( projection_proj4string("utm20N")  )) {
  # Ensure input is a data.table
  setDT(df)

  # Sort the data chronologically to ensure geom_path connects points in the correct order
  setorder(df, tagid, timestamp)

  # Create a clearer categorical variable for the legend based on the 'tag' flag
  df[, event_type := fifelse(tag == 0, "Initial Mark", "Recapture")]
 
  corners = data.frame(lon=c(-66.1, -56.6), lat=c(42.8,47.4))
   
  bb = point_to_bbox( corners, plot_crs=plot_crs )

  additional_features = read_write_fast("c:/home/jae/projects/bstm/docs/movement/data/snowcrab_mapping_features_ggplot.rdz")
  
  df = df[which(is.finite(lon+lat)),]

  df = st_as_sf( df, coords= c("lon", "lat") )
  st_crs(df) =  st_crs( projection_proj4string("lonlat_wgs84") )

  df = st_transform(df, plot_crs )  # redundant .. in case input data is another projection
  # --- ADD THESE 3 LINES ---
  coords <- sf::st_coordinates(df)
  df$lon <- coords[, "X"]
  df$lat <- coords[, "Y"]

  # setDF(df)
p = ggplot(data = df) +
    
    # Map x and y locally here
    geom_path(aes(x = lon, y = lat, group = as.factor(tagid), color = as.factor(tagid)), 
              alpha = 0.6, linewidth = 0.8, show.legend = FALSE) +

    # Map x and y locally here as well
    geom_point(aes(x = lon, y = lat, shape = event_type, color = as.factor(tagid)), 
               size = 2.5) +

    scale_shape_manual(values = c("Initial Mark" = 17, "Recapture" = 16)) +

    additional_features +
    
    coord_sf(xlim =bb$x, ylim =bb$y, expand = FALSE, crs=plot_crs ) +  #

    theme_minimal() +
    labs(
      title = "Animal Movement Trajectories",
      x = "Longitude (UTM)", # Note: since you projected to UTM, these are no longer degrees
      y = "Latitude (UTM)",
      shape = "Event Type"
    ) +
    guides(color = "none")

  return(p) 
}
