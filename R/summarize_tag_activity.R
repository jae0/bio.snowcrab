summarize_tag_activity <- function(df) {
  # Create a local data.table so we don't accidentally modify the original input
  dt <- as.data.table(df)
  
  # Ensure chronological order per tag
  setorder(dt, tagid, timestamp)
  
  # Convert biological features to numeric, suppressing warnings for "na" strings
  suppressWarnings({
    dt[, cw_num := as.numeric(cw)]
    dt[, cc_num := as.numeric(cc)]
  })
  
  # Calculate summary statistics per tagid
  summary_dt <- dt[, .(
    
    # 1. Duration of activity in days
    duration_days = as.numeric(difftime(max(timestamp, na.rm = TRUE), 
                                        min(timestamp, na.rm = TRUE), 
                                        units = "days")),
    
    # 2. Biological changes (Last non-NA value minus First non-NA value)
    cw_change = if(length(na.omit(cw_num)) > 1) {
      last(na.omit(cw_num)) - first(na.omit(cw_num))
    } else { NA_real_ },
    
    cc_change = if(length(na.omit(cc_num)) > 1) {
      last(na.omit(cc_num)) - first(na.omit(cc_num))
    } else { NA_real_ },
    
    # 3. Total distance travelled (in meters)
    total_dist_m = {
      # Grab coordinates and drop any NA locations
      coords <- cbind(lon, lat)
      coords <- coords[complete.cases(coords), , drop = FALSE]
      
      if (nrow(coords) > 1) {
        # Create a spatial path (LINESTRING) and compute its length
        # Assuming input coordinates are standard lon/lat WGS84 (EPSG 4326)
        path <- sf::st_sfc(sf::st_linestring(coords), crs = 4326)
        as.numeric(sf::st_length(path))
      } else {
        0 # Distance is 0 if there is only 1 point or none
      }
    },
    
    # 4. (Optional) Helpful metadata
    n_points = .N
    
  ), by = tagid]
  
  return(summary_dt)
}