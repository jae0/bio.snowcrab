filter_dead = function(dt, time_threshold_days = 30, dist_threshold_meters = 50){

    # 1. Sort chronologically by tag
    setorder(dt, tagid, timestamp)

    # 2. Identify the final known location and time for each tag
    dt[, `:=`(
      last_lon = last(lon),
      last_lat = last(lat),
      last_time = last(timestamp)
    ), by = tagid]

    # 3. Calculate spatial distance (meters) and time difference (days) from each ping to the final ping
    dt[, dist_to_last := distGeo(
      matrix(c(lon, lat), ncol = 2),
      matrix(c(last_lon, last_lat), ncol = 2)
    )]
    dt[, days_to_last := as.numeric(difftime(last_time, timestamp, units = "days"))]

    # 4. Sort reverse-chronologically to build the terminal cluster
    setorder(dt, tagid, -timestamp)

    # cummax() calculates the maximum distance the tag will travel from this point until its battery dies
    dt[, max_dist_to_end := cummax(dist_to_last), by = tagid]

    # 5. Re-sort chronologically
    setorder(dt, tagid, timestamp)

    # 6. Flag records that belong to a terminal stationary period exceeding the time threshold
    dt[, is_dead := {

      # True if the animal never leaves the distance threshold from this point forward
      in_terminal_cluster <- max_dist_to_end <= dist_threshold_meters

      # Calculate the total duration of this continuous non-moving period
      cluster_duration <- if (any(in_terminal_cluster, na.rm = TRUE)) {
        max(days_to_last[which(in_terminal_cluster)])
      } else {
        0
      }

      # Flag as dead if it is in the terminal cluster AND the cluster duration is long enough
      in_terminal_cluster & (cluster_duration >= time_threshold_days)

    }, by = tagid]

    # Clean up intermediate calculation columns to keep the dataframe tidy
    dt[, c("last_lon", "last_lat", "last_time", "dist_to_last", "days_to_last", "max_dist_to_end") := NULL]

  return(dt)
}
