  snowcrab_landings_db = function ( threshold.year = 0 ) {
    # glue historical data with marfis data

    recent = logbook.db( DS="logbook" ) # modern data:  mass in kg  -- must use all data and not postionally filtered data to get accurate totals
    # message( "Note:: Fishing 'yr' for CFA 4X has been set to starting year:: 2001-2002 -> 2001, etc.")
    recent = recent[recent$yr >= threshold.year, ]
    recent$subarea = recent$cfa0
    recent$cfa0 = NULL

    historical = snowcrab_historical_data()
    
    historical = historical[,names(recent)]
    all = rbind(historical, recent)

    all = all[all$yr > threshold.year, ]
    return(all)
  }
