  snowcrab_landings_db = function ( ) {
    # glue historical data with marfis data

    threshold.year = 2004

    b = logbook.db( DS="logbook" ) # modern data:  mass in kg  -- must use all data and not postionally filtered data to get accurate totals
    # message( "Note:: Fishing 'yr' for CFA 4X has been set to starting year:: 2001-2002 -> 2001, etc.")
    b = b[b$yr >= threshold.year, ]

    a = snowcrab_historical_data()
    a = a[a$yr < threshold.year,]
    d = a[,names(b)]
    c = rbind(d, b)

    return(c)
  }
