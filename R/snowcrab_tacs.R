snowcrab_tacs = function() {

    tacs = snowcrab_historical_data()[, c("yr", "cfa", "nlicences", "tac.tons") ]
    tacs[tacs$yr < 2005, ]
    # this needs to go into a database .. TODO
    # copied by had from old docs
    tacs = rbind( tacs,
        c(2005, "cfa4x", 9, 337.6),
        c(2006, "cfa4x", 9, 337.6),
        c(2007, "cfa4x", 9, 230), 
        c(2008, "cfa4x", 9, 230),
        c(2009, "cfa4x", 9, 230), 
        c(2010, "cfa4x", 9, 346),
        c(2011, "cfa4x", 9, 346),
        c(2012, "cfa4x", 9, 263),
        c(2013, "cfa4x", 9, 80), 
        c(2014, "cfa4x", 9, 80), 
        c(2015, "cfa4x", 9, 150),
        c(2016, "cfa4x", 9, 80),
        c(2017, "cfa4x", 9, 110),
        c(2018, "cfa4x", 9, 0),
        c(2019, "cfa4x", 9, 55),
        c(2020, "cfa4x", 9, NA),  ## lookup
        c(2021, "cfa4x", 9, NA),  ## lookup

        c(2005, "cfanorth", 78, 566),
        c(2006, "cfanorth", 78, 487),
        c(2007, "cfanorth", 78, 244), 
        c(2008, "cfanorth", 78, 244),
        c(2009, "cfanorth", 78, 576), 
        c(2010, "cfanorth", 78, 576),
        c(2011, "cfanorth", 78, 534),
        c(2012, "cfanorth", 78, 603),
        c(2013, "cfanorth", 78, 783), 
        c(2014, "cfanorth", 78, 783), 
        c(2015, "cfanorth", 78, 620),
        c(2016, "cfanorth", 78, 286),
        c(2017, "cfanorth", 78, 825),
        c(2018, "cfanorth", 78, 786),
        c(2019, "cfanorth", 78, 631),
        c(2020, "cfanorth", 78, NA), ## lookup
        c(2021, "cfanorth", 78, NA), ## lookup

        c(2005, "cfasouth", 114, 6353),
        c(2006, "cfasouth", 114, 4510),
        c(2007, "cfasouth", 115, 4950), 
        c(2008, "cfasouth", 115, 8316),
        c(2009, "cfasouth", 116, 10800), 
        c(2010, "cfasouth", 116, 13200),
        c(2011, "cfasouth", 116, 12120),
        c(2012, "cfasouth", 116, 11707),
        c(2013, "cfasouth", 116, 11311), 
        c(2014, "cfasouth", 116, 11311), 
        c(2015, "cfasouth", 116, 11311),
        c(2016, "cfasouth", 116, 9614),
        c(2017, "cfasouth", 116, 6730),
        c(2018, "cfasouth", 116, 6057),
        c(2019, "cfasouth", 116, 6663),
        c(2020, "cfasouth", 116, NA),  ## lookup
        c(2021, "cfasouth", 116, NA)  ## lookup

    )

    names(tacs) =c("yr", "region", "Licenses", "TAC") 
    tacs$yr = as.numeric(tacs$yr)
    tacs$Licenses = as.numeric(tacs$Licenses)
    tacs$TAC = as.numeric(tacs$TAC)

    return(tacs)
}