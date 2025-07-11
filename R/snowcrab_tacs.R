snowcrab_tacs = function( vn="region" ) {
    
    # vn can be: "region" or "subareas" or "area"

    if (!"aegis" %in% .packages()) require("aegis")

    CA = read_write_fast( file.path( project.datadirectory("bio.snowcrab", "data", "CA", "CA_db"), "CA.rdz" ) )
    
    tac.db = CA$tac.db
    names(tac.db) = tolower(names(tac.db))
    
    tac.db$yr = as.numeric( tac.db$yr )

    tac.db$subarea = tac.db$area 
    tac.db$subarea = tolower(tac.db$subarea) 
    tac.db$subarea = gsub(" ", "", tac.db$subarea)
    tac.db$subarea = gsub("4x", "cfa4x", tac.db$subarea)
    tac.db$subarea = gsub("millbrook", "", tac.db$subarea) 
    tac.db$subarea = gsub("n-ens", "cfanorth", tac.db$subarea) 
    
    tac.db$region = NA
    tac.db$region[which(grepl(23, tac.db$area) | grepl(24, tac.db$area) | grepl(23, tac.db$area))] = "cfasouth"
    tac.db$region[which(grepl("4X", tac.db$area))] = "cfa4x"
    tac.db$region[which(grepl("N-ENS", tac.db$area))] = "cfanorth"
    
    setDT(tac.db)
    tac.db = tac.db[, .(tac=sum( as.numeric(tac), na.rm=TRUE)), by=c(vn, "yr") ]
  
    lic.hist = CA$lic.hist
    names(lic.hist) = c('yr', 'area', 'count')
    lic.hist$area = gsub(" ", "", lic.hist$area)
    lic.hist$area = gsub(" ", "", lic.hist$area)
    lic.hist$region = lic.hist$area 
    lic.hist$subarea = lic.hist$area  # there is no discrimination for 23 vs 24 in this .. needs to be added  <<< fixme
     
    lic.db = CA$lic.db
    names(lic.db) = c("count", "area", "yr")
    lic.db$yr = as.numeric(lic.db$yr  )
    lic.db$count = as.numeric(lic.db$count  ) 
    lic.db$region = lic.db$area
    
    lic.db$region = NA
    lic.db$region[which(grepl(23, lic.db$area) | grepl(24, lic.db$area) | grepl(23, lic.db$area))] = "cfasouth"
    lic.db$region[which(grepl("4X", lic.db$area))] = "cfa4x"
    lic.db$region[which(grepl("N-ENS", lic.db$area))] = "cfanorth"
    
    lic.db$subarea = lic.db$area 
    lic.db$subarea = tolower(lic.db$subarea) 
    lic.db$subarea = gsub(" ", "", lic.db$subarea)
    lic.db$subarea = gsub("4x", "cfa4x", lic.db$subarea)
    lic.db$subarea = gsub("millbrook", "", lic.db$subarea) 
    lic.db$subarea = gsub("n-ens", "cfanorth", lic.db$subarea) 
 
    setDT(lic.db)
    lic.db = lic.db[, .(count=sum( as.numeric(count), na.rm=TRUE)), by=c(vn, "yr") ]
    lic.hist = lic.hist[, names(lic.db) ]
    lic.new = rbind(lic.db, lic.hist[which(lic.hist$yr < min(lic.db$yr)), ])

    tacs.new = merge(as.data.frame(lic.new), as.data.frame(tac.db), by=c(vn, 'yr'))
    
    tacs.new = tacs.new[, c(vn, "yr", "count", "tac")] 
   
    tacs.new$yr = as.character(tacs.new$yr)

    names(tacs.new) =c(vn, "yr", "Licenses", "TAC") 
    tacs.new$yr = as.numeric(tacs.new$yr)
    tacs.new$Licenses = as.numeric(tacs.new$Licenses)
    tacs.new$TAC = as.numeric(tacs.new$TAC)

    message( "Much has been moved to oracle tables. TAC's still need to be manually written via CA.insertTAC(cfa23.tac=????, cfa24.tac=????, Millbrook.tac = ????, nens.tac= ????, xxxx.tac= ????, year=????)
 Licenses numbers depend on CDD values being entered. Currently via CA.insertCDDcsv() but looking for cleaner solution to directly reference source tables. If 14.TAC_CA_OPPS is current then this should return correctly ")

    return(tacs.new)
}
