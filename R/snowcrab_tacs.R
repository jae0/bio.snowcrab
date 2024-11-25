snowcrab_tacs = function() {
if(!"aegis" %in% .packages())require("aegis")

    tac.db = CA.getTable("TAC")
    
    names(tac.db) = tolower(names(tac.db))
    tac.db$area[which(grepl(23, tac.db$area) | grepl(24, tac.db$area) | grepl(23, tac.db$area))] = "cfasouth"
    tac.db$area[which(grepl("4X", tac.db$area))] = "cfa4x"
    tac.db$area[which(grepl("N-ENS", tac.db$area))] = "cfanorth"
    
    setDT(tac.db)
    tac.db = tac.db[, .(tac=sum( as.numeric(tac), na.rm=TRUE)), by=.(area, yr) ]
    tac.db$yr = as.numeric( tac.db$yr )
  
    lic.hist = CA.getTable("NUMLICENSE")
    names(lic.hist) = c('yr', 'area', 'count')
    lic.hist$area = str_replace_all(lic.hist$area, " ", "")
   
     lic.db = CA.getTable("LIC_SUMMARY")
    names(lic.db) = c("count", "area", "yr")
    lic.db$area[which(grepl(23, lic.db$area) | grepl(24, lic.db$area) | grepl(23, lic.db$area))] = "cfasouth"
    lic.db$area[which(grepl("4X", lic.db$area))] = "cfa4x"
    lic.db$area[which(grepl("N-ENS", lic.db$area))] = "cfanorth"
    
    setDT(lic.db)
    lic.db$yr = as.numeric(lic.db$yr  )
    lic.db$count = as.numeric(lic.db$count  )
    lic.db$area = gsub(" ", "", lic.db$area  )
    lic.db = lic.db[, .(count=sum( as.numeric(count), na.rm=TRUE)), by=.(area, yr) ]
    lic.new = rbind(lic.db, lic.hist[which(lic.hist$yr < min(lic.db$yr)),])

    tacs.new = merge(as.data.frame(lic.new), as.data.frame(tac.db), by=c('area', 'yr'))
    
    tacs.new = tacs.new[, c("area", "yr", "count", "tac")] 
   
    tacs.new$yr = as.character(tacs.new$yr)

    names(tacs.new) =c("region", "yr", "Licenses", "TAC") 
    tacs.new$yr = as.numeric(tacs.new$yr)
    tacs.new$Licenses = as.numeric(tacs.new$Licenses)
    tacs.new$TAC = as.numeric(tacs.new$TAC)

    message( "Much has been moved to oracle tables. TAC's still need to be manually written via CA.insertTAC(cfa23.tac=????, cfa24.tac=????, Millbrook.tac = ????, nens.tac= ????, xxxx.tac= ????, year=????)
 Licenses numbers depend on CDD values being entered. Currently via CA.insertCDDcsv() but looking for cleaner solution to directly refence source tables. If 14.TAC_CA_OPPS is current then this should return correctly ")

    return(tacs.new)
}
