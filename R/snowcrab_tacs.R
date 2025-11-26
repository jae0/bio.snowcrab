

snowcrab_tacs = function(vn="region") {
  # vn can be: "region" or "subareas" or "area", currently only supports region
  if (!"aegis" %in% .packages()) require("aegis")
  
    if(vn == "region"){
      tac2.db = CA.getTable("TAC_REGION")
      tac2.db$REGION[which(tac2.db$REGION=="4X")] = "cfa4x"
      tac2.db$REGION[which(tac2.db$REGION=="SENS")] = "cfasouth"
      tac2.db$REGION[which(tac2.db$REGION=="NENS")] = "cfanorth"
      names(tac2.db) = c("area", "yr", "tac", "Licenses")
    }
    if(vn == "subarea"){
      tac2.db = CA.getTable("TAC_SUBAREA")
      tac2.db$SUBAREA[which(tac2.db$SUBAREA=="4X")] = "cfa4x"
      tac2.db$SUBAREA[which(tac2.db$SUBAREA=="4XE")] = "cfa4x"
      tac2.db$SUBAREA[which(tac2.db$SUBAREA=="4XW")] = "cfa4xW"
      tac2.db$SUBAREA[which(tac2.db$SUBAREA=="CFA23")] = "cfa23"
      tac2.db$SUBAREA[which(tac2.db$SUBAREA=="CFA24")] = "cfa24"
      tac2.db$SUBAREA[which(tac2.db$SUBAREA=="NENS")] = "cfanorth"
      names(tac2.db) = c("area", "yr", "tac", "Licenses")
    }
    
    tacs.new = tac2.db[, c("yr", "area", "Licenses", "tac")] 
    names(tacs.new) =c("yr", vn, "Licenses", "TAC") 
    tacs.new$Licenses = as.numeric(tacs.new$Licenses)
    tacs.new$TAC = as.numeric(tacs.new$TAC)
   
    message( "Much has been moved to oracle tables. TAC's now written via a call to web_fisheriesdata_update(). Licenses numbers are now assumed standard but can be canged in the input variables.")
    
    return(tacs.new)
  }
  

