#### BC - Need to modify script to have one function with if conditions based on a DS variable instead of
#### multiple function. 

##Must insert budget every year!!!!!!!!!!!
CA.insertVariables = function(price_per_mt = 1100, budget, year){
  CA.proj.path = file.path(project.datadirectory("bio.snowcrab","data", "CA","CA_db"))
  con <- dbConnect(RSQLite::SQLite(), file.path(CA.proj.path, "ca.db")) 
  
  statement = paste("INSERT OR REPLACE INTO VariableLookup(yr, price_per_mt, budget) VALUES(",year,",",price_per_mt,",", budget,");", sep = "")
  rs = dbSendQuery(con, statement)
  dbHasCompleted(rs)
  dbClearResult(rs)
  dbDisconnect(con, file.path(CA))
  
  }
##Must insert TAC values each year!!!!!!!!!!!
CA.insertTAC = function(cfa23.tac, cfa24.tac, Millbrook.tac = 250, nens.tac, xxxx.tac, year){
  CA.proj.path = file.path(project.datadirectory("bio.snowcrab","data", "CA","CA_db"))
  con <- dbConnect(RSQLite::SQLite(), file.path(CA.proj.path, "ca.db")) 
  
  statement = paste("SELECT * FROM TAC where area = 'CFA 23' AND yr = ",year, ";", sep = "")
  rs = dbSendQuery(con, statement)
  d1 <- dbFetch(rs)    
  dbHasCompleted(rs)
  dbClearResult(rs)
  if(nrow(d1) == 0){
    statement = paste('INSERT INTO TAC (yr, area, tac) VALUES(
                      "',year,'", "CFA 23", "', cfa23.tac,'");', sep = "")
    rs = dbSendQuery(con, statement)
    dbHasCompleted(rs)
    dbClearResult(rs)
  }
  statement = paste("SELECT * FROM TAC where area = 'CFA 24' AND yr = ",year, ";", sep = "")
  rs = dbSendQuery(con, statement)
  d1 <- dbFetch(rs)      # extract data in chunks of 10 rows
  dbHasCompleted(rs)
  dbClearResult(rs)
  if(nrow(d1) == 0){
    statement = paste('INSERT INTO TAC (yr, area, tac) VALUES(
                      "',year,'", "CFA 24", "', cfa24.tac,'");', sep = "")
    rs = dbSendQuery(con, statement)
    dbHasCompleted(rs)
    dbClearResult(rs)
  }
  
  statement = paste("SELECT * FROM TAC where area = 'CFA 24 Millbrook' AND yr = ",year, ";", sep = "")
  rs = dbSendQuery(con, statement)
  d1 <- dbFetch(rs)      # extract data in chunks of 10 rows
  dbHasCompleted(rs)
  dbClearResult(rs)
  if(nrow(d1) == 0){
    statement = paste('INSERT INTO TAC (yr, area, tac) VALUES(
                      "',year,'", "CFA 24 Millbrook", "', Millbrook.tac,'");', sep = "")
    rs = dbSendQuery(con, statement)
    dbHasCompleted(rs)
    dbClearResult(rs)
  }
  statement = paste("SELECT * FROM TAC where area = 'N-ENS' AND yr = ",year, ";", sep = "")
  rs = dbSendQuery(con, statement)
  d1 <- dbFetch(rs)      # extract data in chunks of 10 rows
  dbHasCompleted(rs)
  dbClearResult(rs)
  if(nrow(d1) == 0){
    statement = paste('INSERT INTO TAC (yr, area, tac) VALUES(
                      "',year,'", "N-ENS", "', nens.tac,'");', sep = "")
    rs = dbSendQuery(con, statement)
    dbHasCompleted(rs)
    dbClearResult(rs)
  }   
  
  statement = paste("SELECT * FROM TAC where area = '4X' AND yr = ",year, ";", sep = "")
  rs = dbSendQuery(con, statement)
  d1 <- dbFetch(rs)      # extract data in chunks of 10 rows
  dbHasCompleted(rs)
  dbClearResult(rs)
  if(nrow(d1) == 0){
    statement = paste('INSERT INTO TAC (yr, area, tac) VALUES(
                      "',year,'", "4X", "', xxxx.tac,'");', sep = "")
    rs = dbSendQuery(con, statement)
    dbHasCompleted(rs)
    dbClearResult(rs)
    
  }   
  
  dbDisconnect(con, file.path(CA))
  
}

## Must load yearly CDD data This data is received as excell sheet. Save each tab as csv to CDDcsv folderand call
## this ffunction
CA.insertCDDcsv = function() {
  CA.proj.path = file.path(project.datadirectory("bio.snowcrab","data", "CA","CA_db"))
  lof = list.files(file.path(CA.proj.path, "CDDcsv"), full.names=T)
  con <- dbConnect(RSQLite::SQLite(), file.path(CA.proj.path, "ca.db"))
  
  
  for(i in 1:length(lof)){
    if(grepl("all", basename(lof[i]))){
      da = read.csv(lof[i])
      da = da[,-1]
      year = unique(da$Quota.Year)
      da = da[,-2]
      names(da) = c("LICENCE", "AREA", "INITIALofTAC", "LICENCEHOLDER")
      dx = split(da, da$AREA)
      for(j in 1:length(dx)){
        dasub = as.data.frame(dx[j])
        names(dasub) = c("LICENCE", "AREA", "INITIALofTAC", "LICENCEHOLDER")
        dasub$CURRENTOFTAC = dasub$INITIALofTAC
        dasub$LICENCE_UOF = "UNKNOWN"
        print(paste(basename(lof[i]), unique(dasub$AREA), sep = "  "))
        print(paste(names(dasub)[3], " has sum of percents of: ", sum(as.numeric(as.character(dasub$INITIALofTAC))), sep=""))
        print("")
      
        dasub <- dasub[order(dasub$LICENCEHOLDER, dasub$LICENCE),]
        
        for(k in 1:nrow(dasub)){
        statement = paste("SELECT * FROM SnowCrabCDDpercents where licence = '",dasub$LICENCE[k],"' AND yr = ",year, ";", sep = "")
        
        rs = dbSendQuery(con, statement)
        
        
        d1 <- dbFetch(rs)      # extract data in chunks of 10 rows
        dbHasCompleted(rs)
        dbClearResult(rs)
        if(nrow(d1) == 0){
          statement = paste('INSERT INTO SnowCrabCDDpercents (licence, licence_holder, initial_percent_tac, current_percent_tac, licence_uof, initial_percent_tac_uof, current_percent_tac_uof, yr) VALUES(
                            "',as.character(dasub$LICENCE[k]),'","',  dasub$LICENCEHOLDER[k],'","',dasub$INITIALofTAC[k],'","',dasub$CURRENTOFTAC[k],'","',dasub$LICENCE_UOF[k],'","',dasub$INITIALofTAC[k],'","',dasub$CURRENTOFTAC[k],'","',year,'");', sep = "")
          
          rs = dbSendQuery(con, statement)
          dbHasCompleted(rs)
          dbClearResult(rs)
        }
        else{
          
        }
        
        
        }
        
        
      }
    }  
                     
    if(!grepl("uof", basename(lof[i])) && !grepl("all", basename(lof[i]))){
      da = read.csv(lof[i])
      
      da2 = read.csv(gsub("\\.csv", "uof.csv", lof[i]))

      names(da) = gsub("QTA", "TAC",names(da))
      da = da[names(da)[grepl("LICENCE", names(da)) | grepl("TAC", names(da))]]
      names(da) = gsub("\\.", "",names(da))
      
      ## 2017 had no current of TAC column, database structure built to capture this information
      ## incase it was important. For years without current of TAC will assume same as initial of TAC. 
      if(length(names(da)) == 3){
        names(da) = c("LICENCE", "LICENCEHOLDER", "INITIALofTAC")
        da$CURRENTOFTAC = da$INITIALofTAC
        
        
        
        
      }
      
      names(da) = c("LICENCE", "LICENCEHOLDER", "INITIALofTAC", "CURRENTOFTAC")
      
 
      
      ind = which(da[,1] == "")
      if(length(ind) > 0) da = da[-ind,]
      ind = which(da[,2] == "")
      if(length(ind) > 0) da = da[-ind,]
      
      
      names(da2) = gsub("QTA", "TAC",names(da2))
      da2 = da2[names(da2)[grepl("LICENCE", names(da2)) | grepl("TAC", names(da2))]]
      names(da2) = gsub("\\.", "",names(da2))
      
      ## 2017 had no current of TAC column, database structure built to capture this information
      ## incase it was important. For years without current of TAC will assume same as initial of TAC. 
      if(length(names(da2)) == 3){
        names(da2) = c("LICENCE", "LICENCEHOLDER", "INITIALofTAC")
        da2$CURRENTOFTAC = da2$INITIALofTAC
      }
      
      names(da2) = c("LICENCE", "LICENCEHOLDER", "INITIALofTAC", "CURRENTOFTAC")

      ind = which(da2[,1] == "")
      if(length(ind) > 0) da2 = da2[-ind,]
      ind = which(da2[,2] == "")
      if(length(ind) > 0) da2 = da2[-ind,]
      
      print(basename(lof[i]))
      print(paste(names(da)[3], " has sum of percents of: ", sum(as.numeric(da[,3])), sep=""))
      print(paste(names(da)[4], " has sum of percents of: ", sum(as.numeric(da[,4])), sep=""))
      print("")
      
      
      print(basename(gsub("\\.csv", "uof.csv", lof[i])))
      print(paste(names(da2)[3], " has sum of percents of: ", sum(as.numeric(da2[,3])), sep=""))
      print(paste(names(da2)[4], " has sum of percents of: ", sum(as.numeric(da2[,4])), sep=""))
      print("")
      
      
      da <- da[order(da$LICENCEHOLDER, da$LICENCE),]
      da2 <- da2[order(da2$LICENCEHOLDER, da2$LICENCE),]
      if(any(as.character(da2$LICENCEHOLDER) != as.character(da$LICENCEHOLDER))){
        print("There is a problem with the CDD data")
      }
      else{
        year = substr(basename(lof[i]), 1, 4)
        dax = cbind(da, da2$LICENCE, da2$INITIALofTAC, da2$CURRENTOFTAC)
        names(dax) = c("licence", "licence_holder",  "initial_percent_tac",  "current_percent_tac",	"licence_uof", "initial_percent_tac_uof",	"current_percent_tac_uof")
        
        for(j in 1:nrow(dax)){
          lic = da$LICENCE[j]
          if(grepl("ACRQ", lic)){
            lic = unlist(strsplit(as.character(lic), "\\("))[1]
          }
          statement = paste("SELECT * FROM SnowCrabCDDpercents where licence = '",lic,"' AND yr = ",year, ";", sep = "")
          
          rs = dbSendQuery(con, statement)
          
          
          d1 <- dbFetch(rs)      # extract data in chunks of 10 rows
          dbHasCompleted(rs)
          dbClearResult(rs)
          if(nrow(d1) == 0){
            statement = paste('INSERT INTO SnowCrabCDDpercents (licence, licence_holder, initial_percent_tac, current_percent_tac, licence_uof, initial_percent_tac_uof, current_percent_tac_uof, yr) VALUES(
                              "',as.character(lic),'","',dax$licence_holder[j],'","',dax$initial_percent_tac[j],'","',dax$current_percent_tac[j],'","',dax$licence_uof[j],'","',dax$initial_percent_tac_uof[j],'","',dax$current_percent_tac_uof[j],'","',year,'");', sep = "")
            
            rs = dbSendQuery(con, statement)
            dbHasCompleted(rs)
            dbClearResult(rs)
          }
          else{
        
          }
          
        }
      }
      
      
    }
  }
  
  dbDisconnect(con, file.path(CA))
  res = CA.checkMatch()
  if(length(res$a) != 0){
    warning("The above licence holders have no matching affilation information in the affilitation table:")
  print(res$a)
    }
  if(length(res$b) != 0){
    warning("The above affilations have no entries in the percents table:")
 print(res$b)
  }
  }

CA.writeIndividualContributions = function(year){
  CA.proj.path = file.path(project.datadirectory("bio.snowcrab","data", "CA","CA_db"))
  
  d1 = CA.getTable("CalcIndividuals", year)
  
  if(!dir.exists(file.path(CA.proj.path,"reports", year)))
    dir.create(file.path(CA.proj.path,"reports", year))
  write.csv(d1, file.path(CA.proj.path, "reports",year, "individuals.csv"))
  
  names(d1) = gsub(" ", "", names(d1))
  print("INDIVIDUAL CONTRIBUTIONS CHECK")
  print(paste("Sum of percents: ", sum(d1$percent_total), sep = ""))
  print(paste("Total CA contributions: ", sum(d1$ca_contribution), sep = ""))
  print(paste("Total Use of Fish Required: ", sum(d1$uof_required_mt), sep = ""))
  
}

CA.writePartnerContributions = function(year){
  CA.proj.path = file.path(project.datadirectory("bio.snowcrab","data", "CA","CA_db"))
  d1 = CA.getTable("CalcPartner", year)
  
  if(!dir.exists(file.path(CA.proj.path,"reports", year)))
    dir.create(file.path(CA.proj.path,"reports", year))
  write.csv(d1, file.path(CA.proj.path, "reports",year, "partner.csv"))
  
  names(d1) = gsub(" ", "", names(d1))
  print("PARTNER CONTRIBUTIONS CHECK")
  print(paste("Sum of percents: ", sum(d1$sum_of_percent), sep = ""))
  print(paste("Total CA contributions: ", sum(d1$sum_ca_contribution), sep = ""))
  print(paste("Total Use of Fish Required: ", sum(d1$sum_uof_required_mt), sep = ""))
  print("")
}

CA.writeAreaContributions = function(year){
  CA.proj.path = file.path(project.datadirectory("bio.snowcrab","data", "CA","CA_db"))
  d1 = CA.getTable("CalcArea", year)
  
  if(!dir.exists(file.path(CA.proj.path,"reports", year)))
    dir.create(file.path(CA.proj.path,"reports", year))
  write.csv(d1, file.path(CA.proj.path, "reports", year, "area.csv"))
  
  names(d1) = gsub(" ", "", names(d1))
  print("AREA CONTRIBUTIONS CHECK")
  print(paste("Sum of percents: ", sum(d1$sum_of_total_tac_percents), sep = ""))
  print(paste("Total CA contributions: ", sum(d1$sum_of_ca_contributions), sep = ""))
  print(paste("Total Use of Fish Required: ", sum(d1$sum_of_uof_required_mt), sep = ""))
  print("")
}


CA.removeTAC = function(year){
  CA.proj.path = file.path(project.datadirectory("bio.snowcrab","data", "CA","CA_db"))
  con <- dbConnect(RSQLite::SQLite(), file.path(CA.proj.path, "ca.db")) 
  
  statement = paste("DELETE FROM TAC where yr = ",year, ";", sep = "")
  rs = dbSendQuery(con, statement)
  dbHasCompleted(rs)
  dbClearResult(rs)
  dbDisconnect(con, file.path(CA))
  
}

CA.getTable = function(tablename, sub.year = NULL){ 
  CA.proj.path = file.path(project.datadirectory("bio.snowcrab","data", "CA","CA_db"))
  con <- dbConnect(RSQLite::SQLite(), file.path(CA.proj.path, "ca.db")) 
  if(is.null(sub.year))
    statement = paste("SELECT * FROM ", tablename, ";", sep = "")
  else
    statement = paste("SELECT * FROM ", tablename, " where yr = '", sub.year, "';", sep = "")
  
  rs = dbSendQuery(con, statement)
  d1 <- dbFetch(rs)      # extract data in chunks of 10 rows
  dbHasCompleted(rs)
  dbClearResult(rs)
  dbDisconnect(con, file.path(CA))
  return(d1)
}

CA.checkMatch = function(){ 
  CA.proj.path = file.path(project.datadirectory("bio.snowcrab","data", "CA","CA_db"))
  con <- dbConnect(RSQLite::SQLite(), file.path(CA.proj.path, "ca.db")) 
  statement = paste("SELECT licence_holder FROM Individuals;", sep = "")
  statement2 = paste("SELECT licence_holder FROM SnowCrabCDDpercents;", sep = "")
 
  rs = dbSendQuery(con, statement)
  d1 <- dbFetch(rs)      # extract data in chunks of 10 rows
  dbHasCompleted(rs)
  dbClearResult(rs)
  
  rs = dbSendQuery(con, statement2)
  d2 <- dbFetch(rs)      # extract data in chunks of 10 rows
  dbHasCompleted(rs)
  dbClearResult(rs)
  dbDisconnect(con, file.path(CA))
  
  d2 = unique(d2)
  d1 = unique(d1)
  x = NULL
  x$a = setdiff(d2$licence_holder, d1$licence_holder) 
  x$b = setdiff(d1$licence_holder, d2$licence_holder)
  return(x)
   }
