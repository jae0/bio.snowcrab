
##Must insert budget every year!!!!!!!!!!!
CA.insertVariables = function(price_per_mt = 2970, budget, year){
  CA.proj.path = file.path(project.datadirectory("bio.snowcrab","data", "CA","CA_db"))
  
  con = ROracle::dbConnect(DBI::dbDriver("Oracle"),dbname=oracle.snowcrab.server , username=oracle.snowcrab.user, password=oracle.snowcrab.password, believeNRows=F)
  #con <- dbConnect(RSQLite::SQLite(), file.path(CA.proj.path, "ca.db")) 
  
  statement = paste("INSERT INTO VariableLookup(yr, price_per_mt, budget) VALUES(",year,",",price_per_mt,",", budget,")", sep = "")
  rs = ROracle::dbSendQuery(con, statement)
  #dbHasCompleted(rs)
  #dbClearResult(rs)
  #dbDisconnect(con, file.path(CA))
  ROracle::dbDisconnect(con)
  }
##Must insert TAC values each year!!!!!!!!!!!
CA.insertTAC = function(cfa23.tac, cfa24.tac, Millbrook.tac = 250, nens.tac, xxxx.tac, year){
  CA.proj.path = file.path(project.datadirectory("bio.snowcrab","data", "CA","CA_db"))

   #con <- dbConnect(RSQLite::SQLite(), file.path(CA.proj.path, "ca.db")) 
  con = ROracle::dbConnect(DBI::dbDriver("Oracle"),dbname=oracle.snowcrab.server , username=oracle.snowcrab.user, password=oracle.snowcrab.password, believeNRows=F)
  
  statement = paste("SELECT * FROM TAC where area = 'CFA 23' AND yr = ",year, "", sep = "")
#  rs = dbSendQuery(con, statement)
  rs = ROracle::dbGetQuery(con, statement )
 # d1 <- dbFetch(rs)    
 # dbHasCompleted(rs)
 # dbClearResult(rs)
  if(nrow(rs) == 0){
    statement = paste("INSERT INTO TAC(YR, AREA, TAC) VALUES('", year,"', '","CFA 23","', '", cfa23.tac,"')", sep = "")
    #rs = dbSendQuery(con, statement)
    rs = ROracle::dbSendQuery(con, statement)
    #dbHasCompleted(rs)
   # dbClearResult(rs)
  }
  statement = paste("SELECT * FROM TAC where area = 'CFA 24' AND yr = ",year, "", sep = "")
 # rs = dbSendQuery(con, statement)
  rs = ROracle::dbGetQuery(con, statement )
  #d1 <- dbFetch(rs)      # extract data in chunks of 10 rows
  #dbHasCompleted(rs)
  #dbClearResult(rs)
  if(nrow(rs) == 0){
    statement = paste("INSERT INTO TAC(YR, AREA, TAC) VALUES('", year,"', '","CFA 24","', '", cfa24.tac,"')", sep = "")
   
    #rs = dbSendQuery(con, statement)
    rs = ROracle::dbSendQuery(con, statement)
    #dbHasCompleted(rs)
    #dbClearResult(rs)
  }
  
  statement = paste("SELECT * FROM TAC where area = 'CFA 24 Millbrook' AND yr = ",year, "", sep = "")
  rs = ROracle::dbGetQuery(con, statement )
 # rs = dbSendQuery(con, statement)
  #d1 <- dbFetch(rs)      # extract data in chunks of 10 rows
 # dbHasCompleted(rs)
 # dbClearResult(rs)
  if(nrow(rs) == 0){
    statement = paste("INSERT INTO TAC(YR, AREA, TAC) VALUES('", year,"', '","CFA 24 Millbrook","', '", Millbrook.tac,"')", sep = "")
    
    rs = ROracle::dbSendQuery(con, statement)
    #rs = dbSendQuery(con, statement)
    #dbHasCompleted(rs)
    #dbClearResult(rs)
  }
  statement = paste("SELECT * FROM TAC where area = 'N-ENS' AND yr = ",year, "", sep = "")
  rs = ROracle::dbGetQuery(con, statement )
  #rs = dbSendQuery(con, statement)
  #d1 <- dbFetch(rs)      # extract data in chunks of 10 rows
  #dbHasCompleted(rs)
  #dbClearResult(rs)
  if(nrow(rs) == 0){
    statement = paste("INSERT INTO TAC(YR, AREA, TAC) VALUES('", year,"', '","N-ENS","', '", nens.tac,"')", sep = "")
    rs = ROracle::dbSendQuery(con, statement)
   # rs = dbSendQuery(con, statement)
   # dbHasCompleted(rs)
   # dbClearResult(rs)
  }   
  
  statement = paste("SELECT * FROM TAC where area = '4X' AND yr = ",year, "", sep = "")
  rs = ROracle::dbGetQuery(con, statement )
 # rs = dbSendQuery(con, statement)
 # d1 <- dbFetch(rs)      # extract data in chunks of 10 rows
# dbHasCompleted(rs)
 # dbClearResult(rs)
  if(nrow(rs) == 0){
    statement = paste("INSERT INTO TAC(YR, AREA, TAC) VALUES('", year,"', '","4X","', '", xxxx.tac,"')", sep = "")
    
    rs = ROracle::dbSendQuery(con, statement)
    # rs = dbSendQuery(con, statement)
   # dbHasCompleted(rs)
   # dbClearResult(rs)
    
  }   
  ROracle::dbDisconnect(con)
  #dbDisconnect(con, file.path(CA))
  
}

## Must load yearly CDD data This data is received as excel sheet. Save each tab as csv to CDDcsv folderand call
## this function
CA.insertCDDcsv = function(yr = NULL) {
  CA.proj.path = file.path(project.datadirectory("bio.snowcrab","data", "CA","CA_db"))
  lof = list.files(file.path(CA.proj.path, "CDDcsv"), full.names=T)
  con = ROracle::dbConnect(DBI::dbDriver("Oracle"),dbname=oracle.snowcrab.server , username=oracle.snowcrab.user, password=oracle.snowcrab.password, believeNRows=F)
  
  if(!is.null(yr)){
    lof = lof[grepl(as.character(yr), lof)]
  }

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
        statement = paste("SELECT * FROM SnowCrabCDDpercents where RTRIM(licence) = '",dasub$LICENCE[k],"' AND yr = ",year, "", sep = "")
        rs = ROracle::dbGetQuery(con, statement)
  
        if(nrow(rs) == 0){
           statement = paste("INSERT INTO SnowCrabCDDpercents(LICENCE, LICENCE_HOLDER, INITIAL_PERCENT_TAC, CURRENT_PERCENT_TAC, LICENCE_UOF, INITIAL_PERCENT_TAC_UOF, CURRENT_PERCENT_TAC_UOF, YR) VALUES('",as.character(dasub$LICENCE[k]),"','",  str_replace_all(dasub$LICENCEHOLDER[k], "'", ""),"','",dasub$INITIALofTAC[k],"','",dasub$CURRENTOFTAC[k],"','",dasub$LICENCE_UOF[k],"','",dasub$INITIALofTAC[k],"','",dasub$CURRENTOFTAC[k],"','",year,"')", sep = "")
           rs = ROracle::dbSendQuery(con, statement)
        }
        else{}
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
          statement = paste("SELECT * FROM SnowCrabCDDpercents where licence = '",lic,"' AND yr = ",year, "", sep = "")
          rs = ROracle::dbGetQuery(con, statement)
          
          if(nrow(rs) == 0){
            statement = paste("INSERT INTO SnowCrabCDDpercents(LICENCE, LICENCE_HOLDER, INITIAL_PERCENT_TAC, CURRENT_PERCENT_TAC, LICENCE_UOF, INITIAL_PERCENT_TAC_UOF, CURRENT_PERCENT_TAC_UOF, YR) VALUES('",as.character(lic),"','",  str_replace_all(dax$licence_holder[j], "'", ""),"','",dax$initial_percent_tac[j],"','",dax$current_percent_tac[j],"','",dax$licence_uof[j],"','",dax$initial_percent_tac_uof[j],"','",dax$current_percent_tac_uof[j],"','",year,"')", sep = "")
            rs = ROracle::dbSendQuery(con, statement)
          }
          else{}
          }
      }
    }
  }
  ROracle::dbDisconnect(con)
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
    print(paste("Sum of percents: ", sum(d1$PERCENT_TOTAL), sep = ""))
    print(paste("Total CA contributions: ", sum(d1$CA_CONTRIBUTION), sep = ""))
    print(paste("Total Use of Fish Required: ", sum(d1$UOF_REQUIRED_MT), sep = ""))
}

CA.writePartnerContributions = function(year){
  CA.proj.path = file.path(project.datadirectory("bio.snowcrab","data", "CA","CA_db"))
  d1 = CA.getTable("CalcPartner", year)
  
  if(!dir.exists(file.path(CA.proj.path,"reports", year)))
    dir.create(file.path(CA.proj.path,"reports", year))
    write.csv(d1, file.path(CA.proj.path, "reports",year, "partner.csv"))
  
    names(d1) = gsub(" ", "", names(d1))
    print("PARTNER CONTRIBUTIONS CHECK")
    print(paste("Sum of percents: ", sum(d1$SUM_OF_PERCENT), sep = ""))
    print(paste("Total CA contributions: ", sum(d1$SUM_CA_CONTRIBUTION), sep = ""))
    print(paste("Total Use of Fish Required: ", sum(d1$SUM_UOF_REQUIRED_MT), sep = ""))
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
  print(paste("Sum of percents: ", sum(d1$SUM_OF_TOTAL_TAC_PERCENTS), sep = ""))
  print(paste("Total CA contributions: ", sum(d1$SUM_OF_CA_CONTRIBUTIONS), sep = ""))
  print(paste("Total Use of Fish Required: ", sum(d1$SUM_OF_UOF_REQUIRED_MT), sep = ""))
  print("")
}


CA.removeTAC = function(year){
  CA.proj.path = file.path(project.datadirectory("bio.snowcrab","data", "CA","CA_db"))
  con = ROracle::dbConnect(DBI::dbDriver("Oracle"),dbname=oracle.snowcrab.server , username=oracle.snowcrab.user, password=oracle.snowcrab.password, believeNRows=F)
  
  statement = paste("DELETE FROM TAC where yr = ",year, "", sep = "")
  rs = ROracle::dbSendQuery(con, statement )
  
  ROracle::dbDisconnect(con)
}

CA.getTable = function(tablename, sub.year = NULL){ 
  CA.proj.path = file.path(project.datadirectory("bio.snowcrab","data", "CA","CA_db"))
  con = ROracle::dbConnect(DBI::dbDriver("Oracle"),dbname=oracle.snowcrab.server , username=oracle.snowcrab.user, password=oracle.snowcrab.password, believeNRows=F)
  
  if(is.null(sub.year))
    statement = paste("SELECT * FROM ", tablename, "", sep = "")
  else
    statement = paste("SELECT * FROM ", tablename, " where yr = '", sub.year, "'", sep = "")
  
  rs = ROracle::dbGetQuery(con, statement)
  ROracle::dbDisconnect(con)
  return(rs)
}

CA.checkMatch = function(){ 
  CA.proj.path = file.path(project.datadirectory("bio.snowcrab","data", "CA","CA_db"))
  con = ROracle::dbConnect(DBI::dbDriver("Oracle"),dbname=oracle.snowcrab.server , username=oracle.snowcrab.user, password=oracle.snowcrab.password, believeNRows=F)
  
  statement = paste("SELECT licence_holder FROM Individuals", sep = "")
  statement2 = paste("SELECT licence_holder FROM SnowCrabCDDpercents", sep = "")
 
  rs = ROracle::dbGetQuery(con, statement)
  rs2 = ROracle::dbGetQuery(con, statement2)
  ROracle::dbDisconnect(con)
  
  d2 = unique(rs2)
  d1 = unique(rs)
  x = NULL
  x$a = setdiff(d2$licence_holder, d1$licence_holder) 
  x$b = setdiff(d1$licence_holder, d2$licence_holder)
  return(x)
}
