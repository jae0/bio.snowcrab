web_fisheriesdata_update = function(lic.4X=9, lic.23=62, lic.24=53, lic.NENS=79){

  ###### Function call downloads quotareport for local storage. 
  ###### This also assume a set number of licences. You can 
  ###### change in the variables passed id this ever changes

  if (!require("rvest")) install.packages("rvest", dependencies = TRUE)
  if (!require("dplyr")) install.packages("dplyr", dependencies = TRUE)

  library(rvest)
  library(dplyr)
  
  # NOTE CA.* functions reside in R/ca_functions

  fd = CA.getTable("fisheriesdata")
  yr.max.p = max(p$yrs)
  yr.max.fd = max(fd$YEAR)
  
  while(yr.max.fd <= yr.max.p){
    link.add = paste('https://inter-j02.dfo-mpo.gc.ca/mqr/quotareports/snowcrab?rptyear=', yr.max.fd, '&rptnote=false&lang=en', sep = '')
    ht.link = read_html(link.add)
    tables <- ht.link %>% html_table(fill = TRUE)
    quota.rep = tables[[1]] %>% as_tibble()
    for(i in 1:nrow(quota.rep)){
      con = ROracle::dbConnect(DBI::dbDriver("Oracle"),dbname=oracle.snowcrab.server , username=oracle.snowcrab.user, password=oracle.snowcrab.password, believeNRows=F)
      
      if(quota.rep$FLEET[i] != ""){
        if(trimws(quota.rep$FLEET[i]) == "CFA 23"){
          probe = paste("select * from fisheriesdata where AREA = 'CFA23' and YEAR = ", yr.max.fd, sep = "")
          rs = ROracle::dbGetQuery(con, probe)
          if(nrow(rs) == 0){
            sql = paste("insert into fisheriesdata (YEAR, REGION, AREA, SUBAREA, LICENSES, TAC, LANDINGS) values (", yr.max.fd, ", 'SENS', 'CFA23', 'CFA23', ", lic.23, ", ",str_replace(quota.rep$`QUOTA LIMIT`[i], ",", ""), ", ", as.numeric(str_replace(quota.rep$`CTD TOTAL`[i], ",", "")), ")") 
          }else{
            sql = paste("UPDATE fisheriesdata 
            set TAC = ", str_replace(quota.rep$`QUOTA LIMIT`[i], ",", ""), ", LANDINGS = ", as.numeric(str_replace(quota.rep$`CTD TOTAL`[i], ",", ""))," 
            WHERE AREA = 'CFA23' and YEAR = ", yr.max.fd, sep = "")
          }
          rq = ROracle::dbSendQuery(con, sql)
        }
        if(trimws(quota.rep$FLEET[i]) == "CFA 24E"){
          probe = paste("select * from fisheriesdata where AREA = 'CFA24' and YEAR = ", yr.max.fd, sep = "")
          rs = ROracle::dbGetQuery(con, probe)
          if(nrow(rs) == 0){
            sql = paste("insert into fisheriesdata (YEAR, REGION, AREA, SUBAREA, LICENSES, TAC, LANDINGS) values (", yr.max.fd, ", 'SENS', 'CFA24', 'CFA24', ", lic.24, ", ",str_replace(quota.rep$`QUOTA LIMIT`[i], ",", ""), ", ", as.numeric(str_replace(quota.rep$`CTD TOTAL`[i], ",", "")), ")") 
          }else{
            sql = paste("UPDATE fisheriesdata 
            set TAC = ", str_replace(quota.rep$`QUOTA LIMIT`[i], ",", ""), ", LANDINGS = ", as.numeric(str_replace(quota.rep$`CTD TOTAL`[i], ",", ""))," 
            WHERE AREA = 'CFA24' and YEAR = ", yr.max.fd, sep = "")
          }
          rq = ROracle::dbSendQuery(con, sql)
        }
        if(trimws(quota.rep$FLEET[i]) == "CFA 20-22"){
          probe = paste("select * from fisheriesdata where AREA = 'NENS' and YEAR = ", yr.max.fd, sep = "")
          rs = ROracle::dbGetQuery(con, probe)
          if(nrow(rs) == 0){
            sql = paste("insert into fisheriesdata (YEAR, REGION, AREA, SUBAREA, LICENSES, TAC, LANDINGS) values (", yr.max.fd, ", 'NENS', 'NENS', 'NENS', ", lic.NENS, ", ",str_replace(quota.rep$`QUOTA LIMIT`[i], ",", ""), ", ", as.numeric(str_replace(quota.rep$`CTD TOTAL`[i], ",", "")), ")") 
          }else{
            sql = paste("UPDATE fisheriesdata 
            set TAC = ", str_replace(quota.rep$`QUOTA LIMIT`[i], ",", ""), ", LANDINGS = ", as.numeric(str_replace(quota.rep$`CTD TOTAL`[i], ",", ""))," 
            WHERE AREA = 'NENS' and YEAR = ", yr.max.fd, sep = "")
          }
          rq = ROracle::dbSendQuery(con, sql)
        }
        if(trimws(quota.rep$FLEET[i]) == "CFA 4XE"){
          probe = paste("select * from fisheriesdata where AREA = '4XE' and YEAR = ", yr.max.fd, sep = "")
          rs = ROracle::dbGetQuery(con, probe)
          if(nrow(rs) == 0){
            sql = paste("insert into fisheriesdata (YEAR, REGION, AREA, SUBAREA, LICENSES, TAC, LANDINGS) values (", yr.max.fd, ", '4X', '4XE', '4XE', ", lic.4X, ", ",str_replace(quota.rep$`QUOTA LIMIT`[i], ",", ""), ", ", as.numeric(str_replace(quota.rep$`CTD TOTAL`[i], ",", "")), ")") 
          }else{
            sql = paste("UPDATE fisheriesdata 
            set TAC = ", str_replace(quota.rep$`QUOTA LIMIT`[i], ",", ""), ", LANDINGS = ", as.numeric(str_replace(quota.rep$`CTD TOTAL`[i], ",", ""))," 
            WHERE AREA = '4XE' and YEAR = ", yr.max.fd, sep = "")
          }
          rq = ROracle::dbSendQuery(con, sql)
        }
        if(trimws(quota.rep$FLEET[i]) == "CFA 4XW"){
          probe = paste("select * from fisheriesdata where AREA = '4XW' and YEAR = ", yr.max.fd, sep = "")
          rs = ROracle::dbGetQuery(con, probe)
          if(nrow(rs) == 0){
            sql = paste("insert into fisheriesdata (YEAR, REGION, AREA, SUBAREA, LICENSES, TAC, LANDINGS) values (", yr.max.fd, ", '4X', '4XW', '4XW', ", lic.4X, ", ",str_replace(quota.rep$`QUOTA LIMIT`[i], ",", ""), ", ", as.numeric(str_replace(quota.rep$`CTD TOTAL`[i], ",", "")), ")") 
          }else{
            sql = paste("UPDATE fisheriesdata 
            set TAC = ", str_replace(quota.rep$`QUOTA LIMIT`[i], ",", ""), ", LANDINGS = ", as.numeric(str_replace(quota.rep$`CTD TOTAL`[i], ",", ""))," 
            WHERE AREA = '4XW' and YEAR = ", yr.max.fd, sep = "")
          }
          rq = ROracle::dbSendQuery(con, sql)
        }
        if(trimws(quota.rep$FLEET[i]) == "CFA 24W"){
          probe = paste("select * from fisheriesdata where AREA = '4X' and YEAR = ", yr.max.fd, sep = "")
          rs = ROracle::dbGetQuery(con, probe)
          if(nrow(rs) == 0){
            sql = paste("insert into fisheriesdata (YEAR, REGION, AREA, SUBAREA, LICENSES, TAC, LANDINGS) values (", yr.max.fd, ", '4X', '4X', '4X', ", lic.4X, ", ",str_replace(quota.rep$`QUOTA LIMIT`[i], ",", ""), ", ", as.numeric(str_replace(quota.rep$`CTD TOTAL`[i], ",", "")), ")") 
          }else{
            sql = paste("UPDATE fisheriesdata 
            set TAC = ", str_replace(quota.rep$`QUOTA LIMIT`[i], ",", ""), ", LANDINGS = ", as.numeric(str_replace(quota.rep$`CTD TOTAL`[i], ",", ""))," 
            WHERE AREA = '4X' and YEAR = ", yr.max.fd, sep = "")
          }
          rq = ROracle::dbSendQuery(con, sql)
        }
        
      }
      ROracle::dbDisconnect(con)
    }
    
    fd = CA.getTable("fisheriesdata")
    yr.max.fd = max(fd$YEAR)+1
  }
  return("Update to CA data base Completed")
}
