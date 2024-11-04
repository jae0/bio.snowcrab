suitable.spread <- function (lat, lon, yr, yr.spread){
  spread = NA
  
  #get all spread stats

  allstats = netmind.db (DS="stats", quiet =T)
  
  #remove by survey year current year and everything but recent years
  allstats$yr = year(allstats$t0)
  allstats$month = month(allstats$t0)
  allstats$yr[which(allstats$month<5)]=allstats$yr[which(allstats$month<5)]-1
  allstats = allstats[which(allstats$yr > (yr-7)),]
  allstats = allstats[which(allstats$yr < yr),]
  
  #Find average spread for recent years to apply a yearly offset
  otheryears = aggregate(allstats$spread, list(allstats$yr), FUN=mean, na.rm=T)
  spreaddiff = yr.spread - mean(otheryears$x)
  
  #Find recent spreads at this location
  spread = mean(allstats[which(((abs(lat - allstats$slat)) < .05) & ((abs(lon - allstats$slon)) < .05)),]$spread, na.rm = T)
  
  #Apply yearly offset
  spread = spread + spreaddiff
  
  return(spread*1000)
} 