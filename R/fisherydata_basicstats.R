  

fisherydata_basicstats = function(year.assessment=year(Sys.Date())) {

  FD = fisherydata_summary()

  year_previous = year.assessment  - 1
  
  l_nens = round(FD$landings[which(FD$region=="cfanorth" & FD$yr==year.assessment)])
  l_sens = round(FD$landings[which(FD$region=="cfasouth" & FD$yr==year.assessment)])
  l_4x = round(FD$landings[which(FD$region=="cfa4x" & FD$yr==year.assessment)])
  
  l_nens_p = round(FD$landings[which(FD$region=="cfanorth" & FD$yr==year_previous)])
  l_sens_p = round(FD$landings[which(FD$region=="cfasouth" & FD$yr==year_previous)])
  l_4x_p = round(FD$landings[which(FD$region=="cfa4x" & FD$yr==year_previous)])

  dt_l_nens = round((l_nens - l_nens_p)  / l_nens_p *100 )
  dt_l_sens = round((l_sens - l_sens_p)  / l_sens_p *100 )
  dt_l_4x = round((l_4x - l_4x_p)  / l_4x_p *100 )

  c_nens = round(FD$cpue[which(FD$region=="cfanorth" & FD$yr==year.assessment)])
  c_sens = round(FD$cpue[which(FD$region=="cfasouth" & FD$yr==year.assessment)])
  c_4x = round(FD$cpue[which(FD$region=="cfa4x" & FD$yr==year.assessment)])

  c_nens_p = round(FD$cpue[which(FD$region=="cfanorth" & FD$yr==year_previous)])
  c_sens_p = round(FD$cpue[which(FD$region=="cfasouth" & FD$yr==year_previous)])
  c_4x_p = round(FD$cpue[which(FD$region=="cfa4x" & FD$yr==year_previous)])

  dt_c_nens = round(( c_nens - c_nens_p ) /c_nens_p * 100 )
  dt_c_sens = round(( c_sens - c_sens_p ) /c_sens_p * 100 )
  dt_c_4x = round(( c_4x - c_4x_p ) /c_4x_p * 100 )

  dt = as.data.frame( FD[ which(FD$yr %in% c(year.assessment - c(0:10))),] )
  n = ( which(dt$region=="cfanorth") )
  s = ( which(dt$region=="cfasouth") )
  x = ( which(dt$region=="cfa4x") )
  dt =  dt[,c("region", "yr", "Licenses", "TAC", "landings", "effort", "cpue")] 
  names(dt) = c("Region", "Year", "Licenses", "TAC", "Landings", "Effort", "CPUE") 
  rownames(dt) = NULL
  
  tac_nens = FD$TAC[which(FD$yr==year.assessment & FD$region=="cfanorth")]
  tac_sens = FD$TAC[which(FD$yr==year.assessment & FD$region=="cfasouth")]
  tac_4x = FD$TAC[which(FD$yr==year.assessment & FD$region=="cfa4x")] # 4x is refered by start year
  tac_4x_p = FD$TAC[which(FD$yr==year_previous & FD$region=="cfa4x")] # 4x is refered by start year


  # carapace condition .. could have used the stuff in 02_fisures and table but it was convoluted
  # should replace that mess with this eventually 
  male = 0
  odb = observer.db("odb")
  setDT(odb)
  odb = odb[ which( odb$sex==male & odb$cw >= 95 & odb$cw < 170 & odb$prodcd_id=="0" & is.finite(odb$shell) ) ,]  # commerical sized crab only
  odb$region = NA
  for ( reg in c("cfanorth", "cfasouth", "cfa4x")) {
    r = polygon_inside(x=odb, region=aegis.polygons::polygon_internal_code(reg), planar=FALSE)
    if (length(r)> 0)  odb$region[r] = reg
  }
  CC = odb[ !is.na(odb$region), .N, by=.(region, fishyr, shell) ]
  CC[, total:=sum(N, na.rm=TRUE), by=.(region, fishyr)]
  CC$percent = round(CC$N / CC$total, 2) * 100
  cc_soft_nens = CC[region=="cfanorth" & fishyr==year.assessment & shell %in% c(1,2), sum(percent)]
  cc_soft_sens = CC[region=="cfasouth" & fishyr==year.assessment & shell %in% c(1,2), sum(percent)]
  cc_soft_4x = CC[region=="cfa4x" & fishyr==year.assessment & shell %in% c(1,2), sum(percent)]
  cc_soft_nens_p = CC[region=="cfanorth" & fishyr==year_previous & shell %in% c(1,2), sum(percent)]
  cc_soft_sens_p = CC[region=="cfasouth" & fishyr==year_previous & shell %in% c(1,2), sum(percent)]
  cc_soft_4x_p = CC[region=="cfa4x" & fishyr==year_previous & shell %in% c(1,2), sum(percent)]


  FS = list(
    l_nens = l_nens,
    l_sens = l_sens,
    l_4x = l_4x,
    l_nens_p = l_nens_p,
    l_sens_p = l_sens_p,
    l_4x_p = l_4x_p,
    dt_l_nens = dt_l_nens,
    dt_l_sens = dt_l_sens,
    dt_l_4x = dt_l_4x,
    c_nens = c_nens,
    c_sens = c_sens,
    c_4x = c_4x, 
    c_nens_p = c_nens_p,
    c_sens_p = c_sens_p,
    c_4x_p = c_4x_p, 
    dt_c_nens = dt_c_nens,
    dt_c_sens = dt_c_sens,
    dt_c_4x = dt_c_4x,
    dt = dt,
    tac_nens = tac_nens,
    tac_sens = tac_sens, 
    tac_4x = tac_4x,
    tac_4x_p = tac_4x_p,
    CC = CC,
    cc_soft_nens = cc_soft_nens, 
    cc_soft_sens = cc_soft_sens, 
    cc_soft_4x = cc_soft_4x, 
    cc_soft_nens_p = cc_soft_nens_p, 
    cc_soft_sens_p = cc_soft_sens_p, 
    cc_soft_4x_p = cc_soft_4x_p
  )
  return(FS)

}
