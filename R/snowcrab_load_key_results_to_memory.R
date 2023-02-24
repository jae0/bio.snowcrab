
snowcrab_load_key_results_to_memory = function( year.assessment=2022, envir = parent.frame()) {
    # function to bring in key fishery stats and assessment results and make available in memory 
    # primary usage is for Rmarkdown documents
  
  year_previous = year.assessment - 1

  p = bio.snowcrab::load.environment( year.assessment=year.assessment )

  SCD = project.datadirectory("bio.snowcrab")
  

  FD = fisherydata_summary()

  l_nens = round(FD$landings[which(FD$region=="cfanorth" & FD$yr==year.assessment)], 1)
  l_sens = round(FD$landings[which(FD$region=="cfasouth" & FD$yr==year.assessment)], 1)
  l_4x = round(FD$landings[which(FD$region=="cfa4x" & FD$yr==year.assessment)], 1)
  
  l_nens_p = round(FD$landings[which(FD$region=="cfanorth" & FD$yr==year_previous)], 1)
  l_sens_p = round(FD$landings[which(FD$region=="cfasouth" & FD$yr==year_previous)], 1 )
  l_4x_p = round(FD$landings[which(FD$region=="cfa4x" & FD$yr==year_previous)], 1)

  dt_l_nens = round((l_nens - l_nens_p)  / l_nens_p *100, 1 )
  dt_l_sens = round((l_sens - l_sens_p)  / l_sens_p *100, 1 )
  dt_l_4x = round((l_4x - l_4x_p)  / l_4x_p *100, 1 )
  
  e_nens = round(FD$effort[which(FD$region=="cfanorth" & FD$yr==year.assessment)], 3)
  e_sens = round(FD$effort[which(FD$region=="cfasouth" & FD$yr==year.assessment)], 3)
  e_4x = round(FD$effort[which(FD$region=="cfa4x" & FD$yr==year.assessment)], 3)

  e_nens_p = round(FD$effort[which(FD$region=="cfanorth" & FD$yr==year_previous)], 3)
  e_sens_p = round(FD$effort[which(FD$region=="cfasouth" & FD$yr==year_previous)], 3)
  e_4x_p = round(FD$effort[which(FD$region=="cfa4x" & FD$yr==year_previous)], 3)

  dt_e_nens = round(( e_nens - e_nens_p ) /e_nens_p * 100, 1 )
  dt_e_sens = round(( e_sens - e_sens_p ) /e_sens_p * 100, 1 )
  dt_e_4x = round(( e_4x - e_4x_p ) /e_4x_p * 100, 1 )

  c_nens = round(FD$cpue[which(FD$region=="cfanorth" & FD$yr==year.assessment)], 2)
  c_sens = round(FD$cpue[which(FD$region=="cfasouth" & FD$yr==year.assessment)], 2)
  c_4x = round(FD$cpue[which(FD$region=="cfa4x" & FD$yr==year.assessment)], 2)

  c_nens_p = round(FD$cpue[which(FD$region=="cfanorth" & FD$yr==year_previous)], 2)
  c_sens_p = round(FD$cpue[which(FD$region=="cfasouth" & FD$yr==year_previous)], 2)
  c_4x_p = round(FD$cpue[which(FD$region=="cfa4x" & FD$yr==year_previous)], 2)

  dt_c_nens = round(( c_nens - c_nens_p ) /c_nens_p * 100, 1 )
  dt_c_sens = round(( c_sens - c_sens_p ) /c_sens_p * 100, 1 )
  dt_c_4x = round(( c_4x - c_4x_p ) /c_4x_p * 100, 1 )

  dt = as.data.frame( FD[ which(FD$yr %in% c(year.assessment - c(0:10))),] )
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
  CC$percent = round(CC$N / CC$total, 3) * 100
  cc_soft_nens = CC[region=="cfanorth" & fishyr==year.assessment & shell %in% c(1,2), sum(percent)]
  cc_soft_sens = CC[region=="cfasouth" & fishyr==year.assessment & shell %in% c(1,2), sum(percent)]
  cc_soft_4x = CC[region=="cfa4x" & fishyr==year.assessment & shell %in% c(1,2), sum(percent)]
  cc_soft_nens_p = CC[region=="cfanorth" & fishyr==year_previous & shell %in% c(1,2), sum(percent)]
  cc_soft_sens_p = CC[region=="cfasouth" & fishyr==year_previous & shell %in% c(1,2), sum(percent)]
  cc_soft_4x_p = CC[region=="cfa4x" & fishyr==year_previous & shell %in% c(1,2), sum(percent)]

  
  # method = "size_structured_dde_normalized"
  method = "logistic_discrete_historical"

  loc = file.path(SCD, "fishery_model", year.assessment, method )

  b1north = read.csv( file.path(loc, "results_turing_cfanorth_bio_fishing.csv") )
  b1south = read.csv( file.path(loc, "results_turing_cfasouth_bio_fishing.csv") )
  b14x = read.csv( file.path(loc, "results_turing_cfa4x_bio_fishing.csv") )

  # fsnorth = read.csv( file.path(loc, "results_turing_cfanorth_summary.csv") )
  # fssouth = read.csv( file.path(loc, "results_turing_cfasouth_summary.csv") )
  # fs4x = read.csv( file.path(loc, "results_turing_cfa4x_summary.csv") )

  # fmnorth = read.csv( file.path(loc, "results_turing_cfanorth_fm.csv") )
  # fmsouth = read.csv( file.path(loc, "results_turing_cfasouth_fm.csv") )
  # fm4x = read.csv( file.path(loc, "results_turing_cfa4x_fm.csv") )

  if (method == "size_structured_dde_normalized") {
    # b2north = read.csv( file.path(loc, "results_turing_cfanorth_bio_nofishing.csv") )
    # b2south = read.csv( file.path(loc, "results_turing_cfasouth_bio_nofishing.csv") )
    # b24x = read.csv( file.path(loc, "results_turing_cfa4x_bio_nofishing.csv") )
  }

  t1 = which(p$yrs == p$year.assessment -1 )
  t0 = which(p$yrs == p$year.assessment )

  B_north = rowMeans(b1north, na.rm=TRUE )
  B_south = rowMeans(b1south, na.rm=TRUE )
  B_4x = rowMeans(b14x, na.rm=TRUE )

  B_north_sd = apply(b1north, 1, sd, na.rm=TRUE )
  B_south_sd = apply(b1south, 1, sd, na.rm=TRUE )
  B_4x_sd = apply(b14x, 1, sd, na.rm=TRUE )
 

  method = "logistic_discrete_historical"
  loc = file.path(SCD, "fishery_model", year.assessment, method )
  fmnorth = read.csv( file.path(loc, "results_turing_cfanorth_fm.csv") )
  fmsouth = read.csv( file.path(loc, "results_turing_cfasouth_fm.csv") )
  fm4x = read.csv( file.path(loc, "results_turing_cfa4x_fm.csv") )
  t1 = which(p$yrs == p$year.assessment -1 )
  t0 = which(p$yrs == p$year.assessment )
  FM_north = rowMeans(fmnorth, na.rm=TRUE )
  FM_south = rowMeans(fmsouth, na.rm=TRUE )
  FM_4x = rowMeans(fm4x, na.rm=TRUE )

  # method = "size_structured_dde_normalized"
  method = "logistic_discrete_historical"

  loc = file.path(SCD, "fishery_model", year.assessment, method )

  fsnorth = read.csv( file.path(loc, "results_turing_cfanorth_summary.csv") )
  fssouth = read.csv( file.path(loc, "results_turing_cfasouth_summary.csv") )
  fs4x = read.csv( file.path(loc, "results_turing_cfa4x_summary.csv") )

  t1 = which(p$yrs == p$year.assessment -1 )
  t0 = which(p$yrs == p$year.assessment )

  inorth = fsnorth[which(fsnorth$parameters=="K"),]
  isouth = fssouth[which(fssouth$parameters=="K"),]
  i4x = fs4x[which(fs4x$parameters=="K"),]
  
  K_north = round(inorth[["mean"]], 2 )
  K_south = round(isouth[["mean"]], 2 )
  K_4x = round(i4x[["mean"]], 2 )

  K_north_sd = round(inorth[["std"]], 2 )
  K_south_sd = round(isouth[["std"]], 2 )
  K_4x_sd = round(i4x[["std"]], 2 )

  jnorth = fsnorth[which(fsnorth$parameters=="r"),]
  jsouth = fssouth[which(fssouth$parameters=="r"),]
  j4x = fs4x[which(fs4x$parameters=="r"),]
  
  r_north = round(jnorth[["mean"]], 2 )
  r_south = round(jsouth[["mean"]], 2 )
  r_4x = round(j4x[["mean"]], 2 )

  r_north_sd = round(inorth[["std"]], 2 )
  r_south_sd = round(isouth[["std"]], 2 )
  r_4x_sd = round(i4x[["std"]], 2 )

  method = "size_structured_dde_normalized"
  loc = file.path(SCD, "fishery_model", year.assessment, method )

  ddefsnorth = read.csv( file.path(loc, "results_turing_cfanorth_summary.csv") )
  ddefssouth = read.csv( file.path(loc, "results_turing_cfasouth_summary.csv") )
  ddefs4x = read.csv( file.path(loc, "results_turing_cfa4x_summary.csv") )

  fnsumm = file.path( SCD, "modelled", 
        "1999_present_fb", "fishery_model_results", "turing1", "biodyn_number_size_struct.RData" )
  load(fnsumm)

  mw_keep = c(-4:0) + nrow(Y) # last five years
  mw_north = mean( Y[mw_keep,"mw_cfanorth_M0"] )
  mw_south = mean( Y[mw_keep,"mw_cfasouth_M0"] )
  mw_4x = mean( Y[mw_keep,"mw_cfa4x_M0"] )

  t1 = which(p$yrs == p$year.assessment -1 )
  t0 = which(p$yrs == p$year.assessment )

  ddenorth = ddefsnorth[which(ddefsnorth$parameters=="K[1]"),]
  ddesouth = ddefssouth[which(ddefssouth$parameters=="K[1]"),]
  dde4x = ddefs4x[which(ddefs4x$parameters=="K[1]"),]
  
  Kdde_north = round(as.numeric(ddenorth[["mean"]]), 2 ) /10^6 * mw_north
  Kdde_south = round(as.numeric(ddesouth[["mean"]]), 2 ) /10^6 * mw_south
  Kdde_4x = round(as.numeric(dde4x[["mean"]]), 2 ) /10^6 * mw_4x

  Kdde_north_sd = round(as.numeric(ddenorth[["std"]]), 2 ) /10^6* mw_north
  Kdde_south_sd = round(as.numeric(ddesouth[["std"]]), 2 ) /10^6 * mw_south
  Kdde_4x_sd = round(as.numeric(dde4x[["std"]]), 2 ) /10^6* mw_4x

  bdde_north = ddefsnorth[which(ddefsnorth$parameters=="b[2]"),]
  bdde_south = ddefssouth[which(ddefssouth$parameters=="b[2]"),]
  bdde_4x = ddefs4x[which(ddefs4x$parameters=="b[2]"),]
  
  bdde2_north = round(as.numeric(bdde_north[["mean"]]), 2 )
  bdde2_south = round(as.numeric(bdde_south[["mean"]]), 2 )
  bdde2_4x = round(as.numeric(bdde_4x[["mean"]]), 2 )

  bdde2_north_sd = round(as.numeric(bdde_north[["std"]]), 2 )
  bdde2_south_sd = round(as.numeric(bdde_south[["std"]]), 2 )
  bdde2_4x_sd = round(as.numeric(bdde_4x[["std"]]), 2 )

  ddefmnorth = read.csv( file.path(loc, "results_turing_cfanorth_fm.csv") )
  ddefmsouth = read.csv( file.path(loc, "results_turing_cfasouth_fm.csv") )
  ddefm4x = read.csv( file.path(loc, "results_turing_cfa4x_fm.csv") )

  ddeFM_north = rowMeans(ddefmnorth, na.rm=TRUE )
  ddeFM_south = rowMeans(ddefmsouth, na.rm=TRUE )
  ddeFM_4x = rowMeans(ddefm4x, na.rm=TRUE )
 
  return( invisible( list2env(as.list(environment()), envir) ) )

}
