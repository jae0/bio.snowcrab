
snowcrab_load_key_results_to_memory = function( 
  year.assessment=2024,
  years_model = 2000:year.assessment, 
  envir = parent.frame(),
  regions = list( region=c("cfanorth", "cfasouth", "cfa4x" )),
  debugging=FALSE,  
  return_as_list=TRUE ) {
    # function to bring in key fishery stats and assessment results and make available in memory 
    # primary usage is for Rmarkdown documents
  
  year_previous = year.assessment - 1

  p = bio.snowcrab::load.environment( year.assessment=year.assessment )

  SCD = project.datadirectory("bio.snowcrab")
  
  FD = fishery_data( regions=regions )  # mass in tonnes
  fda = FD$summary_annual

  l_nens = round(fda$landings[which(fda$region=="cfanorth" & fda$yr==year.assessment)], 1)
  l_sens = round(fda$landings[which(fda$region=="cfasouth" & fda$yr==year.assessment)], 1)
  l_4x = round(fda$landings[which(fda$region=="cfa4x" & fda$yr==year.assessment)], 1)
  
  l_nens_p = round(fda$landings[which(fda$region=="cfanorth" & fda$yr==year_previous)], 1)
  l_sens_p = round(fda$landings[which(fda$region=="cfasouth" & fda$yr==year_previous)], 1 )
  l_4x_p = round(fda$landings[which(fda$region=="cfa4x" & fda$yr==year_previous)], 1)

  dt_l_nens = round((l_nens - l_nens_p)  / l_nens_p *100, 1 )
  dt_l_sens = round((l_sens - l_sens_p)  / l_sens_p *100, 1 )
  dt_l_4x = round((l_4x - l_4x_p)  / l_4x_p *100, 1 )
  
  e_nens = round(fda$effort[which(fda$region=="cfanorth" & fda$yr==year.assessment)], 3)
  e_sens = round(fda$effort[which(fda$region=="cfasouth" & fda$yr==year.assessment)], 3)
  e_4x = round(fda$effort[which(fda$region=="cfa4x" & fda$yr==year.assessment)], 3)

  e_nens_p = round(fda$effort[which(fda$region=="cfanorth" & fda$yr==year_previous)], 3)
  e_sens_p = round(fda$effort[which(fda$region=="cfasouth" & fda$yr==year_previous)], 3)
  e_4x_p = round(fda$effort[which(fda$region=="cfa4x" & fda$yr==year_previous)], 3)

  dt_e_nens = round(( e_nens - e_nens_p ) /e_nens_p * 100, 1 )
  dt_e_sens = round(( e_sens - e_sens_p ) /e_sens_p * 100, 1 )
  dt_e_4x = round(( e_4x - e_4x_p ) /e_4x_p * 100, 1 )

  c_nens = round(fda$cpue[which(fda$region=="cfanorth" & fda$yr==year.assessment)], 2)
  c_sens = round(fda$cpue[which(fda$region=="cfasouth" & fda$yr==year.assessment)], 2)
  c_4x = round(fda$cpue[which(fda$region=="cfa4x" & fda$yr==year.assessment)], 2)

  c_nens_p = round(fda$cpue[which(fda$region=="cfanorth" & fda$yr==year_previous)], 2)
  c_sens_p = round(fda$cpue[which(fda$region=="cfasouth" & fda$yr==year_previous)], 2)
  c_4x_p = round(fda$cpue[which(fda$region=="cfa4x" & fda$yr==year_previous)], 2)

  dt_c_nens = round(( c_nens - c_nens_p ) /c_nens_p * 100, 1 )
  dt_c_sens = round(( c_sens - c_sens_p ) /c_sens_p * 100, 1 )
  dt_c_4x = round(( c_4x - c_4x_p ) /c_4x_p * 100, 1 )

  dt = as.data.frame( fda[ which(fda$yr %in% c(year.assessment - c(0:10))),] )
  dt =  dt[,c("region", "yr", "Licenses", "TAC", "landings", "effort", "cpue")] 
  names(dt) = c("Region", "Year", "Licenses", "TAC", "Landings", "Effort", "CPUE") 
  rownames(dt) = NULL
  
  tac_nens = fda$TAC[which(fda$yr==year.assessment & fda$region=="cfanorth")]
  tac_sens = fda$TAC[which(fda$yr==year.assessment & fda$region=="cfasouth")]
  tac_4x = fda$TAC[which(fda$yr==year.assessment & fda$region=="cfa4x")] # 4x is refered by start year
  tac_4x_p = fda$TAC[which(fda$yr==year_previous & fda$region=="cfa4x")] # 4x is refered by start year
  fda = NULL

  scn = FD$shell_condition
  cc_soft_nens = scn[ region=="cfanorth" & fishyr==year.assessment & shell %in% c(1,2), sum(percent)]
  cc_soft_sens = scn[ region=="cfasouth" & fishyr==year.assessment & shell %in% c(1,2), sum(percent)]
  cc_soft_4x = scn[ region=="cfa4x" & fishyr==year.assessment & shell %in% c(1,2), sum(percent)]
  cc_soft_nens_p = scn[ region=="cfanorth" & fishyr==year_previous & shell %in% c(1,2), sum(percent)]
  cc_soft_sens_p = scn[ region=="cfasouth" & fishyr==year_previous & shell %in% c(1,2), sum(percent)]
  cc_soft_4x_p = scn[ region=="cfa4x" & fishyr==year_previous & shell %in% c(1,2), sum(percent)]
  scn = NULL

  # here mean is used to force result as a scalar
  fob = FD$fraction_observed
  observed_nens = fob[ region=="cfanorth" & yr==year.assessment, mean(observed_landings_pct, na.rm=TRUE) ]
  observed_sens = fob[ region=="cfasouth" & yr==year.assessment, mean(observed_landings_pct, na.rm=TRUE) ]
  observed_4x = fob[ region=="cfa4x" & yr==year.assessment, mean(observed_landings_pct, na.rm=TRUE) ]
  observed_nens_p = fob[ region=="cfanorth" & yr==year_previous, mean(observed_landings_pct, na.rm=TRUE) ]
  observed_sens_p = fob[ region=="cfasouth" & yr==year_previous, mean(observed_landings_pct, na.rm=TRUE) ]
  observed_4x_p = fob[ region=="cfa4x" & yr==year_previous, mean(observed_landings_pct, na.rm=TRUE) ]
  fob = NULL
 
  method = "logistic_discrete_historical"
  loc = file.path(SCD, "fishery_model", year.assessment, method )

  b1north = fread( file.path(loc, "results_turing_cfanorth_bio_fishing.csv"), header=TRUE, sep=";" )
  b1south = fread( file.path(loc, "results_turing_cfasouth_bio_fishing.csv"), header=TRUE, sep=";" )
  b14x = fread( file.path(loc, "results_turing_cfa4x_bio_fishing.csv"), header=TRUE, sep=";" )

  t1 = which(years_model == p$year.assessment -1 )
  t0 = which(years_model == p$year.assessment )

  B_north = rowMeans(b1north, na.rm=TRUE )
  B_south = rowMeans(b1south, na.rm=TRUE )
  B_4x = rowMeans(b14x, na.rm=TRUE )

  B_north_sd = apply(b1north, 1, sd, na.rm=TRUE )
  B_south_sd = apply(b1south, 1, sd, na.rm=TRUE )
  B_4x_sd = apply(b14x, 1, sd, na.rm=TRUE )


  fmnorth = fread( file.path(loc, "results_turing_cfanorth_fm.csv"), header=TRUE, sep=";" )
  fmsouth = fread( file.path(loc, "results_turing_cfasouth_fm.csv"), header=TRUE, sep=";" )
  fm4x = fread( file.path(loc, "results_turing_cfa4x_fm.csv"), header=TRUE, sep=";" )

  FM_north = rowMeans(fmnorth, na.rm=TRUE )
  FM_south = rowMeans(fmsouth, na.rm=TRUE )
  FM_4x = rowMeans(fm4x, na.rm=TRUE )

  FM_north_sd = apply(fmnorth, 1, sd, na.rm=TRUE )
  FM_south_sd = apply(fmsouth, 1, sd, na.rm=TRUE )
  FM_4x_sd = apply(fm4x, 1, sd, na.rm=TRUE )


  fsnorth = fread( file.path(loc, "results_turing_cfanorth_summary.csv"), header=TRUE, sep=";" )
  fssouth = fread( file.path(loc, "results_turing_cfasouth_summary.csv"), header=TRUE, sep=";" )
  fs4x = fread( file.path(loc, "results_turing_cfa4x_summary.csv"), header=TRUE, sep=";" )

  Knorth = fsnorth[which(fsnorth$parameters=="K"),]
  Ksouth = fssouth[which(fssouth$parameters=="K"),]
  K4x = fs4x[which(fs4x$parameters=="K"),]
  
  K_north = round(Knorth[["mean"]], 2 )
  K_south = round(Ksouth[["mean"]], 2 )
  K_4x = round(K4x[["mean"]], 2 )

  K_north_sd = round(Knorth[["std"]], 2 )
  K_south_sd = round(Ksouth[["std"]], 2 )
  K_4x_sd = round(K4x[["std"]], 2 )

  rnorth = fsnorth[which(fsnorth$parameters=="r"),]
  rsouth = fssouth[which(fssouth$parameters=="r"),]
  r4x = fs4x[which(fs4x$parameters=="r"),]

  r_north = round(rnorth[["mean"]], 2 )
  r_south = round(rsouth[["mean"]], 2 )
  r_4x = round(r4x[["mean"]], 2 )

  r_north_sd = round(rnorth[["std"]], 2 )
  r_south_sd = round(rsouth[["std"]], 2 )
  r_4x_sd = round(r4x[["std"]], 2 )


  qnorth = fsnorth[which(fsnorth$parameters=="q1"),]
  qsouth = fssouth[which(fssouth$parameters=="q1"),]
  q4x = fs4x[which(fs4x$parameters=="q1"),]

  q_north = round(qnorth[["mean"]], 2 )
  q_south = round(qsouth[["mean"]], 2 )
  q_4x = round(q4x[["mean"]], 2 )

  q_north_sd = round(qnorth[["std"]], 2 )
  q_south_sd = round(qsouth[["std"]], 2 )
  q_4x_sd = round(q4x[["std"]], 2 )

  if (return_as_list) {
    return( invisible( as.list( environment() ) ) )
  } else {
    return( invisible( list2env(as.list(environment()), envir) ) )
  }

}
