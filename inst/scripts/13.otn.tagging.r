
#### ORPHANED? ... please check what you need and move to functions


# Ocean Tracking Network detections

redo=FALSE
if(redo) {

  project.library('bio.snowcrab')

  options(stringsAsFactors=F)
  t13 = read.csv(file.path(project.datadirectory('bio.snowcrab'),'data','tagging','OTN_cabotline_detections_2013.csv'),header=T)

  t14 =  read.csv(file.path(project.datadirectory('bio.snowcrab'),'data','tagging','OTN_cabotline_detections_2014.csv'),header=T)

  tt = rbind(t13,t14)
  tt1 = tt[which(!tt$sensortype %in% c('release')),]

  aggregate(bottom_depth~tagname+receiver,data=tt1,FUN=mean)


}
