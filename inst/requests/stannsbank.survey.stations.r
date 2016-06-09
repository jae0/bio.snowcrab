

p = bio.snowcrab::initialise.local.environment()

dat = snowcrab.db('set.complete')
pp = read.csv(find.bio.gis('StAnnsMPA.csv'))
require(PBSmapping)

dat$X = dat$lon
dat$Y = dat$lat
dat$EID = 1:nrow(dat)

p1 = findPolys(dat,pp)

da = dat[which(dat$EID %in% p1$EID),]

dal = da[,c('lon','lat','distance','timestamp', names(da)[grep('ms.size',names(da))])]
nn = as.numeric(substr(names(da)[grep('ms.size',names(da))],9,15))
nn = taxonomy.recode(from='spec',tolookup=nn)$vern
names(dal)[5:ncol(dal)] <- nn
 dal$yr = lubridate::year( dal$timestamp )
  dal = dal[which(dal$yr>2003),]

dal = Filter(function(x)!all(is.na(x)),dal)
dal$yr = NULL

#species data
outloc = project.datadirectory( "bio.snowcrab", "R", "requests" )
fn = file.path( outloc, "StAnnsSCSurveySpeciesData.csv" )
write.csv(dal,file=fn, row.names=F)

dal$NSp = apply(dal[,5:ncol(dal)],1,function(x) length(x[!is.na(x)]))

x = dal[,c('lon','lat','distance','timestamp','Nsp')]
xx = x[which( lubridate::year(x$timestamp) > 2003),]

fn = file.path( outloc, "StAnnsSnowCrabSurvey.csv" )
write.csv(xx, fn, row.names=F)



