# Script to extract temperatures from set. complete for Rodger
# March 18th, 2016
# Michelle Greenlaw


p = bio.snowcrab::load.environment()

set = snowcrab.db( DS ="set.clean", p=p )
head(set)
set.t.2015 = set[which(set$yr == '2015'),]
head(set.t.2015)
extract = c('trip', 'set', 'station', 'lon', 'lat', 'timestamp', 'julian', 'yr', 't', 'tsd', 'z', 'zsd' )
set.t.2015 = set.t.2015[, extract]
set.t.2015$z = exp(set.t.2015$z)
set.t.2015$zsd = exp(set.t.2015$zsd)
head(set.t.2015)
outdir = file.path( "~", "tmp" )
fn = file.path( outdir, "SnowCrabTemperature2015.csv" )
write.table(set.t.2015, file=fn, sep = ",", col.names = TRUE, row.names=F)
