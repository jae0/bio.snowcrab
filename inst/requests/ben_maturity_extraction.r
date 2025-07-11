

#p = bio.snowcrab::load.environment()

 setwd( file.path( project.datadirectory("bio.snowcrab"), "output" ) )
 det=read_write_fast("det.georef.rdz")
 set=read_write_fast("set.complete.rdz")

 set = set[, c("trip", "set", "timestamp", "julian", "z", "t" )]
 det = merge( x=det, y=set, by=c("trip", "set"), all.x=T, all.y=F )
 read_write_fast(det, file="det_ben.rdz" )


