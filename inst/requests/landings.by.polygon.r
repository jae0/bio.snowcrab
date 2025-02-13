#area of fishery foot print represented in the emera line

p = bio.snowcrab::load.environment()

a = logbook.db('logbook')
a = a[which(a$cfa=='cfanorth'),]
a = a[complete.cases(a[,c('lon','lat')]),]
a = makePBS(a,polygon=F)
pp = importShapefile(aegis.polygons::polygon_file('emera'))

g = findPolys(a,pp)
g = a[which(a$EID %in% g$EID),]

land = aggregate(landings~year,data=a,FUN=sum)
landr = aggregate(landings~year,data=g,FUN=sum)

#percent of landings from corridor
cbind(2002:2015,landr$landings/land$landings)


#spring versus summer landings
a = logbook.db('logbook')
a = a[which(a$cfa=='cfanorth'),]
a$mon = month(a$date.landed)
