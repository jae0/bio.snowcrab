table_counts = function() {

# Could a table of the survey year, number of stations by CFA, number of crab measured, number of
# temperature recordings, number of ’other species’ recorded, vessel, etc be included as an appendix.
# It is a valuable piece of information to know how much data is going into these analyses and the
# changes over time or if any anomalies occurred.


year.assessment = 2022

p = bio.snowcrab::load.environment( year.assessment=year.assessment )
  
set = snowcrab.db( DS="set.clean", p=p )

cfa4x = polygon_inside(set, "cfa4x")
cfanorth = polygon_inside(set, "cfanorth")
cfasouth = polygon_inside(set, "cfasouth")

set$region = NA
set$region[cfa4x] = "cfa4x"
set$region[cfanorth] = "cfanorth"
set$region[cfasouth] = "cfasouth"

set$region = factor(set$region, levels=regions, labels =region_label)
set = set[(which(!is.na(set$region))), ]
setDT(set)
o = set[,.(Nstation=.N, vessel=unique(vessel), Ntemp=length(is.finite(t))), by=.(yr, region)]   


det = snowcrab.db( DS ="det.georeferenced", p=p )

cfa4x = polygon_inside(det, "cfa4x")
cfanorth = polygon_inside(det, "cfanorth")
cfasouth = polygon_inside(det, "cfasouth")

det$region = NA
det$region[cfa4x] = "cfa4x"
det$region[cfanorth] = "cfanorth"
det$region[cfasouth] = "cfasouth"

det$region = factor(det$region, levels=regions, labels =region_label)
det = det[(which(!is.na(det$region))), ]
setDT(det)
od = det[, .(Ncrab=.N), by=.(yr, region)]


cat = snowcrab.db( DS ="cat.georeferenced", p=p )

cfa4x = polygon_inside(cat, "cfa4x")
cfanorth = polygon_inside(cat, "cfanorth")
cfasouth = polygon_inside(cat, "cfasouth")

cat$region = NA
cat$region[cfa4x] = "cfa4x"
cat$region[cfanorth] = "cfanorth"
cat$region[cfasouth] = "cfasouth"

cat$region = factor(cat$region, levels=regions, labels =region_label)
cat = cat[(which(!is.na(cat$region))), ]
cat = cat[-which(cat$spec==2526),]

setDT(cat)

oc = cat[, .(Nsp = length(unique(spec))), by=.(yr, region)]

out = merge(merge(o, oc), od)
out = out[ order( region, yr ), ]
return(out0)
}