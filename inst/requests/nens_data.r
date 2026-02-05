
# nens data extration for Stephanie Boudreau

# We’re looking for snow crab by size, sex, and maturity (if known), by tow 
# - coordinates, depth, and swept area for N-ENS/4Vn. 
# For the years, we’re starting in 1997 which is when the southern Gulf snow crab survey time series begins. 



require(aegis)  # basic helper tools

yrs = 1997:2025
year.assessment = 2025  # change this as appropriate

p = bio.snowcrab::load.environment( year.assessment=year.assessment )  # set up initial settings

set = snowcrab.db( DS="set.clean", p=p )  # load clean set data
setDT(set)
set = set[, .(trip, set, z)]

det = snowcrab.db( DS="det.georeferenced", p=p )  # merge set.clean
setDT(det)

out = det[ 
  region=="cfanorth" & 
  yr %in% yrs, 
  .(trip, set, crabno, yr, sex, cw, mat, lon, lat, sa) 
]

out = set[out, on=.(trip, set) ]

setnames(out, "z", "depth_m")
setnames(out, "sa", "sweptarea_m2")
setnames(out, "cw", "cw_mm")

out$sex = as.numeric( as.character(out$sex) )
out$mat = as.numeric( as.character(out$mat) )

# re-code variables:
sexes = data.table( sex=c(0, 1 ), sex_code=c("m", "f") )
mats =  data.table( mat=c(0, 1 ), maturity=c("imm", "mat") )

out = sexes[out, on="sex"]
out$sex = NULL
out$mat = NULL

fwrite(out, file="~/tmp/snowcrab_4Vn.csv")
