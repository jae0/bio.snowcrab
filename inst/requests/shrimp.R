
  # ----------------------------------------------------------
  # Shrimp assessment uses this as a recruitment index of snow crab
  # -- N immature male > 56 mm CW in areas 23ab and 24ab
  #  == R2 + R3 + R4

  require(aegis)
  require(data.table)
  require(terra)
  require(sf)
  
  p = bio.snowcrab::load.environment()

  set = snowcrab.db("set.biologicals")

  p$regions.to.model = "cfa.23ab.24ab"  # shrimp area of interest
  p$vars.to.model = "pre.recruit.no"    # size/sex fraction of interest
  p$yrs = p$yrs[ which(p$yrs>1997) ]

  # new method: directly computed averages of core areas
  i = polygon_inside(x=set[, c("plon", "plat")], 
    region= p$regions.to.model, planar=TRUE, proj.type=p$aegis_proj4string_planar_km )

  xs = data.table(set[ i, ])

  out = xs[, .(
    mean = mean(pre.recruit.no, na.rm=TRUE), 
    sd = sd(pre.recruit.no, na.rm=TRUE), 
    n = .N) , 
    by=yr ]
 
  out$se = out$sd/ sqrt(out$n-1)
  out = out[order(yr),]

  outdir  = file.path("~", "tmp")
  dir.create(outdir)
  write.csv ( out, file=file.path(outdir, "bio.snowcrab.recruitment.index.csv"))

  # plot(  out[, c("mean")], ylim=c(0, max(out$mean, na.rm=T)*1.1))
