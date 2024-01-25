
  # Figures   obtained after completion of data assimilation and processing up to the end of "01.snowcrab.r"


  year.assessment = 2023

  p = bio.snowcrab::load.environment( year.assessment=year.assessment )
  
  # ------------------------------------------
  # Time-series: Fisheries landings

  figure.landings.timeseries( yearmax=p$year.assessment, outdir=file.path( p$annual.results,  "timeseries","fishery"), outfile="landings.ts", outfile2="landings.ts.sm" )

  # ------------------------------------------
  # Time-series: Fisheries effort
  figure.effort.timeseries( yearmax=p$year.assessment, outdir=file.path( p$annual.results,"timeseries", "fishery"), outfile="effort.ts", outfile2="effort.ts.sm" )

  # ------------------------------------------
  # Time-series: Fisheries CPUE
  figure.cpue.timeseries( yearmax=p$year.assessment, outdir=file.path( p$annual.results,"timeseries", "fishery"), outfile="cpue.ts", outfile2="cpue.sm.ts" )

  # ------------------------------------------
  # Size frequency distributions, broken down by moult category from at-sea observed data
  figure.observed.size.freq( regions = c("cfanorth", "cfasouth", "cfa4x"), years="all", outdir=file.path( p$annual.results, "figures", "size.freq", "observer")  )


  figure.sizefreq.carapacecondition( X=snowcrab.db( p=p, DS="det.georeferenced" ), cwbr=4, regions=c("cfanorth", "cfasouth", "cfa4x"), 
    outdir=file.path( p$annual.results, "figures", "size.freq", "carapacecondition" )  ) 

 
  # ------------------------------------------
  # Size-frequency distributions of snow crab cw from trawl data, broken down by maturity classes
  # take subset in years
  years = as.character( c(-9:0) + year.assessment )
  regions=c("cfanorth", "cfasouth", "cfa4x")
  outdir=file.path( p$annual.results, "figures", "size.freq", "survey")
  
  M = size_distributions(p=p, toget="simple_direct", xrange=xrange, dx=dx, Y=years )

  # NOTE :: these produce png files (instead of pdfs) change as required.
  # den=arithmetic mean density, denl = geometric mean density  
  plot_histogram_carapace_width( M=M, years=years, regions=regions, plot_sex="female", yvar="denl", 
    outdir=outdir, cols = c("slategray", "gray95" ) )

  plot_histogram_carapace_width( M=M, years=years, regions=regions, plot_sex="male", yvar="denl", 
    outdir=outdir, cols = c("slategray", "gray95" ) )
  

    if (0) {
      # deprecated methods:
      histograms.size.maturity.update( outdir=file.path( p$annual.results, "figures", "size.freq", "survey"),  redo.data=T )
      histograms.size.maturity.single.area( outdir=file.path( p$annual.results, "figures", "size.freq", "survey"),  area='cfa4x',redo.data=T ) #area = cfanorth, cfasouth of cfa4x

      if (oneoff_2022){
        histograms.size.maturity_oneoff(p=p)
        figure.timeseries.survey_oneoff(p=p,
          outdir=file.path(p$annual.results, "timeseries", "survey", "oneoff"), 
          vlab="R0.mass", variables="totmass", plotyears=2004:p$year.assessment) # just R0 to see
      }
    }


  # ------------------------------------------
  # Timeseries of all survey variables
  figure.timeseries.survey(p=p, outdir=file.path(p$annual.results, "timeseries", "survey"), variables="R0.mass", plotyears=2004:p$year.assessment) # just R0 to see
  
  figure.timeseries.survey(p=p, outdir=file.path(p$annual.results, "timeseries", "survey"), variables=c("sexratio.all","sexratio.mat","sexratio.imm"))

  figure.timeseries.survey(p=p, outdir=file.path(p$annual.results, "timeseries", "survey"), variables="cw.male.mat.mean", plotyears=2004:p$year.assessment, backtransform=TRUE) 
  
  figure.timeseries.survey(p=p, outdir=file.path(p$annual.results, "timeseries", "survey"), plotyears=2004:p$year.assessment) # all variables
  figure.timeseries.survey(p=p, outdir=file.path(p$annual.results, "timeseries", "observer"),plotyears=2004:p$year.assessment,type='observer')

  # no longer relevant (incomplete) as temp is now created through temp db. and not in gshyd
  # figure.timeseries.survey(p=p, outdir=file.path(p$annual.results, "timeseries", "survey"), plotyears=2004:p$year.assessment,type='groundfish.t') # groundfish survey temperature
  
  # area-specific figures
  figure_area_based_extraction_from_carstm(DS="temperature", year.assessment )  # can only do done once we have an sppoly for snow crab

  #-----------------------------------------------

  #Timeseries: geometric mean biomass of by-catch from snow crab survey

  # predators and competitors
  #cod, haddock, halibut, plaice, wolfish, thornyskate, smoothskate, winterskate, northernshrimp, jonahcrab, lessertoadcrab
  species = c(10, 11, 30, 40, 201, 50, 2521, 2511, 202, 204, 2211)
  figure.timeseries.bycatch(p=p, species=species, plotyears=2004:p$year.assessment, outdir=file.path(p$annual.results,"timeseries", "survey"), type="no" )

  figure.timeseries.bycatch(p=p, species=species, plotyears=2004:p$year.assessment, outdir=file.path(p$annual.results,"timeseries", "survey"), type="mass" )


  # Naive estimation of by catch: directly from observations
  # Using at-sea-observed catch rate estimates and re-scaling these 
  # to total snow crab fishery effort: 

  # choose area of interest
  region="cfaall"
  region="cfanorth"
  region="cfasouth"
  region="cfa4x"
  
    o = observer.db( DS="bycatch_summary", p=p,  yrs=p$yrs, region=region )   
    
    # efficiency of snow crab capture vs snow crab discard (soft-shelled, sublegal or female crab):
    dev.new(width=8, height=6)
    pl = ggplot( o$eff_summ, aes(x=fishyr, y=loss, ymin=loss-losssd, ymax=loss+losssd) ) +
      geom_pointrange()  + # Vertical line with point in the middle
      geom_errorbar(width = 0.1, col="brown") + # Standard error bars
      geom_point(size = 1.5, col="darkred") +
      labs(x="Year", y="Discard rate of snow crab (At sea observed, by weight)" )
    (pl)

    # compare catch rates
    dev.new(width=10, height=length(o$spec)/2 )  
    plot( o$spec ~ o$bct, xlab = "At sea observed catch rate in snow crab fishery (kg/trap)", ylab="Species", type="p", cex=1.5, pch=19, col="darkorange", xlim=c(0, max(o$bct, na.rm=TRUE)*1.25), yaxt="n" )  
    text( o$bct, o$spec,  labels=o$species, pos=4, srt=0 , cex=1, col="darkslateblue")
    text( max(o$bct, na.rm=TRUE)*0.75, 1, labels=paste( "Snow crab CPUE (At sea obs., mean): ", o$bct_sc, " kg/trap"), col="darkred", cex=1.5 )

    # actual bycatch table .. 
    print( o$bycatch_table )
    # write.csv( o$bycatch_table, file="tmp.csv")  # if you need to save it or export to excel, etc



  # ------------------------------------------
  # Map:  Interpolated mean/geometric mean of various variables in the set data table
  #  p$do.parallel=F
  p$corners = data.frame(plon=c(220, 990), plat=c(4750, 5270) )

  p$mapyears = year.assessment + c(-5:0 )

  outdir = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual" )


  # just for the roadshow
  map.set.information( p=p, outdir=outdir, variables=c('totmass.male.com', 'totmass.female.mat'),mapyears=p$mapyears)

  map.set.information( p=p, outdir=outdir, variables=c('R0.mass'),mapyears=p$mapyears)

  map.set.information( p=p, variables='t',mapyears=p$mapyears,outdir=outdir,log.variable=F,add.zeros=F,theta=100)

  # bycatch (geometric means)
  # predators and competitors
  # cod, haddock, halibut, plaice, wolfish, thornyskate, smoothskate, winterskate, northernshrimp, jonahcrab, lessertoadcrab
  species = c(10, 11, 30, 40, 201, 50, 2521, 2511, 202, 204, 2211)
  bc.vars = c(paste("ms.mass",species,sep='.'),paste("ms.no",species,sep='.'))
  outdir.bc= file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )
  map.set.information( p, variables=bc.vars, mapyears=p$mapyears, outdir=outdir.bc,probs=c(0,0.975)) #


  # all variables (geometric means)
  #map.set.information( p, outdir=outdir) # takes a long time

  # Means
  # variables that shouldn't be logged
  set = snowcrab.db( p=p, DS="set.biologicals")
  variables = bio.snowcrab::snowcrab.variablelist("all.data")
  variables = intersect( variables, names(set) )

  nolog.variables = c("t","z", "julian", variables[grep("cw",variables)])
  map.set.information( p=p, variables=nolog.variables,outdir=outdir,log.variable=F,add.zeros=F,theta=35)

  # logit transform for ratios
  map.set.information( p=p, variables=c("sexratio.all","sexratio.mat","sexratio.imm"),outdir=outdir,log.variable=F,add.zeros=F,theta=40)

  # Geometric Means
  # all except variables that shouldn't be logged
  mass.vars = variables[!variables%in%nolog.variables][grep('mass',variables[!variables%in%nolog.variables])]
  no.vars = variables[!variables%in%nolog.variables][grep('no',variables[!variables%in%nolog.variables])]
  map.set.information( p=p, variables= mass.vars,outdir=outdir)
  map.set.information( p=p, variables= no.vars,outdir=outdir,probs=c(0,0.975))



  # ------------------------------------------
  # Map: Survey locations

  map.survey.locations( p=p, basedir=file.path(p$project.outputdir, "maps", "survey.locations"),  newyear=F, map.method="lattice"  )
  #  map.survey.locations( p, basedir=file.path(p$project.outputdir, "maps", "survey.locations"),  newyear=F, map.method="googleearth"  )

  # ------------------------------------------
  # Map: Observer locations
  map.observer.locations( p=p, basedir=file.path(p$project.outputdir, "maps","observer.locations" ), newyear=F , map.method="lattice"  )

  # ------------------------------------------
  # Map: Logbook recorded locations
  map.logbook.locations( p=p, basedir=file.path(p$project.outputdir, "maps","logbook.locations" ), newyear=F , map.method="lattice"  )



  # ------------------------------------------
  # Map: Logbook data
  map.fisheries.data( 
    outdir=file.path( p$project.outputdir, "maps", "logbook","snowcrab","annual" ), 
    probs=c(0,0.975),
    plot_crs=st_crs( p$aegis_proj4string_planar_km ),
    outformat="png"
  )




 
  # ------------------
  # counts of stations in each area

  # check towquality .. this should always == 1
  set = snowcrab.db(p=p, DS="set.clean")
  if (length( unique( set$towquality) ) != 1 ) print("error -- not good tows")

  out = data.frame(yr=sort( unique(set$yr )) )
  for (reg in c("cfaall", "cfanorth", "cfasouth","cfa4x"  ) ) {
    d = polygon_inside(set[,c("lon","lat")], reg)
    e = as.data.frame( xtabs(~yr, data=set[d,])  )
    names(e) = c("yr", reg)
    e$yr = as.numeric(as.character(e$yr) )
    out = merge(out, e, by="yr", all=T)
  }
  print(out)

  plot.new()
  year = p$year.assessment
  setdata = set[ which(set$yr==year),]
  # row numbers:
  N = polygon_inside(setdata[,c("lon","lat")], "cfanorth")
  S = polygon_inside(setdata[,c("lon","lat")], "cfasouth")
  X = polygon_inside(setdata[,c("lon","lat")], "cfa4x")
  plot(setdata$lon, setdata$lat)
  points(setdata$lon[N], setdata$lat[N],col="red",pch=20)
  points(setdata$lon[S], setdata$lat[S],col="blue",pch=20)
  points(setdata$lon[X], setdata$lat[X],col="black",pch=20)




# -------------------------------------------------------------------------------------
# deprecated ?

  # % mat calculations:
  #as above but

  loc = file.path(sc.R, "size.data")
  dir.create(path=loc, recursive=T, showWarnings=F)
  outfilename = paste( c("mi", "mm", "fi", "fm"), "rdata", sep=".")
  outfile = file.path(loc, paste(outfilename))
  for (f in  outfile) load(f)


  f.i = f.imm[which( rownames(f.imm)%in% sids ) ,]
  f.i.means = apply(X=f.i, MARGIN=2, FUN=mean)
  f.m = f.mat[which( rownames(f.mat)%in% sids ) ,]
  f.m.means = apply(X=f.m, MARGIN=2, FUN=mean)

  toplot = rbind(f.m.means, f.i.means)

  ii = as.data.frame(t(toplot))
  ii$cw = as.numeric(rownames(ii))
  ii$pmat = ii[,1]/ (ii[,1]+ii[,2]) * 100

  plot(ii$cw, ii$pmat)
  abline(h=50)

  str(ii)
  #`data.frame':   70 obs. of  4 variables:
  # $ f.m.means: num  0 0 0 0 0 ...
  # $ f.i.means: num   2.80  6.19 20.05 24.29 74.11 ...
  # $ cw       : num  12 14 16 18 20 22 24 26 28 30 ...
  # $ pmat     : num  0 0 0 0 0 ...







  # ----------------------------------------   NOT USED ____________
  #  Carapace condition from trawl data  < 95mm CW  ... not kriged .. simple proportions

  det0 = snowcrab.db( p=p, DS="det.georeferenced" )
  det0$fishyr = det0$yr  ## the counting routine expectes this variable

  det = det0[ which( det0$cw < 95 ) ,]  # commerical sized crab only
  years = sort( unique( det$yr ) )

  res = NULL
  for (r in p$regions) {
    for (y in years) {
      out = proportion.cc (det, region=r, year=y)
      res = rbind( res, cbind( r, y, t(out)) )
    }}

  cnames = c("region", "fishyr", c(1:5), "ntot")
  colnames(res) = cnames
  print(res)
  res = as.data.frame(res)

  for (i in cnames[-1]) res[,i] = as.numeric(as.character((res[,i])))
  (res)





  ##########################################################################
  ###############################  Retired figures #########################



  # ------------------------------------------
  # Map: Scotian Shelf with CFA lines and labels  .. using gmt
  # this is the basemap from map.r which is then post-labelled in sodipodi
  #  p$outdir = file.path(p$annual.results,"figures")
  #  p$outfile.basename = file.path(p$outdir, "map.CFAs")
  #  map.basemap.with.cfa.lines( p, conversions=c("ps2png")  )


  # ------------------------------------------
  # Habitat usage comparisons (bivariate) ... requires the full "set.rdata" database and "logbook.dZ.rdata" database
  # habitat.usage( usevar="totno.all", covariate="depth", outdir = file.path(p$annual.results, "habitat.templates") )

  # ------------------------------------------
  # Habitat usage comparisons (bivariate) ... requires the full "set.rdata" database and "logbook.dZ.rdata" database
  #habitat.usage( usevar="totno.all", covariate="temperature", outdir = file.path(p$annual.results, "habitat.templates") )

  # ------------------------------------------
  # Habitat usage comparisons (bivariate) ... requires the full "set.rdata" database and "logbook.dZ.rdata" database
  # habitat.usage( usevar="totno.all", covariate="bottom.slope", outdir = file.path(p$annual.results, "habitat.templates") )

  # ------------------------------------------
  # Habitat usage comparisons (bivariate) ... requires the full "set.rdata" database and "logbook.dZ.rdata" database
  #habitat.usage( usevar="totno.all", covariate="bottom.curvature", outdir = file.path(p$annual.results, "habitat.templates") )

  # ------------------------------------------
  # Habitat usage comparisons (bivariate) ... requires the full "set.rdata" database and "logbook.dZ.rdata" database
  #habitat.usage( usevar="totno.all", covariate="substrate", outdir = file.path(p$annual.results, "habitat.templates") )

  # ------------------------------------------
  # Timeseries: Larval brachyura from the SSIP data
  ##figure.timeseries.larvae( outdir=file.path(p$project.outputdir, "timeseries", "larvae") )

  # ------------------------------------------
  # Growth as a a function of instar for Scotian Shelf snow crab
  figure.growth.instar( outdir=file.path(p$project.outputdir, "growth") )


  # ------------------------------------------
  # Map: Larval distributions from the Scotian Shelf Ichtyoplankton Program data
  map.larvae( p=p, outdir=file.path(p$project.outputdir, "maps", "larvae"), conversions=conversions )


  # ------------------------------------------
  # Map: Spatial representation of maturity patterns of snow crab
  #MG Not sure we use these maps either, check with Adam and Jae
  # map.maturity( p, outdir=file.path(p$project.outputdir, "maps", "maturity"), newyear=T )
