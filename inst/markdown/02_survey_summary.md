---
title: "Snow crab trawl survey -- figures and tables"
author:
  - name: Snow crab group
    # orcid: 0000-0003-3632-5723 
    # email: jae.choi@dfo-mpo.gc.ca
    # email: choi.jae.seok@gmail.com
    # corresponding: true
    affiliation: 
      - name: Bedford Institute of Oceanography, Fisheries and Oceans Canada
        city: Dartmouth
        state: NS
        url: www.bio.gc.ca
date: 2024-08-17
keywords: 
  - snow crab trawl survey results
  - basic tables and figures
abstract: |
  Snow crab demographic structure and environmental conditions. 
toc: true
number-sections: true
highlight-style: pygments
bibliography: media/references.bib  
# csl: media/canadian-journal-of-fisheries-and-aquatic-sciences.csl  # see https://www.zotero.org/styles for more
license: "CC BY"
copyright: 
  holder: Jae S. Choi
  year: 2024
citation: 
  container-title: https://github.com/jae0/bio.snowcrab/
  doi: NA
funding: "The snow crab scientific survey was funded by the snow crab fishers of Maritimes Region of Atlantic Canada."
editor:
  render-on-save: false
format:
  html: 
    code-fold: true
    html-math-method: katex
    embed-resources: true
  pdf:
    pdf-engine: lualatex
  docx: default 
  beamer:
    pdf-engine: lualatex
---



# Snow crab trawl survey -- figures and tables

## Preamble

This is a markdown document. It can be viewed in formatted form via:
  
  - a web browser open the webpage: [02_survey_summary.md](https://github.com/jae0/bio.snowcrab/tree/master/inst/markdown/02_survey_summary.md)   

  - a web browser open the local file directly: [02_survey_summary.md](../markdown/02_survey_summary.md) (you might need to install a browser add-in), or 
  
  - a text editor with markdown awareness (e.g. Rstudio, VSCode, etc. ).


As this document uses the Quarto and Rmarkdown dialect of Markdown, you can  create reports, by running one of the following in a shell command prompt or using Rstudio. HTML is the most versatile at the moment: 


```shell
 
# {via Quarto}
make quarto FN=02_survey_summary YR=2024 SOURCE=~/bio/bio.snowcrab/inst/markdown WK=~/bio.data/bio.snowcrab/assessments DOCEXTENSION=html 

# {via Rmarkdown}
make rmarkdown FN=02_survey_summary YR=2024 SOURCE=~/bio/bio.snowcrab/inst/markdown WK=~/bio.data/bio.snowcrab/assessments  DOCTYPE=pdf_document DOCEXTENSION=pdf 

# {via pandoc}
make pdf FN=02_survey_summary  

```

Or, see the Makefile and alter defaults to your needs. Note as it tries to be viewable in all dialects, formatting is not perfect. 
  
Note: Quarto method does not pass params easily. So you must adjust "params" in yaml at the top of this fle, or use to quarto command such as:

```shell
quarto ... -P year.assessment:$(YR) -P media_loc:$(MEDIA) 
```

[See here for more YAML options.](https://quarto.org/docs/output-formats/all-formats.html)
 
  

## Set up environment

  - Ensure the **year.assessment** is correct

  - Ensure that the data pulls of survey data stored on the ISSDB is complete and assimilated to the end of [01_snowcrab_data.md](https://github.com/jae0/bio.snowcrab/tree/master/inst/markdown/01_snowcrab-data.md).
  
 


```{r}
#| eval: true
#| output: false
#| label: setup

require(aegis)

# Get data and format based upon parameters:
year.assessment = 2024
p = bio.snowcrab::load.environment( year.assessment=year.assessment )
yrs = 1996:year.assessment # redo all years

# loadfunctions( "aegis")
# loadfunctions( "bio.snowcrab")  # in case of local edits

SCD = project.datadirectory("bio.snowcrab")

media_loc = project.codedirectory("bio.snowcrab", "inst", "markdown", "media")

sn_env = snowcrab_load_key_results_to_memory( year.assessment, debugging=FALSE, return_as_list=TRUE  ) 

attach(sn_env)
  s
```


## Overview of locations and counts


### Map: Survey locations
```{r}
#| eval: true
#| output: false
#| label: map-survey-locations
#| fig-cap: "Survey station locations"

map.survey.locations( p=p, basedir=file.path(p$project.outputdir, "maps", "survey.locations"),  newyear=F, map.method="lattice"  )

# map.survey.locations( p, basedir=file.path(p$project.outputdir, "maps", "survey.locations"),  newyear=F, map.method="googleearth"  )
```


### Counts of stations in each area

```{r}
#| eval: true
#| output: true
#| label: table-survey-station-count
#| tbl-cap: "Survey station counts"

set = snowcrab.db(p=p, DS="set.clean")
setDT(set)
# check towquality .. this should always == 1
if (length( unique( set$towquality) ) != 1 ) print("error -- not good tows")
set$region = NA
for (reg in c( "cfanorth", "cfasouth", "cfa4x"  ) ) {
  d = polygon_inside(set[,c("lon","lat")], reg)
  set$region[d] = reg 
}
out = dcast( set[, .(N=.N), by=.(region, yr)], yr~region, value.var="N", fill=0, drop=FALSE, na.rm=TRUE )
out[,Total:=sum(cfanorth,cfasouth,cfa4x, na.rm=TRUE)]
out = out[, .(yr, cfanorth, cfasouth, cfa4x)]
names(out) = c("Year", "NENS", "SENS", "4X")
gt::gt(out) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
  summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
  footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
  row_group.padding = gt::px(1))
```



```r
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
```


 

## Carapace condition from trawl data \>= 95mm CW

```{r}
#| eval: true
#| output: false

det = snowcrab.db( p=p, DS="det.georeferenced" )
setDT(det)
det$fishyr = det$yr  ## the counting routine expectes this variable
det = det[ cw >= 95 ,]  # commerical sized crab only
years = sort( unique( det$yr ) )
det$region = NA
for ( reg in regions) {
  r = polygon_inside(x = det, region = aegis.polygons::polygon_internal_code(reg), planar=FALSE)
  det$region[r] = reg
}

```

NENS:

```{r}
#| eval: true
#| output: true
#| label: table-survey-nens-comm
#| tbl-cap: "Distribution of NENS survey: males less than 95 mm CW by year and shell condition."

resN = dcast( det[ region=="cfanorth" & !is.na(shell), .(N=.N), by=.(fishyr, shell) ], fishyr  ~ shell, value.var="N", fill=0, drop=FALSE, na.rm=TRUE )
names(resN) = c("Year", "CC1", "CC2", "CC3", "CC4", "CC5" )
resN$Total = rowSums( resN[, 2:6 ], na.rm=TRUE)
resN[, 2:6 ] = round(resN[, 2:6 ] / resN$Total * 100, digits=2)
gt::gt(resN) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
  summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
  footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
  row_group.padding = gt::px(1))
```

SENS:

```{r}
#| eval: true
#| output: true
#| label: table-survey-sens-comm
#| tbl-cap: "Distribution of SENS survey: males less than 95 mm CW by year and shell condition."
resS = dcast( det[ region=="cfasouth" & !is.na(shell), .(N=.N), by=.(fishyr, shell) ], fishyr  ~ shell, value.var="N", fill=0, drop=FALSE, na.rm=TRUE )
names(resS) = c("Year", "CC1", "CC2", "CC3", "CC4", "CC5" )
resS$Total = rowSums( resS[, 2:6 ], na.rm=TRUE)
resS[, 2:6 ] = round(resS[, 2:6 ] / resS$Total * 100, digits=2)
gt::gt(resS) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
  summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
  footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
  row_group.padding = gt::px(1))
```

4X:

```{r}
#| eval: true
#| output: true
#| label: table-survey-4X-comm
#| tbl-cap: "Distribution of 4X survey: males less than 95 mm CW by year and shell condition."

resX = dcast( det[ region=="cfa4x" & !is.na(shell), .(N=.N), by=.(fishyr, shell) ], fishyr  ~ shell, value.var="N", fill=0, drop=FALSE, na.rm=TRUE )
names(resX) = c("Year", "CC1", "CC2", "CC3", "CC4", "CC5" )
resX$Total = rowSums( resX[, 2:6 ], na.rm=TRUE)
resX[, 2:6 ] = round(resX[, 2:6 ] / resX$Total * 100, digits=2)
gt::gt(resX) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
  summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
  footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
  row_group.padding = gt::px(1))
```

 

## Size-frequency distributions of snow crab cw from trawl data, broken down by maturity classes


```r
# take subset in years
years = as.character( c(-9:0) + year.assessment )
# years = as.character(2004:2013)
# years = as.character(1996:2003)


regions=c("cfanorth", "cfasouth", "cfa4x")
# outdir=file.path( p$annual.results, "figures", "size.freq", "survey_1996_2003" )
# outdir=file.path( p$annual.results, "figures", "size.freq", "survey_2004_2013" )
outdir=file.path( p$annual.results, "figures", "size.freq", "survey" )



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


```
 



## Timeseries of all survey variables

```r
  figure.timeseries.survey(p=p, outdir=file.path(p$annual.results, "timeseries", "survey"), variables="R0.mass", plotyears=2004:p$year.assessment) # just R0 to see
  
  figure.timeseries.survey(p=p, outdir=file.path(p$annual.results, "timeseries", "survey"), variables=c("sexratio.all","sexratio.mat","sexratio.imm"))

  figure.timeseries.survey(p=p, outdir=file.path(p$annual.results, "timeseries", "survey"), variables="cw.male.mat.mean", plotyears=2004:p$year.assessment, backtransform=TRUE) 
  
  figure.timeseries.survey(p=p, outdir=file.path(p$annual.results, "timeseries", "survey"), plotyears=2004:p$year.assessment) # all variables
  figure.timeseries.survey(p=p, outdir=file.path(p$annual.results, "timeseries", "observer"),plotyears=2004:p$year.assessment,type='observer')

  # no longer relevant (incomplete) as temp is now created through temp db. and not in gshyd
  # figure.timeseries.survey(p=p, outdir=file.path(p$annual.results, "timeseries", "survey"), plotyears=2004:p$year.assessment,type='groundfish.t') # groundfish survey temperature
  
  # area-specific figures
  figure_area_based_extraction_from_carstm(DS="temperature", year.assessment )  # can only do done once we have an sppoly for snow crab
```



## Trawl survey catch of other species

### Predators 

cod, haddock, halibut, plaice, wolfish, thornyskate, smoothskate, winterskate

```r
  species = c(10, 11, 30, 40, 50, 201, 202, 204 )

  # by mean number
  figure.timeseries.bycatch(p=p, species=species, plotyears=2004:p$year.assessment, outdir=file.path(p$annual.results,"timeseries", "survey"), type="no" )

  # by mean mass
  figure.timeseries.bycatch(p=p, species=species, plotyears=2004:p$year.assessment, outdir=file.path(p$annual.results,"timeseries", "survey"), type="mass" )

```


### Competitors

northernshrimp, jonahcrab, lessertoadcrab

```r
  species = c( 2521, 2511, 2211)

  # by mean number
  figure.timeseries.bycatch(p=p, species=species, plotyears=2004:p$year.assessment, outdir=file.path(p$annual.results,"timeseries", "survey"), type="no" )

  # by mean mass
  figure.timeseries.bycatch(p=p, species=species, plotyears=2004:p$year.assessment, outdir=file.path(p$annual.results,"timeseries", "survey"), type="mass" )

```



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





  ########################################################
  ########################  Retired figures ###################



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

  res = maturity_region_year(p)  # timeseries of maturity


## References

<!--

deprecated = TRUE
if (deprecated) {
  # the following are deprecated methods as of 2024 (JC), here only for reference .. most functionality is nowin the Quarto/Rmarkdown above. 
  

  # Tables obtained after completion of data assimilation and processing up to the end of "01.snowcrab.r"
 
  year.assessment = 2024

  p = bio.snowcrab::load.environment( year.assessment=year.assessment )
 

  require(gridExtra)
  library("xtable")
  library("R2HTML")

  odb0 = observer.db("odb")
  regions = c("cfanorth", "cfasouth", "cfa4x")
  nregions = length(regions)

  #------------------------------------------------
  #Fisheries statistics per region
  tabledir = file.path(project.datadirectory("bio.snowcrab"), "data", "fisheries")
  outtabledir= file.path(project.datadirectory("bio.snowcrab"), "assessments", p$year.assessment, "tables", "logbook")
  if(!dir.exists(tabledir)) dir.create(tabledir, recursive =T)
  if(!dir.exists(outtabledir)) dir.create(outtabledir, recursive =T)

  setwd(tabledir)

  NFS <- xtable(read.csv("NENS_FisherySummary.csv"))
  SFS <- xtable(read.csv("SENS_FisherySummary.csv"))
  Fx <- xtable(read.csv("4x_FisherySummary.csv"))

  setwd(outtabledir)
  print.xtable(NFS, type="latex", file="NENS_FisherySummary.tex")
  print.xtable(NFS, type="html", file="NENS_FisherySummary.html")

  print.xtable(SFS, type="latex", file="SENS_FisherySummary.tex")
  print.xtable(SFS, type="html", file="SENS_FisherySummary.html")

  print.xtable(Fx, type="latex", file="4x_FisherySummary.tex")
  print.xtable(Fx, type="html", file="4x_FisherySummary.html")

  #regions = c("cfaall")
  #regions = c("cfanorth", "cfasouth", "cfa4x")
  regions = c("cfanorth", "cfa23", "cfa24", "cfa4x")
  l = NULL
  for (r in regions) {
    res = get.fishery.stats.by.region( Reg=r) #need to add the TACs per year and number of licences
    #round the landings to ton
    #round CPUE to no decimal places
    #round the effort to per x1000 trap hauls
    print(r)
    print(res)
  }


  # ----------------------------------------
  #  Carapace condition from observed data  < 95mm CW

  outtabledir= file.path(project.datadirectory("bio.snowcrab"), "assessments", p$year.assessment, "tables", "observer")
  dir.create(outtabledir, recursive=TRUE)

  odb = odb0
  odb = odb[ which( odb$cw < 95 & odb$prodcd_id=="0" ) ,]
  regions = c("cfanorth", "cfasouth", "cfa4x")
  nregions = length(regions)
  years = sort( unique( odb$fishyr ) )

  res = NULL
  for (r in p$regions) {
    for (y in years) {
      out = proportion.cc (odb, region=r, year=y)
      res = rbind( res, cbind( r, y, t(out)) )
    }}

  cnames = c("region", "fishyr", c(1:5), "ntot")
  colnames(res) = cnames
  print(res)
  res = as.data.frame(res)
  res[is.na(res)] <- NA
  ct <- c("CC1", "CC2", "CC3", "CC4", "CC5", "Total")

  setwd(outtabledir)

  Rn= res[res$region=="cfanorth", 3:8]
  print(Rn)
  rownames(Rn) = years
  colnames(Rn) = ct
  print.xtable(Rn, type="latex", file="table.CC.small.north.obs.tex")
  HTML(Rn, file="table.CC.Small.north.obs.html")

  Rs= res[res$region=="cfasouth", 3:8]
  rownames(Rs) = years
  colnames(Rs) = ct
  print.xtable(Rs, type="latex", file="table.CC.small.south.obs.tex")
  HTML(Rs, file="table.CC.small.south.obs.html")

  Rx= res[res$region=="cfa4x", 3:8]
  rownames(Rx) = years
  colnames(Rx) = ct
  print.xtable(Rs, type="latex", file="table.CC.small.4x.obs.tex")
  HTML(Rx, file="table.CC.small.4x.obs.html")


  # ----------------------------------------
  #  Carapace condition from observed data >=95mm CW
  odb = odb0
  odb = odb[ which( odb$cw >= 95 & odb$cw < 170 & odb$prodcd_id=="0" ) ,]  # commerical sized crab only
  years = sort( unique( odb$fishyr ) )

  # get proportion by cc
  regions = c("cfanorth", "cfasouth", "cfa4x")
  years = sort( unique( odb$fishyr ) )

  res = NULL
  for (r in regions) {
    for (y in years) {
      out = proportion.cc (odb, region=r, year=y)
      res = rbind( res, cbind( r, y, t(out)) )
    }}

  cnames = c("region", "fishyr", c(1:5), "ntot")
  colnames(res) = cnames
  print(res)
  res = as.data.frame(res)
  res[is.na(res)] <- NA
  #  for (i in cnames[-1]) res[,i] = as.numeric(as.character((res[,i])))
  setwd(outtabledir)

  ct <- c("CC1", "CC2", "CC3", "CC4", "CC5")
  Rn = res[res$region=="cfanorth", 3:7]
  #Rn = as.matrix( res[ which(res$region=="cfanorth") , as.character(c(1:5)) ] )
  rownames(Rn) = years
  colnames(Rn) = ct
  print.xtable(Rn, type="latex", file="table.CC.large.north.obs.tex")
  HTML(Rn, file="table.CC.large.north.obs.html")

  #Rs = as.matrix( res[ which(res$region=="cfasouth") , as.character(c(1:5)) ] )
  Rs = res[res$region=="cfasouth", 3:7]
  rownames(Rs) = years
  colnames(Rs) = ct
  print.xtable(Rs, type="latex", file="table.CC.large.south.obs.tex")
  HTML(Rs, file="table.CC.large.south.obs.html")

  #Rx = as.matrix( res[ which(res$region=="cfa4x") , as.character(c(1:5)) ] )
  Rx = res[res$region=="cfa4x", 3:7]
  rownames(Rx) = years
  colnames(Rx) = ct
  print.xtable(Rx, type="latex", file="table.CC.large.4x.obs.tex")
  HTML(Rx, file="table.CC.large.4x.obs.html")

  # ----------------------------------------
  #  Percent soft from observed data

  odb = odb0
  odb = odb[ which( odb$cw > 95 & odb$cw < 170 & odb$prodcd_id=="0" ) ,]  # commercial crab
  years = sort( unique( odb$fishyr ) )


  res = NULL
  for (r in p$regions) {
    for (y in years) {
      out = proportion.soft (odb, region=r, year=y)
      res = rbind( res, cbind( r, y, t(out)) )
    }}

  cnames = c("region", "fishyr", "pr.soft", "nsoft", "ntot")
  colnames(res) = cnames
  print(res)
  res = as.data.frame(res)

  for (i in cnames[-1]) res[,i] = as.numeric(as.character((res[,i])))

  Rn = as.matrix( res[ which(res$region=="cfanorth") , ] )
  rownames(Rn) = years
  HTML(Rn, file="table.proportion.soft.north.obs.html")

  Rs = as.matrix( res[ which(res$region=="cfasouth" ),   ] )
  rownames(Rs) = years
  HTML(Rs, file="table.proportion.soft.south.obs.html")

  Rx = as.matrix( res[ which(res$region=="cfa4x") , ] )
  rownames(Rx) = years
  HTML(Rx, file="table.proportion.soft.4x.obs.html")



  # instars of interest: 11 and 12

  # growth increment (assumming average weight in the midpoint of each increment)
  growth.11.to.12 =  predict.mass.g.from.CW.mm( mean(CW.interval.male(12)) ) - predict.mass.g.from.CW.mm (mean(CW.interval.male(11)) )

  # = 419 g
  #  12to13 = ~450



  # Table of proportion discarded
  odb = observer.db("odb")
  regions = c("cfanorth", "cfasouth", "cfa4x")
  years = sort( unique( odb$fishyr ) )
  out = NULL
  for (r in regions) {
    for (y in years) {
      res = proportion.legal (odb, region=r, year=y)
      out = rbind(out, cbind( r, y, res[1], res[2], res[3] ) )
    } }
  out

  HTML(out, file="table.proportion.discarded.html")


  # ---------------------------------------- USED
  #  Carapace condition from trawl data  >= 95mm CW  ... not kriged .. simple proportions

  det0 = snowcrab.db( p=p, DS="det.georeferenced" )
  det0$fishyr = det0$yr  ## the counting routine expectes this variable

  det = det0[ which( det0$cw >= 95 ) ,]  # commerical sized crab only
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

  HTML(res, file="table.CC.large.survey.html")

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
 
  HTML(out, file="table.tow_counts_survey.html")



  # % mat calculations: deprecated  ... size analysis is now in 01_snowcrab.R

  # loc = file.path(sc.R, "size.data")
  # dir.create(path=loc, recursive=T, showWarnings=F)
  # outfilename = paste( c("mi", "mm", "fi", "fm"), "rdata", sep=".")
  # outfile = file.path(loc, paste(outfilename))
  # for (f in  outfile) load(f)


  # f.i = f.imm[which( rownames(f.imm)%in% sids ) ,]
  # f.i.means = apply(X=f.i, MARGIN=2, FUN=mean)
  # f.m = f.mat[which( rownames(f.mat)%in% sids ) ,]
  # f.m.means = apply(X=f.m, MARGIN=2, FUN=mean)

  # toplot = rbind(f.m.means, f.i.means)

  # ii = as.data.frame(t(toplot))
  # ii$cw = as.numeric(rownames(ii))
  # ii$pmat = ii[,1]/ (ii[,1]+ii[,2]) * 100

  # plot(ii$cw, ii$pmat)
  # abline(h=50)

  # str(ii)
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



}
-->
