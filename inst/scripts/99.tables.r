---
title: "Snow crab tables"
author: "Jae S. Choi"
toc: true
number-sections: true
highlight-style: pygments
editor:
  render-on-save: true
format:
  html: 
    code-fold: true
    html-math-method: katex
  pdf:
    pdf-engine: lualatex
    geometry: 
      - top=30mm
      - left=30mm
  docx: default 
---


<!--Copy this file to a work directory (e.g., bio.data/bio.snowcrab/assessments/) and Quarto from run from there...-->

# Set up environment

First set up environment. Data comes from:

- At sea observations of fishery (ISSDB)
- Dockside monitoring (Marfis)
- Snow crab survey 

These are mostly imported and formatted in [01_snowcrab.R](01_snowcrab.R)


```{r}
#| eval: true
#| output: false

    # Get data and format based upon parameters:
    year.assessment = 2023
    p = bio.snowcrab::load.environment( year.assessment=year.assessment )
    
    # loadfunctions( "aegis")
    # loadfunctions( "bio.snowcrab")  # in case of local edits

    # require(ggplot2)
    # require(data.table)
 
    require(gt)  # table formatting
 
    outtabledir = file.path( p$annual.results, "tables" )
    
    years = as.character(1996: year.assessment)

    regions = c("cfanorth", "cfasouth", "cfa4x")
    nregions = length(regions)
  
    FD = fishery_data()  # mass in tonnes
    fda = FD$summary_annual
    dt = as.data.frame( fda[ which(fda$yr %in% c(year.assessment - c(0:10))),] )
    dt =  dt[,c("region", "yr", "Licenses", "TAC", "landings", "effort", "cpue")] 
    names(dt) = c("Region", "Year", "Licenses", "TAC", "Landings", "Effort", "CPUE") 
    rownames(dt) = NULL
  
    odb0 = setDT(observer.db("odb"))
    odb0$region = NA

    for ( reg in regions) {
      r = polygon_inside(x = odb0, region = aegis.polygons::polygon_internal_code(reg), planar=FALSE)
      odb0$region[r] = reg
    }


```


## Fishery statistics from at sea observations

NENS:

```{r}
#| eval: true
#| output: false
#| label: table-fishery-nens-perf
#| tbl-cap: "Fishery performance statistics in NENS. Units are: TACs and Landings (tons, t), Effort ($\\times 10^3$ trap hauls, th) and CPUE (kg/th)."

ii = which(dt$Region=="cfanorth")
oo = dt[ii, c("Year", "Licenses", "TAC", "Landings", "Effort", "CPUE")] 
gt::gt(oo) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
```


SENS:

```{r}
#| eval: true
#| output: false
#| label: table-fishery-sens-perf
#| tbl-cap: "Fishery performance statistics in SENS. Units are: TACs and Landings (tons, t), Effort ($\\times 10^3$ trap hauls, th) and CPUE (kg/th)."
ii = which(dt$Region=="cfasouth")
oo = dt[ii, c("Year", "Licenses", "TAC", "Landings", "Effort", "CPUE")] 
gt::gt(oo) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
```

4X:

```{r}
#| eval: true
#| output: false
#| label: table-fishery-4x-perf
#| tbl-cap: "Fishery performance statistics in 4X. Units are: TACs and Landings (tons, t), Effort ($\\times 10^3$ trap hauls, th) and CPUE (kg/th). There were no landings or TACs in 2018/2019 due to indications of low abundance. The 2022 season is ongoing."
ii = which(dt$Region=="cfa4x")
oo = dt[ii, c("Year", "Licenses", "TAC", "Landings", "Effort", "CPUE")] 
gt::gt(oo) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
```
 
## At sea observed data

### Carapace condition from observed data  < 95mm CW

```{r}
#| eval: true
#| output: false
  odb = odb0[ cw < 95 & prodcd_id==0 & shell %in% c(1:5) & region %in% regions & sex==0, ]  # male
```

NENS:

```{r}
#| eval: true
#| output: false
#| label: table-fishery-nens-sublegal
#| tbl-cap: "Fishery performance statistics in NENS. Distribution of at sea observations of males less than 95 mm CW by year and shell condition."
resN = dcast( odb[ region=="cfanorth", .(N=.N), by=.(fishyr, shell) ], fishyr  ~ shell, value.var="N", fill=0, drop=FALSE, na.rm=TRUE )
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
#| output: false
#| label: table-fishery-sens-sublegal
#| tbl-cap: "Fishery performance statistics in SENS. Distribution of at sea observations of males less than 95 mm CW by year and shell condition."
resS = dcast( odb[ region=="cfasouth", .(N=.N), by=.(fishyr, shell) ], fishyr  ~ shell, value.var="N", fill=0, drop=FALSE, na.rm=TRUE )
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
#| output: false
#| label: table-fishery-4x-sublegal
#| tbl-cap: "Fishery performance statistics in 4X. Distribution of at sea observations of males less than 95 mm CW by year and shell condition."
resX = dcast( odb[ region=="cfa4x", .(N=.N), by=.(fishyr, shell) ], fishyr  ~ shell, value.var="N", fill=0, drop=FALSE, na.rm=TRUE )
names(resX) = c("Year", "CC1", "CC2", "CC3", "CC4", "CC5" )
resX$Total = rowSums( resX[, 2:6 ], na.rm=TRUE)
resX[, 2:6 ] = round(resX[, 2:6 ] / resX$Total * 100, digits=2)
gt::gt(resS) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
```


### Carapace condition from observed data  >= 95mm CW

```{r}
#| eval: true
#| output: false
odb = odb0[ cw >= 95 & cw < 170  & prodcd_id==0 & shell %in% c(1:5) & region %in% regions & sex==0, ]  # male
```

NENS:

```{r}
#| eval: true
#| output: false
#| label: table-fishery-nens-comm
#| tbl-cap: "Fishery performance statistics in NENS. Distribution of at sea observations of males less than 95 mm CW by year and shell condition."
resN = dcast( odb[ region=="cfanorth", .(N=.N), by=.(fishyr, shell) ], fishyr  ~ shell, value.var="N", fill=0, drop=FALSE, na.rm=TRUE )
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
#| output: false
#| label: table-fishery-sens-comm
#| tbl-cap: "Fishery performance statistics in SENS. Distribution of at sea observations of males less than 95 mm CW by year and shell condition."
resS = dcast( odb[ region=="cfasouth", .(N=.N), by=.(fishyr, shell) ], fishyr  ~ shell, value.var="N", fill=0, drop=FALSE, na.rm=TRUE )
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
#| output: false
#| label: table-fishery-4x-comm
#| tbl-cap: "Fishery performance statistics in 4X. Distribution of at sea observations of males less than 95 mm CW by year and shell condition."
resX = dcast( odb[ region=="cfa4x", .(N=.N), by=.(fishyr, shell) ], fishyr  ~ shell, value.var="N", fill=0, drop=FALSE, na.rm=TRUE )
names(resX) = c("Year", "CC1", "CC2", "CC3", "CC4", "CC5" )
resX$Total = rowSums( resX[, 2:6 ], na.rm=TRUE)
resX[, 2:6 ] = round(resX[, 2:6 ] / resX$Total * 100, digits=2)
gt::gt(resS) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
```

   
   
###  Percent soft from observed data

There are two possible definitions: 

- carapace conditions 1 and 2
- durometer < 68


```{r}
#| eval: true
#| output: false
odb = odb0[ cw >= 95 & cw < 170  & prodcd_id==0 & shell %in% c(1:5) & region %in% regions & sex==0, ]  # male
scn = setDT(FD$shell_condition)

```

NENS (durometer-based):

```{r}
#| eval: true
#| output: false
#| label: table-fishery-nens-soft-durometer
#| tbl-cap: "Fishery performance statistics in NENS. Distribution of at sea observations of males soft-shelled (durometer)."
softN  = odb[ region=="cfanorth" & durometer <  68, .(Soft=.N), by=.(fishyr ) ] 
totalN = odb[ region=="cfanorth", .(Total=.N), by=.(fishyr) ] 
resN = softN[totalN, on="fishyr"]
resN = resN[, .(Year=fishyr, Soft=round(Soft/Total*100,2), Total=Total) ]  
gt::gt(resN) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
```

SENS (durometer-based):

```{r}
#| eval: true
#| output: false
#| label: table-fishery-sens-soft-durometer
#| tbl-cap: "Fishery performance statistics in SENS. Distribution of at sea observations of males soft-shelled (durometer)."
softS  = odb[ region=="cfanorth" & durometer <  68, .(Soft=.N), by=.(fishyr ) ] 
totalS = odb[ region=="cfanorth", .(Total=.N), by=.(fishyr) ] 
resS = softS[totalS, on="fishyr"]
resS = resS[, .(Year=fishyr, Soft=round(Soft/Total*100,2), Total=Total) ]  
gt::gt(resS) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
```

4X (durometer-based):

```{r}
#| eval: true
#| output: false
#| label: table-fishery-4x-soft-durometer
#| tbl-cap: "Fishery performance statistics in 4X. Distribution of at sea observations of males soft-shelled (durometer)."
softX  = odb[ region=="cfanorth" & durometer <  68, .(Soft=.N), by=.(fishyr ) ] 
totalX = odb[ region=="cfanorth", .(Total=.N), by=.(fishyr) ] 
resX = softX[totalX, on="fishyr"]
resX = resX[, .(Year=fishyr, Soft=round(Soft/Total*100,2), Total=Total) ]  
gt::gt(resX) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
```

   
   
NENS (shell-condition-based):

```{r}
#| eval: true
#| output: false
#| label: table-fishery-nens-soft-shell
#| tbl-cap: "Fishery performance statistics in NENS. Distribution of at sea observations of males soft-shelled (durometer)."
resN = scn[ region=="cfanorth" & shell %in% c(1,2), .(Soft=sum(percent), Total=unique(total)[1]), by=.(fishyr)]
gt::gt(resN) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
```

SENS (shell-condition-based):

```{r}
#| eval: true
#| output: false
#| label: table-fishery-sens-soft-shell
#| tbl-cap: "Fishery performance statistics in SENS. Distribution of at sea observations of males soft-shelled (durometer)."
resS = scn[ region=="cfasouth" & shell %in% c(1,2), .(Soft=sum(percent), Total=unique(total)[1]), by=.(fishyr)]
gt::gt(resS) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
```

4X (shell-condition-based):

```{r}
#| eval: true
#| output: false
#| label: table-fishery-4x-soft-shell
#| tbl-cap: "Fishery performance statistics in 4X. Distribution of at sea observations of males soft-shelled (durometer)."
resX = scn[ region=="cfa4x" & shell %in% c(1,2), .(Soft=sum(percent), Total=unique(total)[1]), by=.(fishyr)]
gt::gt(resX) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
```

   
  

<!--
# instars of interest: 11 and 12

# growth increment (assumming average weight in the midpoint of each increment)
growth.11.to.12 =  predict.mass.g.from.CW.mm( mean(CW.interval.male(12)) ) - predict.mass.g.from.CW.mm (mean(CW.interval.male(11)) )
(growth.11.to.12)
# = 419 g
#  12to13 = ~450
-->

  region="cfaall"
  o = observer.db( DS="bycatch_summary", p=p,  yrs=p$yrs, region=region )   
```

Compare discard rates Maritimes:

```{r}
#| fig-dpi: 144
#| fig-height: 4
  pl = ggplot( o$eff_summ, aes(x=fishyr, y=loss, ymin=loss-losssd, ymax=loss+losssd) ) +
 
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




```

 


## Naive estimation: directly from observations


deprecated = TRUE
if (deprecated) {
  # the following are deprecated methods as of 2024 (JC), here only for reference .. most functionality is nowin the Quarto/Rmarkdown above. 
  

  # Tables obtained after completion of data assimilation and processing up to the end of "01.snowcrab.r"
 
  year.assessment = 2023

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

  # regions = list( region=c("cfanorth", "cfasouth", "cfa4x"))
  regions = list(subarea=c("cfanorth", "cfa23", "cfa24", "cfa4x"))
  res = fishery_data( toget="summary_annual", regions=regions )
  
  l = NULL
  for (r in unlist(regions)) {
    #round the landings to ton
    #round CPUE to no decimal places
    #round the effort to per x1000 trap hauls
    print(r)
    vn = names(regions)
    print(res[res[[vn]]==r, ])
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
 