

require(spsUtil)
quietly = spsUtil::quiet

require(flextable)
require(gt)  # table formatting
require(ggplot2)
require(data.table) # for speed
require(lubridate)
require(stringr) 

require(janitor)

require(aegis)  # basic helper tools
require(bio.taxonomy)  # handle species codes
require(bio.snowcrab)

require(spsUtil)
quietly = spsUtil::quiet

require(MBA)

require(aegis)  # basic helper tools

year.assessment = 2025  # change this as appropriate

p = bio.snowcrab::load.environment( year.assessment=year.assessment )  # set up initial settings
p$yrs_observer = c(p$year.assessment + c(-4:0))  # if you change this change yrs_observer 


# compute items of interest
odb = observer.db( DS="odb", p=p ) # 3 minutes
obc = observer.db( DS="bycatch_clean_data", p=p, yrs=p$yrs ) # 3 minutes

# observer data
odb0 = setDT(odb)


# mau is variable name in logbook.db()
mau = "region"  # generic default
mau = "subarea"  # cfanoth, 23, 24, 4x

maus = management_areal_units( mau=mau )  



odb0[[mau]] = NA
for ( reg in maus[["internal"]]) {
    r = polygon_inside(x = odb0, region =reg, planar=FALSE)
    odb0[[mau]][r] = reg
}

odb0[fishyr>2020, length(unique(trip)), by=.(region,fishyr)]


# bycatch summaries
BC = list()
for ( reg in c(maus[["internal"]], "cfaall")) {
    oo = observer.db( DS="bycatch_summary", p=p,  yrs=p$yrs, region=reg  )  # using polygon test
    if (is.null(oo)) { 
    message( "bycatch data is likely broken .. check" ) 
    }
    BC[[reg]] = oo
    BC[[reg]]$bycatch_table_effort[ BC[[reg]]$bycatch_table_effort==0 ] = NA
    BC[[reg]]$bycatch_table_effort[ is.na(BC[[reg]]$bycatch_table_effort) ] = "."
    BC[[reg]]$bycatch_table_catch[ BC[[reg]]$bycatch_table_catch==0 ] = NA
    BC[[reg]]$bycatch_table_catch[ is.na(BC[[reg]]$bycatch_table_catch) ] = "."
    
}


#### At-sea-observed effort: no. trips
# "Number of at-sea observed trips."

oo = dcast( odb0[ fishyr>=2004,.(N=length(unique(tripset))), by=c(mau, "fishyr")], 
  fishyr ~ get(mau), value.var="N", fill=0, drop=FALSE, na.rm=TRUE )
if ( "NA" %in% names(oo) ) oo$"NA" = NULL
keep = c("fishyr", maus[["internal"]])
oo = oo[,..keep]
names(oo) = c("Year", maus[["labels"]] )
oo$Total = rowSums( oo[, 2:maus[["n"]] ], na.rm=TRUE)

gt::gt(oo) |> 
  gt::tab_header( title = md("ASO number of trips")) |> 
  gt::tab_options(table.font.size = 14, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))



#| tbl-cap: "Number of at-sea observed trap hauls."

odb0$th = paste(odb0$tripset, odb0$lat, odb0$lon)  ## <<< NOTE: needs a re-think 
oo = dcast( odb0[ fishyr>=2004,.(N=length(unique(th))), by=c(mau, "fishyr")], 
  fishyr ~ get(mau), value.var="N", fill=0, drop=FALSE, na.rm=TRUE )
if ( "NA" %in% names(oo) ) oo$"NA" = NULL
keep = c("fishyr", maus[["internal"]])
oo = oo[,..keep]
names(oo) = c("Year", maus[["labels"]] )
oo$Total = rowSums( oo[, 2:maus[["n"]] ], na.rm=TRUE)

gt::gt(oo) |> 
  gt::tab_header( title = md("ASO number of (effort) trap hauls")) |> 
  gt::tab_options(table.font.size = 14, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))





#| tbl-cap: "Total weight of at-sea observed crab (kg)." 

oo = dcast( odb0[ fishyr>=2004,.(N=sum( mass, na.rm=TRUE)/1000 ), by=c(mau, "fishyr")], 
  fishyr ~ get(mau), value.var="N", fill=0, drop=FALSE, na.rm=TRUE )
if ( "NA" %in% names(oo) ) oo$"NA" = NULL

keep = c("fishyr", maus[["internal"]])
oo = oo[,..keep]
names(oo) = c("Year", maus[["labels"]] )

oo[, 2:maus[["n"]]] = round( oo[, 2:maus[["n"]]] )
oo$Total = rowSums( oo[, 2:maus[["n"]] ], na.rm=TRUE)

gt::gt(oo) |> 
  gt::tab_header( title = md("ASO Total weight of at-sea observed crab (kg)")) |> 
  gt::tab_options(table.font.size = 14, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))





### remaining is not used (yet) --- ignore





#| label: tbl-observer-softgt95
#| eval: true
#| output: true
#| echo: false
#| results: asis
#| layout-ncol: 1
#| tbl-cap: "Soft-shell incidence. There are two possible definitions of soft-shelled crab: (D) based on durometer measurements < 68 on the hardness scale; and (CC) based upon classification as carapace conditions 1 and 2." 
#| tbl-subcap: 
#|   - "(a)"
#|   - "(b)"
#|   - "(c)"
#|   - "(d)"

odb = odb0[ cw >= 95 & cw < 170  & prodcd_id==0 & shell %in% c(1:5) & get(mau) %in% maus[["internal"]] & sex==0, ]  # male
shell_condition = odb[ !is.na(get(mau)), .N, by=c(mau, "fishyr", "shell") ]
shell_condition[, total:=sum(N, na.rm=TRUE), by=c(mau, "fishyr")]
shell_condition$percent = round(shell_condition$N / shell_condition$total, 3) * 100
shell_condition$Year = shell_condition$fishyr
 
for (r in 1:maus[["n"]]){
  reg = maus[["internal"]][r]
  REG = maus[["labels"]][r]
  cat( REG, "\n")

  soft  = odb[ get(mau)==reg & durometer <  68, .(Soft=.N), by=.(fishyr ) ] 
  total = odb[ get(mau)==reg & is.finite(durometer) , .(Total=.N), by=.(fishyr) ] 
  oo = soft[total, on="fishyr"]
  oo = oo[, .(Year=fishyr, Soft=round(Soft/Total*100, 1), Total=Total) ]  
  scond = shell_condition[ get(mau)==reg & shell %in% c(1,2), .(SoftSC=sum(percent), TotalSC=unique(total)[1]), by=.(Year)]
  oo = oo[scond, on="Year"]

  names(oo) = c( "Year", "Soft (D)", "Total (D)", "Soft (CC)", "Total (CC)" )
  out = gt::gt(oo) |> gt::tab_options(table.font.size = 14, data_row.padding = gt::px(1), 
      summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
      footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
      row_group.padding = gt::px(1))
  print(out)
  cat("\n\n")
} 

```
 


#### Total discards of snow crab, by weight

```{r}
#| label: tbl-fishery-discard-total 
#| eval: true
#| output: true
#| echo: false
#| results: asis
#| layout-ncol: 2
#| tbl-cap: "Average by-catch discard rate by weight observed (kg/trap haul; and standard deviation, SD)." 
#| tbl-subcap: 
#|   - "(a)"
#|   - "(b)"
#|   - "(c)"
#|   - "(d)"


for (r in 1:maus[["n"]]){
  reg = maus[["internal"]][r]
  REG = maus[["labels"]][r]
  cat( REG, "\n")
 
  o = BC[[reg]]   
  oo = o$eff_summ[ order(fishyr), ]
  names(oo) = c("Year", "Discards", "SD")
  oo$Discards = round( oo$Discards*100, 1)
  oo$SD = round( oo$SD*100, 1)
  
  names(oo) = c("Year", "Discard mean", "Discard SD")
  out = gt::gt(oo) |> gt::tab_options(table.font.size = 14, data_row.padding = gt::px(1), 
      summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
      footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
      row_group.padding = gt::px(1))
  print(out)
  cat("\n\n")
}

```

$~$


```{r}
#| label: fig-discard_maritimes
#| fig-cap: "At sea observed rate of snow crab discards relative to total catch (discard + kept), all Maritimes."
#| fig-dpi: 144
#| fig-height: 6
 
oo = NULL
for (reg in maus[["internal"]]) {
  o = BC[[reg]][["eff_summ"]]
  o$region = reg
  oo = rbind(oo, o)
}
oo$fishyr = jitter(oo$fishyr, amount=0.02)
 
pl = ggplot( oo, aes(x=fishyr, y=discard_rate, ymin=discard_rate-discard_rate_sd, ymax=discard_rate+discard_rate_sd, fill=region, colour=region) ) +
  geom_line( alpha=0.9, linewidth=1 ) +
  geom_point(aes(shape=region), size=5, alpha=0.7 )+
  geom_pointrange()  + # Vertical line with point in the middle
  geom_errorbar(width = 0.1, col="brown") + # Standard error bars
  geom_point(size = 1.5, col="darkred") +
  scale_colour_manual(values=maus[["color_map"]]) +
  scale_fill_manual(values=maus[["color_map"]]) +
  scale_shape_manual(values = maus[["shapes"]]) +
  theme_light( base_size = 18) + 
  labs(x="Year", y="Discard rate of snow crab (Observed, by weight)" ) + 
  theme( legend.position="inside", legend.position.inside=c(0.2, 0.8), legend.title=element_blank() )
  
(pl)

```

#### Bycatch of non-target species

General approach: estimate bycatch from at sea observed data and project
onto marfis data

##### Bycatch of non-target species: estimates based on fisheries **effort**

Bycatch comes from at-sea-observed effort and catch. Rescale these
naively to total snow crab fishery effort.

```{r}
#| label: tbl-fishery-discard-effort
#| eval: true
#| output: true
#| echo: false
#| results: asis
#| layout-ncol: 1
#| tbl-cap: "By-catch (kg) estimated from fisheries effort. Dots indicate low values. Where species exist in a list but there is no data, this indicates some historical bycatch. The overall average is only for the years shown." 
#| tbl-subcap: 
#|   - "(a)"
#|   - "(b)"
#|   - "(c)"
#|   - "(d)"

for (r in 1:maus[["n"]]){
  reg = maus[["internal"]][r]
  REG = maus[["labels"]][r]
  cat( REG, "\n")
  o = BC[[reg]]   
  oo = o$bycatch_table_effort
  oo[ oo==0 ] = NA
  oo[ is.na(oo) ] = "."
  out = gt::gt(oo) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
  print(out)
  cat("\n\n")
}


```
 
##### Bycatch of non-target species: estimates based on fisheries **catch**

Bycatch comes from at-sea-observed effort and catch. Rescale these
naively to total snow crab fishery catch.

```{r}
#| label: tbl-fishery-discard-catch 
#| eval: true
#| output: true
#| echo: false
#| results: asis
#| layout-ncol: 1
#| tbl-cap: "By-catch (kg) estimated from landings. Dots indicate low values. Where species exist in a list but there is no data, this indicates some historical bycatch. The overall average is only for the years shown." 
#| tbl-subcap: 
#|   - "(a)"
#|   - "(b)"
#|   - "(c)"
#|   - "(d)"

for (r in 1:maus[["n"]]){
  reg = maus[["internal"]][r]
  REG = maus[["labels"]][r]
  cat( REG, "\n")
 
  o = BC[[reg]]   
  oo = o$bycatch_table_catch 
  oo[ oo==0 ] = NA
  oo[ is.na(oo) ] = "."
  
  out = gt::gt(oo) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
  print(out)
  cat("\n\n")

}


```