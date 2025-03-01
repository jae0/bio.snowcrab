---
title: "Assessment of Scotian Shelf Snow Crab"
subtitle: "SAR"
keywords: 
  - snow crab stock status assessment 
abstract: |
  Snow crab stock status assessment.
fontsize: 12pt
metadata-files:
  - _metadata.yml
params:
  year_assessment: 2024
  year_start: 1999
  data_loc:  "~/bio.data/bio.snowcrab"
  sens: 1
  debugging: FALSE
  model_variation: logistic_discrete_historical

---


 
<!--

make quarto FN=snowcrab_sar.md DOCTYPE=html  PARAMS="-P year_assessment:2024"  --directory=~/bio/bio.snowcrab/inst/markdown

NOTE: this is not a full document but rather the basis for adding to a SAR/FSAR ... just use for values and figures and tables

-->



```{r}
#| eval: true
#| output: false
#| echo: false
#| label: setup

  require(knitr)

  knitr::opts_chunk$set(
    root.dir = data_root,
    echo = FALSE,
    out.width="6.2in",
    fig.retina = 2,
    dpi=192
  )

  # things to load into memory (in next step) via _load_results.qmd
  toget = c( "fishery_results", "fishery_model", "ecosystem" )  

```


<!-- 
# _load_results.qmd contains instructions to load data 
#  this is a shared R-script to boot strap and provide a consistent data interface
-->

{{< include _load_results.qmd >}}  

 

 
# SUMMARY


- Fishing effort in `r year_assessment` were `r e_nens` $\times 10^3$ trap hauls in N-ENS, `r e_sens` $\times 10^3$ trap hauls in S-ENS and  `r e_4x` $\times 10^3$ trap hauls in 4X. This represents a  change of `r dt_e_nens` %, `r dt_e_sens` % and  `r dt_e_4x` %, respectively, relative to the previous year. 

- Landings in `r year_assessment` were `r l_nens` t in N-ENS; `r l_sens` t in S-ENS; and `r l_4x` t in CFA 4X (season ongoing), representing a change of `r dt_l_nens`%, `r dt_l_sens`% and  `r dt_l_4x`%, respectively, relative to `r year_previous`. Total Allowable Catches (TACs) for `r year_assessment` were `r tac_nens` t, `r tac_sens` t and `r tac_4x` t in N-ENS, S-ENS and 4X, respectively. 

- Non-standardized fishery catch rates in `r year_assessment` were `r c_nens`, `r c_sens`  and `r c_4x` kg/trap haul in N-ENS, S-ENS and 4X, respectively. This represents a change of respectively, `r dt_c_nens` %, `r dt_c_sens` % and  `r dt_c_4x` % (season ongoing) relative to the previous year. Though the spatial extent of exploitation was smaller, many of the exploited area show elevated catch rates .

- Commercial catches of soft-shelled (newly moulted) Snow Crab for `r year_assessment` were `r cc_soft_nens`%, `r cc_soft_sens`% and `r cc_soft_4x`%, respectively, in N-ENS S-ENS and 4X (season ongoing). In `r year_previous`, it was `r cc_soft_nens_p`%, `r cc_soft_sens_p`% and `r cc_soft_4x_p`%, respectively. Sensecent (CC5) crab levels were higher in N- and S-ENS than in previous years, though low and inconsistent sampling effort makes this uncertain.

- Bycatch of non-target species is low (<< 1% of total catch) in all Snow Crab fishing areas; however, the last reliable estimate comes from 2017. 

- Little to no recruitment is expected for the next 1-3 years in N-ENS. In S-ENS, continued moderate levels of recruitment are expected. In 4X, low to moderate levels of recruitment are expected for 2 years. Egg and larval production is expected to be high in the next year in all areas except N-ENS. 

- In N-ENS, the modelled biomass (pre-fishery) of Snow Crab in `r year_assessment` was `r round(B_north[t0], 2)` t, relative to `r round(B_north[t1], 2)` t in `r year_previous`. In S-ENS, the modelled biomass (pre-fishery) was `r round(B_south[t0], 2)` t, relative to `r round(B_south[t1], 2)` t in `r year_previous`. In 4X, the `r year_assessment`-`r year_assessment+1` season's modelled biomass (pre-fishery) was `r round(B_4x[t0], 2)` t, relative to `r round(B_4x[t1], 2)` t in the `r year_previous`-`r year_assessment` season. 

- Average bottom temperatures observed in the 2022 Snow Crab survey were near or above historical highs in all areas. Average viable habitat surface area has declined to historical lows in 2022.

- Based on stomach sampling, Atlantic Halibut, Atlantic Wolffish, Thorny Skate, and other skate species appear to be the predominant predators of Snow Crab on the Scotian Shelf. Overall, higher predation mortality seems likely in N- and S-ENS and lower in 4X. Co-occurring species such as shrimp densities have declined, possibly due to large-scaled environmental change.

- In N-ENS, though recruitment continues at low levels, a gap in future recruitment to the fishery is expected for the next 1-3 years in N-ENS. N-ENS still exists in the healthy zone and so still has flexibility in approach. However, fishing mortality may have been higher than desirable. A more conservative harvest strategy would enable bridging this coming recruitment gap. A reduced TAC is prudent. 

- In S-ENS, recruitment to the fishery is likely to continue at a moderate rate for the upcoming season. The S-ENS stock remains in the healthy zone. Exploitation rates derived from the fishery model have been declining in recent years. S-ENS has considerable flexibility in approach. All TAC directions are viable. 

- In 4X, low to moderate levels of recruitment are expected for 2 years. 4X exists in the "cautious zone". The area is also in the southern-most extent of Snow Crab distribution in the North Atlantic and viable habitat has been depressed for many years. A reduced TAC is prudent. 





# INTRODUCTION
 
```{r crab-image, out.width='45%', echo=FALSE, fig.align = 'center', fig.cap = 'A mature male Snow Crab (*Chionoecetes opilio*).' }
include_graphics( file.path( media_loc, "snowcrab_image.png" ) )
# \@ref(fig:crab-image)
```


```{r sse-map, out.width='60%', echo=FALSE, fig.align = 'center', fig.cap = 'Map of the SSE and Crab Fishing Areas (CFAs). S-ENS is further subdivided for management into areas 23 to the northeast and 24 to the southwest.' }
include_graphics( file.path (media_loc, "area_map.png" ) )
# \@ref(fig:sse-map)
```

 
Snow Crab (*Chionoecetes opilio*; Figure \@ref(fig:crab-image)) are a circumpolar, subarctic species. In the Scotian Shelf Ecosystem (SSE; Figure \@ref(fig:sse-map)), habitat preference is generally for soft mud bottoms, at depths from 60 to 300 m and temperatures from -1 to 6${}^\circ$ C. The SSE represents the southern-most part of this distribution and so most influenced by environmental variability. More detailed information with regards to Snow Crab life history, habitat requirements and spatiotemporal distributions of different life stages can be found in Choi et al. (2022) and references therein.

Management of Snow Crab in the Scotian Shelf Snow Crab fishery are inherently precautionary:

- The spawning stock biomass is not fished and so *completely protected*.
- Conservative exploitation strategies have generally been the norm since the mid-2000s.
- Spatial refugia (Marine Protected Areas, continental slopes, western inshore CFA 24 and 4X).
- Temporal refugia (fishing seasons)
- Biological refugia (most life stages are protected, especially immature and soft-shelled, females)
- Attention to minimizing effects upon other species by implementing bycatch reduction measures (season timing, biodegradeable mesh, area closures) 
- Collaborative management which prioritizes scientific, fisheries-independent, evidence based decision making.


# Fishery performance


```{r table-fishery-nens, echo=FALSE }
ii = which(dt$Region=="cfanorth")
oo = dt[ii, c("Year", "Licenses", "TAC", "Landings", "Effort", "CPUE")] 
kable( oo, format="simple", row.names=FALSE, align="cccccc",
caption = "Fishery performance statistics in N-ENS. Units are: TACs and Landings (tons, t), Effort ($\\times 10^3$ trap hauls, th) and CPUE (kg/th).")
# \@ref(tab:table-fishery-nens)
```

 

```{r table-fishery-sens, echo=FALSE }
ii = which(dt$Region=="cfasouth")
oo = dt[ii, c("Year", "Licenses", "TAC", "Landings", "Effort", "CPUE")] 
kable( oo, format="simple", row.names=FALSE, align="cccccc",
caption = "Fishery performance statistics in S-ENS. Units are: TACs and Landings (tons, t), Effort ($\\times 10^3$ trap hauls, th) and CPUE (kg/th).")
# \@ref(tab:table-fishery-nens)
```
 


```{r table-fishery-4x, echo=FALSE }
ii = which(dt$Region=="cfa4x")
oo = dt[ii,c("Year", "Licenses", "TAC", "Landings", "Effort", "CPUE")]
kable(oo, format="simple", row.names=FALSE, align="cccccc",
caption = "Fishery performance statistics in 4X. Units are: TACs and Landings (tons, t), Effort ($\\times 10^3$ trap hauls, th) and CPUE (kg/th). There were no landings or TACs in 2018/2019 due to indications of low abundance. The 2022 season is ongoing.")
# \@ref(tab:table-fishery-4x)
```


## Fishing effort

Fishing effort in `r year_assessment` was, respectively, `r e_nens` $\times 10^3$, `r e_sens` $\times 10^3$ and `r e_4x` $\times 10^3$ trap hauls in N-ENS, S-ENS and 4X. Relative to the previous year, these represent changes of `r dt_e_nens` %, `r dt_e_sens` % and  `r dt_e_4x` %, respectively (Tables \@ref(tab:table-fishery-nens), \@ref(tab:table-fishery-sens), \@ref(tab:table-fishery-4x), Figure \@ref(fig:effort-timeseries)). Fishing effort was consistent between `r year_assessment` and `r year_previous` in terms of  spatial distribution. In S-ENS, there was, however, a spatial contraction to inshore areas and away from the area 23-24 boundary (Figure \@ref(fig:effort-map)). This is presumed to be related to, in part, warm water incursions into the area (see Ecosystem considerations, below).


```{r effort-timeseries, echo=FALSE, out.width='60%', fig.align='center', fig.cap = 'Temporal variations in fishing effort $\\times 10^3$ trap hauls.' }
fn1=file.path( data_loc, "assessments", year_assessment, "timeseries", "fishery",   "effort.ts.png" )
knitr::include_graphics( fn1 ) 
# \@ref(fig:effort-timeseries)  
```

```{r effort-map, echo=FALSE, out.width='45%', fig.show='hold',  fig.align='center', fig.cap = 'Snow Crab fishing effort from fisheries logbook data for previous and current years (no $\\times 10^3$ per 10 km X 10 km grid).' }
loc0= file.path( data_loc, "output", "maps", "logbook", "snowcrab", "annual", "effort" )
fn1 = file.path( loc0, paste( "effort", year_previous,   "png", sep=".") ) 
fn2 = file.path( loc0, paste( "effort", year_assessment, "png", sep=".") ) 
include_graphics(  c(fn1, fn2) )
#  \@ref(fig:landings-map) 
```
 


## Fishery landings and TACs

Landings across time are shown in Figure \@ref(fig:landings-timeseries). In `r year_assessment`, they were `r l_nens`, `r l_sens` and `r l_4x` t, in N-ENS, S-ENS and 4X (season ongoing), respectively. Relative to `r year_previous`, they represent changes of `r dt_l_nens`%, `r dt_l_sens`% and  `r dt_l_4x`%, respectively (Tables \@ref(tab:table-fishery-nens), \@ref(tab:table-fishery-sens), \@ref(tab:table-fishery-4x)). Total Allowable Catches (TACs) for `r year_assessment` were `r tac_nens` t, `r tac_sens` t and `r tac_4x` t in N-ENS, S-ENS and 4X, respectively. 


Additional carry forward allowance was implemented by DFO Fisheries Management of up to 25% of the 2020 quota to the 2021 season for N-ENS and S-ENS Snow Crab fishery. These were a response to COVID-19 related uncertainties with safe fishing activity. The carry-forward amounts were: 11.2 t in N-ENS, and 217.4 t in S-ENS (Tables \@ref(tab:table-fishery-nens), \@ref(tab:table-fishery-sens)). In 2022, landings in all areas were below respective TACs.


The landings in N-ENS for 2022 and 2021 were similar in their spatial patterns (Figure \@ref(fig:landings-map)). In S-ENS, landings, as with fishing effort, shifted slightly inshore and away from the area 23-24 boundary (Figure \@ref(fig:landings-map)). There were no landings on the continental slope areas of S-ENS in 2022; it continues to serve as a "reserve" for Snow Crab from fishing. The landings in 4X for 2022 as with 2021, were primarily in the area just south of Sambro, bordering onto area 24. In N-ENS, most landings occured in the spring. 
  

```{r landings-timeseries, echo=FALSE, out.width='60%',  fig.align='center', fig.cap = 'Landings (t) of Snow Crab on the SSE. For 4X, the year refers to the starting year of the season. Inset is a closeup view of the timeseries for N-ENS and 4X.'}
include_graphics( file.path( data_loc, "assessments", year_assessment, "timeseries", "fishery",   "landings.ts.png" ) )
# \@ref(fig:landings-timeseries)
``` 


```{r landings-map, echo=FALSE, out.width='45%', fig.show='hold', fig.align='center', fig.cap = 'Snow Crab landings from fisheries logbook data for previous and current years (tons per 10 km x 10 km grid).' }
loc0= file.path( data_loc, "output", "maps", "logbook", "snowcrab", "annual", "landings" )
fn1 = file.path( loc0, paste( "landings", year_previous,   "png", sep=".") ) 
fn2 = file.path( loc0, paste( "landings", year_assessment, "png", sep=".") ) 
knitr::include_graphics( c(fn1, fn2 ) )
#  \@ref(fig:landings-map)  
```

\clearpage


## Fishery catch rates

Non-standardized fishery catch rates in `r year_assessment` were `r c_nens`, `r c_sens`  and `r c_4x` kg/trap haul in N-ENS, S-ENS and 4X, respectively. This represents a change of respectively, `r dt_c_nens` %, `r dt_c_sens` % and  `r dt_c_4x` % (season ongoing) relative to the previous year (Tables \@ref(tab:table-fishery-nens), \@ref(tab:table-fishery-sens), \@ref(tab:table-fishery-4x), Figures \@ref(fig:cpue-timeseries)). Though the spatial extent of exploitation was smaller, many of the exploited areas show elevated catch rates (Figure \@ref(fig:cpue-map)).


```{r cpue-timeseries, echo=FALSE, out.width='60%', fig.align='center', fig.cap = 'Temporal variations in crude catch rates of Snow Crab (kg/trap haul).'}
include_graphics( file.path( data_loc, "assessments", year_assessment, "timeseries", "fishery",   "cpue.ts.png" ) ) 
# \@ref(fig:cpue-timeseries)  
```

 

```{r cpue-map, echo=FALSE, out.width='45%', fig.show='hold', fig.align='center', fig.cap = 'Snow Crab crude catch rates on the Scotian Shelf for previous and current years. Units are kg/trap haul per 10 km x 10 km grid.' }
loc0= file.path( data_loc, "output", "maps", "logbook", "snowcrab", "annual", "cpue" )
fn1 = file.path( loc0, paste( "cpue", year_previous,   "png", sep=".") ) 
fn2 = file.path( loc0, paste( "cpue", year_assessment, "png", sep=".") ) 
knitr::include_graphics( c(fn1, fn2 ) )
# \@ref(fig:cpue-map)  
```


## At-Sea-Observed information

Carapace condition of the fished component is determined from At-Sea-Observed catches. In 2021, both N-ENS and 4X were not sampled by At-Sea-Observers. In 2022, 4X was not sampled by At-Sea-Observers. Estimates of carapace condition since 2020 are unreliable as they represent only small areas of the fishing grounds and short time periods relative to the whole fishing season (Figure \@ref(fig:observer-locations-map)). 

```{r observer-locations-map, out.width='45%', fig.show='hold', fig.align='center', fig.cap= 'Snow Crab At-sea-observer locations.' }
loc = file.path( data_loc, "output", "maps", "observer.locations" )
yrsplot = year_assessment + c(0:-4)
fn4 = file.path( loc, paste( "observer.locations", yrsplot[4], "png", sep=".") )
fn3 = file.path( loc, paste( "observer.locations", yrsplot[3], "png", sep=".") )
fn2 = file.path( loc, paste( "observer.locations", yrsplot[2], "png", sep=".") )
fn1 = file.path( loc, paste( "observer.locations", yrsplot[1], "png", sep=".") )
include_graphics( c( fn4, fn3, fn2, fn1) )
# \@ref(fig:observer-locations-map)  
#*Figure XXX. Snow Crab observer locations.*
```

In the exploited fraction of Snow Crab, Carapace Condition (CC) is an index of the approximate time since the last molt, and so describes the relative development and subsequent decay of the carapace. CC1 signifies a newly molted crab, soft-shelled, with no epibiont (e.g., barnacles) growth. CC2 crab have begun to harden, but is still considered to be soft and of no commercial value. CC3 and CC4 represent crab preferred by industry. The oldest carapace condition (CC5) signifies extensive shell decay with no expectation of survival into the next year.

Commercial catches of soft-shelled crab were `r cc_soft_nens`% (low sampling), `r cc_soft_sens`% (low sampling) and `r cc_soft_4x`% (no sampling; season ongoing) in N-ENS, S-ENS and 4X, respectively for `r year_assessment`. In `r year_previous`, it was `r cc_soft_nens_p`% (no sampling), `r cc_soft_sens_p`% (low sampling) and `r cc_soft_4x_p`% (no sampling), respectively. Generally, higher soft-shell indicates incoming recruitment to the fishery and their  handling potential and unnecessary handling/discard mortality.

In 2022, CC5 crab levels were higher in N- and S-ENS than in previous years, though again low and inconsistent sampling effort makes this uncertain.

Bycatch in the Snow Crab fishery is also monitored from the At-Sea-Observed catches. However, due to a paupacy of consistent data, very little can be said of bycatch after 2017. Historically, bycatch in the Snow Crab fishery has been minimal (Tables 4, 5) with increasing levels as a function of increasing water temperature: bycatch is higher in warmer conditions, primarily other Crustacea (crab and lobster). Low bycatch has been attributed to trap design (top entry conical traps), the large mesh size (5.25 inches, knot to knot) and the passive nature of the gear (Hebert et al., 2001).  


*Table 4. Bycatch (kg) estimates from the N-ENS and S-ENS Snow Crab fishery. Estimates are extrapolated from At-sea-observed bycatch and biomass of catch (bycatch = [observed biomass of bycatch species / observed landings of Snow Crab] X total landings of Snow Crab). Reliable species specific data beyond 2017 are currently unavailable.*

| Species  |  2015   | 2016  | 2017  |  
| :---------------------:    | :-------:  | :-------: |:-------: |  
|Rock Crab | 19  | 0  | 0 | 
|Cod | 187 | 84 | 353 |
|Jonah Crab | 19 | 854 | 0 |
|Northern Stone Crab | 0 | 670 | 18 |
|Toad Crab | 0 | 84 | 35 |
|Soft Coral | 0 | 0 | 18 |
|Basket Star | 0 | 0 | 18 |
|Sea Urchin | 0 | 33 | 18 |
|Sand Dollars | 0 | 17 | 0 |
|Purple Starfish | 0 | 0 | 35 |
|Sea Cucumbers | 19 | 50 | 495 |
|Whelk | 0 | 17 | 0 |
|Winter Flounder | 0 | 0 | 35 |
|Eelpout | 0 | 0 | 35 |
|Redfish | 75 | 50 | 247 |
|Sea Raven | 37 | 33 | 0 |
|Skate | 0 | 67 | 18 |
|Northern Wolffish | 112 | 17 | 0 |
|Spotted Wolffish | 0 | 0 | 194 |
|Striped Wolffish | 149 | 100 | 371 |
|Total Bycatch | 617  | 2076 | 1890 |




*Table 5. Bycatch (kg) estimates from the 4X Snow Crab fishery. Estimates are extrapolated from At-sea-observed bycatch and biomass of catch (bycatch = [observed biomass of bycatch species / observed landings of Snow Crab] X total landings of Snow Crab). Reliable species specific data beyond 2017 are currently unavailable.*

| Species  |  2015   | 2016  | 2017  |  
| :---------------------:    | :-------:  | :-------: |:-------: |  
|American Lobster | 98 |  48  | 55 |  
|Cod  | 0  | 16  | 0 |  
|Jonah Crab  | 0  | 16  | 14  | 
|Rock Crab  | 0  | 0  | 14 |  
|Lumpfish  | 11  | 0  | 0  |  
|Northern Stone Crab  | 130  | 81 |  82 |   
|Redfish  | 0  | 0  | 14 | 
|Sea Raven  | 239  | 0  | 41  |  
|Total Bycatch  | 478  | 161  | 219  |  



## Stock status

Survey catch rates are are confounded by numerous factors that vary across space and time. Adjustment for these confounding influences were completed with Conditional AutoRegressive SpatioTemporal Models (*CARSTM*; Choi 2020, Choi et al. 2022 and references therein).  The survey was not conducted in 2020 due to Covid-19 related health and safety uncertainties. In 2022, the survey did not complete due to mechanical issues with the survey vessel; inshore areas of S-ENS were most affected (Figure \@ref(fig:survey-locations-map)). 


```{r survey-locations-map, out.width='45%', fig.show='hold', fig.align='center', fig.cap= 'Snow Crab survey locations.' }
loc = file.path( data_loc, "output", "maps", "survey.locations" )
yrsplot = setdiff( year_assessment + c(0:-9), 2020)
fn6 = file.path( loc, paste( "survey.locations", yrsplot[6], "png", sep=".") )
fn5 = file.path( loc, paste( "survey.locations", yrsplot[5], "png", sep=".") )
fn4 = file.path( loc, paste( "survey.locations", yrsplot[4], "png", sep=".") )
fn3 = file.path( loc, paste( "survey.locations", yrsplot[3], "png", sep=".") )
fn2 = file.path( loc, paste( "survey.locations", yrsplot[2], "png", sep=".") )
fn1 = file.path( loc, paste( "survey.locations", yrsplot[1], "png", sep=".") )
include_graphics( c(fn2, fn1) )
# \@ref(fig:survey-locations-map)  
```



## Recruitment

Based on size-frequency histograms of the male Snow Crab population, little to no recruitment is expected for the next 1-3 years in N-ENS (Figure \@ref(fig:sizefeq-male)). In S-ENS, continued moderate levels of recruitment are expected. In 4X, low to moderate levels of recruitment are expected for 2 years. 


```{r sizefeq-male, out.width='90%', echo=FALSE, fig.align='center', fig.cap = 'Size-frequency (areal density; no/km$^2$) histograms by carapace width of male Snow Crab. The vertical line represents the legal size (95 mm). Immature animals are shown with light coloured bars, mature with dark.'}
include_graphics(  file.path( data_loc, "assessments", year_assessment, "figures", "size.freq", "survey",  "male.denl.png" )  )
# \@ref(fig:sizefeq-male)
```
  

## Reproduction

In all areas, there was substantial and continued recruitment of female crab into the mature (egg-bearing) stage of the population from 2016-2022 (Figure \@ref(fig:sizefeq-female)). However, in N-ENS, a decline in numerical densities of both the mature and adolescent components was observed in 2022. Egg and larval production is expected to be moderate to high in the next year in all areas except N-ENS. 

 

```{r sizefeq-female, out.width='90%', echo=FALSE, fig.align='center', fig.cap = 'Size-frequency (areal density; no/km$^2$) histograms by carapace width of female Snow Crab. Immature animals are shown with light coloured bars, mature with dark.'}
include_graphics(  file.path( data_loc, "assessments", year_assessment, "figures", "size.freq", "survey",  "female.denl.png" )  )
# \@ref(fig:sizefeq-female)
```
 

```{r fmat-map, echo=FALSE, out.width='45%', fig.show='hold', fig.align='center', fig.cap= 'Mature female density log$_{10}$(no/km$^2$)  from the Snow Crab survey.' }
loc = file.path( data_loc, "output", "maps", "survey", "snowcrab", "annual", "totno.female.mat" )
yrsplot = setdiff( year_assessment + c(0:-9), 2020)
fn4 = file.path( loc, paste( "totno.female.mat", yrsplot[4], "png", sep=".") )
fn3 = file.path( loc, paste( "totno.female.mat", yrsplot[3], "png", sep=".") )
fn2 = file.path( loc, paste( "totno.female.mat", yrsplot[2], "png", sep=".") )
fn1 = file.path( loc, paste( "totno.female.mat", yrsplot[1], "png", sep=".") )
include_graphics( c(   fn2, fn1) )
# \@ref(fig:fmat-map)  
```


## Sex ratios

The sex ratios (proportion female) of the mature component is particularly important as unbalanced ratios can impact reproductive success (Figure \@ref(fig:sexratio-mature)). In the ESS, there is generally a lack of females. The exception being in inshore areas and areas with high bottom slopes (Figure \@ref(fig:sexratio-map)). A decline in sex ratios has been observed since 2017 in N-ENS (Figure \@ref(fig:sexratio-mature)). In S-ENS the sex ratio increased from 20% in 2021 to just under 35% in 2022. In 4X, mature sex ratios are more stable and balanced and currently near the 50% level.  




```{r sexratio-mature, out.width='60%', echo=FALSE, fig.align='center', fig.cap = 'Timeseries of sex ratios (proportion female) of mature Snow Crab.'  }
include_graphics( file.path( data_loc, "assessments", year_assessment, "timeseries", "survey", "sexratio.mat.png") )
# \@ref(fig:sexratio-mature)
```
 
  


```{r sexratio-map, echo=FALSE, out.width='45%', fig.show='hold', fig.align='center', fig.cap= 'Map of sex ratios (proportion female) of mature Snow Crab.'}
yrsplot = setdiff( year_assessment + c(0:-4), 2020)
loc = file.path( data_loc, "output", "maps", "survey", "snowcrab", "annual", "sexratio.mat" )
fn4 = file.path( loc, paste( "sexratio.mat", yrsplot[4], "png", sep=".") )
fn3 = file.path( loc, paste( "sexratio.mat", yrsplot[3], "png", sep=".") )
fn2 = file.path( loc, paste( "sexratio.mat", yrsplot[2], "png", sep=".") )
fn1 = file.path( loc, paste( "sexratio.mat", yrsplot[1], "png", sep=".") )
include_graphics( c( fn2, fn1  ) )
# \@ref(fig:sexratio-map)  
```


## Biomass Density

The fishable component is defined as Snow Crab that are male, mature, and larger than 95 mm CW. The crude, unadjusted, geometric mean fishable **biomass density** (per unit swept area by the trawl) are shown in Figures \@ref(fig:fbGMTS), and \@ref(fig:fbgeomean-map). A peak in crude biomass densities was observed in 2009 to 2014 and has since been declining in all areas. Note that, high and low biomass density areas fluctuate with time (Figure \@ref(fig:fbgeomean-map)). Biomass density, however, does not equate to total biomass as the areas occupied by crab can contract, expand and shift with environmental conditions and ecosystem change. 


```{r fbGMTS, out.width='60%', echo=FALSE, fig.align='center', fig.cap = 'The crude, unadjusted geometric mean fishable biomass density (t/km$^2$) from the Snow Crab survey. Error bars represent 95\\% Confidence Intervals. Note the absence of data in 2020. Prior to 2004, surveys were conducted in the Spring.'}
fn = file.path(data_loc,'assessments',year_assessment,'timeseries','survey','R0.mass.png')
include_graphics( c(fn) )
#\@ref(fig:fbGMTS)
```


```{r fbgeomean-map, echo=FALSE, out.width='45%', fig.show='hold', fig.align='center', fig.cap= 'Snow Crab survey fishable component biomass density (t/km$^2$). Note, there is no data in 2020.' }
loc = file.path( data_loc, 'output', 'maps', 'survey', 'snowcrab', 'annual', 'R0.mass')
yrsplot =  setdiff(year_assessment + c(0:-9), 2020 ) 
fn6 = file.path( loc, paste( 'R0.mass', yrsplot[6], 'png', sep='.') )
fn5 = file.path( loc, paste( 'R0.mass', yrsplot[5], 'png', sep='.') )
fn4 = file.path( loc, paste( 'R0.mass', yrsplot[4], 'png', sep='.') )
fn3 = file.path( loc, paste( 'R0.mass', yrsplot[3], 'png', sep='.') )
fn2 = file.path( loc, paste( 'R0.mass', yrsplot[2], 'png', sep='.') )
fn1 = file.path( loc, paste( 'R0.mass', yrsplot[1], 'png', sep='.') )
include_graphics( c(fn2, fn1) )
# \@ref(fig:fbgeomean-map)  
```


## Biomass Index

The fishable **biomass index** (statistically adjusted for covariates and autocorrelation; Figure \@ref(fig:fbindex-map)) was computed using conditional auto-regressive spatio-temporal models (Choi 2020). This approach models Snow Crab numerical abundance and mean size with environmental (depth, substrate, temperature) and biological factors (species composition) as covariates. Upon aggregation we see that the overall biomass has had several cycles (Figure \@ref(fig:fbindex-timeseries)). Further, the biomass index model infers the spatiotemporal distribution of the biomass density of the fishable component from the covariates measured in that year (Figure \@ref(fig:fbindex-map)); note also the aggregate timeseries with elevated uncertainty for the 2020 estimate (Figure \@ref(fig:fbindex-timeseries)). 



```{r fbindex-map, echo=FALSE, out.width='45%', fig.show='hold', fig.align='center', fig.cap= 'Biomass index log10(t/km$^2$)  predicted from the Snow Crab survey.' }
loc = file.path( data_loc, 'modelled', 'default_fb', 'predicted_biomass_densities' )
yrsplot =  year_assessment + c(0:-10)
fn10 = file.path( loc, paste( 'biomass', yrsplot[10], 'png', sep='.') )
fn9 = file.path( loc, paste( 'biomass', yrsplot[9], 'png', sep='.') )
fn8 = file.path( loc, paste( 'biomass', yrsplot[8], 'png', sep='.') )
fn7 = file.path( loc, paste( 'biomass', yrsplot[7], 'png', sep='.') )
fn6 = file.path( loc, paste( 'biomass', yrsplot[6], 'png', sep='.') )
fn5 = file.path( loc, paste( 'biomass', yrsplot[5], 'png', sep='.') )
fn4 = file.path( loc, paste( 'biomass', yrsplot[4], 'png', sep='.') )
fn3 = file.path( loc, paste( 'biomass', yrsplot[3], 'png', sep='.') )
fn2 = file.path( loc, paste( 'biomass', yrsplot[2], 'png', sep='.') )
fn1 = file.path( loc, paste( 'biomass', yrsplot[1], 'png', sep='.') )
include_graphics( c(  fn2, fn1) )
# \@ref(fig:fbindex-map)  
#*Figure XXX. Biomass index log${}_{10}$(t/km${}^2$) predicted from the Snow Crab survey. Note there was no survey in 2020.*
```





```{r fbindex-timeseries, out.width='60%', echo=FALSE, fig.align='center', fig.cap = 'The fishable biomass index (t) predicted from the Snow Crab survey. Error bars represent Bayesian 95\\% Credible Intervals. Note large errors in 2020 when there was no survey.' }
include_graphics( file.path( data_loc, 'modelled', 'default_fb', 'aggregated_biomass_timeseries', 'biomass_M0.png') )
# \@ref(fig:fbindex-timeseries)
```



The magnitudes of the biomass index are optimistically high as the spatial expansion uses areal units with large surface areas, on average, much larger than the patchiness of Snow Crab distributions (Choi 2020).  As such, it should only be seen as a spatially and temporally consistent relative index of abundance. The spatial distribution of the biomass index has been consistent over the past six years, with a peak in overall biomass index in 2019 and 2020 (Figure \@ref(fig:fbindex-map)). Since then, a reduction in the biomass index was observed throughout the region. Contraction of spatial range in 4X and the western parts of S-ENS were evident in 2021 to 2022. Upon aggregation, the predicted biomass index declined marginally in all areas (Figure \@ref(fig:fbindex-timeseries)).


## Modelled Biomass

The biomass index along with fishery removals are used to fit a logistic biomass dynamics model ("Model 1") to determine fishable **modelled biomass** (Figure \@ref(fig:logisticPredictions)) and relevant biological reference points (i.e., carrying capacity and fishing mortality at maximum sustainable yield, or F~MSY~). In N-ENS, the modelled biomass (pre-fishery) of Snow Crab in `r year_assessment` was `r round(B_north[t0], 2)` t, relative to `r round(B_north[t1], 2)` t in `r year_previous`. In S-ENS, the `r year_assessment` modelled biomass (pre-fishery) was `r round(B_south[t0], 2)` t, relative to `r round(B_south[t1], 2)` t in `r year_previous`. In 4X, the modelled biomass (pre-fishery) for the `r year_assessment`-`r year_assessment+1` season was `r round(B_4x[t0], 2)` t, relative to `r round(B_4x[t1], 2)` t for the `r year_previous`-`r year_assessment` season. 
 

```{r logisticPredictions, out.width='32%', echo=FALSE, fig.show='hold', fig.align='center', fig.cap = 'Fishable posterior mean modelled biomass (pre-fishery; kt) are shown in dark orange for N-ENS, S-ENS and 4X (left, middle and right). Light orange are posterior samples of modelled biomass (pre-fishery; kt) to illustrate the variability of the predictions. The biomass index (post-fishery, except prior to 2004) after model adjustment by the model catchability coefficient is in gray.' } 
loc = file.path( data_loc, 'fishery_model', year_assessment, 'logistic_discrete_historical' )
fn1 = file.path( loc, 'plot_predictions_cfanorth.png' ) 
fn2 = file.path( loc, 'plot_predictions_cfasouth.png' ) 
fn3 = file.path( loc, 'plot_predictions_cfa4x.png' ) 
knitr::include_graphics(c(fn1, fn2, fn3) )
# \@ref(fig:logisticPredictions)
```



 

## Fishing Mortality

In N-ENS, the `r year_assessment` fishing mortality is estimated to have been `r round(FM_north[t0],3)` (annual exploitation rate of `r round(100*(exp(FM_north[t0])-1),2)`%), up from the  `r year_previous` rate of `r round(FM_north[t1],3)` (annual exploitation rate of `r round(100*(exp(FM_north[t1])-1),1)`%; Figure \@ref(fig:logisticFishingMortality)).

In S-ENS, the `r year_assessment` fishing mortality is estimated to have been `r round(FM_south[t0],3)` (annual exploitation rate of `r round(100*(exp(FM_south[t0])-1),1)`%), decreasing marginally from the `r year_previous` rate of `r round(FM_south[t1],3)` (annual exploitation rate of `r round(100*(exp(FM_south[t1])-1),1)`%; Figure \@ref(fig:logisticFishingMortality)). Localized exploitation rates are likely higher, as not all areas for which biomass is estimated are fished (e.g., continental slope areas and western, inshore areas of CFA 24).

In 4X, the `r year_assessment`-`r year_assessment+1` season (ongoing), fishing mortality is estimated to have been `r round(FM_4x[t0],3)` (annual exploitation rate of `r round(100*(exp(FM_4x[t0])-1),1)`%), decreasing from the `r year_assessment-1`-`r year_assessment` season rate of `r round(FM_4x[t1],3)` (annual exploitation rate of `r round(100*(exp(FM_4x[t1])-1),1)`%; Figure \@ref(fig:logisticFishingMortality)). Localized exploitation rates are likely higher, as not all areas for which biomass is estimated are fished. 



```{r logisticFishingMortality, out.width='32%', echo=FALSE,  fig.show='hold', fig.align='center', fig.cap = 'Time-series of modelled instantaneous fishing mortality for N-ENS (left), S-ENS (middle), and 4X (right). Samples of the posterior densities are presented, with the darkest line being the mean.' }
  odir = file.path( data_loc, "fishery_model", year_assessment, "logistic_discrete_historical" )
  fn1 = file.path( odir, "plot_fishing_mortality_cfanorth.png" ) 
  fn2 = file.path( odir, "plot_fishing_mortality_cfasouth.png" ) 
  fn3 = file.path( odir, "plot_fishing_mortality_cfa4x.png" ) 
include_graphics(c(fn1, fn2, fn3) )
# \@ref(fig:logisticFishingMortality)
```


## Reference Points
 
Reference points are used to guide harvest strategies (Canada Gazette 2022; DFO 2013; Figures \@ref(fig:logistic-hcr), \@ref(fig:logistic-hcr)). Lower and Upper Stock Reference points are 25% and 50% of carrying capacity which delineate "critical", "cautious" and "healthy" zones. The Upper Removal Reference point is the exploitation rate that the fishery tries to stay below; it is defined in terms of the fishing mortality associated with Maximum Sustainable Yield (FMSY). In Model 1, FMSY = $r /2$. As $r \approx 1$ for snow crab, FMSY $\approx 0.5$ is expected. 

The operational target exploitation changes depending upon the "zone" in which a population lands. When in the "healthy" zone, the rule of thumb has been to keep annual exploitation rates between 10% to 30% of the available biomass ($F = 0.11, 0.36$, respectively). In the "cautious" zone, the rule of thumb has been to keep annual exploitation rates between 0% to 20% ($F = 0, 0.22$, respectively). In the "critical" zone, fishery closure is considered until recovery is observed, where recovery indicates at a minimum, modelled biomass > LSR. Other biological and ecosystem considerations such as recruitment, spawning stock (female) biomass, size structure, sex ratios and environmental and ecosystem conditions, provide additional guidance within each range.


*Table 6. Reference points from the logistic biomass dynamics fishery model (Model 1): K is Carrying capacity (kt); and r is Intrinsic rate of increase (non-dimensional). Note that FMSY (fishing mortality associated with 'Maximum Sustainable Yield') is r/2. Similarly, BMSY (biomass associated with 'Maximum Sustainable Yield') is K/2. SD is posterior Standard deviation.*


|       |  $K$ [SD] | $r$ [SD] |
| :---: |    :----:   | :---: |
| N-ENS | `r K_north` [`r K_north_sd`] | `r r_north` [`r r_north_sd`] |
| S-ENS | `r K_south` [`r K_south_sd`] | `r r_south` [`r r_south_sd`] |
| 4X    | `r K_4x`   [`r K_4x_sd`] | `r r_4x` [`r r_4x_sd`] |


The current estimates of key Reference Points from Model 1 are shown in Table 6 and Figure \@ref(fig:ReferencePoints). The related PA thresholds can be computed as:

  - Lower Stock Reference (LSR): $K/4$
  - Upper Stock Reference (USR): $K/2$
  - Upper Removal Reference (URR): keep fishing mortality below *FMSY* = $r/2$

The current state of the fishable components and the above landmarks from Model 1 (Figure \@ref(fig:logistic-hcr)) suggests that:  

  - N-ENS is in the "healthy" zone
  - S-ENS is in the "healthy" zone
  - 4X is in the "cautious" zone



```{r ReferencePoints, out.width='60%', echo=FALSE,   fig.align='center', fig.cap = 'Harvest control rules for the Scotian Shelf Snow Crab fisheries.' }
include_graphics( file.path( media_loc, 'harvest_control_rules.png') ) 
# \@ref(fig:ReferencePoints)
```
 


```{r logistic-hcr, out.width='32%', echo=FALSE, fig.show='hold', fig.align='center', fig.cap = 'Reference Points (fishing mortality and modelled biomass) for N-ENS (left), S-ENS (middle), and 4X (right), derived from the biomass dynmaics model. The large yellow dot indicates most recent year.' }
  odir = file.path( data_loc, 'fishery_model', year_assessment, 'logistic_discrete_historical' )
  fn1 = file.path( odir, 'plot_hcr_cfanorth.png' ) 
  fn2 = file.path( odir, 'plot_hcr_cfasouth.png' ) 
  fn3 = file.path( odir, 'plot_hcr_cfa4x.png' ) 
  include_graphics(c(fn1, fn2, fn3) )
#  \@ref(fig:logistic-hcr)
```


It should be noted that using these parameters assume that the population dynamics are well described by the fishery model. This is, of course, **not true**. For example, the observation of fisheries landings is assumed to be known without error. This is not true as illegal and unreported exploitation occurs. These and other unaccounted factors can easily bias parameter estimates. As such, **caution is required in using these reference points.** Other contextual reference points must be used in conjunction:

  - Strength of recruitment (short-term, long-term)
  - Strength of spawning stock (females)
  - Ecosystem variability (predator and prey trends and distributions) within norms)
  - Habitat viability within norms
  - Availability of spatial and temporal refugia within norms 
 

# Ecosystem Considerations

## Bottom Temperature

A general warming trend has been observed in the Snow Crab survey since the early 1990s on the Scotian Shelf (Choi et al. 2022). Temperatures are more stable in N-ENS than S-ENS; 4X exhibits the most erratic and highest annual mean bottom temperatures.  The average temperature is found to have increased well beyond the $7^\circ$C threshold in 4X. N-ENS and S-ENS also continued to experience historical highs in bottom temperature and elevated spatial variability of bottom temperatures (Figures \@ref(fig:bottom-temperatures), \@ref(fig:bottom-temperatures-map)). 
 


```{r bottom-temperatures, out.width='60%', echo=FALSE, fig.align='center', fig.cap = 'Temporal variations in bottom temperature estimated from a historical analysis of temperature data. Red horizontal line is at $7^\\circ$C. Presented are 95\\% Credible Intervals of spatial variability in temperature at each time slice, after adjustment for spatiotemporal autocorrelation' }
include_graphics( file.path( data_loc, 'assessments', year_assessment, 'timeseries', 'temperature_bottom.png') )
# \@ref(fig:bottom-temperatures)
```


```{r bottom-temperatures-map, out.width='45%', echo=FALSE, fig.show='hold', fig.align='center', fig.cap = 'Spatial variations in bottom temperature estimated from a historical analysis of temperature data for 1 September.' }

loc = file.path( data_root, 'aegis', 'temperature', 'modelled', 'default', 'maps' )
yrsplot =  year_assessment + c(0:-10)
fn10 = file.path( loc, paste( 'predictions.',  yrsplot[10], '.0.75',  '.png', sep='') )
fn9  = file.path( loc, paste( 'predictions.',  yrsplot[9],  '.0.75',  '.png', sep='') )
fn8  = file.path( loc, paste( 'predictions.',  yrsplot[8],  '.0.75',  '.png', sep='') )
fn7  = file.path( loc, paste( 'predictions.',  yrsplot[7],  '.0.75',  '.png', sep='') )
fn6  = file.path( loc, paste( 'predictions.',  yrsplot[6],  '.0.75',  '.png', sep='') )
fn5  = file.path( loc, paste( 'predictions.',  yrsplot[5],  '.0.75',  '.png', sep='') )
fn4  = file.path( loc, paste( 'predictions.',  yrsplot[4],  '.0.75',  '.png', sep='') )
fn3  = file.path( loc, paste( 'predictions.',  yrsplot[3],  '.0.75',  '.png', sep='') )
fn2  = file.path( loc, paste( 'predictions.',  yrsplot[2],  '.0.75',  '.png', sep='') )
fn1  = file.path( loc, paste( 'predictions.',  yrsplot[1],  '.0.75',  '.png', sep='') )

include_graphics( c( fn2, fn1) )
# \@ref(fig:bottom-temperatures-map)
```


 
## Viable habitat

Snow Crab being cold water stenotherms, stability of environmental conditions is critical for their survival. The Maritimes Region being at the confluence of many oceanic currents renders the area highly variable. Rapid climate change and uncertainty exacerbates this situation. The viable habitat estimated for each area across time has shown some variations (Figures \@ref(fig:fb-habitat-timeseries), \@ref(fig:fb-habitat-map)) in the historical record. As can be seen, 4X showed a significantly lower average viable habitat levels relative to the N-ENS and S-ENS.  A peak in average probability of observing fishable snow crab ("viable habitat") was observed in 2010 for 4X, 2011 for N-ENS and 2012 for S-ENS. Since 2015, the average viable habitat has declined to historical lows and remained so in 2022. 
 

```{r fb-habitat-timeseries, out.width='60%', echo=FALSE,   fig.align='center', fig.cap = 'Habitat viability (probability; fishable Snow Crab). Means and 95\\% Credible Intervals are presented.' }
loc = file.path( data_loc, 'modelled', 'default_fb', 'aggregated_habitat_timeseries' )
include_graphics( file.path( loc, 'habitat_M0.png') )
# \@ref(fig:fb-habitat-timeseries)
```
  
 

```{r fb-habitat-map, out.width='45%', fig.show='hold', fig.align='center', fig.cap= 'Habitat viability (probability; fishable Snow Crab).' }
loc = file.path( data_loc, 'modelled', 'default_fb', 'predicted_habitat' )
yrsplot =  year_assessment + c(0:-10)
fn10 = file.path( loc, paste( 'habitat.', yrsplot[10], '.png', sep='') )
fn9 = file.path( loc, paste( 'habitat.', yrsplot[9], '.png', sep='') )
fn8 = file.path( loc, paste( 'habitat.', yrsplot[8], '.png', sep='') )
fn7 = file.path( loc, paste( 'habitat.', yrsplot[7], '.png', sep='') )
fn6 = file.path( loc, paste( 'habitat.', yrsplot[6], '.png', sep='') )
fn5 = file.path( loc, paste( 'habitat.', yrsplot[5], '.png', sep='') )
fn4 = file.path( loc, paste( 'habitat.', yrsplot[4], '.png', sep='') )
fn3 = file.path( loc, paste( 'habitat.', yrsplot[3], '.png', sep='') )
fn2 = file.path( loc, paste( 'habitat.', yrsplot[2], '.png', sep='') )
fn1 = file.path( loc, paste( 'habitat.', yrsplot[1], '.png', sep='') )
include_graphics( c(  fn2, fn1) )
# \@ref(fig:fb-habitat-map)  
```



## Predators, preys, competitors

Being long-lived, the influence of predators can be significant. Especially important are predators of the smaller immature and female Snow Crab. Increasing predation not only lowers the abundance and recruitment, it can also reduce the reproductive potential of Snow Crab and therefore long-term population cycles. N-ENS and S-ENS are well known to have skewed sex ratios with few mature females for extended periods of time, quite possibly linked to differential predation mortality (mature females being much smaller and full of fat-rich gonads and eggs). 

Based on stomach sampling, Atlantic Halibut, Atlantic Wolffish, Thorny Skate, and other skate species are known predators of Snow Crab. Large Atlantic Halibut with mature female Snow Crab in their stomachs have been reported. Anecdotal information of some seals having fed upon Snow Crab are also known. 

Some of these predators (e.g., Halibut; DFO 2018) have significantly increased in abundance in the Region. However, for all, the abundance and encounter rates in areas overlapping with snow crab habitat what is more important, but this is not known. We do know from the bycatch in the Snow Crab survey that there are elevated areal **densities** with many snow crab trawl samples. This means that encounter rates will also likely increase and so too potentially predation mortality. However, high density does not equate to high abundance nor high total predation mortality; but this remains unknown and requires further analysis. The following presents information of areal density of co-occuring species; they are potential predators, competitors and prey. 

 
Atlantic Halibut densities have increased rapidly since 2010 (Figure \@ref(fig:halibut-timeseries); DFO 2018) on Snow Crab grounds. Most of these increases were towards The Gully, Slope Edge and near Sable Island (\@ref(fig:halibut-map)).

Thorny skate densities have been increasing as well (Figure \@ref(fig:thornyskate-timeseries)), especially in N-ENS and along the margins of Banquereau Bank (Figure \@ref(fig:thornyskate-map)). A minor decline in the densities have been seen in 4X.

Striped Atlantic Wolffish densities have been high, though declining in N-ENS since 2007 (Figure \@ref(fig:Wolffish-timeseries)). Highest densities were towards the Laurentian Channel (Figure \@ref(fig:Wolffish-map)).

Northern shrimp co-occur as they share similar habitat preferences and are also potential prey items of Snow Crab. Their numerical densities have declined after a peak in 2011, especially in S-ENS (Figures \@ref(fig:Shrimp-timeseries), \@ref(fig:Shrimp-map)). 

Lesser toad crab is a co-occurring species and potential competitor. Their numbers have declined to low levels throughout, after a peak in densities in 2007 and 2016 in N-ENS (Figures \@ref(fig:lessertoadcrab-timeseries), \@ref(fig:lessertoadcrab-map)).

Overall, higher predation mortality seems likely in N- and S-ENS and lower in 4X. Further co-occurring species have declined, possibly due to large-scaled habitat variations.



```{r halibut-timeseries, out.width='50%', echo=FALSE,   fig.align='center', fig.cap = 'Atlantic Halibut crude, unadjusted geometric mean numerical density (no/km$^2$) from annual Snow Crab survey. Error bars are 95\\%  Confidence Intervals.' }
include_graphics( file.path( data_loc, 'assessments', year_assessment, 'timeseries', 'survey', 'ms.no.30.png') )
# \@ref(fig:halibut-timeseries)
```
 

```{r halibut-map, out.width='45%', fig.show='hold', fig.align='center', fig.cap= 'Halibut density log10(no/km$^2$) from the Snow Crab survey.' }
loc = file.path( data_loc, 'output', 'maps', 'survey', 'snowcrab', 'annual', 'bycatch', 'ms.no.30' )
yrsplot = setdiff( year_assessment + c(0:-9), 2020)
fn4 = file.path( loc, paste( 'ms.no.30', yrsplot[4], 'png', sep='.') )
fn3 = file.path( loc, paste( 'ms.no.30', yrsplot[3], 'png', sep='.') )
fn2 = file.path( loc, paste( 'ms.no.30', yrsplot[2], 'png', sep='.') )
fn1 = file.path( loc, paste( 'ms.no.30', yrsplot[1], 'png', sep='.') )
include_graphics( c( fn2, fn1) )
# \@ref(fig:halibut-map)  
```


```{r thornyskate-timeseries, out.width='50%', echo=FALSE,  fig.align='center', fig.cap = 'Thorny Skate crude, unadjusted geometric mean numerical density (no/km$^2$) from annual Snow Crab survey. Error bars are 95\\%  Confidence Intervals.'}
include_graphics( file.path( data_loc, 'assessments', year_assessment, 'timeseries', 'survey', 'ms.no.201.png') )
# \@ref(fig:thornyskate-timeseries)
```

 

```{r thornyskate-map, out.width='45%', fig.show='hold', fig.align='center', fig.cap= 'Thorny skate density log10(no/km$^2$) from the Snow Crab survey.' }
loc = file.path( data_loc, 'output', 'maps', 'survey', 'snowcrab', 'annual', 'bycatch', 'ms.no.201' )
yrsplot = setdiff( year_assessment + c(0:-9), 2020)
fn4 = file.path( loc, paste( 'ms.no.201', yrsplot[4], 'png', sep='.') )
fn3 = file.path( loc, paste( 'ms.no.201', yrsplot[3], 'png', sep='.') )
fn2 = file.path( loc, paste( 'ms.no.201', yrsplot[2], 'png', sep='.') )
fn1 = file.path( loc, paste( 'ms.no.201', yrsplot[1], 'png', sep='.') )
include_graphics( c(  fn2, fn1) )
# \@ref(fig:thornyskate-map)  
```




```{r Wolffish-timeseries, out.width='50%', echo=FALSE,   fig.align='center', fig.cap = 'Striped Atlantic Wolffish crude, unadjusted geometric mean numerical density (no/km$^2$) from annual Snow Crab survey. Error bars are 95\\%  Confidence Intervals.' }
include_graphics( file.path( data_loc, 'assessments', year_assessment, 'timeseries', 'survey', 'ms.no.50.png') )
# \@ref(fig:Wolffish-timeseries)
```
 

```{r Wolffish-map, out.width='45%', fig.show='hold', fig.align='center', fig.cap= 'Striped Atlantic Wolffish density log10(no/km$^2$) from the Snow Crab survey.' }
loc = file.path( data_loc, 'output', 'maps', 'survey', 'snowcrab', 'annual', 'bycatch', 'ms.no.50' )
yrsplot = setdiff( year_assessment + c(0:-9), 2020)
fn4 = file.path( loc, paste( 'ms.no.50', yrsplot[4], 'png', sep='.') )
fn3 = file.path( loc, paste( 'ms.no.50', yrsplot[3], 'png', sep='.') )
fn2 = file.path( loc, paste( 'ms.no.50', yrsplot[2], 'png', sep='.') )
fn1 = file.path( loc, paste( 'ms.no.50', yrsplot[1], 'png', sep='.') )
include_graphics( c( fn2, fn1) )
# \@ref(fig:Wolffish-map)  
```


```{r Shrimp-timeseries, out.width='50%', echo=FALSE,   fig.align='center', fig.cap = 'Northern Shrimp crude, unadjusted geometric mean numerical density (n/$km^2$) from annual Snow Crab survey. Error bars are 95\\%  Confidence Intervals.' }
include_graphics( file.path( data_loc, 'assessments', year_assessment, 'timeseries', 'survey', 'ms.no.2211.png') )
# \@ref(fig:Shrimp-timeseries)
```
 

```{r Shrimp-map, out.width='45%', fig.show='hold', fig.align='center', fig.cap= 'Northern Shrimp density log10(no/km$^2$) from the Snow Crab survey.' }
loc = file.path( data_loc, 'output', 'maps', 'survey', 'snowcrab', 'annual', 'bycatch', 'ms.no.2211' )
yrsplot = setdiff( year_assessment + c(0:-9), 2020)
fn4 = file.path( loc, paste( 'ms.no.2211', yrsplot[4], 'png', sep='.') )
fn3 = file.path( loc, paste( 'ms.no.2211', yrsplot[3], 'png', sep='.') )
fn2 = file.path( loc, paste( 'ms.no.2211', yrsplot[2], 'png', sep='.') )
fn1 = file.path( loc, paste( 'ms.no.2211', yrsplot[1], 'png', sep='.') )
include_graphics( c( fn2, fn1) )
# \@ref(fig:Shrimp-map)  
```

 
```{r lessertoadcrab-timeseries, out.width='50%', echo=FALSE,   fig.align='center', fig.cap = 'Lesser Toad Crab crude, unadjusted geometric mean numerical density (no/km$^2$) from annual Snow Crab survey. Error bars are 95\\%  Confidence Intervals.' }
include_graphics( file.path( data_loc, 'assessments', year_assessment, 'timeseries', 'survey', 'ms.no.2521.png') )
# \@ref(fig:lessertoadcrab-timeseries)
```
 

```{r lessertoadcrab-map, out.width='45%', fig.show='hold', fig.align='center', fig.cap= 'Lesser Toad Crab density log10(no/km$^2$) from the Snow Crab survey.' }
loc = file.path( data_loc, 'output', 'maps', 'survey', 'snowcrab', 'annual', 'bycatch', 'ms.no.2521' )
yrsplot = setdiff( year_assessment + c(0:-9), 2020)
fn4 = file.path( loc, paste( 'ms.no.2521', yrsplot[4], 'png', sep='.') )
fn3 = file.path( loc, paste( 'ms.no.2521', yrsplot[3], 'png', sep='.') )
fn2 = file.path( loc, paste( 'ms.no.2521', yrsplot[2], 'png', sep='.') )
fn1 = file.path( loc, paste( 'ms.no.2521', yrsplot[1], 'png', sep='.') )
include_graphics( c(  fn2, fn1) )
# \@ref(fig:lessertoadcrab-map)  
```



## Other sources of uncertainty
 
* Bycatch of other species in the Snow Crab fishery cannot be reliably computed as they are derived from At-Sea-Observed data. 

* Bycatch of Snow Crab in other fisheries remains an area requiring attention. Anecdotal information from fishers suggest illegal retention of sublegal and female Snow Crab as bycatch and their use as bait. Illegal removals are also known to occur; however, the scale of such activities are not known.

* Oil and gas exploration (especially the use of seismics), development and exploitation (pollution) continues to be issues as their long-term effects remain unknown, though they are known to be generally deleterious to most life.

* Undersea, high-voltage cables is another source of concern. Many marine organisms are known to be sensitive to electromagnetic fields: fish in particular through their lateral line canals use it to detect prey. As such, predators and prey distributions can be affected causing indirect effects upon Snow Crab. Direct effects are also possible through trenching, but long term effects remain unknown.  

* Marine Protected Areas (MPAs) continue to be developed (e.g., Canada Gazette 2016). The presence of a refuge from fishing activities is potentially positive for Snow Crab. However, positive effects upon other organisms (predators or prey) can have counter-balancing effects. The overall long-term effects of the MPAs upon Snow Crab are unknown.

* Capture of soft-shell Snow Crab is always a concern. Prompt and careful return of immature (small-clawed, non-terminally molted) crab to the water is an important conservation measure that will enhance the 2-3 year productivity of the fishable component.

* Illegal, unreported, and unregulated fishing activities are known to occur. Such activities hinder the application of a precautionary approach to the management of this resource and cause potential bias and uncertainty in Reference Point estimation.

All of the above uncertainties are human induced and unlike the larger scaled climatic and ecosystemic uncertainties, there is some measure of human intervention possible. 

To remain adaptive in the face of these and other as yet unknown uncertainties including the ecosystemic, climatic uncertainties of which we are already aware, is the true challenge. It requires a balance of both resilience and robustness (see, Choi and Patten 2001). The Precautionary Approach as practiced by the snow crab fishers in Maritimes Region represents a unique model of such an adaptive approach. It values qualitity information and the communication and discussion of approaches in a distributed, collective information network that goes well beyond the simplistic and potentially maladaptive rubric of a purely "reference-points" based approach. Continued vigilence is necessary.


## CONCLUSIONS AND ADVICE

The ESS ecosystem is still experiencing a lot of volatility driven by rapid ecosystem and climatic variations. Under such conditions, it is prudent to be careful. Further, the overall indications of population status suggest that Snow Crab are still able to persist under extreme conditions if they are episodic, albeit, with some shifts in spatial distribution towards cooler and deeper waters. 

The modelled solutions represent a few of many possible views of the dynamics of snow crab. Over-emphasis of any one of these modelled solutions and associated Reference Points in determining a strategy for fisheries management is **not prudent** and certainly not precautionary. 


- In N-ENS, though recruitment continues at low levels, a gap in future recruitment to the fishery is expected for the next 1-3 years in N-ENS. N-ENS still exists in the healthy zone and so still has flexibility in approach. However, fishing mortality may have been higher than desirable. A more conservative harvest strategy would enable bridging this coming recruitment gap. A reduced TAC is prudent. 

- In S-ENS, recruitment to the fishery is likely to continue at a moderate rate for the upcoming season. The S-ENS stock remains in the healthy zone. Exploitation rates derived from the fishery model have been declining in recent years. S-ENS has considerable flexibility in approach. All TAC directions are viable. 

- In 4X, low to moderate levels of recruitment are expected for 2 years. 4X exists in the "cautious zone". The area is also in the southern-most extent of Snow Crab distribution in the North Atlantic and viable habitat has been depressed for many years. A reduced TAC is prudent. 




# REFERENCES

Canada Gazette. 2022. [Regulations Amending the Fishery (General) Regulations. Part II, Volume 156, Number 8.](https://www.gazette.gc.ca/rp-pr/p2/2022/2022-04-13/html/sor-dors73-eng.html)

Canada Gazette. 2016. St. Anns Bank Marine Protected Area Regulations. Canada Gazette, Part I, Vol 150, Issue 51: 4143-4149.
Choi, J.S. 2020. A Framework for the assessment of Snow Crab (*Chioneocete opilio*) in Maritimes Region (NAFO Div 4VWX) . DFO Can. Sci. Advis. Sec. Res. Doc. 2020/nnn. v + xxx p.

Choi, J.S. 2022. Reconstructing the Decline of Atlantic Cod with the Help of Environmental Variability in the Scotian Shelf of Canada. bioRxiv. https://doi.org/10.1101/2022.05.05.490753.

Choi, J. S., and B. C. Patten. 2001. Sustainable Development: Lessons from the Paradox of Enrichment. Ecosystem Health 7: 163–77.

Choi, Jae S., B. Cameron, K. Christie, A. Glass, and E. MacEachern. 2022. Temperature and Depth Dependence of the Spatial Distribution of Snow Crab. bioRxiv. https://doi.org/10.1101/2022.12.20.520893.

Choi, Jae S. 2023. A Multi-Stage, Delay Differential Model of Snow Crab Population Dynamics in the Scotian Shelf of Atlantic Canada. bioRxiv. https://doi.org/10.1101/2023.02.13.528296.

DFO. 2013. [Integrated Fisheries Management Plan for Eastern Nova Scotia and 4X Snow Crab (*Chionoecetes Opillio*.)](http://www.dfo-mpo.gc.ca/fm-gp/peches-fisheries/ifmp-gmp/snow-crab-neige/snow-crab-neiges2013-eng.htm)

DFO. 2018. Stock Status Update of Atlantic Halibut (Hippoglossus hippoglossus) on the Scotian Shelf and Southern Grand Banks in NAFO Divisions 3NOPs4VWX5Zc. DFO Can. Sci. Advis. Sec. Sci. Resp. 2018/022.

Hebert M, Miron G, Moriyasu M, Vienneau R, and DeGrace P. Efficiency and ghost fishing of Snow Crab (Chionoecetes opilio) traps in the Gulf of St. Lawrence. Fish Res. 2001; 52(3): 143-153. 10.1016/S0165-7836(00)00259-9   

```{r conclude, include=FALSE}
  detach(sn_env)
```
