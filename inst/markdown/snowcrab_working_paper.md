---
title: "Snow Crab status, Scotian Shelf"
subtitle: "Working Paper"

keywords: 
  - snow crab stock status assessment 
abstract: |
  Snow crab stock status assessment.

metadata-files:
  - _metadata.yml

params:
  year_assessment: 2025
  year_start: 1999
  data_loc:  "~/bio.data/bio.snowcrab"
  sens: 1
  debugging: FALSE
  model_variation: logistic_discrete_historical
  todo: [fishery_results,survey,ecosystem,fishery_model]

---


<!-- 

make quarto FN=snowcrab_working_paper.md DOCTYPE=html  PARAMS="-P year_assessment:2025 -P todo:[fishery_results,survey,ecosystem,fishery_model]"  --directory=~/bio/bio.snowcrab/inst/markdown
   

NOTE: this is not a full document but rather the basis for adding to a RES DOC ... just use for values and figures and tables

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
 
  
```


<!-- 
# _load_results.qmd contains instructions in "todo" 
#  this is a shared R-script to boot strap and provide a consistent data interface
-->

{{< include _load_results.qmd >}}  


# NOTE: this is just a template to get most of the results in one place and not a finalized document 

 
# ABSTRACT

Overall, fishing effort was much reduced in 2022, especially in S-ENS and 4X. TACs were mostly caught, with the exception of 4X. Catch rates were at or near all-time highs in all areas, except for 4X. At-Sea-Observed information was sparse and so fishery bycatch and snow crab size distributions were not evaluated.

Survey indices suggest a downtrend in numerical densities of Snow Crab recruiting to fishable sizes, especially in N-ENS. Similarly, declines in numerical densities of mature females, especially in N-ENS, suggest potential declines in egg production and long term recruitment. Mature snow crab sex ratios remain close to balanced levels, except again in N-ENS. Crude geometric mean biomass densities of the fishable component have also declined in all areas, and for S-ENS and 4X, it was at an historical low. Adjustment by spacetime index modelling suggest that the biomass index is also low, but habitat variability, especially as they relate to bottom temperature variability was very significant. 

Modelled inference (Model 1) closely follows the biomass index and suggests that the overall biomass of fishable Snow Crab has marginally declined over the last year and that fishing mortality has declined and are well below FMSY, with the exception of N-ENS. Model 1 places N- and S-ENS in the "healthy" zone and 4X in the "cautious" zone. All regions have flexibility. However, due to a lower level of recruitment in N-ENS and potentially elevated predation mortality, care should be taken in 2025 until recruitment returns.


# SUMMARY


- Fishing effort in `r year_assessment` were `r ey[["cfasnorth"]]` $\times 10^3$ trap hauls in N-ENS, `r ey[["cfasouth"]]` $\times 10^3$ trap hauls in S-ENS and  `r ey[["cfa4x"]]` $\times 10^3$ trap hauls in 4X. This represents a change of `r edt[["cfanorth"]]` %, `r edt[["cfasouth"]]` % and  `r edt[["cfa4x"]]` %, respectively, relative to the previous year. 

- Landings in `r year_assessment` were `r ly[["cfanorth"]]` t in N-ENS; `r ly[["cfasouth"]]` t in S-ENS; and `r ly[["cfa4x"]]` t in CFA 4X (season ongoing), representing a change of `r ldt[["cfanorth"]]`%, `r ldt[["cfasouth"]]`% and  `r ldt[["cfa4x"]]`%, respectively, relative to `r year_previous`. Total Allowable Catches (TACs) for `r year_assessment` were `r tacy[["cfanorth"]]` t, `r tacy[["cfasouth"]]` t and `r tacy[["cfa4x"]]` t in N-ENS, S-ENS and 4X, respectively. 

- Non-standardized fishery catch rates in `r year_assessment` were `r cy[["cfanorth"]]`, `r cy[["cfasouth"]]fasouth"]]`  and `r cy[["cfa4x"]]` kg/trap haul in N-ENS, S-ENS and 4X, respectively. This represents a change of respectively, `r cdt[["cfanorth"]]` %, `r cdt[["cfasouth"]]` % and  `r cdt[["cfa4x"]]` % (season ongoing) relative to the previous year. Though the spatial extent of exploitation was dispersed, with many of the exploited areas showing lower catch rates.

- Commercial catches of soft-shelled (newly moulted) Snow Crab for `r year_assessment` were `r ccsy[["cfanorth"]]rth"]]rth"]]`%, `r ccsy[["cfasouth"]]`% and `r ccsy[["cfa4x"]]`%, respectively, in N-ENS S-ENS and 4X (season ongoing). In `r year_previous`, it was `r ccsp[["cfanorth"]]`%, `r ccsp[["cfasouth"]]`% and `r ccsp[["cfa4x"]]`%, respectively. Sensecent (CC5) crab levels were higher in N- and S-ENS than in previous years, though low and inconsistent sampling effort makes this uncertain.

- Bycatch of non-target species continuesto be low (<< 1% of total catch) in all Snow Crab fishing areas except 4X which was not observed. 

- Little to no recruitment is expected for the next 1-3 years in N-ENS. In S-ENS, continued recruitment is expected. In 4X, low levels of recruitment are expected. Egg and larval production is expected to be high in the next year in all areas except N-ENS. 

- In N-ENS, the modelled biomass (pre-fishery) of Snow Crab in `r year_assessment` was `r round(B_north[t0], 2)` t, relative to `r round(B_north[t1], 2)` t in `r year_previous`. In S-ENS, the modelled biomass (pre-fishery) was `r round(B_south[t0], 2)` t, relative to `r round(B_south[t1], 2)` t in `r year_previous`. In 4X, the `r year_assessment`-`r year_assessment+1` season's modelled biomass (pre-fishery) was `r round(B_4x[t0], 2)` t, relative to `r round(B_4x[t1], 2)` t in the `r year_previous`-`r year_assessment` season. 

- Average bottom temperatures observed in the Snow Crab survey have declined to historical lows since the historical highs in 2022. Average viable habitat surface area has however declined to near historical lows since peaking in the early 2010s.

- In N-ENS, though recruitment continues at low levels, a gap in future recruitment to the fishery is expected for the next 1-3 years in N-ENS. N-ENS still exists in the healthy zone and so still has flexibility in approach. However, fishing mortality may have been higher than desirable. A more conservative harvest strategy would enable bridging this coming recruitment gap. A reduced TAC is prudent. 

- In S-ENS, recruitment to the fishery is likely to continue at a moderate rate for the upcoming season. The S-ENS stock remains in the healthy zone. Exploitation rates derived from the fishery model have been declining in recent years. S-ENS has considerable flexibility in approach. All TAC directions are viable. 

- In 4X, low to moderate levels of recruitment are expected for 2 years. 4X exists in the "cautious zone". The area is also in the southern-most extent of Snow Crab distribution in the North Atlantic and viable habitat has been depressed for many years. A reduced TAC is prudent. 




## Terms of reference 

To provide an assessment of the Maritimes Region Snow crab fishery (Figure -@fig-management-areas):

  1. Status of Snow Crabs
  2. Relative abundance and exploitation rates
  3. Evaluate consequences of different harvest levels upon abundance and exploitation rates
  4. Bycatch in survey

Amendments of the Fisheries Act in 2022 (*Fish Stock Provisions*), have encoded management approaches and requirements relative to biological reference points. 

The above are addressed through directed surveys funded by a Collaborative Agreement with the Snow Crab Fishers in the Maritimes Region of Canada.


## Fishery performance

- A summary for each main fishing area is shown in -@tbl-fishery-performance
- Biological data of fished component obtained from At-Sea-Observer data base (figure -@fig-map-observer-locations)
- Catch data obtained from dockside monitoring of fishery logbooks in the Marfis data base (100% coverage)

- In `r year_assessment`, they were `r ly[["cfanorth"]]`, `r ly[["cfasouth"]]` and `r ly[["cfa4x"]]` t, in N-ENS, S-ENS and 4X (season ongoing), respectively. 

- Relative to `r year_previous`, they represent changes of `r ldt[["cfanorth"]]`%, `r ldt[["cfasouth"]]`% and  `r ldt[["cfa4x"]]`%, respectively. 

- Total Allowable Catches (TACs) for `r year_assessment` were `r tacy[["cfanorth"]]` t, `r tacy[["cfasouth"]]` t and `r tacy[["cfa4x"]]` t in N-ENS, S-ENS and 4X, respectively. 

Non-standardized fishery catch rates in `r year_assessment` were `r cy[["cfanorth"]]`, `r cy[["cfasouth"]]`  and `r cy[["cfa4x"]]` kg/trap haul in N-ENS, S-ENS and 4X, respectively; `r cdt[["cfanorth"]]` %, `r cdt[["cfasouth"]]` % and  `r cdt[["cfa4x"]]` % (season ongoing) relative to the previous year 

## Effort

- N-ENS: increased in 2024 relative to the previous year and spatially more dispersed. 
- S-ENS: increased in 2024 relative to the previous year and spatially more dispersed.
- CFA 4X: decreased in 2024-2025 relative to the previous year and spatially contracted.
- Reference: Figures -@fig-effort-timeseries,  -@fig-effort-map
 

## Landings

- N-ENS: declined marginally and spatially more dispersed 
- S-ENS: declined marginally and spatially more dispersed 
- CFA 4X: only a small fraction of TACs landed and fishery is ongoing
- Reference: Figure -@fig-landings-timeseries, -@fig-landings-map


## Catch rates

- N-ENS: decreased to 2019 levels, declines inshore
- S-ENS: decreased to 2017 levels, declines through, especially CFA 24
- CFA 4X: increased to 2021 levels, but only a small fraction of TACs landed in Sambro area; fishery is ongoing
- Reference: -@fig-cpue-timeseries, -@fig-cpue-map
    

## Discard of snow crab

- N-ENS: increased since 2022, but within historical range; coverage was 1.3% of landings
- S-ENS: increased since 2022, but within historical range; coverage was 2.4% of landings
- CFA 4X: No data, fishery is ongoing 
- Reference: -@fig-discard_maritimes


## Discard of soft-shelf crab

- Generally, higher soft-shell indicates incoming recruitment to the fishery and their handling potential and unnecessary handling/discard mortality.

- Commercial catches of soft-shelled crab were `r ccsy[["cfanorth"]]`% (low sampling), `r ccsy[["cfasouth"]]`% (low sampling) and `r ccsy[["cfa4x"]]`% (no sampling; season ongoing) in N-ENS, S-ENS and 4X, respectively for `r year_assessment`. 

- In `r year_previous`, it was `r ccsp[["cfanorth"]]`% (no sampling), `r ccsp[["cfasouth"]]`% (low sampling) and `r ccsp[["cfa4x"]]`% (no sampling), respectively. 


- N-ENS: increased since 2022, but within historical range
- S-ENS: increased since 2022, but within historical range
- CFA 4X: No data, fishery is ongoing
- Reference: -@fig-observed-softshell-map, -@tbl-observer-softgt95
 

## Discard of non-target species ("Bycatch")

- N-ENS: average of 0.02% of landings; primarily other Crustacea (crab and lobster)  
- S-ENS: average of 0.03% of landings; primarily other Crustacea (crab and lobster)
- CFA 4X: average of 0.87% of landings; primarily other Crustacea (crab and lobster)  
- Reference: -@fig-map-observer-locations, -@tbl-fishery-discard-effort


## Bitter crab disease

- Low level background infection rate of < 0.1% of the surveyed crab (-@tbl-bcd)
- Spatial distribution is widespread and usually in shallower locations (-@fig-bcd-map) 

## Survey indices

- Details of survey design and index modelling can be found in (Christie et al. 2024; Choi 2020;  -@fig-survey-locations-map )
- Numerous vessel and captain changes have occured over the years
- 2004: transition from a Spring to Fall survey and mandate transfer from Moncton to BIO
- 2020: due to Covid-19 related to health and safety uncertainties, survey was not conducted 
- 2022: due to mechanical issues with the survey vessel, inshore areas of S-ENS were partially sampled  
- 2024: due to financial issues, areas associated with Marine Protected Areas (St Ann's, Gully) were not sampled
  - N-ENS: 58 stations completed 
  - S-ENS: 282 stations completed
  - CFA 4X: 24 stations completed



## Growth patterns inferred from size modes

- Reference: -@fig-growth-modes


## Recruitment: males

- N-ENS: recruitment gap noted in 2021 has been bridged (year-class centered at 68 mm CW in 2023 in now at 90 mm CW) but subsequent year-classes are weak suggesting high natural mortality (predation?); little internal recruitment is expected for the next 3-4 years
- S-ENS: strong and stable size structure suggestive of a stable size/age structure; recruitment expected in 2025 and elevated soft-shell capture will require care
- CFA 4X: strong years classes (30 mm CW in 2019) were progressing but in 2024 has been lost due to elevated natural morality (predation and high temperatures?); minimal internal recruitment expected 
- Reference: -@fig-sizefeq-male
  

## Recruitment: females

- N-ENS: mature fraction has senesced with minimal recruitment into mature phase; next large year class is at 25 mm CW
- S-ENS: strong and stable mature fraction; steady recruitment into mature phase expected in 2025  
- CFA 4X: mature fraction has been strong for much of the historical record; recruitment into mature phase is expected to be low and the female population to senesce
- Reference: -@fig-sizefeq-female
 

## Reproduction

- N-ENS: densities (~potential egg production) has been low since 2021 and near historical lows; with high densities in inshore and Glace Bay hole
- S-ENS: densities (~potential egg production) declined marginally, especially offshore 
- CFA 4X: densities (~potential egg production) has continued to decline since 2021 and is near historial lows;highest densities close to Sambro and Lunenburg
- Reference: -@fig-totno-female-mat-timeseries, -@fig-totno-female-mat-map

 
## Sex ratios

- N-ENS: sex ratio of mature crab are female have increased marginally; with increases in near shore and in shallower areas
- S-ENS: sex ratios have declined marginally; with declines occurring in Missaine and inshore 24
- CFA 4X: sex ratios have declned but still high; peaks in roseway and further downstream
- Reference: -@fig-sexratio-mat-timeseries, -@fig-sexratio-mat-map

The sex ratios (proportion female) of the mature component is particularly important as locally unbalanced ratios can impact upon encounter rates and ultimately reproductive success. In the ESS, there is generally a lack of females, in contrast to, for example, the Gulf of St-Lawrence where the reverse is often the case. Higher sex ratios are usually found in inshore and bottom slope areas. A decline in sex ratios has been observed since 2017 in N-ENS. In S-ENS the sex ratio increased from 20% in 2021 to just under 35% in 2022. In 4X, mature sex ratios are more balanced and currently near the 50% level.   



## Fishable biomass density

- N-ENS: densities have declined; declines mostly inshore
- S-ENS: densities have increased; increases mostly offshore
- CFA 4X: densities have declined to historical lows; peak only in Mahonr Bay area (south of Lunenburg)
- Reference: -@fig-R0-timeseries, -@fig-R0-map


The fishable component is defined as Snow Crab that are male, mature, and larger than 95 mm CW. Fishable biomass density is the geometric mean biomass per unit swept area by the trawl. A peak in biomass density was observed in 2009 to 2014 and has since been declining in all areas. Note that high and low biomass density areas fluctuate with time. Biomass density, however, does not equate to total biomass as the areas occupied by crab can contract, expand and shift with environmental conditions and ecosystem variability. 


## Bottom Temperatures

- N-ENS: declined to historical lows since historical highs in 2022
- S-ENS: declined to historical lows since historical highs in 2022
- CFA 4X: declined to historical lows since historical highs in 2022
- Reference: -@fig-bottom-temperatures-timeseries, -@fig-figures-temperature-bottom-map

Average bottom temperatures **observed** in the 2022 Snow Crab survey were near or above historical highs in all areas. A general warming trend has been observed in the Snow Crab survey since the early 1990s on the Scotian Shelf (Choi et al. 2022). Temperatures are more stable in N-ENS than S-ENS; 4X exhibits the most erratic and highest annual mean bottom temperatures. Of particular note, the observed temperatures in the 2022 Snow Crab survey for S-ENS increased well above the average. This is expected to be a source of bias and certainty for S-ENS abundance predictions as it is outside the range normally encountered by the statistical models. Furthermore, the Groundfish surveys did not operate over Snow Crab grounds in 2020 and 2022 and so temperature data is also very sparse for the area of interest.

Upon aggregation and modelling of historical temperature data, the average temperature is found to have increased well beyond the $7^\circ$C threshold in 4X. N-ENS and S-ENS also continued to experience historical highs in bottom temperature and elevated spatial variability of bottom temperatures (). In particular, the spike observed in S-ENS (above) was less extreme and so likely due to a short term influx of warm water. 

Overall, since 1999, there has been a persistent spatial gradient of almost $15^\circ$C in bottom temperatures in the Maritimes Region. This large gradient in a spatially complex and temporally dynamic area makes assessment a particular challenge for a stenothermic organism such as Snow Crab (Choi et al. 2022). 


## Interspecific interactions

- Predators, preys, competitors

  - N-ENS:
  - S-ENS:
  - CFA 4X:
  - Reference: 

Being long-lived, the influence of predators can be significant. Especially important are predators of the smaller immature and female Snow Crab. Increasing predation not only lowers the abundance and recruitment, it can also reduce the reproductive potential of Snow Crab and therefore long-term population cycles. N-ENS and S-ENS are well known to have skewed sex ratios with few mature females for extended periods of time, quite possibly linked to differential predation mortality (mature females being much smaller and full of fat-rich gonads and eggs). 

Based on stomach sampling, Atlantic Halibut, Atlantic Wolffish, Thorny Skate, and other skate species are known predators of Snow Crab. Large Atlantic Halibut with mature female Snow Crab in their stomachs have been reported. Anecdotal information of some seals having fed upon Snow Crab are also known. 

Some of these predators (e.g., Halibut; DFO 2018) have significantly increased in abundance in the Region. However, for all, the abundance and encounter rates in areas overlapping with snow crab habitat is more important, but this is not known. We do know from the bycatch in the Snow Crab survey that there are elevated areal **densities** with many snow crab trawl samples. This means that encounter rates will also likely increase and so too potentially predation mortality. However, high density does not equate to high abundance nor high total predation mortality; but this remains unknown and requires further analysis. The following presents information of areal density of co-occuring species; they are potential predators, competitors and prey. 

Atlantic Halibut densities have increased rapidly since 2010 (-@fig-halibut-timeseries; DFO 2018) on Snow Crab grounds. Most of these increases were towards The Gully, Slope Edge and near Sable Island (-@fig-halibut-map).

Thorny skate densities have been increasing as well (-@fig-thornyskate-timeseries), especially in N-ENS and along the margins of Banquereau Bank (-@fig-thornyskate-map). A minor decline in the densities have been seen in 4X.

Striped Atlantic Wolffish densities have been high, though declining in N-ENS since 2007 (-@fig-stripatlwolffish-timeseries). Highest densities were towards the Laurentian Channel (-@fig-stripatlwolffish-map).

Northern shrimp co-occur as they share similar habitat preferences and are also potential prey items of Snow Crab. Their numerical densities have declined after a peak in 2011, especially in S-ENS (-@fig-northernshrimp-timeseries, -@fig-northernshrimp-map). 

- Lesser toad crab is a co-occurring species and potential competitor. Their numbers have declined to low levels throughout, after a peak in densities in 2007 and 2016 in N-ENS (-@fig-lyrecrab-timeseries, -@fig-lyrecrab-map).

- Overall, higher predation mortality seems likely in N- and S-ENS and lower in 4X. Further shrimp with similar habitat preferences have declined, possibly due to large-scaled habitat variations and predation.





## Viable habitat

- Snow Crab being cold water stenotherms, stability of environmental conditions is critical for their survival. However, it is not just temperature that is in control, but rather, where the temperature is good in combination with substrate, depths, co-occuring species, etc. The Maritimes Region being at the confluence of many oceanic currents renders the area highly variable. Rapid specis composition change and climate change and uncertainty exacerbates this situation. The viable habitat estimated for each area across time has shown some variations in the historical record. 

- 4X showed a significantly lower average viable habitat levels relative to the N-ENS and S-ENS.  

- A peak in average probability of observing fishable snow crab ("viable habitat") was observed in 2010 for 4X, 2011 for N-ENS and 2012 for S-ENS. Since 2015, the average viable habitat has declined to historical lows and remained so in 2022. 

- N-ENS: Pr has declined to near historical lows (peaks in 2011 and 2020); peaks in inner trench and GBH
- S-ENS: Pr has declined to near historical lows (peaks in 2012); peaks near Sable and missaine
- CFA 4X: declined to near historical lows (peaks in 2010); peak near Lunenburg
- Reference: -@fig-fb-habitat-timeseries, -@fig-fb-habitat-map

 

## Fishable biomass index

- N-ENS: decline to historical lows; decline throughout area
- S-ENS: increased marginally; increase offshore (marginally, north of sable)
- CFA 4X: increased marginally; area south of sambro and lunenburg
- Reference:  -@fig-fbindex-timeseries, -@fig-fbindex-map

The fishable **biomass index** ( ) are optimistically high as the spatial expansion uses areal units with large surface areas, larger than the patchiness of Snow Crab distributions (spatial autocorrelation length is <20 km, on average; Choi 2020). As such, it should only be seen as a spatially and temporally comparable relative index of abundance.  

The spatial distribution of the biomass index has been consistent over the past six years, with a peak in overall biomass index in 2019 and 2020. Since then, a  reduction was observed throughout the region, with the exception of the core areas. A contraction of spatial range in 4X and the western parts of S-ENS were also evident in 2021 to 2022. Upon aggregation, the biomass index declined marginally in all areas  


## Posterior estimates of fishable biomass

- N-ENS: marginal increase in prefishery biomass
- S-ENS: marginal increase in prefishery biomass
- CFA 4X: marginal increase in prefishery biomass
- Reference: -@fig-logisticPredictions

The biomass index along with fishery removals are used to fit a *Logistic Biomass Dynamics Model* to determine fishable **modelled biomass** (biomass estimated from the fisheries model) and relevant biological reference points (i.e., carrying capacity and fishing mortality at maximum sustainable yield, or F~MSY~). In N-ENS, the modelled biomass (pre-fishery) of Snow Crab in `r year_assessment` was `r round(B_north[t0], 2)` kt, relative to `r round(B_north[t1], 2)` kt in `r year_previous`. In S-ENS, the `r year_assessment` modelled biomass (pre-fishery) was `r round(B_south[t0], 2)` kt, relative to `r round(B_south[t1], 2)` kt in `r year_previous`. In 4X, the modelled biomass (pre-fishery) for the `r year_assessment`-`r year_assessment+1` season was `r round(B_4x[t0], 2)` kt, relative to `r round(B_4x[t1], 2)` kt for the `r year_previous`-`r year_assessment` season. 
  

## Posterior estimates of fishing mortality

- N-ENS: increased
- S-ENS: decreased
- CFA 4X: decreased
- Reference: -@fig-logisticFishingMortality 

In N-ENS, the `r year_assessment` fishing mortality is estimated to have been `r round(FM_north[t0],3)` (annual exploitation rate of `r round(100*(exp(FM_north[t0])-1),2)`%), up from the  `r year_previous` rate of `r round(FM_north[t1],3)` (annual exploitation rate of `r round(100*(exp(FM_north[t1])-1),1)`%.

In S-ENS, the `r year_assessment` fishing mortality is estimated to have been `r round(FM_south[t0],3)` (annual exploitation rate of `r round(100*(exp(FM_south[t0])-1),1)`%), decreasing marginally from the `r year_previous` rate of `r round(FM_south[t1],3)` (annual exploitation rate of `r round(100*(exp(FM_south[t1])-1),1)`%. Localized exploitation rates are likely higher, as not all areas for which biomass is estimated are fished (e.g., continental slope areas and western, inshore areas of CFA 24).

In 4X, the `r year_assessment`-`r year_assessment+1` season (ongoing), fishing mortality is estimated to have been `r round(FM_4x[t0],3)` (annual exploitation rate of `r round(100*(exp(FM_4x[t0])-1),1)`%), decreasing from the `r year_assessment-1`-`r year_assessment` season rate of `r round(FM_4x[t1],3)` (annual exploitation rate of `r round(100*(exp(FM_4x[t1])-1),1)`%. Localized exploitation rates are likely higher, as not all areas for which biomass is estimated are fished. 


## Reference Points

- N-ENS: In cautious zone; FM is high; 
- S-ENS: In healthy zone
- CFA 4X: In critical zone
- Reference: -@fig-logistic-hcr

Reference points are used to guide harvest strategies (Canada Gazette 2022; DFO 2013). More specifically, the state of the fishery relative to the following "Reference Points" are assessed. In terms of modelled biomass, Lower and Upper Stock Reference are 25% and 50% of carrying capacity which delineate "critical", "cautious" and "healthy" zones. In terms of exploitation rates, the Upper Removal Reference is the exploitation rate that we try not to go beyond; it is defined in term of the fishing mortality associated with Maximum Sustainable Yield (FMSY). In the biomass dynamics model, FMSY = $r /2$. As $r \approx 1$ for snow crab, FMSY $\approx 0.5$ is expected. 

The operational target exploitation changes depending upon the "zone" in which a population lands. When in the "healthy" zone, the rule of thumb has been to keep annual exploitation rates between 10% to 32% of the available biomass ($F = 0.11, 0.36$, respectively). In the "cautious" zone, the rule of thumb has been to keep annual exploitation rates between 0% to 20% ($F = 0, 0.22$, respectively). In the "critical" zone, fishery closure is considered until recovery is observed, where recovery indicates at a minimum, modelled biomass > LSR. Other biological and ecosystem considerations such as recruitment, spawning stock (female) biomass, size structure, sex ratios and environmental and ecosystem conditions, provide additional guidance within each range.

Model 1 estimates key Reference Points as shown in -@tbl-reference-points and -@fig-ReferencePoints. The related PA thresholds can be computed as:

  - Lower Stock Reference (LSR): $K/4$
  - Upper Stock Reference (USR): $K/2$
  - Upper Removal Reference (URR): keep fishing mortality below *FMSY* = $r/2$

Model 1 suggests the current state of the fishable components to be (-@fig-logistic-hcr):  

  - N-ENS is in the "healthy" zone
  - S-ENS is in the "healthy" zone
  - 4X is in the "cautious" zone

It should be emphasized that using these parameters assumes that the population dynamics are well described by the fishery model. This is, of course, **not true**. For example, the observation of fisheries landings is assumed to be known without error. This is not true as illegal and unreported exploitation occurs. These and other unaccounted factors (recruitment strength, environmental variability, predation intensity, disease) can easily bias parameter estimates. As such, **caution is required in using these reference points.** Other contextual reference points must be used in conjunction:

  - Strength of recruitment (short-term, long-term)
  - Strength of spawning stock (females)
  - Ecosystem variability (predator and prey trends and distributions) within norms
  - Habitat viability within norms
  - Availability of spatial and temporal refugia within norms 

We turn to some these additional factors in the next section.

 

## Other sources of uncertainty
 
- Capture of soft-shell Snow Crab is always a concern. Prompt and careful return of immature (small-clawed, non-terminally molted) crab to the water is an important conservation measure that will enhance the 2-3 year productivity of the fishable component.

- Bycatch of Snow Crab in other fisheries remains an area requiring attention. Anecdotal information from fishers suggest illegal retention of sublegal and female Snow Crab as bycatch and their use as bait. Illegal removals are also known to occur; however, the scale of such activities are not known.

- Illegal, unreported, and unregulated fishing activities are known to occur. Such activities hinder the application of a precautionary approach to the management of this resource and cause potential bias and uncertainty in Reference Point estimation.

- Marine Protected Areas (MPAs) continue to be developed (e.g., Canada Gazette 2016). The presence of a refuge from fishing activities is potentially a direct positive effect upon Snow Crab. However, positive effects upon other organisms (predators or prey) can have counter-balancing effects. The overall long-term effects of the MPAs upon Snow Crab are unknown.



# Conclusions  

- N-ENS:
- S-ENS:
- CFA 4X:
- Reference: 

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





## Tables
 
```{r}
#| label: tbl-fishery-performance
#| echo: false
#| results: asis
#| eval: true
#| output: true
#| tbl-cap: "Fishery performance statistics. Note: 4X years represent the starting year and currently ongoing."
#| tbl-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

for (r in 1:nregions) {
  reg = regions[r]
  REG = reg_labels[r]
  cat(REG, "\n")
  oo = dt[ which(dt$Region==reg), c("Year", "Licenses", "TAC", "Landings", "Effort", "CPUE")] 

  names(oo) = c( "Year", "Licenses", "TAC (t)", "Landings (t)", "Effort (1000 th)", "CPUE (kg/th)" )

  out = gt::gt(oo) |> gt::tab_options(table.font.size = 14, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
  print(out)
  cat("\n\n")
} 

```
 
&nbsp;  $~$  <br /> 


```{r}
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

odb = odb0[ cw >= 95 & cw < 170  & prodcd_id==0 & shell %in% c(1:5) & get(vnr) %in% regions & sex==0, ]  # male
shell_condition = odb[ !is.na(get(vnr)), .N, by=c(vnr, "fishyr", "shell") ]
shell_condition[, total:=sum(N, na.rm=TRUE), by=c(vnr, "fishyr")]
shell_condition$percent = round(shell_condition$N / shell_condition$total, 3) * 100
shell_condition$Year = shell_condition$fishyr
 
for (r in 1:nregions){
  reg = regions[r]
  REG = reg_labels[r]
  cat( REG, "\n")

  soft  = odb[ get(vnr)==reg & durometer <  68, .(Soft=.N), by=.(fishyr ) ] 
  total = odb[ get(vnr)==reg & is.finite(durometer) , .(Total=.N), by=.(fishyr) ] 
  oo = soft[total, on="fishyr"]
  oo = oo[, .(Year=fishyr, Soft=round(Soft/Total*100, 1), Total=Total) ]  
  scond = shell_condition[ get(vnr)==reg & shell %in% c(1,2), .(SoftSC=sum(percent), TotalSC=unique(total)[1]), by=.(Year)]
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

&nbsp;  $~$  <br /> 


```{r}
#| label: tbl-fishery-discard-effort
#| eval: true
#| output: true
#| echo: false
#| results: asis
#| layout-ncol: 1
#| tbl-cap: "Bycatch (kg) estimated from fisheries effort. Dots indicate low values. Where species exist in a list but there is no data, this indicates some historical bycatch. The overall average is only for the years shown." 
#| tbl-subcap: 
#|   - "(a)"
#|   - "(b)"
#|   - "(c)"
#|   - "(d)"

for (r in 1:nregions){
  reg = regions[r]
  REG = reg_labels[r]
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

&nbsp;  $~$  <br /> 


 

```{r}
#| label: tbl-bcd
#| echo: false
#| eval: true
#| output: true
#| tbl-cap: "[Bitter Crab Disease in Maritimes Region](https://www.dfo-mpo.gc.ca/science/aah-saa/diseases-maladies/hematcb-eng.html) is an infection from a dinoflagellate (*Hematodinium*) that causes muscle degeneration. They are widespread (Alaska, NW Atlantic, Greenland) and usually found in warm-water, physiologically stressful conditions. In the Maritimes, it seems to be a low level background infection, found everywhere in the fishing grounds."
#| fig.show: hold 
#| fig-dpi: 144
#| fig-height: 10

include_graphics( file.path( data_loc, "output", "bcd.png") )

```
 

&nbsp;  $~$  <br /> 


|   | N-ENS | S-ENS | 4X |
|----- | ----- | ----- | ----- |
| |  |  |  |
|q       | `r round(q_north, 3)` (`r round(q_north_sd, 3)`) | `r round(q_south, 3)` (`r round(q_south_sd, 3)`) | `r round(q_4x, 3)` (`r round(q_4x_sd, 3)`) |
|r       | `r round(r_north, 3)` (`r round(r_north_sd, 3)`) | `r round(r_south, 3)` (`r round(r_south_sd, 3)`) | `r round(r_4x, 3)` (`r round(r_4x_sd, 3)`) |
|K       | `r round(K_north, 2)` (`r round(K_north_sd, 2)`) | `r round(K_south, 2)` (`r round(K_south_sd, 2)`) | `r round(K_4x, 2)` (`r round(K_4x_sd, 2)`) |
|Prefishery Biomass   | `r round(B_north[t0], 2)` (`r round(B_north_sd[t0], 2)`) | `r round(B_south[t0], 2)`  (`r round(B_south_sd[t0], 2)`) | `r round(B_4x[t0], 2)`  (`r round(B_4x_sd[t0], 2)`)  |
|Fishing Mortality    | `r round(FM_north[t0], 3)` (`r round(FM_north_sd[t0], 3)`) | `r round(FM_south[t0], 3)` (`r round(FM_south_sd[t0], 3)`) | `r round(FM_4x[t0], 3)` (`r round(FM_4x_sd[t0], 3)`) |

: Reference points from the logistic biomass dynamics fishery model. K is Carrying capacity (kt); and r is Intrinsic rate of increase (non-dimensional). Note that FMSY (fishing mortality associated with 'Maximum Sustainable Yield') is r/2. Similarly, BMSY (biomass associated with 'Maximum Sustainable Yield') is K/2. SD is posterior Standard deviations.* {#tbl-reference-points}

&nbsp;  $~$  <br /> 


```{r}
#| label: tbl-predators
#| echo: false
#| eval: true
#| output: true
#| tbl-cap: "Main predators based upon frequency of occuence of snow crab in finfish stomach samples, unadjusted for sampling effort."

gt::gt(counts[1:11,]) 

```



## Figures


```{r}
#| label: fig-management-areas
#| eval: true
#| echo: false 
#| output: true
#| fig-cap: "The Scotian Shelf (NW Atlantic Ocean; NAFO Div. 4VWX). Shown are isobaths and major bathymetric features. Managed Crab Fishing Areas (CFAs; divided by dashed lines) include: NENS, SENS, 4X. SENS is further subdivided (dotted line) into 23 (NW) and 24 (SE)."
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold 

fn = file.path( media_loc, "snowcrab_cfas.png" )
knitr::include_graphics( fn ) 

```
 

 
&nbsp;  $~$  <br /> 


```{r}
#| label: fig-map-observer-locations
#| eval: true 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| echo: false 
#| layout-ncol: 2
#| fig-cap: "Snow Crab At-sea-observer locations."
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

loc = file.path( data_loc, "output", "maps", "observer.locations" )
yrsplot = year_assessment + c(0:-3) 
 
fns = paste( "observer.locations", yrsplot, "png", sep="." ) 
fn = file.path( loc, fns ) 

include_graphics( fn )

```
   


```{r}
#| label: fig-effort-timeseries
#| eval: true 
#| output: true
#| fig-dpi: 144
#| fig-height: 8
#| fig-cap: "Temporal variations in fishing effort."
 
if (params$sens==1) {
  ts_dir = file.path( data_loc, "assessments", year_assessment, "timeseries", "fishery" )
} else if (params$sens==2) {
  ts_dir = file.path( data_loc, "assessments", year_assessment, "timeseries", "fishery", "split" )
}

include_graphics( file.path( ts_dir, "effort.ts.png" ) )
 

```

&nbsp;  $~$  <br /> 


```{r}
#| label: fig-effort-map
#| eval: true 
#| output: true
#| fig-dpi: 144
#| fig-height: 8
#| echo: false 
#| layout-ncol: 2
#| fig-cap: "Snow Crab fishing effort from fisheries logbook data for previous and current years (no X 10$^3$ per 10 km X 10 km grid)."
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""
 
loc0 = file.path( data_loc, "output", "maps", "logbook", "snowcrab", "annual", "effort" )
yrsplot = year_assessment + c(0:-3)
fn = file.path( loc0, paste( "effort", yrsplot, "png", sep=".") ) 

include_graphics( fn ) 

```
 
&nbsp;  $~$  <br /> 

```{r}
#| label: fig-landings-timeseries
#| eval: true 
#| output: true
#| fig-cap: "Landings (t) of Snow Crab on the SSE. For 4X, the year refers to the starting year of the season. Inset is a closeup view of the timeseries for N-ENS and 4X."
#| fig-dpi: 144
#| fig-height: 8

if (params$sens==1) {
  ts_dir = file.path( data_loc, "assessments", year_assessment, "timeseries", "fishery" )
} else if (params$sens==2) {
  ts_dir = file.path( data_loc, "assessments", year_assessment, "timeseries", "fishery", "split" )
}

include_graphics( file.path( ts_dir, "landings.ts.png" ) )
```

&nbsp;  $~$  <br /> 


```{r}
#| label: fig-landings-map
#| eval: true 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| echo: false 
#| layout-ncol: 2
#| fig-cap: "Snow Crab landings from fisheries logbook data for previous and current years (tons per 10 km x 10 km grid)."
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

loc0 = file.path( data_loc, "output", "maps", "logbook", "snowcrab", "annual", "landings" )
yrsplot = year_assessment + c(0:-3)
fn = file.path( loc0, paste( "landings", yrsplot, "png", sep=".") ) 

include_graphics( fn ) 

```


&nbsp;  $~$  <br /> 


```{r}
#| label: fig-cpue-timeseries
#| eval: true 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig-cap: "Temporal variations in crude catch rates of Snow Crab (kg/trap haul)."
 
if (params$sens==1) {
  ts_dir = file.path( data_loc, "assessments", year_assessment, "timeseries", "fishery" )
} else if (params$sens==2) {
  ts_dir = file.path( data_loc, "assessments", year_assessment, "timeseries", "fishery", "split" )
}

include_graphics( file.path( ts_dir, "cpue.ts.png" ) ) 
```

&nbsp;  $~$  <br /> 


```{r}
#| label: fig-cpue-map
#| eval: true 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| echo: false 
#| layout-ncol: 2
#| fig-cap: "Snow Crab crude catch rates on the Scotian Shelf for previous and current years. Units are kg/trap haul per 10 km x 10 km grid."
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

loc0= file.path( data_loc, "output", "maps", "logbook", "snowcrab", "annual", "cpue" )
yrsplot = year_assessment + c(0:-3)
fn = file.path( loc0, paste( "cpue", yrsplot, "png", sep=".") ) 

include_graphics( fn ) 
 
```  

&nbsp;  $~$  <br /> 



```{r}
#| label: fig-discard_maritimes
#| fig-cap: "At sea observed rate of snow crab discards relative to total catch (discard + kept), all Maritimes."
#| fig-dpi: 144
#| fig-height: 6
 
oo = NULL
for (reg in regions) {
  o = BC[[reg]][["eff_summ"]]
  o$region = reg
  oo = rbind(oo, o)
}
oo$fishyr = jitter(oo$fishyr, amount=0.02)

color_map = c("#E69F00", "#56B4E9",  "#CC79A7" , "#D55E00", "#F0E442")[1:length(regions)]
shapes = c(15, 17, 19, 21, 23)[1:length(regions)]

pl = ggplot( oo, aes(x=fishyr, y=discard_rate, ymin=discard_rate-discard_rate_sd, ymax=discard_rate+discard_rate_sd, fill=region, colour=region) ) +
  geom_line( alpha=0.9, linewidth=1 ) +
  geom_point(aes(shape=region), size=5, alpha=0.7 )+
  geom_pointrange()  + # Vertical line with point in the middle
  geom_errorbar(width = 0.1, col="brown") + # Standard error bars
  geom_point(size = 1.5, col="darkred") +
  scale_colour_manual(values=color_map) +
  scale_fill_manual(values=color_map) +
  scale_shape_manual(values = shapes) +
  theme_light( base_size = 18) + 
  labs(x="Year", y="Discard rate of snow crab (Observed, by weight)" ) + 
  theme( legend.position="inside", legend.position.inside=c(0.2, 0.8), legend.title=element_blank() )
  
(pl)

```

&nbsp;  $~$  <br /> 



```{r}
#| label: fig-observed-softshell-map
#| eval: true 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| echo: false 
#| layout-ncol: 2
#| fig-cap: "Map of observed soft shell locations."
#| fig-subcap:  
#|   - "N-ENS"
#|   - "CFA23"
#|   - "CFA24" 

fns = c( 
  "nens_soft_crab_positions_68.png",
  "cfa23_soft_crab_positions_68.png",
  "cfa24_soft_crab_positions_68.png"
)
 
include_graphics( file.path( media_supplementary, fns )  ) 

```
 

&nbsp;  $~$  <br /> 



```{r}
#| label: fig-bcd-map
#| eval: true
#| echo: false 
#| output: true
#| fig-cap: "Bitter crab disease observations since 2008"
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
  
fns = file.path( media_loc, c(
  "BCD_map.png"  
) )

knitr::include_graphics( fns ) 

``` 


&nbsp;  $~$  <br /> 



```{r}
#| label: fig-survey-locations-map 
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| echo: false 
#| layout-ncol: 2
#| fig-cap: "Survey locations."
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

loc = file.path( data_loc, "output", "maps", "survey.locations" )
years = year_assessment + c(0:-3)
fn = check_file_exists( file.path( loc, paste( "survey.locations", years, "png", sep=".") ))
include_graphics( fn )
```
 

&nbsp;  $~$  <br /> 

  

```{r}
#| label: fig-growth-modes
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: "Modes identified from survey data using *Kernel Density Estmation* of local moving data windows and associated inferred growth derived from *Kernel Mixture Models*. Legend: sex|maturity|instar"
#| fig-subcap: 
#|   - "Female modes"
#|   - "Male modes"
#|   - "Female growth trajectory"
#|   - "Male growth trajectory"

fns = file.path( media_loc, c(
  "density_f_imodes.png", 
  "density_m_imodes.png", 
  "plot_growth_female.png",
  "plot_growth_male.png"
))

include_graphics( fns )
```
 
 
&nbsp;  $~$  <br /> 


```{r}
#| label: fig-sizefeq-male
#| eval: true
#| output: true
#| fig-cap: "Size-frequency (areal density; no/km$^2$) histograms by carapace width of male Snow Crab. The vertical line represents the legal size (95 mm). Immature animals are shown with light coloured bars, mature with dark."
#| fig-dpi: 144
#| fig-height: 10

if (params$sens==1) {
  sf_outdir = file.path( p$annual.results, "figures", "size.freq", "survey", "period1")
} else if (params$sens==2) {
  sf_outdir = file.path( p$annual.results, "figures", "size.freq", "survey_split", "period1")
}

include_graphics( file.path( sf_outdir,  "male.denl.png" ) )

```
 

 

&nbsp;  $~$  <br /> 


 
```{r}
#| label: fig-sizefeq-female
#| eval: true
#| output: true
#| fig-cap: "Size-frequency (areal density; no/km$^2$) histograms by carapace width of female Snow Crab. Immature animals are shown with light coloured bars, mature with dark."
#| fig-dpi: 144
#| fig-height: 10


if (params$sens==1) {
  sf_outdir = file.path( p$annual.results, "figures", "size.freq", "survey", "period1")
} else if (params$sens==2) {
  sf_outdir = file.path( p$annual.results, "figures", "size.freq", "survey_split", "period1")
}

fn = file.path( sf_outdir, "female.denl.png" )

include_graphics( fn )
```

 


&nbsp;  $~$  <br /> 




```{r}
#| label: fig-totno-female-mat-timeseries
#| eval: true
#| output: true
#| fig-cap: "The crude, unadjusted geometric mean of mature female density log$_{10}$(no/km$^2$) from the Snow Crab survey. Error bars represent 95\\% Confidence Intervals. Note the absence of data in 2020. Prior to 2004, surveys were conducted in the Spring."
#| fig-dpi: 144
#| fig-height: 4 
 
 
if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

fn = file.path( ts_outdir, paste("totno.female.mat", "png", sep=".") )
include_graphics( fn )

```

&nbsp;  $~$  <br /> 


```{r}
#| label: fig-totno-female-mat-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: "Mature female density log$_{10}$(no/km$^2$) from the Snow Crab survey."
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""
 
map_outdir = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual" )
map_years  = year_assessment + c(0:-3)
  
fn = check_file_exists( file.path( 
  map_outdir, "totno.female.mat", paste( "totno.female.mat", map_years, "png", sep="." )  
) )

include_graphics( fn )
```
 
&nbsp;  $~$  <br /> 

```{r}
#| label: fig-sexratio-mat-timeseries
#| eval: true
#| output: true
#| fig-cap: "The crude, unadjusted geometric mean of sex ratios (proportion female) of mature Snow Crab. Error bars represent 95\\% Confidence Intervals. Note the absence of data in 2020. Prior to 2004, surveys were conducted in the Spring."
#| fig-dpi: 144
#| fig-height: 4 
 
if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

fn = file.path( ts_outdir, paste("sexratio.mat", "png", sep=".") )
include_graphics( fn )

```

&nbsp;  $~$  <br /> 




```{r}
#| label: fig-sexratio-mat-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: "Snow Crab survey sex ratios (proportion female) of mature Snow Crab."
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""
 

map_outdir = file.path( p$project.outputdir, "maps", "survey", "snowcrab", "annual" )
map_years  = year_assessment + c(0:-3)
  
fn = check_file_exists( file.path( 
  map_outdir, "sexratio.mat", paste( "sexratio.mat", map_years, "png", sep="." )  
) )

include_graphics( fn )
```
 
&nbsp;  $~$  <br /> 
 
```{r}
#| label: fig-R0-timeseries
#| eval: true
#| output: true
#| fig-cap: "The crude, unadjusted geometric mean fishable biomass density log$_{10}$(t/km$^2$) from the Snow Crab survey. Error bars represent 95\\% Confidence Intervals. Note the absence of data in 2020. Prior to 2004, surveys were conducted in the Spring."
#| fig-dpi: 144
#| fig-height: 4 
 
if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

fn = file.path( ts_outdir, paste("R0.mass", "png", sep=".") )
include_graphics( fn )

```

&nbsp;  $~$  <br /> 



```{r}
#| label: fig-R0-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: "Snow Crab survey fishable component biomass density log$_{10}$(t/km$^2$)."
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""
 

map_outdir = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual" )
map_years  = year_assessment + c(0:-3)
  
fn = check_file_exists( file.path( 
  map_outdir, "R0.mass", paste( "R0.mass", map_years, "png", sep="." )  
) )

include_graphics( fn )
```


 

```{r}
#| label: fig-bottom-temperatures-timeseries
#| eval: true
#| echo: false 
#| output: true
#| fig-cap: "Bottom temperatures"
#| fig-subcap: 
#|   - "Annual variations in bottom temperature observed during the Snow Crab survey. The horizontal (black) line indicates the long-term, median temperature within each subarea. Error bars represent standard errors."
#|   - "Posterior densities of predicted average bottom temperatures from an analysis of historical temperature data using [carstm](https://github.com/jae0/carstm). Red horizontal line is at $7^\\circ$C."
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold 

tloc = file.path( data_loc, "assessments", year_assessment, "timeseries"  )

fns = c( 
  file.path("survey", "t.png"), 
  "temperature_bottom.png" 
)

knitr::include_graphics( file.path( tloc, fns) )

``` 

&nbsp;  $~$  <br /> 


```{r}
#| label: fig-figures-temperature-bottom-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: "Bottom temperature ($^\\circ$C) observed during the Snow Crab survey."
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""
 
map_outdir = file.path( p$project.outputdir, "maps", "survey", "snowcrab", "annual" )
map_years  = year_assessment + c(0:-3)
  
fn = check_file_exists( file.path( 
  map_outdir, "t", paste( "t", map_years, "png", sep="." )  
) )

include_graphics( fn )
```
 


&nbsp;  $~$  <br /> 





```{r}
#| label: fig-speciescomposition-biplot
#| eval: true
#| echo: false 
#| output: true
#| fig.show: hold 
#| fig-dpi: 144
#| fig-height: 6
#| fig-cap: "Species ordination (PCA: eigenanalysis of correlation matrices). PC1 is associatd with bottom temperatures. PC2 is associated with depth. Snow crab is shown as an orange dot."
     
xlab = paste("PC1 (", pca$variance_percent[1], "%)", sep="" )
ylab = paste("PC2 (", pca$variance_percent[2], "%)", sep="" )
plot( PC2 ~ PC1, pcadata, type="n", xlab=xlab, ylab=ylab )
text( PC2 ~ PC1, labels=vern, data=pcadata, cex=0.7, col="slateblue"  )
i = grep("Snow crab", pcadata$vern, ignore.case=TRUE)
points( PC2 ~ PC1, pcadata[i,], pch=19, cex=3.0, col="darkorange" )

```

 
```{r}
#| label: fig-speciescomposition-timeseries
#| eval: true
#| echo: false 
#| output: true
#| fig.show: hold 
#| fig-dpi: 144
#| fig-height: 6
#| fig-cap: "Species composition in time. Primary gradient (PC1) is related to bottom temperatures; second (PC2) to depth. Groundfish surveys were not conducted in 2020 and 2022 in the snow crab domain. Snow crab surveys were not conducted in 2020, and incomplete in 2022 in S-ENS."
#| fig-subcap: 
#|   - "Mean annual PC1 score."
#|   - "Mean annual PC2 score."
 
pc = c(1, 2)

spc_loc = file.path( data_root, "aegis", "speciescomposition", "modelled", "default" )
fnpr = file.path( spc_loc, "figures", paste("pca", pc, "_time.png", sep="" ) )
knitr::include_graphics( fnpr ) 

``` 
 
```{r}
#| label: fig-speciescomposition-map-pc1
#| eval: true
#| echo: false 
#| output: true
#| fig.show: hold 
#| fig-dpi: 144
#| fig-height: 6
#| fig-cap: "Species composition (PC1) in space. Groundfish surveys were not conducted in 2020 and 2022 in the snow crab domain. Snow crab surveys were not conducted in 2020, and incomplete in 2022 in S-ENS."
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""
   
yrsplot =  year_assessment + c(0:-3)
 
vn = "pca1"

spc_loc = file.path( data_root, "aegis", "speciescomposition", "modelled", "default", "maps" )
fns = file.path( spc_loc, paste( vn, "predictions", yrsplot, "png", sep=".") )

knitr::include_graphics( fns ) 

``` 
 
```{r}
#| label: fig-speciescomposition-map-pc2
#| eval: true
#| echo: false 
#| output: true
#| fig.show: hold 
#| fig-dpi: 144
#| fig-height: 6
#| fig-cap: "Species composition (PC2) in space. Groundfish surveys were not conducted in 2020 and 2022 in the snow crab domain. Snow crab surveys were not conducted in 2020, and incomplete in 2022 in S-ENS."
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""
   
yrsplot =  year_assessment + c(0:-3)
 
vn = "pca2"

spc_loc = file.path( data_root, "aegis", "speciescomposition", "modelled", "default", "maps" )
fns = file.path( spc_loc, paste( vn, "predictions", yrsplot, "png", sep=".") )

knitr::include_graphics( fns ) 

``` 
 

```{r}
#| label: fig-predator-biplot
#| eval: true
#| echo: false 
#| output: true
#| fig.show: hold 
#| fig-dpi: 144
#| fig-height: 6
#| fig-cap: "Main predators of snow crab on Scotian Shelf of Atlantic Canada (1999-2020). Relative location of snow crab predators (green) in the species composition ordination. Snow crab in orange. Of 58,287 finfish stomach samples, 159 had snow crab (0.28%). There is no information on snow crab diet in the database."
#| fig-subcap: 
#|   - "Total random effect (space and space-time) - PC1"
#|   - "Total random effect (space and space-time) - PC1"
#|   - "Mean annual PC1 score."
#|   - "Mean annual PC2 score."
#| layout-ncol: 1

# potential predators:
lookup= c( "cod", "halibut", "sculpin", "skate", "plaice", "hake", "wolffish", "atlantic cod", "atlantic halibut", "longhorn sculpin", "thorny skate", "striped atlantic wolffish", "haddock", "american plaice", "smooth skate", "winter skate", "white hake", "shorthorn sculpin", "eelpout newfoundland", "squirrel or red hake", "sea raven", "ocean pout", "barndoor skate", "Squid" )

xlab = paste("PC1 (", pca$variance_percent[1], "%)", sep="" )
ylab = paste("PC2 (", pca$variance_percent[2], "%)", sep="" )
plot( PC2 ~ PC1, pcadata, type="n", xlab=xlab, ylab=ylab )
text( PC2 ~ PC1, labels=vern, data=pcadata, cex=0.75, col="slategrey"  )
i = grep("Snow crab", pcadata$vern, ignore.case=TRUE)
points( PC2 ~ PC1, pcadata[i,], pch=19, cex=3.0, col="darkorange" )
j = NULL
for (k in lookup) j = c(j, grep( k, pcadata$vern, ignore.case=TRUE))    
j = unique(j)
points( PC2 ~ PC1, pcadata[j,], pch=19, cex=2.0, col="lightgreen" )
text( PC2 ~ PC1, labels=vern, data=pcadata[j,], cex=0.75, col="darkgreen"  )

```



```{r}
#| label: fig-diet-biplot
#| eval: true
#| echo: false 
#| output: true
#| fig.show: hold 
#| fig-dpi: 144
#| fig-height: 6
#| fig-cap: "Relative location of snow crab prey (green) in the species composition ordination. Snow crab in orange. Most of the potental prey are found to the right of snow crab (i.e. colder-water species) at a variety of depths."
#| fig-subcap: 
#|   - "Total random effect (space and space-time) - PC1"
#|   - "Total random effect (space and space-time) - PC1"
#|   - "Mean annual PC1 score."
#|   - "Mean annual PC2 score."
#| layout-ncol: 1

# potential food items:
lookup= c( "echinoderm", "polychaete", "maldane", "nereis", "shrimp", "pandalus", "rock crab", "toad crab", "lesser toad crab", "quahog", "artica islandica", "mollusc", "mytilus", "modiolus", "hiatella", "starfish", "sea anemone", "brittle star", "sea star", "sea anemone", "ophiura", "ophiopholis", "edwardsia", "metridium", "euphasid" )

xlab = paste("PC1 (", pca$variance_percent[1], "%)", sep="" )
ylab = paste("PC2 (", pca$variance_percent[2], "%)", sep="" )
plot( PC2 ~ PC1, pcadata, type="n", xlab=xlab, ylab=ylab )
text( PC2 ~ PC1, labels=vern, data=pcadata, cex=0.75, col="slategrey"  )
i = grep("Snow crab", pcadata$vern, ignore.case=TRUE)
points( PC2 ~ PC1, pcadata[i,], pch=19, cex=3.0, col="darkorange" )
j = NULL
for (k in lookup) j = c(j, grep( k, pcadata$vern, ignore.case=TRUE))    
j = unique(j)
points( PC2 ~ PC1, pcadata[j,], pch=19, cex=2.0, col="lightgreen" )
text( PC2 ~ PC1, labels=vern, data=pcadata[j,], cex=0.75, col="darkgreen"  )

```


 

```{r}
#| label: fig-competitor-biplot
#| eval: true
#| echo: false 
#| output: true
#| fig.show: hold 
#| fig-dpi: 144
#| fig-height: 6
#| fig-cap: "Potential competitors of snow crab on Scotian Shelf of Atlantic Canada (1999-2020). Relative location of snow crab predators (green) in the species composition ordination. Snow crab in orange."
#| fig-subcap: 
#|   - "Total random effect (space and space-time) - PC1"
#|   - "Total random effect (space and space-time) - PC1"
#|   - "Mean annual PC1 score."
#|   - "Mean annual PC2 score."
#| layout-ncol: 1

# potential predators:

lookup= c( "pandalus", "Jonah Crab", "Atlantic Rock Crab" , "Toad Crab", "Hyas Coarctatus", "Northern Stone Crab" )   # add more here: most are not direct compeitors as they have slightly different depth/temp preferences


xlab = paste("PC1 (", pca$variance_percent[1], "%)", sep="" )
ylab = paste("PC2 (", pca$variance_percent[2], "%)", sep="" )
plot( PC2 ~ PC1, pcadata, type="n", xlab=xlab, ylab=ylab )
text( PC2 ~ PC1, labels=vern, data=pcadata, cex=0.75, col="slategrey"  )
i = grep("Snow crab", pcadata$vern, ignore.case=TRUE)
points( PC2 ~ PC1, pcadata[i,], pch=19, cex=3.0, col="darkorange" )
j = NULL
for (k in lookup) j = c(j, grep( k, pcadata$vern, ignore.case=TRUE))    
j = unique(j)
points( PC2 ~ PC1, pcadata[j,], pch=19, cex=2.0, col="lightgreen" )
text( PC2 ~ PC1, labels=vern, data=pcadata[j,], cex=0.75, col="darkgreen"  )

```



```{r}
#| label: fig-fb-habitat-timeseries
#| eval: true
#| echo: false 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "Habitat viability (probability; fishable Snow Crab). Means and 95\\% Credible Intervals are presented."

loc = file.path( data_loc, "modelled", "default_fb", "aggregated_habitat_timeseries" )
include_graphics( file.path( loc, "habitat_M0.png") )

```


```{r}
#| label: fig-fb-habitat-map
#| eval: true
#| echo: false 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "Habitat viability (probability; fishable Snow Crab)."
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""
#|   - ""
#|   - ""
#|   - ""
#|   - ""
#|   - ""
#|   - ""
#| layout: [[100], [100], [50,50], [50,50], [50,50], [50,50]]
  
loc = file.path( data_loc, "modelled", "default_fb", "predicted_habitat" )
vn = "habitat."
yrsplot =  year_assessment + c(0:-9)

fns = file.path( loc, paste( vn, yrsplot, ".png", sep="") )
include_graphics( fns )
 
```

 

&nbsp;  $~$  <br /> 


```{r}
#| label: fig-fbindex-timeseries
#| eval: true
#| echo: false 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "The fishable biomass index (t) predicted by CARSTM of Snow Crab survey densities. Error bars represent Bayesian 95\\% Credible Intervals. Note large errors in 2020 when there was no survey."

fn = file.path( data_loc, "modelled", "default_fb", "aggregated_biomass_timeseries" , "biomass_M0.png")
include_graphics( fn )

```

&nbsp;  $~$  <br /> 


```{r}
#| label: fig-fbindex-map
#| eval: true
#| echo: false 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "Biomass index log~10(t/km$^2$) predicted from the Snow Crab survey."
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""
#|   - ""
#|   - ""
#|   - ""
#|   - ""
#|   - ""
#|   - ""
#| layout: [[100], [100], [50,50], [50,50], [50,50], [50,50]]
  
loc = file.path( data_loc, "modelled", "default_fb", "predicted_biomass_densities" )
yrsplot =  year_assessment + c(0:-9)

fns = file.path( loc, paste( "biomass", yrsplot, "png", sep=".") )
include_graphics( fns )
```
 
&nbsp;  $~$  <br /> 


```{r}
#| label: fig-logisticPredictions
#| eval: true
#| echo: false
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "Fishable, posterior mean modelled biomass (pre-fishery; kt) are shown in dark orange. Light orange are posterior samples of modelled biomass (pre-fishery; kt) to illustrate the variability of the predictions. The biomass index (post-fishery, except prior to 2004) after model adjustment by the model catchability coefficient is in gray."
#| fig-subcap:
#|   - "N-ENS"
#|   - "S-ENS"
#|   - "4X"

loc = file.path( data_loc, "fishery_model", year_assessment, "logistic_discrete_historical" )
fns = file.path( loc, c(
  "plot_predictions_cfanorth.png",
  "plot_predictions_cfasouth.png",
  "plot_predictions_cfa4x.png"
) )

include_graphics( fns )

```
 
```{r}
#| label: fig-logisticPredictions2
#| eval: true
#| echo: false
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "Fishable, posterior mean modelled biomass (post-fishery; kt) are shown in dark orange. Light orange are posterior samples of modelled biomass (post-fishery; kt) to illustrate the variability of the predictions. The biomass index (post-fishery, except prior to 2004) after model adjustment by the model catchability coefficient is in gray."
#| fig-subcap:
#|   - "N-ENS"
#|   - "S-ENS"
#|   - "4X"

loc = file.path( data_loc, "fishery_model", year_assessment, "logistic_discrete_historical" )
fns = file.path( loc, c(
  "plot_predictions_postfishery_cfanorth.png",
  "plot_predictions_postfishery_cfasouth.png",
  "plot_predictions_postfishery_cfa4x.png"
) )

include_graphics( fns )

```
 
 
```{r}
#| label: fig-logisticFishingMortality
#| eval: true
#| echo: false
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "Time-series of modelled instantaneous fishing mortality from Model 1, for N-ENS (left), S-ENS (middle), and 4X (right). Samples of the posterior densities are presented, with the darkest line being the mean."
#| fig-subcap:
#|   - "N-ENS"
#|   - "S-ENS"
#|   - "4X"

odir = file.path( fishery_model_results, year_assessment, "logistic_discrete_historical" )
fns = file.path( odir, c(
  "plot_fishing_mortality_cfanorth.png",
  "plot_fishing_mortality_cfasouth.png",
  "plot_fishing_mortality_cfa4x.png"
))

include_graphics( fns )

```



```{r}
#| label: fig-ReferencePoints
#| eval: true
#| echo: false
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "Harvest control rules for the Scotian Shelf Snow Crab fisheries."

fn = file.path( media_loc, "harvest_control_rules.png")
include_graphics( fn )

```
 


```{r}
#| label: fig-logistic-hcr
#| eval: true
#| echo: false
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "Reference Points (fishing mortality and modelled biomass) from the Fishery Model, for N-ENS (left), S-ENS (middle), and 4X (right). The large yellow dot indicates most recent year and the 95\\% CI. Not: the model does not account for illegal and unreported landings, and interspecific interactions. Prefishery."
#| fig-subcap:
#|   - "N-ENS"
#|   - "S-ENS"
#|   - "4X"

odir = file.path( fishery_model_results, year_assessment, "logistic_discrete_historical" )

fns = file.path( odir, c(
  "plot_hcr_cfanorth.png" ,
  "plot_hcr_cfasouth.png",
  "plot_hcr_cfa4x.png"
) )

include_graphics( fns )

```
 

```{r}
#| label: fig-logistic-hcr2
#| eval: true
#| echo: false
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "Reference Points (fishing mortality and modelled biomass) from the Fishery Model, for N-ENS (left), S-ENS (middle), and 4X (right). The large yellow dot indicates most recent year and the 95\\% CI. Not: the model does not account for illegal and unreported landings, and interspecific interactions. Post-fishery."
#| fig-subcap:
#|   - "N-ENS"
#|   - "S-ENS"
#|   - "4X"

odir = file.path( fishery_model_results, year_assessment, "logistic_discrete_historical" )

fns = file.path( odir, c(
  "plot_hcr_postfishery_cfanorth.png" ,
  "plot_hcr_postfishery_cfasouth.png",
  "plot_hcr_postfishery_cfa4x.png"
) )

include_graphics( fns )

```
 



# Appendices

## Rapid climate change

```{r}
#| label: fig-ocean-currents
#| eval: true
#| echo: false 
#| output: true
#| fig-cap: "Ocean currents in the Maritimes.  Source: [DFO](https://www.dfo-mpo.gc.ca/oceans/publications/soto-rceo/2018/atlantic-ecosystems-ecosystemes-atlantiques/index-eng.html)"
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold 

fn = file.path( media_loc, c(
  "maritimes_currents.png" 
) )
knitr::include_graphics( fn ) 

```

&nbsp;  $~$  <br /> 



```{r}
#| label: fig-rapid-climate-change
#| eval: true
#| echo: false 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold 
#| fig-cap: "Global surface (2 meter) air temperature.  Source: [The Crisis Report](https://richardcrim.substack.com/p/the-crisis-report-99) and [James E. Hansen](https://www.columbia.edu/~jeh1/mailings/2024/ICJ.PressBriefing.09December2024.pdf). Note 2023 was and El Nino year."
#| fig-subcap: 
#|   - "Anomalies relative to pre-industrial baseline."
#|   - "Seasonal variations by year."
#|   - "Annual temperatures."

fn = file.path( media_loc, c(
  "gst_anomaly.png",
  "gst_seasonal.png",
  "gst_ts.png"
) )

knitr::include_graphics( fn ) 

```

&nbsp;  $~$  <br /> 

  
```{r}
#| label: fig-ocean-productivity-chla
#| eval: true
#| echo: false 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold 
#| fig-cap: "Chlorphyll-a in the NW Atlantic.  Source: [Copernicus Marine Service](https://marine.copernicus.eu/access-data/ocean-monitoring-indicators/chlorophyll-and-primary-production)"
#| fig-subcap: 
#|   - "Full timeseries of surface Chl-a estimated from satellite imagery."
#|   - "Trends in surface Chl-a over time."

fn = file.path( media_loc, c(
  "copernicus_chla.png",
  "copernicus_chla_map.png"
) )
knitr::include_graphics( fn ) 

```

&nbsp;  $~$  <br /> 


## Fishery model diagnostics


```{r}
#| label: fig-logistic-surplus-production
#| results: asis
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig-cap: "Surplus (Shaeffer) production. Each year is represented by a different colour."
#| fig-subcap:
#|   - "N-ENS"
#|   - "S-ENS"
#|   - "4X"

fns = file.path( fm_loc, paste("plot_surplus_production_", regions, ".png", sep="" ) ) 

include_graphics( fns )
   
``` 




```{r}
#| label: fig-logistic-prior-posterior-K
#| results: asis
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig-cap: "Prior-posterior comparisons of carrying capacity (K; kt)."
#| fig-subcap:
#|   - "N-ENS"
#|   - "S-ENS"
#|   - "4X"
   
fns = file.path( fm_loc, paste("plot_prior_K_", regions, ".png", sep="" ) ) 

include_graphics( fns )
   

``` 



```{r}
#| label: fig-logistic-prior-posterior-r
#| results: asis
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig-cap: "Prior-posterior comparisons of th iIntrinsic rate of biomass increase (r)."
#| fig-subcap:
#|   - "N-ENS"
#|   - "S-ENS"
#|   - "4X"

  
fns = file.path( fm_loc, paste("plot_prior_r_", regions, ".png", sep="" ) ) 

include_graphics( fns )

``` 

 
```{r}
#| label: fig-logistic-prior-posterior-q
#| results: asis
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig-cap: "Prior-posterior comparisons of the catchability coefficient (q)."
#| fig-subcap:
#|   - "N-ENS"
#|   - "S-ENS"
#|   - "4X"

  
fns = file.path( fm_loc, paste("plot_prior_q1_", regions, ".png", sep="" ) ) 

include_graphics( fns )

``` 



```{r}
#| label: fig-logistic-prior-posterior-obserror
#| results: asis
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig-cap: "Prior-posterior comparisons of Observation error (kt)."
#| fig-subcap:
#|   - "N-ENS"
#|   - "S-ENS"
#|   - "4X"

fns = file.path( fm_loc, paste("plot_prior_bosd_", regions, ".png", sep="" ) ) 

include_graphics( fns )

``` 
 


```{r}
#| label: fig-logistic-prior-posterior-processerror
#| results: asis
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig-cap: "Prior-posterior comparisons of Model process error (kt)."
#| fig-subcap:
#|   - "N-ENS"
#|   - "S-ENS"
#|   - "4X"
 
fns = file.path( fm_loc, paste("plot_prior_bpsd_", regions, ".png", sep="" ) ) 

include_graphics( fns )

``` 
 
  
```{r}
#| label: fig-logistic-state-space
#| results: asis
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig-cap: "State space (kt): year vs year+1."
#| fig-subcap:
#|   - "N-ENS"
#|   - "S-ENS"
#|   - "4X"

fns = file.path( fm_loc, paste("plot_state_space_", regions, ".png", sep="" ) ) 

include_graphics( fns )

``` 
 


 

&nbsp;  $~$  <br /> 

## Movement

 
```{r}
#| label: fig-movement-tracks
#| eval: true
#| echo: false 
#| output: true
#| fig-cap: "Snow Crab movement"
#| fig-subcap: 
#|   - "Tracks from 1996-2004"
#|   - "Tracks from 2004 - present"
#|   - "Distance between mark and recapture (km)"
#|   - "Minimum speed for each mark-recapture event (km/month)"
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| layout: [[100], [100], [50,50] ]
  
fns = file.path( media_loc, c(
  "movement0.png", 
  "movement.png" ,
  "snowcrab_movement_distances.png", 
  "snowcrab_movement_rates.png" 
) )

knitr::include_graphics( fns ) 

``` 

## Habitat

&nbsp;  $~$  <br /> 

```{r}
#| label: fig-temp-depth
#| eval: true
#| echo: false 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "Habitat preferences associated with depth and temperature."
 
fn2=file.path( media_loc, "viable_habitat_depth_temp.png" )
knitr::include_graphics( c( fn2 ) ) 

```
 

```{r}
#| label: fig-viable-habitat-persistent
#| eval: true
#| echo: false 
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| fig-cap: "Persistent habitat, independent of temperature, time, etc."

fn1 = file.path( media_loc, "viable_habitat.png" ) 
knitr::include_graphics( c(fn1  ) ) 

```

## Co-occurring species

```{r}
#| label: fig-atlcod-timeseries
#| eval: true
#| output: true
#| fig-cap: "Mean density of Atlantic cod log$_{10}$(no/km$^2$) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

species_predator = 10

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```

&nbsp;  $~$  <br /> 



```{r}
#| label: fig-atlcod-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: Atlantic cod, density; log$_{10}$(no/km$^2$). 
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""
 
map_years  = year_assessment + c(0:-3)
  
species_predator = 10
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```

```{r}
#| label: fig-haddock-timeseries
#| eval: true
#| output: true
#| fig-cap: "Mean density of Haddock log$_{10}$(no/km$^2$) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}
species_predator = 11

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```

&nbsp;  $~$  <br /> 



```{r}
#| label: fig-haddock-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: Haddock, density; log$_{10}$(no/km$^2$). 
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""
 
map_years  = year_assessment + c(0:-3)
  
species_predator = 11
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```


&nbsp;  $~$  <br /> 


  
```{r}
#| label: fig-halibut-timeseries
#| eval: true
#| output: true
#| fig-cap: "Mean density of Halibut log$_{10}$(no/km$^2$) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

species_predator = 30

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```

&nbsp;  $~$  <br /> 


```{r}
#| label: fig-halibut-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: Halibut, density; log$_{10}$(no/km$^2$). 
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""
 
map_years  = year_assessment + c(0:-3)
  
species_predator = 30
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```


  
```{r}
#| label: fig-amerplaice-timeseries
#| eval: true
#| output: true
#| fig-cap: "Mean density of American plaice log$_{10}$(no/km$^2$) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

species_predator = 40

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```

&nbsp;  $~$  <br /> 



```{r}
#| label: fig-amerplaice-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: American plaice, density; log$_{10}$(no/km$^2$). 
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

map_years  = year_assessment + c(0:-3)
  
species_predator = 40
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```
 
  
```{r}
#| label: fig-stripatlwolffish-timeseries
#| eval: true
#| output: true
#| fig-cap: "Mean density of Striped Atlantic wolffish log$_{10}$(no/km$^2$) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

species_predator = 50

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```

&nbsp;  $~$  <br /> 



```{r}
#| label: fig-stripatlwolffish-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: Striped Atlantic wolffish, density; log$_{10}$(no/km$^2$). 
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

map_years  = year_assessment + c(0:-3)
  
species_predator = 50
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```
 

&nbsp;  $~$  <br /> 


  
```{r}
#| label: fig-thornyskate-timeseries
#| eval: true
#| output: true
#| fig-cap: "Mean density of Thorny skate log$_{10}$(no/km$^2$) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

species_predator = 201

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```

&nbsp;  $~$  <br /> 



```{r}
#| label: fig-thornyskate-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: Thorny skate, density; log$_{10}$(no/km$^2$). 
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

map_years  = year_assessment + c(0:-3)
  
species_predator = 201
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```
 


&nbsp;  $~$  <br /> 





```{r}
#| label: fig-northernshrimp-timeseries
#| eval: true
#| output: true
#| fig-cap: "Mean density of Northern shrimp log$_{10}$(no/km$^2$) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

species_predator = 2211

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```

&nbsp;  $~$  <br /> 



```{r}
#| label: fig-northernshrimp-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: Northern shrimp, density; log$_{10}$(no/km$^2$). 
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

map_years  = year_assessment + c(0:-3)
  
species_predator = 2211
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```
 

&nbsp;  $~$  <br /> 


```{r}
#| label: fig-jonahcrab-timeseries
#| eval: true
#| output: true
#| fig-cap: "Mean density of Jonah crab log$_{10}$(no/km$^2$) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

species_predator = 2511

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```

&nbsp;  $~$  <br /> 



```{r}
#| label: fig-jonahcrab-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: Jonah crab, density; log$_{10}$(no/km$^2$). 
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

map_years  = year_assessment + c(0:-3)
  
species_predator = 2511
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```
 
&nbsp;  $~$  <br /> 



```{r}
#| label: fig-lyrecrab-timeseries
#| eval: true
#| output: true
#| fig-cap: "Mean density of Arctic Lyre crab (Lesser toad crab) log$_{10}$(no/km$^2$) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

species_predator = 2521

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```

&nbsp;  $~$  <br /> 



```{r}
#| label: fig-lyrecrab-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: Arctic Lyre crab, density; log$_{10}$(no/km$^2$). 
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

map_years  = year_assessment + c(0:-3)
  
species_predator = 2521
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```
 
&nbsp;  $~$  <br /> 


```{r}
#| label: fig-nstonecrab-timeseries
#| eval: true
#| output: true
#| fig-cap: "Mean density of Northern stone crab log$_{10}$(no/km$^2$) from surveys with 95\\% Confidence Intervals."
#| fig-dpi: 144
#| fig-height: 4 

if (params$sens==1) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey")
} else if (params$sens==2) {
  ts_outdir = file.path( p$annual.results, "timeseries", "survey", "split")
}

species_predator = 2523

bc_vars = paste("ms.no", species_predator, sep='.')
fn = file.path( ts_outdir, paste(bc_vars, "png", sep=".") )
include_graphics( fn )

```

&nbsp;  $~$  <br /> 



```{r}
#| label: fig-nstonecrab-map
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: Northern stone crab, density; log$_{10}$(no/km$^2$). 
#| fig-subcap: 
#|   - ""
#|   - ""
#|   - ""
#|   - ""

map_years  = year_assessment + c(0:-3)
  
species_predator = 2523
bc_vars = paste("ms.no", species_predator, sep='.')
outdir_bc = file.path( p$project.outputdir, "maps", "survey", "snowcrab","annual", "bycatch" )

fn = check_file_exists( file.path( outdir_bc, bc_vars, paste(bc_vars, map_years, "png", sep=".") ) )
include_graphics( fn )
    
```
 
