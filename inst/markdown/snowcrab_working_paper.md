---
title: "Snow Crab status, Scotian Shelf 2023"
subtitle: "Working Paper"
author: "Jae S. Choi"
date: "`r format(Sys.time(), '%Y-%m-%d')`"
always_allow_html: true
output:
  # bookdown::word_document2:
  bookdown::pdf_document2:
    fig_caption: yes
    toc: true
    toc_depth: 2
    latex_engine: pdflatex
    number_sections: true
    keep_tex: true
#   # bookdown::pdf_document2:
#   #   toc: true
#   #   fig_caption: yes
header-includes: 
  - \usepackage{graphicx}
  - \usepackage{float}
  # - \usepackage{subfig}
  #bookdown::html_document2:
  #bookdown::word_document2:
  #bookdown::pdf_document2:
  # word_document:  
  #  reference_docx: media/RES2021-eng.docx
  # csl: media/csas.csl
  # csl is citation style download from: https://www.zotero.org/styles
  # csas.csl copied from "csasdown"
# bibliography: media/snowcrab.bib
fontsize: 12pt
params:
  year.assessment: 2023
  media_loc: "media"
  debugging: FALSE 
---


<!-- Preamble

This is a Markdown document ... To create HTML or PDF, etc, run: 


  make quarto FN=snowcrab_working_paper YR=2023 SOURCE=~/projects/bio.snowcrab/inst/markdown WK=~/bio.data/bio.snowcrab/assessments  DOCEXTENSION=html  # {via Quarto}

  make rmarkdown FN=snowcrab_working_paper YR=2023 SOURCE=~/projects/bio.snowcrab/inst/markdown WK=~/bio.data/bio.snowcrab/assessments DOCTYPE=bookdown::pdf_document2 DOCEXTENSION=pdf # {via Rmarkdown}

  make pdf FN=snowcrab_working_paper  # {via pandoc}

Alter year and directories to reflect setup or copy Makefile and alter defaults to your needs.
  
 
 
-->



```{r setup, include=FALSE}
  require(knitr)
  knitr::opts_chunk$set(
    root.dir = data_root,
    echo = FALSE,
    out.width="6.2in",
#     dev.args = list(type = "cairo"),
    fig.retina = 2,
    dpi=192
  )

  # inits and data loading (front load all required data)

  require(aegis)
  
  year.assessment = params$year.assessment
  year_previous = year.assessment - 1
  p = bio.snowcrab::load.environment( year.assessment=year.assessment )
  SCD = project.datadirectory("bio.snowcrab")
  media_loc = params$media_loc
  
  # fishery_model_results = file.path( "/home", "jae", "projects", "dynamical_model", "snowcrab", "outputs" )
  fishery_model_results = file.path( SCD, "fishery_model" )

  sn_env = snowcrab_load_key_results_to_memory( year.assessment, debugging=params$debugging, return_as_list=TRUE  ) 

  attach(sn_env)

  # predator diet data
  diet_data_dir = file.path( SCD, "data", "diets" )
  require(data.table) # for speed
  require(lubridate)
  require(stringr) 
  require(gt)  # table formatting
  library(janitor)
  require(ggplot2)
  require(aegis) # map-related 
  require(bio.taxonomy)  # handle species codes

  # assimilate the CSV data tables:
  # diet = get_feeding_data( diet_data_dir, redo=TRUE )  # if there is a data update
  diet = get_feeding_data( diet_data_dir, redo=FALSE )
  tx = taxa_to_code("snow crab")  
  # matching codes are 
  #  spec    tsn                  tx                   vern tx_index
  #1  528 172379        BENTHODESMUS           BENTHODESMUS     1659
  #2 2522  98427        CHIONOECETES SPIDER QUEEN SNOW UNID      728
  #3 2526  98428 CHIONOECETES OPILIO        SNOW CRAB QUEEN      729
  # 2 and 3 are correct

  snowcrab_predators = diet[ preyspeccd %in% c(2522, 2526), ]  # n=159 oservations out of a total of 58287 observations in db (=0.28% of all data)
  snowcrab_predators$Species = code_to_taxa(snowcrab_predators$spec)$vern
  snowcrab_predators$Predator = factor(snowcrab_predators$Species)
  
  counts = snowcrab_predators[ , .(Frequency=.N), by=.(Species)]
  setorderv(counts, "Frequency", order=-1)
  
  # species composition
  psp = speciescomposition_parameters( yrs=p$yrs, carstm_model_label="default" )
  pca = speciescomposition_db( DS="pca", p=psp )  

  pcadata = as.data.frame( pca$loadings )
  pcadata$vern = stringr::str_to_title( taxonomy.recode( from="spec", to="taxa", tolookup=rownames( pcadata ) )$vern )

  # bycatch summaries
  o_cfaall = observer.db( DS="bycatch_summary", p=p,  yrs=p$yrs, region="cfaall" )
  o_cfanorth = observer.db( DS="bycatch_summary", p=p,  yrs=p$yrs, region="cfanorth" )   
  o_cfasouth = observer.db( DS="bycatch_summary", p=p,  yrs=p$yrs, region="cfasouth" )   
  o_cfa4x = observer.db( DS="bycatch_summary", p=p,  yrs=p$yrs, region="cfa4x" )   

```
 

 
 
# ABSTRACT

In the Scotian Shelf Ecosystem (SSE), Snow Crab (*Chionoecetes opilio*) has been a dominant macro-invertebrate since the decline of the groundfish in the 1990s. They are mostly observed in deep, soft-bottom substrates ranging from 60 to 300 m and at temperatures generally less than 6${}^\circ$ C. The SSE Snow Crab are in the southern-most extreme of their spatial distribution in the Northwest Atlantic Ocean and vulnerable to climate variability. 

The assessment of 4VWX Snow Crab status (TOR #1; see Terms of reference, below) is based on fishery independent surveys with a focus upon indicators of abundance, reproductive potential, recruitment, and exploitation rates (TOR #2). Robust, Bayesian, Hurdle-type spatiotemporal models are used to incorporate habitat viability attributable to ecosystem influences (depth, species composition and bottom temperature variations, etc.) and so account for biases due to deficiencies of sample design and incomplete surveys. Models of fishery dynamics are used to infer historical abundance and a new model for forward projections is used to help evaluate the consequences of different harvest levels (TOR #3). Further, we highlight some of the bycatch in the survey that highlights potentially elevated predation-related mortality in some areas (TOR #4).  

Commercial catch rates and other fishery statistics are reported. Overall, fishing effort was much reduced in 2022, especially in S-ENS and 4X. TACs were mostly caught, with the exception of 4X. Catch rates were at or near all-time highs in all areas, except for 4X. At-Sea-Observed information was sparse and so fishery bycatch and snow crab size distributions were not evaluated.

Survey indices suggest a downtrend in numerical densities of Snow Crab recruiting to fishable sizes, especially in N-ENS. Similarly, declines in numerical densities of mature females, especially in N-ENS, suggest potential declines in egg production and long term recruitment. Mature snow crab sex ratios remain close to balanced levels, except again in N-ENS. Crude geometric mean biomass densities of the fishable component have also declined in all areas, and for S-ENS and 4X, it was at an historical low. Adjustment by spacetime index modelling suggest that the biomass index is also low, but habitat variability, especially as they relate to bottom temperature variability was very significant. 

Modelled inference (Model 1) closely follows the biomass index and suggests that the overall biomass of fishable Snow Crab has marginally declined over the last year and that fishing mortality has declined and are well below FMSY, with the exception of N-ENS. Model 1 places N- and S-ENS in the "healthy" zone and 4X in the "cautious" zone. All regions have flexibility. However, due to a lower level of recruitment in N-ENS and potentially elevated predation mortality, care should be taken in 2023 until recruitment returns.


# SUMMARY


- Fishing effort in `r year.assessment` were `r e_nens` $\times 10^3$ trap hauls in N-ENS, `r e_sens` $\times 10^3$ trap hauls in S-ENS and  `r e_4x` $\times 10^3$ trap hauls in 4X. This represents a  change of `r dt_e_nens` %, `r dt_e_sens` % and  `r dt_e_4x` %, respectively, relative to the previous year. 

- Landings in `r year.assessment` were `r l_nens` t in N-ENS; `r l_sens` t in S-ENS; and `r l_4x` t in CFA 4X (season ongoing), representing a change of `r dt_l_nens`%, `r dt_l_sens`% and  `r dt_l_4x`%, respectively, relative to `r year_previous`. Total Allowable Catches (TACs) for `r year.assessment` were `r tac_nens` t, `r tac_sens` t and `r tac_4x` t in N-ENS, S-ENS and 4X, respectively. 

- Non-standardized fishery catch rates in `r year.assessment` were `r c_nens`, `r c_sens`  and `r c_4x` kg/trap haul in N-ENS, S-ENS and 4X, respectively. This represents a change of respectively, `r dt_c_nens` %, `r dt_c_sens` % and  `r dt_c_4x` % (season ongoing) relative to the previous year. Though the spatial extent of exploitation was smaller, many of the exploited area show elevated catch rates.

- Commercial catches of soft-shelled (newly moulted) Snow Crab for `r year.assessment` were `r cc_soft_nens`%, `r cc_soft_sens`% and `r cc_soft_4x`%, respectively, in N-ENS S-ENS and 4X (season ongoing). In `r year_previous`, it was `r cc_soft_nens_p`%, `r cc_soft_sens_p`% and `r cc_soft_4x_p`%, respectively. Sensecent (CC5) crab levels were higher in N- and S-ENS than in previous years, though low and inconsistent sampling effort makes this uncertain.

- Bycatch of non-target species is low (<< 1% of total catch) in all Snow Crab fishing areas; however, the last reliable estimate comes from 2017. 

- Little to no recruitment is expected for the next 1-3 years in N-ENS. In S-ENS, continued moderate levels of recruitment are expected. In 4X, low to moderate levels of recruitment are expected for 2 years. Egg and larval production is expected to be high in the next year in all areas except N-ENS. 

- In N-ENS, the modelled biomass (pre-fishery) of Snow Crab in `r year.assessment` was `r round(B_north[t0], 2)` t, relative to `r round(B_north[t1], 2)` t in `r year_previous`. In S-ENS, the modelled biomass (pre-fishery) was `r round(B_south[t0], 2)` t, relative to `r round(B_south[t1], 2)` t in `r year_previous`. In 4X, the `r year.assessment`-`r year.assessment+1` season's modelled biomass (pre-fishery) was `r round(B_4x[t0], 2)` t, relative to `r round(B_4x[t1], 2)` t in the `r year_previous`-`r year.assessment` season. 

- Average bottom temperatures observed in the 2022 Snow Crab survey were near or above historical highs in all areas. Average viable habitat surface area has declined to historical lows in 2022.

- Based on stomach sampling, Atlantic Halibut, Atlantic Wolffish, Thorny Skate, and other skate species appear to be the predominant predators of Snow Crab on the Scotian Shelf. Overall, higher predation mortality seems likely in N- and S-ENS and lower in 4X. Co-occurring species such as shrimp densities have declined, possibly due to large-scaled environmental change.

- In N-ENS, though recruitment continues at low levels, a gap in future recruitment to the fishery is expected for the next 1-3 years in N-ENS. N-ENS still exists in the healthy zone and so still has flexibility in approach. However, fishing mortality may have been higher than desirable. A more conservative harvest strategy would enable bridging this coming recruitment gap. A reduced TAC is prudent. 

- In S-ENS, recruitment to the fishery is likely to continue at a moderate rate for the upcoming season. The S-ENS stock remains in the healthy zone. Exploitation rates derived from the fishery model have been declining in recent years. S-ENS has considerable flexibility in approach. All TAC directions are viable. 

- In 4X, low to moderate levels of recruitment are expected for 2 years. 4X exists in the "cautious zone". The area is also in the southern-most extent of Snow Crab distribution in the North Atlantic and viable habitat has been depressed for many years. A reduced TAC is prudent. 




# BACKGROUND

## Terms of Reference (TOR)

To provide an assessment of the:

  1. status of Snow Crabs
  2. relative abundance and exploitation rates
  3. evaluate consequences of different harvest levels upon abundance and exploitation rates
  4. bycatch in survey


Amendments of the Fisheries Act in 2022 (*Fish Stock Provisions*), have encoded management approaches and requirements relative to biological reference points. In the SSE Snow Crab Fishery, these biological reference points were originally defined heuristically with a  *phenomenological* model (a descriptive, model without mechanism, herein "Model 1"; Choi 2023), with all participants to the management process being aware that it was simplistic at best and that it served as supplemental information to help provide context rather than be the tool used for management. As there is now a strong legal requirement that such reference points become a primary management tool, there is a need to define these reference point with greater care and evaluate their robustness and utility with other models. 

## Ecology

Snow Crab (*Chionoecetes opilio*; Figure \@ref(fig:crab-image)) are a circumpolar, subarctic species. In the Northwest Atlantic Ocean, they are found from northern Labrador to near the Gulf of Maine. The Scotian Shelf Ecosystem (SSE) represents the southern-most part of this distribution and so most influenced by environmental variability. In the SSE, habitat preference is generally for soft mud bottoms, at depths from 60 to 300 m and temperatures from -1 to 6${}^\circ$ C. Temperatures greater than 7${}^\circ$ C are known to be metabolically detrimental. Mature females and immature crab are found in slightly more complex habitats with shelter at marginally shallower depths (Choi et al. 2022). Their primary food items show ontogenetic shift. In their bottom dwelling phase, they rely upon: shrimp, fish, starfish, sea urchins, worms, detritus, large zooplankton, other crabs, molluscs, sea snails, and sea anemones. Their predators are primarily humans, Atlantic Halibut, Thorny Skate, Wolfish, Atlantic Cod, seals, American Plaice, squids, and other crabs. Crab in the size range of 3 to 30 mm carapace width (CW) are particularly vulnerable to predation, as are soft-shelled crab in the spring molting season. Cannibalism of immature crab by females have also been observed. More detailed information with regards to snow habitat requirements and spatiotemporal distributions of different life stages in the SSE are presented in Choi et al. (2022).


## Precautionary Approach 

A precautionary approach is a guiding principle for fisheries harvesting in Canada ("Fish Stock Provisions" Amendments of the Fisheries Act, 2022). As such it is worthwhile emphasizing the many existing measures and fishing practices in the SSE Snow Crab fishery that are inherently precautionary that are not always appreciated by the general public:

- **The spawning stock (breeding females)** is  *completely protected*.

- **Conservative exploitation strategies** have generally been the norm since the mid-2000s.

- **Spatial refugia** are associated with Marine Protected Areas (MPAs), along the continental slopes, and much of the western inshore portion of CFA 24 and 4X.

- **Temporal refugia** associated with fishing seasons. 

- **Biological refugia** -- most life stages are protected. 

  - No removal of female crab, the "spawning stock". None.

  - Immature and soft-shelled (newly-molted) Snow Crab are easily damaged. They are not harvested, and handling mortality is minimized via spring harvesting, voluntary area closures, and at-sea observer monitoring of soft-shell incidence, helping to maximize the potential yield per animal to the biomass. 

  - Reproductive potential of spawning stock is not disrupted. Most removals of males occur after mating and sub-legal mature crab (able to reproduce) are never removed. Their relative numbers (sex ratios) are monitored.

- **Collaborative management** 

  - Fishers prioritize scientific, fisheries-independent, evidence-based decision making. They participate in temperature monitoring and tagging studies. They also fund directly the majority of the science work required to collect fisheries-independent data through Collaborative Agreements. They fully understand the intrinsic value of the long-term sustainability of Snow Crab and its connection to their own long term socioeconomic success. 
  
  - Fishers have had many years of observational experience and are well informed about the biology and ecology of Snow Crab and associated species. As such, they represent a distributed knowledge network that embodies a great depth of traditional and historical knowledge of the SSE. Their participation in scientific reviews continues to keep them up-to-date with emerging issues and able to be active agents in an **adaptive system** of fishery management.   

  - Management of the fishery was initially based on effort controls (season, license, trap limits) from 1982 to 1993 with harvest occurring from June-November, of hard-shelled males larger than 95 mm Carapace Width (CW). Additional management measures were introduced from 1994 to 1999: Individual Boat Quotas (IBQs), Total Allowable Catches (TACs), Individual Transferable Quotas (ITQs), 100% dockside monitoring, mandatory logbooks and at-sea monitoring by certified observers and satellite-based Vessel Monitoring Systems (VMS); and adoption of biodegradeable mesh to reduce ghost-fishing of lost gear.

  - Spring fishing efforts were voluntarily implemented in N-ENS and S-ENS to avoid soft-shelled  crab handling mortality, in the mid-2000s. The fishery in area 4X having a fall to winter season, spanning calendar years, have managed to avoid unnecessarily soft-shell captures from its inception in the late 1990s. Voluntary measures to land only fully mature crab, easily identified by large-claws were also agreed upon in the late-2000s to increase the longevity of the sexually mature population, and reduce the number of snow crab individuals removed and were based upon conservation and precautionary principles. 



# FISHERY STATISTICS

## Fishing effort

Fishing effort in `r year.assessment` was `r e_nens` $\times 10^3$, `r e_sens` $\times 10^3$ and `r e_4x` $\times 10^3$ trap hauls in N-ENS, S-ENS and 4X, respectively. Relative to the previous year, these represent changes of `r dt_e_nens` %, `r dt_e_sens` % and  `r dt_e_4x` %, respectively (Tables \@ref(tab:table-fishery-nens), \@ref(tab:table-fishery-sens), \@ref(tab:table-fishery-4x), Figure \@ref(fig:effort-timeseries)). Fishing effort was consistent between `r year.assessment` and `r year_previous` in terms of  spatial distribution. In S-ENS, there was, however, a minor spatial contraction to inshore areas and away from the area 23-24 boundary (Figure \@ref(fig:effort-map)). This is presumed to be related to, in part, warm water incursions into the area (see Ecosystem considerations, below).


## Fishery landings and TACs

Landings across time are shown in Figure \@ref(fig:landings-timeseries). In `r year.assessment`, they were `r l_nens`, `r l_sens` and `r l_4x` t, in N-ENS, S-ENS and 4X (season ongoing), respectively. Relative to `r year_previous`, they represent changes of `r dt_l_nens`%, `r dt_l_sens`% and  `r dt_l_4x`%, respectively (Tables \@ref(tab:table-fishery-nens), \@ref(tab:table-fishery-sens), \@ref(tab:table-fishery-4x)). Total Allowable Catches (TACs) for `r year.assessment` were `r tac_nens` t, `r tac_sens` t and `r tac_4x` t in N-ENS, S-ENS and 4X, respectively. 


Additional carry forward allowance was implemented by DFO Fisheries Management of up to 25% of the 2020 quota to the 2021 season for N-ENS and S-ENS Snow Crab fishery. These were a response to COVID-19 related uncertainties with safe fishing activity. The carry-forward amounts were: 11.2 t in N-ENS, and 217.4 t in S-ENS (Tables \@ref(tab:table-fishery-nens), \@ref(tab:table-fishery-sens)). In 2022, landings in all areas were below respective TACs.
 

The landings in N-ENS for 2022 and 2021 were similar in their spatial patterns (Figure \@ref(fig:sse-map)). In S-ENS, landings, as with fishing effort, shifted slightly inshore and away from the area 23-24 boundary (Figure \@ref(fig:landings-map)). There were no landings on the continental slope areas of S-ENS in 2022 which continues to serve as a "reserve" for Snow Crab from fishing. The landings in 4X for 2022 as with 2021, were primarily in the area just south of Sambro, bordering onto area 24. In N-ENS, most landings occured in the spring. 
  

## Fishery catch rates

Non-standardized fishery catch rates in `r year.assessment` were `r c_nens`, `r c_sens`  and `r c_4x` kg/trap haul in N-ENS, S-ENS and 4X, respectively. This represents a change of respectively, `r dt_c_nens` %, `r dt_c_sens` % and  `r dt_c_4x` % (season ongoing) relative to the previous year (Tables \@ref(tab:table-fishery-nens), \@ref(tab:table-fishery-sens), \@ref(tab:table-fishery-4x), Figures \@ref(fig:cpue-timeseries)). Though the spatial extent of exploitation was smaller, many of the exploited area show elevated catch rates (Figure \@ref(fig:cpue-map)).


## At-Sea-Observed information

Carapace condition of the fished component is determined from At-Sea-Observed catches. In 2021, both N-ENS and 4X were not sampled by At-Sea-Observers. In 2022, 4X was not sampled by At-Sea-Observers. Estimates of carapace condition since 2020 are unreliable as they represent only small areas of the fishing grounds and short time periods relative to the whole fishing season (Figure \@ref(fig:observer-locations-map)). 

In the exploited fraction of Snow Crab, Carapace Condition (CC) is an index of the approximate time since the last molt so describes the relative development and subsequent decay of the carapace. CC1 signifies a newly molted crab, soft-shelled, with no epibiont (e.g., barnacles) growth. CC2 crab have begun to harden, but is still considered to be soft and of no commercial value. CC3 and CC4 represent ideal commercial crab. The oldest carapace condition (CC5) signifies extensive shell decay with no expectation of survival into the next year.

Commercial catches of soft-shelled crab were `r cc_soft_nens`% (low sampling), `r cc_soft_sens`% (low sampling) and `r cc_soft_4x`% (no sampling; season ongoing) in N-ENS, S-ENS and 4X, respectively for `r year.assessment`. In `r year_previous`, it was `r cc_soft_nens_p`% (no sampling), `r cc_soft_sens_p`% (low sampling) and `r cc_soft_4x_p`% (no sampling), respectively. Generally, higher soft-shell indicates incoming recruitment to the fishery and their  handling potential and unnecessary handling/discard mortality.

In 2022, CC5 crab levels were higher in N- and S-ENS than in previous years, though again low and inconsistent sampling effort (space, time aliasing) makes this uncertain.

Bycatch in the Snow Crab fishery is also monitored from the At-Sea-Observed catches. Due to a paupacy of consistent data, very little can be said of bycatch after 2017. Historically, bycatch in the Snow Crab fishery has been minimal (Tables 4, 5) with increasing levels as a function of increasing water temperature: there are higher by catch in warmer conditions, primarily other Crustacea (crab and lobster), especially when viable habitat is low and animals with divergent habitat requirements find themselves next to each other. Low bycatch has been attributed to trap design (top entry conical traps), the large mesh size (5.25 inches, knot to knot) and the passive nature of the gear (Hebert et al., 2001).  

 

# SURVEY INDICES
 
Survey catch rates are confounded by numerous factors that vary across space and time. This is because distributions of Snow Crab and variables that influence these distributions vary across space and time, while survey effort does not (they are mostly fixed stations, in time and space). Survey catch rates depend upon seasonality, bottom temperatures, predator distributions, food availability, reproductive behavior, substrate/shelter availability, relative occurrence of soft and immature crab, species co-occurrence, vessel/captain experience, gear configuration, ambient currents, etc. These factors are taken into account where possible with *CARSTM* (Conditional AutoRegressive SpatioTemporal Models; Choi 2020, Choi et al. 2022; and references therein), a robust extension of Bayesian, Hurdle-type generalized linear random effects models to spatially and temporally connected units.  

The survey was not conducted in 2020 due to Covid-19 related to health and safety uncertainties. In 2022, the survey did not complete due to mechanical issues with the survey vessel. Inshore areas of S-ENS were most affected (Figure \@ref(fig:survey-locations-map)). 


## Size and carapace condtion

Geometric mean size of the mature male component of Snow Crab has varied between 83 to 108 mm CW (Figure \@ref(fig:meansize-male-mat)) in the historical record. Changes in size can be caused by  poor environmental conditions that encourage early maturation, size-selective predation, loss of the largest individuals due to fishing and the start or end of a recruitment pulse. N-ENS and 4X have seen the greatest volatility. S-ENS has been stable. In 2022 it has declined in S-ENS (though this may be an result of incomplete surveys) and 4X (death of the largest crab possibly due to warm water incursions) while increasing in N-ENS (low recruitment).  

The carapace condition of mature male crab captured by the Snow Crab survey are shown in Figure (\@ref(fig:sizefeq-male-survey-cc)). The relative distributions are mostly comparable across time. However, there has been a slight increase in the proportion of CC5 crab in 2022 for S-ENS. An increase was also seen in At-sea-observed fishery data and so could be indicative of accelerated ageing associated with high temperature conditions (see below). However, the incomplete sampling in the inshore areas of S-ENS and in the At-Sea-Observed data makes inference uncertain. 
 

## Recruitment

Quantitative determination of recruitment levels into the fishable component is confounded by a number of factors. These include terminal molt (the timing offset of molting in spring and the survey in the fall), the inability to age crab, the inability to predict the age that male crab will terminally molt and the focus of the survey upon the fishable component which results in a biased under-representation of recruitment. As habitat requirements/preferences of adolescent and larval components are divergent to those of the fishable component (Choi et al 2022), their estimation is challenging. 

Based on size-frequency histograms of the male Snow Crab population, little to no recruitment is expected for the next 1-3 years in N-ENS (Figure \@ref(fig:sizefeq-male)). In S-ENS, continued moderate levels of recruitment are expected. In 4X, low to moderate levels of recruitment are expected for 2 years. 

## Reproduction

In all areas, there was substantial and continued recruitment of female crab into the mature (egg-bearing) segment of the population from 2016-2022 (Figure \@ref(fig:sizefeq-female)). However, in N-ENS for 2022, a decline in numerical densities was observed as well as low densities of adolescent females. Egg and larval production is expected to be moderate to high in the next year in all areas except N-ENS. Mature female abundance (Figure \@ref(fig:fmat-timeseries)) distributions are heterogeneous and often in shallower areas.  

 
## Sex ratios

The sex ratios (proportion female) of the mature component is particularly important as locally unbalanced ratios can impact upon encounter rates and ultimately reproductive success (Figure \@ref(fig:sexratio-mature)). In the ESS, there is generally a lack of females, in contrast to, for example, the Gulf of St-Lawrence where the reverse is often the case. Higher sex ratios are usually found in inshore and bottom slope areas (Figure \@ref(fig:sexratio-map)). A decline in sex ratios has been observed since 2017 in N-ENS (Figure \@ref(fig:sexratio-mature)). In S-ENS the sex ratio increased from 20% in 2021 to just under 35% in 2022. In 4X, mature sex ratios are more balanced and currently near the 50% level.  



## Biomass Density

The fishable component is defined as Snow Crab that are male, mature, and larger than 95 mm CW. The crude **biomass density** (that is, unadjusted, geometric mean fishable per unit swept area by the trawl) of the fishable component, sampled by the Snow Crab survey is shown in Figures \@ref(fig:fbGMTS), and \@ref(fig:fbgeomean-map). A peak in crude biomass density was observed in 2009 to 2014 and has since been declining in all areas. Note that high and low biomass density areas fluctuate with time (Figure \@ref(fig:fbgeomean-map)). Biomass density, however, does not equate to total biomass as the areas occupied by crab can contract, expand and shift with environmental conditions and ecosystem variability. 


## Biomass Index

The fishable **biomass index** is a statistically modelled estimate after adjustment for ecosystem indices and spatiotemporal autocorrelation; Figure \@ref(fig:fbindex-map)). More specifically, it was computed using conditional auto-regressive spatio-temporal models (Choi 2020) of probability of observation, numerical abundance of positive valued occurrence and mean size with environmental covariates of depth, substrate, bottom temperature, two axes of species composition ordinations. Further, the biomass index model infers (that is, *imputes*) the spatiotemporal distribution of the biomass density of the fishable component (Figure \@ref(fig:fbindex-map)). The locations with no sampling from 2022 are also imputed. Upon aggregation of the solutions, we see that the overall biomass has had several cycles (Figure \@ref(fig:fbindex-timeseries)). Note, however, the more elevated uncertainty for the imputed 2020 and partially imputed 2022 estimates (Figure \@ref(fig:fbindex-timeseries)).

The magnitudes of the biomass index are optimistically high as the spatial expansion uses areal units with large surface areas, larger than the patchiness of Snow Crab distributions (spatial autocorrelation length is <20 km, on average; Choi 2020). As such, it should only be seen as a spatially and temporally comparable relative index of abundance. As the size of areal units decreases to scales smaller than the spatial autocorrelation length, the solutions will converge to reality; but the costs in time and resources of achieving this are prohibitive as sampling intensity required to achieve such a result would require at least another two orders of magnitude increase in sampling effort (in time and space).

The spatial distribution of the biomass index has been consistent over the past six years, with a peak in overall biomass index in 2019 and 2020 (Figure \@ref(fig:fbindex-map)). Since then, a  reduction was observed throughout the region, with the exception of the core areas. A contraction of spatial range in 4X and the western parts of S-ENS were also evident in 2021 to 2022. Upon aggregation, the biomass index declined marginally in all areas (Figure \@ref(fig:fbindex-timeseries)).


## Modelled Biomass

The biomass index along with fishery removals are used to fit a *Logistic Biomass Dynamics Model* (Model 1; Choi 2023) to determine fishable **modelled biomass** (biomass estimated from the fisheries model; Figure \@ref(fig:logisticPredictions)) and relevant biological reference points (i.e., carrying capacity and fishing mortality at maximum sustainable yield, or F~MSY~). In N-ENS, the modelled biomass (pre-fishery) of Snow Crab in `r year.assessment` was `r round(B_north[t0], 2)` t, relative to `r round(B_north[t1], 2)` t in `r year_previous`. In S-ENS, the `r year.assessment` modelled biomass (pre-fishery) was `r round(B_south[t0], 2)` t, relative to `r round(B_south[t1], 2)` t in `r year_previous`. In 4X, the modelled biomass (pre-fishery) for the `r year.assessment`-`r year.assessment+1` season was `r round(B_4x[t0], 2)` t, relative to `r round(B_4x[t1], 2)` t for the `r year_previous`-`r year.assessment` season. 
  

## Fishing Mortality


In N-ENS, the `r year.assessment` fishing mortality is estimated to have been `r round(FM_north[t0],3)` (annual exploitation rate of `r round(100*(exp(FM_north[t0])-1),2)`%), up from the  `r year_previous` rate of `r round(FM_north[t1],3)` (annual exploitation rate of `r round(100*(exp(FM_north[t1])-1),1)`%; Figure \@ref(fig:logisticFishingMortality)).

In S-ENS, the `r year.assessment` fishing mortality is estimated to have been `r round(FM_south[t0],3)` (annual exploitation rate of `r round(100*(exp(FM_south[t0])-1),1)`%), decreasing marginally from the `r year_previous` rate of `r round(FM_south[t1],3)` (annual exploitation rate of `r round(100*(exp(FM_south[t1])-1),1)`%; Figure \@ref(fig:logisticFishingMortality)). Localized exploitation rates are likely higher, as not all areas for which biomass is estimated are fished (e.g., continental slope areas and western, inshore areas of CFA 24).

In 4X, the `r year.assessment`-`r year.assessment+1` season (ongoing), fishing mortality is estimated to have been `r round(FM_4x[t0],3)` (annual exploitation rate of `r round(100*(exp(FM_4x[t0])-1),1)`%), decreasing from the `r year.assessment-1`-`r year.assessment` season rate of `r round(FM_4x[t1],3)` (annual exploitation rate of `r round(100*(exp(FM_4x[t1])-1),1)`%; Figure \@ref(fig:logisticFishingMortality)). Localized exploitation rates are likely higher, as not all areas for which biomass is estimated are fished. 


## Reference Points
 
Reference points are used to guide harvest strategies (Canada Gazette 2022; DFO 2013; Figures \@ref(fig:logistic-hcr)). More specifically, the state of the fishery relative to the following "Reference Points" are assessed. In terms of modelled biomass, Lower and Upper Stock Reference are 25% and 50% of carrying capacity which delineate "critical", "cautious" and "healthy" zones. In terms of exploitation rates, the Upper Removal Reference is the exploitation rate that we try not to go beyond; it is defined in term of the fishing mortality associated with Maximum Sustainable Yield (FMSY). In the biomass dynamics model, FMSY = $r /2$. As $r \approx 1$ for snow crab, FMSY $\approx 0.5$ is expected. 

The operational target exploitation changes depending upon the "zone" in which a population lands. When in the "healthy" zone, the rule of thumb has been to keep annual exploitation rates between 10% to 32% of the available biomass ($F = 0.11, 0.36$, respectively). In the "cautious" zone, the rule of thumb has been to keep annual exploitation rates between 0% to 20% ($F = 0, 0.22$, respectively). In the "critical" zone, fishery closure is considered until recovery is observed, where recovery indicates at a minimum, modelled biomass > LSR. Other biological and ecosystem considerations such as recruitment, spawning stock (female) biomass, size structure, sex ratios and environmental and ecosystem conditions, provide additional guidance within each range.

Model 1 estimates key Reference Points as shown in Table 6 and Figure \@ref(fig:ReferencePoints). The related PA thresholds can be computed as:

  - Lower Stock Reference (LSR): $K/4$
  - Upper Stock Reference (USR): $K/2$
  - Upper Removal Reference (URR): keep fishing mortality below *FMSY* = $r/2$

Model 1 suggests the current state of the fishable components to be (Figure \@ref(fig:logistic-hcr)):  

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


# Ecosystem Considerations

## Bottom Temperature

Average bottom temperatures **observed** in the 2022 Snow Crab survey were near or above historical highs in all areas (Figure \@ref(fig:bottom-temperatures-survey)). A general warming trend has been observed in the Snow Crab survey since the early 1990s on the Scotian Shelf (Choi et al. 2022). Temperatures are more stable in N-ENS than S-ENS; 4X exhibits the most erratic and highest annual mean bottom temperatures. Of particular note, the observed temperatures in the 2022 Snow Crab survey for S-ENS increased well above the average. This is expected to be a source of bias and certainty for S-ENS abundance predictions as it is outside the range normally encountered by the statistical models. Furthermore, the Groundfish surveys did not operate over Snow Crab grounds in 2020 and 2022 and so temperature data is also very sparse for the area of interest.

Upon aggregation and modelling of historical temperature data, the average temperature is found to have increased well beyond the $7^\circ$C threshold in 4X. N-ENS and S-ENS also continued to experience historical highs in bottom temperature and elevated spatial variability of bottom temperatures (Figure \@ref(fig:bottom-temperatures)). In particular, the spike observed in S-ENS (above) was less extreme and so likely due to a short term influx of warm water. 

Overall, since 1999, there has been a persistent spatial gradient of almost $15^\circ$C in bottom temperatures in the Maritimes Region (Figure \@ref(fig:bottom-temperatures-spatialeffect)). This large gradient in a spatially complex and temporally dynamic area makes assessment a particular challenge for a stenothermic organism such as Snow Crab (Choi et al. 2022). 


## Viable habitat

Snow Crab being cold water stenotherms, stability of environmental conditions is critical for their survival. The Maritimes Region being at the confluence of many oceanic currents renders the area highly variable. Rapid climate change and uncertainty exacerbates this situation. The viable habitat  estimated for each area across time has shown some variations (Figures \@ref(fig:fb-habitat-timeseries), \@ref(fig:fb-habitat-map)) in the historical record. As can be seen, 4X showed a significantly lower average viable habitat levels relative to the N-ENS and S-ENS.  A peak in average probability of observing fishable snow crab ("viable habitat") was observed in 2010 for 4X, 2011 for N-ENS and 2012 for S-ENS. Since 2015, the average viable habitat has declined to historical lows and remained so in 2022. 
 

## Viable habitat in a stage-dependent model

The Model 1 representations of Snow Crab populations are imperfect at best. It is a naive phenomenological model that fits a pattern rather than biological processes (see Choi 2023). Amongst its' difficulties are:

  - Projections into the future are always optimistic as rapid growth to carrying capacity is explicitly assumed in the model structure. For example, each area is expected to increase in biomass in the absence of fishing, based upon the 1-year projection of the logistic model  (Figure \@ref(fig:logisticPredictions)). Indeed, the more that a population is depressed closer to zero, the faster this recovery is expected to be. This of course is because it ignores mechanism: the presence or absence of recruitment in any given year. It ignores the presence or absence of reproductive females for long term production. It ignores the need of a species to maintain a dominance hierarchy and that it will not be replaced by a competitor. It ignores the possibility of multiple equilibria where other ecological balances are possible.

  - Discounting/ignoring other biological processes will bias the magnitude of rate processes estimated by the model. Thus, external disturbances due to extreme weather events, dynamics of predator and prey species, disease, undocumented removals, etc., being ignored by the model will have an effect upon estimates of Reference Points. They will generally result in higher model estimates of the intrinsic rate of increase that tries to account for these additional mortality factors. FMSY and potentially carrying capacity and so biomass-based reference points will be higher as a consequence, and so result in more optimistic assessments of status than might be reasonable.

  - Currently, Model 1's status assessment of all areas seem overly optimistic. This is especially the case in 4X where we have supporting information of very high natural mortality rates associated with predation and extreme environmental variability, through the loss of strong year classes and weak recruitment into the fishery. The same issues also exist for N-ENS, where strong adolescent crab year classes seem to disappear at a rate faster that might be expected, possibly due to elevated predation mortality in the last decade. Similarly, S-ENS also has had  very strong year classes diminishing rapidly before entry into the fishable component.    

## Predators, preys, competitors

Being long-lived, the influence of predators can be significant. Especially important are predators of the smaller immature and female Snow Crab. Increasing predation not only lowers the abundance and recruitment, it can also reduce the reproductive potential of Snow Crab and therefore long-term population cycles. N-ENS and S-ENS are well known to have skewed sex ratios with few mature females for extended periods of time, quite possibly linked to differential predation mortality (mature females being much smaller and full of fat-rich gonads and eggs). 

Based on stomach sampling, Atlantic Halibut, Atlantic Wolffish, Thorny Skate, and other skate species are known predators of Snow Crab. Large Atlantic Halibut with mature female Snow Crab in their stomachs have been reported. Anecdotal information of some seals having fed upon Snow Crab are also known. 

Some of these predators (e.g., Halibut; DFO 2018) have significantly increased in abundance in the Region. However, for all, the abundance and encounter rates in areas overlapping with snow crab habitat is more important, but this is not known. We do know from the bycatch in the Snow Crab survey that there are elevated areal **densities** with many snow crab trawl samples. This means that encounter rates will also likely increase and so too potentially predation mortality. However, high density does not equate to high abundance nor high total predation mortality; but this remains unknown and requires further analysis. The following presents information of areal density of co-occuring species; they are potential predators, competitors and prey. 

Atlantic Halibut densities have increased rapidly since 2010 (Figure \@ref(fig:halibut-timeseries); DFO 2018) on Snow Crab grounds. Most of these increases were towards The Gully, Slope Edge and near Sable Island (\@ref(fig:halibut-map)).

Thorny skate densities have been increasing as well (Figure \@ref(fig:thornyskate-timeseries)), especially in N-ENS and along the margins of Banquereau Bank (Figure \@ref(fig:thornyskate-map)). A minor decline in the densities have been seen in 4X.

Striped Atlantic Wolffish densities have been high, though declining in N-ENS since 2007 (Figure \@ref(fig:Wolffish-timeseries)). Highest densities were towards the Laurentian Channel (Figure \@ref(fig:Wolffish-map)).

Northern shrimp co-occur as they share similar habitat preferences and are also potential prey items of Snow Crab. Their numerical densities have declined after a peak in 2011, especially in S-ENS (Figures \@ref(fig:Shrimp-timeseries), \@ref(fig:Shrimp-map)). 

Lesser toad crab is a co-occurring species and potential competitor. Their numbers have declined to low levels throughout, after a peak in densities in 2007 and 2016 in N-ENS (Figures \@ref(fig:lessertoadcrab-timeseries), \@ref(fig:lessertoadcrab-map)).

Overall, higher predation mortality seems likely in N- and S-ENS and lower in 4X. Further shrimp with similar habitat preferences have declined, possibly due to large-scaled habitat variations and predation.


## Other sources of uncertainty
 
* Bycatch of other species in the Snow Crab fishery cannot be reliably computed as they are derived from At-Sea-Observed data. 

* Bycatch of Snow Crab in other fisheries remains an area requiring attention. Anecdotal information from fishers suggest illegal retention of sublegal and female Snow Crab as bycatch and their use as bait. Illegal removals are also known to occur; however, the scale of such activities are not known.

* Oil and gas exploration (especially the use of seismics), development and exploitation (pollution) continues to be issues as their long-term effects remain unknown, though they are known to be generally deleterious to most life.

* Undersea, high-voltage cables is another source of concern. Many marine organisms are known to be sensitive to electromagnetic fields: fish in particular through their lateral line canals use it to detect prey. As such, predators and prey distributions can be affected causing indirect effects upon Snow Crab. Direct effects are also possible through trenching, but long term effects remain unknown.  

* Marine Protected Areas (MPAs) continue to be developed (e.g., Canada Gazette 2016). The presence of a refuge from fishing activities is potentially positive for Snow Crab. However, positive effects upon other organisms (predators or prey) can have counter-balancing effects. The overall long-term effects of the MPAs upon Snow Crab are unknown.

* Capture of soft-shell Snow Crab is always a concern. Prompt and careful return of immature (small-clawed, non-terminally molted) crab to the water is an important conservation measure that will enhance the 2-3 year productivity of the fishable component.

* Illegal, unreported, and unregulated fishing activities are known to occur. Such activities hinder the application of a precautionary approach to the management of this resource and cause potential bias and uncertainty in Reference Point estimation.

All of the above uncertainties are human induced and unlike the larger scaled climatic and ecosystemic uncertainties, there is some measure of human intervention possible. To remain adaptive in the face of these and other as yet unknown uncertainties including the climatic and ecosystemic uncertainties of which we are already aware, is the true challenge. It requires a balance of both resilience and robustness (see, Choi and Patten 2001). The Precautionary Approach as practiced by the snow crab fishers in Maritimes Region represents a unique model of such an adaptive approach. It values qualitity information and the communication and discussion of approaches in a distributed, collective information network that goes well beyond the simplistic and potentially maladaptive rubric of a purely "reference-points" based approach. Continued vigilence is necessary.


# CONCLUSIONS  

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






\newpage

# Tables




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


*Table 4. Bycatch (kg) estimates from the N-ENS and S-ENS Snow Crab fishery. Estimates are extrapolated from At-Sea-Observed bycatch and biomass of catch (bycatch = [observed biomass of bycatch species / observed landings of Snow Crab] X total landings of Snow Crab). Reliable species specific data beyond 2017 are currently unavailable.*

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




*Table 5. Bycatch (kg) estimates from the 4X Snow Crab fishery. Estimates are extrapolated from At-Sea-Observed bycatch and biomass of catch (bycatch = [observed biomass of bycatch species / observed landings of Snow Crab] X total landings of Snow Crab). Reliable species specific data beyond 2017 are currently unavailable.*

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



<!-- &nbsp; -->


*Table 6. Reference points from the logistic biomass dynamics fishery model (Model 1): K is Carrying capacity (kt); and r is Intrinsic rate of increase (non-dimensional). Note that FMSY (fishing mortality associated with 'Maximum Sustainable Yield') is r/2. Similarly, BMSY (biomass associated with 'Maximum Sustainable Yield') is K/2. SD is posterior Standard deviations.*


|       |  $K$ [SD] | $r$ [SD] |
| :---: |    :----:   | :---: |
| N-ENS | `r K_north` [`r K_north_sd`] | `r r_north` [`r r_north_sd`] |
| S-ENS | `r K_south` [`r K_south_sd`] | `r r_south` [`r r_south_sd`] |
| 4X    | `r K_4x`   [`r K_4x_sd`] | `r r_4x` [`r r_4x_sd`] |



\clearpage

\newpage

# Figures

```{r crab-image, out.width='45%', echo=FALSE, fig.align = 'center', fig.cap = 'A mature male Snow Crab (*Chionoecetes opilio*).' }
include_graphics( file.path( params$media_loc, "snowcrab_image.png" ) )
# \@ref(fig:crab-image)
```





```{r sse-map, out.width='60%', echo=FALSE, fig.align = 'center', fig.cap = 'Map of the SSE and Crab Fishing Areas (CFAs). S-ENS is further subdivided for management into areas 23 to the northeast and 24 to the southwest.' }
include_graphics( file.path (params$media_loc, "area_map.png" ) )
# \@ref(fig:sse-map)
```

 

```{r effort-timeseries, echo=FALSE, out.width='75%', fig.align='center', fig.cap = 'Temporal variations in fishing effort $\\times 10^3$ trap hauls.' }
fn1=file.path( SCD, "assessments", year.assessment, "timeseries", "fishery",   "effort.ts.pdf" )
knitr::include_graphics( fn1 ) 
# \@ref(fig:effort-timeseries)  
```





```{r effort-map, echo=FALSE, out.width='48%', fig.show='hold',  fig.align='center', fig.cap = 'Snow Crab fishing effort from fisheries logbook data for previous and current years. Units are No. $\\times 10^3$ per (10 km X 10 km) grid.' }
loc0= file.path( SCD, "output", "maps", "logbook", "snowcrab", "annual", "effort" )
fn1 = file.path( loc0, paste( "effort", year_previous,   "png", sep=".") ) 
fn2 = file.path( loc0, paste( "effort", year.assessment, "png", sep=".") ) 
include_graphics(  c(fn1, fn2) )
#  \@ref(fig:landings-map) 
```
 


```{r landings-timeseries, echo=FALSE, out.width='75%',  fig.align='center', fig.cap = 'Landings (t) of Snow Crab on the SSE. For 4X, the year refers to the starting year of the season.  Inset is a closeup view of the timeseries for N-ENS and 4X.'}
include_graphics( file.path( SCD, "assessments", year.assessment, "timeseries", "fishery",   "landings.ts.pdf" ) )
# \@ref(fig:landings-timeseries)
``` 


```{r landings-map, echo=FALSE, out.width='48%', fig.show='hold', fig.align='center', fig.cap = 'Snow Crab landings from fisheries logbook data for previous and current years. Units are tons per 10 km x 10 km grid.' }
loc0= file.path( SCD, "output", "maps", "logbook", "snowcrab", "annual", "landings" )
fn1 = file.path( loc0, paste( "landings", year_previous,   "png", sep=".") ) 
fn2 = file.path( loc0, paste( "landings", year.assessment, "png", sep=".") ) 
knitr::include_graphics( c(fn1, fn2 ) )
#  \@ref(fig:landings-map)  
```


 

```{r cpue-timeseries, echo=FALSE, out.width='75%', fig.align='center', fig.cap = 'Temporal variations in crude catch rates of Snow Crab, expressed as kg per trap haul.'}
include_graphics( file.path( SCD, "assessments", year.assessment, "timeseries", "fishery",   "cpue.ts.pdf" ) ) 
# \@ref(fig:cpue-timeseries)  
```

 

```{r cpue-map, echo=FALSE, out.width='48%', fig.show='hold', fig.align='center', fig.cap = 'Snow Crab crude catch rates on the Scotian Shelf for previous and current years. Units are kg/trap haul per 10 km x 10 km grid.' }
loc0= file.path( SCD, "output", "maps", "logbook", "snowcrab", "annual", "cpue" )
fn1 = file.path( loc0, paste( "cpue", year_previous,   "png", sep=".") ) 
fn2 = file.path( loc0, paste( "cpue", year.assessment, "png", sep=".") ) 
knitr::include_graphics( c(fn1, fn2 ) )
# \@ref(fig:cpue-map)  
```


 
```{r observer-locations-map, out.width='45%', fig.show='hold', fig.align='center', fig.cap= 'Snow Crab At-sea-observer locations.' }
loc = file.path( SCD, "output", "maps", "observer.locations" )
yrsplot = year.assessment + c(0:-4)
fn4 = file.path( loc, paste( "observer.locations", yrsplot[4], "png", sep=".") )
fn3 = file.path( loc, paste( "observer.locations", yrsplot[3], "png", sep=".") )
fn2 = file.path( loc, paste( "observer.locations", yrsplot[2], "png", sep=".") )
fn1 = file.path( loc, paste( "observer.locations", yrsplot[1], "png", sep=".") )
include_graphics( c( fn4, fn3, fn2, fn1) )
# \@ref(fig:observer-locations-map)  
#*Figure XXX. Snow Crab observer locations.*
```

 
```{r observer-CC, echo=FALSE, out.width='45%', fig.show='hold', fig.align='center', fig.cap = 'Size frequency distribution of Snow Crab sampled by At-sea-observers, broken down by Carapace Condition (CC). For 4X, the year refers to the starting year of the season; the current season is ongoing. Vertical lines indicate 95 mm Carapace Width, the minimum legal commercial size.' }
  loc = file.path( SCD, "assessments", year.assessment, "figures", "size.freq", "observer")
  fn1 = file.path( loc, paste( "size.freqcfanorth", (year_previous), ".pdf", sep="" ) )
  fn2 = file.path( loc, paste( "size.freqcfanorth", (year.assessment  ), ".pdf", sep="" ) )
  fn3 = file.path( loc, paste( "size.freqcfasouth", (year_previous), ".pdf", sep="" ) )
  fn4 = file.path( loc, paste( "size.freqcfasouth", (year.assessment  ), ".pdf", sep="" ) )
  fn5 = file.path( loc, paste( "size.freqcfa4x", (year_previous), ".pdf", sep="" ) )
  fn6 = file.path( loc, paste( "size.freqcfa4x", (year.assessment  ), ".pdf", sep="" ) )
  include_graphics(  c(fn1, fn2, fn3, fn4, fn5, fn6) )
# \@ref(fig:observer-CC)  
# *Figure XXX. Size frequency distribution of Snow Crab sampled by At-Sea-Observers, broken down by Carapace Condition (CC). Left side are for `r year_previous` and right side for `r year.assessment`. Top row is N-ENS, middle row is S-ENS, and bottom row is 4X. For 4X, the year refers to the starting year of the season; the current season is ongoing. Vertical lines indicate 95 mm Carapace Width, the minimum legal commercial size.*
```

  

\clearpage


```{r survey-locations-map, out.width='45%', fig.show='hold', fig.align='center', fig.cap= 'Snow Crab survey locations.' }
loc = file.path( SCD, "output", "maps", "survey.locations" )
yrsplot = setdiff( year.assessment + c(0:-9), 2020)
fn6 = file.path( loc, paste( "survey.locations", yrsplot[6], "png", sep=".") )
fn5 = file.path( loc, paste( "survey.locations", yrsplot[5], "png", sep=".") )
fn4 = file.path( loc, paste( "survey.locations", yrsplot[4], "png", sep=".") )
fn3 = file.path( loc, paste( "survey.locations", yrsplot[3], "png", sep=".") )
fn2 = file.path( loc, paste( "survey.locations", yrsplot[2], "png", sep=".") )
fn1 = file.path( loc, paste( "survey.locations", yrsplot[1], "png", sep=".") )
include_graphics( c(fn6, fn5, fn4, fn3, fn2, fn1) )
# \@ref(fig:survey-locations-map)  
# *Figure XXX. Snow Crab survey locations.*
```


```{r meansize-male-mat, out.width='60%', echo=FALSE, fig.align='center', fig.cap = 'Mean size of mature male Snow Crab log10(CW; mm) from surveys with 95\\% Confidence Intervals'}
include_graphics(  file.path( SCD, "assessments", year.assessment, "timeseries", "survey", "cw.male.mat.mean.pdf" )  )
# \@ref(fig:meansize-male-mat)
```
  

```{r sizefeq-male-survey-cc, out.width='32%', fig.show='hold', echo=FALSE, fig.align='center', fig.cap = 'Size-frequency of mature male Snow Crab by carapace width (mm) and carapace condition from surveys. Rows are years and columns as N-ENS (left), S-ENS(middle and 4X(right).'}
  odir = file.path( SCD, "assessments", year.assessment, "figures", "size.freq", "carapacecondition" )
  fn1 = file.path( odir, "sizefreq.cfanorth.2019.pdf" ) 
  fn2 = file.path( odir, "sizefreq.cfasouth.2019.pdf" ) 
  fn3 = file.path( odir, "sizefreq.cfa4x.2019.pdf" ) 
  fn4 = file.path( odir, "sizefreq.cfanorth.2021.pdf" ) 
  fn5 = file.path( odir, "sizefreq.cfasouth.2021.pdf" ) 
  fn6 = file.path( odir, "sizefreq.cfa4x.2021.pdf" ) 
  fn7 = file.path( odir, "sizefreq.cfanorth.2022.pdf" ) 
  fn8 = file.path( odir, "sizefreq.cfasouth.2022.pdf" ) 
  fn9 = file.path( odir, "sizefreq.cfa4x.2022.pdf" ) 
  include_graphics(c(fn1, fn2, fn3, fn4, fn5, fn6, fn7, fn8, fn9) )
# \@ref(fig:sizefeq-male-survey-cc)
```
 

 

 
```{r sizefeq-male, out.width='90%', echo=FALSE, fig.align='center', fig.cap = 'Size-frequency (areal density; no/km$^2$) histograms by carapace width of male Snow Crab. The vertical line represents the legal size (95 mm). Immature animals are shown with light coloured bars, mature with dark.'}
include_graphics(  file.path( SCD, "assessments", year.assessment, "figures", "size.freq", "survey",  "male.denl.png" )  )
# \@ref(fig:sizefeq-male)
```
  

```{r sizefeq-female, out.width='90%', echo=FALSE, fig.align='center', fig.cap = 'Size-frequency (areal density; no/km$^2$) histograms by carapace width of female Snow Crab. Immature animals are shown with light coloured bars, mature with dark.'}
include_graphics(  file.path( SCD, "assessments", year.assessment, "figures", "size.freq", "survey",  "female.denl.png" )  )
# \@ref(fig:sizefeq-female)
```
 
\clearpage


```{r fmat-timeseries, out.width='60%', echo=FALSE, fig.align='center', fig.cap = 'Mature female density log$_{10}$(no/km$^2$)  from the Snow Crab survey.'  }
include_graphics( file.path( SCD, "assessments", year.assessment, "timeseries", "survey", "totno.female.mat.pdf") )
# \@ref(fig:fmat-timeseries)
```
 

```{r fmat-map, echo=FALSE, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Mature female density log$_{10}$(no/km$^2$)  from the Snow Crab survey.' }
loc = file.path( SCD, "output", "maps", "survey", "snowcrab", "annual", "totno.female.mat" )
yrsplot = setdiff( year.assessment + c(0:-9), 2020)
fn4 = file.path( loc, paste( "totno.female.mat", yrsplot[4], "png", sep=".") )
fn3 = file.path( loc, paste( "totno.female.mat", yrsplot[3], "png", sep=".") )
fn2 = file.path( loc, paste( "totno.female.mat", yrsplot[2], "png", sep=".") )
fn1 = file.path( loc, paste( "totno.female.mat", yrsplot[1], "png", sep=".") )
include_graphics( c( fn3, fn2, fn1) )
# \@ref(fig:fmat-map)  
```



```{r sexratio-mature, out.width='60%', echo=FALSE, fig.align='center', fig.cap = 'Timeseries of sex ratios (proportion female) of mature Snow Crab.'  }
include_graphics( file.path( SCD, "assessments", year.assessment, "timeseries", "survey", "sexratio.mat.pdf") )
# \@ref(fig:sexratio-mature)
```
 
  


```{r sexratio-map, echo=FALSE, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Map of sex ratios (proportion female) of mature Snow Crab.'}
yrsplot = setdiff( year.assessment + c(0:-4), 2020)
loc = file.path( SCD, "output", "maps", "survey", "snowcrab", "annual", "sexratio.mat" )
fn4 = file.path( loc, paste( "sexratio.mat", yrsplot[4], "png", sep=".") )
fn3 = file.path( loc, paste( "sexratio.mat", yrsplot[3], "png", sep=".") )
fn2 = file.path( loc, paste( "sexratio.mat", yrsplot[2], "png", sep=".") )
fn1 = file.path( loc, paste( "sexratio.mat", yrsplot[1], "png", sep=".") )
include_graphics( c( fn3, fn2, fn1 ) )
# \@ref(fig:sexratio-map) 
```


\clearpage


```{r fbGMTS, out.width='60%', echo=FALSE, fig.align='center', fig.cap = 'The crude, unadjusted geometric mean fishable biomass density log~10(t/km$^2$) from the Snow Crab survey. Error bars represent 95\\% Confidence Intervals. Note the absence of data in 2020. Prior to 2004, surveys were conducted in the Spring.'}
fn = file.path(SCD,'assessments',year.assessment,'timeseries','survey','R0.mass.pdf')
include_graphics( c(fn) )
#\@ref(fig:fbGMTS)
```


```{r fbgeomean-map, echo=FALSE, out.width='45%', fig.show='hold', fig.align='center', fig.cap= 'Snow Crab survey fishable component biomass density log~10(t/km$^2$). Note, there is no data in 2020.' }
loc = file.path( SCD, 'output', 'maps', 'survey', 'snowcrab', 'annual', 'R0.mass')
yrsplot =  setdiff(year.assessment + c(0:-9), 2020 ) 
fn6 = file.path( loc, paste( 'R0.mass', yrsplot[6], 'png', sep='.') )
fn5 = file.path( loc, paste( 'R0.mass', yrsplot[5], 'png', sep='.') )
fn4 = file.path( loc, paste( 'R0.mass', yrsplot[4], 'png', sep='.') )
fn3 = file.path( loc, paste( 'R0.mass', yrsplot[3], 'png', sep='.') )
fn2 = file.path( loc, paste( 'R0.mass', yrsplot[2], 'png', sep='.') )
fn1 = file.path( loc, paste( 'R0.mass', yrsplot[1], 'png', sep='.') )
include_graphics( c(fn6, fn5, fn4, fn3, fn2, fn1) )
# \@ref(fig:fbgeomean-map)  
```
 

```{r fbindex-map, echo=FALSE, out.width='45%', fig.show='hold', fig.align='center', fig.cap= 'Biomass index log~10(t/km$^2$)  predicted from the Snow Crab survey.' }
loc = file.path( SCD, 'modelled', 'default_fb', 'predicted_biomass_densities' )
yrsplot =  year.assessment + c(0:-10)
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
include_graphics( c( fn8, fn7, fn6, fn5, fn4, fn3, fn2, fn1) )
```





```{r fbindex-timeseries, out.width='80%', echo=FALSE, fig.align='center', fig.cap = 'The fishable biomass index (t) predicted by CARSTM of Snow Crab survey densities. Error bars represent Bayesian 95\\% Credible Intervals. Note large errors in 2020 when there was no survey.' }
include_graphics( file.path( SCD, 'modelled', 'default_fb', 'aggregated_biomass_timeseries',  'biomass_M0.png') )
# \@ref(fig:fbindex-timeseries)
```



```{r logisticPredictions, out.width='32%', echo=FALSE, fig.show='hold', fig.align='center', fig.cap = 'Model 1 fishable, posterior mean modelled biomass (pre-fishery; kt) are shown in dark orange for N-ENS, S-ENS and 4X (left, middle and right). Light orange are posterior samples of modelled biomass (pre-fishery; kt) to illustrate the variability of the predictions. The biomass index (post-fishery, except prior to 2004) after model adjustment by the model catchability coefficient is in gray.' } 
loc = file.path( SCD, 'fishery_model', year.assessment, 'logistic_discrete_historical' )
fn1 = file.path( loc, 'plot_predictions_cfanorth.pdf' ) 
fn2 = file.path( loc, 'plot_predictions_cfasouth.pdf' ) 
fn3 = file.path( loc, 'plot_predictions_cfa4x.pdf' ) 
knitr::include_graphics(c(fn1, fn2, fn3) )
# \@ref(fig:logisticPredictions)
```



 

```{r logisticFishingMortality, out.width='32%', echo=FALSE,  fig.show='hold', fig.align='center', fig.cap = 'Time-series of modelled instantaneous fishing mortality from Model 1, for N-ENS (left), S-ENS (middle), and 4X (right). Samples of the posterior densities are presented, with the darkest line being the mean.' }
  odir = file.path( SCD, "fishery_model", year.assessment, "logistic_discrete_historical" )
  fn1 = file.path( odir, "plot_fishing_mortality_cfanorth.pdf" ) 
  fn2 = file.path( odir, "plot_fishing_mortality_cfasouth.pdf" ) 
  fn3 = file.path( odir, "plot_fishing_mortality_cfa4x.pdf" ) 
include_graphics(c(fn1, fn2, fn3) )
# \@ref(fig:logisticFishingMortality)
```



```{r ReferencePoints, out.width='60%', echo=FALSE,   fig.align='center', fig.cap = 'Harvest control rules for the Scotian Shelf Snow Crab fisheries.  N-ENS (left), S-ENS (middle), and 4X (right).' }
include_graphics( file.path( params$media_loc, 'harvest_control_rules.png') ) 
# \@ref(fig:ReferencePoints)
```
 
 

```{r logistic-hcr, out.width='32%', echo=FALSE, fig.show='hold', fig.align='center', fig.cap = 'Reference Points (fishing mortality and modelled biomass) from Model 1, for N-ENS (left), S-ENS (middle), and 4X (right). The large yellow dot indicates most recent year.' }
  odir = file.path( SCD, 'fishery_model', year.assessment, 'logistic_discrete_historical' )
  fn1 = file.path( odir, 'plot_hcr_cfanorth.pdf' ) 
  fn2 = file.path( odir, 'plot_hcr_cfasouth.pdf' ) 
  fn3 = file.path( odir, 'plot_hcr_cfa4x.pdf' ) 
  include_graphics(c(fn1, fn2, fn3) )
#  \@ref(fig:logistic-hcr)
```
   
\clearpage


```{r bottom-temperatures-survey, out.width='60%', echo=FALSE, fig.align='center', fig.cap = 'Annual variations in bottom temperature observed during the Snow Crab survey. The horizontal (black) line indicates the long-term, median temperature within each subarea. Error bars represent standard errors.' }
include_graphics( file.path( SCD, 'assessments', year.assessment, 'timeseries', 'survey', 't.pdf') )
# \@ref(fig:bottom-temperatures-survey)
```
 
 

```{r bottom-temperatures, out.width='80%', echo=FALSE, fig.align='center', fig.cap = 'Temporal variations in bottom temperature estimated from a historical analysis of temperature data. Red horizontal line is at $7^\\circ$C. Presented are 95\\% Credible Intervals of spatial variability in temperature at each time slice, after adjustment for spatiotemporal autocorrelation' }
include_graphics( file.path( SCD, 'assessments', year.assessment, 'timeseries', 'temperature_bottom.pdf') )
# \@ref(fig:bottom-temperatures)
```
 

```{r bottom-temperatures-spatialeffect, out.width='200%', echo=FALSE, fig.align='center', fig.cap = 'Persistent spatial effect of bottom temperature, relative to the overall mean, after adjustment for spatiotemporal variability and autocorrelations. Time period from 1999 to present.' }
loc = file.path( data_root, 'aegis', 'temperature', 'modelled', 'default', 'maps' )
include_graphics( file.path( loc, 'space_re_total.png') )
# \@ref(fig:bottom-temperatures-spatialeffect)
```
 
 

```{r bottom-temperatures-map, out.width='45%', echo=FALSE, fig.show='hold', fig.align='center', fig.cap = 'Spatial variations in bottom temperature estimated from a historical analysis of temperature data for 1 September.' }

loc = file.path( data_root, 'aegis', 'temperature', 'modelled', 'default', 'maps' )
yrsplot =  year.assessment + c(0:-10)
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
include_graphics( c( fn6, fn5, fn4, fn3, fn2, fn1) )
# \@ref(fig:bottom-temperatures-map)
# *Spatial variations in bottom temperature estimated from a historical analysis of temperature data for 1 September.*
```




\clearpage


```{r fb-habitat-timeseries, out.width='80%', echo=FALSE,   fig.align='center', fig.cap = 'Habitat viability (probability; fishable Snow Crab). Means and 95\\% Credible Intervals are presented.' }
loc = file.path( SCD, 'modelled', 'default_fb', 'aggregated_habitat_timeseries' )
include_graphics( file.path( loc, 'habitat_M0.png') )
# \@ref(fig:fb-habitat-timeseries)
```
 

 

```{r fb-habitat-map, out.width='45%', fig.show='hold', fig.align='center', fig.cap= 'Habitat viability (probability; fishable Snow Crab).' }
loc = file.path( SCD, 'modelled', 'default_fb', 'maps' )
yrsplot =  year.assessment + c(0:-10)
fn10 = file.path( loc, paste( 'Predicted_presence_absence.predictions.f', yrsplot[10], '.png', sep='') )
fn9 = file.path( loc, paste( 'Predicted_presence_absence.predictions.', yrsplot[9], '.png', sep='') )
fn8 = file.path( loc, paste( 'Predicted_presence_absence.predictions.', yrsplot[8], '.png', sep='') )
fn7 = file.path( loc, paste( 'Predicted_presence_absence.predictions.', yrsplot[7], '.png', sep='') )
fn6 = file.path( loc, paste( 'Predicted_presence_absence.predictions.', yrsplot[6], '.png', sep='') )
fn5 = file.path( loc, paste( 'Predicted_presence_absence.predictions.', yrsplot[5], '.png', sep='') )
fn4 = file.path( loc, paste( 'Predicted_presence_absence.predictions.', yrsplot[4], '.png', sep='') )
fn3 = file.path( loc, paste( 'Predicted_presence_absence.predictions.', yrsplot[3], '.png', sep='') )
fn2 = file.path( loc, paste( 'Predicted_presence_absence.predictions.', yrsplot[2], '.png', sep='') )
fn1 = file.path( loc, paste( 'Predicted_presence_absence.predictions.', yrsplot[1], '.png', sep='') )
include_graphics( c( fn8, fn7, fn6, fn5, fn4, fn3, fn2, fn1) )
# \@ref(fig:fb-habitat-map)  
# *Figure XXX. Habitat viability (probability; fishable Snow Crab)* 
```


\clearpage


```{r halibut-timeseries, out.width='50%', echo=FALSE,   fig.align='center', fig.cap = 'Atlantic Halibut crude, unadjusted geometric mean numerical density (no/km$^2$) from annual Snow Crab survey. Error bars are 95\\%  Confidence Intervals.' }
include_graphics( file.path( SCD, 'assessments', year.assessment, 'timeseries', 'survey', 'ms.no.30.pdf') )
# \@ref(fig:halibut-timeseries)
```
 

```{r halibut-map, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Halibut density log10(no/km$^2$) from the Snow Crab survey.' }
loc = file.path( SCD, 'output', 'maps', 'survey', 'snowcrab', 'annual', 'bycatch', 'ms.no.30' )
yrsplot = setdiff( year.assessment + c(0:-9), 2020)
fn4 = file.path( loc, paste( 'ms.no.30', yrsplot[4], 'png', sep='.') )
fn3 = file.path( loc, paste( 'ms.no.30', yrsplot[3], 'png', sep='.') )
fn2 = file.path( loc, paste( 'ms.no.30', yrsplot[2], 'png', sep='.') )
fn1 = file.path( loc, paste( 'ms.no.30', yrsplot[1], 'png', sep='.') )
include_graphics( c(  fn3, fn2, fn1) )
# \@ref(fig:halibut-map)  
```



\clearpage


```{r thornyskate-timeseries, out.width='60%', echo=FALSE,  fig.align='center', fig.cap = 'Thorny Skate crude, unadjusted geometric mean numerical density (no/km$^2$) from annual Snow Crab survey. Error bars are 95\\%  Confidence Intervals.'}
include_graphics( file.path( SCD, 'assessments', year.assessment, 'timeseries', 'survey', 'ms.no.201.pdf') )
# \@ref(fig:thornyskate-timeseries)
```

 

```{r thornyskate-map, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Thorny skate density log10(no/km$^2$) from the Snow Crab survey.' }
loc = file.path( SCD, 'output', 'maps', 'survey', 'snowcrab', 'annual', 'bycatch', 'ms.no.201' )
yrsplot = setdiff( year.assessment + c(0:-9), 2020)
fn4 = file.path( loc, paste( 'ms.no.201', yrsplot[4], 'png', sep='.') )
fn3 = file.path( loc, paste( 'ms.no.201', yrsplot[3], 'png', sep='.') )
fn2 = file.path( loc, paste( 'ms.no.201', yrsplot[2], 'png', sep='.') )
fn1 = file.path( loc, paste( 'ms.no.201', yrsplot[1], 'png', sep='.') )
include_graphics( c(  fn3, fn2, fn1) )
# \@ref(fig:thornyskate-map)  
```


\clearpage


```{r Wolffish-timeseries, out.width='60%', echo=FALSE,   fig.align='center', fig.cap = 'Striped Atlantic Wolffish crude, unadjusted geometric mean numerical density (no/km$^2$) from annual Snow Crab survey. Error bars are 95\\%  Confidence Intervals.' }
include_graphics( file.path( SCD, 'assessments', year.assessment, 'timeseries', 'survey', 'ms.no.50.pdf') )
# \@ref(fig:Wolffish-timeseries)
```
 

```{r Wolffish-map, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Striped Atlantic Wolffish density log10(no/km$^2$) from the Snow Crab survey.' }
loc = file.path( SCD, 'output', 'maps', 'survey', 'snowcrab', 'annual', 'bycatch', 'ms.no.50' )
yrsplot = setdiff( year.assessment + c(0:-9), 2020)
fn4 = file.path( loc, paste( 'ms.no.50', yrsplot[4], 'png', sep='.') )
fn3 = file.path( loc, paste( 'ms.no.50', yrsplot[3], 'png', sep='.') )
fn2 = file.path( loc, paste( 'ms.no.50', yrsplot[2], 'png', sep='.') )
fn1 = file.path( loc, paste( 'ms.no.50', yrsplot[1], 'png', sep='.') )
include_graphics( c( fn3, fn2, fn1) )
# \@ref(fig:Wolffish-map)  
```



\clearpage


```{r Shrimp-timeseries, out.width='60%', echo=FALSE,   fig.align='center', fig.cap = 'Northern Shrimp crude, unadjusted geometric mean numerical density (n/$km^2$) from annual Snow Crab survey. Error bars are 95\\%  Confidence Intervals.' }
include_graphics( file.path( SCD, 'assessments', year.assessment, 'timeseries', 'survey', 'ms.no.2211.pdf') )
# \@ref(fig:Shrimp-timeseries)
```
 

```{r Shrimp-map, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Northern Shrimp density log10(no/km$^2$) from the Snow Crab survey.' }
loc = file.path( SCD, 'output', 'maps', 'survey', 'snowcrab', 'annual', 'bycatch', 'ms.no.2211' )
yrsplot = setdiff( year.assessment + c(0:-9), 2020)
fn4 = file.path( loc, paste( 'ms.no.2211', yrsplot[4], 'png', sep='.') )
fn3 = file.path( loc, paste( 'ms.no.2211', yrsplot[3], 'png', sep='.') )
fn2 = file.path( loc, paste( 'ms.no.2211', yrsplot[2], 'png', sep='.') )
fn1 = file.path( loc, paste( 'ms.no.2211', yrsplot[1], 'png', sep='.') )
include_graphics( c( fn3, fn2, fn1) )
# \@ref(fig:Shrimp-map)  
```

 

\clearpage


```{r lessertoadcrab-timeseries, out.width='60%', echo=FALSE,   fig.align='center', fig.cap = 'Lesser Toad Crab crude, unadjusted geometric mean numerical density (no/km$^2$) from annual Snow Crab survey. Error bars are 95\\%  Confidence Intervals.' }
include_graphics( file.path( SCD, 'assessments', year.assessment, 'timeseries', 'survey', 'ms.no.2521.pdf') )
# \@ref(fig:lessertoadcrab-timeseries)
```
 

```{r lessertoadcrab-map, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Lesser Toad Crab density log10(no/km$^2$) from the Snow Crab survey.' }
loc = file.path( SCD, 'output', 'maps', 'survey', 'snowcrab', 'annual', 'bycatch', 'ms.no.2521' )
yrsplot = setdiff( year.assessment + c(0:-9), 2020)
fn4 = file.path( loc, paste( 'ms.no.2521', yrsplot[4], 'png', sep='.') )
fn3 = file.path( loc, paste( 'ms.no.2521', yrsplot[3], 'png', sep='.') )
fn2 = file.path( loc, paste( 'ms.no.2521', yrsplot[2], 'png', sep='.') )
fn1 = file.path( loc, paste( 'ms.no.2521', yrsplot[1], 'png', sep='.') )
include_graphics( c( fn3, fn2, fn1) )
# \@ref(fig:lessertoadcrab-map)  
```




# Appendix 

## Counts

```
     yr region Nstation         vessel Ntemp Nsp Ncrab
 1: 1996  S-ENS       24         OPILIO    24  34  2972
 2: 1997  N-ENS       26 GLACE BAY LADY    26  40  3216
 3: 1997  S-ENS      124 GLACE BAY LADY   124  46 11365
 4: 1998  N-ENS       56 MARCO BRITTANY    56  47  3881
 5: 1998  S-ENS      157 MARCO BRITTANY   157  56 13272
 6: 1999  N-ENS       59 MARCO BRITTANY    59  44  4091
 7: 1999  S-ENS      215 MARCO BRITTANY   215  55 10305
 8: 2000  N-ENS       69 MARCO BRITTANY    69  47  1968
 9: 2000  S-ENS      253 MARCO BRITTANY   253  60 11718
10: 2001  N-ENS       99 MARCO BRITTANY    99  60  1995
11: 2001  N-ENS       99   MARCO MICHEL    99  60  1995
12: 2001  S-ENS      257   MARCO MICHEL   257  71 10674
13: 2001  S-ENS      257 MARCO BRITTANY   257  71 10674
14: 2001  S-ENS      257  DEN C. MARTIN   257  71 10674
15: 2001     4X        3  DEN C. MARTIN     3  16    20
16: 2002  N-ENS      177   MARCO MICHEL   177  60  2539
17: 2002  S-ENS      231   MARCO MICHEL   231  63  6917
18: 2002     4X       21   MARCO MICHEL    21  33   110
19: 2003  N-ENS       79   MARCO MICHEL    79  47  2035
20: 2003  S-ENS      274   MARCO MICHEL   274  65  8368
21: 2003     4X        3   MARCO MICHEL     3  12     5
22: 2004  N-ENS       61    GENTLE LADY    61  68  1850
23: 2004  S-ENS      292    GENTLE LADY   292 100 13877
24: 2004     4X       26    GENTLE LADY    26  48   246
25: 2005  N-ENS       63    GENTLE LADY    63  76  4284
26: 2005  S-ENS      294    GENTLE LADY   294 111 18736
27: 2005     4X       32    GENTLE LADY    32  51   782
28: 2006  N-ENS       63    GENTLE LADY    63  65  8676
29: 2006  S-ENS      281    GENTLE LADY   281  90 26337
30: 2006     4X       30    GENTLE LADY    30  56  2083
31: 2007  N-ENS       65    GENTLE LADY    65  74  9027
32: 2007  S-ENS      284    GENTLE LADY   284 108 23573
33: 2007     4X       29    GENTLE LADY    29  57  1320
34: 2008  N-ENS       71    GENTLE LADY    71  82  7953
35: 2008  S-ENS      300    GENTLE LADY   300 117 22278
36: 2008     4X       34    GENTLE LADY    34  72  1084
37: 2009  N-ENS       71    GENTLE LADY    71  81  4360
38: 2009  S-ENS      303    GENTLE LADY   303 114 22430
39: 2009     4X       33    GENTLE LADY    33  72  2377
40: 2010  N-ENS       70    GENTLE LADY    70  80  4339
41: 2010  S-ENS      306    GENTLE LADY   306 102 23671
42: 2010     4X       31    GENTLE LADY    31  65  3188
43: 2011  N-ENS       72    GENTLE LADY    72  78  3141
44: 2011  S-ENS      310    GENTLE LADY   310 100 20464
45: 2011     4X       32    GENTLE LADY    32  65  1247
46: 2012  N-ENS       71    GENTLE LADY    71  75  2081
47: 2012  S-ENS      304    GENTLE LADY   304  97 17144
48: 2012     4X       32    GENTLE LADY    32  63  1130
49: 2013  N-ENS       71    GENTLE LADY    71  72  3262
50: 2013  S-ENS      304    GENTLE LADY   304 106 18003
51: 2013     4X       33    GENTLE LADY    33  67   450
52: 2014  N-ENS       65    MISS JESSIE    65  83  4399
53: 2014  S-ENS      281    MISS JESSIE   281 100 18768
54: 2014     4X       13    MISS JESSIE    13  59   345
55: 2015  N-ENS       71    MISS JESSIE    71  90  4226
56: 2015  S-ENS      308    MISS JESSIE   308 112 14967
57: 2015     4X       34    MISS JESSIE    34  66   556
58: 2016  N-ENS       70    MISS JESSIE    70  79  4605
59: 2016  S-ENS      307    MISS JESSIE   307 119 15620
60: 2016     4X       34    MISS JESSIE    34  63   469
61: 2017  N-ENS       70    MISS JESSIE    70  83  3637
62: 2017  S-ENS      266    MISS JESSIE   266 109 11573
63: 2017     4X       28    MISS JESSIE    28  59   243
64: 2018  N-ENS       66    MISS JESSIE    66  76  3212
65: 2018  S-ENS      298    MISS JESSIE   298 111 11635
66: 2018     4X       20    MISS JESSIE    20  51   583
67: 2019  N-ENS       68    MISS JESSIE    68  80  5262
68: 2019  S-ENS      342    MISS JESSIE   342 108 21602
69: 2019     4X       20    MISS JESSIE    20  49   971
70: 2021  N-ENS       68    MISS JESSIE    68  81  3717
71: 2021  S-ENS      293    MISS JESSIE   293 131 21411
72: 2021     4X       19    MISS JESSIE    19  50   778
73: 2022  N-ENS       68 R.S.JOURNEY II    68 114  1466
74: 2022  S-ENS      214 R.S.JOURNEY II   214 140 16203
75: 2022     4X       20 R.S.JOURNEY II    20  69   726
```