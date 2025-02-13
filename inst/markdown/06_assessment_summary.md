---
title: "Snow crab stock assessment"
author:
  - name: 
      given: Snow Crab Unit
      family: DFO Science
    # orcid: 0000-0003-3632-5723 
    # email: jae.choi@dfo-mpo.gc.ca
    # email: choi.jae.seok@gmail.com
    # corresponding: true
    affiliation: 
      - name: Bedford Institute of Oceanography, Fisheries and Oceans Canada
        city: Dartmouth
        state: NS
        # url: www.bio.gc.ca
date: last-modified
date-format: "YYYY-MM-D"
keywords: 
  - snow crab stock status assessment 
abstract: |
  Snow crab stock status assessment.
toc: true
toc-depth: 4
number-sections: true
highlight-style: pygments
# bibliography: media/references.bib  
# csl: media/canadian-journal-of-fisheries-and-aquatic-sciences.csl  # see https://www.zotero.org/styles for more
# license: "CC BY"
copyright: 
  holder: snow-crab-unit
  year: 2024
# citation: 
#  container-title: https://github.com/jae0/bio.snowcrab/
#  doi: NA
funding: "The snow crab scientific survey was funded by the snow crab fishers of Maritimes Region of Atlantic Canada."
editor:
  render-on-save: false
execute:
  echo: true
format:
  html: 
    code-fold: true
    code-overflow: wrap
    html-math-method: katex
    self-contained: true
    embed-resources: true
params:
  year_assessment: 2024
  year_start: 1999
  data_loc:  "~/bio.data/bio.snowcrab" 
  sens: 1
  debugging: FALSE
  model_variation: logistic_discrete_historical
---


 
<!-- 

## Preamble

Summary 4 of 4 -- This file is designed to be an HTML document that describes and summarizes the assessment of stock status. 
 

make quarto FN=06_assessment_summary.md YR=2024 DATADIR=~/bio.data/bio.snowcrab DOCTYPE=html PARAMS="-P year_assessment:2024" --directory=~/bio/bio.snowcrab/inst/markdown



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
  # dev.args = list(type = "cairo"),
  fig.retina = 2,
  dpi=192
)

require(spsUtil)
quietly = spsUtil::quiet

require(ggplot2)
require(MBA)
require(gt)  # table formatting
require(aegis)  # basic helper tools
    
loadfunctions( "aegis")
loadfunctions( "bio.snowcrab")  # in case of local edits

year_assessment = params$year_assessment
model_variation = params$model_variation

data_loc= params$data_loc
media_loc = file.path( params$media_loc, "media" )


p = load.environment( year.assessment=year_assessment )
p$corners = data.frame(plon=c(220, 990), plat=c(4750, 5270) )
p$mapyears = year_assessment + c(-5:0 )   # default in case not specified

year_previous = year_assessment - 1

p$fishery_model_years = 2000:year_assessment   ### NOTE these years need to be consistent with years in the 
 
lregions = list(region=c("cfanorth", "cfasouth", "cfa4x"))
reg_labels = c("N-ENS", "S-ENS", "CFA 4X")  # formatted for label

if (params$sens==2) {
  lregions = list(subarea=c("cfanorth", "cfa23",  "cfa24", "cfa4x"))
  reg_labels = c("CFA 20-22", "CFA 23", "CFA 24", "CFA 4X")  # formatted for label
}

regions = unlist(lregions)
nregions = length(regions)

# directories
outtabledir = file.path( p$annual.results, "tables" )

# fishery_model_results = file.path( "/home", "jae", "projects", "dynamical_model", "snowcrab", "outputs" )
fishery_model_results = file.path( data_loc, "fishery_model" )
   
fm_loc = file.path( data_loc, 'fishery_model', year_assessment, model_variation )
  

# as modelled years in fishery model can differ from iput data years, make sure  "years_model" is correct
p$fishery_model_years = 2000:year_assessment
sn_env = snowcrab_load_key_results_to_memory( year_assessment, years_model=p$fishery_model_years, return_as_list=TRUE  ) 

attach(sn_env)
 

```


&nbsp;  $~$  <br /> 


# Stock status

- Survey timing changed from Spring to Autumn in 2004.

- No survey in 2020 (Covid-19) and incomplete 2022 (mechanical issues).

  - Inshore areas of S-ENS were most affected.
  - N-ENS and CFA 4X were not affected.


## Size structure

Factors: early maturation, size-selective predation, fishing of largest individuals, start or end of a recruitment pulse, timing of survey.
 
## Recruitment

- Little to no recruitment is expected for the next 1-3 years in N-ENS.
- Moderate levels of recruitment are expected in S-ENS.
- Low to moderate levels of recruitment are expected for 2 years in 4X.


## Reproduction

- All areas had recruitment of female crab into the mature (egg-bearing) segment of the population from 2016-2022.
- In N-ENS for 2022, a decline in numerical densities, and low densities of adolescent females.
- Egg and larval production is expected to be moderate to high in the next year in all areas except N-ENS.
 
## Reproduction (mature females)
 
Distributions are heterogeneous and often in shallower areas.
 
## Sex ratios (proportion female, mature)

- Mostly male-dominated: larger size may be protective against predation?
- Imbalance indicates differential mortality: predation, competition and fishing
- In 4X, sex ratios are balanced.
 

## Modelled Biomass (pre-fishery)

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
 
 


N-ENS: `{r} round(B_north[t0], 2)` t in `{r} year_assessment`

- `{r} round(B_north[t1], 2)` t in `{r} year_previous`.

S-ENS: `{r} round(B_south[t0], 2)` t in `{r} year_assessment`

- `{r} round(B_south[t1], 2)` t in `{r} year_previous`.

4X:  `{r} round(B_4x[t0], 2)` t in `{r} year_assessment`-`{r} year_assessment+1`

- `{r} round(B_4x[t1], 2)` t for the `{r} year_previous`-`{r} year_assessment` season.

 
## Fishing Mortality



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


N-ENS: `{r} round(FM_north[t0],3)` (annual exploitation rate of `{r} round(100*(exp(FM_north[t0])-1),2)`%) in `{r} year_assessment`

  - Up from the  `{r} year_previous` rate of `{r} round(FM_north[t1],3)` (annual exploitation rate of `{r} round(100*(exp(FM_north[t1])-1),1)`%)

S-ENS: `{r} round(FM_south[t0],3)` (annual exploitation rate of `{r} round(100*(exp(FM_south[t0])-1),1)`%) in `{r} year_assessment`

  - Decreasing marginally from the `{r} year_previous` rate of `{r} round(FM_south[t1],3)` (annual exploitation rate of `{r} round(100*(exp(FM_south[t1])-1),1)`%)

4X: `{r} round(FM_4x[t0],3)` (annual exploitation rate of `{r} round(100*(exp(FM_4x[t0])-1),1)`%) in `{r} year_assessment`-`{r} year_assessment+1` season

  - Decreasing from the `{r} year_assessment-1`-`{r} year_assessment` season rate of `{r} round(FM_4x[t1],3)` (annual exploitation rate of `{r} round(100*(exp(FM_4x[t1])-1),1)`%)

Localized exploitation rates are likely higher, as not all areas for which biomass is estimated are fished.





## Reference Points

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
#| fig-cap: "Reference Points (fishing mortality and modelled biomass) from the Fishery Model, for N-ENS (left), S-ENS (middle), and 4X (right). The large yellow dot indicates most recent year and the 95\\% CI. Not: the model does not account for illegal and unreported landings, and interspecific interactions."
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
 

|   | N-ENS | S-ENS | 4X |
|----- | ----- | ----- | ----- |
| |  |  |  |
|q       | `r round(q_north, 3)` (`r round(q_north_sd, 3)`) | `r round(q_south, 3)` (`r round(q_south_sd, 3)`) | `r round(q_4x, 3)` (`r round(q_4x_sd, 3)`) |
|r       | `r round(r_north, 3)` (`r round(r_north_sd, 3)`) | `r round(r_south, 3)` (`r round(r_south_sd, 3)`) | `r round(r_4x, 3)` (`r round(r_4x_sd, 3)`) |
|K       | `r round(K_north, 2)` (`r round(K_north_sd, 2)`) | `r round(K_south, 2)` (`r round(K_south_sd, 2)`) | `r round(K_4x, 2)` (`r round(K_4x_sd, 2)`) |
|Prefishery Biomass   | `r round(B_north[t0], 2)` (`r round(B_north_sd[t0], 2)`) | `r round(B_south[t0], 2)`  (`r round(B_south_sd[t0], 2)`) | `r round(B_4x[t0], 2)`  (`r round(B_4x_sd[t0], 2)`)  |
|Fishing Mortality    | `r round(FM_north[t0], 3)` (`r round(FM_north_sd[t0], 3)`) | `r round(FM_south[t0], 3)` (`r round(FM_south_sd[t0], 3)`) | `r round(FM_4x[t0], 3)` (`r round(FM_4x_sd[t0], 3)`) |

Note: Values in parentheses are Posterior standard deviations.


----------- The text needs to be updated ------------------

- N-ENS is in the "cautious" zone
- S-ENS is in the "healthy" zone
- 4X is in the "critical" zone



# Conclusions

The ESS ecosystem is still experiencing a lot of volatility and prudence is wise:

- Rapid ecosystem change (groundfish collapse in 1980s and 1990s)
- Rapid climatic change: variability high, very warm period from 2012-2024
- Snow crab:
  - Can persist if extreme conditions are episodic by shifting spatial distributions
  - 4X was possibly too long and did not have the habitat space
  - Climate variability can induce juxtaposition of warmer water species (prey and predators) with snow crab
  - Some shifts in spatial distribution towards cooler and deeper waters

Modelled solutions:

- A few of many possible views
- Over-emphasis of any one view and associated Reference Points is **not precautionary**.


## Conclusions: N-ENS

- Viable habitat stable near historical mean and marginally improved relative to 2022.
- Female larval production peaked in 2017 and has since declined.
- Recruitment continues at low levels, a noticible gap in recruitment to the fishery persists.
- Predation mortality (Cod, Halibut, Skate) is likely high and causing reduced recruitment to the fishery.
- Fishing mortality has been increasing since 2019.
- Fishable biomass continues a decreasing trend.
- PA template suggest "healthy zone".
- A more careful harvest strategy would enable bridging this coming recruitment gap.
- A reduced TAC is prudent.

In N-ENS, though recruitment continues at low levels,
a gap in future recruitment to the fishery is expected for the next 1-3 years


 bridging this coming recruitment gap. A reduced TAC is prudent.


## Conclusions: S-ENS

- Viable habitat marginally improved relative to a historic low in 2022.
- Female larval production peaked in 2023/2024.
- Recruitment to the fishery from a strong year class has begun, full entry by 2025.
- Predation mortality (Halibut, Plaice, Sculpin, Skate) is likely high and causing reduced recruitment to the fishery.
- Fishing mortality has been increasing since 2019.
- Fishable biomass is at a historical low.
- PA template suggest "healthy zone".

- Flexiblity in harvest strategy exists due to strong recuitment.
- A status quo TAC is prudent until confirmation of recruitment.

recruitment to the fishery is likely to continue at a moderate rate for the upcoming season


## Conclusions: 4X

- Viable habitat has been at historical lows since 2015.
- Female larval production peaked in 2022.
- Recruitment to the fishery in the near-term is not evident. There is a pulse 4-5 years away.
- Predation mortality (Halibut) and competition with other Crustacea is likely high and causing reduced recruitment to the fishery.
- Fishing mortality has declined to negligible levels since a peak in 2019.
- Fishable biomass is at a historical low.
- PA template suggest "critical zone".
- Stock status is extremely poor.
- Fishery closure until recovery is prudent.

   low to moderate levels of recruitment are expected for 2 years.
   4X exists in the "cautious zone".

   habitat has been depressed for many years. A reduced TAC is prudent.




  

## Surplus production ("Shaeffer form") {.c}

```{r}
#| label: fig-logistic-surplus-production
#| results: asis
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig-cap: "Surplus production. Each year is represented by a different colour."


for (r in 1:nregions){
  reg = regions[r]
  REG = reg_labels[r] 
  cat("#### ", REG, "\n")
  fn = file.path( fm_loc, paste("plot_surplus_production_", reg, ".png", sep="" ) ) 
  show_image(check_file_exists(fn)) 
  cat("\n\n")
} 
  
``` 





 
## Prior-posterior comparsons: Carrying capacity (K; kt) {.c}

```{r}
#| label: fig-logistic-prior-posterior-K
#| results: asis
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig-cap: "Prior-posterior comparisons of carrying capacity (K; kt)."

for (r in 1:nregions){
  reg = regions[r]
  REG = reg_labels[r] 
  cat("#### ", REG, "\n")
  fn = file.path( fm_loc, paste("plot_prior_K_", reg, ".png", sep="" ) ) 
  show_image(check_file_exists(fn)) 
  cat("\n\n")
} 
  

``` 




## Prior-posterior comparsons: Intrinsic rate of biomass increase (r) {.c}

```{r}
#| label: fig-logistic-prior-posterior-r
#| results: asis
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig-cap: "Prior-posterior comparisons of th iIntrinsic rate of biomass increase (r)."


for (r in 1:nregions){
  reg = regions[r]
  REG = reg_labels[r] 
  cat("#### ", REG, "\n")
  fn = file.path( fm_loc, paste("plot_prior_r_", reg, ".png", sep="" ) ) 
  show_image(check_file_exists(fn)) 
  cat("\n\n")
} 
  

``` 



 
## Prior-posterior comparsons: Catchability coefficient (q) {.c}

```{r}
#| label: fig-logistic-prior-posterior-q
#| results: asis
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig-cap: "Prior-posterior comparisons of the catchability coefficient (q)."


for (r in 1:nregions){
  reg = regions[r]
  REG = reg_labels[r] 
  cat("#### ", REG, "\n")
  fn = file.path( fm_loc, paste("plot_prior_q1_", reg, ".png", sep="" ) ) 
  show_image(check_file_exists(fn)) 
  cat("\n\n")
} 
  

``` 




 
## Prior-posterior comparsons: Observation error  {.c}

```{r}
#| label: fig-logistic-prior-posterior-obserror
#| results: asis
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig-cap: "Prior-posterior comparisons of Observation error (kt)."
 
for (r in 1:nregions){
  reg = regions[r]
  REG = reg_labels[r] 
  cat("#### ", REG, "\n")
  fn = file.path( fm_loc, paste("plot_prior_bosd_", reg, ".png", sep="" ) ) 
  show_image(check_file_exists(fn)) 
  cat("\n\n")
} 
  

``` 
 

 
## Prior-posterior comparsons: Model process error {.c}

```{r}
#| label: fig-logistic-prior-posterior-processerror
#| results: asis
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig-cap: "Prior-posterior comparisons of Model process error (kt)."
 

for (r in 1:nregions){
  reg = regions[r]
  REG = reg_labels[r] 
  cat("#### ", REG, "\n")
  fn = file.path( fm_loc, paste("plot_prior_bpsd_", reg, ".png", sep="" ) ) 
  show_image(check_file_exists(fn)) 
  cat("\n\n")
} 
  

``` 
 
 
## State space {.c}

```{r}
#| label: fig-logistic-state-space
#| results: asis
#| echo: false
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4
#| fig-cap: "State space (kt): year vs year+1."
 

for (r in 1:nregions){
  reg = regions[r]
  REG = reg_labels[r] 
  cat("#### ", REG, "\n")
  fn = file.path( fm_loc, paste("plot_state_space_", reg, ".png", sep="" ) ) 
  show_image(check_file_exists(fn)) 
  cat("\n\n")
} 
  

``` 
 

