---
title: "Snow crab fishery model"
keywords: 
  - snow crab fishery model results 
abstract: |
  Snow crab fishery model summary.
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

Summary 4 of 4 -- This file is designed to be an HTML document that describes and summarizes the assessment of stock status. 
 
make quarto FN=06_fishery_model_summary.md YR=2024 DATADIR=~/bio.data/bio.snowcrab DOCTYPE=html PARAMS="-P year_assessment:2024" --directory=~/bio/bio.snowcrab/inst/markdown

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
  toget = c( "fishery_results", "fishery_model" )  

```


<!-- 
# _load_results.qmd contains instructions to load data 
#  this is a shared R-script to boot strap and provide a consistent data interface
-->

{{< include _load_results.qmd >}}  



## Fishery model


$$\begin{aligned}
b_{t+1} & \sim N\left({{{b_{t}+r}b_{t}}{{({1-b_{t}})}-\mathit{Fishing}}_{t}K^{-1},\sigma_{p}}\right) \\
Y_{t} & \sim N\left({q^{-1}\mathit{Kb}_{t},\sigma_{o}}\right) \\
r & \sim N\left({1.0,0.1}\right) \\
K & \sim N\left({\kappa,0.1\cdot\kappa}\right) \\
q & \sim N\left({1.0,0.1}\right) \\
\end{aligned}$$

$\kappa={\lbrack{5.0,60,1.25}\rbrack}$

${b=B}K^{-1}$,  a real unobservable (latent) process

- MCMC (Julia/Turing):
  - 4 chains, NUTS sampler (warmup: 10000, retain: 10000)
  - target_acceptance_rate, max_depth, init_ϵ = 0.75, 9, 0.05



## Surplus (Shaeffer) production

```{r}
#| label: fig-logistic-surplus-production
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

fns = paste("plot_surplus_production_", regions, ".png", sep="" )
fns = file.path( fm_loc, fns ) 

include_graphics( fns )
   
``` 


## Carrying capacity 

```{r}
#| label: fig-logistic-prior-posterior-K
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

## Intrinsic rate of increase

```{r}
#| label: fig-logistic-prior-posterior-r
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

## Catchability coefficient
 
```{r}
#| label: fig-logistic-prior-posterior-q
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


## Observation error

```{r}
#| label: fig-logistic-prior-posterior-obserror
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
 
## Model process error

```{r}
#| label: fig-logistic-prior-posterior-processerror
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

## State space 
  
```{r}
#| label: fig-logistic-state-space
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


----------- The text will need to be updated ------------------

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




  
