---
title: "Snow Crab Assessment"
subtitle: "Life history and Methodological considerations"
author:
  - name: 
      given: Snow Crab Unit, DFO Science
    # family: Choi
    # orcid: 0000-0003-3632-5723 
    # email: jae.choi@dfo-mpo.gc.ca
    # email: choi.jae.seok@gmail.com
    # corresponding: true
    affiliation: 
      - name: Bedford Institute of Oceanography, Fisheries and Oceans Canada
        city: Dartmouth
        state: NS
        # url: www.bio.gc.ca
# date: "`r format(Sys.time(), '%d %B, %Y')`"
date: last-modified
date-format: "YYYY-MM-D"
toc: false
toc-depth: 4
number-sections: false
highlight-style: pygments
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
  docx:
    toc: true
    toc_depth: 2
    number-sections: false
    highlight-style: github
    tbl-cap-location: top
    reference-doc: media/csas_template.docx
    # reference-doc: media/RES2024-eng.docx
  revealjs:
    theme: moon    
    self-contained: true
    scrollable: true
    smaller: true
    dpi: 144
  html: 
    code-fold: true
    code-overflow: wrap
    html-math-method: katex
    self-contained: true
    embed-resources: true
  pdf:
    pdf-engine: lualatex
  beamer:
    theme: "metropolis"
    colortheme: "seagull"
    fonttheme: "professionalfonts"
    fig_caption: yes
    latex_engine: lualatex 
    keep_tex: true 
# classoption: 
#  - aspectratio=169 #16:9 wide
#  - t  # top align
# header-includes: 
#   - \usepackage{graphicx}
 #  - \usepackage[font={scriptsize}, labelfont={bf}]{caption}
 #  - \usepackage{longtable}
 # - \usepackage{booktabs}
 # - \usepackage{caption}
 # - \usepackage{float}
 # - \usepackage{multicol}
  # - \usepackage{subfig}
  # - \newcommand{\btiny}{\begin{tiny}}
  # - \newcommand{\etiny}{\end{tiny}}
#revealjs-plugins:
#  - revealjs-text-resizer
params:
  year_assessment: 2024
  year_start: 1999
  data_loc:  "~/bio.data/bio.snowcrab"
  sens: 1
  model_variation: logistic_discrete_historical
  debugging: FALSE
--- 


<!-- Preamble
  
  make quarto FN=snowcrab_presentation_methods.md DOCTYPE=revealjs  PARAMS="-P year_assessment:2024"  --directory=~/bio/bio.snowcrab/inst/markdown

  make quarto FN=snowcrab_presentation_methods.md DOCTYPE=html   PARAMS="-P year_assessment:2024"  --directory=~/bio/bio.snowcrab/inst/markdown
      
  make quarto FN=snowcrab_presentation_methods.md DOCTYPE=beamer  PARAMS="-P year_assessment:2024"  --directory=~/bio/bio.snowcrab/inst/markdown 

-->
 


<!-- Set up R-environment -->

```{r}
#| label: setup
#| eval: true 
#| output: false
#| echo: false

  require(knitr)

  opts_chunk$set(
    root.dir = data_root,
    echo = FALSE,
    out.width="6.2in",
    # dev.args = list(type = "cairo"),
    fig.retina = 2,
    dpi=192
  )
  
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

  require(ggplot2)
  require(MBA)
  require(aegis)  # basic helper tools

  loadfunctions( "aegis")
  loadfunctions( "bio.snowcrab")  # in case of local edits

  year_assessment = params$year_assessment
  
  year_start = params$year_start

  year_previous = year_assessment - 1

  model_variation = params$model_variation
  
 
data_loc= params$data_loc
media_loc = file.path( params$media_loc, "media" )

 
#### params and directories

  p = load.environment( year.assessment=year_assessment )  
  
  p$corners = data.frame(plon=c(220, 990), plat=c(4750, 5270) )

  p$mapyears = year_assessment + c(-5:0 )   # default in case not specified
  
  years = as.character(1996: year_assessment)
  yrs_observer = year_assessment + c(0:-4)
  
  fishery_model_results = file.path( data_loc, "fishery_model" )
 
  # note copied "roadshow figures" temporaily here ... figure creation should be be assimilated TODO
  media_supplementary = file.path( data_loc, "assessments",  year_assessment, "media_supplementary")

  outtabledir = file.path( p$annual.results, "tables" ) 

  lregions = list(region=c("cfanorth", "cfasouth", "cfa4x"))
  reg_labels = c("N-ENS", "S-ENS", "CFA 4X")  # formatted for label

  if (params$sens==2) {
    lregions = list(subarea=c("cfanorth", "cfa23",  "cfa24", "cfa4x"))
    reg_labels = c("CFA 20-22", "CFA 23", "CFA 24", "CFA 4X")  # formatted for label
  }

  vnr = names(lregions)
  regions = unname( unlist(lregions) )
  
  nregions = length(regions)
  
 
```
  


## Life history {.c}

```{r}
#| label: fig-photos-snowcrab-1
#| eval: true
#| echo: false 
#| output: true
#| fig-cap: "Snow Crab male, mature benthic form."
#| fig-subcap: 
#|   - ""
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| layout-ncol: 1

fns = file.path( media_loc, "snowcrab_male.png" )

include_graphics( fns ) 

```

&nbsp;  $~$  <br /> 

## {#Photographs data-menu-title="Photographs"} 

```{r}
#| label: fig-photos-snowcrab-2
#| eval: true
#| echo: false 
#| output: true
#| fig-cap: "Snow Crab images. Note sexual dimorphism."
#| fig-subcap: 
#|   - "Pelagic zoea"
#|   - "Mating pair - note sexual dimorphism, with a smaller female. Note epibiont growth on female which suggests it is an older female (multiparous)."
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| layout-ncol: 2

fns = file.path( media_loc, c(
  "snowcrab_zoea.png", 
  "snowcrab_male_and_female.png"
) )

include_graphics( fns ) 

```

&nbsp;  $~$  <br /> 
 


## Life history stages
 
```{r}
#| label: fig-lifehistory
#| eval: true
#| echo: false 
#| output: true
#| fig-cap: "Life history patterns of snow crab and approximate timing of the main life history stages of snow crab and size (carapace width; CW mm) and instar (Roman numerals). Size and timings are specific to the area of study and vary with environmental conditions, food availability and genetic variability."
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold

fn1=file.path( media_loc, "life_history.png" )
include_graphics( fn1 ) 

```

## Male growth stanzas 
 
```{r}
#| label: fig-growth-stanzas
#| eval: true
#| echo: false 
#| output: true
#| fig-cap: "The growth stanzas of the male component and decision paths to maturity and terminal moult. Black ellipses indicate terminally molted animals."
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
 
fn1=file.path( media_loc, "life_history_male.png" )
include_graphics( fn1 ) 

```

 
## Growth patterns inferred from size modes

```{r}
#| label: fig-growth-modes
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: "Modes (CW, ln mm) identified from survey data using *Kernel Density Estmation* of local moving data windows. Legend: sex|maturity|instar"
#| fig-subcap: 
#|   - "Female modes"
#|   - "Male modes" 


fns = file.path( media_loc, c(
  "density_f_imodes.png", 
  "density_m_imodes.png"
))

include_graphics( fns )
```
 
  
## Growth patterns inferred from modes
 
```{r}
#| label: fig-growth-modes-growth
#| eval: true
#| output: true
#| fig-dpi: 144
#| fig-height: 4 
#| echo: false 
#| layout-ncol: 2
#| fig-cap: "Inferred growth derived from *Kernel Mixture Models* (priors)."
#| fig-subcap: 
#|   - "Female growth trajectory"
#|   - "Male growth trajectory"

fns = file.path( media_loc, c(
  "plot_growth_female.png",
  "plot_growth_male.png"
))

include_graphics( fns )
```
 
 
&nbsp;  $~$  <br /> 

 

## Notable traits

  - Sexual dimorphism 
  - Pelagic in larval stages; benthic in pre-adolescent and adult stages
  - Biannual, annual molts depending upon size/age and environmental conditions
  - Terminal molt to maturity and survives for up to 5 years
  - Total life span: up to 15 years. 
  - Stenothermic (narrow) temperature requirements < 6 C
  - Ontogenetic shifts in habitat preferences: warmer, complex substrates to mud.  
  - Cannibalism of immature crab by mature female snow crab are known 
  - Primiparous female (57.4 mm CW) produces between 35,000 to 46,000 eggs
  - Multiparous females more fecund (>100,000 eggs) 
  - Eggs extruded February and April and brooded for up to two years (> 80% in SSE follow an annual cycle) and hatched/released from April to June.
   

## Life history: Spatially and temporally complex
  - 14+ instars/life stages with different preferences; 
  - Multiply prey and predators (ontogenetic shift); many ecosystem changes 
  - Growth: biannual (up to instar 5), annual, biennial (skip moulters)
  - Reproduction: annual and biennial; strong sexual selection (dimorphism) and parental investment for up to 2 years; lags of up to 10 years
  - Mortality: pulsed due to moult vulnerability, cannibalism
  - Movement: seasonal, annual, ontogenetic 
  - Evolution: sexual dimorphism, parental investment of brood, stenothermic
  - Population structure integrates all these effects ("Cumulative Effects")

## Population structure integrates everything
  - Spatial structures: scales > 100's km; local crab holes ( 10's to 100's m ) 
    - Confluence of major oceanic currents (temperature, salinity, nutrients, primary and secondary production)
  - Temporal processes: scales > 10 yrs; daily, tidal currents,  
    - 10-15 y life cycle, seasonal, annual, multiyear, decadal, evolutionary, etc. 
  - Biological (individual, age classes, population, etc.): cross multiple space, time and organismal scales 
    - Ontogeny (developmental change and habitat requirements)
    - Life history (habitat preferences)
    - Inter-specific interactions (feeding relationships change with biological stage)
    - Interactions with humans

## Multiple processes at multiple scales 

  - Populations integrate these processes  
  - **Induces** auto-correlations at multiple spatial-temporal-organisational scales
  - We perceive this as "autocorrelation"
  - First law of geography: "everything is related to everything else, but near things are more related than distant things" (Tobler 1970) 
  - "Nonstationary" processes (mean and variance not stable in space and time)
  - Ignoring correlated errors causes inappropriate parameter estimates


## Observations error (precision, accuracy) 

  - Pragmatic compromise between costs vs information gain
  - SSE snow crab domain with ~400 stations and ~109,120  $\text{km}^{2}$ 
  - Each station represents 273 $\text{km}^{2}$ , but actually samples ~0.0039 $\text{km}^{2}$ 
    - A factor of 1:70,000 (spatially)
  - Each sample is about 5 minutes but represents a whole year (525,600 minutes)
    - A factor of 1:105,120
  - In space-time: 
    - A factor of ~ 1:10 billion (note: 1 Angstrom = 10$^{-10}$m) .. the same scale as an atom! 
    - A few samples will do if it is spatially homogenous ("well mixed" = constant forcings, no ecosystem variability, no spatiotemporal autocorrelation)
  - Spatial bias/aliasing
  - Temporal bias/aliasing ("Spring" before 2004)
  

## Solutions: Experimental design

  - Naive (Ignore everything) "Random" Sampling / Fixed Station Sampling (IID)
  - Stratified (Hope this works) "Random" Sampling -- depth is all that matters (IID within, IID between)
    - Ignore spatiotemporal structure 
    - Ignore Ecosystem variability
    - Assume **no autocorrelated errors**
  - Geostatistical design (gridded, pseudo-random)
    - Address spatial autocorrelation as a stationary process
    - Abundance is NOT first and second order stationary and non-Gaussian (distributional bias)
    - Ignores time (temporal discretization/aliasing)
    - Ignores ecosystem variability
    - Variogram solutions unstable when population size declines

## Solutions: Experimental design 2

  - Hierarchical process temporal and covariate as an external forcing and spatial as a local process (UKED)
    - Variogram instability across time/stage when abundance declined and distribution was spotty
    - NOT second-order stationary
  - Mosaic process:
    - Hierarchical first and second-order **non-stationary** model (modular)
    - Global Hurdle/INLA/GAM/GLM 
    - Local Random spatial effects (residuals) Matern (SPDE; FFT)
    - Independent space and time process
    - Slow
  - Model-based inference: Bayesian hierarchical mixed-effects model
    - Conditional Autoregressive Space-Time Models (CARSTM, separable: Kronecker product space BYM2, time AR1)
    - Hurdle/Poisson and Binomial number process, with Gaussian weight process


## Latent ecological process

- Real (latent, unobserved) ecological processes are generaly spatiotemporal processes  

- Observations from samples/surveys are used to infer the real (latent, unobserved) state

- Stock assessments, almost always, focus only upon a purely *temporal process* by integrating the latent spatiotemporal process  
 
- Experimental design assumes samples are random or random-stratified, almost always, as a purely *spatial process*

- This is a problem: 
  - temporal aliasing is likely (as it is usually ignored)
  - spatial aliasing is likely when sampling environments that cannot be randomly sampled 
    - many locations are not directly accessible to trawls (rocks and bedrock)
    - easier to sample locations (softer, gravel or mud substrates) coincide with preferred snow habitats

 
## Snow crab survey locations {.c}
   
```{r survey-locations-map, out.width='60%', fig.show='hold', fig.align='center', fig.cap= 'Snow Crab survey locations. No survey in 2020 (Covid-19) and incomplete 2022 (mechanical issues).' }
loc = file.path( data_loc, "output", "maps", "survey.locations" )
yrsplot = setdiff( year_assessment + c(0:-9), 2020)
fn6 = file.path( loc, paste( "survey.locations", yrsplot[6], "png", sep=".") )
fn5 = file.path( loc, paste( "survey.locations", yrsplot[5], "png", sep=".") )
fn4 = file.path( loc, paste( "survey.locations", yrsplot[4], "png", sep=".") )
fn3 = file.path( loc, paste( "survey.locations", yrsplot[3], "png", sep=".") )
fn2 = file.path( loc, paste( "survey.locations", yrsplot[2], "png", sep=".") )
fn1 = file.path( loc, paste( "survey.locations", yrsplot[1], "png", sep=".") )
include_graphics( c( fn2, fn1) )
# \@ref(fig:survey-locations-map)  
``` 

 
## {}

```{r surveydomain, out.width='100%', fig.show='hold', fig.align='center', fig.cap= 'Snow Crab survey prediction grid.' }
fn1 = file.path( media_loc, "carstm_prediction_domain.png" )
include_graphics( c( fn1) )
``` 

 
## Spatial clustering  
 
```{r}
#| label: fig-aggregation
#| eval: true
#| echo: false 
#| output: true
#| fig-cap: "Spider crab tend to cluster/aggregate in space."
#| fig-subcap: 
#|   - "Australian spider crab \\emph{Leptomithrax gaimardii} aggregation for moulting and migration."
#|   - "Alaska red king crab \\emph{Paralithodes camtschaticus} aggregation in Alaska for egg release, migrations."
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
 
fns = file.path( media_loc, c( 
  "australian_leptomithrax_gaimardii.png", 
  "kingcrab_aggregation.png" 
) )

include_graphics( fns ) 

```

- *Other* crab species show "Mounding" for protection from predation (larval, moulting and females).

- Narrow habitat preferences force them to move and cluster when environment is poor.

## Spatial clustering 2 
 

```{r}
#| label: fig-aggregation2
#| eval: true
#| echo: false 
#| output: true
#| fig-cap: "High density locations of Snow Crab, approximately 1 per square meter."
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
 
fn = file.path(p$project.outputdir, "maps", "map_highdensity_locations.png" )

include_graphics( fn ) 

```


## Summary
 


:::: {.columns}

::: {.column width="50%"}

- Sexual dimorphism
- Ontogenetic shifts (changes with age and stage) in habitat preferences: 
  - Pelagic in larval stages, with diurnal veritical migration  
  - Benthic in pre-adolescent and adult stages
  - Complex substrates for small snow crab to muddy 
  - Adult females and smaller immature benthic crab have a preference for temperatures < 0 Celcius
  - Adult males have a preference for < 4 Celcius  
- Biannual, annual and biennial molts depending upon size/age and environmental conditions
- Terminal molt to maturity and survives for up to 5 years thereafter
:::

::: {.column width="50%"}
 
- Total life span: up to 13 years (females) and 15 years (males)
- Cannibalism of immature crab by mature female snow crab are known 
- Fecundity:
  - Primiparous - 57.4 mm CW can produce between 35,000 to 46,000 eggs
  - Multiparous - more than 100,000 eggs  
- Eggs extruded February and April and brooded for up to two years 
  - more than 80% in the Area of Interest follow an annual cycle and hatched/released from April to June.
- Movement and connectivity
  - pelagic and zoea stages -- no information: vertical migration related to temperature and light
  - adult females -- no information: possibly related to food availability, temperature and predator avoidance 
  - adult males is long-tailed in distribution and possibly related to mate finding and food availablity

:::

::::
 

## Sampling bias: depth
```{r bias-depth, out.width='40%', fig.show='hold', fig.align='center', fig.cap= 'Comparison of depths in snow crab domain ("predictions") and survey locations ("observations"). Bias: survey trawls preferentially sample deeper locations. Red line is overall average.' }

include_graphics( file.path( framework_loc, "bias_depth.png" ) )
```
 
## Sampling bias: substrate grain size

```{r bias-substrate, out.width='40%', fig.show='hold', fig.align='center', fig.cap= 'Comparison of substrate grain size (log; mm) in snow crab domain ("predictions") and survey locations ("observations"). No bias observed. Red line is overall average.' }

include_graphics( file.path( framework_loc, "bias_substrate.png" ) ) 
``` 

## Sampling bias: bottom temperature

```{r bias-temperature, out.width='40%', fig.show='hold', fig.align='center', fig.cap= 'Comparison of bottom temperature (degrees Celcius) in snow crab domain ("predictions") and survey locations ("observations"). Bias: surveys preferentially sample cold water bottoms. Red line is overall average.' }

include_graphics( file.path( framework_loc, "bias_temp.png" ) ) 
``` 
 
## Sampling bias: species composition 1 (temperature related)

```{r bias-pc1, out.width='40%', fig.show='hold', fig.align='center', fig.cap= 'Comparison of species composition gradient (PC1) in snow crab domain ("predictions") and survey locations ("observations"). Bias: surveys preferentially sample cold water species. Red line is overall average.' } 

include_graphics( file.path( framework_loc, "bias_pca1.png" ) ) 
``` 
  

## Sampling bias: species composition 2 (depth related)

```{r bias-pc2, out.width='40%', fig.show='hold', fig.align='center', fig.cap= 'Comparison of species composition gradient (PC2) in snow crab domain ("predictions") and survey locations ("observations"). Bias: surveys preferentially sample deeper species. Red line is overall average.' }

include_graphics( file.path( framework_loc, "bias_pca2.png" ) ) 
``` 

## Sampling bias: summary

Indications of sampling bias relative to spatial domain of snow crab:

- deeper and more favourable locations are being sampled
- statistically overlapping but biologically important
 

## Generalized linear models (GLM)


$S={S_{1},...,S_{K}}$ is a set of $k=1,...,K$ non-overlapping (areal) units.

$\boldsymbol{Y}=(y_{1},...,y_{K})$ are observations on $S$, then:


$$\begin{aligned}
Y & \sim f(y|\Omega)
g(\mu) &=\boldsymbol{x}^{T}\boldsymbol{\beta}+\boldsymbol{O}+\boldsymbol{\varepsilon}.
\end{aligned}$$

- $\mu = \text{E}(Y)$  is the expected value of $\boldsymbol{Y}$

- $f(\cdot)$ indicates an exponential function  

- $\Omega$ is the set of the parameters of the function $f(\cdot)$ 

## Generalized linear models (GLM) 2



- $g(\cdot) = f^{-1}(\cdot)$ is the linearizing link function 
 
- $\boldsymbol{x=}(x_{kv})\in\Re^{K\times V}$ is the matrix of covariates

- $\boldsymbol{\beta}$ are the $V$ covariate parameters with MVN prior, with mean $\mu_{\beta}$ and diagonal variance matrix $\Sigma_{\beta}$
  
- $\boldsymbol{O=}(o_{1},...,o_{K}\boldsymbol{)}$ are offsets, if any
 
- $\boldsymbol{\varepsilon}=(\varepsilon_{1},...,\varepsilon_{K})$ are residual errors, if any




##  Generalized linear models (GLM) ...
 
For each distributional family:
 
$Y\sim\text{Normal}(\mu,\sigma^{2})$ 

  - $\mu=\boldsymbol{x}^{T}\boldsymbol{\beta}+\boldsymbol{O}+\boldsymbol{\varepsilon}$

$Y\sim\text{Poisson}(\mu)$, and

  - $\text{ln}(\mu)=\boldsymbol{x}^{T}\boldsymbol{\beta}+\boldsymbol{O}+\boldsymbol{\varepsilon}$.

$Y\sim\text{Binomial}(\eta,\theta)$ 

  - $\text{ln}(\theta/(1-\theta))=\boldsymbol{x}^{T}\boldsymbol{\beta}+\boldsymbol{O}+\boldsymbol{\varepsilon}$,
  - $\eta$ is the vector of number of trials 
  - $\theta$ the vector of probabilities of success in each trial 
 

##  Generalized linear models (GLM) ... 2
 

**Stratified random sampling** assumes: 

$$\begin{aligned}
\varepsilon_{s} & \sim\text{N}(0,{}^{\varepsilon}\!\sigma_{s}^{2})
\end{aligned}$$

i.e., IID errors in space (time is usually ignored)

  
## Autocorrelated spatial errors 

 \small

**Kriging** extends $\varepsilon$ to spatial contraints ("variograms"). 

$$\begin{aligned}
Y_{t} & \sim\text{MVN}(\boldsymbol{\mu}_{t},\boldsymbol{\Sigma}_{t})\\
g(\mu_{t}) & =\boldsymbol{x_{t}}^{T}\boldsymbol{\beta_{t}}+\boldsymbol{\omega}_{t}+\boldsymbol{\varepsilon}_{t}\\
\varepsilon_{t} & \sim N(0,{}^{\varepsilon}\!\sigma_{t}^{2})\\
\omega_{t} & \sim\text{GP}(\boldsymbol{0},C(s_{t},s_{t}';^{\omega}\!\theta_{t}))\\
\boldsymbol{\Sigma_{t}} & =\left[C(\text{s}_{it},\text{s}_{jt};^{\omega}\!\theta_{t})\right]_{i,j=1}^{K}+{}^{\varepsilon}\!\sigma_{t}^{2}I_{S} \\
C(h)_{\text{Mat\'{e}rn}} & = ;^{\omega}\sigma^{2}\frac{1}{2^{\nu-1}\Gamma(\nu)}(\sqrt{2\nu}h/\phi)^{\nu}\ K_{\nu}(\sqrt{2\nu}h/\phi).
\end{aligned}$$
 
\footnotesize
-   ignores temporal structure, focus upon spatial-only process

-   assumes first and second order stationarity (incorrect)

-   problems when abundance declines and spatial autocorrelation changes
    across time or becomes unstable or impossible to estimable

\normalsize


<!--
## Spatial autocorrelation
 
Spatial autocorrelation function describes the proportion
of the spatial variance $C(h=0)=;^{\omega}\sigma^{2}$ decrease as distance
increases $h$.

Spatial covariance function $C(h)$ scaled by the total variance $C(0)$:

$\rho(h)=C(h)/C(0)$. 

The spatial covariance function:

$C(h)=C(s,s';\theta)$ 

expresses the tendency of observations closer together to be more similar to each other than those further away (Tobler's First law). 

Here, the distance, 

$h=\parallel s-s'\parallel$,

where $\parallel\cdot\parallel$ indicates a spatial norm which in $d=2$ spatial dimensions is simply the Euclidean distance:

$h=(\Delta\text{northing}^{2}+\Delta\text{easting}^{2})^{1/2}$. 




## Spatial autocorrelation 1

Commonly used forms include:

$$\begin{matrix}C(h)_{\text{Spherical}} & = & \left\{ \begin{matrix}^{\omega}\sigma^{2}(1-\frac{3}{2}h/\phi+\frac{1}{2}(h/\phi)^{3}); & 0<h<=\phi\\
0; & h>\phi,
\end{matrix}\right.\\
C(h)_{\text{Exponential}} & = & ^{\omega}\sigma^{2}e^{-h/\phi},\\
C(h)_{\text{Gaussian}} & = & ^{\omega}\sigma^{2}e^{-(h/\phi)^{2}},\\
C(h)_{\text{Powered exponential}} & = & ^{\omega}\sigma^{2}e^{-|h/\phi|^{p}},\\
C(h)_{\text{Mat\'{e}rn}} & = & ^{\omega}\sigma^{2}\frac{1}{2^{\nu-1}\Gamma(\nu)}(\sqrt{2\nu}h/\phi)^{\nu}\ K_{\nu}(\sqrt{2\nu}h/\phi).
\end{matrix}$$


## Spatial autocorrelation 2

At zero distance,
$C(0)=\text{Cov}(Y_{s},Y_{s})=\text{Var}(Y_{s})={}^{\varepsilon}\sigma^{2}+^{\omega}\sigma^{2}$
(*i.e.*, global spatial variance, also called the sill), where
$^{\varepsilon}\sigma$ is the non-spatial, unstructured error (also
called the nugget), $^{\omega}\sigma$ is the spatially structured error
(also called the partial sill), and
$^{\omega}\theta=\{\phi,\nu,p,\ldots\}$ are function-specific parameters
including $\phi$ the *range* parameter. $\Gamma(\cdot)$ is the Gamma
function and $K_{\nu}(\cdot)$ is the Bessel function of the second kind
with smoothness $\nu$. The Matérn covariance function is frequently used
in the more recent literature as the shape of this function is more
flexible and known to be connected to a Gaussian spatial random process. The *stmv* approach  permits a local
estimation of these autocorrelation parameters and the spatial scale.
 



-->

  
<!-- 


## STMV
 
```{r stmv-concept, out.width='40%', fig.show='hold', fig.align='center', fig.cap= 'Spatial distribution of data (blue dots) overlaid by a statistical grid. The $m$ nodes represent the centers of each local subdomain $S_{m}$ which extends to a distance (right-facing arrows; solid squares) that varies depending upon the underlying spatial variability of the data, defined as the distance at which the spatial autocorrelation drops to some small value (e.g., $\rho_{s}=0.1$). Data within this distance and parameters obtained from the local analysis are, under the assumption of second order stationarity' } 

include_graphics( file.path( framework_loc, "stmv_concept.png" ) ) 
```
 


## Random field-based approximations using SPDE (INLA)
 
$$\begin{matrix}Y & \sim & \text{N}(\mu_{st},^{\varphi}\!\sigma^{2})\\
g(\mu_{st}) & = & \boldsymbol{x}^{T}\boldsymbol{\beta}_{st}+\varphi_{st},\\
\varphi_{st} & \sim & \text{N}(0,^{\varphi}\!\sigma^{2}).
\end{matrix}$$

$$\begin{matrix}\varphi_{mt}^{*} & = & \omega_{mt}+\varepsilon_{mt},\\
\omega_{mt} & \sim & \text{GP}(0,C(\mathbf{s},\mathbf{s}';\mathbf{^{\boldsymbol{\omega}}\theta}_{mt}=\{\nu_{mt},\phi_{mt},^{\omega}\sigma_{mt}\})),\\
\varepsilon_{mt} & \sim & \text{Normal}(0,{}^{\varepsilon}\sigma_{m}^{2}).
\end{matrix}$$
 
SPDE representation of a random field is efficiently computed from:

$$\frac{\partial}{\partial t}(\kappa(s)^{2}-\Delta)^{\alpha/2}(\tau(\mathbf{s})x(\mathbf{s},\mathbf{t}))=\mathcal{W}(\mathbf{s},\mathbf{t}),\;(\mathbf{s},\mathbf{t})\in\Omega\times\mathbb{R}$$

-->


## autocorrelation

```{r autocorrelation, out.width='40%', fig.show='hold', fig.align='center', fig.cap= '' }

include_graphics( file.path( framework_loc, "autocorrelation.png" ) ) 
``` 
\tiny

Matérn autocorrelation function, $\rho(h)=C(h)/C(0)$, the covariance function $C(h)$ scaled by the total variance $C(0)$, for two values of $\nu$ (dark lines). As $\nu$ increases $(\nu=100)$, it approaches the Gaussian curve (upper dark curve on the left side) while at smaller values $(\nu=0.5)$ the curve is exponential (lower dark curve on the left side). This flexibility has made it a popular choice in geostatistics. The associated semivariograms (scaled to unit variance) $\gamma(h)$ are shown in light stippled lines. Spatial scale is defined heuristically as the distance $h$ at which the autocorrelation falls to a low value (dashed horizontal line). The semivariance (also called a semivariogram), $\gamma(h)$, is more commonly used in the geostatistics literature, and is simply the covariance function $C(h)$ reflected on the horizontal axis of the global variance $C(0)$ such that $\gamma(h)=C(0)-C(h)=\frac{1}{2}\ \text{Var}[Y_{s}-Y_{s}']=^{\omega}\sigma^{2}[1-\rho(h)]$.

\normalsize

 
## CARSTM: Conditional Autoregressive Space-Time Models (via INLA)

\small

$$\begin{aligned}
Y &\sim f(\mu,{}^{\varepsilon}\!\sigma^{2})\\
g(\mu) &=\boldsymbol{x}^{T}\boldsymbol{\beta}+\boldsymbol{O}+\phi+\boldsymbol{\varepsilon}\\
\varepsilon &\sim\text{N}(0,^{\varepsilon}\!\sigma^{2})\\
\phi &\sim N(\boldsymbol{0},Q^{-1})\\
Q &=[\tau(D-\alpha W)]^{-1}\\
p(\phi_{i}|\phi_{i\sim j},\tau) &=N(\alpha\sum_{i\sim j}^{K}w_{ij}\phi_{j},\tau^{-1})
\end{aligned}$$

-   Bayesian Hierarchical models of time \| space (INLA)

-   Spatial autocorrelation modelled as a local autocorrelation
    (CAR/BYM)

-   Temporal autocorrelation modelled as a local AR1 (ARIMA)

-   Ecosystem variability enters as covariates (smooths)

\normalsize


## CARSTM: Temperature Posterior Predictive Check

```{r temp_ppc, out.width='50%', fig.show='hold', fig.align='center', fig.cap= '' }
loc = file.path( data_root, "aegis", "temperature", "modelled", "default" )
include_graphics( file.path( loc, "posterior_predictive_check.png" ) ) 
``` 

## CARSTM: Species composition Posterior Predictive Check

```{r pca_ppc, out.width='40%', fig.show='hold', fig.align='center', fig.cap= '' }
loc = file.path( data_root, "aegis", "speciescomposition", "modelled", "default" )
fn1 = file.path( loc, "pca1_posterior_predictive_check.png" )
fn2 = file.path( loc, "pca2_posterior_predictive_check.png" )
fn3 = file.path( loc, "pca3_posterior_predictive_check.png" )
include_graphics( c(fn1, fn2 ) ) 
``` 
 
## CARSTM: Snow crab Posterior Predictive Check

```{r snowcrab_ppc, out.width='32%', fig.show='hold', fig.align='center', fig.cap= '' }
loc = file.path( data_root, "bio.snowcrab", "modelled", "default_fb" )
fn1 = file.path( loc, "totno_posterior_predictive_check.png" )
fn2 = file.path( loc, "meansize_posterior_predictive_check.png" )
fn3 = file.path( loc, "pa_posterior_predictive_check.png" )
include_graphics( c(fn1, fn2, fn3) ) 
``` 
 



## Fishery Model 

\small

$$\begin{aligned}
b_{t+1} & \sim N\left({{{b_{t}+r}b_{t}}{{({1-b_{t}})}-\mathit{Fishing}}_{t}K^{-1},\sigma_{p}}\right) \\
Y_{t} & \sim N\left({q^{-1}\mathit{Kb}_{t},\sigma_{o}}\right) \\
r & \sim N\left({1.0,0.1}\right) \\
K & \sim N\left({\kappa,0.1\cdot\kappa}\right) \\
q & \sim N\left({1.0,0.1}\right) \\
\end{aligned}$$

$\kappa={\lbrack{5.0,60,1.25}\rbrack}$

${b=B}K^{-1}$,  a real unobservable (latent) process

\normalsize  
 
## Prior and posterior comparisons: K 

```{r pppcK, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Prior (light blue) and posterior (purple) distributions of K.' }
loc = file.path(data_root, 'bio.snowcrab', 'fishery_model', year_assessment, 'logistic_discrete_historical' )
fn1 = file.path(loc, 'plot_prior_K_cfanorth.png')
fn2 = file.path(loc, 'plot_prior_K_cfasouth.png')
fn3 = file.path(loc, 'plot_prior_K_cfa4x.png')
# include_graphics( c(fn1)  ) 
include_graphics( c(fn1, fn2, fn3)  ) 
``` 



## Prior and posterior comparisons: r
 
```{r pppcr, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Prior (light blue) and posterior (purple) distributions of r.' }
loc = file.path(data_root, 'bio.snowcrab', 'fishery_model', year_assessment, 'logistic_discrete_historical' )
fn1 = file.path(loc, 'plot_prior_r_cfanorth.png')
fn2 = file.path(loc, 'plot_prior_r_cfasouth.png')
fn3 = file.path(loc, 'plot_prior_r_cfa4x.png')
include_graphics( c(fn1, fn2, fn3)  ) 
``` 

## Prior and posterior comparisons: q
 
```{r pppcq, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Prior (light blue) and posterior (purple) distributions of q.' }
loc = file.path(data_root, 'bio.snowcrab', 'fishery_model', year_assessment, 'logistic_discrete_historical' )
fn1 = file.path(loc, 'plot_prior_q1_cfanorth.png')
fn2 = file.path(loc, 'plot_prior_q1_cfasouth.png')
fn3 = file.path(loc, 'plot_prior_q1_cfa4x.png')
include_graphics( c(fn1, fn2, fn3)  )  
``` 

## Prior and posterior comparisons: observation error
 
```{r pppcbosd, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Prior (light blue) and posterior (purple) distributions of observation error.' }
loc = file.path(data_root, 'bio.snowcrab', 'fishery_model', year_assessment, 'logistic_discrete_historical' )
fn1 = file.path(loc, 'plot_prior_bosd_cfanorth.png')
fn2 = file.path(loc, 'plot_prior_bosd_cfasouth.png')
fn3 = file.path(loc, 'plot_prior_bosd_cfa4x.png') 
include_graphics( c( fn1, fn2, fn3)  )  
``` 

## Prior and posterior comparisons: process error
 
```{r pppcbpsd, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Prior (light blue) and posterior (purple) distributions of proess error.' }
loc = file.path(data_root, 'bio.snowcrab', 'fishery_model', year_assessment, 'logistic_discrete_historical' )
fn1 = file.path(loc, 'plot_prior_bpsd_cfanorth.png')
fn2 = file.path(loc, 'plot_prior_bpsd_cfasouth.png')
fn3 = file.path(loc, 'plot_prior_bpsd_cfa4x.png') 
include_graphics( c(fn1, fn2, fn3)  )   
``` 





## End
 
 

<!--  NOTES:

## References


Banerjee, S., Carlin, B. P., and Gelfand, A. E.. 2004. Hierarchical Modeling and Analysis for Spatial Data. Monographs on Statistics and Applied Probability. Chapman and Hall/CRC.

Besag, Julian. 1974. Spatial interaction and the statistical analysis of lattice systems. Journal of the Royal Statistical Society Series B (Methodological) 1974: 192-236.

Canada Gazette. 2022. [Regulations Amending the Fishery (General) Regulations. Part II, Volume 156, Number 8.](https://www.gazette.gc.ca/rp-pr/p2/2022/2022-04-13/html/sor-dors73-eng.html)

Canada Gazette. 2016. St. Anns Bank Marine Protected Area Regulations. Canada Gazette, Part I, Vol 150, Issue 51: 4143-4149.

Choi, J.S. 2023. A Framework for the assessment of Snow Crab (*Chioneocete opilio*) in Maritimes Region (NAFO Div 4VWX) . DFO Can. Sci. Advis. Sec. Res. Doc. 2023/nnn. v + xxx p.

Choi, J.S. 2022. Reconstructing the Decline of Atlantic Cod with the Help of Environmental Variability in the Scotian Shelf of Canada. bioRxiv. https://doi.org/10.1101/2022.05.05.490753.
 
Choi, Jae S., B. Cameron, K. Christie, A. Glass, and E. MacEachern. 2022. Temperature and Depth Dependence of the Spatial Distribution of Snow Crab. bioRxiv. https://doi.org/10.1101/2022.12.20.520893.


Choi, Jae S. 2023. A Multi-Stage, Delay Differential Model of Snow Crab Population Dynamics in the Scotian Shelf of Atlantic Canada. bioRxiv. https://doi.org/10.1101/2023.02.13.528296.
s
DFO. 2013. [Integrated Fisheries Management Plan for Eastern Nova Scotia and 4X Snow Crab (*Chionoecetes Opillio*.)](http://www.dfo-mpo.gc.ca/fm-gp/peches-fisheries/ifmp-gmp/snow-crab-neige/snow-crab-neiges2013-eng.htm)
 
Riebler, A., Sørbye, S.H., Simpson D., and Rue, H. 2016. An intuitive Bayesian spatial model for disease mapping that accounts for scaling. Statistical methods in medical research 25: 1145-1165.

Simpson, D., Rue, H., Riebler, A., Martins, T.G., and Sørbye, SH. 2017. Penalising Model Component Complexity: A Principled, Practical Approach to Constructing Priors. Statist. Sci. 32: 1-28.
 
-->  

 
