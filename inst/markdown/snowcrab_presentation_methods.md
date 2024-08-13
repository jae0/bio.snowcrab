---
title: "Snow Crab, Scotian Shelf, Canada (NAFO Div. 4VWX) in 2023"
subtitle: "Life history and Methodological considerations"
author: "Snow Crab Group"
# author: "Jae S. Choi"
# footnote: "jae.choi@dfo-mpo.gc.ca"
institute: "Bedford Institute of Oceanography, DFO Science"
# date: "`r format(Sys.time(), '%d %B, %Y')`"
output:
  beamer_presentation:
    theme: "metropolis"
    colortheme: "seagull"
    fonttheme: "professionalfonts"
    fig_caption: yes
    # latex_engine: pdflatex
    latex_engine: lualatex 
    keep_tex: true
classoption: 
  - aspectratio=169 #16:9 wide
  - t  # top align
header-includes: 
  - \usepackage{graphicx}
  - \usepackage[font={scriptsize}, labelfont={bf}]{caption}
  # - \usepackage{float}
  # - \usepackage{subfig}
  # - \newcommand{\btiny}{\begin{tiny}}
  # - \newcommand{\etiny}{\end{tiny}}
params:
  year.assessment: 2023
  media_loc: "media"
  debugging: FALSE 
--- 


<!-- Preamble


This is a Markdown document ... To create HTML or PDF, etc, run: 


  make quarto FN=snow_crab_presentation YR=2023 SOURCE=~/projects/bio.snowcrab/inst/markdown WK=~/bio.data/bio.snowcrab/assessments # {via Quarto}

  make rmarkdown FN=snow_crab_presentation YR=2023 SOURCE=~/projects/bio.snowcrab/inst/markdown WK=~/bio.data/bio.snowcrab/assessments {via Rmarkdown}

  make pdf FN=snow_crab_presentation  # {via pandoc}

Alter year and directories to reflect setup or copy Makefile and alter defaults to your needs.
  
  

Huge > huge > LARGE > Large > large > normalsize > small > footnotesize > scriptsize > tiny


::: columns 
:::: column 

::::
 
:::: column

::::
:::


-->



<!-- Set up R-environment -->

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

  media_loc = params$media_loc
  year.assessment = params$year.assessment
  year_previous = year.assessment - 1
  p = bio.snowcrab::load.environment( year.assessment=year.assessment )
  SCD = project.datadirectory("bio.snowcrab")
  
    

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
  require(janitor)
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
  


# Life history {.c}

## Photos {.c}
::: columns 
:::: column 
```{r photos, echo=FALSE, out.width='90%', fig.align='center', fig.show='hold',  fig.cap = 'Snow Crab pelagic Zoea.' }
fn1=file.path( media_loc, "snowcrab_zoea.png" )
knitr::include_graphics( c(fn1 ) ) 
# \@ref(fig:photos)  
```
::::
:::: column
```{r photos2, echo=FALSE, out.width='90%', fig.align='center', fig.show='hold',  fig.cap = 'Snow Crab male and mating pair. Note sexual dimorphism.' }
fn2=file.path( media_loc, "snowcrab_male.png" )
fn3=file.path( media_loc, "snowcrab_male_and_female.png" )
knitr::include_graphics( c( fn2, fn3) ) 
# \@ref(fig:photos)  
```
::::
:::


## Life history stages{.c}
 
```{r lifehistory, echo=FALSE, out.width='95%', fig.align='center', fig.cap = 'Life history patterns of snow crab and approximate timing of the main life history stages of snow crab and size (carapace width; CW mm) and instar (Roman numerals). Size and timings are specific to the area of study and vary with environmental conditions, food availability and genetic variability.' }
loc = file.path( Sys.getenv("HOME"), "projects", "dynamical_model", "snowcrab", "media" )
fn1=file.path( loc, "life_history.png" )
knitr::include_graphics( fn1 ) 
# \@ref(fig:lifehistory)  
```

## Growth stanzas {.c}
 
```{r lifehistory_male, echo=FALSE, out.width='55%', fig.align='center', fig.cap = 'The growth stanzas of the male component and decision paths to maturity and terminal moult. Black ellipses indicate terminally molted animals.' }
loc = file.path( Sys.getenv("HOME"), "projects", "dynamical_model", "snowcrab", "media" )
fn1=file.path( loc, "life_history_male.png" )
knitr::include_graphics( fn1 ) 
# \@ref(fig:lifehistory_male)  
```

## Growth modes{.c}

```{r growth_modes1, echo=FALSE, out.width='100%', fig.align='center', fig.cap = 'Modal analysis.' }
loc = file.path( Sys.getenv("HOME"), "bio.data", "bio.snowcrab", "output" )
fn1=file.path( loc, "size_structure", "growth_summary1.png" )
knitr::include_graphics( c(fn1) ) 
```
 
## Growth modes 2 {.c}

```{r growth_modes2, echo=FALSE, out.width='30%', fig.align='center', fig.cap = 'Modal analysis.' }
loc = file.path( Sys.getenv("HOME"), "bio.data", "bio.snowcrab", "output" )
fn1=file.path( loc, "size_structure", "growth_summary2.png" )
knitr::include_graphics( c(fn1) ) 
```

 
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
  
::: columns 
:::: column 
```{r survey-locations-map, out.width='60%', fig.show='hold', fig.align='center', fig.cap= 'Snow Crab survey locations. No survey in 2020 (Covid-19) and incomplete 2022 (mechanical issues).' }
loc = file.path( SCD, "output", "maps", "survey.locations" )
yrsplot = setdiff( year.assessment + c(0:-9), 2020)
fn6 = file.path( loc, paste( "survey.locations", yrsplot[6], "png", sep=".") )
fn5 = file.path( loc, paste( "survey.locations", yrsplot[5], "png", sep=".") )
fn4 = file.path( loc, paste( "survey.locations", yrsplot[4], "png", sep=".") )
fn3 = file.path( loc, paste( "survey.locations", yrsplot[3], "png", sep=".") )
fn2 = file.path( loc, paste( "survey.locations", yrsplot[2], "png", sep=".") )
fn1 = file.path( loc, paste( "survey.locations", yrsplot[1], "png", sep=".") )
knitr::include_graphics( c( fn2, fn1) )
# \@ref(fig:survey-locations-map)  
```
::::
:::: column

```{r surveydomain, out.width='100%', fig.show='hold', fig.align='center', fig.cap= 'Snow Crab survey prediction grid.' }
fn1 = file.path( media_loc, "carstm_prediction_domain.png" )
knitr::include_graphics( c( fn1) )
```
::::
:::


## Clustering  {.c}
 
```{r aggregation, echo=FALSE, out.width='40%', fig.align='center', fig.show='hold',  fig.cap = 'Australian spider crab \\emph{Leptomithrax gaimardii} aggregation for moulting and migration and Alaska red king crab \\emph{Paralithodes camtschaticus} aggregation in Alaska for egg release, migrations.' }
fn1=file.path( media_loc, "australian_leptomithrax_gaimardii.png" )
fn2=file.path( media_loc, "kingcrab_aggregation.png" ) 
knitr::include_graphics( c(fn1, fn2 ) ) 
# \@ref(fig:aggregation)  
```

- *Other* crab species show "Mounding" for protection from predation (larval, moulting and females)

- Narrow habitat preferences force them to move and cluster when environment is poor
 

## Clustering2  {.c}
 
```{r clustering, echo=FALSE, out.width='60%', fig.align='center', fig.show='hold',  fig.cap = 'High density locations of Snow Crab, approximately 1 per square meter.'}
fn = file.path(p$project.outputdir, "maps", "map_highdensity_locations.png" )
knitr::include_graphics( fn ) 
# \@ref(fig:aggregation)  
if (0) {
    # high density locations directly from databases
    M = snowcrab.db( DS="set.complete", p=p ) 
    setDT(M)
    i = which(M$totno.all > 2.5*10^5)
    H = M[i, .( plon, plat, towquality, dist, distance, surfacearea, vessel, yr, z, julian, no.male.all, no.female.all, cw.mean, totno.all, totno.male.imm, totno.male.mat, totno.female.imm, totno.female.mat, totno.female.primiparous, totno.female.multiparous, totno.female.berried)]
    H$log10density = log10(H$totno.all)
    library(ggplot2)
    cst = coastline_db( p=p, project_to=st_crs(pg) ) 
    isodepths = c(100, 200, 300)
    isob = isobath_db( DS="isobath", depths=isodepths, project_to=st_crs(pg))
    isob$level = as.factor( isob$level)
    plt = ggplot() +
      geom_sf( data=cst, show.legend=FALSE ) +
      geom_sf( data=isob, aes( alpha=0.1, fill=level), lwd=0.1, show.legend=FALSE) +
    geom_point(data=H, aes(x=plon, y=plat, colour=log10density), size=5) +
      coord_sf(xlim = c(270, 940 ), ylim = c(4780, 5200 )) +
    theme(legend.position="inside", legend.position.inside=c(0.08, 0.8)) 
    png(filename=fn, width=1000,height=600, res=144)
      (plt)
    dev.off()
}

```

## Sampling bias: depth
```{r bias-depth, out.width='40%', fig.show='hold', fig.align='center', fig.cap= 'Comparison of depths in snow crab domain ("predictions") and survey locations ("observations"). Bias: survey trawls preferentially sample deeper locations. Red line is overall average.' }
loc = file.path( homedir, "projects", "snowcrabframework" )
knitr::include_graphics( file.path( loc, "bias_depth.png" ) )
```
 
## Sampling bias: substrate grain size

```{r bias-substrate, out.width='40%', fig.show='hold', fig.align='center', fig.cap= 'Comparison of substrate grain size (log; mm) in snow crab domain ("predictions") and survey locations ("observations"). No bias observed. Red line is overall average.' }
loc = file.path( homedir, "projects", "snowcrabframework" )
knitr::include_graphics( file.path( loc, "bias_substrate.png" ) ) 
``` 

## Sampling bias: bottom temperature

```{r bias-temperature, out.width='40%', fig.show='hold', fig.align='center', fig.cap= 'Comparison of bottom temperature (degrees Celcius) in snow crab domain ("predictions") and survey locations ("observations"). Bias: surveys preferentially sample cold water bottoms. Red line is overall average.' }
loc = file.path( homedir, "projects", "snowcrabframework" )
knitr::include_graphics( file.path( loc, "bias_temp.png" ) ) 
``` 
 
## Sampling bias: species composition 1 (temperature related)

```{r bias-pc1, out.width='40%', fig.show='hold', fig.align='center', fig.cap= 'Comparison of species composition gradient (PC1) in snow crab domain ("predictions") and survey locations ("observations"). Bias: surveys preferentially sample cold water species. Red line is overall average.' }
loc = file.path( homedir, "projects", "snowcrabframework" )
knitr::include_graphics( file.path( loc, "bias_pca1.png" ) ) 
``` 
  

## Sampling bias: species composition 2 (depth related)

```{r bias-pc2, out.width='40%', fig.show='hold', fig.align='center', fig.cap= 'Comparison of species composition gradient (PC2) in snow crab domain ("predictions") and survey locations ("observations"). Bias: surveys preferentially sample deeper species. Red line is overall average.' }
loc = file.path( homedir, "projects", "snowcrabframework" )
knitr::include_graphics( file.path( loc, "bias_pca2.png" ) ) 
``` 

## Sampling bias: summary

Indications of sampling bias relative to spatial domain of snow crab:

- deeper and more favourable locations are being sampled
- statistically overlapping but biologically important
 

## Generalized linear models (GLM)

\small

::: columns 
:::: column 

$S={S_{1},...,S_{K}}$ is a set of $k=1,...,K$ non-overlapping (areal) units.

$\boldsymbol{Y}=(y_{1},...,y_{K})$ are observations on $S$, then:


$$\begin{aligned}
Y & \sim f(y|\Omega)
g(\mu) &=\boldsymbol{x}^{T}\boldsymbol{\beta}+\boldsymbol{O}+\boldsymbol{\varepsilon}.
\end{aligned}$$

- $\mu = \text{E}(Y)$  is the expected value of $\boldsymbol{Y}$

- $f(\cdot)$ indicates an exponential function  

- $\Omega$ is the set of the parameters of the function $f(\cdot)$ 

::::
 
:::: column

- $g(\cdot) = f^{-1}(\cdot)$ is the linearizing link function 
 
- $\boldsymbol{x=}(x_{kv})\in\Re^{K\times V}$ is the matrix of covariates

- $\boldsymbol{\beta}$ are the $V$ covariate parameters with MVN prior, with mean $\mu_{\beta}$ and diagonal variance matrix $\Sigma_{\beta}$
  
- $\boldsymbol{O=}(o_{1},...,o_{K}\boldsymbol{)}$ are offsets, if any
 
- $\boldsymbol{\varepsilon}=(\varepsilon_{1},...,\varepsilon_{K})$ are residual errors, if any

::::
:::

\normalsize




##  Generalized linear models (GLM) ...


\small
::: columns 
:::: column 
For each distributional family:
\vspace{2mm}

$Y\sim\text{Normal}(\mu,\sigma^{2})$ 

  - $\mu=\boldsymbol{x}^{T}\boldsymbol{\beta}+\boldsymbol{O}+\boldsymbol{\varepsilon}$

$Y\sim\text{Poisson}(\mu)$, and

  - $\text{ln}(\mu)=\boldsymbol{x}^{T}\boldsymbol{\beta}+\boldsymbol{O}+\boldsymbol{\varepsilon}$.

$Y\sim\text{Binomial}(\eta,\theta)$ 

  - $\text{ln}(\theta/(1-\theta))=\boldsymbol{x}^{T}\boldsymbol{\beta}+\boldsymbol{O}+\boldsymbol{\varepsilon}$,
  - $\eta$ is the vector of number of trials 
  - $\theta$ the vector of probabilities of success in each trial 
::::
 
:::: column

**Stratified random sampling** assumes: 

$$\begin{aligned}
\varepsilon_{s} & \sim\text{N}(0,{}^{\varepsilon}\!\sigma_{s}^{2})
\end{aligned}$$

i.e., IID errors in space (time is usually ignored)

::::
:::


\normalsize

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
loc = file.path( homedir, "projects", "snowcrabframework" )
knitr::include_graphics( file.path( loc, "stmv_concept.png" ) ) 
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
loc = file.path( homedir, "projects", "snowcrabframework" )
knitr::include_graphics( file.path( loc, "autocorrelation.png" ) ) 
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
knitr::include_graphics( file.path( loc, "posterior_predictive_check.png" ) ) 
``` 

## CARSTM: Species composition Posterior Predictive Check

```{r pca_ppc, out.width='40%', fig.show='hold', fig.align='center', fig.cap= '' }
loc = file.path( data_root, "aegis", "speciescomposition", "modelled", "default" )
fn1 = file.path( loc, "pca1_posterior_predictive_check.png" )
fn2 = file.path( loc, "pca2_posterior_predictive_check.png" )
fn3 = file.path( loc, "pca3_posterior_predictive_check.png" )
knitr::include_graphics( c(fn1, fn2 ) ) 
``` 
 
## CARSTM: Snow crab Posterior Predictive Check

```{r snowcrab_ppc, out.width='32%', fig.show='hold', fig.align='center', fig.cap= '' }
loc = file.path( data_root, "bio.snowcrab", "modelled", "default_fb" )
fn1 = file.path( loc, "totno_posterior_predictive_check.png" )
fn2 = file.path( loc, "meansize_posterior_predictive_check.png" )
fn3 = file.path( loc, "pa_posterior_predictive_check.png" )
knitr::include_graphics( c(fn1, fn2, fn3) ) 
``` 
 


# Model Diagnostics

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
loc = file.path(data_root, 'bio.snowcrab', 'fishery_model', '2023', 'logistic_discrete_historical' )
fn1 = file.path(loc, 'plot_prior_K_cfanorth.png')
fn2 = file.path(loc, 'plot_prior_K_cfasouth.png')
fn3 = file.path(loc, 'plot_prior_K_cfa4x.png')
# knitr::include_graphics( c(fn1)  ) 
knitr::include_graphics( c(fn1, fn2, fn3)  ) 
``` 



## Prior and posterior comparisons: r
 
```{r pppcr, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Prior (light blue) and posterior (purple) distributions of r.' }
loc = file.path(data_root, 'bio.snowcrab', 'fishery_model', '2023', 'logistic_discrete_historical' )
fn1 = file.path(loc, 'plot_prior_r_cfanorth.png')
fn2 = file.path(loc, 'plot_prior_r_cfasouth.png')
fn3 = file.path(loc, 'plot_prior_r_cfa4x.png')
knitr::include_graphics( c(fn1, fn2, fn3)  ) 
``` 

## Prior and posterior comparisons: q
 
```{r pppcq, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Prior (light blue) and posterior (purple) distributions of q.' }
loc = file.path(data_root, 'bio.snowcrab', 'fishery_model', '2023', 'logistic_discrete_historical' )
fn1 = file.path(loc, 'plot_prior_q1_cfanorth.png')
fn2 = file.path(loc, 'plot_prior_q1_cfasouth.png')
fn3 = file.path(loc, 'plot_prior_q1_cfa4x.png')
knitr::include_graphics( c(fn1, fn2, fn3)  )  
``` 

## Prior and posterior comparisons: observaton error
 
```{r pppcbosd, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Prior (light blue) and posterior (purple) distributions of observation error.' }
loc = file.path(data_root, 'bio.snowcrab', 'fishery_model', '2023', 'logistic_discrete_historical' )
fn1 = file.path(loc, 'plot_prior_bosd_cfanorth.png')
fn2 = file.path(loc, 'plot_prior_bosd_cfasouth.png')
fn3 = file.path(loc, 'plot_prior_bosd_cfa4x.png') 
knitr::include_graphics( c( fn1, fn2, fn3)  )  
``` 
## Prior and posterior comparisons: process error
 
```{r pppcbpsd, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Prior (light blue) and posterior (purple) distributions of proess error.' }
loc = file.path(data_root, 'bio.snowcrab', 'fishery_model', '2023', 'logistic_discrete_historical' )
fn1 = file.path(loc, 'plot_prior_bpsd_cfanorth.png')
fn2 = file.path(loc, 'plot_prior_bpsd_cfasouth.png')
fn3 = file.path(loc, 'plot_prior_bpsd_cfa4x.png') 
knitr::include_graphics( c(fn1, fn2, fn3)  )   
``` 





# End
 


## Fishery Model Biomass (pre-fishery){.c}

Results using less informative priors

\begin{tiny}
```{r logisticPredictions, out.width='32%', echo=FALSE, fig.show='hold', fig.align='center', fig.cap = 'Model 1 fishable, posterior mean modelled biomass (pre-fishery; kt) are shown in dark orange for N-ENS, S-ENS and 4X (left, middle and right). Light orange are posterior samples of modelled biomass (pre-fishery; kt) to illustrate the variability of the predictions. The biomass index (post-fishery, except prior to 2004) after model adjustment by the model catchability coefficient is in gray.' } 
# loc = file.path( SCD, 'fishery_model', year.assessment, 'logistic_discrete_historical' )
loc = file.path( homedir, "projects", "dynamical_model", "snowcrab", "outputs_alt",  year.assessment, "logistic_discrete_historical" )
fn1 = file.path( loc, 'plot_predictions_cfanorth.png' ) 
fn2 = file.path( loc, 'plot_predictions_cfasouth.png' ) 
fn3 = file.path( loc, 'plot_predictions_cfa4x.png' ) 
include_graphics(c(fn1, fn2, fn3) )
# \@ref(fig:logisticPredictions)
```
\end{tiny}

 

## Fishing Mortality {.c}

```{r logisticFishingMortality, out.width='32%', echo=FALSE,  fig.show='hold', fig.align='center', fig.cap = 'Time-series of modelled instantaneous fishing mortality from Model 1, for N-ENS (left), S-ENS (middle), and 4X (right). Samples of the posterior densities are presented, with the darkest line being the mean.' }
  # odir = file.path( fishery_model_results, year.assessment, "logistic_discrete_historical" )
  odir = file.path( homedir, "projects", "dynamical_model", "snowcrab", "outputs_alt",  year.assessment, "logistic_discrete_historical" )
  fn1 = file.path( odir, "plot_fishing_mortality_cfanorth.png" ) 
  fn2 = file.path( odir, "plot_fishing_mortality_cfasouth.png" ) 
  fn3 = file.path( odir, "plot_fishing_mortality_cfa4x.png" ) 
include_graphics(c(fn1, fn2, fn3) )
# \@ref(fig:logisticFishingMortality)
```


<!-- 
## Fishery Model Summary  {.c}

|   | N-ENS | S-ENS | 4X |  
|----- | ----- | ----- | ----- |
| |  |  |  |
|q       | `r round(q_north, 3)` (`r round(q_north_sd, 3)`) | `r round(q_south, 3)` (`r round(q_south_sd, 3)`) | `r round(q_4x, 3)` (`r round(q_4x_sd, 3)`) |
|r       | `r round(r_north, 3)` (`r round(r_north_sd, 3)`) | `r round(r_south, 3)` (`r round(r_south_sd, 3)`) | `r round(r_4x, 3)` (`r round(r_4x_sd, 3)`) |
|K       | `r round(K_north, 2)` (`r round(K_north_sd, 2)`) | `r round(K_south, 2)` (`r round(K_south_sd, 2)`) | `r round(K_4x, 2)` (`r round(K_4x_sd, 2)`) |
|Prefishery Biomass   | `r round(B_north[t0], 2)` (`r round(B_north_sd[t0], 2)`) | `r round(B_south[t0], 2)`  (`r round(B_south_sd[t0], 2)`) | `r round(B_4x[t0], 2)`  (`r round(B_4x_sd[t0], 2)`)  |
|Fishing Mortality    | `r round(FM_north[t0], 3)` (`r round(FM_north_sd[t0], 3)`) | `r round(FM_south[t0], 3)` (`r round(FM_south_sd[t0], 3)`) | `r round(FM_4x[t0], 3)` (`r round(FM_4x_sd[t0], 3)`) |


\tiny
Note: Values in parentheses are Posterior standard deviations.
\normalsize

-->


## Reference Points ... {.c}

```{r logistic-hcr, out.width='29%', echo=FALSE, fig.show='hold', fig.align='center', fig.cap = 'Reference Points (fishing mortality and modelled biomass) from the Fishery Model, for N-ENS (left), S-ENS (middle), and 4X (right). The large yellow dot indicates most recent year and the 95\\% CI. Not: the model does not account for illegal and unreported landings, and interspecific interactions.' }
  # odir = file.path( fishery_model_results, year.assessment, "logistic_discrete_historical" )
  odir = file.path( homedir, "projects", "dynamical_model", "snowcrab", "outputs_alt",  year.assessment, "logistic_discrete_historical" )
  fn1 = file.path( odir, 'plot_hcr_cfanorth.png' ) 
  fn2 = file.path( odir, 'plot_hcr_cfasouth.png' ) 
  fn3 = file.path( odir, 'plot_hcr_cfa4x.png' ) 
  include_graphics(c(fn1, fn2, fn3) )
#  \@ref(fig:logistic-hcr)
```
 

## Prior and posterior comparisons: K 

```{r pppcKtest, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Prior (light blue) and posterior (purple) distributions of K.' }
# loc = file.path(data_root, 'bio.snowcrab', 'fishery_model', '2023', 'logistic_discrete_historical' )
loc = file.path( homedir, "projects", "dynamical_model", "snowcrab", "outputs_alt",  year.assessment, "logistic_discrete_historical" )
fn1 = file.path(loc, 'plot_prior_K_cfanorth.png')
fn2 = file.path(loc, 'plot_prior_K_cfasouth.png')
fn3 = file.path(loc, 'plot_prior_K_cfa4x.png')
# knitr::include_graphics( c(fn1)  ) 
knitr::include_graphics( c(fn1, fn2, fn3)  ) 
``` 



## Prior and posterior comparisons: r
 
```{r pppcrtest, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Prior (light blue) and posterior (purple) distributions of r.' }
# loc = file.path(data_root, 'bio.snowcrab', 'fishery_model', '2023', 'logistic_discrete_historical' )
loc = file.path( homedir, "projects", "dynamical_model", "snowcrab", "outputs_alt",  year.assessment, "logistic_discrete_historical" )
fn1 = file.path(loc, 'plot_prior_r_cfanorth.png')
fn2 = file.path(loc, 'plot_prior_r_cfasouth.png')
fn3 = file.path(loc, 'plot_prior_r_cfa4x.png')
knitr::include_graphics( c(fn1, fn2, fn3)  ) 
``` 

## Prior and posterior comparisons: q
 
```{r pppcqtest, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Prior (light blue) and posterior (purple) distributions of q.' }
# loc = file.path(data_root, 'bio.snowcrab', 'fishery_model', '2023', 'logistic_discrete_historical' )
loc = file.path( homedir, "projects", "dynamical_model", "snowcrab", "outputs_alt",  year.assessment, "logistic_discrete_historical" )
fn1 = file.path(loc, 'plot_prior_q1_cfanorth.png')
fn2 = file.path(loc, 'plot_prior_q1_cfasouth.png')
fn3 = file.path(loc, 'plot_prior_q1_cfa4x.png')
knitr::include_graphics( c(fn1, fn2, fn3)  )  
``` 

## Prior and posterior comparisons: observaton error
 
```{r pppcbosdtest, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Prior (light blue) and posterior (purple) distributions of observation error.' }
# loc = file.path(data_root, 'bio.snowcrab', 'fishery_model', '2023', 'logistic_discrete_historical' )
loc = file.path( homedir, "projects", "dynamical_model", "snowcrab", "outputs_alt",  year.assessment, "logistic_discrete_historical" )
fn1 = file.path(loc, 'plot_prior_bosd_cfanorth.png')
fn2 = file.path(loc, 'plot_prior_bosd_cfasouth.png')
fn3 = file.path(loc, 'plot_prior_bosd_cfa4x.png') 
knitr::include_graphics( c( fn1, fn2, fn3)  )  
``` 
## Prior and posterior comparisons: process error
 
```{r pppcbpsdtest, out.width='32%', fig.show='hold', fig.align='center', fig.cap= 'Prior (light blue) and posterior (purple) distributions of proess error.' }
# loc = file.path(data_root, 'bio.snowcrab', 'fishery_model', '2023', 'logistic_discrete_historical' )
loc = file.path( homedir, "projects", "dynamical_model", "snowcrab", "outputs_alt",  year.assessment, "logistic_discrete_historical" )
fn1 = file.path(loc, 'plot_prior_bpsd_cfanorth.png')
fn2 = file.path(loc, 'plot_prior_bpsd_cfasouth.png')
fn3 = file.path(loc, 'plot_prior_bpsd_cfa4x.png') 
knitr::include_graphics( c(fn1, fn2, fn3)  )   
``` 



<!--  NOTES:

## References

\begin{tiny}


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


\end{tiny}

-->  

 
