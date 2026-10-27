# 05_snowcrab_fishery_model.md

## Purpose

Prepare aggregated time-series data for inference using a Bayesian fishery model. Currently, the default model is the discrete logistic formulation. It is labeled, "logistic_historical" in the following. There are other variations possible. See the calling code.

---

## Prepare data for discrete logistic models and run it for a chosen area

This step requires:the following to have been completed:

- [bio.snowcrab/inst/markdown/03_biomass_index_carstm.md](https://github.com/jae0/bio.snowcrab/inst/markdown/03_biomass_index_carstm.md)

```R

require(bio.snowcrab)

year.assessment = 2025  # change this as appropriate

p = bio.snowcrab::load.environment( year.assessment=year.assessment )  # set up initial settings

library(coda) # For handling mcmc.list objects
library(nimble)
 
model_variation = "logistic_discrete_historical"   
outputs_dir = file.path( p$datadir, "fishery_model", string(year_assessment), model_variation )
  
mkpath( outputs_dir )

year_assessment = 2025
aulab = "cfanorth"
 

# choose:
is_test = FALSE
# is_test = TRUE

niter = if (is_test) 1000 else 50000
nburnin = if (is_test) 500 else 10000  
nchains = if (is_test) 1 else 4
nthin = if (is_test) 1 else 50
truncate_max_years = if (is_test) 1 else 5

do_retrospective = TRUE
do_plots = TRUE
outformat = "png"

# 1. Prepare data and parameter setup
setup = prepare_fishery_model_data(
  year_assessment = year_assessment,
  aulab = aulab,
  outputs_dir = outputs_dir,
  model_variation = "logistic_discrete_historical"
)

# 2. Build model inputs using factored model library
model_inputs = prepare_nimble_model_inputs(setup, model_type = "gaussian")

nimble_model = nimbleModel(
  code = model_inputs$code,
  constants = model_inputs$constants,
  data = model_inputs$data,
  dimensions = model_inputs$dimensions,
  inits = model_inputs$inits,
  buildDerivs = model_inputs$buildDerivs,
  check = TRUE
)

C_nimble_model = compileNimble(nimble_model)

# 3. Configure MCMC with WAIC enabled
mcmc_conf = configureMCMC(
  nimble_model,
  monitors = model_inputs$monitors,
  enableWAIC = TRUE
)
mcmc_builder = buildMCMC(mcmc_conf)
C_mcmc_builder = compileNimble(mcmc_builder, project = nimble_model)

# 4. Run MCMC sampling
cat("\nRunning MCMC sampling...\n")
mcmc_out = runMCMC(
  C_mcmc_builder,
  niter = niter,
  nburnin = nburnin,
  nchains = nchains,
  thin = nthin,
  progressBar = TRUE,
  samplesAsCodaMCMC = TRUE,
  WAIC = TRUE
)

samples = if (is.list(mcmc_out) && "samples" %in% names(mcmc_out)) {
  mcmc_out$samples
} else {
  mcmc_out
}

waic_res = NULL
if (is.list(mcmc_out) && "WAIC" %in% names(mcmc_out)) {
  waic_res = mcmc_out$WAIC
  cat(sprintf("\nWAIC: %.2f (pWAIC: %.2f)\n", waic_res$WAIC, waic_res$pWAIC))
}

# 5. Prior sampling and posterior predictive check

n_prior = nrow(samples)
n_ppc = n_prior

prefix = if (is_test) paste0(aulab, "_gaussian_test_") else paste0(aulab, "_gaussian_")

prior = sample_priors_nimble(setup, model_type = "gaussian", n_samples = n_prior)

ppc = posterior_predictive_check(
  samples = samples,
  setup = setup,
  model_type = "gaussian",
  n_sims = n_ppc,
  output_dir = setup$outputs_dir,
  prefix = prefix,
  outformat = outformat
)

ppc$fit    # Pearson r, r^2, RMSE, bias, coverage
ppc$plot   # observed vs predicted

# 6. Post-processing, figures, and zone classification
output_prefix = if (is_test) "results_nimble_gaussian_test" else "results_nimble_gaussian"
post_process_and_save_fishery_model(
  samples = as.mcmc.list(samples),
  setup = setup,
  prior = prior,
  do_plots = do_plots,
  do_hcr_probabilities = TRUE,
  outformat = outformat,
  output_prefix = output_prefix
)

# 7. Retrospective analysis (optional)
retrospective_results = NULL
if (do_retrospective && truncate_max_years > 0) {
  cat(sprintf("\nStarting retrospective analysis (truncating up to %d years)...\n",
              truncate_max_years))
  retro_iter = if (is_test) 500 else 2000
  retro_burnin = if (is_test) 200 else 1000

  retrospective_mcmc_settings = list(
    niter = retro_iter,
    nburnin = retro_burnin,
    nchains = 2,
    progressBar = TRUE,
    samplesAsCodaMCMC = TRUE
  )

  retrospective_results = do_retrospective_analysis(
    nimble_code = model_inputs$code,
    PM_base = setup$PM,
    nimble_constants_base = model_inputs$constants,
    model_data_base = model_inputs$data,
    initial_values_base = model_inputs$inits,
    monitor_params = c("K", "r", "q1"),
    truncate_max_years = truncate_max_years,
    rho_target_param = "K",
    mcmc_settings = retrospective_mcmc_settings,
    buildDerivs = model_inputs$buildDerivs
  )

  plot_retrospective_analysis_results(
    retrospective_results,
    setup,
    "K",
    outformat = outformat
  )
}


# or, more compactly: 

# run_biomass_dynamics_gaussian( year_assessment = year_assessment, aulab = aulab, is_test = is_test )

# run_biomass_dynamics_lognormal( year_assessment = year_assessment, aulab = aulab, is_test = is_test_mode )

# run_delay_difference( year_assessment = year_assessment, aulab = aulab, is_test = is_test_mode ) 


```

## Data

### Survey Timing & Catch Penalty Logic

The survey timing on the Scotian Shelf changed in 2004:
- **Pre-2004**: Spring survey (prior to fishery removals). Removals occur after the survey, so survey index reflects unexploited spring biomass ($C_{pen} = 0$).
- **2004**: Transition year. Penalty is averaged: $C_{pen} = (C_t + C_{t-1}) / 2$.
- **Post-2004**: Fall survey (after fishery removals). Removals have already taken place ($C_{pen} = C_t$).
- **CFA 4X**: Consistently fall survey throughout the entire time series ($yeartransition = 0$, $C_{pen} = C_t$).

The function `compute_removed_penalty(PM)` computes this vector safely across all survey observation years without indexing or vector recycling artifacts.


## Biomass Dynamics (Logistic Discrete) Model for Snow Crab


**Logistic Surplus Production**: Parsimonious 2-parameter biomass production skeleton ($r$, $K$). Appropriate when recruitment or age-structure indices are unavailable or uninformative. Offers closed-form biological reference points ($B_{msy} = K/2$, $F_{msy} = r/2$).


### State Variable

Rather than tracking absolute biomass $B_t$ (kt), the model tracks **relative
biomass** (depletion):

$$m_t = \frac{B_t}{K}$$

where $K$ is the carrying capacity. $m_t \in (0, \infty)$ but is biologically
expected in $(0, 1]$.

---

### Process Model: Discrete Logistic Growth

The deterministic skeleton follows the Schaefer surplus-production equation in
discrete time. The expected relative biomass at time $t$ is:

$$
\hat{m}_t = m_{t-1} + r \cdot m_{t-1} \cdot (1 - m_{t-1}) - \frac{C_{t-1}}{K}
$$

where:
- $r$: intrinsic rate of population growth.
- $K$: carrying capacity (kt).
- $C_{t-1}$: total catch (removals, kt) in year $t-1$.

In the Gaussian variant $\hat{m}_t$ is not floored; the truncated normal is
defined for any real mean, and a floor would bias $m_t$ upward in years where
$C_{t-1}/K$ exceeds surplus production. The lognormal variant uses
$\max(mlim_{low}, \hat{m}_t)$ only because $\log$ requires a positive argument.

#### Catch-feasibility constraints (Gaussian variant)

The state $m_t$ is bounded by the truncation, but $\hat{m}_t$ is not. A value
$\hat{m}_t \le mlim_{low}$ means the recorded catch exceeded the stock plus its
production for that draw of $(K, r, m_{t-1})$. Truncation alone does not
penalise this: as $\mu = \hat{m}_t$ falls below $L = mlim_{low}$, the
renormalised density tends to $\lambda e^{-\lambda (m_t - L)}$ with
$\lambda = (L - \mu)/\sigma^2$, which places $m_t$ at the bound and assigns no
cost to the infeasible mean. Because catches were taken, such parameter
combinations are excluded explicitly:

$$
\hat{m}_t > mlim_{low} \quad (t = 2, \dots, n_T), \qquad
m_{t} K > C_{\mathrm{pen},t} \quad (\text{survey years}).
$$

Each inequality is a `dconstraint` node observed as 1 (`feasible_proc`,
`feasible_fall`; built by `feasibility_constraint_data()`), so the posterior
density is zero outside the feasible region. Initial values must satisfy both
constraints. Retrospective fits regenerate the constraint vectors for the
truncated series.

#### Process Error — Gaussian variant

$$
m_t \sim \mathrm{TruncatedNormal}\!\left(\mu = \hat{m}_t,\; \sigma = bpsd,\; mlim_{low},\; mlim_{high}\right)
$$

Here $bpsd$ is an **absolute SD** on the depletion scale, identical in
interpretation to the Julia/Turing implementation. The bounds $[mlim_{low}, mlim_{high}] = [0.01, 0.99]$
prevent biologically implausible population collapse or runaway states, matching Julia's constraint space.

#### Process Error — Lognormal variant

$$
m_t \sim \mathrm{TruncatedLogNormal}\!\left(\mu_\ell,\; \sigma_\ell = bpsd,\; mlim_{low},\; mlim_{high}\right)
\quad \text{where} \quad
\mu_\ell = \log(\max(mlim_{low}, \hat{m}_t)) - \tfrac{1}{2}\, bpsd^2
$$

Here $bpsd$ is the **log-scale SD** ($\sigma_\ell$). The bias correction
$-\tfrac{1}{2}\sigma_\ell^2$ ensures the median of $m_t$ equals
$\hat{m}_t$. The process is bounded on $[mlim_{low}, mlim_{high}]$ to match the Julia reference domain.
A prior mean of $bpsd = 0.2$ implies a 95% multiplicative interval of
approximately $[0.67,\,1.49]\hat{m}_t$.

---

### Observation Model

The survey biomass index $S_{obs,t}$ is linked to the depletion state via
catchability $q_1$. A portion of the annual catch has already been taken at the
time of the survey; the removal penalty $C_{\mathrm{pen},t}$ (kt) is divided
by $K$ to convert it to the depletion scale before subtraction:

$$
\mu_{S,t} = \frac{m_t - C_{\mathrm{pen},t}/K}{q_1}
$$

Both $m_t$ and $C_{\mathrm{pen},t}/K$ are dimensionless, so $q_1$ is
**dimensionless** and $O(1)$. The prior $q_1 \sim \mathrm{TruncNormal}(1, 0.5)$
is therefore physically meaningful for all subareas regardless of their $K$
scale (cfanorth $K\approx 5$ kt, cfasouth $K\approx 65$ kt). This formulation
is identical to the Julia/Turing reference implementation.

#### Likelihood — Gaussian and Lognormal variants

$$
S_{obs,t} \sim \mathrm{Normal}\!\left(\mu_{S,t},\; \sigma = bosd\right)
$$

$bosd$ is an **absolute SD** on the normalised $[0,1]$ index scale.
Observation error is additive Gaussian noise around the predicted post-harvest depletion.

---

### Prior Assumptions

#### Gaussian variant
| Parameter                | Prior                                                                         | Interpretation                                                                   |
| --------------------------| -------------------------------------------------------------------------------| ----------------------------------------------------------------------------------|
| $K$                      | $\mathrm{TruncatedNormal}(\mu_K, \sigma_K, 0.01 K_0, \infty)$                 | carrying capacity (kt)                                                           |
| $r$                      | $\mathrm{TruncatedNormal}(\mu_r, \sigma_r, 0.01, \infty)$                     | intrinsic growth rate                                                            |
| $q_1$                    | $\mathrm{TruncatedNormal}(\mu_{q1}, \sigma_{q1}, 0.01, \infty)$               | dimensionless catchability; $q_1=1$ means index $\approx$ post-harvest depletion |
| $bpsd$ (process abs. SD) | $\mathrm{TruncatedNormal}(\mu_{bpsd}, \sigma_{bpsd}, 0.01, \infty)$           | SD of $m_t$ residuals                                                            |
| $bosd$ (obs abs. SD)     | $\mathrm{TruncatedNormal}(\mu_{bosd}, \sigma_{bosd}, 0.01, \infty)$           | SD of survey index residuals on $[0,1]$ scale                                    |
| $m_t$                    | $\mathrm{TruncatedNormal}(\hat{m}_t, \sigma_{bpsd}, mlim_{low}, mlim_{high})$ | relative depletion                                                               |

#### Lognormal variant
| Parameter                     | Prior                                                                                |
| -------------------------------| --------------------------------------------------------------------------------------|
| $K$                           | $\mathrm{TruncatedLogNormal}(\mu_{\ell,K}, \sigma_{\ell,K}, 0.01 K_0, \infty)$       |
| $r$                           | $\mathrm{TruncatedLogNormal}(\mu_{\ell,r}, \sigma_{\ell,r}, 0.01, \infty)$           |
| $q_1$                         | $\mathrm{TruncatedLogNormal}(\mu_{\ell,q1}, \sigma_{\ell,q1}, 0.01, \infty)$         |
| $bpsd$ (process log-SD)       | $\mathrm{TruncatedLogNormal}(\mu_{\ell,bpsd}, \sigma_{\ell,bpsd}, 0.01, \infty)$     |
| $bosd$ (obs SD, survey scale) | $\mathrm{TruncatedLogNormal}(\mu_{\ell,bosd}, \sigma_{\ell,bosd}, 0.01, \infty)$     |
| $m_t$                         | $\mathrm{TruncatedLogNormal}(\mu_\ell, \sigma_{\ell,bpsd}, mlim_{low}, mlim_{high})$ |

The lognormal parameters $(\mu_\ell, \sigma_\ell)$ are computed from the intended
mean $\mu$ and SD $\sigma$ of the prior via the moment-matching identities:

$$
\sigma_\ell = \sqrt{\log\!\left(1 + \left(\frac{\sigma}{\mu}\right)^2\right)},
\qquad
\mu_\ell = \log(\mu) - \tfrac{1}{2}\sigma_\ell^2
$$

---

### Derived Reference Points

From $r$ and $K$, MSY reference points are derived analytically (Schaefer model):

$$
MSY = \frac{r K}{4}, \qquad
B_{MSY} = \frac{K}{2}, \qquad
F_{MSY} = \frac{r}{2}
$$

The Harvest Control Rule (HCR) classifies fall biomass $B_{\mathrm{fall}} = m_t K - C_t$
relative to $B_{MSY}$ thresholds:

| Zone     | Condition                         |
| ----------| -----------------------------------|
| Critical | $B_{\mathrm{fall}} \le K/4$       |
| Cautious | $K/4 < B_{\mathrm{fall}} \le K/2$ |
| Healthy  | $B_{\mathrm{fall}} > K/2$         |

---

### Exploitation Rate and Fishing Mortality

Posterior exploitation rates ($HR$) and instantaneous fishing mortality ($F$) are
derived from the posterior biomass trajectories and observed removals:

$$
HR_t = \frac{C_t}{B_t}, \qquad F_t = -\log(1 - HR_t)
$$

In accordance with analytical standards, $HR_t$ is not artificially clamped or
truncated. Values $HR_t \ge 1$ indicate that observed removals exceeded posterior
biomass in that sample and year; when encountered during diagnostics, an
informative warning is emitted to flag potential model misspecification or
data anomalies rather than masking the condition via numerical truncation.

---

### Prior Sampling for Diagnostics

Prior predictive samples are generated by sampling from the corresponding
truncated-normal (Gaussian variant) or lognormal (lognormal variant)
distributions used in the model through `sample_priors_nimble()`. When
parameter bounds are defined in test configurations, identical truncated bounds
are applied during prior sampling to ensure that the empirical prior display
strictly corresponds to the model specification.

---

### Retrospective Analysis

To assess model consistency, years of terminal data are iteratively removed and
the model is re-estimated. **Mohn's Rho** ($\rho$) measures systematic
retrospective bias:

$$
\rho = \frac{1}{n} \sum_{i=1}^{n} \frac{\theta_i - \theta_{\mathrm{full}}}{\theta_{\mathrm{full}}}
$$

where $\theta_i$ is the parameter estimate with $i$ years truncated and
$\theta_{\mathrm{full}}$ is the full-model estimate. Values $|\rho| > 0.2$
generally indicate problematic retrospective patterns (Mohn 1999).

---

### Key References

- **Schaefer, M. B.** (1954). Some aspects of the dynamics of populations important to the management of the commercial marine fisheries. *Bulletin of the Inter-American Tropical Tuna Commission*, 1(2), 25–56.
- **Meyer, R., & Millar, R. B.** (1999). BUGS in Bayesian stock assessments. *Canadian Journal of Fisheries and Aquatic Sciences*, 56(6), 1078–1087. https://doi.org/10.1139/f99-043
- **Millar, R. B., & Meyer, R.** (2000). Non-linear state space modelling of fisheries biomass dynamics by using Metropolis-Hastings within-Gibbs sampling. *Journal of the Royal Statistical Society: Series C (Applied Statistics)*, 49(3), 327–342. https://doi.org/10.1111/1467-9876.00195
- **Mohn, R.** (1999). The retrospective problem in stock assessment: An investigation of methods for computing bias. *ICES Journal of Marine Science*, 56(4), 473–488. https://doi.org/10.1006/jmsc.1999.0486
- **Gelman, A., & Rubin, D. B.** (1992). Inference from iterative simulation using multiple sequences. *Statistical Science*, 7(4), 457–472. https://doi.org/10.1214/ss/1177011136
- **Watanabe, S.** (2010). Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory. *Journal of Machine Learning Research*, 11, 3571–3594.
- **DFO.** (2009). *A Fishery Decision-Making Framework Incorporating the Precautionary Approach*. Fisheries and Oceans Canada, Science Advisory Report 2009/005.
- **Choi, J. S.** (2020). Assessment of Eastern Nova Scotia (Area 23, 24) and 4X Snow Crab (*Chionoecetes opilio*). *DFO Canadian Science Advisory Secretariat Research Document* 2020/nnn.
- **de Valpine, P., Turek, D., Paciorek, C. J., Anderson-Bergman, C., Temple Lang, D., & Bodik, R.** (2017). Programming with models: writing statistical algorithms for general model structures with NIMBLE. *Journal of Computational and Graphical Statistics*, 26(2), 403–413. https://doi.org/10.1080/10618600.2016.1172487



## Delay-Difference Model for Snow Crab

**Deriso-Schnute Delay-Difference (`delay_difference_run.R`)**: Explicitly separates natural mortality ($M$), individual body growth ($\rho$), and recruit size ($w_k$), modeling recruitment variability via an AR(1) process. Utility when pre-recruit survey data or size composition data allow disentangling growth from cohort pulses. 



### Process Model: Deriso-Schnute Delay-Difference

Unlike the discrete logistic model which relies solely on the intrinsic rate of growth ($r$) and carrying capacity ($K$), the delay-difference model unpacks population dynamics into survival, somatic growth, and recruitment.

#### Absolute Biomass Equation
The classic Deriso-Schnute delay-difference model for absolute biomass ($B_t$) follows the structure:
$$
B_t = (1 + \rho) e^{-M} B_{t-1} - \rho e^{-2M} B_{t-2} - (1 + \rho) e^{-M} C_{t-1} + \rho e^{-2M} C_{t-2} + w_k R_t
$$

Where:
- $M$: Instantaneous natural mortality rate.
- $\rho$: Brody growth coefficient (defining somatic growth of individuals).
- $C_t$: Catch/removals at time $t$.
- $w_k$: Average weight of new recruits entering the exploitable population.
- $R_t$: Number of recruits entering the population at time $t$.

#### Relative Biomass Scaling
To maintain numerical stability and consistency with the previous logistic models, the model tracks **relative biomass** ($m_t = B_t / K$). By dividing the absolute equation by carrying capacity ($K$), the core process model becomes:

$$
m_t = (1 + \rho) e^{-M} m_{t-1} - \rho e^{-2M} m_{t-2} - \frac{(1 + \rho) e^{-M} C_{t-1}}{K} + \frac{\rho e^{-2M} C_{t-2}}{K} + \frac{w_k R_t}{K}
$$

**Process Error**: The relative biomass is assumed to follow a log-normal distribution with process standard deviation `bpsd`.
$$
m_t \sim Lognormal(\log(m_{process\_mean}) - 0.5 \times bpsd^2, bpsd)
$$

### Recruitment Modeling

Recruitment $R_t$ is highly variable and often non-stationary. Rather than estimating independent values for each year (which would overparameterize the model), it is constrained by an **Autoregressive (AR1) Random Walk** in log-space.

#### State Equation for Recruitment
$$
\log(R_t) \sim Normal(\mu_R + \phi_R (\log(R_{t-1}) - \mu_R), \sigma_R)
$$

Where:
- $\mu_R$: The baseline average (drift) of log-recruitment.
- $\phi_R$: The AR(1) coefficient, representing temporal autocorrelation in recruitment success. Constrained to $[-0.99, 0.99]$ to ensure stationarity.
- $\sigma_R$: The standard deviation (process error) of the random walk.

### Observation Models

The latent states ($m_t$ and $R_t$) are anchored by observation indices.

#### Survey Biomass Likelihood
The observed survey biomass ($S_{obs, t}$) is a function of the relative biomass,
catchability $q_1$, and the catch penalty $C_{\mathrm{pen},t}$ (kt) scaled by $K$:

$$
S_{obs, t} \sim Normal\left(\frac{m_t - C_{\mathrm{pen},t}/K}{q_1},\; bosd\right)
$$

Both $m_t$ and $C_{\mathrm{pen},t}/K$ are dimensionless; $q_1$ is therefore
dimensionless and $O(1)$, consistent with the logistic model formulation.
$bosd$ is an absolute SD on the normalised $[0,1]$ survey-index scale.

#### Recruitment Likelihood
To properly estimate the recruitment random walk, an observed recruitment index (`R_observed`) is provided. The observation likelihood links the true unobserved recruits $R_t$ to the index using a catchability/scaling parameter $q_R$:

$$
R_{observed, t} \sim Lognormal\left(\log(q_R \times R_t) - 0.5 \times sd_{R\_obs}^2, sd_{R\_obs}\right)
$$

### Prior Assumptions

Several critical parameters have prior distributions defined to inform their estimation:
- **$M$ (Natural Mortality)**: $M \sim \mathrm{LogNormal}(0.2, 0.05)$. Baseline natural mortality ~0.2 yr$^{-1}$; the model may deviate if data warrants.
- **$\rho$ (Brody Growth Coefficient)**: $\rho \sim \mathrm{LogNormal}(0.8, 0.1)$. Slow somatic growth decay typical of longer-lived crustaceans.
- **$w_k$ (Recruit Weight)**: $w_k \sim \mathrm{LogNormal}(0.5, 0.1)$. Mean weight of incoming recruits.
- **$\phi_R$ (AR1 Coefficient)**: $\phi_R \sim \mathrm{Uniform}(-0.99, 0.99)$. Stationarity-ensuring uninformative prior on recruitment autocorrelation.
- **$q_R$ (Recruitment Catchability)**: $q_R \sim \mathrm{LogNormal}(0, 1)$. Scales the latent $R_t$ to the observed recruitment index.
- **$sd_{R\_obs}$ (Recruitment Obs. SD)**: $sd_{R\_obs} \sim \mathrm{LogNormal}(0, 0.5)$. Observation error on the recruitment index.
- **$K, q_1, bpsd, bosd$**: Retain priors from the logistic implementation.

All `LogNormal` priors are parameterised by `(meanlog, sdlog)` per Nimble
convention; the intended natural-scale mean $\mu$ and SD $\sigma$ are
converted via moment-matching:
$\sigma_\ell = \sqrt{\log(1 + (\sigma/\mu)^2)}$,
$\mu_\ell = \log(\mu) - \tfrac{1}{2}\sigma_\ell^2$.

### Retrospective Analysis

The delay difference model incorporates a **retrospective analysis** to identify systematic biases (such as retrospective patterns where the model consistently over- or under-estimates biomass as new data is added).

#### Methodology
1. The model is first fit to the full available dataset.
2. The dataset is then iteratively truncated by removing the terminal year of data (e.g., up to `truncate_years_count = 5` years).
3. The model is refitted to each truncated dataset.

#### Mohn's Rho Calculation
For the target parameter (typically $K$), Mohn's $\rho$ is calculated to quantify the retrospective pattern:
$$
\rho = \frac{1}{n} \sum_{i=1}^{n} \frac{\theta_{i} - \theta_{full}}{\theta_{full}}
$$
Where:
- $\theta_i$ is the parameter estimate from the truncated model (ending in year $i$).
- $\theta_{full}$ is the parameter estimate from the full model for the corresponding year.

A large positive or negative Mohn's $\rho$ indicates a strong retrospective bias.

---

### Reference Points and Harvest Control Rules

Unlike the Schaefer surplus production model, the Deriso-Schnute delay-difference formulation has no closed-form MSY equation. Equilibrium biomass $B_{eq}(F)$ and equilibrium yield $Y_{eq}(F)$ under constant instantaneous fishing mortality $F$ are derived analytically:

$$
B_{eq}(F) = \frac{w_k R_{eq}}{1 - (1+\rho)e^{-(M+F)} + \rho e^{-(2M+F)}}
$$

$$
Y_{eq}(F) = (1 - e^{-F}) B_{eq}(F)
$$****

The biological reference point $F_{0.1}$ is computed numerically as the fishing mortality where the slope of the equilibrium yield curve drops to 10% of its value at the origin:

$$
\left.\frac{\partial Y_{eq}}{\partial F}\right|_{F=F_{0.1}} = 0.1 \times \left.\frac{\partial Y_{eq}}{\partial F}\right|_{F=0}
$$

Reference biomass is defined as $B_{msy} = B_{eq}(F_{0.1})$ and $\text{MSY} = Y_{eq}(F_{0.1})$. Stock status zones follow the Precautionary Approach framework:
- **Critical Zone**: $B_t \le B_{msy} / 2$
- **Cautious Zone**: $B_{msy} / 2 < B_t \le B_{msy}$
- **Healthy Zone**: $B_t > B_{msy}$

This is computed for each MCMC posterior draw via `delay_difference_reference_points_R()`.


### Key References

- **Deriso, R. B.** (1980). Harvesting strategies and parameter estimation for an age-structured model. *Canadian Journal of Fisheries and Aquatic Sciences*, 37(2), 268–282. https://doi.org/10.1139/f80-034
- **Schnute, J.** (1985). A general theory for analysis of catch and effort data. *Canadian Journal of Fisheries and Aquatic Sciences*, 42(3), 414–429. https://doi.org/10.1139/f85-056
- **Hilborn, R., & Walters, C. J.** (1992). *Quantitative Fisheries Stock Assessment: Choice, Dynamics and Uncertainty*. Chapman and Hall, New York.
- **Quinn, T. J., & Deriso, R. B.** (1999). *Quantitative Fish Dynamics*. Oxford University Press, New York.
- **Mohn, R.** (1999). The retrospective problem in stock assessment: An investigation of methods for computing bias. *ICES Journal of Marine Science*, 56(4), 473–488. https://doi.org/10.1006/jmsc.1999.0486
- **DFO.** (2009). *A Fishery Decision-Making Framework Incorporating the Precautionary Approach*. Fisheries and Oceans Canada, Science Advisory Report 2009/005.
- **de Valpine, P., Turek, D., Paciorek, C. J., Anderson-Bergman, C., Temple Lang, D., & Bodik, R.** (2017). Programming with models: writing statistical algorithms for general model structures with NIMBLE. *Journal of Computational and Graphical Statistics*, 26(2), 403–413. https://doi.org/10.1080/10618600.2016.1172487


```
   

```
  aulab ="cfanorth"
  include( joinpath(project_directory, "logistic_discrete.jl" ) )

  aulab ="cfasouth"
  include( joinpath(project_directory, "logistic_discrete.jl" ) )

  aulab ="cfa4x"
  include( joinpath(project_directory, "logistic_discrete.jl" ) )


```


## Save fishery model results for quick load 

```{r}
#| eval: false
#| output: false

# reload parameters   
require(aegis)  # basic helper tools
year.assessment = 2025  # change this as appropriate
p = bio.snowcrab::load.environment( year.assessment=year.assessment )  # set up initial settings


year_assessment = min(p$yrs)
year_start = max(p$yrs)
data_loc = p$datadir
mau = "region"
model_variation = "logistic_discrete_historical"
  

# update new results for fast access
snowcrab_load_key_results_to_memory( p=p, todo="fishery_model", mau=mau, redo=TRUE ) 

```


[And that is it. Now we continue to look at some of the main results](06_assessment_summary.md).

 