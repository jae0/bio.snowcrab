Historical approach

The modelling approach that has been used in the past was a simple biomass dynamics model to empirically estimate two biological parameters, the intrinsic rate of increase (r) and carrying capacity (K). These are the minimal parameters required to estimate the biological reference points that help gauge the status of a population in what is currently referred to as a Precautionary Approach (*add* *references*). Maximum likelihood and least squares approaches are frequently used for estimation, however, when used with snow crab, these methods have been found to be numerically unstable and not operational. For snow crab, a Bayesian approach to parameter estimation was used where informative priors for these parameters helped to facilitate stable estimation (*add* *references*). The inference engines use were JAGS (Just Another Gibbs Sampler; *reference*) prior to 20XX and thereafter STAN (reference), due to use of the more efficient NUTS sampler in the latter.

The estimation approach was to use a discrete form of the biomass dynamics ("process") model, combined with an observation model to relate process means to observed means. The model specification is as follows:

$b_{t+1}=r \cdot b_t(1-b_t/K) $

$r \sim N(1, 0.1)$

$K \sim N(\mu, 0.1 * \mu)$

for each areal unit. The prior mean of the carrying capacity are based upon prior historical analyses when estimation was based upon various kinds of kriging, GAM-based interpolations and stmv-based procedures (*add references).*

The observation model was an extremely simple one where observations ($Y_t$) were related to the latent biomass at a given year ($t$), by a linear (scalar) coefficient  $(q)$, often referred to as a "catchability" coefficient and Gaussian observation error ($\sigma_o $):

$Y_t \sim N(q \cdot b_t, \sigma_o)$

and a latent process's mean and error:

$b_t \sim N(b_t^*, \sigma_p)$

The strength of this approach is its simplicity; but his simplicity is also its weakness.

In particular, the logistic for is a phenomenological formulation, that is a model without mechanism. Biomass is to create biomass at a given rate (r) to a limit of K. Structurally, the model amounts to a second-order Taylor-series expansion of the variability of dynamics and is in essence a heuristic smoother/interpolation operator.

In the case of the biomass of the fishable component of snow crab, it does not create the fishable biomass in the next, but rather is a composite of remaining biomass, recruitment, fishing and other sources of mortality as well as growth.  These rates are not directly related to fishable biomass in the previous year, except perhaps fishery mortality. Ultimately, it is a heuristic smoothing operator without directly treating the mechanisms.  To make predictions more useful, and biological parameters more meaningful in the context of the Precautionary approach, we must identify how these parameters relate to the biological mechanisms involved to explain why total biomass should have a stable upper limit or some rate of growth. This mechanistic understanding is possible only if we shift the focus to numbers and numerical density and average weight and moulting and mortality. That is, a mechanistic model of numbers permits this understanding.

For such a purpose, we consider the following six-compartment model:

  F8 -> F 
  F8 -> M4
  M4 -> M3
  M3 -> M2
  M2 -> M1
  M1 -> M0

The F indicate mature females of instar 9+, determined from size range and maturity. There is approximately an 8+ year period required for these females to produce the next generation of instar 9 females. This is represented by the notation 'F8'. A similar amount of time is required to produce males of instar 9 (denoted M4). Subsequently, instar 9 males transition to instar 10 (M3) and then to instar 11 (M2) and instar 12 (M1). M0 are considered fishable biomass (mature, > 95mm CW). A (small) fraction of instar 11 crab (M2) will be large enough to be considered legally fishable whereas all instar 12 (M1) and greater will be considered fishable size. M1, therefore, represents a composite group of all crab that have entered fishable size but are still immature. Some fraction of this group will mature into the fishable component (M0) in the next year. Others will continue to moult to instar 13 and higher, but for our purposes, they will be considered M0 only if they are morphologically mature. The numeric indices for males, therefore, indicate an approximate time in years before entry into the fishable component. The first two transition rates are akin to a birth rate (b) where the bounds in the rate are determined by egg production and survival to instar 9 (approximately 40 to 60 mm CW). The latter four are first order transition rates (v) that are bound by definition to be in the range (0 to 1). 

Each compartment has associated second order death rates, analogous to the logistic formulation, by decreasing to a first order process when numbers $s$ are are low but as a factor of $(1 - s/K)$. However,   