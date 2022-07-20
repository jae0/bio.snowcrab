Historical approach

The modelling approach that has been used in the past was a simple biomass dynamics model to empirically estimate two biological parameters, the intrinsic rate of increase (r) and carrying capacity (K). These are the minimal parameters required to estimate the biological reference points that help guage the status of a population in what is currently referred to  as a Precautionary Approach (*add* *references*). Maximum likelihood and least squares approaches are frequently used for estimation, however, when used with snow crab, these methods have been found to be numerically unstable and not operational. For snow crab, a Bayesian approach to parameter estimation was used where informative priors for these parameters helped to facillitate stable estimation (*add* *references*). The inference engines use were JAGS (Just Another Gibbs Sampler; *reference*) prior to 20XX and thereafter STAN (reference), due to use of the more efficient NUTS sampler in the latter.

The estimation approach was to use a discrete form of the biomass dynamics ("process") model, combined with an observation model to relate process means to observed means. The model specification is as follows:

$b_{t+1}=r \cdot b_t(1-b_t/K) $

$r \sim N(1, 0.1)$

$K \sim N(\mu, 0.1 * \mu)$

for each areal unit. The prior mean of the carrying capacity are based upon prior historical analyses when estimation was based upon various kinds of kriging, GAM-based interpolations and stmv-based procedures (*add references).*

The obervation model was an extremely simple one where observations ($Y_t$) were  related to the mean of the latent process model at a given year ($t$), by a linear (scalar) coefficient  $(q)$, often reffered to as a of "catchability" coefficient and Gaussian observation error ($\sigma_o $):

$Y_t \sim N(q \cdot b_t , \sigma_o)$

and a latent process's mean and error:

$b_t \sim N(b_t^*, \sigma_p)$

The strength of this approach is its simplicty; but his simplicity is also its weakness.


In particular, the loigsitc for is a phenomelogical formulation, that is a model without mechanism. Biomass is to create biomass at a given rate (r) to a limit of K. Structurally, the model amounts to a second-order Taylor-series expansion of the variability of dynamics and is in essence a heuristic smoother/interpolation operator. 

In the case of snow crab, and specifically the fishable component, the biomass of the fishable component at one time slice does not create the fishable biomass in the next, but rather is a composite of remaining biomass, recruitment, fishing and other sources of mortality as well as growth.  These rates are not directly related to fishable biomass in the previous year, except perhaps fishery mortlaity. Ultimately, it is a heuristic shortcut to mechanisms without directly treating the mechanisms.  To make predictions more useful,   and biological parameters more meaningful, we must identify how these parameters relate to the biological mechanisms involved to explain why total biomass should have a stable upper limit or some rate of growth. This mechanistic understanding is possible if we shift the focus to numbers and numerical density and average weight and moulting and mortality. That is, a mechanistic model of numbers permits this understanding.
