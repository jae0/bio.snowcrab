
# -----------------------------------------------------------------------------
# Refactored pipeline functions
# -----------------------------------------------------------------------------

prepare_fishery_model_data = function(
  year_assessment = 2025,
  aulab = "cfanorth",
  model_variation = "logistic_discrete_historical",
  inputs_dir = NULL,
  outputs_dir = NULL
) {
  require(aegis)  # basic helper tools
  require(terra) # function "geom" in snowcrab_parameters() .. remove it?

  p = bio.snowcrab::load.environment( year.assessment=year_assessment )  # set up initial settings

  pkgs = c(
    "nimble", "nimbleHMC", "ggplot2", "dplyr", "tidyr", "igraph", "pracma",
    "data.table","coda", "numDeriv"
  )
  for (pk in pkgs) { if (!require(pk, character.only = TRUE)) install.packages(pk, deps=TRUE) }
  for (pk in pkgs) { library(pk, character.only = TRUE) }

  yrs = 2000:year_assessment

  bio_data_directory = project.datadirectory( "bio.snowcrab", "modelled", "default_fb")

  if (is.null(inputs_dir)) {
    inputs_dir = bio_data_directory
  }

  if (is.null(outputs_dir)) {
    outputs_dir = project.datadirectory( "bio.snowcrab", "fishery_model", year_assessment, model_variation )
  }

  o = readRDS(  file.path(inputs_dir, "biodyn_biomass.rds") )

  # Data preparation (translated from Julia code)
  Y_df = o$Y[o$Y$yrs %in% yrs, ]  # index of abundance (kt)
  removals_df = o$L[o$L$yrs %in% yrs, ]  # landings (kt)

  nT = length(yrs)
  nP = 0  # number of predictions into future (with no fishing)
  nM = nP + nT  # total number of prediction years
  dt = 1  # time resolution of solutions .. discrete annual so 1 year
  nS = 1 # no state variables, not used
  no_digits = 3  # time floating point rounding
  smallnumber = 1.0e-6 # floating point value sufficient to assume 0 valued

  # "survey index"
  S_obs = Y_df[[aulab]]  # select subarea inputs for modelling
  S_obs = S_obs / max(S_obs, na.rm = TRUE) # Normalise

  # id index for each subarea
  ki = switch(
    aulab,
    "cfanorth" = 1,
    "cfasouth" = 2,
    "cfa4x"    = 3,
    0  # default
  )

  Kmu = c(5.0, 65.0, 1.5) # vector of guess of magnitude of carrying capacity (kt)
  kmu = Kmu[ki]

  survey_time = Y_df$yrs    # time of observations for survey .. for plotting

  # adjust survey time to be approximately in spring or autumn
  spring_idx = which(survey_time < 2004)
  fall_idx   = which(survey_time >= 2004)
  survey_time[spring_idx] = survey_time[spring_idx] + 5.0/12.0  # arbitrary ..  "spring"
  survey_time[fall_idx]   = survey_time[fall_idx] + 10.0/12.0   # time of survey in "fall" ; all of 4X sep-mar gets called "fall" for simplicity in plot
  survey_time = round(survey_time, digits=no_digits)    # time of observations for survey

  removed_val = removals_df[[aulab]]

  # predictions ("m") are "spring" .. arbitrarily setting to 4/12
  predtime =  4.0/ 12.0
  fall_time_const = 11.9 / 12.0

  prediction_time = floor( seq(min(yrs), max(yrs)+nP, by = 1) ) + round(predtime/dt, digits=no_digits)   # april (m== spring)
  yrs_pred_report_idx = which((prediction_time >= 1999.0) & (prediction_time <= (as.numeric(year_assessment)+nP +1.0)))
  prediction_time_ss = prediction_time[yrs_pred_report_idx]

  fall_time = floor( seq(min(yrs), max(yrs)+nP, by = 1) ) + round(fall_time_const/dt, digits=no_digits)   # april (m== spring)
  yrs_pred_report_fall_idx = which((fall_time >= 1999.0) & (fall_time <= (as.numeric(year_assessment)+nP +1.0)))
  fall_time_ss = fall_time[yrs_pred_report_fall_idx]

  iok = which(!is.na(S_obs))  # index of data locations

  # basic parameters for "logistic_discrete_historical"
  PM = list(
    yrs = yrs,
    nS = nS,
    nT = nT,
    nP = nP,
    nM = nM,
    K = c(kmu, kmu * 0.2 ), # Vector for Nimble's dnorm parameters
    r = c(1.0, 0.5 ),
    q1 = c(1.0, 0.5 ),
    bpsd = c(0.2, 0.1 ),
    bosd = c(0.2, 0.1 ),
    mlim = c(0.01, 0.99 ), # bounds on relative depletion m_t (matches Julia reference)
    m0 = c(0.8, 0.2 ),
    removed = removed_val,
    S = S_obs,
    iok = iok,
    yeartransition = 0
  )

  if (aulab == "cfanorth" | aulab == "cfasouth") {
    xxyy = which(PM$yrs == 2004)
    if (length(xxyy) > 0) {
      PM$yeartransition = xxyy[1]
    }
  }

  # Define initial values for the sampler
  initial_values = list(
    K = kmu,
    r = 1.0,
    q1 = 1.0,
    bpsd = 0.1,
    bosd = PM$bosd[1],   # prior mean; survey-index scale (not scaled by K)
    m = rep(0.8, nM)
  )

  # Data exclusively holds observed or dynamic vectors used for likelihood/loop calculations
  removed_penalty_for_S_observed = compute_removed_penalty(PM)

  nimble_data = list(
    S_observed = S_obs[iok],
    removed = PM$removed
  )

  if (model_variation == "delay_difference") {
    if (exists("R_observed_kt", o)) {
      nimble_data$R_observed = o$R_observed_kt
    } else {
      stop("recruitment data not found")
    }
  }

  return(list(
    PM = PM,
    initial_values = initial_values,
    nimble_data = nimble_data,
    removed_penalty_for_S_observed = removed_penalty_for_S_observed,
    survey_time = survey_time,
    prediction_time_ss = prediction_time_ss,
    fall_time_ss = fall_time_ss,
    outputs_dir = outputs_dir,
    inputs_dir = inputs_dir,
    S_obs = S_obs,
    yrs = yrs,
    year_assessment = year_assessment,
    aulab = aulab,
    model_variation = model_variation
  ))
}


#' Data vectors for catch-feasibility constraints
#'
#' @description
#' Observed removals impose hard limits on the logistic state-space model.
#' The expected relative biomass after removals and the post-harvest biomass at
#' the time of the survey must both be positive:
#' \deqn{\hat{m}_t = m_{t-1} + r m_{t-1}(1 - m_{t-1}) - C_{t-1}/K > mlim_{low},
#'   \quad t = 2, \dots, n_T}
#' \deqn{m_{i_j} - C_{pen,j}/K > 0, \quad j = 1, \dots, n_S}
#' Each inequality is a \code{dconstraint} node observed as 1, so parameter
#' combinations implying removals larger than the stock have zero posterior
#' density instead of being clamped or absorbed by the truncation bound.
#'
#' @param nT Integer number of historical years (\eqn{n_T \ge 2}).
#' @param n_observed_S Integer number of survey observations (\eqn{n_S \ge 1}).
#' @return A list with \code{feasible_proc} (length \eqn{n_T - 1}) and
#'   \code{feasible_fall} (length \eqn{n_S}), all elements equal to 1.
feasibility_constraint_data = function(nT, n_observed_S) {
  if (nT < 2) stop("feasibility_constraint_data: nT must be >= 2.")
  if (n_observed_S < 1) stop("feasibility_constraint_data: n_observed_S must be >= 1.")
  list(
    feasible_proc = rep(1, nT - 1),
    feasible_fall = rep(1, n_observed_S)
  )
}


#' Prepare constants, data, dimensions, and inits for a Nimble model
#'
#' @description
#' Systematically constructs model components for Gaussian, Lognormal, or
#' Delay-Difference fishery models based on the unified \code{setup} object.
#'
#' @param setup List returned by \code{prepare_fishery_model_data()}.
#' @param model_type Character string: "lognormal", "gaussian", or "delay_difference". 
#' @param seed Integer seed for reproducibility of mock recruitment data.
#' @return A list containing \code{code}, \code{constants}, \code{data}, \code{dimensions},
#'   \code{inits}, \code{monitors}, and \code{buildDerivs}.
prepare_nimble_model_inputs = function(
    setup,
    model_type = c("lognormal", "gaussian", "delay_difference"), 
    seed = 42
) {
  
  set.seed(seed)

  model_type = match.arg(model_type)
  PM = setup$PM
  nimble_data = setup$nimble_data
  initial_values = setup$initial_values
  removed_penalty = setup$removed_penalty_for_S_observed

  # Base constants common across all models
  base_constants = list(
    nT = as.integer(PM$nT),
    nM = as.integer(PM$nM),
    m_idx_for_S = as.integer(PM$iok),
    n_observed_S = length(PM$iok),
    removed_penalty_for_S_observed = as.numeric(removed_penalty)
  )

  if (model_type == "gaussian") {

    #' state-space Schaefer surplus production model: truncated normal priors,
    #' process error bounded by \eqn{[mlim_{low}, mlim_{high}]}, and Gaussian observation error:
    #' \deqn{m_t = m_{t-1} + r m_{t-1}(1 - m_{t-1}) - \frac{C_{t-1}}{K}}
    #' \deqn{m_t \sim \text{TruncNormal}(m_t, \sigma_{proc}, mlim_{low}, mlim_{high})}
    #' \deqn{S_j \sim \text{Normal}\left(\frac{m_{i_j} - C_{pen,j}/K}{q_1}, \sigma_{obs}\right)}
    #' with catch-feasibility constraints (see \code{feasibility_constraint_data}):
    #' \deqn{\hat{m}_t > mlim_{low}, \qquad m_{i_j} K > C_{pen,j}}
    
    code = nimbleCode({
      
      # note: Nimble models require use of "->"

      # Priors: normal, truncated below only (identical to the Turing reference)
      K ~ T(dnorm(mean = mean_K, sd = sd_K), lower_K, Inf)
      r ~ T(dnorm(mean = mean_r, sd = sd_r), lower_par, Inf)
      q1 ~ T(dnorm(mean = mean_q1, sd = sd_q1), lower_par, Inf)
      bpsd ~ T(dnorm(mean = mean_bpsd_prior, sd = sd_bpsd_prior), lower_par, Inf)
      bosd ~ T(dnorm(mean = mean_bosd_prior, sd = sd_bosd_prior), lower_par, Inf)

      # Initial relative biomass m[1] bounded by [mlim_low, mlim_high]
      m[1] ~ T(dnorm(mean = mean_m0, sd = sd_m0), mlim_low, mlim_high)

      # Process model for historical observations. The expected value is not
      # floored: the truncated normal is defined for any real mean, and a
      # floor would bias m upward whenever C_{t-1}/K exceeds surplus production.
      # Feasibility: removals cannot exceed the stock plus its production.
      for (i in 2:nT) {
        mean_m_process[i] <- m[i - 1] + r * m[i - 1] * (1.0 - m[i - 1]) -
                             removed[i - 1] / K
        feasible_proc[i - 1] ~ dconstraint(mean_m_process[i] > mlim_low)
        m[i] ~ T(dnorm(mean = mean_m_process[i], sd = bpsd), mlim_low, mlim_high)
      }

      # Prediction loop into future (no fishery removals)
      if (nM > nT) {
        for (i in (nT + 1):nM) {
          mean_m_predict[i] <- m[i - 1] + r * m[i - 1] * (1.0 - m[i - 1])
          m[i] ~ T(dnorm(mean = mean_m_predict[i], sd = bpsd), mlim_low, mlim_high)
        }
      }

      # Observation likelihood on relative survey scale (no floor on mu_S).
      # Feasibility: post-harvest biomass at survey time must be positive.
      for (j in 1:n_observed_S) {
        feasible_fall[j] ~ dconstraint(
          m[m_idx_for_S[j]] * K > removed_penalty_for_S_observed[j]
        )
        mu_S[j] <- (m[m_idx_for_S[j]] - removed_penalty_for_S_observed[j] / K) / q1
        S_observed[j] ~ dnorm(mu_S[j], sd = bosd)
      }
    })
    
    buildDerivs = TRUE
    monitors = c("K", "r", "q1", "bpsd", "bosd", "m")

    if (is.null(PM$mlim) || length(PM$mlim) != 2 || PM$mlim[1] >= PM$mlim[2]) {
      stop("PM$mlim must be a length-2 vector c(lower, upper) with lower < upper.")
    }

    constants = c(base_constants, list(
      sd_K = PM$K[2], mean_K = PM$K[1], lower_K = PM$K[1] * 0.01,
      sd_r = PM$r[2], mean_r = PM$r[1],
      sd_q1 = PM$q1[2], mean_q1 = PM$q1[1],
      sd_bpsd_prior = PM$bpsd[2], mean_bpsd_prior = PM$bpsd[1],
      sd_bosd_prior = PM$bosd[2], mean_bosd_prior = PM$bosd[1],
      sd_m0 = PM$m0[2], mean_m0 = PM$m0[1],
      lower_par = 0.01,
      mlim_low = PM$mlim[1],
      mlim_high = PM$mlim[2]
    ))

    data = c(
      nimble_data,
      feasibility_constraint_data(PM$nT, length(PM$iok))
    )
    dimensions = list(
      S_observed = length(data$S_observed),
      removed = length(data$removed),
      feasible_proc = length(data$feasible_proc),
      feasible_fall = length(data$feasible_fall)
    )
    inits = initial_values

  } else if (model_type == "lognormal") {

    #' State-space Schaefer surplus production model with lognormal
    #' priors bounded by reference thresholds, lognormal process error
    #' bounded on the relative depletion domain \eqn{[mlim_{low}, mlim_{high}]},
    #' and Gaussian observation error:
    #' \deqn{m_t \sim \text{TruncLogNormal}\left(\ln(\max(mlim_{low}, \hat{m}_t)) - \frac{1}{2}\sigma_{proc}^2, \sigma_{proc}, mlim_{low}, mlim_{high}\right)}
    #' \deqn{S_j \sim \text{Normal}\left(\frac{m_{i_j} - C_{pen,j}/K}{q_1}, \sigma_{obs}\right)}
    
    code = nimbleCode({

      # note: Nimble models require use of "->"

      # Priors using lognormal with lower bounds
      K ~ T(dlnorm(meanlog = meanlog_K, sdlog = sdlog_K), lower_K, upper_K)
      r ~ T(dlnorm(meanlog = meanlog_r, sdlog = sdlog_r), 0.01, 3.0)
      q1 ~ T(dlnorm(meanlog = meanlog_q1, sdlog = sdlog_q1), 0.01, 3.0)
      bpsd ~ T(dlnorm(meanlog = meanlog_bpsd_prior, sdlog = sdlog_bpsd_prior), 0.01, Inf)
      bosd ~ T(dlnorm(meanlog = meanlog_bosd_prior, sdlog = sdlog_bosd_prior), 0.01, Inf)

      # Initial relative biomass m[1] bounded by [mlim_low, mlim_high]
      m[1] ~ T(dlnorm(meanlog = meanlog_m0, sdlog = sdlog_m0), mlim_low, mlim_high)

      # Process model for historical observations: bounded by [mlim_low, mlim_high]
      for (i in 2:nT) {
        meanlog_m_process[i] <- log(max(
          mlim_low,
          m[i - 1] + r * m[i - 1] * (1.0 - m[i - 1]) - removed[i - 1] / K
        )) - 0.5 * pow(bpsd, 2)
        m[i] ~ T(dlnorm(meanlog = meanlog_m_process[i], sdlog = bpsd), mlim_low, mlim_high)
      }

      # Prediction loop into future (no fishery removals)
      if (nM > nT) {
        for (i in (nT + 1):nM) {
          meanlog_m_predict[i] <- log(max(
            mlim_low,
            m[i - 1] + r * m[i - 1] * (1.0 - m[i - 1])
          )) - 0.5 * pow(bpsd, 2)
          m[i] ~ T(dlnorm(meanlog = meanlog_m_predict[i], sdlog = bpsd), mlim_low, mlim_high)
        }
      }

      # Observation likelihood on relative survey scale
      for (j in 1:n_observed_S) {
        mu_S[j] <- (m[m_idx_for_S[j]] - removed_penalty_for_S_observed[j] / K) / q1
        S_observed[j] ~ dnorm(mu_S[j], sd = bosd)
      }
    })

    buildDerivs = TRUE
    monitors = c("K", "r", "q1", "bpsd", "bosd", "m")

    K_p = get_lognorm_params(PM$K[1], PM$K[2])
    r_p = get_lognorm_params(PM$r[1], PM$r[2])
    q1_p = get_lognorm_params(PM$q1[1], PM$q1[2])
    bpsd_p = get_lognorm_params(PM$bpsd[1], PM$bpsd[2])
    bosd_p = get_lognorm_params(PM$bosd[1], PM$bosd[2])
    m0_p = get_lognorm_params(PM$m0[1], PM$m0[2])

    constants = c(base_constants, list(
      sdlog_K = K_p$sdlog, meanlog_K = K_p$meanlog, lower_K = PM$K[1] * 0.01, upper_K = PM$K[1] * 5.0,
      sdlog_r = r_p$sdlog, meanlog_r = r_p$meanlog,  lower_r = 0.1, upper_r = 3.0,  
      sdlog_q1 = q1_p$sdlog, meanlog_q1 = q1_p$meanlog,  lower_q1= 0.1, upper_q1 = 3.0,
      sdlog_bpsd_prior = bpsd_p$sdlog, meanlog_bpsd_prior = bpsd_p$meanlog,
      sdlog_bosd_prior = bosd_p$sdlog, meanlog_bosd_prior = bosd_p$meanlog,
      sdlog_m0 = m0_p$sdlog, meanlog_m0 = m0_p$meanlog,
      mlim_low = PM$mlim[1],
      mlim_high = PM$mlim[2]
    ))

    data = nimble_data
    dimensions = list(
      S_observed = length(data$S_observed),
      removed = length(data$removed)
    )
    inits = initial_values

  } else if (model_type == "delay_difference") {
        
    #' state-space Deriso-Schnute delay-difference model with
    #' recruitment dynamics modeled as an AR(1) latent process:
    #' \deqn{m_t = (1+\rho)e^{-M} m_{t-1} - \rho e^{-2M} m_{t-2} - \frac{(1+\rho)e^{-M}C_{t-1}}{K} + \frac{\rho e^{-2M}C_{t-2}}{K} + \frac{w_k R_t}{K}}

    code = nimbleCode({

      # note: Nimble models require use of "->"

      # Priors on core scaling and variance parameters
      K ~ dlnorm(meanlog = meanlog_K, sdlog = sdlog_K)
      q1 ~ dlnorm(meanlog = meanlog_q1, sdlog = sdlog_q1)
      bpsd ~ dlnorm(meanlog = meanlog_bpsd_prior, sdlog = sdlog_bpsd_prior)
      bosd ~ dlnorm(meanlog = meanlog_bosd_prior, sdlog = sdlog_bosd_prior)

      # Delay difference biological parameters
      M ~ dlnorm(meanlog = meanlog_M, sdlog = sdlog_M)
      rho ~ dlnorm(meanlog = meanlog_rho, sdlog = sdlog_rho)
      wk ~ dlnorm(meanlog = meanlog_wk, sdlog = sdlog_wk)

      # Recruitment process priors
      sigmaR ~ dlnorm(meanlog = meanlog_sigmaR, sdlog = sdlog_sigmaR)
      phi_R ~ dunif(-0.99, 0.99)
      mu_R ~ dnorm(mean = mean_mu_R, sd = sd_mu_R)
      q_R ~ dlnorm(0, sdlog = 1)
      sd_R_obs ~ dlnorm(0, sdlog = 0.5)

      # Initialize relative biomass for the two delay lag years
      m[1] ~ dlnorm(meanlog = meanlog_m0, sdlog = sdlog_m0)
      m[2] ~ dlnorm(meanlog = meanlog_m0, sdlog = sdlog_m0)

      # Initialize recruitment AR(1) state
      logR_process[1] ~ dnorm(mean = mu_R, sd = sigmaR)
      R[1] <- exp(logR_process[1])

      logR_process[2] ~ dnorm(mean = mu_R + phi_R * (logR_process[1] - mu_R), sd = sigmaR)
      R[2] <- exp(logR_process[2])

      # Process model starting from t = 3
      for (i in 3:nT) {
        mean_logR[i] <- mu_R + phi_R * (logR_process[i - 1] - mu_R)
        logR_process[i] ~ dnorm(mean = mean_logR[i], sd = sigmaR)
        R[i] <- exp(logR_process[i])

        mean_m_process[i] <- (1 + rho) * exp(-M) * m[i - 1] -
                            rho * exp(-2 * M) * m[i - 2] -
                            ((1 + rho) * exp(-M) * removed[i - 1]) / K +
                            (rho * exp(-2 * M) * removed[i - 2]) / K +
                            (wk * R[i]) / K

        meanlog_m_process[i] <- log(max(1e-9, mean_m_process[i])) - 0.5 * pow(bpsd, 2)
        m[i] ~ dlnorm(meanlog = meanlog_m_process[i], sdlog = bpsd)
      }

      # Prediction loop into future (no removals)
      if (nM > nT) {
        for (i in (nT + 1):nM) {
          mean_logR[i] <- mu_R + phi_R * (logR_process[i - 1] - mu_R)
          logR_process[i] ~ dnorm(mean = mean_logR[i], sd = sigmaR)
          R[i] <- exp(logR_process[i])

          mean_m_predict[i] <- (1 + rho) * exp(-M) * m[i - 1] -
                              rho * exp(-2 * M) * m[i - 2] +
                              (wk * R[i]) / K

          meanlog_m_predict[i] <- log(max(1e-9, mean_m_predict[i])) - 0.5 * pow(bpsd, 2)
          m[i] ~ dlnorm(meanlog = meanlog_m_predict[i], sdlog = bpsd)
        }
      }

      # Survey observation likelihood
      for (j in 1:n_observed_S) {
        mu_S[j] <- (m[m_idx_for_S[j]] - removed_penalty_for_S_observed[j] / K) / q1
        S_observed[j] ~ dnorm(mu_S[j], sd = bosd)
      }

      # Recruitment index observation likelihood
      for (t in 1:nT) {
        meanlog_R_obs[t] <- log(q_R * R[t]) - 0.5 * pow(sd_R_obs, 2)
        R_observed[t] ~ dlnorm(meanlog = meanlog_R_obs[t], sdlog = sd_R_obs)
      }
    })
    
    buildDerivs = FALSE
    monitors = c("K", "q1", "bpsd", "bosd", "M", "rho", "wk", "sigmaR", "phi_R", "mu_R", "q_R", "m", "R")

    K_p = get_lognorm_params(PM$K[1], PM$K[2])
    q1_p = get_lognorm_params(PM$q1[1], PM$q1[2])
    bpsd_p = get_lognorm_params(PM$bpsd[1], PM$bpsd[2])
    bosd_p = get_lognorm_params(PM$bosd[1], PM$bosd[2])
    M_p = get_lognorm_params(0.2, 0.05)
    rho_p = get_lognorm_params(0.8, 0.1)
    wk_p = get_lognorm_params(0.5, 0.1)
    sigmaR_p = get_lognorm_params(0.3, 0.05)
    m0_p = get_lognorm_params(0.5, 0.1)

    # recruitment index 
    data = nimble_data
    n_obs_years = PM$nT
    mu_true = mean(log(data$R_observed))
 
    constants = c(base_constants, list(
      sdlog_K = K_p$sdlog, meanlog_K = K_p$meanlog,
      sdlog_q1 = q1_p$sdlog, meanlog_q1 = q1_p$meanlog,
      sdlog_bpsd_prior = bpsd_p$sdlog, meanlog_bpsd_prior = bpsd_p$meanlog,
      sdlog_bosd_prior = bosd_p$sdlog, meanlog_bosd_prior = bosd_p$meanlog,
      sdlog_M = M_p$sdlog, meanlog_M = M_p$meanlog,
      sdlog_rho = rho_p$sdlog, meanlog_rho = rho_p$meanlog,
      sdlog_wk = wk_p$sdlog, meanlog_wk = wk_p$meanlog,
      sdlog_sigmaR = sigmaR_p$sdlog, meanlog_sigmaR = sigmaR_p$meanlog,
      sdlog_m0 = m0_p$sdlog, meanlog_m0 = m0_p$meanlog,
      mean_mu_R = mu_true,
      sd_mu_R = 1.0
    ))

    dimensions = list(
      S_observed = length(data$S_observed),
      removed = length(data$removed),
      R_observed = n_obs_years,
      m = PM$nM,
      R = PM$nM,
      logR_process = PM$nM
    )

    inits = list(
      K = PM$K[1],
      q1 = PM$q1[1],
      bpsd = PM$bpsd[1],
      bosd = PM$bosd[1],
      M = 0.2,
      rho = 0.8,
      wk = 0.5,
      sigmaR = 0.3,
      phi_R = 0.5,
      mu_R = mu_true,
      q_R = 1.0,
      sd_R_obs = 0.2,
      m = rep(0.5, PM$nM),
      logR_process = rep(mu_true, PM$nM)
    )
  }

  return(list(
    code = code,
    constants = constants,
    data = data,
    dimensions = dimensions,
    inits = inits,
    monitors = monitors,
    buildDerivs = buildDerivs
  ))
}



#' Sample prior distributions for Nimble biomass dynamics models
#'
#' @description
#' Generates independent Monte Carlo samples directly from prior distributions
#' specified for model parameters.
#'
#' @param setup List returned by \code{prepare_fishery_model_data()}.
#' @param model_type Character string: "lognormal", "gaussian", or "delay_difference".
#' @param n_samples Integer number of prior samples to draw (default 10000).
#' @return A \code{coda::mcmc.list} object containing the prior samples.
sample_priors_nimble = function(
    setup,
    model_type = c("lognormal", "gaussian", "delay_difference"),
    n_samples = 10000
) {
  model_type = match.arg(model_type)
  PM = setup$PM

  if (model_type == "gaussian") {
    rtruncnorm = function(n, mean, sd, lower = 0.01) {
      res = numeric(0)
      while (length(res) < n) {
        s = rnorm(n * 2, mean, sd)
        res = c(res, s[s > lower])
      }
      return(res[1:n])
    }

    prior_mat = cbind(
      K = rtruncnorm(n_samples, PM$K[1], PM$K[2], lower = PM$K[1] * 0.01),
      r = rtruncnorm(n_samples, PM$r[1], PM$r[2], lower = 0.01),
      q1 = rtruncnorm(n_samples, PM$q1[1], PM$q1[2], lower = 0.01),
      bpsd = rtruncnorm(n_samples, PM$bpsd[1], PM$bpsd[2], lower = 0.01),
      bosd = rtruncnorm(n_samples, PM$bosd[1], PM$bosd[2], lower = 0.01)
    )

  } else if (model_type == "lognormal") {
    rtruncrlnorm = function(n, meanlog, sdlog, lower = 0.01) {
      res = numeric(0)
      while (length(res) < n) {
        s = rlnorm(n * 2, meanlog, sdlog)
        res = c(res, s[s > lower])
      }
      return(res[1:n])
    }

    K_p = get_lognorm_params(PM$K[1], PM$K[2])
    r_p = get_lognorm_params(PM$r[1], PM$r[2])
    q1_p = get_lognorm_params(PM$q1[1], PM$q1[2])
    bpsd_p = get_lognorm_params(PM$bpsd[1], PM$bpsd[2])
    bosd_p = get_lognorm_params(PM$bosd[1], PM$bosd[2])

    prior_mat = cbind(
      K = rtruncrlnorm(n_samples, K_p$meanlog, K_p$sdlog, lower = PM$K[1] * 0.01),
      r = rtruncrlnorm(n_samples, r_p$meanlog, r_p$sdlog, lower = 0.01),
      q1 = rtruncrlnorm(n_samples, q1_p$meanlog, q1_p$sdlog, lower = 0.01),
      bpsd = rtruncrlnorm(n_samples, bpsd_p$meanlog, bpsd_p$sdlog, lower = 0.01),
      bosd = rtruncrlnorm(n_samples, bosd_p$meanlog, bosd_p$sdlog, lower = 0.01)
    )

  } else if (model_type == "delay_difference") {
    
    K_p = get_lognorm_params(PM$K[1], PM$K[2])
    q1_p = get_lognorm_params(PM$q1[1], PM$q1[2])
    bpsd_p = get_lognorm_params(PM$bpsd[1], PM$bpsd[2])
    bosd_p = get_lognorm_params(PM$bosd[1], PM$bosd[2])
    M_p = get_lognorm_params(0.2, 0.05)
    rho_p = get_lognorm_params(0.8, 0.1)
    wk_p = get_lognorm_params(0.5, 0.1)
    sigmaR_p = get_lognorm_params(0.3, 0.05)
    mu_true = log(PM$K[1] * 0.05)

    prior_mat = cbind(
      K = rlnorm(n_samples, K_p$meanlog, K_p$sdlog),
      q1 = rlnorm(n_samples, q1_p$meanlog, q1_p$sdlog),
      bpsd = rlnorm(n_samples, bpsd_p$meanlog, bpsd_p$sdlog),
      bosd = rlnorm(n_samples, bosd_p$meanlog, bosd_p$sdlog),
      M = rlnorm(n_samples, M_p$meanlog, M_p$sdlog),
      rho = rlnorm(n_samples, rho_p$meanlog, rho_p$sdlog),
      wk = rlnorm(n_samples, wk_p$meanlog, wk_p$sdlog),
      sigmaR = rlnorm(n_samples, sigmaR_p$meanlog, sigmaR_p$sdlog),
      phi_R = runif(n_samples, -0.99, 0.99),
      mu_R = rnorm(n_samples, mu_true, 1.0),
      q_R = rlnorm(n_samples, 0, 1.0)
    )
  }

  return(coda::mcmc.list(coda::mcmc(prior_mat)))
}


 
#' Compute lognormal parameters from arithmetic mean and standard deviation
#'
#' @description
#' Converts arithmetic mean \eqn{\mu} and standard deviation \eqn{\sigma} to
#' the log-scale mean (\eqn{\mu_{\log}}) and standard deviation (\eqn{\sigma_{\log}})
#' via method-of-moments:
#' \deqn{\sigma_{\log} = \sqrt{\ln\left(1 + \frac{\sigma^2}{\mu^2}\right)}}
#' \deqn{\mu_{\log} = \ln(\mu) - \frac{1}{2}\sigma_{\log}^2}
#'
#' @param mu Numeric arithmetic mean (\eqn{\mu > 0}).
#' @param sd Numeric arithmetic standard deviation (\eqn{\sigma > 0}).
#' @return A list containing \code{meanlog} and \code{sdlog}.
get_lognorm_params = function(mu, sd) {
  if (mu <= 0 || sd <= 0) {
    stop("Mean (mu) and standard deviation (sd) must be strictly positive.")
  }
  sdlog = sqrt(log(1.0 + (sd / mu)^2))
  meanlog = log(mu) - 0.5 * sdlog^2
  return(list(meanlog = meanlog, sdlog = sdlog))
}




#' Posterior predictive check (PPC) for fitted biomass dynamics models
#'
#' @description
#' Simulates replicated survey index observations \eqn{S^{rep}} given parameter
#' draws from the joint posterior distribution, and compares replicated data
#' distributions against observed data \eqn{S^{obs}}.
#'
#' For each draw \eqn{d} the expected index is
#' \deqn{\mu_{j}^{(d)} = \frac{m_{i_j}^{(d)} - C_{pen,j}/K^{(d)}}{q_1^{(d)}}, \qquad
#'   S_j^{rep,(d)} \sim N(\mu_j^{(d)}, bosd^{(d)})}
#' The predicted index is the posterior mean \eqn{\hat{S}_j = E[\mu_j]}.
#' Agreement between \eqn{S^{obs}} and \eqn{\hat{S}} is summarised by:
#' \itemize{
#'   \item Pearson correlation \eqn{r = \mathrm{cor}(S^{obs}, \hat{S})} and \eqn{r^2};
#'   \item the posterior distribution of \eqn{r^{(d)} = \mathrm{cor}(S^{obs}, \mu^{(d)})},
#'     reported as a median and 95\% credible interval;
#'   \item \eqn{RMSE = \sqrt{n^{-1}\sum_j (S_j^{obs} - \hat{S}_j)^2}} and
#'     bias \eqn{= n^{-1}\sum_j (\hat{S}_j - S_j^{obs})}.
#' }
#' Bayesian p-values \eqn{P(T(S^{rep}) \ge T(S^{obs}))} are computed for the
#' mean, SD and maximum (Gelman et al. 2013, ch. 6).
#'
#' @param samples A \code{coda::mcmc.list} object from MCMC sampling.
#' @param setup List returned by \code{prepare_fishery_model_data()}.
#' @param model_type Character string: "lognormal", "gaussian", or "delay_difference".
#' @param n_sims Integer number of posterior predictive draws (default 500).
#' @param output_dir Character path to save PPC table and plot (optional).
#' @param prefix Character filename prefix.
#' @param outformat Character plot file format passed to \code{ggsave} (default "png").
#' @return A list containing \code{S_obs}, \code{S_rep}, \code{mu_S},
#'   \code{summary} (test statistics and Bayesian p-values), \code{fit}
#'   (correlation, \eqn{r^2}, RMSE, bias), \code{r_draws}, \code{predictions}
#'   (per-observation data frame) and \code{plot} (ggplot object).
posterior_predictive_check = function(
    samples,
    setup,
    model_type = c("lognormal", "gaussian", "delay_difference"),
    n_sims = 500,
    output_dir = NULL,
    prefix = "",
    outformat = "png"
) {

  model_type = match.arg(model_type)
  PM = setup$PM
  S_obs = setup$S_obs[PM$iok]
  n_obs = length(S_obs)
  penalty = setup$removed_penalty_for_S_observed

  if (n_obs < 3) {
    stop("posterior_predictive_check: at least 3 observations are required for correlation.")
  }

  mcmc_mat = as.matrix(samples)
  n_total = nrow(mcmc_mat)
  draws = sample(seq_len(n_total), min(n_sims, n_total))
  n_used = length(draws)

  m_cols = paste0("m[", PM$iok, "]")
  missing_cols = setdiff(c("K", "q1", "bosd", m_cols), colnames(mcmc_mat))
  if (length(missing_cols) > 0) {
    stop("posterior_predictive_check: samples lack columns: ",
         paste(missing_cols, collapse = ", "))
  }

  # Expected index (n_obs x n_used) and replicated observations
  K_d = mcmc_mat[draws, "K"]
  q1_d = mcmc_mat[draws, "q1"]
  bosd_d = mcmc_mat[draws, "bosd"]
  m_d = t(mcmc_mat[draws, m_cols, drop = FALSE])

  mu_S = (m_d - outer(penalty, 1 / K_d)) / matrix(q1_d, n_obs, n_used, byrow = TRUE)
  S_rep = mu_S + matrix(rnorm(n_obs * n_used), n_obs, n_used) *
    matrix(bosd_d, n_obs, n_used, byrow = TRUE)

  # Test statistics: Mean, SD, Max
  obs_mean = mean(S_obs)
  obs_sd = sd(S_obs)
  obs_max = max(S_obs)

  rep_means = apply(S_rep, 2, mean)
  rep_sds = apply(S_rep, 2, sd)
  rep_maxs = apply(S_rep, 2, max)

  # Bayesian p-values: P(T(S_rep) >= T(S_obs))
  p_mean = mean(rep_means >= obs_mean)
  p_sd = mean(rep_sds >= obs_sd)
  p_max = mean(rep_maxs >= obs_max)

  ppc_summary = data.frame(
    statistic = c("Mean", "SD", "Max"),
    observed = c(obs_mean, obs_sd, obs_max),
    rep_mean = c(mean(rep_means), mean(rep_sds), mean(rep_maxs)),
    rep_q025 = c(quantile(rep_means, 0.025), quantile(rep_sds, 0.025), quantile(rep_maxs, 0.025)),
    rep_q975 = c(quantile(rep_means, 0.975), quantile(rep_sds, 0.975), quantile(rep_maxs, 0.975)),
    bayesian_p_value = round(c(p_mean, p_sd, p_max), 3)
  )

  # Observed vs predicted: Pearson correlation and error statistics
  S_pred = rowMeans(mu_S)
  pred_lo = apply(S_rep, 1, quantile, probs = 0.025)
  pred_hi = apply(S_rep, 1, quantile, probs = 0.975)

  r_hat = cor(S_obs, S_pred, method = "pearson")
  r_test = cor.test(S_obs, S_pred, method = "pearson")
  r_draws = apply(mu_S, 2, function(mu) cor(S_obs, mu, method = "pearson"))
  r_q = quantile(r_draws, c(0.025, 0.5, 0.975))
  rmse = sqrt(mean((S_obs - S_pred)^2))
  bias = mean(S_pred - S_obs)
  coverage = mean(S_obs >= pred_lo & S_obs <= pred_hi)

  fit = data.frame(
    n = n_obs,
    pearson_r = r_hat,
    r_squared = r_hat^2,
    r_freq_q025 = r_test$conf.int[1],
    r_freq_q975 = r_test$conf.int[2],
    r_post_median = r_q[2],
    r_post_q025 = r_q[1],
    r_post_q975 = r_q[3],
    rmse = rmse,
    bias = bias,
    coverage_95 = coverage,
    row.names = NULL
  )

  cat("\n--- Posterior Predictive Check (PPC) Summary ---\n")
  print(ppc_summary, row.names = FALSE)
  cat("\n--- Observed vs Predicted ---\n")
  print(fit, row.names = FALSE, digits = 3)

  predictions = data.frame(
    year = PM$yrs[PM$iok],
    observed = S_obs,
    predicted = S_pred,
    pred_q025 = pred_lo,
    pred_q975 = pred_hi
  )

  stats_label = paste(
    sprintf("n = %d", n_obs),
    sprintf("Pearson r = %.3f", r_hat),
    sprintf("posterior r: %.3f [%.3f, %.3f]", r_q[2], r_q[1], r_q[3]),
    sprintf("r\u00b2 = %.3f", r_hat^2),
    sprintf("RMSE = %.3f", rmse),
    sprintf("bias = %.3f", bias),
    sprintf("95%% PI coverage = %.2f", coverage),
    sprintf("p(mean) = %.2f, p(SD) = %.2f, p(max) = %.2f", p_mean, p_sd, p_max),
    sep = "\n"
  )

  axis_rng = range(c(predictions$observed, predictions$pred_q025, predictions$pred_q975))

  pl = ggplot2::ggplot(predictions, ggplot2::aes(x = observed, y = predicted)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = pred_q025, ymax = pred_q975),
      width = 0, colour = "steelblue", alpha = 0.4
    ) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                         colour = "firebrick", linewidth = 0.7) +
    ggplot2::geom_point(colour = "steelblue4", size = 2.5) +
    ggplot2::geom_text(ggplot2::aes(label = year), size = 2.5, hjust = -0.3,
                       colour = "grey30") +
    ggplot2::annotate(
      "label", x = axis_rng[1], y = axis_rng[2], label = stats_label,
      hjust = 0, vjust = 1, size = 3.2,
      fill = "white", alpha = 0.85
    ) +
    ggplot2::coord_equal(xlim = axis_rng, ylim = axis_rng) +
    ggplot2::labs(
      x = "Observed survey index",
      y = "Predicted survey index (posterior mean, 95% predictive interval)",
      title = paste0("Posterior predictive check: ", setup$aulab, " (", model_type, ")")
    ) +
    ggplot2::theme_bw(base_size = 12)

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    fn = file.path(output_dir, paste0(prefix, "posterior_predictive_check.csv"))
    write.table(ppc_summary, file = fn, sep = ";", row.names = FALSE)
    fn_fit = file.path(output_dir, paste0(prefix, "posterior_predictive_fit.csv"))
    write.table(fit, file = fn_fit, sep = ";", row.names = FALSE)
    fn_pred = file.path(output_dir, paste0(prefix, "posterior_predictive_predictions.csv"))
    write.table(predictions, file = fn_pred, sep = ";", row.names = FALSE)
    fn_pl = file.path(output_dir, paste0(prefix, "posterior_predictive_xy.", outformat))
    ggplot2::ggsave(fn_pl, pl, width = 7, height = 7, dpi = 150)
    cat(paste("PPC outputs saved to:", output_dir, "\n"))
  }

  return(list(
    S_obs = S_obs,
    S_rep = S_rep,
    mu_S = mu_S,
    summary = ppc_summary,
    fit = fit,
    r_draws = r_draws,
    predictions = predictions,
    plot = pl
  ))
}



namedlist_to_dataframe = function(...) {
  args = list(...)
  names_args = names(args)
  if (is.null(names_args)) {
    stop("Arguments to namedlist_to_dataframe must be named.")
  }

  # Convert named list to a data frame using base expand.grid
  df = do.call(expand.grid, args)
  return(df)
}
 

# fishery_model_predictions
# 'res' is expected to be a coda::mcmc.list object from nimble
# 'nM' should be available in the global environment or passed as an argument
fishery_model_predictions = function(res, n_sample = -1, nM) {
  if (!inherits(res, "mcmc.list")) {
    stop("'res' must be an 'mcmc.list' object from the coda package.")
  }

  nchains = length(res) # Number of chains
  nsims_per_chain = nrow(res[[1]]) # Number of iterations per chain

  if (n_sample == -1) {
    n_sample = nchains * nsims_per_chain
    oo = namedlist_to_dataframe(sims = 1:nsims_per_chain, chains = 1:nchains)
  } else {
    # Randomly sample iterations and chains
    sim_indices = sample(1:nsims_per_chain, n_sample, replace = TRUE)
    chain_indices = sample(1:nchains, n_sample, replace = TRUE)
    oo = data.frame(sims = sim_indices, chains = chain_indices)
  }

  md = matrix(0, nrow = nM, ncol = n_sample)
  mb = matrix(0, nrow = nM, ncol = n_sample)

  for (z in 1:n_sample) {
    j = oo$sims[z]   # iteration index within chain
    l = oo$chains[z] # chain index

    current_chain = res[[l]]

    # Extract K for current iteration and chain
    current_K = current_chain[j, "K"]

    for (i in 1:nM) {
      # Parameter name for m[i]
      m_param_name = paste0("m[", i, "]")
      md[i, z] = current_chain[j, m_param_name]
      mb[i, z] = md[i, z] * current_K
    }
  }

  # Returning a list equivalent to Julia's tuple
  return(list(md, NULL, mb, NULL, NULL, NULL))
}



# logistic_discrete_reference_points
logistic_discrete_reference_points = function(r, K) {
  # K is already on the natural scale in the Nimble model (sampled via dlnorm)
  msy   = r * K / 4.0      # maximum height of of the latent productivity (yield)
  bmsy  = K / 2.0          # biomass at MSY
  fmsy  = 2.0 * msy / K    # fishing mortality at MSY (which simplifies to r / 2.0)
  return(list(msy, bmsy, fmsy))
}

# fishing_mortality_instantaneous
fishing_mortality_instantaneous = function(removed, abundance) {
  return(-log(1.0 - (removed / abundance)))
}



# Expects 'removed' as a vector (length nT or nM),
# 'bio' as a matrix (e.g., nM x n_sample) representing spring biomass,
# and 'survey_time' to define the time dimension (though not directly used for indexing here).
fishery_model_mortality = function(removed, bio) {
  # The Julia code implies 'fb' is a matrix matching the time dimension of 'removed'
  # and the sample dimension of 'bio'.
  # Assuming 'bio' (e.g., mb from fishery_model_predictions) has dimensions (nM x n_sample)
  # and 'removed' is a vector of length 'nT' (or 'nM').
  # We need to broadcast 'removed' across the 'n_sample' dimension.

  # Ensure 'bio' has enough rows for the 'removed' vector
  if (nrow(bio) < length(removed)) {
    stop("Number of rows in 'bio' must be at least the length of 'removed'.")
  }

  fb = bio[1:length(removed), , drop = FALSE] # Select rows of bio corresponding to 'removed' length

  # Create a matrix of 'removed' values by replicating it across columns
  Fkt = matrix(removed, nrow = nrow(fb), ncol = ncol(fb), byrow = FALSE)

  HR = Fkt / fb   # exploitation rate = removals / spring biomass (element-wise)

  # Check for analytical anomaly where removals exceed standing spring biomass
  if (any(HR > 1.0, na.rm = TRUE)) {
    warning("Harvest rate HR > 1.0 detected: removals exceed estimated spring standing biomass.")
  }

  # Instantaneous F: -log(1 - HR), safeguarding against log(0) or negative arguments
  FM = -1.0 * log(1.0 - pmin(HR, 0.99999))

  # Returning a list equivalent to Julia's tuple
  return(list(Fkt, HR, FM))
}


# 'res' is expected to be a coda::mcmc.list object
# 'Sai' is expected to be a vector of observed survey indices (length nT)
# Note: with the relative-scale likelihood S ≈ m/q1, the estimated biomass is
# B = m*K = S*q1*K. This formula is therefore unchanged, but q1 is now
# dimensionless (O(1)) rather than in kt. Catch penalty is ignored here
# (approximation for the survey overlay plot only).
abundance_from_index = function(Sai, res) {
  mcmc_matrix = as.matrix(res)
  K_samples = mcmc_matrix[, "K"]
  q1_samples = mcmc_matrix[, "q1"]

  # S_m[i, j] = Sai[i] * q1_samples[j] * K_samples[j]
  S_m = outer(Sai, K_samples * q1_samples, FUN = "*")

  return(S_m)
}

#' Compute catch penalty for observed survey biomass
#'
#' @description
#' Calculates the catch deduction vector corresponding to observed survey
#' times, accounting for the spring-to-fall survey timing transition in 2004.
#' In spring-surveyed years (before yeartransition), surveys occur prior to the
#' main fishery removals, so the penalty is 0. In the transition year, the
#' penalty is the average of current and previous year removals:
#'   \deqn{C_{pen} = \frac{C_t + C_{t-1}}{2}}
#' In fall-surveyed years, the survey occurs after fishery removals, so the
#' full removal is subtracted:
#'   \deqn{C_{pen} = C_t}
#' When yeartransition == 0 (e.g., CFA 4X), surveys are consistently in fall,
#' so the full removal is used for all observations.
#'
#' @param PM A list containing model parameters and data:
#'   \itemize{
#'     \item \code{iok}: Integer vector of indices (1 to nT) where survey observations exist.
#'     \item \code{removed}: Numeric vector of fishery removals of length nT.
#'     \item \code{yeartransition}: Integer year index of survey timing transition (0 if none).
#'   }
#' @return Numeric vector of length \code{length(PM$iok)} containing catch penalties.
compute_removed_penalty = function(PM) {
  n_obs = length(PM$iok)
  penalty = numeric(n_obs)

  if (is.null(PM$yeartransition) || PM$yeartransition == 0) {
    penalty = as.numeric(PM$removed[PM$iok])
  } else {
    for (j in seq_len(n_obs)) {
      i = PM$iok[j]
      if (i < PM$yeartransition) {
        penalty[j] = 0.0
      } else if (i == PM$yeartransition) {
        prev_rem = if (i > 1) PM$removed[i - 1] else 0.0
        penalty[j] = (prev_rem + PM$removed[i]) / 2.0
      } else {
        penalty[j] = PM$removed[i]
      }
    }
  }

  return(as.numeric(penalty))
}

#' Compute equilibrium reference points for Deriso-Schnute delay-difference model
#'
#' @description
#' Calculates equilibrium biomass, yield, and fisheries management reference
#' points for the Deriso-Schnute delay-difference model. In equilibrium, relative
#' biomass \eqn{m_{eq}(F)} under instantaneous fishing mortality \eqn{F} satisfies:
#' \deqn{m_{eq}(F) = \frac{w_k R_{eq} / K}{1 - (1+\rho)e^{-(M+F)} + \rho e^{-(2M+F)}}}
#' with equilibrium yield \eqn{Y_{eq}(F) = (1 - e^{-F}) m_{eq}(F) K}.
#' Reference points (\eqn{F_{0.1}}, \eqn{B_{msy}}, \eqn{MSY}) are computed by
#' solving for the fishing mortality where marginal yield equals 10% of the
#' slope at the origin:
#' \deqn{\left.\frac{\partial Y_{eq}}{\partial F}\right|_{F=F_{0.1}} = 0.1 \times \left.\frac{\partial Y_{eq}}{\partial F}\right|_{F=0}}
#'
#' @param M Numeric natural mortality rate (\eqn{M > 0}).
#' @param rho Numeric Brody growth coefficient (\eqn{0 < \rho < 1}).
#' @param wk Numeric weight of recruits (\eqn{w_k > 0}).
#' @param K Numeric carrying capacity scaling parameter (\eqn{K > 0}).
#' @param Req Numeric mean recruitment (\eqn{R_{eq} > 0}).
#' @return A list containing \code{msy}, \code{bmsy}, \code{fmsy}, \code{b0}, and \code{msy_hr}.
delay_difference_reference_points = function(M, rho, wk, K, Req) {
  n_draws = max(length(M), length(rho), length(wk), length(K), length(Req))
  if (length(M) < n_draws) M = rep_len(M, n_draws)
  if (length(rho) < n_draws) rho = rep_len(rho, n_draws)
  if (length(wk) < n_draws) wk = rep_len(wk, n_draws)
  if (length(K) < n_draws) K = rep_len(K, n_draws)
  if (length(Req) < n_draws) Req = rep_len(Req, n_draws)

  calc_single = function(m_val, rho_val, wk_val, k_val, req_val) {
    if (any(is.na(c(m_val, rho_val, wk_val, k_val, req_val)))) {
      return(c(msy = NA_real_, bmsy = NA_real_, fmsy = NA_real_, b0 = NA_real_))
    }
    denom0 = 1.0 - (1.0 + rho_val) * exp(-m_val) + rho_val * exp(-2.0 * m_val)
    if (denom0 <= 0) {
      return(c(msy = NA_real_, bmsy = NA_real_, fmsy = NA_real_, b0 = NA_real_))
    }
    b0 = (wk_val * req_val) / denom0

    y_eq = function(f_val) {
      d = 1.0 - (1.0 + rho_val) * exp(-(m_val + f_val)) + rho_val * exp(-(2.0 * m_val + f_val))
      if (d <= 0) return(0.0)
      (1.0 - exp(-f_val)) * (wk_val * req_val) / d
    }

    h = 1e-5
    s0 = (y_eq(h) - y_eq(0.0)) / h
    if (s0 <= 0) {
      return(c(msy = NA_real_, bmsy = NA_real_, fmsy = NA_real_, b0 = b0))
    }

    target = function(f_val) {
      dy = (y_eq(f_val + h) - y_eq(f_val - h)) / (2.0 * h)
      dy - 0.1 * s0
    }

    f01 = tryCatch({
      uniroot(target, c(1e-4, 3.0))$root
    }, error = function(e) {
      m_val
    })

    msy_val = y_eq(f01)
    bmsy_val = msy_val / (1.0 - exp(-f01))

    return(c(msy = msy_val, bmsy = bmsy_val, fmsy = f01, b0 = b0))
  }

  res_mat = matrix(NA_real_, nrow = n_draws, ncol = 4)
  colnames(res_mat) = c("msy", "bmsy", "fmsy", "b0")

  for (i in seq_len(n_draws)) {
    res_mat[i, ] = calc_single(M[i], rho[i], wk[i], K[i], Req[i])
  }

  msy = res_mat[, "msy"]
  bmsy = res_mat[, "bmsy"]
  fmsy = res_mat[, "fmsy"]
  b0 = res_mat[, "b0"]
  msy_hr = 1.0 - exp(-fmsy)

  return(list(msy = msy, bmsy = bmsy, fmsy = fmsy, b0 = b0, msy_hr = msy_hr))
}

#' Convergence diagnostics for MCMC samples
#'
#' @description
#' Evaluates MCMC chain convergence and mixing efficiency using Gelman-Rubin
#' potential scale reduction factor (\eqn{\hat{R}}), effective sample size (ESS),
#' and parameter trace inspections. Parameters exceeding convergence thresholds
#' are flagged with informative warnings.
#'
#' @param samples A \code{coda::mcmc.list} object containing MCMC chains.
#' @param threshold_rhat Numeric threshold for Gelman-Rubin \eqn{\hat{R}} (default 1.01).
#' @param threshold_ess Numeric threshold for minimum effective sample size (default 400).
#' @param params Character vector of parameter names to evaluate. If NULL, evaluates all non-vector parameters.
#' @param output_dir Character path to directory for saving diagnostic summaries. If NULL, does not save to disk.
#' @param prefix Character prefix for output filenames (default "").
#' @return A data frame containing parameter names, Point est. \eqn{\hat{R}}, Upper CI \eqn{\hat{R}}, ESS, and convergence flag.
compute_convergence_diagnostics = function(
    samples,
    threshold_rhat = 1.01,
    threshold_ess = 400,
    params = NULL,
    output_dir = NULL,
    prefix = ""
) {
  if (!inherits(samples, "mcmc.list")) {
    stop("'samples' must be an 'mcmc.list' object from coda.")
  }

  all_varnames = coda::varnames(samples)
  if (is.null(params)) {
    params = all_varnames[!grepl("\\[", all_varnames)]
    if (length(params) == 0) params = all_varnames
  } else {
    params = intersect(params, all_varnames)
  }

  if (length(params) == 0) {
    warning("No matching parameters found for convergence diagnostics.")
    return(NULL)
  }

  samples_sub = samples[, params, drop = FALSE]
  n_chains = length(samples_sub)

  rhat_point = rep(NA_real_, length(params))
  rhat_upper = rep(NA_real_, length(params))

  if (n_chains >= 2) {
    tryCatch({
      gd = coda::gelman.diag(samples_sub, autoburnin = FALSE, multivariate = FALSE)
      psrf = gd$psrf
      rhat_point = psrf[, "Point est."]
      rhat_upper = psrf[, "Upper C.I."]
    }, error = function(e) {
      warning(paste("Gelman-Rubin calculation encountered an issue:", e$message))
      for (k in seq_along(params)) {
        tryCatch({
          single_gd = coda::gelman.diag(samples_sub[, params[k], drop = FALSE], autoburnin = FALSE)
          rhat_point[k] = single_gd$psrf[1, "Point est."]
          rhat_upper[k] = single_gd$psrf[1, "Upper C.I."]
        }, error = function(err) NULL)
      }
    })
  } else {
    warning("Gelman-Rubin diagnostic requires at least 2 chains; skipping R-hat.")
  }

  ess = rep(NA_real_, length(params))
  tryCatch({
    ess = coda::effectiveSize(samples_sub)
  }, error = function(e) {
    warning(paste("Effective sample size calculation error:", e$message))
  })

  diag_df = data.frame(
    parameter = params,
    rhat = round(as.numeric(rhat_point), 3),
    rhat_upper_ci = round(as.numeric(rhat_upper), 3),
    ess = round(as.numeric(ess), 1),
    converged = TRUE,
    stringsAsFactors = FALSE
  )

  if (n_chains >= 2) {
    rhat_failed = which(diag_df$rhat > threshold_rhat | is.na(diag_df$rhat))
    diag_df$converged[rhat_failed] = FALSE
  }
  ess_failed = which(diag_df$ess < threshold_ess | is.na(diag_df$ess))
  diag_df$converged[ess_failed] = FALSE

  cat("\n--- MCMC Convergence Diagnostics Summary ---\n")
  print(diag_df, row.names = FALSE)

  failed_params = diag_df$parameter[!diag_df$converged]
  if (length(failed_params) > 0) {
    cat(sprintf("\n[WARNING] %d parameter(s) failed convergence criteria (R-hat > %.2f or ESS < %d):\n",
                length(failed_params), threshold_rhat, threshold_ess))
    cat(paste("  -", failed_params, collapse = "\n"), "\n")
  } else {
    cat(sprintf("\n[OK] All evaluated parameters meet convergence criteria (R-hat <= %.2f, ESS >= %d).\n",
                threshold_rhat, threshold_ess))
  }

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    fn_csv = file.path(output_dir, paste0(prefix, "convergence_diagnostics.csv"))
    write.table(diag_df, file = fn_csv, sep = ";", row.names = FALSE)
    cat(paste("Diagnostics saved to:", fn_csv, "\n"))
  }

  return(diag_df)
}


# probability_pa
# 'res' is coda::mcmc.list, 'bio' is matrix (nM x n_sample),
# 'HR' and 'FM' are matrices, 'yr' is scalar.
probability_pa = function(res, bio, HR, FM, yr, PM, yrs, prediction_time_ss) {
  # Find index for the given year
  yri = which(PM$yrs == yr)
  if (length(yri) == 0) stop(paste("Year", yr, "not found in PM$yrs."))

  # Extract spring biomass for the relevant year
  # bio is assumed to be nM x n_sample, rows correspond to time points in PM$yrs
  fb_spring_vec = as.vector(bio[yri, ])

  # Get removed value for the specific year 'yr' and subtract from spring biomass
  removed_for_yr = PM$removed[yri]
  if (is.na(removed_for_yr)) {
    removed_for_yr = 0
  }
  fb_fall_vec = fb_spring_vec - removed_for_yr

  # Extract K and r samples from the MCMC results
  mcmc_matrix = as.matrix(res)
  K_samples = mcmc_matrix[, "K"]

  if ("r" %in% colnames(mcmc_matrix)) {
    r_samples = mcmc_matrix[, "r"]
    ref_points = logistic_discrete_reference_points(r_samples, K_samples)
    msy = unlist(ref_points[[1]])
    bmsy = unlist(ref_points[[2]])
    fmsy = unlist(ref_points[[3]])
    msy_hr = 1 - exp(-fmsy)
  } else if (all(c("M", "rho", "wk") %in% colnames(mcmc_matrix))) {
    M_s = mcmc_matrix[, "M"]
    rho_s = mcmc_matrix[, "rho"]
    wk_s = mcmc_matrix[, "wk"]
    Req_s = if ("mu_R" %in% colnames(mcmc_matrix)) exp(mcmc_matrix[, "mu_R"]) else K_samples * 0.05
    ref_points = delay_difference_reference_points(M_s, rho_s, wk_s, K_samples, Req_s)
    msy = ref_points$msy
    bmsy = ref_points$bmsy
    fmsy = ref_points$fmsy
    msy_hr = ref_points$msy_hr
    r_samples = rep(NA_real_, length(K_samples))
  } else {
    r_samples = rep(NA_real_, length(K_samples))
    msy = rep(NA_real_, length(K_samples))
    bmsy = rep(NA_real_, length(K_samples))
    fmsy = rep(NA_real_, length(K_samples))
    msy_hr = rep(NA_real_, length(K_samples))
  }

  ndat = length(K_samples)
  cls = numeric(3) # Initialize counts for zones: critical, cautious, healthy

  # Classify fall biomass into zones based on LRP (bmsy/2) and USR (bmsy)
  n_valid = 0
  for (i in 1:ndat) {
    if (is.na(fb_fall_vec[i])) next
    lrp = if (!is.na(bmsy[i])) bmsy[i] / 2.0 else (if (!is.na(K_samples[i])) K_samples[i] / 4.0 else NA_real_)
    usr = if (!is.na(bmsy[i])) bmsy[i] else (if (!is.na(K_samples[i])) K_samples[i] / 2.0 else NA_real_)
    if (is.na(lrp) || is.na(usr)) next

    n_valid = n_valid + 1
    if (fb_fall_vec[i] <= lrp) {
      cls[1] = cls[1] + 1 # Critical
    } else if (fb_fall_vec[i] <= usr) {
      cls[2] = cls[2] + 1 # Cautious
    } else {
      cls[3] = cls[3] + 1 # Healthy
    }
  }

  probs = if (n_valid > 0) round(cls / n_valid, digits = 3) else rep(NA_real_, 3) # Probabilities for each zone

  # Helper function to calculate mean, sd, and quantiles
  calc_stats = function(x) {
    round(c(mean(x, na.rm = TRUE), sd(x, na.rm = TRUE), quantile(x, 0.025, na.rm = TRUE), quantile(x, 0.975, na.rm = TRUE)), digits = 3)
  }

  # Calculate statistics for various variables
  fb_fall_stats = calc_stats(fb_fall_vec)
  fb_spring_stats = calc_stats(fb_spring_vec)

  # Extract FM and HR for the specific year and flatten to a vector
  if (yri <= nrow(FM)) {
    fm_vec = as.vector(FM[yri, ])  # FM is a matrix (time x samples), select row 'yri'
    fm_stats = calc_stats(fm_vec)

    hr_vec = as.vector(HR[yri, ]) # HR is a matrix (time x samples), select row 'yri'
    hr_stats = calc_stats(hr_vec)

    hrM_vec = 1 - exp(-fm_vec)  # HR derived from FM
    hrM_stats = calc_stats(hrM_vec)
  } else {
    fm_stats = rep(NA, 4)
    hr_stats = rep(NA, 4)
    hrM_stats = rep(NA, 4)
  }

  msy_hr_stats = calc_stats(msy_hr)
  fmsy_stats = calc_stats(fmsy)
  r_stats = calc_stats(r_samples) # r_samples is already a vector of samples

  # Define row names for the data frame
  variable_names = c("mean", "sd", "q0.025", "q0.975")
  zone_names = c("critical", "cautious", "healthy", NA) # NA for the 4th row of the dataframe for zone

  # Create and return the data frame
  toreturn = data.frame(
    variable = variable_names,
    fb_fall = fb_fall_stats,
    fb_spring = fb_spring_stats,
    fm = fm_stats,
    hr = hr_stats,
    hrM = hrM_stats,
    fmsy = fmsy_stats,
    msy_hr = msy_hr_stats,
    r = r_stats,
    zone = zone_names,
    probability_zone = c(probs, NA) # Combine probabilities and NA for the 4th zone
  )

  return(toreturn)
}



plot_prior_posterior = function(vn, prior_samples, posterior_samples, bw = 0.02) {
  prior_mat = as.matrix(prior_samples)
  posterior_mat = as.matrix(posterior_samples)

  has_prior = vn %in% colnames(prior_mat)

  if (has_prior) {
    prior_vec = as.vector(prior_mat[, vn])
    posterior_vec = as.vector(posterior_mat[, vn])

    # Combine into a data frame for ggplot
    plot_df = data.frame(
      value = c(prior_vec, posterior_vec),
      Type = factor(c(rep("Prior", length(prior_vec)), rep("Posterior", length(posterior_vec))))
    )

    # Create the density plot using ggplot2
    pl = ggplot(plot_df, aes(x = value, fill = Type, color = Type)) +
      geom_density(alpha = 0.5, bw = bw, linewidth = 0) +
      scale_fill_manual(values = c("Prior" = "slateblue", "Posterior" = "purple")) +
      scale_color_manual(values = c("Prior" = "slateblue", "Posterior" = "purple")) +
      labs(
        x = vn,
        y = "Density",
        title = paste("Prior vs. Posterior Density for", vn)
      ) +
      theme_bw() +
      theme(
        legend.position = "top",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14)
      )
  } else {
    warning(paste("Variable", vn, "not found in prior_samples. Plotting posterior only."))
    posterior_vec = as.vector(posterior_mat[, vn])

    plot_df = data.frame(
      value = posterior_vec,
      Type = factor(rep("Posterior", length(posterior_vec)))
    )

    pl = ggplot(plot_df, aes(x = value, fill = Type, color = Type)) +
      geom_density(alpha = 0.5, bw = bw, linewidth = 0) +
      scale_fill_manual(values = c("Posterior" = "purple")) +
      scale_color_manual(values = c("Posterior" = "purple")) +
      labs(
        x = vn,
        y = "Density",
        title = paste("Posterior Density for", vn)
      ) +
      theme_bw() +
      theme(
        legend.position = "top",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14)
      )
  }

  return(pl)
}

fishery_model_plot = function(
    toplot = c("fishing", "survey"),
    n_sample = 250,
    res,
    bio, # bio is assumed to be the 'mb' (unscaled biomass) from fishery_model_predictions
    FM,
    HR,
    S, # original survey data (vector)
    PM, # PM list with removed, etc.
    prediction_time_ss,
    fall_time_ss,
    survey_time,
    yrs,
    alphav = 0.075,
    pl_in = NULL,
    time_range = c(floor(min(prediction_time_ss, na.rm=TRUE)) - 1.0, ceiling(max(prediction_time_ss, na.rm=TRUE)) + 0.5)
) {

  nsims = ncol(bio) # Number of simulations/samples in bio matrix
  # Adjust n_sample if it's greater than nsims
  n_sample_actual = min(n_sample, nsims)
  ss = sample(1:nsims, n_sample_actual, replace = FALSE) # sample indices for plotting individual traces

  # Initialize ggplot object if not provided
  if (is.null(pl_in)) {
    pl = ggplot() + theme_bw() + labs(y = "", x = "Year")
  } else {
    pl = pl_in
  }

  if ("trace" %in% toplot) {
    warning("trace is not valid for a discrete model")
  }

  if ("nofishing" %in% toplot) {
    warning("nofishing not implemented")
  }

  # --- Plot Spring Biomass (with fishing) ---
  if ("fishing" %in% toplot) {
    g = bio # [yr, sim]

    # Prepare data for individual simulation lines
    df_sims_fishing = as.data.frame(g[, ss])
    colnames(df_sims_fishing) = paste0("sim", 1:n_sample_actual)
    df_sims_fishing$time = prediction_time_ss
    df_sims_fishing_long = df_sims_fishing %>% pivot_longer(cols = starts_with("sim"), names_to = "sim_id", values_to = "biomass")

    # Prepare data for mean line
    mean_g = rowMeans(g, na.rm = TRUE)
    df_mean_fishing = data.frame(time = prediction_time_ss, biomass = mean_g)

    pl = pl +
      geom_line(data = df_sims_fishing_long, aes(x = time, y = biomass, group = sim_id), alpha = alphav, color = "orange") +
      geom_line(data = df_mean_fishing, aes(x = time, y = biomass), color = "darkorange", linewidth = 1.5) +
      labs(title = "Spring Biomass (with Fishing)") +
      coord_cartesian(
        xlim = time_range,
        ylim = c(0, quantile(as.vector(g), 0.99, na.rm = TRUE) * 1.05)
      ) +
      theme(legend.position = "none", axis.text = element_text(size = 14))
  }

  # --- Plot Fall Biomass (with fishing) ---
  if ("fall" %in% toplot) {
    # PM$removed is aligned to PM$yrs (nM rows); bio has been pre-subsetted
    # to length(prediction_time_ss) rows by the caller. Use tail() to align.
    removed_plot = tail(PM$removed, nrow(bio))

    if (length(removed_plot) != nrow(bio)) {
      stop("Length of removed_plot must match number of rows in bio for 'fall' plot.")
    }
    g_fall = bio - matrix(removed_plot, nrow = nrow(bio), ncol = ncol(bio), byrow = FALSE)

    # Prepare data for individual simulation lines
    df_sims_fall = as.data.frame(g_fall[, ss])
    colnames(df_sims_fall) = paste0("sim", 1:n_sample_actual)
    df_sims_fall$time = fall_time_ss # Julia used fall_time_ss here
    df_sims_fall_long = df_sims_fall %>% pivot_longer(cols = starts_with("sim"), names_to = "sim_id", values_to = "biomass")

    # Prepare data for mean line
    mean_g_fall = rowMeans(g_fall, na.rm = TRUE)
    df_mean_fall = data.frame(time = fall_time_ss, biomass = mean_g_fall)

    pl = pl +
      geom_line(data = df_sims_fall_long, aes(x = time, y = biomass, group = sim_id), alpha = alphav, color = "orange") +
      geom_line(data = df_mean_fall, aes(x = time, y = biomass), color = "darkorange", linewidth = 1.5) +
      labs(title = "Fall Biomass (with Fishing)") +
      coord_cartesian(
        xlim = time_range,
        ylim = c(0, quantile(as.vector(g_fall), 0.99, na.rm = TRUE) * 1.05)
      ) +
      theme(legend.position = "none", axis.text = element_text(size = 14))
  }

  if ("footprint" %in% toplot) {
    warning("footprint not implemented")
  }

  # --- Plot Survey Index ---
  if ("survey" %in% toplot) {
    # S_m = abundance_from_index(S, res) # returns a matrix (length(S) x n_total_samples)
    # S_K = rowMeans(S_m, na.rm = TRUE) # average by year

    # Assuming S_m and S_K are already calculated or can be calculated from global context
    # If S_m and S_K are not globally available, you might need to pass them or calculate here.
    S_m_val = abundance_from_index(S, res) # Call the R translation
    S_K = rowMeans(S_m_val, na.rm = TRUE)

    df_survey = data.frame(time = survey_time, S_K = S_K)

    pl = pl +
      geom_line(data = df_survey, aes(x = time, y = S_K), color = "gray", linewidth = 1) +
      geom_point(data = df_survey, aes(x = time, y = S_K), color = "darkgray", size = 2) +
      labs(title = "Survey Index") +
      coord_cartesian(xlim = time_range) +
      theme(legend.position = "none", axis.text = element_text(size = 14))
  }

  # --- Plot Fishing Mortality ---
  if ("fishing_mortality" %in% toplot) {
    FM_mean = rowMeans(FM, na.rm = TRUE)
    FM_mean[is.na(FM_mean)] = 0 # Julia's FMmean[isnan.(FMmean)] .= zero(eltype(FM))
    ub = quantile(as.vector(FM), 0.99, na.rm = TRUE) * 1.05 # Julia uses quantile on flattened array

    df_sims_fm = as.data.frame(FM[, ss])
    colnames(df_sims_fm) = paste0("sim", 1:n_sample_actual)
    df_sims_fm$time = prediction_time_ss
    df_sims_fm_long = df_sims_fm %>% pivot_longer(cols = starts_with("sim"), names_to = "sim_id", values_to = "fm")

    df_mean_fm = data.frame(time = prediction_time_ss, fm = FM_mean)

    pl = pl +
      geom_line(data = df_sims_fm_long, aes(x = time, y = fm, group = sim_id), alpha = 0.02, color = "lightslateblue") +
      geom_line(data = df_mean_fm, aes(x = time, y = fm), color = "slateblue", linewidth = 1.5) +
      labs(title = "Fishing Mortality") +
      coord_cartesian(xlim = time_range, ylim = c(0, ub)) +
      theme(legend.position = "none", axis.text = element_text(size = 14))
  }

  # --- Plot Exploitation Rate ---
  if ("exploitation_rate" %in% toplot) {
    HR_mean = rowMeans(HR, na.rm = TRUE)
    HR_mean[is.na(HR_mean)] = 0 # Julia's HRmean[isnan.(HRmean)] .= zero(eltype(HR))
    ub = quantile(as.vector(HR), 0.99, na.rm = TRUE) * 1.05

    df_sims_hr = as.data.frame(HR[, ss])
    colnames(df_sims_hr) = paste0("sim", 1:n_sample_actual)
    df_sims_hr$time = prediction_time_ss
    df_sims_hr_long = df_sims_hr %>% pivot_longer(cols = starts_with("sim"), names_to = "sim_id", values_to = "hr")

    df_mean_hr = data.frame(time = prediction_time_ss, hr = HR_mean)

    pl = pl +
      geom_line(data = df_sims_hr_long, aes(x = time, y = hr, group = sim_id), alpha = 0.02, color = "lightslateblue") +
      geom_line(data = df_mean_hr, aes(x = time, y = hr), color = "slateblue", linewidth = 1.5) +
      labs(title = "Exploitation Rate") +
      coord_cartesian(xlim = time_range, ylim = c(0, ub)) +
      theme(legend.position = "none", axis.text = element_text(size = 14))
  }

  if ("fishing_mortality_vs_footprint" %in% toplot) {
    warning("footprint not implemented")
  }

  if ("harvest_control_rule_footprint" %in% toplot) {
    warning("footprint not implemented")
  }

  # --- Harvest Control Rule Plots ---
  # These are complex as they plot multiple lines/points with reference lines.
  # They also share a lot of logic, so I'll try to generalize.

  # Extract parameter samples from MCMC results
  mcmc_matrix = as.matrix(res)
  K_samples = mcmc_matrix[, "K"]
  q1_samples = mcmc_matrix[, "q1"]

  # Reference points calculation for logistic or delay-difference models
  if ("r" %in% colnames(mcmc_matrix)) {
    r_samples = mcmc_matrix[, "r"]
    ref_points = logistic_discrete_reference_points(r_samples, K_samples)
    fmsy_F_ref = unlist(ref_points[[3]])
    fmsy_hr_ref = 1.0 - exp(-fmsy_F_ref)
    bmsy_ref = unlist(ref_points[[2]])
  } else if (all(c("M", "rho", "wk") %in% colnames(mcmc_matrix))) {
    M_s = mcmc_matrix[, "M"]
    rho_s = mcmc_matrix[, "rho"]
    wk_s = mcmc_matrix[, "wk"]
    Req_s = if ("mu_R" %in% colnames(mcmc_matrix)) exp(mcmc_matrix[, "mu_R"]) else K_samples * 0.05
    ref_points = delay_difference_reference_points(M_s, rho_s, wk_s, K_samples, Req_s)
    fmsy_F_ref = ref_points$fmsy
    fmsy_hr_ref = ref_points$msy_hr
    bmsy_ref = ref_points$bmsy
  } else {
    fmsy_F_ref = rep(NA_real_, length(K_samples))
    fmsy_hr_ref = rep(NA_real_, length(K_samples))
    bmsy_ref = K_samples / 2.0
  }
  lrp_ref = bmsy_ref / 2.0

  if (any(c("harvest_control_rule", "harvest_control_rule_F",
            "harvest_control_rule_fall", "harvest_control_rule_fall_F") %in% toplot)) {

    # Common HCR elements with model-appropriate reference points
    df_hcr_ref = data.frame(
      fmsy_hr = fmsy_hr_ref,
      fmsy_F = fmsy_F_ref,
      K = K_samples,
      bmsy = bmsy_ref,
      lrp = lrp_ref
    )

    # Base ggplot for HCR plots
    pl_hcr_base = ggplot() + theme_bw() +
      theme(legend.position = "none", axis.text = element_text(size = 14)) +

      # Posterior draw reference lines
      geom_vline(data = df_hcr_ref[ss, ], aes(xintercept = K), alpha = 0.05, color = "limegreen") +
      geom_vline(data = df_hcr_ref[ss, ], aes(xintercept = bmsy), alpha = 0.05, color = "darkkhaki") +
      geom_vline(data = df_hcr_ref[ss, ], aes(xintercept = lrp), alpha = 0.05, color = "darkred") +

      # Carrying capacity K
      geom_vline(xintercept = mean(df_hcr_ref$K, na.rm = TRUE), color = "chartreuse4", linewidth = 1.5, alpha = 0.6) +
      geom_vline(xintercept = quantile(df_hcr_ref$K, 0.975, na.rm = TRUE), color = "chartreuse4", linewidth = 0.8, linetype = "dashed", alpha = 0.5) +
      geom_vline(xintercept = quantile(df_hcr_ref$K, 0.025, na.rm = TRUE), color = "chartreuse4", linewidth = 0.8, linetype = "dashed", alpha = 0.5) +

      # Upper Stock Reference (B_msy)
      geom_vline(xintercept = mean(df_hcr_ref$bmsy, na.rm = TRUE), color = "darkkhaki", linewidth = 1.5, alpha = 0.6) +
      geom_vline(xintercept = quantile(df_hcr_ref$bmsy, 0.975, na.rm = TRUE), color = "darkkhaki", linewidth = 0.8, linetype = "dashed", alpha = 0.5) +
      geom_vline(xintercept = quantile(df_hcr_ref$bmsy, 0.025, na.rm = TRUE), color = "darkkhaki", linewidth = 0.8, linetype = "dashed", alpha = 0.5) +

      # Limit Reference Point (LRP = B_msy / 2)
      geom_vline(xintercept = mean(df_hcr_ref$lrp, na.rm = TRUE), color = "darkred", linewidth = 1.5, alpha = 0.6) +
      geom_vline(xintercept = quantile(df_hcr_ref$lrp, 0.975, na.rm = TRUE), color = "darkred", linewidth = 0.8, linetype = "dashed", alpha = 0.5) +
      geom_vline(xintercept = quantile(df_hcr_ref$lrp, 0.025, na.rm = TRUE), color = "darkred", linewidth = 0.8, linetype = "dashed", alpha = 0.5)


    # Prepare common time series annotations (mimic Julia's `colours` and `text.()`)
    nt = length(prediction_time_ss)
    # Julia's `ColorSchemes.tab20c` can be approximated by `scale_color_viridis_d` or similar
    # or a custom palette. For simplicity, let's use a default ggplot color mapping for now.
    df_hcr_time_labels = data.frame(time = prediction_time_ss, label = as.character(trunc(prediction_time_ss)))


    if ("harvest_control_rule" %in% toplot) { # HR-based, Spring Biomass
      fb_hcr = bio[1:length(prediction_time_ss), ] # Spring biomass
      fb_mean = apply(fb_hcr, 1, median, na.rm = TRUE)
      hr_mean = apply(HR[1:length(prediction_time_ss), ], 1, median, na.rm = TRUE)

      # Data for last  HR vs. fb (points)
      fbbb = c(quantile(fb_hcr[nt,], 0.025, na.rm = TRUE), quantile(fb_hcr[nt,], 0.975, na.rm = TRUE))
      HRbb = c(quantile(HR[nt,], 0.025, na.rm = TRUE), quantile(HR[nt,], 0.975, na.rm = TRUE))
      df_hcr_last_year_dist = data.frame(fb = fb_hcr[nt,], hr = HR[nt,])
      df_hcr_last_year_quantiles = data.frame(
        fb_mean = fb_mean[nt], hr_mean = hr_mean[nt],
        fb_low = fbbb[1], fb_high = fbbb[2],
        hr_low = HRbb[1], hr_high = HRbb[2]
      )

      df_hcr_ts_mean = data.frame(fb = fb_mean, hr = hr_mean, time = prediction_time_ss)

      ub = max(quantile(K_samples, 0.9, na.rm = TRUE), quantile(fb_mean, 0.975, na.rm = TRUE)) # xlim upper bound
      uby = max(quantile(fmsy_hr_ref, 0.975, na.rm = TRUE), quantile(as.vector(HR), 0.975, na.rm = TRUE)) # ylim upper bound

      pl_hcr = pl_hcr_base +
        geom_hline(data = df_hcr_ref[ss, ], aes(yintercept = fmsy_hr), alpha = 0.01, color = "lightgray") + # Specific hline for HR-based
        geom_hline(yintercept = mean(df_hcr_ref$fmsy_hr, na.rm = TRUE), color = "darkgray", linewidth = 1.5, alpha = 0.6) +
        geom_hline(yintercept = quantile(df_hcr_ref$fmsy_hr, 0.975, na.rm = TRUE), color = "gray", linewidth = 0.8, linetype = "dashed", alpha = 0.5) +
        geom_hline(yintercept = quantile(df_hcr_ref$fmsy_hr, 0.025, na.rm = TRUE), color = "gray", linewidth = 0.8, linetype = "dashed", alpha = 0.5) +

        geom_point(data = df_hcr_last_year_dist, aes(x = fb, y = hr), alpha = 0.01, color = "magenta", size = 1.5) +
        geom_segment(data = df_hcr_last_year_quantiles, aes(x = fb_low, xend = fb_high, y = hr_mean, yend = hr_mean), color = "magenta", linewidth = 1) +
        geom_segment(data = df_hcr_last_year_quantiles, aes(x = fb_mean, xend = fb_mean, y = hr_low, yend = hr_high), color = "magenta", linewidth = 1) +
        geom_point(data = df_hcr_ts_mean[nt,], aes(x = fb, y = hr), alpha = 0.9, color = "gold", size = 4, stroke = 1) +
        geom_path(data = df_hcr_ts_mean, aes(x = fb, y = hr), color = "slateblue", linewidth = 1) +
        geom_point(data = df_hcr_ts_mean, aes(x = fb, y = hr, color = factor(time)), size = 2) +
        geom_text(data = df_hcr_ts_mean, aes(x = fb + 0.051, y = hr - 0.0025, label = as.character(trunc(time)), color = factor(time)), size = 2.5, hjust = 0, vjust = 0.5) +
        scale_color_viridis_d(option = "A", direction = 1) + # Use a discrete color scale for time points
        labs(x = "Biomass (Spring)", y = "Exploitation Rate") +
        coord_cartesian(xlim = c(0, ub), ylim = c(0, uby))

      pl = pl_hcr # Update main plot with HCR plot
    }

    if ("harvest_control_rule_F" %in% toplot) { # F-based, Spring Biomass
      fb_hcr = bio[1:length(prediction_time_ss), ] # Spring biomass
      fb_mean = apply(fb_hcr, 1, median, na.rm = TRUE)
      fm_mean = apply(FM[1:length(prediction_time_ss), ], 1, median, na.rm = TRUE)

      # Data for last  FM vs. fb (points)
      fbbb = c(quantile(fb_hcr[nt,], 0.025, na.rm = TRUE), quantile(fb_hcr[nt,], 0.975, na.rm = TRUE))
      FMbb = c(quantile(FM[nt,], 0.025, na.rm = TRUE), quantile(FM[nt,], 0.975, na.rm = TRUE))
      df_hcr_last_year_dist_F = data.frame(fb = fb_hcr[nt,], fm = FM[nt,])
      df_hcr_last_year_quantiles_F = data.frame(
        fb_mean = fb_mean[nt], fm_mean = fm_mean[nt],
        fb_low = fbbb[1], fb_high = fbbb[2],
        fm_low = FMbb[1], fm_high = FMbb[2]
      )

      df_hcr_ts_mean_F = data.frame(fb = fb_mean, fm = fm_mean, time = prediction_time_ss)

      ub = max(quantile(K_samples, 0.9, na.rm = TRUE), quantile(fb_mean, 0.975, na.rm = TRUE))
      uby = max(quantile(fmsy_F_ref, 0.9, na.rm = TRUE), quantile(as.vector(FM), 0.95, na.rm = TRUE))

      pl_hcr_F = pl_hcr_base + # Start with the base HCR plot
        geom_hline(data = df_hcr_ref[ss, ], aes(yintercept = fmsy_F), alpha = 0.01, color = "lightgray") + # for F-based
        geom_hline(yintercept = mean(df_hcr_ref$fmsy_F, na.rm = TRUE), color = "darkgray", linewidth = 1.5, alpha = 0.6) +
        geom_hline(yintercept = quantile(df_hcr_ref$fmsy_F, 0.975, na.rm = TRUE), color = "gray", linewidth = 0.8, linetype = "dashed", alpha = 0.5) +
        geom_hline(yintercept = quantile(df_hcr_ref$fmsy_F, 0.025, na.rm = TRUE), color = "gray", linewidth = 0.8, linetype = "dashed", alpha = 0.5) +
        geom_point(data = df_hcr_last_year_dist_F, aes(x = fb, y = fm), alpha = 0.01, color = "magenta", size = 1.5) +
        geom_segment(data = df_hcr_last_year_quantiles_F, aes(x = fb_low, xend = fb_high, y = fm_mean, yend = fm_mean), color = "magenta", linewidth = 1) +
        geom_segment(data = df_hcr_last_year_quantiles_F, aes(x = fb_mean, xend = fb_mean, y = fm_low, yend = fm_high), color = "magenta", linewidth = 1) +
        geom_point(data = df_hcr_ts_mean_F[nt,], aes(x = fb, y = fm), alpha = 0.9, color = "gold", size = 4, stroke = 1) +
        geom_path(data = df_hcr_ts_mean_F, aes(x = fb, y = fm), color = "slateblue", linewidth = 1) +
        geom_point(data = df_hcr_ts_mean_F, aes(x = fb, y = fm, color = factor(time)), size = 2) +
        geom_text(data = df_hcr_ts_mean_F, aes(x = fb + 0.051, y = fm - 0.0025, label = as.character(trunc(time)), color = factor(time)), size = 2.5, hjust = 0, vjust = 0.5) +
        scale_color_viridis_d(option = "A", direction = 1) +
        labs(x = "Biomass (Spring)", y = "Fishing Mortality") +
        coord_cartesian(xlim = c(0, ub), ylim = c(0, uby))
      pl = pl_hcr_F # Update main plot with HCR_F plot
    }

    if ("harvest_control_rule_fall" %in% toplot) { # HR-based, Fall Biomass
      removed_plot = tail(PM$removed, length(prediction_time_ss))
      fb_hcr_fall = bio[1:length(prediction_time_ss), ] - matrix(removed_plot, nrow = length(prediction_time_ss), ncol = nsims, byrow = FALSE)
      fb_mean_fall = apply(fb_hcr_fall, 1, median, na.rm = TRUE)
      hr_mean_fall = apply(HR[1:length(prediction_time_ss), ], 1, median, na.rm = TRUE)

      fbbb_fall = c(quantile(fb_hcr_fall[nt,], 0.025, na.rm = TRUE), quantile(fb_hcr_fall[nt,], 0.975, na.rm = TRUE))
      HRbb_fall = c(quantile(HR[nt,], 0.025, na.rm = TRUE), quantile(HR[nt,], 0.975, na.rm = TRUE))
      df_hcr_last_year_dist_fall = data.frame(fb = fb_hcr_fall[nt,], hr = HR[nt,])
      df_hcr_last_year_quantiles_fall = data.frame(
        fb_mean = fb_mean_fall[nt], hr_mean = hr_mean_fall[nt],
        fb_low = fbbb_fall[1], fb_high = fbbb_fall[2],
        hr_low = HRbb_fall[1], hr_high = HRbb_fall[2]
      )

      df_hcr_ts_mean_fall = data.frame(fb = fb_mean_fall, hr = hr_mean_fall, time = prediction_time_ss)

      ub = max(quantile(K_samples, 0.9, na.rm = TRUE), quantile(fb_mean_fall, 0.975, na.rm = TRUE))
      uby = max(quantile(fmsy_hr_ref, 0.9, na.rm = TRUE), quantile(as.vector(HR), 0.95, na.rm = TRUE))

      pl_hcr_fall = pl_hcr_base +
        geom_hline(data = df_hcr_ref[ss, ], aes(yintercept = fmsy_hr), alpha = 0.01, color = "lightgray") + # Specific hline for HR-based
        geom_hline(yintercept = mean(df_hcr_ref$fmsy_hr, na.rm = TRUE), color = "darkgray", linewidth = 1.5, alpha = 0.6) +
        geom_hline(yintercept = quantile(df_hcr_ref$fmsy_hr, 0.975, na.rm = TRUE), color = "gray", linewidth = 0.8, linetype = "dashed", alpha = 0.5) +
        geom_hline(yintercept = quantile(df_hcr_ref$fmsy_hr, 0.025, na.rm = TRUE), color = "gray", linewidth = 0.8, linetype = "dashed", alpha = 0.5) +

        geom_point(data = df_hcr_last_year_dist_fall, aes(x = fb, y = hr), alpha = 0.01, color = "magenta", size = 1.5) +
        geom_segment(data = df_hcr_last_year_quantiles_fall, aes(x = fb_low, xend = fb_high, y = hr_mean, yend = hr_mean), color = "magenta", linewidth = 1) +
        geom_segment(data = df_hcr_last_year_quantiles_fall, aes(x = fb_mean, xend = fb_mean, y = hr_low, yend = hr_high), color = "magenta", linewidth = 1) +
        geom_point(data = df_hcr_ts_mean_fall[nt,], aes(x = fb, y = hr), alpha = 0.9, color = "gold", size = 4, stroke = 1) +
        geom_path(data = df_hcr_ts_mean_fall, aes(x = fb, y = hr), color = "slateblue", linewidth = 1) +
        geom_point(data = df_hcr_ts_mean_fall, aes(x = fb, y = hr, color = factor(time)), size = 2) +
        geom_text(data = df_hcr_ts_mean_fall, aes(x = fb + 0.051, y = hr - 0.0025, label = as.character(trunc(time)), color = factor(time)), size = 2.5, hjust = 0, vjust = 0.5) +
        scale_color_viridis_d(option = "A", direction = 1) +
        labs(x = "Biomass (Fall)", y = "Exploitation Rate") +
        coord_cartesian(xlim = c(0, ub), ylim = c(0, uby))
      pl = pl_hcr_fall # Update main plot with HCR_fall plot
    }

    if ("harvest_control_rule_fall_F" %in% toplot) { # F-based, Fall Biomass
      removed_plot_fall_F = tail(PM$removed, length(prediction_time_ss))
      fb_hcr_fall = bio[1:length(prediction_time_ss), ] - matrix(removed_plot_fall_F, nrow = length(prediction_time_ss), ncol = nsims, byrow = FALSE)
      fb_mean_fall = apply(fb_hcr_fall, 1, median, na.rm = TRUE)
      fm_mean_fall = apply(FM[1:length(prediction_time_ss), ], 1, median, na.rm = TRUE)

      fbbb_fall = c(quantile(fb_hcr_fall[nt,], 0.025, na.rm = TRUE), quantile(fb_hcr_fall[nt,], 0.975, na.rm = TRUE))
      FMbb_fall = c(quantile(FM[nt,], 0.025, na.rm = TRUE), quantile(FM[nt,], 0.975, na.rm = TRUE))
      df_hcr_last_year_dist_fall_F = data.frame(fb = fb_hcr_fall[nt,], fm = FM[nt,])
      df_hcr_last_year_quantiles_fall_F = data.frame(
        fb_mean = fb_mean_fall[nt], fm_mean = fm_mean_fall[nt],
        fb_low = fbbb_fall[1], fb_high = fbbb_fall[2],
        fm_low = FMbb_fall[1], fm_high = FMbb_fall[2]
      )

      df_hcr_ts_mean_fall_F = data.frame(fb = fb_mean_fall, fm = fm_mean_fall, time = prediction_time_ss)

      ub = max(quantile(K_samples, 0.9, na.rm = TRUE), quantile(fb_mean_fall, 0.975, na.rm = TRUE))
      uby = max(quantile(fmsy_F_ref, 0.9, na.rm = TRUE), quantile(as.vector(FM), 0.95, na.rm = TRUE))

      pl_hcr_fall_F = pl_hcr_base +
        geom_point(data = df_hcr_last_year_dist_fall_F, aes(x = fb, y = fm), alpha = 0.01, color = "magenta", size = 1.5) +
        geom_segment(data = df_hcr_last_year_quantiles_fall_F, aes(x = fb_low, xend = fb_high, y = fm_mean, yend = fm_mean), color = "magenta", linewidth = 1) +
        geom_segment(data = df_hcr_last_year_quantiles_fall_F, aes(x = fb_mean, xend = fb_mean, y = fm_low, yend = fm_high), color = "magenta", linewidth = 1) +
        geom_point(data = df_hcr_ts_mean_fall_F[nt,], aes(x = fb, y = fm), alpha = 0.9, color = "gold", size = 4, stroke = 1) +
        geom_path(data = df_hcr_ts_mean_fall_F, aes(x = fb, y = fm), color = "slateblue", linewidth = 1) +
        geom_point(data = df_hcr_ts_mean_fall_F, aes(x = fb, y = fm, color = factor(time)), size = 2) +
        geom_text(data = df_hcr_ts_mean_fall_F, aes(x = fb + 0.051, y = fm - 0.0025, label = as.character(trunc(time)), color = factor(time)), size = 2.5, hjust = 0, vjust = 0.5) +
        scale_color_viridis_d(option = "A", direction = 1) +
        labs(x = "Biomass (Fall)", y = "Fishing Mortality") +
        coord_cartesian(xlim = c(0, ub), ylim = c(0, uby))
      pl = pl_hcr_fall_F # Update main plot with HCR_fall_F plot
    }
  }

  # --- Prediction Space ---
  if ("prediction_space" %in% toplot) {
    K_pred = K_samples[ss]
    q1_pred = q1_samples[ss]

    # bio has been pre-subsetted to length(prediction_time_ss) rows by the caller.
    # Use nrow(bio) for indexing to stay consistent with whatever subset was passed.
    t1_idx = 2:nrow(bio)
    t0_idx = 1:(nrow(bio) - 1)

    # b1 is the predicted survey index: (B / K) / q1
    b1 = (bio[t1_idx, ss] / matrix(K_pred, nrow = length(t1_idx), ncol = n_sample_actual, byrow = TRUE)) / matrix(q1_pred, nrow = length(t1_idx), ncol = n_sample_actual, byrow = TRUE)

    # b0 is the observed survey index S (already on the index scale)
    b0 = matrix(S[t0_idx], nrow = length(t0_idx), ncol = n_sample_actual, byrow = FALSE)

    df_pred_space = data.frame(b0 = as.vector(b0), b1 = as.vector(b1))
    pl = pl +
      geom_point(data = df_pred_space, aes(x = b0, y = b1), alpha = 0.2, color = "darkgray", size = 1.5) +
      labs(x = expression(paste("Observed Index S(t)")), y = expression(paste("Predicted Index S(t+1)")), title = "Prediction Space") +
      theme(legend.position = "none", axis.text = element_text(size = 14))
  }

  # --- State Space ---
  if ("state_space" %in% toplot) {
    K_state = K_samples[ss]
    q1_state = q1_samples[ss]

    t1_idx = 2:nrow(bio)
    t0_idx = 1:(nrow(bio) - 1)

    # b_scaled is the predicted survey index: (B / K) / q1
    b_scaled = (bio[, ss] / matrix(K_state, nrow = nrow(bio), ncol = n_sample_actual, byrow = TRUE)) / matrix(q1_state, nrow = nrow(bio), ncol = n_sample_actual, byrow = TRUE)
    b1_state = b_scaled[t1_idx, ]
    b0_state = b_scaled[t0_idx, ]

    df_state_space = data.frame(b0 = as.vector(b0_state), b1 = as.vector(b1_state))
    pl = pl +
      geom_point(data = df_state_space, aes(x = b0, y = b1), alpha = 0.2, color = "darkgray", size = 1.5) +
      labs(x = expression(paste("Predicted Index S(t)")), y = expression(paste("Predicted Index S(t+1)")), title = "State Space") +
      theme(legend.position = "none", axis.text = element_text(size = 14))
  }

  # --- Surplus Production ---
  if ("surplus_production" %in% toplot) {
    # r_samples is set above only when "r" is in the MCMC output columns.
    # Use a flag to avoid relying on exists() which can pick up stale session objects.
    has_r = "r" %in% colnames(mcmc_matrix)
    if (!has_r) {
      warning("'r' not found in MCMC samples. Skipping surplus_production plot (only valid for logistic models).")
    } else {
      K_sp = K_samples[ss]
      r_sp = r_samples[ss]
      q1_sp = q1_samples[ss]

      # b_scaled_sp is the predicted survey index = (B/K) * q1 = m * q1,
      # consistent with Julia's b = bio[:, ss] ./ K' .* q1'
      b_scaled_sp = (bio[, ss] / matrix(K_sp, nrow = nrow(bio), ncol = n_sample_actual, byrow = TRUE)) * matrix(q1_sp, nrow = nrow(bio), ncol = n_sample_actual, byrow = TRUE)
      sp = matrix(r_sp, nrow = nrow(b_scaled_sp), ncol = n_sample_actual, byrow = TRUE) * b_scaled_sp * (1.0 - b_scaled_sp)

      # Flatten for ggplot
      df_surplus_prod = data.frame(
        biomass = as.vector(b_scaled_sp),
        surplus_prod = as.vector(sp),
        time = rep(prediction_time_ss, times = n_sample_actual) # Aligned with column-major matrix unrolling
      )

      pl = pl +
        geom_point(data = df_surplus_prod, aes(x = biomass, y = surplus_prod, color = factor(time)), alpha = 0.3, size = 1.5) +
        labs(x = "Biomass", y = "Surplus Production", title = "Surplus Production") +
        scale_color_viridis_d(option = "A", direction = 1, name = "Year") +
        theme(legend.position = "right", axis.text = element_text(size = 14))
    }
  }

  return(pl)
}


do_retrospective_analysis = function(
    nimble_code,
    PM_base,              # The base PM list containing survey indices, years, etc.
    nimble_constants_base,# The base flattened constants list for nimbleModel
    model_data_base,      # The list of data (containing S_observed and removed)
    initial_values_base,  # The list of initial values
    monitor_params,       # Vector of parameter names to monitor in MCMC
    truncate_max_years,   # Maximum number of years to truncate (analyzes from 1 to this max)
    rho_target_param,     # Parameter name (string) for Mohn's rho calculation (e.g., "K")
    mcmc_settings = list(niter = 2000, nburnin = 1000, nchains = 2, progressBar = TRUE, samplesAsCodaMCMC = TRUE),
    buildDerivs = TRUE
) {

  cat("\n--- Starting Retrospective Analysis ---\n")

  # 1. Run the Full Model to get the reference estimate
  cat("\nRunning full model for reference estimates...\n")
  nimble_model_full = nimbleModel(
    code = nimble_code,
    constants = nimble_constants_base,
    data = model_data_base,
    dimensions = lapply(model_data_base, length),
    inits = initial_values_base,
    buildDerivs = buildDerivs,
    check = TRUE
  )
  C_nimble_model_full = compileNimble(nimble_model_full)
  mcmc_conf_full = configureMCMC(nimble_model_full, monitors = unique(c(monitor_params, "m")))
  mcmc_builder_full = buildMCMC(mcmc_conf_full)
  C_mcmc_builder_full = compileNimble(mcmc_builder_full, project = nimble_model_full)
  full_samples = runMCMC(C_mcmc_builder_full,
                          niter = mcmc_settings$niter,
                          nburnin = mcmc_settings$nburnin,
                          nchains = mcmc_settings$nchains,
                          progressBar = mcmc_settings$progressBar,
                          samplesAsCodaMCMC = mcmc_settings$samplesAsCodaMCMC)

  # Extract the reference estimates (mean and 95% CI) for all monitored parameters
  full_param_means = sapply(monitor_params, function(p) mean(as.matrix(full_samples)[, p]))
  full_param_quantiles = sapply(monitor_params, function(p) quantile(as.matrix(full_samples)[, p], probs = c(0.025, 0.975)))
  full_param_lower_ci = full_param_quantiles[1, ]
  full_param_upper_ci = full_param_quantiles[2, ]
  cat(paste0("Full model mean of ", rho_target_param, ": ", round(full_param_means[rho_target_param], 3), "\n"))

  # Store retrospective results for all monitored parameters
  retrospective_samples_list = list()
  retrospective_param_means_matrix = matrix(NA, nrow = truncate_max_years, ncol = length(monitor_params),
                                             dimnames = list(NULL, monitor_params))
  retrospective_param_lower_ci_matrix = matrix(NA, nrow = truncate_max_years, ncol = length(monitor_params),
                                                dimnames = list(NULL, monitor_params))
  retrospective_param_upper_ci_matrix = matrix(NA, nrow = truncate_max_years, ncol = length(monitor_params),
                                                dimnames = list(NULL, monitor_params))

  # 2. Loop for truncating years
  for (k in 1:truncate_max_years) {
    cat(paste0("\nRunning retrospective model: truncating ", k, " year(s) of data...\n"))

    # Deep copy base constants and data
    modified_PM = PM_base
    modified_data = model_data_base
    modified_initial_values = initial_values_base
    modified_constants = nimble_constants_base

    # Modify PM for retrospective analysis
    current_nT = PM_base$nT - k
    if (current_nT <= 0) {
      warning(paste0("Cannot truncate ", k, " years: not enough data. Stopping retrospective analysis."))
      break
    }

    modified_PM$yrs = PM_base$yrs[1:current_nT]
    modified_PM$removed = PM_base$removed[1:current_nT]
    modified_PM$S = PM_base$S[1:current_nT]
    modified_PM$nT = current_nT
    modified_PM$nM = PM_base$nP + current_nT
    modified_PM$iok = which(!is.na(modified_PM$S))

    # Adjust yeartransition if it falls outside the new data range
    if (modified_PM$yeartransition > current_nT && modified_PM$yeartransition != 0) {
      modified_PM$yeartransition = 0
      cat(paste0("  Warning: PM$yeartransition (", PM_base$yeartransition, ") is outside new data range, resetting to 0.\n"))
    }

    # Rebuild flattened nimble constants
    modified_constants$nT = as.integer(modified_PM$nT)
    modified_constants$nM = as.integer(modified_PM$nM)
    modified_constants$m_idx_for_S = as.integer(modified_PM$iok)
    modified_constants$n_observed_S = length(modified_PM$iok)

    removed_penalty_for_S_observed = compute_removed_penalty(modified_PM)
    modified_constants$removed_penalty_for_S_observed = as.numeric(removed_penalty_for_S_observed)

    # Update model_data
    modified_data$S_observed = modified_PM$S[modified_PM$iok]
    modified_data$removed = modified_PM$removed

    # Regenerate catch-feasibility constraint vectors when the model uses them
    feasibility_names = c("feasible_proc", "feasible_fall")
    if (any(feasibility_names %in% names(modified_data))) {
      fd = feasibility_constraint_data(modified_PM$nT, length(modified_PM$iok))
      modified_data[feasibility_names] = fd[feasibility_names]
    }

    # Truncate any additional nT-length data vectors (e.g. R_observed for
    # the delay-difference model). If a data vector has length == original
    # nT, truncate it to current_nT.
    original_nT = PM_base$nT
    for (dname in names(modified_data)) {
      if (dname %in% c("S_observed", "removed", feasibility_names)) next
      if (length(modified_data[[dname]]) == original_nT) {
        modified_data[[dname]] = model_data_base[[dname]][1:current_nT]
      }
    }

    # Update initial_values for 'm' to match new nM
    # 'm' is relative biomass (proportion of K), so initialize to 0.5
    modified_initial_values$m = rep(0.5, modified_PM$nM)

    # Truncate logR_process and R initial values if present (delay-difference)
    if (!is.null(modified_initial_values$logR_process)) {
      modified_initial_values$logR_process = rep(
        modified_initial_values$logR_process[1], modified_PM$nM
      )
    }

    # Build dimensions dynamically from modified_data
    modified_dims = lapply(modified_data, length)

    # Run Retrospective Model
    retrospective_nimble_model = nimbleModel(
      code = nimble_code,
      constants = modified_constants,
      data = modified_data,
      dimensions = modified_dims,
      inits = modified_initial_values,
      buildDerivs = buildDerivs,
      check = TRUE
    )
    C_retrospective_nimble_model = compileNimble(retrospective_nimble_model)

    retrospective_mcmc_conf = configureMCMC(retrospective_nimble_model, monitors = unique(c(monitor_params, "m")))
    retrospective_mcmc_builder = buildMCMC(retrospective_mcmc_conf)
    C_retrospective_mcmc_builder = compileNimble(retrospective_mcmc_builder, project = retrospective_nimble_model) # Use compiled model

    retrospective_samples = runMCMC(C_retrospective_mcmc_builder,
                                     niter = mcmc_settings$niter,
                                     nburnin = mcmc_settings$nburnin,
                                     nchains = mcmc_settings$nchains,
                                     progressBar = mcmc_settings$progressBar,
                                     samplesAsCodaMCMC = mcmc_settings$samplesAsCodaMCMC)

    retrospective_samples_list[[k]] = retrospective_samples

    current_means = sapply(monitor_params, function(p) mean(as.matrix(retrospective_samples)[, p]))
    current_quantiles = sapply(monitor_params, function(p) quantile(as.matrix(retrospective_samples)[, p], probs = c(0.025, 0.975)))

    retrospective_param_means_matrix[k, ] = current_means
    retrospective_param_lower_ci_matrix[k, ] = current_quantiles[1, ]
    retrospective_param_upper_ci_matrix[k, ] = current_quantiles[2, ]

    cat(paste0("  Retrospective model mean of ", rho_target_param, " (", k, " years truncated): ", round(retrospective_param_means_matrix[k, rho_target_param], 3), "\n"))
  }

  # 3. Calculate Mohn's Rho for all monitored parameters
  cat("\n--- Calculating Mohn's Rho ---\n")
  # Remove any NAs if the loop broke early
  retrospective_param_means_matrix = retrospective_param_means_matrix[1:length(retrospective_samples_list), , drop = FALSE]
  retrospective_param_lower_ci_matrix = retrospective_param_lower_ci_matrix[1:length(retrospective_samples_list), , drop = FALSE]
  retrospective_param_upper_ci_matrix = retrospective_param_upper_ci_matrix[1:length(retrospective_samples_list), , drop = FALSE]

  mohn_rho_all = sapply(monitor_params, function(p) {
    if (nrow(retrospective_param_means_matrix) > 0) {
      current_full_mean = full_param_means[p]
      current_retrospective_means = retrospective_param_means_matrix[, p]
      rho_values_per_param = (current_retrospective_means - current_full_mean) / current_full_mean
      mean(rho_values_per_param, na.rm = TRUE)
    } else {
      NA
    }
  })

  if (length(mohn_rho_all) > 0) {
    cat(paste0("Mohn's Rho for ", rho_target_param, ": ", round(mohn_rho_all[rho_target_param], 4), "\n"))
  } else {
    cat("Not enough retrospective runs to calculate Mohn's Rho.\n")
  }

  cat("\n--- Retrospective Analysis Complete ---\n")
  return(list(
    full_model_samples = full_samples,
    retrospective_model_samples = retrospective_samples_list,
    full_param_means = full_param_means,
    full_param_lower_ci = full_param_lower_ci, # New
    full_param_upper_ci = full_param_upper_ci, # New
    retrospective_param_means = retrospective_param_means_matrix,
    retrospective_param_lower_ci = retrospective_param_lower_ci_matrix, # New
    retrospective_param_upper_ci = retrospective_param_upper_ci_matrix, # New
    mohn_rho = mohn_rho_all
  ))
}



post_process_and_save_fishery_model = function(
    samples, setup, prior = NULL, do_plots = TRUE, do_hcr_probabilities = TRUE, outformat = "png",
    outputs_dir = setup$outputs_dir,
    output_prefix = "results_nimble",
    run_diagnostics = TRUE
) {
  # Ensure outputs directory exists
  dir.create(outputs_dir, recursive = TRUE, showWarnings = FALSE)

  # Extract setup variables
  PM = setup$PM
  aulab = setup$aulab
  year_assessment = setup$year_assessment
  prediction_time_ss = setup$prediction_time_ss
  fall_time_ss = setup$fall_time_ss
  survey_time = setup$survey_time
  yrs = setup$yrs
  S_obs = setup$S_obs

  # Pad PM$removed to length of yrs with 0s (assume 0 catch for future/unobserved years to compute probabilities)
  if (length(PM$removed) < length(yrs)) {
    PM$removed = c(PM$removed, rep(0, length(yrs) - length(PM$removed)))
  }
  PM$removed[is.na(PM$removed)] = 0

  cat("\n\n--- Starting Post-Processing ---\n")

  # 1. Define File Paths for Output
  res_fn = file.path(outputs_dir, paste0(output_prefix, "_", aulab, ".rds"))
  prior_fn = file.path(outputs_dir, paste0(output_prefix, "_", aulab, "_prior", ".rds"))
  summary_fn = file.path(outputs_dir, paste0(output_prefix, "_", aulab, "_summary", ".csv"))
  bio_fn1 = file.path(outputs_dir, paste0(output_prefix, "_", aulab, "_bio_fishing", ".csv"))
  fm_fn = file.path(outputs_dir, paste0(output_prefix, "_", aulab, "_fm", ".csv"))
  fn_ypa = file.path(outputs_dir, paste0(output_prefix, "_", aulab, "_hcr_probabilities", ".csv"))

  # 2. Save Raw MCMC Outputs
  saveRDS(samples, file = res_fn)
  cat(paste0("Posterior samples saved to: ", res_fn, "\n"))

  if (!is.null(prior) && inherits(prior, "mcmc.list")) {
    saveRDS(prior, file = prior_fn)
    cat(paste0("Prior samples saved to: ", prior_fn, "\n"))
  } else {
    cat("No 'prior' object provided; skipping prior save.\n")
  }

  # Run convergence diagnostics
  if (run_diagnostics) {
    cat("\nRunning convergence diagnostics...\n")
    compute_convergence_diagnostics(
      samples = samples,
      output_dir = outputs_dir,
      prefix = paste0(aulab, "_")
    )
  }

  # 3. Core Calculations (Predictions and Mortality)
  cat("\nCalculating predictions and mortality metrics...\n")
  prediction_results = fishery_model_predictions(res = samples, nM = PM$nM)
  md_R = prediction_results[[1]]
  mb_R = prediction_results[[3]]

  fm_mortality_results = fishery_model_mortality(removed = PM$removed, bio = mb_R)
  HR_R = fm_mortality_results[[2]]
  FM_R = fm_mortality_results[[3]]

  # Ensure FM and HR match nrow(mb_R) length by padding with NAs
  if (nrow(FM_R) < nrow(mb_R)) {
    pad_rows = nrow(mb_R) - nrow(FM_R)
    FM_R = rbind(FM_R, matrix(NA, nrow=pad_rows, ncol=ncol(FM_R)))
    HR_R = rbind(HR_R, matrix(NA, nrow=pad_rows, ncol=ncol(HR_R)))
  }

  # 4. Save Summaries and Derived Data to CSV
  cat("\nSaving summary and model outputs to CSV files...\n")

  # Summary of Posterior Samples
  summary_df = as.data.frame(summary(samples)$statistics)
  write.table(summary_df, file = summary_fn, sep = ";", row.names = TRUE)
  cat(paste0("Summary saved to: ", summary_fn, "\n"))

  # Compute year-range indices once (used for both CSV export and plotting).
  # Selects the last length(prediction_time_ss) rows of mb_R, aligning
  # row indices to the years in prediction_time_ss (1999 onward).
  idx_plot = (nrow(mb_R) - length(prediction_time_ss) + 1):nrow(mb_R)

  # Fall Biomass (bio_post) — use the same year-range slice as the plots.
  mb_R_subset = mb_R[idx_plot, , drop = FALSE]
  removed_subset = PM$removed[idx_plot]
  bio_post = mb_R_subset - matrix(removed_subset, nrow = nrow(mb_R_subset), ncol = ncol(mb_R_subset), byrow = FALSE)
  write.table(as.data.frame(bio_post), file = bio_fn1, sep = ";", row.names = FALSE)
  cat(paste0("Fall biomass (bio_post) saved to: ", bio_fn1, "\n"))

  # Fishing Mortality (FM)
  write.table(as.data.frame(FM_R), file = fm_fn, sep = ";", row.names = FALSE)
  cat(paste0("Fishing mortality (FM) saved to: ", fm_fn, "\n"))

  # HCR Probabilities
  if (do_hcr_probabilities) {
    # Pass the full mb_R so that yri (which indexes into PM$yrs) aligns with bio rows.
    ypa = probability_pa(
      res = samples, bio = mb_R, HR = HR_R, FM = FM_R,
      yr = year_assessment,
      PM = PM, yrs = yrs, prediction_time_ss = prediction_time_ss
    )
    write.table(ypa, file = fn_ypa, sep = ";", row.names = FALSE)
    cat(paste0("HCR probabilities saved to: ", fn_ypa, "\n"))

    cat(paste0("\nProbabilities by Zone: ", aulab, " (year: ", year_assessment, ")\n"))
    print(ypa)

    ypa1 = probability_pa(
      res = samples, bio = mb_R, HR = HR_R, FM = FM_R,
      yr = year_assessment - 1,
      PM = PM, yrs = yrs, prediction_time_ss = prediction_time_ss
    )
    cat(paste0("\nProbabilities by Zone: ", aulab, " (year: ", year_assessment - 1, ")\n"))
    print(ypa1)
  }

  # 5. Visualizations
  if (do_plots) {
    cat(paste0("\nGenerating and saving plots to: ", outputs_dir, "\n"))

    # Use idx_plot (computed above) for plot subsetting
    mb_R_plot = mb_R[idx_plot, , drop = FALSE]
    FM_R_plot = FM_R[idx_plot, , drop = FALSE]
    HR_R_plot = HR_R[idx_plot, , drop = FALSE]

    # Biomass Trajectories
    pl_fishing = fishery_model_plot(toplot=c("survey", "fishing"), res=samples, bio=mb_R_plot, FM=FM_R_plot, HR=HR_R_plot, S=S_obs, PM=PM, prediction_time_ss=prediction_time_ss, fall_time_ss=fall_time_ss, survey_time=survey_time, yrs=yrs)
    ggsave(filename=file.path(outputs_dir, paste0("plot_predictions_", aulab, ".", outformat)), plot=pl_fishing, device=outformat)

    pl_fall = fishery_model_plot(toplot=c("survey", "fall"), res=samples, bio=mb_R_plot, FM=FM_R_plot, HR=HR_R_plot, S=S_obs, PM=PM, prediction_time_ss=prediction_time_ss, fall_time_ss=fall_time_ss, survey_time=survey_time, yrs=yrs)
    ggsave(filename=file.path(outputs_dir, paste0("plot_predictions_fall_", aulab, ".", outformat)), plot=pl_fall, device=outformat)

    # Mortality & Exploitation
    pl_fm = fishery_model_plot(toplot="fishing_mortality", res=samples, bio=mb_R_plot, FM=FM_R_plot, HR=HR_R_plot, S=S_obs, PM=PM, prediction_time_ss=prediction_time_ss, fall_time_ss=fall_time_ss, survey_time=survey_time, yrs=yrs)
    ggsave(filename=file.path(outputs_dir, paste0("plot_fishing_mortality_", aulab, ".", outformat)), plot=pl_fm, device=outformat)

    pl_hr = fishery_model_plot(toplot="exploitation_rate", res=samples, bio=mb_R_plot, FM=FM_R_plot, HR=HR_R_plot, S=S_obs, PM=PM, prediction_time_ss=prediction_time_ss, fall_time_ss=fall_time_ss, survey_time=survey_time, yrs=yrs)
    ggsave(filename=file.path(outputs_dir, paste0("plot_exploitation_rate_", aulab, ".", outformat)), plot=pl_hr, device=outformat)

    # Harvest Control Rules
    pl_hcr = fishery_model_plot(toplot="harvest_control_rule", n_sample=1000, res=samples, bio=mb_R_plot, FM=FM_R_plot, HR=HR_R_plot, S=S_obs, PM=PM, prediction_time_ss=prediction_time_ss, fall_time_ss=fall_time_ss, survey_time=survey_time, yrs=yrs)
    ggsave(filename=file.path(outputs_dir, paste0("plot_hcr_", aulab, ".", outformat)), plot=pl_hcr, device=outformat)

    pl_hcr_fall = fishery_model_plot(toplot="harvest_control_rule_fall", n_sample=1000, res=samples, bio=mb_R_plot, FM=FM_R_plot, HR=HR_R_plot, S=S_obs, PM=PM, prediction_time_ss=prediction_time_ss, fall_time_ss=fall_time_ss, survey_time=survey_time, yrs=yrs)
    ggsave(filename=file.path(outputs_dir, paste0("plot_hcr_fall_", aulab, ".", outformat)), plot=pl_hcr_fall, device=outformat)

    # Prior vs Posterior
    if (!is.null(prior) && inherits(prior, "mcmc.list")) {
      prior_vars = coda::varnames(prior)
      post_vars = coda::varnames(samples)
      common_vars = intersect(prior_vars, post_vars)

      for (p_var in c("K", "r", "q1", "bpsd", "bosd", "M", "rho", "wk", "sigmaR")) {
        if (p_var %in% common_vars) {
          bw_val = if (p_var == "K") 0.04 else NULL
          pl_prior = plot_prior_posterior(p_var, prior_samples = prior, posterior_samples = samples, bw = bw_val)
          ggsave(filename = file.path(outputs_dir, paste0("plot_prior_", p_var, "_", aulab, ".", outformat)), plot = pl_prior, device = outformat)
        }
      }
    }

    # Other Diagnostics
    has_r = "r" %in% coda::varnames(samples)
    if (has_r) {
      pl_sp = fishery_model_plot(toplot="surplus_production", n_sample=1000, res=samples, bio=mb_R_plot, FM=FM_R_plot, HR=HR_R_plot, S=S_obs, PM=PM, prediction_time_ss=prediction_time_ss, fall_time_ss=fall_time_ss, survey_time=survey_time, yrs=yrs)
      ggsave(filename=file.path(outputs_dir, paste0("plot_surplus_production_", aulab, ".", outformat)), plot=pl_sp, device=outformat)
    } else {
      cat("Skipping surplus production plot (only valid for models with intrinsic growth rate 'r').\n")
    }

    pl_ss = fishery_model_plot(toplot="state_space", n_sample=1000, res=samples, bio=mb_R_plot, FM=FM_R_plot, HR=HR_R_plot, S=S_obs, PM=PM, prediction_time_ss=prediction_time_ss, fall_time_ss=fall_time_ss, survey_time=survey_time, yrs=yrs)
    ggsave(filename=file.path(outputs_dir, paste0("plot_state_space_", aulab, ".", outformat)), plot=pl_ss, device=outformat)

    cat("Plotting complete.\n")
  }
}

plot_retrospective_analysis_results = function(retrospective_results, setup, rho_target_param, outformat="png", outputs_dir = setup$outputs_dir) {
  library(ggplot2)
  library(tidyr)

  aulab = setup$aulab
  PM = setup$PM

  # Extract values specifically for rho_target_param
  target_full_mean = retrospective_results$full_param_means[rho_target_param]
  target_retro_means = retrospective_results$retrospective_param_means[, rho_target_param]

  # Create a dataframe for plotting
  plot_data = data.frame(
    years_truncated = 0:length(target_retro_means),
    param_mean = c(target_full_mean, target_retro_means)
  )

  # Create the plot
  retrospective_plot = ggplot(plot_data, aes(x = years_truncated, y = param_mean)) +
    geom_point(color = "blue", size = 3) +
    geom_line(color = "blue", linetype = "dashed") +
    geom_hline(yintercept = target_full_mean, linetype = "solid", color = "red", linewidth = 1) +
    annotate("text", x = max(plot_data$years_truncated) * 0.7, y = target_full_mean, label = "Full Model Estimate", color = "red", vjust = -0.5) +
    labs(
      title = paste("Retrospective Analysis of", rho_target_param),
      x = "Years Truncated (0 = Full Model)",
      y = paste("Mean of", rho_target_param)
    ) +
    theme_bw()

  print(retrospective_plot)
  ggsave(filename=file.path(outputs_dir, paste0("plot_retrospective_", rho_target_param, "_", aulab, ".", outformat)), plot=retrospective_plot, device=outformat)

  # Visualize Retrospective Analysis Results for Multiple Parameters
  monitored_params = names(retrospective_results$full_param_means)

  df_full = data.frame(
    years_truncated = 0,
    Parameter = monitored_params,
    param_mean = as.numeric(retrospective_results$full_param_means)
  )

  df_retrospective = as.data.frame(retrospective_results$retrospective_param_means)
  df_retrospective$years_truncated = 1:nrow(df_retrospective)

  df_retrospective_long = df_retrospective %>%
    pivot_longer(cols = all_of(monitored_params), names_to = "Parameter", values_to = "param_mean")

  plot_data = rbind(df_full, df_retrospective_long)

  hline_data = data.frame(
    Parameter = monitored_params,
    yintercept = as.numeric(retrospective_results$full_param_means)
  )

  retrospective_plot_multi = ggplot(plot_data, aes(x = years_truncated, y = param_mean, color = Parameter)) +
    geom_point(size = 3) +
    geom_line(linetype = "dashed") +
    geom_hline(data = hline_data, aes(yintercept = yintercept, color = Parameter), linetype = "solid", linewidth = 1, show.legend = FALSE) +
    geom_text(data = hline_data, aes(x = max(plot_data$years_truncated) * 0.7, y = yintercept, label = paste("Full model", Parameter)),
              vjust = -0.5, hjust = 0, size = 3, show.legend = FALSE) +
    labs(
      title = "Retrospective Analysis of Monitored Parameters",
      x = "Years Truncated (0 = Full Model)",
      y = "Mean Parameter Estimate"
    ) +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank())

  print(retrospective_plot_multi)
  ggsave(filename=file.path(outputs_dir, paste0("plot_retrospective_multi_", aulab, ".", outformat)), plot=retrospective_plot_multi, device=outformat)

  # --- Plot Retrospective Biomass Time Series ---
  cat("\nExtracting and plotting retrospective Biomass time series...\n")

  base_yrs = PM$yrs

  get_median_B = function(samples_mat, nT) {
    K_s = samples_mat[, "K"]
    m_cols = paste0("m[", 1:nT, "]")
    m_s = samples_mat[, m_cols, drop=FALSE]
    B_mat = m_s * K_s
    apply(B_mat, 2, median)
  }

  df_ts_list = list()

  full_mat = as.matrix(retrospective_results$full_model_samples)
  full_B = get_median_B(full_mat, length(base_yrs))
  df_ts_list[[1]] = data.frame(
    Year = base_yrs,
    Biomass = full_B,
    Model = "Full Model"
  )

  for(k in 1:length(retrospective_results$retrospective_model_samples)) {
    retro_mat = as.matrix(retrospective_results$retrospective_model_samples[[k]])
    retro_nT = length(base_yrs) - k
    retro_B = get_median_B(retro_mat, retro_nT)
    df_ts_list[[k+1]] = data.frame(
      Year = base_yrs[1:retro_nT],
      Biomass = retro_B,
      Model = paste0("Truncated Yr ", k)
    )
  }

  df_ts_all = do.call(rbind, df_ts_list)
  df_ts_all$Model = factor(df_ts_all$Model, levels = c("Full Model", paste0("Truncated Yr ", 1:length(retrospective_results$retrospective_model_samples))))

  retro_ts_plot = ggplot(df_ts_all, aes(x = Year, y = Biomass, color = Model, size = Model)) +
    geom_line() +
    scale_color_manual(values = c("Full Model" = "black", setNames(rainbow(length(retrospective_results$retrospective_model_samples)), paste0("Truncated Yr ", 1:length(retrospective_results$retrospective_model_samples))))) +
    scale_size_manual(values = c("Full Model" = 1.5, setNames(rep(0.8, length(retrospective_results$retrospective_model_samples)), paste0("Truncated Yr ", 1:length(retrospective_results$retrospective_model_samples))))) +
    labs(title = "Retrospective Analysis: Biomass Time Series", x = "Year", y = "Median Biomass") +
    theme_bw() +
    theme(legend.position = "right")

  print(retro_ts_plot)
  ggsave(filename=file.path(outputs_dir, paste0("plot_retrospective_biomass_ts_", aulab, ".", outformat)), plot=retro_ts_plot, device=outformat)
}




#' Run state-space truncated Gaussian surplus production model
#'
#' @description
#' Executes data preparation, Nimble model compilation, MCMC sampling, WAIC
#' calculation, convergence diagnostics, prior sampling, posterior predictive
#' checks, post-processing, and retrospective analysis for the state-space
#' truncated Gaussian surplus production model.
#'
#' Supports both production mode (5,000 iterations, 4 chains, 5-year retrospective)
#' and fast test mode (1,000 iterations, 2 chains, 1-year retrospective).
#'
#' @param year_assessment Assessment year (default 2025).
#' @param aulab Subarea label: "cfanorth", "cfasouth", or "cfa4x" (default "cfanorth").
#' @param outputs_dir Character path to output directory.
#' @param is_test Logical; if TRUE, uses reduced MCMC iterations and shortened retrospective.
#' @param niter Total MCMC iterations per chain (default 5000 in production, 1000 in test).
#' @param nburnin Burn-in iterations per chain (default 1000 in production, 500 in test).
#' @param nchains Number of parallel MCMC chains (default 4 in production, 2 in test).
#' @param do_retrospective Logical; whether to execute retrospective analysis (default TRUE).
#' @param truncate_max_years Number of years to truncate in retrospective analysis.
#' @param do_plots Logical; whether to generate and save diagnostic plots (default TRUE).
#' @param outformat Character graphics format (default "png").
#' @return A list containing \code{samples}, \code{setup}, \code{WAIC},
#'   \code{diagnostics}, \code{PPC}, and optionally \code{retrospective}.
run_biomass_dynamics_gaussian = function(
    year_assessment = 2025,
    aulab = "cfanorth",
    outputs_dir = "c:/home/jae/projects/snowcrab_nimble/output",
    is_test = FALSE,
    niter = if (is_test) 1000 else 5000,
    nburnin = if (is_test) 500 else 1000,
    nchains = if (is_test) 2 else 4,
    do_retrospective = TRUE,
    truncate_max_years = if (is_test) 1 else 5,
    do_plots = TRUE,
    outformat = "png"
) {
  mode_label = if (is_test) "TEST (Fast)" else "PRODUCTION"
  cat("=================================================================\n")
  cat(sprintf("Truncated Gaussian Surplus Production Model: %s Mode\n", mode_label))
  cat(sprintf("Area: %s | Year: %d | Chains: %d | Iter: %d | Burnin: %d\n",
              aulab, year_assessment, nchains, niter, nburnin))
  cat("=================================================================\n\n")


if (0) {
  

    setwd("c:/home/jae/projects/snowcrab_nimble/")

    source("biomass_dynamics_functions.R")
    
    year_assessment = 2025
    
    aulab = "cfanorth"
    outputs_dir = "c:/home/jae/projects/snowcrab_nimble/output"
    is_test = FALSE
    niter = if (is_test) 1000 else 5000
    nburnin = if (is_test) 500 else 1000  
    nchains = if (is_test) 2 else 4
    do_retrospective = TRUE
    truncate_max_years = if (is_test) 1 else 5
    do_plots = TRUE
    outformat = "png"
}


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
  n_prior = if (is_test) 2000 else 10000
  n_ppc = if (is_test) 200 else 500
  prefix = if (is_test) paste0(aulab, "_gaussian_test_") else paste0(aulab, "_gaussian_")

  prior = sample_priors_nimble(setup, model_type = "gaussian", n_samples = n_prior)
  ppc = posterior_predictive_check(
    samples = samples,
    setup = setup,
    model_type = "gaussian",
    n_sims = n_ppc,
    output_dir = setup$outputs_dir,
    prefix = prefix
  )

  # 6. Post-processing, figures, and zone classification
  output_prefix = if (is_test) "results_nimble_gaussian_test" else "results_nimble_gaussian"
  post_process_and_save_fishery_model(
    samples = samples,
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

  cat("\nExecution of run_biomass_dynamics_gaussian completed successfully.\n")
  return(list(
    samples = samples,
    setup = setup,
    WAIC = waic_res,
    prior = prior,
    ppc = ppc,
    retrospective = retrospective_results
  ))
}


#' state-space lognormal surplus production model
#'
#' @description
#' Executes data preparation, Nimble model compilation, MCMC sampling, WAIC
#' calculation, convergence diagnostics, prior sampling, posterior predictive
#' checks, post-processing, and retrospective analysis for the state-space lognormal
#' surplus production model.
#'
#' Supports both production mode (5,000 iterations, 4 chains, 5-year retrospective)
#' and fast test mode (1,000 iterations, 2 chains, 1-year retrospective).
#'
#' @param year_assessment Assessment year (default 2025).
#' @param aulab Subarea label: "cfanorth", "cfasouth", or "cfa4x" (default "cfanorth").
#' @param outputs_dir Character path to output directory.
#' @param is_test Logical; if TRUE, uses reduced MCMC iterations and shortened retrospective.
#' @param niter Total MCMC iterations per chain (default 5000 in production, 1000 in test).
#' @param nburnin Burn-in iterations per chain (default 1000 in production, 500 in test).
#' @param nchains Number of parallel MCMC chains (default 4 in production, 2 in test).
#' @param do_retrospective Logical; whether to execute retrospective analysis (default TRUE).
#' @param truncate_max_years Number of years to truncate in retrospective analysis.
#' @param do_plots Logical; whether to generate and save diagnostic plots (default TRUE).
#' @param outformat Character graphics format (default "png").
#' @return A list containing \code{samples}, \code{setup}, \code{WAIC},
#'   \code{diagnostics}, \code{PPC}, and optionally \code{retrospective}.
run_biomass_dynamics_lognormal = function(
    year_assessment = 2025,
    aulab = "cfanorth",
    outputs_dir = "c:/home/jae/projects/snowcrab_nimble/output",
    is_test = FALSE,
    niter = if (is_test) 1000 else 5000,
    nburnin = if (is_test) 500 else 1000,
    nchains = if (is_test) 2 else 4,
    do_retrospective = TRUE,
    truncate_max_years = if (is_test) 1 else 5,
    do_plots = TRUE,
    outformat = "png"
) {
  mode_label = if (is_test) "TEST (Fast)" else "PRODUCTION"
  cat("=================================================================\n")
  cat(sprintf("Lognormal Surplus Production Model: %s Mode\n", mode_label))
  cat(sprintf("Area: %s | Year: %d | Chains: %d | Iter: %d | Burnin: %d\n",
              aulab, year_assessment, nchains, niter, nburnin))
  cat("=================================================================\n\n")

  # 1. Prepare data and parameter setup
  setup = prepare_fishery_model_data(
    year_assessment = year_assessment,
    aulab = aulab,
    outputs_dir = outputs_dir,
    model_variation = "logistic_discrete_historical"
  )

  # 2. Build model inputs using factored model library
  model_inputs = prepare_nimble_model_inputs(setup, model_type = "lognormal")

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
  n_prior = if (is_test) 2000 else 10000
  n_ppc = if (is_test) 200 else 500
  prefix = if (is_test) paste0(aulab, "_test_") else paste0(aulab, "_")

  prior = sample_priors_nimble(setup, model_type = "lognormal", n_samples = n_prior)
  ppc = posterior_predictive_check(
    samples = samples,
    setup = setup,
    model_type = "lognormal",
    n_sims = n_ppc,
    output_dir = setup$outputs_dir,
    prefix = prefix
  )

  # 6. Post-processing, figures, and zone classification
  output_prefix = if (is_test) "results_nimble_test" else "results_nimble"
  post_process_and_save_fishery_model(
    samples = samples,
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

  cat("\nExecution of run_biomass_dynamics_lognormal completed successfully.\n")
  return(list(
    samples = samples,
    setup = setup,
    WAIC = waic_res,
    prior = prior,
    ppc = ppc,
    retrospective = retrospective_results
  ))
}



#' Deriso-Schnute delay-difference biomass dynamics model
#'
#' @description
#' Executes data preparation, Nimble model compilation, MCMC sampling, WAIC
#' calculation, convergence diagnostics, prior sampling, posterior predictive
#' checks, numerical reference point estimation ($F_{0.1}$, $B_{msy}$, $\text{MSY}$),
#' post-processing, and retrospective analysis for the Deriso-Schnute
#' delay-difference model with AR(1) recruitment.
#'
#' Supports both production mode (2,000 iterations, 2 chains, 5-year retrospective)
#' and fast test mode (1,000 iterations, 2 chains, 1-year retrospective).
#'
#' @param year_assessment Assessment year (default 2025).
#' @param aulab Subarea label: "cfanorth", "cfasouth", or "cfa4x" (default "cfanorth").
#' @param outputs_dir Character path to output directory.
#' @param is_test Logical; if TRUE, uses reduced MCMC iterations and shortened retrospective.
#' @param niter Total MCMC iterations per chain (default 2000 in production, 1000 in test).
#' @param nburnin Burn-in iterations per chain (default 1000 in production, 500 in test).
#' @param nchains Number of parallel MCMC chains (default 2).
#' @param do_retrospective Logical; whether to execute retrospective analysis (default TRUE).
#' @param truncate_max_years Number of years to truncate in retrospective analysis.
#' @param do_plots Logical; whether to generate and save diagnostic plots (default TRUE).
#' @param outformat Character graphics format (default "png").
#' @return A list containing \code{samples}, \code{setup}, \code{WAIC},
#'   \code{diagnostics}, \code{PPC}, and optionally \code{retrospective}.
run_delay_difference = function(
    year_assessment = 2025,
    aulab = "cfanorth",
    outputs_dir = "c:/home/jae/projects/snowcrab_nimble/output",
    is_test = FALSE,
    niter = if (is_test) 1000 else 2000,
    nburnin = if (is_test) 500 else 1000,
    nchains = 2,
    do_retrospective = TRUE,
    truncate_max_years = if (is_test) 1 else 5,
    do_plots = TRUE,
    outformat = "png"
) {
  mode_label = if (is_test) "TEST (Fast)" else "PRODUCTION"
  cat("=================================================================\n")
  cat(sprintf("Deriso-Schnute Delay-Difference Model: %s Mode\n", mode_label))
  cat(sprintf("Area: %s | Year: %d | Chains: %d | Iter: %d | Burnin: %d\n",
              aulab, year_assessment, nchains, niter, nburnin))
  cat("=================================================================\n\n")

  # 1. Prepare data and parameter setup
  setup = prepare_fishery_model_data(
    year_assessment = year_assessment,
    aulab = aulab,
    outputs_dir = outputs_dir,
    model_variation = "delay_difference"
  )

  # 2. Build model inputs using factored model library
  model_inputs = prepare_nimble_model_inputs( setup, model_type = "delay_difference" )

  nimble_model = nimbleModel(
    code = model_inputs$code,
    constants = model_inputs$constants,
    data = model_inputs$data,
    dimensions = model_inputs$dimensions,
    inits = model_inputs$inits,
    buildDerivs = model_inputs$buildDerivs, # FALSE for delay-difference
    check = TRUE
  )

  C_nimble_model = compileNimble(nimble_model)

  # 3. Configure MCMC
  mcmc_conf = configureMCMC(
    nimble_model,
    monitors = model_inputs$monitors,
    enableWAIC = TRUE
  )
  mcmc_builder = buildMCMC(mcmc_conf)
  C_mcmc_builder = compileNimble(mcmc_builder, project = nimble_model)

  # 4. Run MCMC sampling
  cat("\nRunning MCMC sampling for Delay-Difference Model...\n")
  mcmc_out = runMCMC(
    C_mcmc_builder,
    niter = niter,
    nburnin = nburnin,
    nchains = nchains,
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
  n_prior = if (is_test) 2000 else 5000
  n_ppc = if (is_test) 200 else 500
  prefix = if (is_test) paste0(aulab, "_delay_diff_test_") else paste0(aulab, "_delay_diff_")

  prior = sample_priors_nimble(setup, model_type = "delay_difference", n_samples = n_prior)
  ppc = posterior_predictive_check(
    samples = samples,
    setup = setup,
    model_type = "delay_difference",
    n_sims = n_ppc,
    output_dir = setup$outputs_dir,
    prefix = prefix
  )

  # 6. Post-processing, figures, and zone classification
  output_prefix = if (is_test) "results_nimble_delay_diff_test" else "results_nimble_delay_diff"
  post_process_and_save_fishery_model(
    samples = as.mcmc.list(samples),
    setup = setup,
    prior = prior,
    do_plots = do_plots,
    do_hcr_probabilities = TRUE,
    outformat = outformat,
    output_prefix = output_prefix
  )

  # 7. Retrospective analysis (buildDerivs = FALSE)
  retrospective_results = NULL
  if (do_retrospective && truncate_max_years > 0) {
    cat(sprintf("\nStarting retrospective analysis (truncating up to %d years)...\n",
                truncate_max_years))
    retro_iter = if (is_test) 500 else 1000
    retro_burnin = if (is_test) 200 else 500

    retrospective_mcmc_settings = list(
      niter = retro_iter,
      nburnin = retro_burnin,
      nchains = 2,
      progressBar = TRUE,
      samplesAsCodaMCMC = TRUE
    )

    monitor_params_for_rho = c("K", "q1", "bpsd", "bosd", "M", "rho", "wk", "sigmaR", "phi_R", "mu_R", "q_R")
    rho_target_param = "K"

    retrospective_results = do_retrospective_analysis(
      nimble_code = model_inputs$code,
      PM_base = setup$PM,
      nimble_constants_base = model_inputs$constants,
      model_data_base = model_inputs$data,
      initial_values_base = model_inputs$inits,
      monitor_params = monitor_params_for_rho,
      truncate_max_years = truncate_max_years,
      rho_target_param = rho_target_param,
      mcmc_settings = retrospective_mcmc_settings,
      buildDerivs = FALSE
    )

    plot_retrospective_analysis_results(
      retrospective_results,
      setup,
      rho_target_param,
      outformat = outformat
    )
  }

  cat("\nExecution of run_delay_difference completed successfully.\n")
  return(list(
    samples = samples,
    setup = setup,
    WAIC = waic_res,
    prior = prior,
    ppc = ppc,
    retrospective = retrospective_results
  ))
}
