#' Default Parameter Settings for Prior Distributions
#'
#' Default parameters for the prior distributions used in the \code{makePriors} function.
#'
#' @format A list with the following components:
#' \describe{
#'   \item{asymptote}{A list with components \code{g1} and \code{g2}, default values for the asymptote parameters.}
#'   \item{threshold}{A list with components \code{min} and \code{max}, default values for the threshold parameters.}
#'   \item{median}{A list with components \code{m1} and \code{m2}, default values for the median parameters.}
#'   \item{first_quartile}{A list with components \code{q1} and \code{q2}, default values for the first quartile parameters.}
#' }
#' @export
prior_params_default <- list(
  asymptote = list(g1 = 1, g2 = 1),
  threshold = list(min = 15, max = 35),
  median = list(m1 = 2, m2 = 2),
  first_quartile = list(q1 = 6, q2 = 3)
)

# TODO DEFAULT RISK PROPORTION DATA
# Decide what to do here/how to consider censoring

#' Default Distribution Data
#'
#' Default data frame structure with row names for use in the \code{makePriors} function.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{age}{Age values (NA for default).}
#'   \item{at_risk}{Count of people at risk (NA for default).}
#' }
#' @export
distribution_data_default <- data.frame(
  row.names = c("min", "first_quartile", "median", "max"),
  age = c(NA, NA, NA, NA),
  at_risk = c(NA, NA, NA, NA)
)
# Checks a prior_params list and expands it to the suffixed keys of the active mode. For each parameter and each active slot, the suffixed entry is kept if present, otherwise the bare entry is copied in; the bare entries are then dropped. Used by makePriors() and penetrance().
expand_prior_params_internal <- function(prior_params, sex_specific) {
  if (!is.list(prior_params)) {
    stop("Error: 'prior_params' must be a list.")
  }

  active_sex_slots <- if (sex_specific) c("male", "female") else "noSex"
  wrong_mode_slots <- if (sex_specific) "noSex" else c("male", "female")
  positive_fields <- list(asymptote = c("g1", "g2"), median = c("m1", "m2"), first_quartile = c("q1", "q2"))
  known_keys <- as.vector(outer(names(prior_params_default), c("", paste0("_", c("noSex", "male", "female"))), paste0))
  expanded_params <- prior_params[setdiff(names(prior_params), known_keys)]

  for (parameter_name in names(prior_params_default)) {
    for (sex_slot in wrong_mode_slots) {
      wrong_key <- paste0(parameter_name, "_", sex_slot)
      if (!is.null(prior_params[[wrong_key]])) {
        stop(paste0("Error: 'prior_params$", wrong_key, "' cannot be used when 'sex_specific' is ", sex_specific))
      }
    }

    for (sex_slot in active_sex_slots) {
      suffixed_key <- paste0(parameter_name, "_", sex_slot)
      key_name <- suffixed_key
      prior_value <- prior_params[[key_name]]
      if (is.null(prior_value)) {
        key_name <- parameter_name
        prior_value <- prior_params[[parameter_name]]
      }
      if (is.null(prior_value)) {
        stop(paste0("Error: 'prior_params' must contain '", suffixed_key, "' or the bare '", parameter_name, "'."))
      }
      if (!is.list(prior_value)) {
        stop(paste0("Error: 'prior_params$", key_name, "' must be a list."))
      }

      if (parameter_name == "threshold") {
        if (!is.numeric(prior_value[["min"]]) || !is.numeric(prior_value[["max"]]) || prior_value[["min"]] >= prior_value[["max"]]) {
          stop(paste0("Error: 'prior_params$", key_name, "$min' must be numeric and less than 'prior_params$", key_name, "$max'."))
        }
      } else {
        for (field_name in positive_fields[[parameter_name]]) {
          field_value <- prior_value[[field_name]]
          if (!is.numeric(field_value) || length(field_value) != 1 || !is.finite(field_value) || field_value <= 0) {
            stop(paste0("Error: 'prior_params$", key_name, "$", field_name, "' must be a positive numeric value."))
          }
        }
      }

      expanded_params[[suffixed_key]] <- prior_value
    }
  }

  return(expanded_params)
}

#' Make Priors
#'
#' This function generates prior distributions based on user input or default parameters.
#' It is designed to aid in the statistical analysis of risk proportions in populations, particularly in the context of cancer research.
#' The distributions are calculated for various statistical metrics such as asymptote, threshold, median, and first quartile.
#'
#' @param data A data frame with row names min, first_quartile, median, and max, and columns age and at_risk (a count of people at risk, not a proportion). If NULL or entirely NA, default parameters are used. If provided, all entries in the age and at_risk columns must be numeric and non-NA, ages must strictly increase from min to max, at_risk must not increase from min to max, and no at_risk count may exceed sample_size.
#' @param sample_size Numeric, positive, the total number of people in the source study. Required when data is provided, since the asymptote is estimated from the fraction of sample_size no longer at risk by the median age.
#' @param ratio List, with components noSex, male, and female, each the odds ratio (OR) or relative risk (RR) used in asymptote parameter calculations for the corresponding baseline. The components must match sex_specific: noSex when FALSE, male and/or female when TRUE. If sex_specific is TRUE, a ratio may be given for one sex only and the other keeps its asymptote prior from another method.
#' @param ratio_is_or Logical, indicating whether the ratios supplied are odds ratios (OR) (TRUE) or relative risks (FALSE).
#' @param ratio_ci_lower List, with components noSex, male, and female, the lower bound of each reported ratio's confidence interval (on the same OR/RR scale as ratio, per ratio_is_or). Required alongside the matching ratio component(s).
#' @param ratio_ci_upper List, with components noSex, male, and female, the upper bound of each reported ratio's confidence interval (on the same OR/RR scale as ratio, per ratio_is_or). Required alongside the matching ratio component(s).
#' @param ci_level Numeric strictly between 0 and 1, the coverage of the confidence intervals given in ratio_ci_lower and ratio_ci_upper. Default is 0.95.
#' @param prior_params List of prior hyperparameters for the beta distributions. Can be supplied bare or with suffix (_male, _female, or _noSex. If NULL, prior_params_default is used.
#' @param baseline_data Baseline risk data, required with ratio. A data frame with Male and Female columns if sex_specific is TRUE, otherwise a numeric vector or single-column data frame.
#' @param lifetime_risk Numeric between 0 and 1, an externally reported lifetime risk used as the asymptote prior's mean when asymptote_method = "lifetime_risk".
#' @param lifetime_risk_cv Numeric, positive. The assumed coefficient of variation (relative SD) of lifetime_risk, used to derive the asymptote prior's concentration when asymptote_method = "lifetime_risk". Default is 0.15.
#' @param asymptote_method Character, one of "default", "lifetime_risk", or "estimate". Selects how the asymptote prior is built for the time interval risk data pathway. "default" uses a flat Beta(1,1). "lifetime_risk" uses the user-supplied lifetime_risk as the Beta mean. "estimate" derives the asymptote from the time interval (KM-curve-style) risk data itself -- not yet implemented. NULL (default) leaves the asymptote prior untouched.
#' @param max_age Numeric, the maximum age to plot the sampled curves over when prior_predictive_check is TRUE. Default is 94, matching penetrance()'s default max_age.
#' @param prior_predictive_check Logical, indicating whether to draw n_draws samples from prior_params and plot the resulting penetrance curves overlaid, as a sanity check on the joint prior specification. Default is FALSE.
#' @param n_draws Numeric, the number of curves to draw and overlay when prior_predictive_check is TRUE. Default is 50.
#' @param sex_specific Logical. TRUE returns the priors separately for each sex, FALSE returns one for both (e.g. asymptote_noSex). Default is FALSE.
#'
#' @details
#' If data is provided (the quartile pathway), the median and first quartile priors are Beta distributions whose means are the reported ages rescaled onto [0, 1] and whose concentrations are the at_risk counts at those ages, capped at 30 so the prior stays only moderately informative. The threshold prior is uniform from 10 years below to 10 years above the reported minimum age, with the lower bound floored at 0.
#'
#' The reported minimum and maximum ages are stored in the result as lower_anchor and upper_anchor, and the sampler rescales the median and first quartile on those anchors so they are scored on the same scale they were elicited on.
#'
#' The asymptote prior is a Beta distribution with mean 2 * (1 - at_risk_median / sample_size) and concentration min(at_risk_median, 30). This follows from the reported median age being the age by which half of the eventually affected carriers are diagnosed, so the fraction diagnosed by that age is half the asymptote. If this estimate falls outside (0, 1), for example because follow-up in the source study was incomplete or heavily censored, an error is thrown.
#'
#' If ratio is supplied (the ratio pathway), only the asymptote prior is set. The baseline lifetime risk is the sum of baseline_data (per sex when sex_specific is TRUE). An odds ratio, with its CI bounds, is first converted to a relative risk. The asymptote prior is then the Beta whose lower, 50%, and upper quantiles (the ci_level interval) are closest to the baseline lifetime risk times the ratio's lower bound, point estimate, and upper bound. The reported ratio is treated as the median, so the prior's mean can differ slightly from the baseline times the ratio. Both Beta shapes are at least 1. An error is thrown if the implied lifetime risk at the upper bound is 1 or more.
#'
#' If sex_specific is TRUE, the priors derived from data are returned for both sexes (e.g. median_male and median_female) and are identical, since a single table carries no sex-specific information. If asymptote_method or ratio is also supplied, the asymptote prior from those pathways overrides the one derived from data.
#'
#' @return A list of the prior parameters for asymptote, threshold, median, and first_quartile, to be passed as prior_params to penetrance(). Every parameter is returned with a suffix (_male, _female, or _noSex).
#'
#' @seealso \code{\link{qbeta}}, \code{\link{runif}}
#' 
#' @export
makePriors <- function(data = NULL,
                      sample_size,
                      ratio = list(noSex = NULL, male = NULL, female = NULL), ratio_is_or,
                      ratio_ci_lower = list(noSex = NULL, male = NULL, female = NULL),
                      ratio_ci_upper = list(noSex = NULL, male = NULL, female = NULL),
                      ci_level = 0.95,
                      prior_params = prior_params_default, baseline_data,
                      asymptote_method = NULL,
                      lifetime_risk = NULL, lifetime_risk_cv = 0.15,
                      max_age = 94, prior_predictive_check = FALSE, n_draws = 50, sex_specific = FALSE) {

  data_supplied <- !(is.null(data) || all(is.na(data)))

  # Sex slots the priors are written for under the active mode
  active_sex_slots <- if (sex_specific) c("male", "female") else "noSex"

  # Every prior is stored under a sex-suffixed key
  prior_key <- function(parameter_name, sex_slot) {
    return(paste0(parameter_name, "_", sex_slot))
  }

  if (!data_supplied) {
    if (is.null(prior_params)) {
      prior_params <- prior_params_default
    }
    prior_params <- expand_prior_params_internal(prior_params, sex_specific)
  }

  # No data, ratio, or asymptote_method supplied: return prior_params as default (after the predictive check, if requested)
  if (!data_supplied && all(sapply(ratio, is.null)) && is.null(asymptote_method) && !prior_predictive_check) {
    return(prior_params)
  }

  # Converts OR to RR given the baseline risk P2:
  # OR = [P1/(1-P1)] / [P2/(1-P2)], solved for P1 and substituted into RR = P1/P2.
  # P1 = OR·P2 / (1 - P2 + OR·P2); where P1 = asymptote, P2 = baseline lifetime risk

  convert_or_to_rr <- function(or, baseline_risk) {
    return(or / (1 - baseline_risk + (or * baseline_risk)))
  }

  # Matches a Beta(g1, g2) to a target mean and variance, solving
  # Var = mean(1-mean)/(concentration+1) for concentration.
  compute_concentration <- function(mean, variance) {
    return((mean * (1 - mean)) / variance - 1)
  }

  # Fit the Beta whose lower, 50%, and upper quantiles (the ci_level interval) are closest, by squared error, to the targets
  # c(lower, point estimate, upper). The reported ratio is treated as the median.
  fit_beta_to_ci <- function(targets) {
    probs <- c((1 - ci_level) / 2, 0.5, 1 - (1 - ci_level) / 2)
    squared_error <- function(shapes) {
      return(sum((qbeta(probs, shapes[1], shapes[2]) - targets)^2))
    }
    fit <- optim(pmax(c(targets[2] * 10, (1 - targets[2]) * 10), 1), squared_error, method = "L-BFGS-B", lower = c(1, 1))
    return(list(g1 = fit$par[1], g2 = fit$par[2]))
  }

  # Checks that a ratio and its CI bounds are usable: all supplied, single finite positive
  # numbers, and ordered lower < ratio < upper.
  # Runs on the values as reported, before any OR-to-RR conversion.
  validate_ratio_ci <- function(ratio_value, ci_lower_value, ci_upper_value, slot) {
    if (is.null(ci_lower_value)) {
      stop(paste0("Error: 'ratio_ci_lower$", slot, "' must be supplied when 'ratio$", slot, "' is supplied."))
    }
    if (is.null(ci_upper_value)) {
      stop(paste0("Error: 'ratio_ci_upper$", slot, "' must be supplied when 'ratio$", slot, "' is supplied."))
    }
    if (!is.numeric(ratio_value) || length(ratio_value) != 1 || !is.finite(ratio_value) || ratio_value <= 0) {
      stop(paste0("Error: 'ratio$", slot, "' must be a single positive numeric value."))
    }
    if (!is.numeric(ci_lower_value) || length(ci_lower_value) != 1 || !is.finite(ci_lower_value) || ci_lower_value <= 0) {
      stop(paste0("Error: 'ratio_ci_lower$", slot, "' must be a single positive numeric value."))
    }
    if (!is.numeric(ci_upper_value) || length(ci_upper_value) != 1 || !is.finite(ci_upper_value) || ci_upper_value <= 0) {
      stop(paste0("Error: 'ratio_ci_upper$", slot, "' must be a single positive numeric value."))
    }
    if (!(ci_lower_value < ratio_value && ratio_value < ci_upper_value)) {
      stop(paste0("Error: 'ratio$", slot, "' must lie strictly between 'ratio_ci_lower$", slot, "' and 'ratio_ci_upper$", slot, "'."))
    }
  }

  normalize_median <- function(stat, min_age, max_age) {
    return((stat - min_age) / (max_age - min_age))
  }

  normalize_first_quartile <- function(stat, min_age, median_age) {
    return((stat - min_age) / (median_age - min_age))
  }

  # cap concentration at 30 to avoid extremely narrow distributions when at_risk is high
  max_concentration <- 30

  compute_parameters_median <- function(stat, at_risk, min_age, max_age) {
    median_norm <- normalize_median(stat, min_age, max_age)
    concentration <- min(at_risk, max_concentration)
    alpha <- median_norm * concentration
    beta <- concentration - alpha
    return(list(m1 = alpha, m2 = beta))
  }

  compute_parameters_quartile <- function(stat, at_risk, min_age, median_age) {
    quartile_norm <- normalize_first_quartile(stat, min_age, median_age)
    concentration <- min(at_risk, max_concentration)
    alpha <- quartile_norm * concentration
    beta <- concentration - alpha
    return(list(q1 = alpha, q2 = beta))
  }

  # Median ages are quantiles among eventually-affected carriers only, so G(median_age) = 0.5
  # and F(median_age) = asymptote * G(median_age), giving asymptote = 2 * F(median_age),
  # with F(median_age) estimated as 1 - at_risk(median_age) / sample_size.
  compute_parameters_asymptote <- function(at_risk, sample_size) {
    asymptote_estimate <- 2 * (1 - at_risk / sample_size)
    if (asymptote_estimate <= 0 || asymptote_estimate >= 1) {
      stop(paste0("Error: the asymptote implied by 'data' (2 * (1 - at_risk_median / sample_size) = ", round(asymptote_estimate, 3),
                  ") must be between 0 and 1. Check 'data$at_risk' and 'sample_size', and whether follow-up in the source study was incomplete or heavily censored."))
    }
    concentration <- min(at_risk, max_concentration)
    alpha <- asymptote_estimate * concentration
    beta <- concentration - alpha
    return(list(g1 = alpha, g2 = beta))
  }

  # Time interval risk data pathway: asymptote (skeleton only -- the
  # threshold/median/first_quartile estimation for this pathway, driven by
  # an age/at_risk table across many ages rather than just the 4 quartile
  # rows, is not yet built. This scaffolds the 3 user-selectable paths for
  # the asymptote step only, so the rest of the pathway can be built against
  # it later.
  compute_asymptote_time_interval <- function(asymptote_method, lifetime_risk, lifetime_risk_cv) {
    if (asymptote_method == "default") {
      # Path 1: no informative asymptote data provided, use flat Beta(1,1).
      return(list(g1 = 1, g2 = 1))
    }

    if (asymptote_method == "lifetime_risk") {
      # Path 2: user supplies a lifetime_risk estimate directly for use as a Beta mean.
      if (is.null(lifetime_risk)) {
        stop("Error: 'lifetime_risk' must be supplied when asymptote_method = 'lifetime_risk'.")
      }
      if (!is.numeric(lifetime_risk) || lifetime_risk <= 0 || lifetime_risk >= 1) {
        stop("Error: 'lifetime_risk' must be a numeric value between 0 and 1.")
      }
      if (!is.numeric(lifetime_risk_cv) || lifetime_risk_cv <= 0) {
        stop("Error: 'lifetime_risk_cv' must be a positive numeric value.")
      }

      # Concentration is derived from an assumed coefficient of variation
      # (relative SD) on lifetime_risk, since a bare point estimate carries no
      # reported CI to work with directly (unlike the ratio pathway). Reuses
      # the same Beta mean/variance moment-matching (compute_concentration)
      # the ratio pathway uses -- only the source of the target variance differs.
      se <- lifetime_risk_cv * lifetime_risk
      concentration <- compute_concentration(lifetime_risk, se^2)

      g1 <- lifetime_risk * concentration
      g2 <- concentration - g1
      return(list(g1 = g1, g2 = g2))
    }

    if (asymptote_method == "estimate") {
      # Path 3: derive the asymptote from the time interval risk data itself.
      # Mechanism not yet decided -- stub for
      # now, matching the compute_concentration_from_ci() stub convention in
      # priorElicitation.R.
      stop("Error: 'estimate' asymptote_method for the time interval risk data pathway is not yet implemented.")
    }

    stop("Error: 'asymptote_method' must be one of 'default', 'lifetime_risk', or 'estimate'.")
  }

  # Quartile pathway
  # The literature's max age is read from data directly rather than assigned to a local,
  # since max_age is already an argument (used by the prior predictive check).
  if (data_supplied) {
    if (!is.data.frame(data) || !all(c("age", "at_risk") %in% names(data))) {
      stop("Error: 'data' must be a data frame with columns 'age' and 'at_risk'.")
    }
    if (!identical(sort(rownames(data)), sort(c("min", "first_quartile", "median", "max")))) {
      stop("Error: 'data' must have row names 'min', 'first_quartile', 'median', and 'max'.")
    }
    if (any(is.na(data$age)) || any(!sapply(data$age, is.numeric))) {
      stop("Missing or non-numeric age entries in the data. Add numeric ages.")
    }
    if (any(is.na(data$at_risk)) || any(!sapply(data$at_risk, is.numeric))) {
      stop("Missing or non-numeric risk entries in the data. Add individuals at risk.")
    }
    if (missing(sample_size) || is.null(sample_size) || !is.numeric(sample_size) || sample_size <= 0) {
      stop("Error: 'sample_size' must be a positive numeric value when 'data' is supplied.")
    }
    if (any(data$at_risk > sample_size)) {
      stop("Error: 'data$at_risk' cannot exceed 'sample_size'.")
    }

    min_age <- data["min", "age"]
    first_quartile_age <- data["first_quartile", "age"]
    median_age <- data["median", "age"]

    if (!(min_age < first_quartile_age && first_quartile_age < median_age && median_age < data["max", "age"])) {
      stop("Error: 'data$age' must be strictly increasing from min to first_quartile to median to max.")
    }

    risk_median <- data["median", "at_risk"]
    risk_first_quartile <- data["first_quartile", "at_risk"]

    # Nobody re-enters the risk set, so the count at risk cannot increase with age
    if (!(data["min", "at_risk"] >= risk_first_quartile && risk_first_quartile >= risk_median && risk_median >= data["max", "at_risk"])) {
      stop("Error: 'data$at_risk' must not increase from min to first_quartile to median to max.")
    }

    res_median <- compute_parameters_median(median_age, risk_median, min_age, data["max", "age"])
    res_first_quartile <- compute_parameters_quartile(first_quartile_age, risk_first_quartile, min_age, median_age)
    res_asymptote <- compute_parameters_asymptote(risk_median, sample_size)

    quartile_priors <- list(
      asymptote = list(g1 = res_asymptote$g1, g2 = res_asymptote$g2),
      threshold = list(min = max(0, min_age - 10), max = min_age + 10),
      median = list(m1 = res_median$m1, m2 = res_median$m2),
      first_quartile = list(q1 = res_first_quartile$q1, q2 = res_first_quartile$q2)
    )

    # The one table is the only information available, so it is used for every slot of the active mode
    prior_params <- list()
    for (parameter_name in names(quartile_priors)) {
      for (sex_slot in active_sex_slots) {
        prior_params[[prior_key(parameter_name, sex_slot)]] <- quartile_priors[[parameter_name]]
      }
    }

    # The median and first quartile Betas were elicited on [min, max] age anchors, so the chain scores them on the same anchors
    prior_params$lower_anchor <- min_age
    prior_params$upper_anchor <- data["max", "age"]
  }

  # Time interval risk data pathway
  if (!is.null(asymptote_method)) {
    asymptote_prior <- compute_asymptote_time_interval(asymptote_method, lifetime_risk, lifetime_risk_cv)
    for (sex_slot in active_sex_slots) {
      prior_params[[prior_key("asymptote", sex_slot)]] <- asymptote_prior
    }
  }

  # Ratio pathway
  if (!all(sapply(ratio, is.null))) {

    if (sex_specific && !is.null(ratio$noSex)) {
      stop("Error: 'ratio$noSex' cannot be used when 'sex_specific' is TRUE.")
    }
    if (!sex_specific && (!is.null(ratio$male) || !is.null(ratio$female))) {
      stop("Error: 'ratio$male' and 'ratio$female' cannot be used when 'sex_specific' is FALSE.")
    }
    validate_baseline_data(baseline_data, sex_specific)

    if (!is.numeric(ci_level) || length(ci_level) != 1 || !is.finite(ci_level) || ci_level <= 0 || ci_level >= 1) {
      stop("Error: 'ci_level' must be a single number between 0 and 1.")
    }

    # Fit the asymptote's prior for each slot that has a ratio
    for (sex_slot in c("male", "female", "noSex")) {
      if (is.null(ratio[[sex_slot]])) next

      validate_ratio_ci(ratio[[sex_slot]], ratio_ci_lower[[sex_slot]], ratio_ci_upper[[sex_slot]], sex_slot)
      baseline_lifetime_risk <- switch(sex_slot,
                                       male = sum(baseline_data$Male),
                                       female = sum(baseline_data$Female),
                                       noSex = sum(baseline_data))

      ratio_values <- c(lower = ratio_ci_lower[[sex_slot]], point = ratio[[sex_slot]], upper = ratio_ci_upper[[sex_slot]])
      if (ratio_is_or == TRUE) {
        ratio_values <- convert_or_to_rr(ratio_values, baseline_lifetime_risk)
      }

      asymptote_targets <- baseline_lifetime_risk * ratio_values
      if (asymptote_targets[["upper"]] >= 1) {
        stop(paste0("Error: the carrier lifetime risk implied by the upper CI bound of 'ratio$", sex_slot, "' (", round(asymptote_targets[["upper"]], 3), ") must be less than 1."))
      }

      prior_params[[prior_key("asymptote", sex_slot)]] <- fit_beta_to_ci(asymptote_targets)
    }

    # If sex_specific = TRUE and one sex lacks a ratio, its asymptote is derived from asymptote_method, data, or prior_params
    if (sex_specific && (!is.null(ratio$male) || !is.null(ratio$female))) {
      if (!is.null(asymptote_method)) {
        source_label <- "'asymptote_method'"
      } else if (data_supplied) {
        source_label <- "'data'"
      } else {
        source_label <- "'prior_params'"
      }

      for (sex_slot in c("male", "female")) {
        if (!is.null(ratio[[sex_slot]])) next
        message("Note: no ratio supplied for ", sex_slot, ", so '", prior_key("asymptote", sex_slot), "' derived from ", source_label, ".")
      }
    }
  }

  # Prior predictive check: draw n_draws samples from prior_params and plot the
  # resulting penetrance curves overlaid, to sanity-check the joint prior.
  if (prior_predictive_check) {
    if (is.null(max_age)) {
      stop("Error: 'max_age' must be supplied when 'prior_predictive_check' is TRUE.")
    }

    x_values <- seq(0, max_age, length.out = max_age + 1)

    weibull_values <- function(alpha, beta, threshold, asymptote) {
      pweibull(x_values - threshold, shape = alpha, scale = beta) * asymptote
    }

    # Step 1: draw n_draws samples from each of the four priors
    draw_params <- function(sex_slot) {
      threshold_params <- prior_params[[prior_key("threshold", sex_slot)]]
      median_params <- prior_params[[prior_key("median", sex_slot)]]
      first_quartile_params <- prior_params[[prior_key("first_quartile", sex_slot)]]
      asymptote_params <- prior_params[[prior_key("asymptote", sex_slot)]]

      threshold_draws <- runif(n_draws, threshold_params$min, threshold_params$max)
      scaled_median_draws <- qbeta(runif(n_draws), median_params$m1, median_params$m2)
      scaled_first_quartile_draws <- qbeta(runif(n_draws), first_quartile_params$q1, first_quartile_params$q2)
      asymptote_draws <- qbeta(runif(n_draws), asymptote_params$g1, asymptote_params$g2)

      # Un-normalize onto actual ages (mirrors the scaling relations in mhChain.R): fixed anchors if stored, otherwise the drawn threshold and max_age
      lower_anchor <- if (is.null(prior_params$lower_anchor)) threshold_draws else prior_params$lower_anchor
      upper_anchor <- if (is.null(prior_params$upper_anchor)) max_age else prior_params$upper_anchor
      median_draws <- lower_anchor + scaled_median_draws * (upper_anchor - lower_anchor)
      first_quartile_draws <- lower_anchor + scaled_first_quartile_draws * (median_draws - lower_anchor)

      list(threshold = threshold_draws, median = median_draws,
           first_quartile = first_quartile_draws, asymptote = asymptote_draws)
    }

    # Step 2: turn each of the n_draws parameter sets into a Weibull curve and plot them
    plot_curves <- function(draws, color, title) {
      valid <- mapply(validate_weibull_parameters, draws$first_quartile, draws$median, draws$threshold, draws$asymptote)
      weibull_params <- calculate_weibull_parameters(draws$median, draws$first_quartile, draws$threshold)
      alphas <- ifelse(valid, weibull_params$alpha, NA_real_)
      betas <- ifelse(valid, weibull_params$beta, NA_real_)

      dist_matrix <- matrix(unlist(mapply(weibull_values, alphas, betas, draws$threshold, draws$asymptote, SIMPLIFY = FALSE)),
                            nrow = length(x_values), byrow = FALSE)

      plot(x_values, dist_matrix[, 1],
           type = "n",
           ylim = c(0, 1),
           xlab = "Age", ylab = "Cumulative Penetrance",
           main = title
      )
      for (i in seq_len(ncol(dist_matrix))) {
        lines(x_values, dist_matrix[, i], col = adjustcolor(color, alpha.f = 0.15))
      }

      mean_curve <- rowMeans(dist_matrix, na.rm = TRUE)
      lines(x_values, mean_curve, col = color, lwd = 2)

      legend("topleft", legend = title, col = color, lty = 1, cex = 0.8)
    }

    if (sex_specific) {
      # Separate images for visual clarity -- overlaying both sexes on one
      # plot makes it hard to judge either curve family on its own.
      plot_curves(draw_params("male"), "blue", "Male")
      plot_curves(draw_params("female"), "red", "Female")
    } else {
      plot_curves(draw_params("noSex"), "green", "Overall")
    }
  }

  return(prior_params)
}