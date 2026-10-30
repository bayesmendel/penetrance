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
#' Make Priors
#'
#' This function generates prior distributions based on user input or default parameters.
#' It is designed to aid in the statistical analysis of risk proportions in populations, particularly in the context of cancer research.
#' The distributions are calculated for various statistical metrics such as asymptote, threshold, median, and first quartile.
#'
#' @param data A data frame with row names min, first_quartile, median, and max, and columns age and at_risk (a count of people at risk, not a proportion). If NULL or entirely NA, default parameters are used. If provided, all entries in the age and at_risk columns must be numeric and non-NA, ages must strictly increase from min to max, and no at_risk count may exceed sample_size.
#' @param sample_size Numeric, positive, the total number of people in the source study. Required when data is provided, since the asymptote is estimated from the fraction of sample_size no longer at risk by the median age.
#' @param ratio List, with components noSex, male, and female, each the odds ratio (OR) or relative risk (RR) used in asymptote parameter calculations for the corresponding baseline. Only the component(s) matching the inferred sex-specific mode need be supplied.
#' @param ratio_is_or Logical, indicating whether the ratios supplied are odds ratios (OR) (TRUE) or relative risks (FALSE).
#' @param ratio_ci_lower List, with components noSex, male, and female, the lower bound of each reported ratio's confidence interval (on the same OR/RR scale as ratio, per ratio_is_or). Required alongside the matching ratio component(s).
#' @param ratio_ci_upper List, with components noSex, male, and female, the upper bound of each reported ratio's confidence interval (on the same OR/RR scale as ratio, per ratio_is_or). Required alongside the matching ratio component(s).
#' @param prior_params List, containing prior parameters for the beta distributions. If NULL, default parameters are used.
#' @param baseline_data Data frame with the baseline risk data.
#' @param lifetime_risk Numeric between 0 and 1, an externally reported lifetime risk used as the asymptote prior's mean when asymptote_method = "lifetime_risk".
#' @param lifetime_risk_cv Numeric, positive. The assumed coefficient of variation (relative SD) of lifetime_risk, used to derive the asymptote prior's concentration when asymptote_method = "lifetime_risk". Default is 0.15.
#' @param asymptote_method Character, one of "default", "lifetime_risk", or "estimate". Selects how the asymptote prior is built for the time interval risk data pathway. "default" uses a flat Beta(1,1). "lifetime_risk" uses the user-supplied lifetime_risk as the Beta mean. "estimate" derives the asymptote from the time interval (KM-curve-style) risk data itself -- not yet implemented. NULL (default) leaves prior_params$asymptote untouched.
#' @param max_age Numeric, the maximum age to plot the sampled curves over when prior_predictive_check is TRUE. Default is 94, matching penetrance()'s default max_age.
#' @param prior_predictive_check Logical, indicating whether to draw n_draws samples from prior_params and plot the resulting penetrance curves overlaid, as a sanity check on the joint prior specification. Default is FALSE.
#' @param n_draws Numeric, the number of curves to draw and overlay when prior_predictive_check is TRUE. Default is 50.
#' @param sex_specific Logical, indicating whether the prior predictive check should draw and overlay male/female curves separately (TRUE) or a single non-sex-specific set of curves (FALSE). Default is FALSE.
#'
#' @details
#' If data is provided (the quartile pathway), the median and first quartile priors are Beta distributions whose means are the reported ages rescaled onto [0, 1] and whose concentrations are the at_risk counts at those ages. The threshold prior is uniform from 10 years below to 10 years above the reported minimum age, with the lower bound floored at 0.
#'
#' The asymptote prior is a Beta distribution with mean 2 * (1 - at_risk_median / sample_size) and concentration at_risk_median. This follows from the reported median age being the age by which half of the eventually affected carriers are diagnosed, so the fraction diagnosed by that age is half the asymptote. If this estimate falls outside (0, 1), for example because follow-up in the source study was incomplete or heavily censored, an error is thrown.
#'
#' If asymptote_method or ratio is also supplied, the asymptote prior from those pathways overrides the one derived from data.
#'
#' @return A list of the prior parameters for asymptote, threshold, median, and first_quartile, to be passed as prior_params to penetrance().
#'
#' @seealso \code{\link{qbeta}}, \code{\link{runif}}
#' 
#' @export
makePriors <- function(data,
                      sample_size,
                      ratio = list(noSex = NULL, male = NULL, female = NULL), ratio_is_or,
                      ratio_ci_lower = list(noSex = NULL, male = NULL, female = NULL),
                      ratio_ci_upper = list(noSex = NULL, male = NULL, female = NULL),
                      prior_params, baseline_data,
                      asymptote_method = NULL,
                      lifetime_risk = NULL, lifetime_risk_cv = 0.15,
                      max_age = 94, prior_predictive_check = FALSE, n_draws = 50, sex_specific = FALSE) {

  # No data or ratio supplied: fall back to the default priors as-is.
  if ((is.null(data) || all(is.na(data))) && all(sapply(ratio, is.null)) && is.null(asymptote_method)) {
    return(prior_params_default)
  }

  # Converts OR to RR given the baseline risk P2:
  # OR = [P1/(1-P1)] / [P2/(1-P2)], solved for P1 and substituted into RR = P1/P2.
  # P1 = OR·P2 / (1 - P2 + OR·P2); where P1 = asymptote, P2 = baseline lifetime risk

  convert_or_to_rr <- function(or, baseline_risk) {
    return(or / (1 - baseline_risk + (or * baseline_risk)))
  }

  # Standard error of ratio on the raw scale, valid when the reported CI is a symmetric, raw-scale Wald interval.
  compute_se_raw <- function(ratio_ci_lower, ratio_ci_upper) {
    return((ratio_ci_upper - ratio_ci_lower) / (2 * 1.96))
  }

  # Standard error of log(ratio), valid when the reported CI is built on the log scale and are asymmetric on the raw scale.
  compute_se_log <- function(ratio_ci_lower, ratio_ci_upper) {
    return((log(ratio_ci_upper) - log(ratio_ci_lower)) / (2 * 1.96))
  }

  # Variance of ratio: a CI symmetric around ratio implies a raw-scale Wald interval (Var = SE_raw^2); an asymmetric CI implies a log-scale interval, propagated to the raw scale via the delta method (Var = ratio^2 * SE_log^2).
  compute_ratio_variance <- function(ratio, ratio_ci_lower, ratio_ci_upper) {
    is_symmetric <- isTRUE(all.equal(ratio - ratio_ci_lower, ratio_ci_upper - ratio))

    if (is_symmetric) {
      se_raw <- compute_se_raw(ratio_ci_lower, ratio_ci_upper)
      return(se_raw^2)
    }

    se_log <- compute_se_log(ratio_ci_lower, ratio_ci_upper)
    return(ratio^2 * se_log^2)
  }

  # Propagates ratio's variance to the asymptote scale, treating baseline as a fixed constant: Var(baseline * ratio) = baseline^2 * Var(ratio).
  compute_asymptote_variance <- function(baseline, ratio_variance) {
    return(baseline^2 * ratio_variance)
  }

  # Matches a Beta(g1, g2) to a target mean and variance, solving
  # Var = mean(1-mean)/(concentration+1) for concentration.
  compute_concentration <- function(mean, variance) {
    return((mean * (1 - mean)) / variance - 1)
  }

  normalize_median <- function(stat, min_age, max_age) {
    return((stat - min_age) / (max_age - min_age))
  }

  normalize_first_quartile <- function(stat, min_age, median_age) {
    return((stat - min_age) / (median_age - min_age))
  }

  compute_parameters_median <- function(stat, at_risk, min_age, max_age) {
    median_norm <- normalize_median(stat, min_age, max_age)
    alpha <- median_norm * at_risk
    beta <- at_risk - alpha
    return(list(m1 = alpha, m2 = beta))
  }

  compute_parameters_quartile <- function(stat, at_risk, min_age, median_age) {
    quartile_norm <- normalize_first_quartile(stat, min_age, median_age)
    alpha <- quartile_norm * at_risk
    beta <- at_risk - alpha
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
    alpha <- asymptote_estimate * at_risk
    beta <- at_risk - alpha
    return(list(g1 = alpha, g2 = beta))
  }

  # Falls back to the shared (bare) prior when a sex/mode-specific one hasn't been set
  get_param <- function(param_name, suffix) {
    suffixed_param <- prior_params[[paste0(param_name, suffix)]]
    if (!is.null(suffixed_param)) return(suffixed_param)
    prior_params[[param_name]]
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
  if (!(is.null(data) || all(is.na(data)))) {
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

    res_median <- compute_parameters_median(median_age, risk_median, min_age, data["max", "age"])
    res_first_quartile <- compute_parameters_quartile(first_quartile_age, risk_first_quartile, min_age, median_age)
    res_asymptote <- compute_parameters_asymptote(risk_median, sample_size)

    prior_params <- list(
      asymptote = list(g1 = res_asymptote$g1, g2 = res_asymptote$g2),
      threshold = list(min = max(0, min_age - 10), max = min_age + 10),
      median = list(m1 = res_median$m1, m2 = res_median$m2),
      first_quartile = list(q1 = res_first_quartile$q1, q2 = res_first_quartile$q2)
    )
  }

  # Time interval risk data pathway
  if (!is.null(asymptote_method)) {
    prior_params$asymptote <- compute_asymptote_time_interval(asymptote_method, lifetime_risk, lifetime_risk_cv)
  }

  # Ratio pathway
  if (!all(sapply(ratio, is.null))) {
    # Compute beta parameters for each prior
    if (!is.null(ratio$male)) {
      SEER_lifetime_male <- sum(baseline_data$Male)
      if (ratio_is_or == TRUE) {
        ratio$male <- convert_or_to_rr(ratio$male, SEER_lifetime_male)
        ratio_ci_lower$male <- convert_or_to_rr(ratio_ci_lower$male, SEER_lifetime_male)
        ratio_ci_upper$male <- convert_or_to_rr(ratio_ci_upper$male, SEER_lifetime_male)
      }

      if (SEER_lifetime_male * ratio$male >= 1) {
        stop("Error: 'SEER_lifetime_male * ratio$male' (the implied male carrier lifetime risk) must be less than 1.")
      }

      mean_male <- SEER_lifetime_male * ratio$male
      ratio_variance_male <- compute_ratio_variance(ratio$male, ratio_ci_lower$male, ratio_ci_upper$male)
      asymptote_variance_male <- compute_asymptote_variance(SEER_lifetime_male, ratio_variance_male)
      concentration_male <- compute_concentration(mean_male, asymptote_variance_male)

      g1_male <- mean_male * concentration_male
      g2_male <- concentration_male - g1_male

      prior_params$asymptote_male <- list(g1 = g1_male, g2 = g2_male)
    }

    if (!is.null(ratio$female)) {
      SEER_lifetime_female <- sum(baseline_data$Female)
      if (ratio_is_or == TRUE) {
        ratio$female <- convert_or_to_rr(ratio$female, SEER_lifetime_female)
        ratio_ci_lower$female <- convert_or_to_rr(ratio_ci_lower$female, SEER_lifetime_female)
        ratio_ci_upper$female <- convert_or_to_rr(ratio_ci_upper$female, SEER_lifetime_female)
      }

      if (SEER_lifetime_female * ratio$female >= 1) {
        stop("Error: 'SEER_lifetime_female * ratio$female' (the implied female carrier lifetime risk) must be less than 1.")
      }

      mean_female <- SEER_lifetime_female * ratio$female
      ratio_variance_female <- compute_ratio_variance(ratio$female, ratio_ci_lower$female, ratio_ci_upper$female)
      asymptote_variance_female <- compute_asymptote_variance(SEER_lifetime_female, ratio_variance_female)
      concentration_female <- compute_concentration(mean_female, asymptote_variance_female)

      g1_female <- mean_female * concentration_female
      g2_female <- concentration_female - g1_female

      prior_params$asymptote_female <- list(g1 = g1_female, g2 = g2_female)
    }

    if (!is.null(ratio$noSex)) {
      SEER_lifetime_noSex <- sum(baseline_data)
      if (ratio_is_or == TRUE) {
        ratio$noSex <- convert_or_to_rr(ratio$noSex, SEER_lifetime_noSex)
        ratio_ci_lower$noSex <- convert_or_to_rr(ratio_ci_lower$noSex, SEER_lifetime_noSex)
        ratio_ci_upper$noSex <- convert_or_to_rr(ratio_ci_upper$noSex, SEER_lifetime_noSex)
      }

      if (SEER_lifetime_noSex * ratio$noSex >= 1) {
        stop("Error: 'SEER_lifetime_noSex * ratio$noSex' (the implied carrier lifetime risk) must be less than 1.")
      }

      mean_noSex <- SEER_lifetime_noSex * ratio$noSex
      ratio_variance_noSex <- compute_ratio_variance(ratio$noSex, ratio_ci_lower$noSex, ratio_ci_upper$noSex)
      asymptote_variance_noSex <- compute_asymptote_variance(SEER_lifetime_noSex, ratio_variance_noSex)
      concentration_noSex <- compute_concentration(mean_noSex, asymptote_variance_noSex)

      g1_noSex <- mean_noSex * concentration_noSex
      g2_noSex <- concentration_noSex - g1_noSex

      prior_params$asymptote_noSex <- list(g1 = g1_noSex, g2 = g2_noSex)
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
    draw_params <- function(suffix) {
      threshold_params <- get_param("threshold", suffix)
      median_params <- get_param("median", suffix)
      first_quartile_params <- get_param("first_quartile", suffix)
      asymptote_params <- get_param("asymptote", suffix)

      threshold_draws <- runif(n_draws, threshold_params$min, threshold_params$max)
      scaled_median_draws <- qbeta(runif(n_draws), median_params$m1, median_params$m2)
      scaled_first_quartile_draws <- qbeta(runif(n_draws), first_quartile_params$q1, first_quartile_params$q2)
      asymptote_draws <- qbeta(runif(n_draws), asymptote_params$g1, asymptote_params$g2)

      # Un-normalize onto actual ages (mirrors the scaling relations in mhChain.R)
      median_draws <- threshold_draws + scaled_median_draws * (max_age - threshold_draws)
      first_quartile_draws <- threshold_draws + scaled_first_quartile_draws * (median_draws - threshold_draws)

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
      plot_curves(draw_params("_male"), "blue", "Male")
      plot_curves(draw_params("_female"), "red", "Female")
    } else {
      plot_curves(draw_params("_noSex"), "green", "Overall")
    }
  }

  return(prior_params)
}