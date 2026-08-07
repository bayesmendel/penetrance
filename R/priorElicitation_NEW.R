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
#' @param data A data frame containing age and risk data. If NULL or entirely NA, default parameters are used. If provided, all entries in the age column must be numeric and non-NA.
#' @param sample_size Numeric, the total sample size used for risk proportion calculations.
#' @param ratio List, with components \code{noSex}, \code{male}, and \code{female}, each the odds ratio (OR) or relative risk (RR) used in asymptote parameter calculations for the corresponding baseline. Only the component(s) matching the inferred sex-specific mode need be supplied.
#' @param ratio_is_or Logical, indicating whether the ratios supplied are odds ratios (OR) (TRUE) or relative risks (FALSE).
#' @param ratio_ci_lower List, with components \code{noSex}, \code{male}, and \code{female}, the lower bound of each reported ratio's confidence interval (on the same OR/RR scale as \code{ratio}, per \code{ratio_is_or}). Required alongside the matching \code{ratio} component(s).
#' @param ratio_ci_upper List, with components \code{noSex}, \code{male}, and \code{female}, the upper bound of each reported ratio's confidence interval (on the same OR/RR scale as \code{ratio}, per \code{ratio_is_or}). Required alongside the matching \code{ratio} component(s).
#' @param prior_params List, containing prior parameters for the beta distributions. If NULL, default parameters are used.
#' @param risk_proportion Data frame, with default proportions of people at risk.
#' @param baseline_data Data frame with the baseline risk data.
#' @param max_age Numeric, the maximum age to plot the sampled curves over when \code{prior_predictive_check} is \code{TRUE}. Default is 94, matching \code{penetrance()}'s default \code{max_age}.
#' @param prior_predictive_check Logical, indicating whether to draw \code{n_draws} samples from \code{prior_params} and plot the resulting penetrance curves overlaid, as a sanity check on the joint prior specification. Default is FALSE.
#' @param n_draws Numeric, the number of curves to draw and overlay when \code{prior_predictive_check} is \code{TRUE}. Default is 50.
#' @param sex_specific Logical, indicating whether the prior predictive check should draw and overlay male/female curves separately (\code{TRUE}) or a single non-sex-specific set of curves (\code{FALSE}). Default is FALSE.
#'
#' @details
#' TODO: Update accordingly
#'
#' @return A list of functions representing the prior distributions for asymptote, threshold, median, and first quartile.
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
                      lifetime_risk = NULL, max_age = 94, prior_predictive_check = FALSE, n_draws = 50, sex_specific = FALSE) {

  # No data or ratio supplied: fall back to the default priors as-is.
  if ((is.null(data) || all(is.na(data))) && all(sapply(ratio, is.null))) {
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

  # Falls back to the shared (bare) prior when a sex/mode-specific one hasn't been set
  get_param <- function(param_name, suffix) {
    suffixed_param <- prior_params[[paste0(param_name, suffix)]]
    if (!is.null(suffixed_param)) return(suffixed_param)
    prior_params[[param_name]]
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