# ----------------------------------------
# penetrance_postprocess.R
# ----------------------------------------
# Utilities to post-process posterior samples produced by the penetrance CRAN package.
#
# These functions convert sex-specific posterior samples (combined across chains)
# into age-specific summaries (posterior mean + credible intervals), which can then
# be plotted with full flexibility (ggplot2 / base R / custom themes).
#
# Expected input structure:
# - combined_chains is typically the output of penetrance and should contain (sex-specific) posterior vectors:
#   - median_male_results, first_quartile_male_results, threshold_male_results, asymptote_male_results
#   - median_female_results, first_quartile_female_results, threshold_female_results, asymptote_female_results
#
# Notes:
# - prob controls the credible interval width (e.g. prob = 0.95 gives 95% CrI).
# - min_age and max_age define the evaluated age grid: seq(min_age, max_age).
# - If you want to apply burn-in after estimation, you can run penetrance with burn_in = 0
#   and apply burn-in manually before calling these functions (e.g., subset posterior vectors).
# ----------------------------------------


# ----------------------------------------
# penetrance_distribution_cum()
# ----------------------------------------
# Compute cumulative penetrance (CDF) summaries from posterior samples.
#
#' @param combined_chains A list of posterior samples combined across chains (sex-specific) obtained from the penetrance function output.
#' @param prob Numeric, probability level for credible intervals (0 < prob < 1).
#' @param min_age Integer, minimum age to evaluate.
#' @param max_age Integer, maximum age to evaluate.
#' @param sex Character, "Male" or "Female".
#'
#' @return A list with:
#'   - age: integer vector seq(min_age, max_age)
#'   - mean_density: posterior mean cumulative penetrance at each age
#'   - ci_lower: lower credible interval bound at each age
#'   - ci_upper: upper credible interval bound at each age
# ----------------------------------------

penetrance_distribution_cum <- function(combined_chains, prob, min_age, max_age, sex) {
  
  if (!sex %in% c("Male","Female")) stop("sex must be 'Male' or 'Female'")
  if (!is.numeric(prob) || prob <= 0 || prob >= 1) stop("prob must be between 0 and 1")
  if (!is.numeric(min_age) || !is.numeric(max_age) || min_age > max_age) stop("min_age must be <= max_age")
  
  if (sex == "Male") {
    params <- calculate_weibull_parameters(
      combined_chains$median_male_results,
      combined_chains$first_quartile_male_results,
      combined_chains$threshold_male_results
    )
    alphas <- params$alpha
    betas <- params$beta
    thresholds <- combined_chains$threshold_male_results
    asymptotes <- combined_chains$asymptote_male_results
  } else {
    params <- calculate_weibull_parameters(
      combined_chains$median_female_results,
      combined_chains$first_quartile_female_results,
      combined_chains$threshold_female_results
    )
    alphas <- params$alpha
    betas <- params$beta
    thresholds <- combined_chains$threshold_female_results
    asymptotes <- combined_chains$asymptote_female_results
  }
  
  x_values <- seq(min_age, max_age)
  
  distributions <- mapply(function(alpha, beta, threshold, asymptote) {
    pweibull(x_values - threshold, shape = alpha, scale = beta) * asymptote
  }, alphas, betas, thresholds, asymptotes, SIMPLIFY = FALSE)
  
  distributions_matrix <- matrix(unlist(distributions), nrow = length(x_values), byrow = FALSE)
  mean_density <- rowMeans(distributions_matrix, na.rm = TRUE)
  ci_lower <- apply(distributions_matrix, 1, quantile, probs = (1 - prob) / 2, na.rm = TRUE)
  ci_upper <- apply(distributions_matrix, 1, quantile, probs = 1 - (1 - prob) / 2, na.rm = TRUE)
  
  return(list(
    age = x_values,
    mean_density = mean_density,
    ci_lower = ci_lower,
    ci_upper = ci_upper
  ))
}

# ----------------------------------------
# penetrance_distribution_den()
# ----------------------------------------
# Compute the density function (PDF) summaries from posterior samples.
#
#' @param combined_chains A list of posterior samples combined across chains (sex-specific) obtained from the penetrance function output.
#' @param prob Numeric, probability level for credible intervals (0 < prob < 1).
#' @param min_age Integer, minimum age to evaluate.
#' @param max_age Integer, maximum age to evaluate.
#' @param sex Character, "Male" or "Female".
#'
#' @return A list with:
#'   - age: integer vector seq(min_age, max_age)
#'   - mean_density: posterior mean penetrance density at each age
#'   - ci_lower: lower credible interval bound at each age
#'   - ci_upper: upper credible interval bound at each age
# ----------------------------------------

penetrance_distribution_den <- function(combined_chains, prob, min_age, max_age, 
                                        sex) {
  
  if (!sex %in% c("Male","Female")) stop("sex must be 'Male' or 'Female'")
  if (!is.numeric(prob) || prob <= 0 || prob >= 1) stop("prob must be between 0 and 1")
  
  if (sex == "Male") {
    params <- calculate_weibull_parameters(
      combined_chains$median_male_results,
      combined_chains$first_quartile_male_results,
      combined_chains$threshold_male_results
    )
    alphas <- params$alpha
    betas <- params$beta
    thresholds <- combined_chains$threshold_male_results
    asymptotes <- combined_chains$asymptote_male_results
  } else {
    params <- calculate_weibull_parameters(
      combined_chains$median_female_results,
      combined_chains$first_quartile_female_results,
      combined_chains$threshold_female_results
    )
    alphas <- params$alpha
    betas <- params$beta
    thresholds <- combined_chains$threshold_female_results
    asymptotes <- combined_chains$asymptote_female_results
  }
  
  x_values <- seq(min_age, max_age)
  
  pdf_distributions <- mapply(function(alpha, beta, threshold, asymptote) {
    dweibull(x_values - threshold, shape = alpha, scale = beta) * asymptote
  }, alphas, betas, thresholds, asymptotes, SIMPLIFY = FALSE)
  
  pdf_matrix <- matrix(unlist(pdf_distributions), nrow = length(x_values), byrow = FALSE)
  
  mean_density <- rowMeans(pdf_matrix, na.rm = TRUE)
  ci_lower <- apply(pdf_matrix, 1, quantile, probs = (1 - prob) / 2, na.rm = TRUE)
  ci_upper <- apply(pdf_matrix, 1, quantile, probs = 1 - (1 - prob) / 2, na.rm = TRUE)
  
  return(list(
    age = x_values,
    mean_density = mean_density,
    ci_lower = ci_lower,
    ci_upper = ci_upper
  ))
  
}