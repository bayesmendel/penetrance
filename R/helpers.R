#' Calculate Weibull Parameters
#'
#' This function calculates the shape (\code{alpha}) and scale (\code{beta}) parameters
#' of a Weibull distribution given the median, first quartile, and delta values.
#'
#' @param given_median The median of the data.
#' @param given_first_quartile The first quartile of the data.
#' @param delta A constant offset value.
#'
#' @return A list containing the calculated Weibull parameters:
#'   \item{alpha}{The shape parameter of the Weibull distribution}
#'   \item{beta}{The scale parameter of the Weibull distribution}
#' 
#' @examples
#' # Calculate Weibull parameters
#' params <- calculate_weibull_parameters(
#'   given_median = 50,
#'   given_first_quartile = 30,
#'   delta = 15
#' )
#' print(params)
#' @export
calculate_weibull_parameters <- function(given_median, given_first_quartile, delta) {
  # Calculate alpha
  alpha <- log(-log(0.5) / -log(0.75)) / log((given_median - delta) / (given_first_quartile - delta))
  
  # Calculate beta using the median (M)
  beta <- (given_median - delta) / (-log(0.5))^(1 / alpha)
  
  return(list(alpha = alpha, beta = beta))
}

#' Validate Weibull Parameters
#'
#' This function validates the given parameters for calculating Weibull distribution.
#'
#' @param given_first_quartile The first quartile of the data.
#' @param given_median The median of the data.
#' @param threshold A constant threshold value.
#' @param asymptote A constant asymptote value (gamma).
#'
#' @return Logical value indicating whether the parameters are valid (TRUE) or not (FALSE)
#' 
#' @examples
#' # Validate parameters
#' is_valid <- validate_weibull_parameters(
#'   given_first_quartile = 30,
#'   given_median = 50,
#'   threshold = 15,
#'   asymptote = 0.8
#' )
#' print(is_valid)
#' @export
validate_weibull_parameters <- function(given_first_quartile, given_median, threshold, asymptote) {
  # Handle NA values
  if (is.na(given_median) || is.na(given_first_quartile) || is.na(threshold) || is.na(asymptote)) {
    return(FALSE)
  }
  
  # Check for negative or zero values
  if (given_median <= 0 || given_first_quartile <= 0 || threshold < 0) {
    return(FALSE)
  }
  
  # Check if asymptote (gamma) is within the valid range (0,1)
  if (asymptote <= 0 || asymptote >= 1) {
    return(FALSE)
  }
  
  # Check if the logarithmic calculations will be valid
  if (given_first_quartile <= threshold || given_median <= threshold) {
    return(FALSE)
  }
  
  # Check if the denominator in the alpha calculation would be zero
  if ((given_first_quartile - threshold) == (given_median - threshold)) {
    return(FALSE)
  }
  
  # Check ordering
  if (given_first_quartile >= given_median) {
    return(FALSE)
  }
  
  # If all checks pass, return TRUE
  return(TRUE)
}

#' Transform Data Frame
#'
#' This function transforms a data frame from the standard format used in PanelPRO
#' into the required format which conforms to the requirements of penetrance (and clipp).
#'
#' @param df The input data frame in the usual PanelPRO format.
#'
#' @return A data frame in the format required for clipp with the following columns:
#'   \item{individual}{ID of the individual}
#'   \item{isProband}{Indicator if the individual is a proband}
#'   \item{family}{Family ID}
#'   \item{mother}{Mother's ID}
#'   \item{father}{Father's ID}
#'   \item{aff}{Affection status}
#'   \item{sex}{Sex (2 for female, 1 for male)}
#'   \item{age}{Age at diagnosis or current age}
#'   \item{geno}{Genotype information}
#' 
#' @examples
#' # Create example data frame
#' df <- data.frame(
#'   ID = 1:2,
#'   PedigreeID = c(1,1),
#'   Sex = c(0,1),
#'   MotherID = c(NA,1),
#'   FatherID = c(NA,NA),
#'   isProband = c(1,0),
#'   CurAge = c(45,20),
#'   isAff = c(1,0),
#'   Age = c(40,NA),
#'   Geno = c(1,0)
#' )
#' 
#' # Transform the data frame
#' transformed_df <- transformDF(df)
#' @export
transformDF <- function(df) {
  # Rename and transform columns
  df$individual <- df$ID
  df$isProband <- df$isProband
  df$family <- df$PedigreeID
  df$mother <- df$MotherID
  df$father <- df$FatherID
  df$aff <- df$isAff
  df$sex <- ifelse(df$Sex == 0, 2, df$Sex) # Convert 0s to 2s (female), keep 1s as is (male)
  
  # Apply row-wise logic to assign age based on aff column
  df$age <- apply(df, 1, function(row) {
    aff <- as.numeric(row[["isAff"]])
    if (aff == 1) {
      return(as.numeric(row[["Age"]]))
    } else {
      return(as.numeric(row[["CurAge"]]))
    }
  })
  
  # Process 'geno' column, checking for both 'geno' and 'Geno'
  if ("Geno" %in% colnames(df)) {
    df$geno <- ifelse(is.na(df$Geno), "", ifelse(df$Geno == 1, "1/2", ifelse(df$Geno == 0, "1/1", df$Geno)))
  } else {
    # Handle case where neither geno nor Geno column is present, maybe default or error
    warning("Neither 'geno' nor 'Geno' column found in the input data frame. 'geno' column in output will be empty.")
    df$geno <- ""
  }
  
  # Select only the necessary columns
  df <- df[c("individual", "isProband", "family", "mother", "father", "aff", "sex", "age", "geno")]
  
  return(df)
}

#' Validate Allele Frequency Input
#'
#' This function validates whether the provided input is a valid allele frequency.
#' Allele frequencies must be numeric, scalar, and within the range [0, 1].
#' The function also provides warnings for edge cases and common mistakes.
#'
#' @param allele_freq The allele frequency value to validate. Should be a numeric value between 0 and 1.
#' @param param_name Character string specifying the parameter name (for error messages). Default is "allele_freq".
#' @param warn_threshold Numeric value above which to issue a warning about unusually high allele frequency. Default is 0.01 (1%).
#'
#' @return Logical value TRUE if the allele frequency is valid, otherwise stops with an error message.
#'
#' @details
#' The function checks:
#' \itemize{
#'   \item{Whether the input is numeric}
#'   \item{Whether the input is a single value (not a vector)}
#'   \item{Whether the value is between 0 and 1 (inclusive)}
#'   \item{Whether the value is unusually high (> warn_threshold), which may indicate the user provided carrier prevalence instead of allele frequency}
#'   \item{Whether the value is exactly 0 or 1, which may not be biologically meaningful}
#' }
#'
#' @examples
#' # Valid allele frequencies
#' validate_allele_freq(0.0001)  # Common for rare variants
#' validate_allele_freq(0.001)
#' validate_allele_freq(0.05)
#'
#' \dontrun{
#' # Invalid inputs (will throw errors)
#' validate_allele_freq("0.001")       # Not numeric
#' validate_allele_freq(c(0.001, 0.002))  # Vector instead of scalar
#' validate_allele_freq(-0.001)       # Negative value
#' validate_allele_freq(1.5)          # Greater than 1
#' validate_allele_freq(NA)           # Missing value
#' 
#' # Valid but will trigger warnings
#' validate_allele_freq(0.02)         # Unusually high (>1%), warning issued
#' validate_allele_freq(0)            # Edge case, warning issued
#' validate_allele_freq(1)            # Edge case, warning issued
#' }
#'
#' @export
validate_allele_freq <- function(allele_freq, param_name = "allele_freq", warn_threshold = 0.01) {
  # Check if the input is missing
  if (missing(allele_freq)) {
    stop(paste0("Error: '", param_name, "' is required but was not provided."))
  }
  
  # Check if the input is NA
  if (any(is.na(allele_freq))) {
    stop(paste0("Error: '", param_name, "' cannot be NA or contain NA values."))
  }
  
  # Check if the input is numeric
  if (!is.numeric(allele_freq)) {
    stop(paste0("Error: '", param_name, "' must be numeric. ",
                "Received type: ", class(allele_freq)[1], ". ",
                "Please provide a numeric value between 0 and 1."))
  }
  
  # Check if the input is a single value (not a vector)
  if (length(allele_freq) != 1) {
    stop(paste0("Error: '", param_name, "' must be a single numeric value, not a vector. ",
                "Received ", length(allele_freq), " values. ",
                "Please provide a single allele frequency."))
  }
  
  # Check if the value is within valid range [0, 1]
  if (allele_freq < 0 || allele_freq > 1) {
    stop(paste0("Error: '", param_name, "' must be between 0 and 1 (inclusive). ",
                "Received: ", allele_freq, ". ",
                "Allele frequencies represent proportions and must be in the range [0, 1]."))
  }
  
  # Check for infinite values
  if (is.infinite(allele_freq)) {
    stop(paste0("Error: '", param_name, "' cannot be infinite."))
  }
  
  # Warning for edge case: exactly 0
  if (allele_freq == 0) {
    warning(paste0("Warning: '", param_name, "' is exactly 0. ",
                   "This implies the allele does not exist in the population, which may not be biologically meaningful. ",
                   "Consider if this is the intended value."))
  }
  
  # Warning for edge case: exactly 1
  if (allele_freq == 1) {
    warning(paste0("Warning: '", param_name, "' is exactly 1 (100%). ",
                   "This implies the allele is fixed in the population (all individuals carry it), which may not be biologically meaningful. ",
                   "Consider if this is the intended value."))
  }
  
  # Warning for unusually high allele frequency
  if (allele_freq > warn_threshold && allele_freq < 1) {
    warning(paste0("Warning: '", param_name, "' is ", allele_freq, " (", allele_freq * 100, "%), which is relatively high for a disease-associated variant. ",
                   "Please verify this is the allele frequency (p) and not the carrier prevalence (approximately 2p). ",
                   "For example, if the carrier prevalence is 2%, the allele frequency should be approximately 1% (0.01)."))
  }
  
  return(TRUE)
}

#' Validate Baseline Penetrance Data
#'
#' This function validates baseline penetrance data to ensure it represents age-specific
#' probabilities rather than cumulative risk. It checks for monotonicity and whether
#' the sum exceeds 1, both of which suggest the data may be cumulative rather than
#' age-specific.
#'
#' @param baseline_data The baseline data to validate. Can be:
#'   - A data frame with 'Male' and 'Female' columns (when sex_specific = TRUE)
#'   - A numeric vector (when sex_specific = FALSE)
#'   - A single-column data frame (when sex_specific = FALSE)
#' @param sex_specific Logical, indicating whether the data is sex-specific. Default is TRUE.
#' @param param_name Character string specifying the parameter name (for messages). Default is "baseline_data".
#' @param tolerance Numeric value for checking strict monotonicity (to account for floating point precision). Default is 1e-10.
#'
#' @return Logical value TRUE if validation passes (with possible warnings), otherwise stops with an error.
#'
#' @details
#' The function performs the following checks:
#' \itemize{
#'   \item{For monotonicity: If values are strictly non-decreasing (monotone increasing), this suggests
#'         cumulative risk rather than age-specific probabilities. A warning is issued.}
#'   \item{For sum > 1: If the sum of all probabilities exceeds 1, this is problematic because
#'         these should be age-specific probabilities. A warning is issued.}
#'   \item{Individual values must be between 0 and 1 (probabilities)}
#'   \item{No NA or infinite values are allowed}
#' }
#'
#' Age-specific baseline risk represents the probability of developing disease at each specific age,
#' while cumulative risk represents the total probability up to that age. For proper penetrance
#' estimation, age-specific (not cumulative) risk should be used.
#'
#' @examples
#' # Valid age-specific data (varies, not monotone)
#' age_specific <- c(0.001, 0.002, 0.003, 0.002, 0.004, 0.003, 0.005)
#' validate_baseline_data(age_specific, sex_specific = FALSE)
#'
#' # Valid sex-specific data
#' baseline_df <- data.frame(
#'   Male = c(0.001, 0.002, 0.001, 0.003),
#'   Female = c(0.002, 0.003, 0.002, 0.004)
#' )
#' validate_baseline_data(baseline_df, sex_specific = TRUE)
#'
#' \dontrun{
#' # Will trigger warnings
#' # Monotone increasing (suggests cumulative risk)
#' cumulative <- c(0.001, 0.002, 0.003, 0.004, 0.005)
#' validate_baseline_data(cumulative, sex_specific = FALSE)
#'
#' # Sum greater than 1
#' high_values <- rep(0.1, 15)  # sum = 1.5
#' validate_baseline_data(high_values, sex_specific = FALSE)
#'
#' # Invalid data
#' invalid_data <- c(0.001, -0.002, 0.003)  # Negative value
#' validate_baseline_data(invalid_data, sex_specific = FALSE)
#' }
#'
#' @export
validate_baseline_data <- function(baseline_data, sex_specific = TRUE, 
                                   param_name = "baseline_data", 
                                   tolerance = 1e-10) {
  
  # Helper function to check if a vector is monotone increasing
  is_monotone_increasing <- function(x, tol = tolerance) {
    if (length(x) <= 1) return(FALSE)
    # Remove NAs for checking
    x_clean <- x[!is.na(x)]
    if (length(x_clean) <= 1) return(FALSE)
    
    # Check if strictly non-decreasing (allowing for small tolerance)
    diffs <- diff(x_clean)
    return(all(diffs >= -tol))
  }
  
  # Helper function to validate a numeric vector of probabilities
  validate_probability_vector <- function(vec, vec_name = "data") {
    # Check for NA values
    if (any(is.na(vec))) {
      stop(paste0("Error: '", param_name, "' (", vec_name, ") contains NA values. ",
                  "All baseline probabilities must be specified."))
    }
    
    # Check for infinite values
    if (any(is.infinite(vec))) {
      stop(paste0("Error: '", param_name, "' (", vec_name, ") contains infinite values."))
    }
    
    # Check if all values are between 0 and 1
    if (any(vec < 0) || any(vec > 1)) {
      stop(paste0("Error: '", param_name, "' (", vec_name, ") must contain only values between 0 and 1. ",
                  "Baseline data should represent probabilities (age-specific risk)."))
    }
    
    # Check if monotone increasing (suggests cumulative risk)
    if (is_monotone_increasing(vec)) {
      warning(paste0("Warning: '", param_name, "' (", vec_name, ") appears to be monotone increasing. ",
                     "This suggests the data may represent CUMULATIVE risk rather than AGE-SPECIFIC risk. ",
                     "For penetrance estimation, please ensure you are using age-specific probabilities, ",
                     "not cumulative probabilities. Age-specific risk can fluctuate with age, ",
                     "while cumulative risk always increases."),
              immediate. = TRUE)
    }
    
    # Check if sum is greater than 1
    vec_sum <- sum(vec, na.rm = TRUE)
    if (vec_sum > 1) {
      warning(paste0("Warning: '", param_name, "' (", vec_name, ") has a sum of ", 
                     round(vec_sum, 4), " which is greater than 1. ",
                     "This strongly suggests the data may represent CUMULATIVE risk rather than AGE-SPECIFIC risk. ",
                     "For age-specific probabilities, the sum typically should not exceed 1. ",
                     "Please verify your data represents age-specific (not cumulative) baseline risk."),
              immediate. = TRUE)
    }
    
    return(TRUE)
  }
  
  # Main validation logic based on sex_specific
  if (sex_specific) {
    # Expect a data frame with Male and Female columns
    if (!is.data.frame(baseline_data)) {
      stop(paste0("Error: '", param_name, "' must be a data frame when sex_specific is TRUE."))
    }
    
    required_cols <- c("Male", "Female")
    if (!all(required_cols %in% colnames(baseline_data))) {
      stop(paste0("Error: '", param_name, "' must have columns named 'Male' and 'Female' when sex_specific is TRUE."))
    }
    
    # Validate each sex-specific column
    validate_probability_vector(baseline_data$Male, "Male")
    validate_probability_vector(baseline_data$Female, "Female")
    
  } else {
    # Not sex_specific: expect a vector or single-column data frame
    if (is.data.frame(baseline_data)) {
      if (ncol(baseline_data) != 1) {
        stop(paste0("Error: '", param_name, "' must be a single-column data frame when sex_specific is FALSE. ",
                    "Received ", ncol(baseline_data), " columns."))
      }
      vec <- baseline_data[[1]]
    } else if (is.vector(baseline_data) && is.numeric(baseline_data)) {
      vec <- baseline_data
    } else {
      stop(paste0("Error: '", param_name, "' must be a numeric vector or single-column data frame when sex_specific is FALSE."))
    }
    
    # Validate the vector
    validate_probability_vector(vec, "values")
  }
  
  return(TRUE)
}