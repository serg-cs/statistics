#' Sample Mean Distribution (Infinite Population / With Replacement)
#'
#' Calculates the density or cumulative probability for the sample mean,
#' assuming an infinite population or sampling with replacement.
#'
#' Under these conditions, the standard deviation of the sample mean (Standard Error)
#' is calculated as \eqn{\sigma / \sqrt{n}}.
#'
#' @param x The sample mean value to evaluate (quantile).
#' @param mu Population mean (numeric scalar).
#' @param sigma Population standard deviation (numeric scalar, must be > 0).
#' @param n Sample size (numeric scalar, must be > 0).
#' @param type Type of calculation: "density" for PDF height, "cumulative" for area (P(X <= x)).
#' @return Probability density or cumulative probability.
#' @export
get_sample_mean_infinite_population <- function(x, mu, sigma, n, type = c("density", "cumulative")) {
  # 1. Validation
  type <- match.arg(type)
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(mu), length(mu) == 1)
  stopifnot(is.numeric(sigma), length(sigma) == 1, sigma > 0)
  stopifnot(is.numeric(n), length(n) == 1, n > 0)
  
  # 2. Calculate Standard Error (SE)
  # For infinite population: SE = sigma / sqrt(n)
  std_error <- sigma / sqrt(n)
  
  # 3. Calculation using Normal Distribution with adjusted SD
  if (type == "density") {
    stats::dnorm(x, mean = mu, sd = std_error)
  } else {
    stats::pnorm(x, mean = mu, sd = std_error)
  }
}

#' Sample Mean Distribution (Finite Population / Without Replacement)
#'
#' Calculates the density or cumulative probability for the sample mean,
#' assuming a finite population of size N sampled without replacement.
#'
#' This applies the Finite Population Correction Factor (FPCF).
#' The Standard Error is calculated as: \eqn{\frac{\sigma}{\sqrt{n}} \times \sqrt{\frac{N-n}{N-1}}}.
#'
#' @param x The sample mean value to evaluate (quantile).
#' @param mu Population mean (numeric scalar).
#' @param sigma Population standard deviation (numeric scalar, must be > 0).
#' @param n Sample size (numeric scalar, must be > 0).
#' @param N Population size (numeric scalar, must be > n).
#' @param type Type of calculation: "density" for PDF height, "cumulative" for area (P(X <= x)).
#' @return Probability density or cumulative probability.
#' @export
get_sample_mean_finite_population <- function(x, mu, sigma, n, N, type = c("density", "cumulative")) {
  # 1. Validation
  type <- match.arg(type)
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(mu), length(mu) == 1)
  stopifnot(is.numeric(sigma), length(sigma) == 1, sigma > 0)
  stopifnot(is.numeric(n), length(n) == 1, n > 0)
  stopifnot(is.numeric(N), length(N) == 1)
  
  # Logical constraint: Sample size cannot exceed Population size in this context
  stopifnot(n <= N)
  
  # 2. Calculate Standard Error (SE) with Correction Factor
  # Base SE
  base_se <- sigma / sqrt(n)
  
  # Finite Population Correction Factor: sqrt((N-n)/(N-1))
  # If N is very large, this factor approaches 1.
  if (N > 1) {
    correction_factor <- sqrt((N - n) / (N - 1))
  } else {
    correction_factor <- 0 # Edge case handling if N=1 (though usually N >> n)
  }
  
  adjusted_se <- base_se * correction_factor
  
  # 3. Calculation
  if (type == "density") {
    stats::dnorm(x, mean = mu, sd = adjusted_se)
  } else {
    stats::pnorm(x, mean = mu, sd = adjusted_se)
  }
}

#' Sample Quasi-Variance Distribution (Chi-Square)
#'
#' Calculates the density or cumulative probability for the Sample Quasi-Variance ($S^2$).
#'
#' Unlike the sample mean, the sample variance does not follow a normal distribution.
#' Instead, the quantity \eqn{\frac{(n-1)S^2}{\sigma^2}} follows a Chi-Square distribution
#' with \eqn{n-1} degrees of freedom.
#'
#' This function automatically handles the transformation from the raw variance value ($S^2$)
#' to the Chi-Square statistic.
#'
#' @param x The sample quasi-variance value to evaluate (must be >= 0).
#' @param sigma Population standard deviation (numeric scalar, must be > 0).
#' @param n Sample size (numeric scalar, must be > 1).
#' @param type Type of calculation: "density" for PDF height, "cumulative" for area (P(S^2 <= x)).
#' @return Probability density or cumulative probability.
#' @export
get_sample_quasi_variance_distribution <- function(x, sigma, n, type = c("density", "cumulative")) {
  # 1. Validation
  type <- match.arg(type)
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(sigma), length(sigma) == 1, sigma > 0)
  stopifnot(is.numeric(n), length(n) == 1, n > 1)
  
  # Variance cannot be negative
  if (any(x < 0)) warning("Provided variance 'x' contains negative values, which are impossible.")
  
  # 2. Calculate Degrees of Freedom
  df <- n - 1
  
  # 3. Transform S^2 to Chi-Square Statistic
  # Formula: Chi2 = (n-1) * S^2 / sigma^2
  chi_sq_stat <- (df * x) / (sigma^2)
  
  # 4. Calculation
  if (type == "density") {
    # For density, we must apply the Chain Rule (Jacobian) because we are 
    # transforming the variable.
    # f_S2(s2) = f_Chi2(y) * |dy/ds2|
    # where y = (n-1)s2/sigma^2  =>  dy/ds2 = (n-1)/sigma^2
    
    jacobian <- df / (sigma^2)
    stats::dchisq(chi_sq_stat, df = df) * jacobian
    
  } else {
    # For Cumulative, P(S^2 <= x) is exactly P(Chi2 <= transformed_x)
    stats::pchisq(chi_sq_stat, df = df)
  }
}

#' Sample Proportion Distribution (Normal Approximation)
#'
#' Calculates the density or cumulative probability for the Sample Proportion ($\hat{p}$),
#' using the Normal Approximation.
#'
#' This approximation is valid when the sample size is large (typically n > 30).
#' The standard error is calculated as \eqn{\sqrt{\frac{p(1-p)}{n}}}.
#'
#' @param x The sample proportion value to evaluate (between 0 and 1).
#' @param p Population proportion (numeric scalar between 0 and 1).
#' @param n Sample size (numeric scalar, must be > 0).
#' @param type Type of calculation: "density" for PDF height, "cumulative" for area (P(p_hat <= x)).
#' @return Probability density or cumulative probability.
#' @export
get_sample_proportion_distribution <- function(x, p, n, type = c("density", "cumulative")) {
  # 1. Validation
  type <- match.arg(type)
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(p), length(p) == 1, p >= 0, p <= 1)
  stopifnot(is.numeric(n), length(n) == 1, n > 0)
  
  # Warning for small sample sizes where Normal Approx might fail
  if (n <= 30) warning("Sample size n <= 30. Normal approximation for proportions may be inaccurate.")
  
  # 2. Calculate Standard Error (SE)
  # Formula: sqrt( p * (1-p) / n )
  std_error <- sqrt((p * (1 - p)) / n)
  
  # 3. Calculation using Normal Distribution
  if (type == "density") {
    stats::dnorm(x, mean = p, sd = std_error)
  } else {
    stats::pnorm(x, mean = p, sd = std_error)
  }
}
#' Variance of the Sample Mean
#'
#' Calculates the variance of the sampling distribution of the sample mean.
#' This value represents how much the sample mean (\eqn{\bar{X}}) is expected to vary
#' from sample to sample.
#'
#' This is NOT the variance of the population ($\sigma^2$), but the variance of the
#' average of \eqn{n} observations: \eqn{Var(\bar{X}) = \frac{\sigma^2}{n}}.
#'
#' If a population size \eqn{N} is provided (finite population), the Finite Population
#' Correction Factor (FPCF) is applied: \eqn{\frac{N-n}{N-1}}.
#'
#' @param sigma_sq Population variance (\eqn{\sigma^2}). Numeric scalar > 0.
#'                 (If you only have standard deviation, square it: sigma^2).
#' @param n Sample size. Numeric scalar > 0.
#' @param N Population size (Optional). If provided, applies the Finite Population Correction.
#'          Default is Inf (infinite population).
#' @return The variance of the sample mean (numeric scalar).
#' @export
get_variance_of_sample_mean <- function(sigma_sq, n, N = Inf) {
  # 1. Validation
  stopifnot(is.numeric(sigma_sq), length(sigma_sq) == 1, sigma_sq > 0)
  stopifnot(is.numeric(n), length(n) == 1, n > 0)
  
  # 2. Base Calculation (Infinite Population)
  # Var(X_bar) = sigma^2 / n
  base_variance <- sigma_sq / n
  
  # 3. Finite Population Correction (if N is finite)
  if (!is.infinite(N)) {
    stopifnot(is.numeric(N), length(N) == 1)
    stopifnot(N >= n) # Population must be larger than sample
    
    # Correction factor for Variance is ((N-n)/(N-1))
    # Note: For Standard Error (sigma), it is the sqrt() of this.
    correction_factor <- (N - n) / (N - 1)
    
    return(base_variance * correction_factor)
  }
  
  # 4. Return result for Infinite Population
  return(base_variance)
}
#' Chi-Squared Goodness of Fit Test
#'
#' Performs a Chi-Squared test to determine if a dataset follows a specific 
#' theoretical distribution.
#'
#' The test statistic is calculated as:
#' \eqn{\hat{\chi}^2 = \sum \frac{(o_i - e_i)^2}{e_i}}[cite: 935].
#'
#' @param observed A numeric vector of observed frequencies ($o_i$).
#' @param expected A numeric vector of expected frequencies ($e_i$).
#' @param lambda The significance level ($\alpha$) for the decision (default 0.05).
#'   If p-value < lambda, the Null Hypothesis is rejected[cite: 563].
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{statistic}: The calculated Chi-squared value.
#'   \item \code{df}: Degrees of freedom (calculated as k - 1)[cite: 936].
#'   \item \code{p_value}: The probability of observing a statistic this extreme.
#'   \item \code{decision}: Text string indicating whether to reject $H_0$.
#' }
#' @note 
#' If parameters were estimated from the sample to generate the expected counts, 
#' the degrees of freedom ideally should be \eqn{k - m - 1}, where m is the number 
#' of estimated parameters[cite: 951]. This function defaults to \eqn{k - 1}.
#'
#' @examples
#' # Example from Slide 74 (Independence Test context)
#' obs <- c(198, 28, 62, 39, 6, 12, 105, 15, 35)
#' exp <- c(196.9, 28.2, 62.7, 38.9, 5.5, 12.4, 106.0, 15.1, 33.7)
#' chi_squared_test(obs, exp, lambda = 0.05)
#' @export
test_chi_squared <- function(observed, expected, lambda = 0.05) {
  
  # 1. Validation
  if (length(observed) != length(expected)) {
    stop("Error: Datasets must have the same length.")
  }
  
  # 2. Normalize Expected Frequencies
  # [cite_start]Ensure Sum(O) == Sum(E) as required for the calculation[cite: 948].
  total_obs <- sum(observed)
  total_exp <- sum(expected)
  
  if (abs(total_obs - total_exp) > 1e-6) {
    expected <- expected * (total_obs / total_exp)
  }
  
  # 3. Calculate Chi-Squared Statistic
  # Formula: Sum( (O - E)^2 / E )
  chi_sq_stat <- sum((observed - expected)^2 / expected)
  
  # 4. Calculate Degrees of Freedom and P-Value
  df <- length(observed) - 1
  p_value <- pchisq(chi_sq_stat, df, lower.tail = FALSE)
  
  # 5. Make Decision based on Lambda (Significance Level alpha)
  reject_null <- p_value < lambda
  decision_text <- ifelse(reject_null, 
                          "Reject Null Hypothesis (Significant Difference)", 
                          "Fail to Reject Null Hypothesis (No Significant Difference)")
  
  # 6. Return Result
  return(list(
    statistic = chi_sq_stat,
    df = df,
    p_value = p_value,
    significance_level_lambda = lambda,
    reject_null = reject_null,
    decision = decision_text
  ))
}

#' Get Expected Counts for Normal Distribution
#'
#' Fits a Normal Distribution \eqn{N(\mu, \sigma)} to the raw data and calculates
#' expected frequencies for goodness-of-fit testing.
#'
#' Parameters are estimated from the sample:
#' \itemize{
#'   \item Mean (\eqn{\mu}) is estimated by \eqn{\bar{X}}[cite: 902].
#'   \item Standard Deviation (\eqn{\sigma}) is estimated by the sample SD ($S$)[cite: 902].
#' }
#'
#' @param raw_data A numeric vector of raw observations.
#' @return A numeric vector of expected counts based on histogram bins.
#' @export
get_expected_normaldistr_counts <- function(raw_data) {
  
  # 1. Clean data and estimate parameters
  data <- na.omit(raw_data)
  n <- length(data)
  mean_est <- mean(data)
  sd_est <- sd(data)
  
  # 2. Determine Bins (using the default 'Sturges' method R uses for histograms)
  h <- hist(data, plot = FALSE)
  breaks <- h$breaks
  n_bins <- length(breaks) - 1
  
  # 3. Calculate Probabilities for each bin
  probs <- numeric(n_bins)
  
  for(i in 1:n_bins) {
    upper <- breaks[i+1]
    lower <- breaks[i]
    
    if (i == 1) {
      # First bin: -Infinity to Upper
      probs[i] <- pnorm(upper, mean = mean_est, sd = sd_est)
    } else if (i == n_bins) {
      # Last bin: Lower to +Infinity
      probs[i] <- 1 - pnorm(lower, mean = mean_est, sd = sd_est)
    } else {
      # Middle bins
      probs[i] <- pnorm(upper, mean = mean_est, sd = sd_est) - 
        pnorm(lower, mean = mean_est, sd = sd_est)
    }
  }
  
  # 4. Return the Expected Frequencies Vector
  return(probs * n)
}

#' Get Expected Counts for Poisson Distribution
#'
#' Fits a Poisson Distribution \eqn{P(\lambda)} to the raw data and calculates
#' expected frequencies.
#'
#' The parameter \eqn{\lambda} is estimated using the sample mean \eqn{\bar{X}}[cite: 900].
#'
#' @param raw_data A numeric vector of raw observations (integers).
#' @return A numeric vector of expected counts. The last bin represents
#'   "observed max value or more" to ensure total probability sums to 1.
#' @examples
#' # Matches logic from Slide 69 (Geiger counter example) [cite: 961]
#' get_expected_poisson_counts(c(0, 1, 1, 2, 2, 2, 3))
#' @export
get_expected_poissondistr_counts <- function(raw_data) {
  
  # 1. Setup
  data <- na.omit(raw_data)
  N <- length(data)
  
  # 2. Estimate Parameter (Lambda)
  lambda_est <- mean(data)
  
  # 3. Determine Range (0 to Max)
  max_val <- max(data)
  
  # 4. Calculate Probabilities
  # For bins 0 to (max-1), use standard Poisson probability
  probs <- dpois(0:(max_val - 1), lambda = lambda_est)
  
  # For the last bin (max_val), use 1 - sum(others) to capture the "tail"
  # [cite_start]This corresponds to "max_val or more" (e.g. Slide 69 uses "6 or more")[cite: 962].
  prob_tail <- 1 - sum(probs)
  
  # Combine
  all_probs <- c(probs, prob_tail)
  
  # 5. Return ONLY the Expected Counts Vector
  return(all_probs * N)
}

#' Get Expected Counts for Binomial Distribution
#'
#' Fits a Binomial Distribution \eqn{B(n, p)} to the raw data and calculates
#' expected frequencies.
#'
#' The parameter \eqn{p} is estimated using the sample proportion[cite: 901].
#'
#' @param raw_data A numeric vector of raw observations.
#' @param size The number of trials (n). If NULL, estimates 'n' as the maximum observed value.
#' @return A numeric vector of expected counts.
#' @export
get_expected_binomialdistr_counts <- function(raw_data, size = NULL) {
  
  # 1. Setup
  data <- na.omit(raw_data)
  N <- length(data)
  
  # 2. Estimate Parameters
  # If size isn't given, assume max observed value is the number of trials
  size_est <- if(is.null(size)) max(data) else size
  prob_est <- mean(data) / size_est
  
  # 3. Calculate Probabilities
  # Binomial is finite (0 to size), so we just calculate all of them directly
  probs <- dbinom(0:size_est, size = size_est, prob = prob_est)
  
  # Normalize slightly just in case of floating point rounding errors
  probs <- probs / sum(probs)
  
  # 4. Return ONLY the Expected Counts Vector
  return(probs * N)
}
#' Confidence Interval for the Mean of a Normal Distribution
#'
#' Calculates the (1-alpha) confidence interval for the population mean mu.
#' Handles three cases based on inputs:
#' 1. Known Variance: Uses Normal distribution (Z).
#' 2. Unknown Variance, Large Sample (n > 30): Uses Normal approximation (Z).
#' 3. Unknown Variance, Small Sample (n <= 30): Uses Student's t-distribution.
#'
#' @param mean Sample mean (x_bar).
#' @param sd Standard deviation (sigma if known, s if unknown).
#' @param n Sample size.
#' @param variance_known Logical. TRUE if 'sd' is the population sigma. FALSE if 'sd' is sample s.
#' @param alpha Significance level (default 0.05 for 95% confidence).
#' @return A numeric vector of length 2: c(lower_bound, upper_bound).
#' @export
get_ci_mean_normal <- function(mean, sd, n, variance_known = FALSE, alpha = 0.05) {
  # Validation
  stopifnot(n > 0, sd > 0, alpha > 0, alpha < 1)
  
  # Determine Critical Value and Standard Error
  if (variance_known) {
    # Case 1: Known Variance -> Z-distribution
    crit_val <- qnorm(1 - alpha / 2)
  } else {
    if (n > 30) {
      # Case 2: Unknown Var, Large Sample -> Z-distribution approx
      crit_val <- qnorm(1 - alpha / 2)
    } else {
      # Case 3: Unknown Var, Small Sample -> T-distribution
      crit_val <- qt(1 - alpha / 2, df = n - 1)
    }
  }
  
  # Margin of Error
  margin_error <- crit_val * (sd / sqrt(n))
  
  # Interval
  return(c(lower = mean - margin_error, upper = mean + margin_error))
}
#' Confidence Interval for the Variance of a Normal Distribution
#'
#' Calculates the (1-alpha) confidence interval for the population variance sigma^2.
#' Uses the Chi-Square distribution.
#' Note: The interval is not symmetric around the sample variance.
#'
#' @param s2 Sample quasi-variance (s^2).
#' @param n Sample size.
#' @param alpha Significance level (default 0.05).
#' @return A numeric vector of length 2: c(lower_bound, upper_bound).
#' @export
get_ci_variance_normal <- function(s2, n, alpha = 0.05) {
  stopifnot(n > 1, s2 >= 0)
  
  df <- n - 1
  
  # Critical Values from Chi-Square
  # Lower limit uses the Upper Tail critical value (denominator is larger)
  # Upper limit uses the Lower Tail critical value (denominator is smaller)
  chi_upper <- qchisq(1 - alpha / 2, df = df)
  chi_lower <- qchisq(alpha / 2, df = df)
  
  lower <- (df * s2) / chi_upper
  upper <- (df * s2) / chi_lower
  
  return(c(lower = lower, upper = upper))
}
#' Confidence Interval for a Binomial Proportion
#'
#' Calculates the (1-alpha) confidence interval for the population proportion p.
#' Assumes large sample size (Normal approximation).
#' @param p_hat Sample proportion (between 0 and 1).
#' @param n Sample size.
#' @param alpha Significance level (default 0.05).
#' @return A numeric vector of length 2.
#' @export
get_ci_proportion_binomial <- function(p_hat, n, alpha = 0.05) {
  stopifnot(p_hat >= 0, p_hat <= 1, n > 0)
  
  # Critical Value (Z)
  z_crit <- qnorm(1 - alpha / 2)
  
  # Standard Error
  se <- sqrt((p_hat * (1 - p_hat)) / n)
  
  margin_error <- z_crit * se
  
  # Clamp results to [0, 1] as probabilities cannot exceed these
  lower <- max(0, p_hat - margin_error)
  upper <- min(1, p_hat + margin_error)
  
  return(c(lower = lower, upper = upper))
}
#' Confidence Interval for Difference of Two Means
#'
#' Calculates the CI for (mu1 - mu2). Handles known/unknown variances, 
#' large/small samples, and equal/unequal variance assumptions.
#' @param mean1,mean2 Sample means.
#' @param sd1,sd2 Standard deviations (sigma or s).
#' @param n1,n2 Sample sizes.
#' @param variance_known Logical. If TRUE, sd1/sd2 are population sigmas.
#' @param var_equal Logical. Only used if variance_known=FALSE and small samples. 
#'        Assumes sigma1^2 = sigma2^2 (Pooled Variance).
#' @param alpha Significance level.
#' @export
get_ci_diff_means <- function(mean1, mean2, sd1, sd2, n1, n2, 
                          variance_known = FALSE, var_equal = FALSE, alpha = 0.05) {
  
  diff_mean <- mean1 - mean2
  
  if (variance_known) {
    # [cite_start]Case 1: Known Variances (Z-dist) [cite: 286]
    crit_val <- qnorm(1 - alpha/2)
    se <- sqrt(sd1^2/n1 + sd2^2/n2)
    
  } else {
    # Unknown Variances
    if (n1 + n2 > 30) {
      # [cite_start]Case 2: Large Samples (Z-dist approx) [cite: 286]
      crit_val <- qnorm(1 - alpha/2)
      se <- sqrt(sd1^2/n1 + sd2^2/n2)
      
    } else {
      # Small Samples (T-dist)
      if (var_equal) {
        # [cite_start]Case 3: Equal Variances (Pooled) [cite: 286, 288]
        df <- n1 + n2 - 2
        crit_val <- qt(1 - alpha/2, df = df)
        
        # Pooled Variance Calculation
        sp_sq <- ((n1 - 1)*sd1^2 + (n2 - 1)*sd2^2) / df
        se <- sqrt(sp_sq) * sqrt(1/n1 + 1/n2)
        
      } else {
        # [cite_start]Case 4: Unequal Variances (Welch) [cite: 290, 291]
        # Welch-Satterthwaite Degrees of Freedom
        num <- (sd1^2/n1 + sd2^2/n2)^2
        den <- ((sd1^2/n1)^2 / (n1 + 1)) + ((sd2^2/n2)^2 / (n2 + 1)) 
        # Note: Slide 21 uses (n+1) in denom? Standard Welch uses (n-1). 
        # [cite_start]Checking Slide 21[cite: 291]: It says n1+1 and n2+1. 
        # CAUTION: Standard R 't.test' uses (n-1). I will follow the Slide exactly.
        f <- (num / den) - 2 
        
        crit_val <- qt(1 - alpha/2, df = f)
        se <- sqrt(sd1^2/n1 + sd2^2/n2)
      }
    }
  }
  
  margin <- crit_val * se
  return(c(lower = diff_mean - margin, upper = diff_mean + margin))
}
#' Confidence Interval for Ratio of Variances
#'
#' Calculates the (1-alpha) CI for the ratio sigma1^2 / sigma2^2.
#' Uses the F-
#' distribution.
#'
#' @param s1_sq,s2_sq Sample quasi-variances (s^2).
#' @param n1,n2 Sample sizes.
#' @param alpha Significance level.
#' @export
get_ci_ratio_variances <- function(s1_sq, s2_sq, n1, n2, alpha = 0.05) {
  stopifnot(s1_sq > 0, s2_sq > 0, n1 > 1, n2 > 1)
  
  ratio <- s1_sq / s2_sq
  df1 <- n1 - 1
  df2 <- n2 - 1
  
  # F Critical Values
  # Formula: Lower = Ratio / F(upper), Upper = Ratio / F(lower)
  f_upper <- qf(1 - alpha/2, df1 = df1, df2 = df2)
  f_lower <- qf(alpha/2, df1 = df1, df2 = df2)
  
  return(c(lower = ratio / f_upper, upper = ratio / f_lower))
}
#' Confidence Interval for Difference of Proportions
#'
#' Calculates the (1-alpha) CI for (p1 - p2).
#' Assumes large samples (Normal approximation).
#' @param p1,p2 Sample proportions.
#' @param n1,n2 Sample sizes.
#' @param alpha Significance level.
#' @export
get_ci_diff_proportions <- function(p1, p2, n1, n2, alpha = 0.05) {
  # Validation for large sample assumption mentioned in slide
  if (n1 + n2 <= 30) warning("Sample sum <= 30. Normal approximation may be inaccurate.")
  
  diff_p <- p1 - p2
  crit_val <- qnorm(1 - alpha/2)
  
  # [cite_start]Standard Error [cite: 303]
  se <- sqrt( (p1*(1-p1)/n1) + (p2*(1-p2)/n2) )
  
  margin <- crit_val * se
  
  return(c(lower = diff_p - margin, upper = diff_p + margin))
}
#' Calculate Chi-Squared Statistic
#'
#' Calculates the test statistic (E) for a contingency table.
#' Formula: Sum( (Observed - Expected)^2 / Expected )
#'
#' @param observed_matrix A numeric matrix of observed counts.
#' @return A numeric scalar representing the calculated Chi-squared statistic.
#' @export
get_matrix_chi_squared_statistic <- function(observed_matrix) {
  
  # 1. Validation
  if(!is.matrix(observed_matrix)) stop("Input must be a matrix.")
  
  # 2. Calculate Expected Values (e_ij)
  row_totals <- rowSums(observed_matrix)
  col_totals <- colSums(observed_matrix)
  n <- sum(observed_matrix)
  
  # Outer product calculates (row_i * col_j) for all combinations
  expected_matrix <- outer(row_totals, col_totals) / n
  
  # 3. Calculate the Statistic (Summation)
  # Formula: Sum over all i, j of [ (o_ij - e_ij)^2 / e_ij ]
  chi_squared_stat <- sum((observed_matrix - expected_matrix)^2 / expected_matrix)
  
  return(chi_squared_stat)
}
#' Calculate Expected Values Matrix
#'
#' Calculates the expected frequencies for a contingency table assuming
#' independence between rows and columns (H0).
#'
#' @param observed_matrix A numeric matrix or table containing observed counts.
#' @return A matrix of the same dimensions containing the expected values.
#' @export
get_matrix_expected_values <- function(observed_matrix) {
  
  # 1. Validation
  if(!is.matrix(observed_matrix) && !is.table(observed_matrix)) {
    stop("Input must be a matrix or a table.")
  }
  
  # 2. Calculate Marginals (Row sums, Col sums, and Total n)
  row_totals <- rowSums(observed_matrix)
  col_totals <- colSums(observed_matrix)
  n <- sum(observed_matrix)
  
  # 3. Calculate Expected Values
  # Formula: (RowTotal * ColTotal) / n
  # We use outer() to multiply every row sum by every col sum efficiently
  expected_matrix <- outer(row_totals, col_totals) / n
  
  # 4. Check the condition from the image (e_ij > 5)
  # The image states: "It must be verified that e_ij > 5"
  if (any(expected_matrix <= 5)) {
    warning("Condition failed: Some expected values are <= 5. ",
            "The Chi-square approximation may be inaccurate.")
  }
  
  return(expected_matrix)
}
#' Get Chi-Squared Critical Value
#'
#' Calculates the critical value (threshold) for the Chi-squared test
#' based on the matrix dimensions and significance level alpha.
#' Corresponds to: Chi2_alpha; (k-1)(m-1)
#'
#' @param observed_matrix A numeric matrix (to determine k and m).
#' @param alpha Significance level (default is 0.05).
#' @return The critical value from the Chi-squared distribution.
#' @export
get_matrix_chi_critical_value <- function(observed_matrix, alpha = 0.05) {
  
  # 1. Get dimensions (k = rows, m = columns)
  k <- nrow(observed_matrix)
  m <- ncol(observed_matrix)
  
  # 2. Calculate Degrees of Freedom
  # Formula: (k - 1) * (m - 1)
  df <- (k - 1) * (m - 1)
  
  # 3. Calculate Critical Value
  # We use lower.tail = FALSE to get the value for the upper tail area (alpha)
  critical_value <- qchisq(alpha, df = df, lower.tail = FALSE)
  
  return(critical_value)
}
