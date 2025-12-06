#' Discrete Binomial Distribution Probability
#'
#' Calculates the probability of getting exactly `x` successes in `n` trials
#' using the Binomial Probability Mass Function.
#'
#' @param n Number of trials (numeric scalar, must be >= x).
#' @param p Probability of success on each trial (between 0 and 1).
#' @param x Number of successes (numeric scalar).
#' @return Probability value.
#' @examples
#' # Probability of getting exactly 5 heads in 10 flips of a fair coin
#' discrete_binomial_distribution(n = 10, p = 0.5, x = 5)
#' @export
discrete_binomial_distribution <- function(n, p, x) {
  # 1. Validation of types and dimensions
  stopifnot(is.numeric(n), length(n) == 1)
  stopifnot(is.numeric(p), length(p) == 1)
  stopifnot(is.numeric(x), length(x) == 1)
  
  # 2. Validation of logical constraints
  stopifnot(p >= 0, p <= 1)      # Probability must be [0,1]
  stopifnot(n >= 0, x >= 0)      # Counts must be non-negative
  stopifnot(n >= x)              # Trials must be >= Successes
  
  # 3. Calculation using robust stats package
  stats::dbinom(x, size = n, prob = p)
}

#' Discrete Geometric Distribution Probability
#'
#' Calculates the probability that the first success occurs specifically on the
#' `k`-th trial. (Corresponds to the distribution definition supported on {1, 2, ...}).
#'
#' @param p Probability of success on each trial (0 < p <= 1).
#' @param k Number of trials until the first success (numeric, must be >= 1).
#' @return Probability value.
#' @examples
#' # Probability that the first success happens exactly on the 3rd trial (p=0.5)
#' discrete_geometric_distribution(p = 0.5, k = 3)
#' @export
discrete_geometric_distribution <- function(p, k) {
  # 1. Validation
  stopifnot(is.numeric(p), length(p) == 1)
  stopifnot(is.numeric(k), length(k) == 1)
  
  # p must be strictly > 0 for geometric distribution (must eventually succeed)
  stopifnot(p > 0, p <= 1)
  stopifnot(k >= 1) 
  
  # 2. Calculation
  # Note: R's stats::dgeom defines 'x' as number of failures *before* success.
  # If 'k' is the trial number of the first success, then there are k-1 failures.
  stats::dgeom(x = k - 1, prob = p)
}

#' Discrete Poisson Distribution Probability
#'
#' Calculates the probability of a given number of events occurring in a fixed
#' interval of time or space.
#'
#' @param lambda Rate parameter (average number of events). Must be non-negative.
#' @param x Number of events (non-negative integer).
#' @return Probability of observing exactly `x` events.
#' @examples
#' # Probability of 3 emails arriving when the average is 5
#' discrete_poisson_distribution(lambda = 5, x = 3)
#' @export
discrete_poisson_distribution <- function(lambda, x) {
  # 1. Validation
  stopifnot(is.numeric(lambda), length(lambda) == 1, lambda >= 0)
  stopifnot(is.numeric(x), length(x) == 1, x >= 0)
  
  # Check for non-integer x (warning or stop depending on strictness, 
  # here we enforce integer-like logic via the stats function, but ensure numeric input)
  
  # 2. Calculation
  # We use stats::dpois instead of manual formula (lambda^x * e^-lambda / x!)
  # to avoid floating point overflow on factorials and loss of precision on 'e'.
  stats::dpois(x, lambda)
}

#' Standard Normal Distribution Density
#'
#' Calculates the probability density function (PDF) for the standard normal
#' distribution (mean = 0, standard deviation = 1).
#'
#' Note: Unlike discrete distributions, for continuous distributions like this,
#' the return value is the "density" at `x`, not the probability of `x` (which is 0).
#'
#' @param x Numeric vector of quantiles (input values).
#' @return Density value(s) at `x`.
#' @examples
#' # Density at the mean (highest point)
#' standard_normal_distribution(0)
#'
#' # Density at +/- 1.96 (approx 95% interval boundaries)
#' standard_normal_distribution(c(-1.96, 1.96))
#' @export
#' Standard Normal Distribution Density (with Auto-Standardization)
#'
#' Calculates the probability density function (PDF) for the standard normal
#' distribution.
#'
#' If a custom `mean` and `sd` are provided, the function automatically
#' standardizes the input `x` to a Z-score \eqn{z = (x - \mu) / \sigma} before
#' calculating the density on the standard curve.
#'
#' @param x Numeric vector of quantiles (raw values).
#' @param mean Mean of the source distribution (defaults to 0).
#' @param sd Standard deviation of the source distribution (defaults to 1, must be > 0).
#' @return Density value(s) of the standard normal distribution at the standardized Z-score of `x`.
#' @examples
#' # Standard case (mean=0, sd=1) -> Density at 0
#' standard_normal_distribution(0)
#'
#' # Custom case: x=115, mean=100, sd=15 -> z=1 -> Density at 1
#' standard_normal_distribution(115, mean = 100, sd = 15)
#' @export
#' Standard Normal Distribution (PDF or CDF)
#'
#' Calculates the probability density (PDF) or cumulative probability (CDF)
#' for the standard normal distribution.
#'
#' If a custom `mean` and `sd` are provided, the function automatically
#' standardizes the input `x` to a Z-score \eqn{z = (x - \mu) / \sigma} before
#' calculating.
#'
#' @param x Numeric vector of quantiles.
#' @param mean Mean of the source distribution (defaults to 0).
#' @param sd Standard deviation of the source distribution (defaults to 1, must be > 0).
#' @param type Type of calculation: "density" for P(X=x) (PDF height),
#'   "cumulative" for P(X<=x) (CDF area).
#' @return Value of density or cumulative probability.
#' @examples
#' # Density at 0 (Height of the curve)
#' standard_normal_distribution(0, type = "density")
#'
#' # Probability that X <= 1.96 (Area to the left, approx 0.975)
#' standard_normal_distribution(1.96, type = "cumulative")
#' @export
standard_normal_distribution <- function(x, mean = 0, sd = 1, type = c("density", "cumulative")) {
  # 1. Validation
  type <- match.arg(type)
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(mean), length(mean) == 1)
  stopifnot(is.numeric(sd), length(sd) == 1, sd > 0)
  
  # 2. Standardization
  # Convert raw x to z-score (z = (x - mean) / sd)
  z <- (x - mean) / sd
  
  # 3. Calculation
  if (type == "density") {
    # P(X = x) context (Density height)
    stats::dnorm(z, mean = 0, sd = 1)
  } else {
    # P(X <= x) context (Cumulative area)
    stats::pnorm(z, mean = 0, sd = 1)
  }
}
#' Student's t-Distribution (PDF or CDF)
#'
#' Calculates the probability density (PDF) or cumulative probability (CDF)
#' for the Student's t-distribution with `df` degrees of freedom.
#'
#' This distribution is used instead of the Normal distribution when sample
#' sizes are small or population variance is unknown. As `df` increases,
#' this distribution converges to the Standard Normal distribution.
#'
#' @param x Numeric vector of quantiles.
#' @param df Degrees of freedom (must be > 0). Generally calculated as n - 1.
#' @param type Type of calculation: "density" for P(X=x) (height),
#'   "cumulative" for P(X<=x) (area to the left).
#' @return Value of density or cumulative probability.
#' @examples
#' # Density at 0 with 10 degrees of freedom
#' student_t_distribution(0, df = 10, type = "density")
#'
#' # Probability that t <= 2.228 with 10 df (approx 0.975)
#' student_t_distribution(2.228, df = 10, type = "cumulative")
#' @export
student_t_distribution <- function(x, df, type = c("density", "cumulative")) {
  # 1. Validation
  type <- match.arg(type)
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(df), length(df) == 1, df > 0)
  
  # 2. Calculation
  if (type == "density") {
    # P(X = x) context
    stats::dt(x, df = df)
  } else {
    # P(X <= x) context
    stats::pt(x, df = df)
  }
}