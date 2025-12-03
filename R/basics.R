#' Mean (average) with NAs ignored
#' @param x Numeric vector.
#' @return Mean value of `x`.
#' @examples
#' mean_na(1:5)
#' @export
mean_na <- function(x) {
    stopifnot(is.numeric(x))
    mean(x, na.rm = TRUE)
}

#' Variance with NAs ignored
#' @param x Numeric vector.
#' @return Sample variance of `x`.
#' @examples
#' variance_na(1:5)
#' @export
variance_na <- function(x) {
    stopifnot(is.numeric(x))
    stats::var(x, na.rm = TRUE)
}

#' Standard deviation with NAs ignored
#' @param x Numeric vector.
#' @return Sample standard deviation of `x`.
#' @examples
#' sd_na(1:5)
#' @export
sd_na <- function(x) {
    stopifnot(is.numeric(x))
    stats::sd(x, na.rm = TRUE)
}

#' Quantiles
#' @param x Numeric vector.
#' @param probs Probabilities (0-1). Defaults to 0.25, 0.5, and 0.75.
#' @return Vector of quantiles.
#' @examples
#' quantiles(1:9)
#' @export
quantiles <- function(x, probs = c(0.25, 0.5, 0.75)) {
    stopifnot(is.numeric(x))
    stats::quantile(x, probs = probs, na.rm = TRUE, names = TRUE)
}

#' Normalization [0,1]
#' @param x Numeric vector.
#' @return Vector rescaled to [0,1]. If min==max, returns zeros.
#' @examples
#' normalize(1:10)
#' @export
normalize <- function(x) {
    stopifnot(is.numeric(x))
    rng <- range(x, na.rm = TRUE)
    if (isTRUE(all.equal(rng[1], rng[2]))) {
        return(ifelse(is.na(x), NA_real_, 0))
    }
    (x - rng[1]) / (rng[2] - rng[1])
}

#' Standardization (z-score)
#' @param x Numeric vector.
#' @return (x - mean)/sd with NAs ignored.
#' @examples
#' standardize(1:10)
#' @export
standardize <- function(x) {
    stopifnot(is.numeric(x))
    m <- mean(x, na.rm = TRUE)
    s <- stats::sd(x, na.rm = TRUE)
    if (isTRUE(all.equal(s, 0))) {
        return(ifelse(is.na(x), NA_real_, 0))
    }
    (x - m) / s
}

#' Quick correlation
#' @param x,y Numeric vectors.
#' @param method "pearson", "spearman", or "kendall".
#' @return Correlation between `x` and `y`.
#' @examples
#' correlation(1:10, 1:10)
#' @export
correlation <- function(x, y, method = c("pearson", "spearman", "kendall")) {
    method <- match.arg(method)
    stopifnot(is.numeric(x), is.numeric(y), length(x) == length(y))
    stats::cor(x, y, use = "complete.obs", method = method)
}

#' Simple linear regression
#' @param x Numeric predictor.
#' @param y Numeric response.
#' @return An `lm` object.
#' @examples
#' set.seed(1)
#' x <- 1:10
#' y <- x + rnorm(10)
#' m <- linear_regression(x, y)
#' summary(m)
#' @export
linear_regression <- function(x, y) {
    stopifnot(is.numeric(x), is.numeric(y), length(x) == length(y))
    stats::lm(y ~ x)
}

#' Quick plot (scatter/line or histogram)
#' @param x Numeric vector.
#' @param y Optional: if provided, makes plot(x,y). Otherwise, hist(x).
#' @param type Plot type for base R `plot` ("p" points, "l" lines, etc.).
#' @examples
#' plot_quick(1:10, (1:10)^2)
#' plot_quick(rnorm(100))
#' @export
plot_quick <- function(x, y = NULL, type = "p") {
    if (is.null(y)) {
        graphics::hist(x, main = "Histogram", xlab = deparse(substitute(x)))
    } else {
        graphics::plot(x, y, type = type, xlab = deparse(substitute(x)), ylab = deparse(substitute(y)))
    }
}

#' Simple residual diagnostics for a linear model
#' @param model An `lm` object.
#' @examples
#' m <- lm(mpg ~ wt, data = mtcars)
#' plot_residuals(m)
#' @export
plot_residuals <- function(model) {
    stopifnot(inherits(model, "lm"))
    op <- graphics::par(mfrow = c(1, 2))
    on.exit(graphics::par(op), add = TRUE)
    graphics::plot(stats::fitted(model), stats::residuals(model),
        xlab = "Fitted", ylab = "Residuals", main = "Residuals vs Fitted"
    )
    graphics::abline(h = 0, lty = 2)
    stats::qqnorm(stats::rstandard(model), main = "QQ of standardized residuals")
    stats::qqline(stats::rstandard(model))
}

#' Example data generator (linear with noise)
#' @param n Number of observations.
#' @param seed Optional seed for reproducibility.
#' @return data.frame with columns x, y.
#' @examples
#' df <- example_data(20, seed = 123)
#' @export
example_data <- function(n = 50, seed = NULL) {
    stopifnot(is.numeric(n), length(n) == 1, n > 1)
    if (!is.null(seed)) set.seed(seed)
    x <- seq_len(n)
    y <- 3 + 1.5 * x + stats::rnorm(n, sd = sd(x) / 2)
    data.frame(x = x, y = y)
}
