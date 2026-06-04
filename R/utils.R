#' Calculate weighted RMSE between loess-smoothed predictions of y1 and y2 vs x
#'
#' @param y1 First vector of values
#' @param y2 Second vector of values
#' @param x Vector of x coordinates
#' @return Weighted RMSE value
#' @family error metrics
#' @export
rmse_trend <- function(y1, y2, x) {
  # Input validation
  if (length(y1) != length(y2) || length(y1) != length(x)) {
    stop("Vectors y1, y2, and x must have the same length")
  }
  if (any(is.na(y1)) || any(is.na(y2)) || any(is.na(x))) {
    stop("Input vectors must not contain NA values")
  }

  # Fit loess models for y1 vs x and y2 vs x
  fit1 <- loess(y1 ~ x, span = 0.75, se = TRUE)
  fit2 <- loess(y2 ~ x, span = 0.75, se = TRUE)

  # Get predicted values and standard errors at x
  pred1 <- predict(fit1, newdata = data.frame(x = x), se = TRUE)
  pred2 <- predict(fit2, newdata = data.frame(x = x), se = TRUE)

  # Extract fitted values
  y1_smooth <- pred1$fit
  y2_smooth <- pred2$fit

  # Extract standard errors (uncertainty)
  se1 <- pred1$se.fit
  se2 <- pred2$se.fit

  # Calculate weights as inverse of combined variance (1 / (se1^2 + se2^2))
  weights <- 1 / (se1^2 + se2^2)

  # Handle potential infinite or NA weights
  weights[is.infinite(weights) | is.na(weights)] <- 1

  # Calculate weighted RMSE
  squared_diff <- (y1_smooth - y2_smooth)^2
  weighted_mse <- sum(weights * squared_diff) / sum(weights)
  rmse <- sqrt(weighted_mse)

  return(rmse)
}

#' Calculate RMSE between two vectors
#'
#' @param y1 First vector of values
#' @param y2 Second vector of values
#' @return RMSE value
#' @family error metrics
#' @export
rmse <- function(y1, y2) {
  sqrt(mean((y1 - y2)^2, na.rm = TRUE))
}

#' Calculate R-squared between two vectors
#'
#' @param y1 First vector of values
#' @param y2 Second vector of values
#' @return R-squared value
#' @family error metrics
#' @export
r2 <- function(y1, y2) {
  cor(y1, y2, use = "pairwise.complete.obs")^2
}

#' Calculate Mean Relative Residual between two vectors
#'
#' Calculates the mean of absolute differences relative to the reference values.
#' MRR = (1/n) * sum(|y1 - y2| / |y1|)
#'
#' @param y1 Reference vector of values
#' @param y2 Comparison vector of values
#' @return Mean Relative Residual value
#' @family error metrics
#' @export
mrr <- function(y1, y2) {
  # Input validation
  if (length(y1) != length(y2)) {
    stop("Vectors y1 and y2 must have the same length")
  }
  if (any(is.na(y1)) || any(is.na(y2))) {
    stop("Input vectors must not contain NA values")
  }

  # Avoid division by zero
  valid_indices <- which(y1 != 0)
  if (length(valid_indices) < length(y1)) {
    warning("Some reference values are zero. Using only non-zero values for calculation.")
    y1 <- y1[valid_indices]
    y2 <- y2[valid_indices]
  }

  # Calculate relative residuals
  relative_residuals <- abs(y1 - y2) / abs(y1)

  # Calculate mean
  mean(relative_residuals, na.rm = TRUE)
}

#' Calculate Mean Relative Residual between loess-smoothed trends
#'
#' Fits LOESS curves to both y1 vs x and y2 vs x, then calculates
#' the mean relative residual between the smoothed curves.
#'
#' @param y1 First vector of values to compare
#' @param y2 Second vector of values to compare
#' @param x Vector of x coordinates for the trends
#' @return Mean Relative Residual between trends
#' @family error metrics
#' @export
mrr_trend <- function(y1, y2, x) {
  # Input validation
  if (length(y1) != length(y2) || length(y1) != length(x)) {
    stop("Vectors y1, y2, and x must have the same length")
  }
  if (any(is.na(y1)) || any(is.na(y2)) || any(is.na(x))) {
    stop("Input vectors must not contain NA values")
  }

  # Fit loess models for y1 vs x and y2 vs x
  fit1 <- loess(y1 ~ x, span = 0.75, se = TRUE)
  fit2 <- loess(y2 ~ x, span = 0.75, se = TRUE)

  # Get predicted values at x
  y1_smooth <- predict(fit1, newdata = data.frame(x = x))
  y2_smooth <- predict(fit2, newdata = data.frame(x = x))

  # Calculate MRR between smoothed curves
  # Avoid division by zero
  valid_indices <- which(y1_smooth != 0)
  if (length(valid_indices) < length(y1_smooth)) {
    warning("Some smoothed reference values are zero. Using only non-zero values for calculation.")
    y1_smooth <- y1_smooth[valid_indices]
    y2_smooth <- y2_smooth[valid_indices]

    # We need to subset the same elements from y2_smooth
    y2_smooth <- y2_smooth[valid_indices]
  }

  # Calculate relative residuals
  relative_residuals <- abs(y1_smooth - y2_smooth) / abs(y1_smooth)

  # Calculate mean
  mean(relative_residuals, na.rm = TRUE)
}
