#' Compare two LOESS fits
#'
#' @param x Vector of x values
#' @param y1 Vector of y1 values
#' @param y2 Vector of y2 values
#' @param span Smoothing parameter for LOESS
#' @param significance_level Significance level for confidence intervals
#' @param grid Method to determine evaluation grid ("actual_x" or "uniform_x")
#' @param x_name Name for x variable in plots
#' @param y1_name Name for y1 variable in plots
#' @param y2_name Name for y2 variable in plots
#' @return List with comparison metrics and plotting data
#' @family loess comparison
#' @export
compare_loess_fits <- function(x, y1, y2, span = 0.75,
                               significance_level = 0.01,
                               grid = "actual_x", # New parameter - can be "actual_x" or "uniform_x"
                               x_name = deparse(substitute(x)),
                               y1_name = deparse(substitute(y1)),
                               y2_name = deparse(substitute(y2))) {
  # Check inputs
  if (length(y1) != length(x) || length(y2) != length(x)) {
    stop("x, y1, and y2 must have the same length")
  }

  # Validate grid parameter
  if (!grid %in% c("actual_x", "uniform_x")) {
    stop("grid parameter must be either 'actual_x' or 'uniform_x'")
  }

  # Clean variable names (remove any "data$" prefix)
  x_name <- gsub(".*\\$", "", x_name)
  y1_name <- gsub(".*\\$", "", y1_name)
  y2_name <- gsub(".*\\$", "", y2_name)

  # Create data frames
  data1 <- data.frame(x = x, y = y1)
  data2 <- data.frame(x = x, y = y2)

  # Determine evaluation points
  if (grid == "actual_x") {
    # Use all actual x values (sorted), including duplicates, to weight regions with more data points
    evaluation_points <- sort(x[!is.na(x)])
  } else if (grid == "uniform_x") {
    # Create a uniform grid with 100 points covering the range of x
    x_range <- range(x, na.rm = TRUE)
    evaluation_points <- seq(x_range[1], x_range[2], length.out = 100)
  } else {
    stop("Invalid value for grid. Must be 'actual_x' or 'uniform_x'.")
  }

  # Fit LOESS models
  loess1 <- loess(y ~ x, data = data1, span = span)
  loess2 <- loess(y ~ x, data = data2, span = span)

  # Predict with standard errors
  pred1 <- predict(loess1, newdata = data.frame(x = evaluation_points), se = TRUE)
  pred2 <- predict(loess2, newdata = data.frame(x = evaluation_points), se = TRUE)

  # Extract fits and standard errors
  fit1 <- pred1$fit
  fit2 <- pred2$fit
  se1 <- pred1$se.fit
  se2 <- pred2$se.fit

  # Calculate difference and propagate error
  # For difference of independent variables: var(A-B) = var(A) + var(B)
  diff <- fit1 - fit2
  diff_var <- se1^2 + se2^2
  diff_se <- sqrt(diff_var)

  # Calculate confidence intervals
  # Using normal approximation: mean +/- z*se
  z_score <- qnorm(1 - significance_level/2)
  lower <- diff - z_score * diff_se
  upper <- diff + z_score * diff_se

  # Determine significance (zero outside CI)
  is_significant <- (lower > 0) | (upper < 0)

  # Calculate percentage of evaluation points with significant difference
  percent_significant <- mean(is_significant, na.rm = TRUE) * 100

  # Calculate coverage (percentage where fits are statistically similar)
  coverage <- 100 - percent_significant

  # For weighted RMSE, use inverse variance as weights
  inv_var_weights <- 1 / diff_var

  # Handle NAs or zeros in weights
  inv_var_weights[is.na(inv_var_weights) | !is.finite(inv_var_weights)] <- 0

  # Normalize weights to sum to 1 (if any non-zero weights exist)
  if (sum(inv_var_weights) > 0) {
    inv_var_weights <- inv_var_weights / sum(inv_var_weights)
  }

  # Calculate weighted RMSE
  weighted_mse <- sum(inv_var_weights * diff^2, na.rm = TRUE)
  weighted_rmse <- sqrt(weighted_mse)

  # Create data frames for plotting (with original variable names)
  # Use the variable names provided or extracted
  fit_data <- data.frame(
    x = evaluation_points,
    y1_fit = fit1,
    y2_fit = fit2,
    y1_se = se1,
    y2_se = se2
  )

  diff_data <- data.frame(
    x = evaluation_points,
    diff = diff,
    lower = lower,
    upper = upper,
    is_significant = is_significant
  )

  # Create original data with original variable names
  original_data <- data.frame(
    x = rep(x, 2),
    y = c(y1, y2),
    series = rep(c(y1_name, y2_name), each = length(x))
  )

  # Return a simplified list with only the essential results
  return(list(
    # Key metrics
    weighted_rmse = weighted_rmse,
    percent_significant = percent_significant,
    coverage = coverage,

    # Data for plotting
    fit_data = fit_data,
    diff_data = diff_data,
    original_data = original_data,

    # Additional information
    span = span,
    significance_level = significance_level,
    grid = grid,

    # Variable names for plotting
    x_name = x_name,
    y1_name = y1_name,
    y2_name = y2_name
  ))
}

#' Plot LOESS comparison
#'
#' @param comparison_result Result from compare_loess_fits
#' @return Either a grid of plots or a list of plots
#' @noRd
plot_loess_comparison <- function(comparison_result) {
  # Check if ggplot2 is available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package ggplot2 is required but not installed")
  }

  # Unpack data
  fit_data <- comparison_result$fit_data
  diff_data <- comparison_result$diff_data
  original_data <- comparison_result$original_data

  # Get variable names
  x_name <- comparison_result$x_name
  y1_name <- comparison_result$y1_name
  y2_name <- comparison_result$y2_name

  # Create annotation text for the measures
  anno_text <- sprintf(
    "Weighted RMSE: %.4f\nCoverage: %.1f%%",
    comparison_result$weighted_rmse,
    comparison_result$coverage
  )

  # Create a common theme for both plots to ensure consistent sizing
  common_theme <- ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "top",
      plot.title = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 9)
    )

  # Plot 1: Original fits with standard error bands
  p1 <- ggplot2::ggplot() +
    # Original data points
    ggplot2::geom_point(data = original_data, ggplot2::aes(x = x, y = y, color = series), alpha = 0.3) +
    # Standard error bands for y1
    ggplot2::geom_ribbon(data = fit_data,
                ggplot2::aes(x = x, ymin = y1_fit - y1_se, ymax = y1_fit + y1_se),
                fill = "blue", alpha = 0.2) +
    # Standard error bands for y2
    ggplot2::geom_ribbon(data = fit_data,
                ggplot2::aes(x = x, ymin = y2_fit - y2_se, ymax = y2_fit + y2_se),
                fill = "red", alpha = 0.2) +
    # LOESS fit line for y1
    ggplot2::geom_line(data = fit_data, ggplot2::aes(x = x, y = y1_fit), color = "blue", size = 1) +
    # LOESS fit line for y2
    ggplot2::geom_line(data = fit_data, ggplot2::aes(x = x, y = y2_fit), color = "red", size = 1) +
    # Colors for the legend
    ggplot2::scale_color_manual(values = c("blue", "red"),
                       labels = c(y1_name, y2_name)) +
    common_theme +
    # Labels using the original variable names
    ggplot2::labs(title = "LOESS Fits Comparison",
         x = x_name, y = "Value",
         color = "") +  # Empty legend title for cleaner look
    # Add the metrics annotation
    ggplot2::annotate("text", x = min(fit_data$x),
             y = max(c(fit_data$y1_fit, fit_data$y2_fit), na.rm = TRUE),
             label = anno_text, hjust = 0, vjust = 1)

  # For prettier difference plot title
  diff_title <- paste("Difference:", y1_name, "-", y2_name, "with Confidence Intervals")

  # Truncate if too long
  if(nchar(diff_title) > 60) {
    diff_title <- paste0(substr(diff_title, 1, 57), "...")
  }

  # Plot 2: Difference with confidence polygon
  p2 <- ggplot2::ggplot(diff_data) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::geom_ribbon(ggplot2::aes(x = x, ymin = lower, ymax = upper),
                fill = "gray80", alpha = 0.5) +
    ggplot2::geom_line(ggplot2::aes(x = x, y = diff), color = "black", size = 1) +
    ggplot2::geom_point(data = subset(diff_data, is_significant),
               ggplot2::aes(x = x, y = diff), color = "red", size = 1.5) +
    common_theme +
    ggplot2::labs(title = diff_title,
         x = x_name,
         y = "Difference")  # Simplified y-axis label

  # Return plots with consistent widths
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    # Use layout to ensure equal widths
    return(gridExtra::grid.arrange(p1, p2,
                                   ncol = 1,
                                   heights = c(1, 1)))
  } else {
    return(list(fits = p1, diff = p2))
  }
}

#' Calculate weighted RMSE between two LOESS fits
#'
#' @param x Vector of predictor values (same for both models)
#' @param y1 Vector of first response values
#' @param y2 Vector of second response values
#' @param span Smoothing parameter for LOESS (default: 0.75)
#' @return Numeric weighted RMSE between the LOESS fits
#' @family loess comparison
#' @export
compare_loess_rmse <- function(x, y1, y2, span = 0.75) {
  # Check inputs
  if (length(x) != length(y1) || length(x) != length(y2)) {
    stop("x, y1, and y2 must have the same length")
  }

  # Create data frames
  data1 <- data.frame(x = x, y = y1)
  data2 <- data.frame(x = x, y = y2)

  # Use sorted x values as evaluation points
  evaluation_points <- sort(x[!is.na(x)])

  # Fit LOESS models
  loess1 <- loess(y ~ x, data = data1, span = span)
  loess2 <- loess(y ~ x, data = data2, span = span)

  # Predict with standard errors
  pred1 <- predict(loess1, newdata = data.frame(x = evaluation_points), se = TRUE)
  pred2 <- predict(loess2, newdata = data.frame(x = evaluation_points), se = TRUE)

  # Extract fits and standard errors
  fit1 <- pred1$fit
  fit2 <- pred2$fit
  se1 <- pred1$se.fit
  se2 <- pred2$se.fit

  # Calculate difference and propagate error
  diff <- fit1 - fit2
  diff_var <- se1^2 + se2^2

  # Calculate inverse variance weights
  inv_var_weights <- 1 / diff_var
  inv_var_weights[is.na(inv_var_weights) | !is.finite(inv_var_weights)] <- 0

  # Normalize weights
  if (sum(inv_var_weights) > 0) {
    inv_var_weights <- inv_var_weights / sum(inv_var_weights)
  }

  # Calculate weighted RMSE
  weighted_mse <- sum(inv_var_weights * diff^2, na.rm = TRUE)
  weighted_rmse <- sqrt(weighted_mse)

  return(weighted_rmse)
}
