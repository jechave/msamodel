# MSA model fitting — maximum-likelihood point estimation
# Maximises the profiled Gaussian log-likelihood (calculate_loglik_lrmsd_i_msa) over
# (a1, log2(a2+1)) and returns a point estimate plus an asymptotic covariance.

#' Maximum-likelihood point fit of the lrmsd_i MSA model
#'
#' Maximises the profiled Gaussian log-likelihood (`calculate_loglik_lrmsd_i_msa()`)
#' over `(a1, a2)` by numerical optimisation, returning a point estimate plus an
#' asymptotic covariance from the Hessian at the optimum. Intended for large proteins
#' and path simulations.
#'
#' The optimiser works in the coordinates `a1` and
#' `log2(a2 + 1)` (so `a2 = 2^(log2(a2+1)) - 1 >= 0`), on the box `a1_range` ×
#' `log2_a2_plus1_range` — the coordinates in which the prior is uniform. The
#' returned covariance `cov` is on the `(a1, log2(a2+1))` scale; the standard error
#' of `a2` is obtained by the delta method (`da2/d(log2(a2+1)) = 2^(log2(a2+1)) * ln 2`).
#'
#' @param spm_pp Preprocessed data from [preprocess_spm()] (must include `site_map`).
#' @param observed_data Tibble with columns `pdb_site` and `lrmsd_i_obs` (the fit
#'   target), as documented for `calculate_loglik_lrmsd_i_msa()`.
#' @param a1_range Length-2 `[min, max]` box bound for `a1`.
#' @param log2_a2_plus1_range Length-2 `[min, max]` box bound for `log2(a2 + 1)`.
#' @param init Optional length-2 numeric start `c(a1, log2(a2+1))`. When `NULL`
#'   (default), a deterministic coarse grid-max of the likelihood over the box is
#'   used as the start (robust against a bad local start; cheap).
#' @param grid_n Number of points per axis for the default grid-max start (ignored
#'   when `init` is supplied).
#' @return A list with the point estimate and asymptotic uncertainty:
#'   \describe{
#'     \item{a1, a2}{Point estimate of the stability (`a1`) and activity (`a2`)
#'       selection strengths, on the natural scale (the paper's `aS`, `aA`).}
#'     \item{logLik}{Profiled Gaussian log-likelihood at the optimum.}
#'     \item{deviance}{Residual deviance `sum(residuals^2)` (= `nobs * sigma_hat^2`).
#'       A goodness-of-fit primitive; consumed by [gof_lrmsd_i_msa_ml()].}
#'     \item{null_deviance}{Deviance of the flat/mean-only null (a constant profile),
#'       `sum((observed - mean(observed))^2)`. With `deviance`, gives
#'       `D2 = 1 - deviance/null_deviance`.}
#'     \item{nobs}{Number of matched (observed) sites the likelihood scores.}
#'     \item{k}{Free-parameter count for AIC/BIC: `a1`, `a2`, and the profiled `sigma`
#'       (`3`), matching the `logLik.lm`/`broom` convention.}
#'     \item{sigma_hat}{Profiled noise scale `sqrt(mean(residuals^2))` at the optimum.}
#'     \item{cov}{2×2 covariance matrix on the `(a1, log2(a2+1))` scale.}
#'     \item{se_a1, se_a2}{Standard errors; `se_a2` via the delta method.}
#'     \item{convergence}{`optim` convergence code (0 = success).}
#'   }
#' @seealso [gof_lrmsd_i_msa_ml()] (goodness-of-fit from this fit),
#'   [predict_lrmsd_i_msa_ml()] (propagate the fit to a banded profile),
#'   `calculate_loglik_lrmsd_i_msa()` (the objective).
#' @family fitting
#' @export
#' @examples
#' \dontrun{
#' pp <- preprocess_spm(znb_spm)
#' ml <- fit_lrmsd_i_msa_ml(pp, znb_profile)
#' c(a1 = ml$a1, a2 = ml$a2)
#' }
fit_lrmsd_i_msa_ml <- function(spm_pp,
                       observed_data,
                       a1_range = c(0, 10),
                       log2_a2_plus1_range = c(0, 13),
                       init = NULL,
                       grid_n = 25) {
  # Validate box bounds -- fail loud.
  if (length(a1_range) != 2 || a1_range[1] >= a1_range[2]) {
    stop("a1_range must be a vector of length 2 with min < max")
  }
  if (length(log2_a2_plus1_range) != 2 ||
      log2_a2_plus1_range[1] >= log2_a2_plus1_range[2]) {
    stop("log2_a2_plus1_range must be a vector of length 2 with min < max")
  }

  lower <- c(a1_range[1], log2_a2_plus1_range[1])
  upper <- c(a1_range[2], log2_a2_plus1_range[2])

  # Negative profiled log-likelihood in (a1, b) coordinates, b = log2(a2 + 1).
  # calculate_loglik_lrmsd_i_msa already mean-centers both profiles, so this IS the
  # shared objective (no separate centering here).
  nll <- function(theta) {
    -calculate_loglik_lrmsd_i_msa(spm_pp, observed_data,
                          a1 = theta[1], a2 = 2^theta[2] - 1)
  }

  # Start: caller-supplied, else a deterministic coarse grid-max over the box.
  if (is.null(init)) {
    a1_grid <- seq(lower[1], upper[1], length.out = grid_n)
    b_grid  <- seq(lower[2], upper[2], length.out = grid_n)
    grid    <- expand.grid(a1 = a1_grid, b = b_grid)
    ll      <- apply(grid, 1L, function(row) -nll(c(row[["a1"]], row[["b"]])))
    init    <- as.numeric(grid[which.max(ll), c("a1", "b")])
  } else {
    if (length(init) != 2) stop("init must be a length-2 numeric c(a1, log2(a2+1))")
    if (init[1] < lower[1] || init[1] > upper[1] ||
        init[2] < lower[2] || init[2] > upper[2]) {
      stop("init must lie within the box [a1_range] x [log2_a2_plus1_range]")
    }
  }

  opt <- optim(init, nll, method = "L-BFGS-B", lower = lower, upper = upper)

  par_fit <- opt$par
  a1_hat  <- par_fit[1]
  b_hat   <- par_fit[2]
  a2_hat  <- 2^b_hat - 1

  # Asymptotic covariance from the Hessian of the NLL at the optimum.
  H <- optimHess(par_fit, nll)
  cov <- tryCatch(solve(H), error = function(e) NULL)
  if (is.null(cov) || any(!is.finite(cov)) || any(diag(cov) <= 0)) {
    stop("Hessian at the ML optimum is singular or not positive-definite; ",
         "cannot form the asymptotic covariance. The likelihood may be flat in ",
         "one direction (e.g. a parameter pinned at a box bound).")
  }
  se <- sqrt(diag(cov))
  se_a1 <- se[1]
  se_b  <- se[2]
  # Delta method: a2 = 2^b - 1, da2/db = 2^b * ln 2.
  se_a2 <- abs(2^b_hat * log(2)) * se_b

  # Profiled noise scale at the optimum (same formula the likelihood uses).
  pred <- calculate_lrmsd_i_msa(spm_pp, a1_hat, a2_hat)
  obs <- observed_data %>%
    dplyr::select(pdb_site, lrmsd_i_obs) %>%
    inner_join(spm_pp$site_map, by = "pdb_site") %>%
    dplyr::select(i, lrmsd_i_obs)
  cmp <- obs %>%
    inner_join(pred, by = "i") %>%
    mutate(
      nlrmsd_i_obs = lrmsd_i_obs - mean(lrmsd_i_obs),
      nlrmsd_i_msa = lrmsd_i_msa - mean(lrmsd_i_msa)
    )
  residuals <- cmp$nlrmsd_i_obs - cmp$nlrmsd_i_msa
  sigma_hat <- sqrt(mean(residuals^2))

  # Goodness-of-fit primitives (flat/mean-only null; k counts sigma). D^2/AIC/BIC are
  # derived from these by gof_lrmsd_i_msa_ml(). deviance = sum(resid^2) (= glm's
  # residual deviance for a Gaussian, = nobs * sigma_hat^2); null_deviance vs a
  # constant profile.
  deviance      <- sum(residuals^2)
  null_deviance <- calculate_null_deviance(cmp$lrmsd_i_obs)

  list(
    a1            = unname(a1_hat),
    a2            = unname(a2_hat),
    logLik        = -opt$value,
    deviance      = deviance,
    null_deviance = null_deviance,
    nobs          = length(residuals),
    k             = 3L,
    sigma_hat     = sigma_hat,
    cov           = cov,
    se_a1         = unname(se_a1),
    se_a2         = unname(se_a2),
    convergence   = opt$convergence
  )
}

#' Maximum-likelihood point fit of the lrmsd_n MSA model (mode form)
#'
#' Mode counterpart of [fit_lrmsd_i_msa_ml()]: maximises the mode-form profiled
#' Gaussian log-likelihood (`calculate_loglik_lrmsd_n_msa()`) over `(a1, a2)` by
#' numerical optimisation, returning a point estimate plus an asymptotic covariance
#' from the Hessian at the optimum. Identical machinery to the site fit; the response
#' index is the mode `n` (no `site_map` / `pdb_site`).
#'
#' The optimiser works in the same coordinates as the site fit: `a1` and
#' `log2(a2 + 1)` (so `a2 = 2^(log2(a2+1)) - 1 >= 0`), on the box `a1_range` ×
#' `log2_a2_plus1_range` — the coordinates in which the prior is uniform. The
#' returned covariance `cov` is on the `(a1, log2(a2+1))` scale; the standard error
#' of `a2` is obtained by the delta method (`da2/d(log2(a2+1)) = 2^(log2(a2+1)) * ln 2`).
#'
#' @param spm_pp_mode Preprocessed data from [preprocess_spm_mode()] (energy_data +
#'   the `dr2_njm` matrix; no `site_map`).
#' @param observed_data Tibble with columns `n` (mode index) and `lrmsd_n_obs` (the
#'   fit target), as documented for `calculate_loglik_lrmsd_n_msa()`.
#' @param a1_range Length-2 `[min, max]` box bound for `a1`.
#' @param log2_a2_plus1_range Length-2 `[min, max]` box bound for `log2(a2 + 1)`.
#' @param init Optional length-2 numeric start `c(a1, log2(a2+1))`. When `NULL`
#'   (default), a deterministic coarse grid-max of the likelihood over the box is
#'   used as the start (robust against a bad local start; cheap).
#' @param grid_n Number of points per axis for the default grid-max start (ignored
#'   when `init` is supplied).
#' @return A list with the point estimate and asymptotic uncertainty, identical in
#'   shape to [fit_lrmsd_i_msa_ml()]:
#'   \describe{
#'     \item{a1, a2}{Point estimate of the stability (`a1`) and activity (`a2`)
#'       selection strengths, on the natural scale (the paper's `aS`, `aA`).}
#'     \item{logLik}{Profiled Gaussian log-likelihood at the optimum.}
#'     \item{deviance}{Residual deviance `sum(residuals^2)` (= `nobs * sigma_hat^2`).
#'       A goodness-of-fit primitive; consumed by [gof_lrmsd_n_msa_ml()].}
#'     \item{null_deviance}{Deviance of the flat/mean-only null (a constant profile),
#'       `sum((observed - mean(observed))^2)`. With `deviance`, gives
#'       `D2 = 1 - deviance/null_deviance`.}
#'     \item{nobs}{Number of matched (observed) modes the likelihood scores.}
#'     \item{k}{Free-parameter count for AIC/BIC: `a1`, `a2`, and the profiled `sigma`
#'       (`3`), matching the `logLik.lm`/`broom` convention.}
#'     \item{sigma_hat}{Profiled noise scale `sqrt(mean(residuals^2))` at the optimum.}
#'     \item{cov}{2×2 covariance matrix on the `(a1, log2(a2+1))` scale.}
#'     \item{se_a1, se_a2}{Standard errors; `se_a2` via the delta method.}
#'     \item{convergence}{`optim` convergence code (0 = success).}
#'   }
#' @seealso [fit_lrmsd_i_msa_ml()] (the site counterpart),
#'   [gof_lrmsd_n_msa_ml()] (goodness-of-fit from this fit),
#'   `calculate_loglik_lrmsd_n_msa()` (the objective).
#' @family fitting
#' @export
#' @examples
#' \dontrun{
#' pp <- preprocess_spm_mode(znb_spm)
#' ml <- fit_lrmsd_n_msa_ml(pp, znb_profile_n)
#' c(a1 = ml$a1, a2 = ml$a2)
#' }
fit_lrmsd_n_msa_ml <- function(spm_pp_mode,
                       observed_data,
                       a1_range = c(0, 10),
                       log2_a2_plus1_range = c(0, 13),
                       init = NULL,
                       grid_n = 25) {
  # Validate box bounds (same contract as fit_lrmsd_i_msa_ml) -- fail loud.
  if (length(a1_range) != 2 || a1_range[1] >= a1_range[2]) {
    stop("a1_range must be a vector of length 2 with min < max")
  }
  if (length(log2_a2_plus1_range) != 2 ||
      log2_a2_plus1_range[1] >= log2_a2_plus1_range[2]) {
    stop("log2_a2_plus1_range must be a vector of length 2 with min < max")
  }

  lower <- c(a1_range[1], log2_a2_plus1_range[1])
  upper <- c(a1_range[2], log2_a2_plus1_range[2])

  # Negative profiled log-likelihood in (a1, b) coordinates, b = log2(a2 + 1).
  nll <- function(theta) {
    -calculate_loglik_lrmsd_n_msa(spm_pp_mode, observed_data,
                          a1 = theta[1], a2 = 2^theta[2] - 1)
  }

  # Start: caller-supplied, else a deterministic coarse grid-max over the box.
  if (is.null(init)) {
    a1_grid <- seq(lower[1], upper[1], length.out = grid_n)
    b_grid  <- seq(lower[2], upper[2], length.out = grid_n)
    grid    <- expand.grid(a1 = a1_grid, b = b_grid)
    ll      <- apply(grid, 1L, function(row) -nll(c(row[["a1"]], row[["b"]])))
    init    <- as.numeric(grid[which.max(ll), c("a1", "b")])
  } else {
    if (length(init) != 2) stop("init must be a length-2 numeric c(a1, log2(a2+1))")
    if (init[1] < lower[1] || init[1] > upper[1] ||
        init[2] < lower[2] || init[2] > upper[2]) {
      stop("init must lie within the box [a1_range] x [log2_a2_plus1_range]")
    }
  }

  opt <- optim(init, nll, method = "L-BFGS-B", lower = lower, upper = upper)

  par_fit <- opt$par
  a1_hat  <- par_fit[1]
  b_hat   <- par_fit[2]
  a2_hat  <- 2^b_hat - 1

  # Asymptotic covariance from the Hessian of the NLL at the optimum.
  H <- optimHess(par_fit, nll)
  cov <- tryCatch(solve(H), error = function(e) NULL)
  if (is.null(cov) || any(!is.finite(cov)) || any(diag(cov) <= 0)) {
    stop("Hessian at the ML optimum is singular or not positive-definite; ",
         "cannot form the asymptotic covariance. The likelihood may be flat in ",
         "one direction (e.g. a parameter pinned at a box bound).")
  }
  se <- sqrt(diag(cov))
  se_a1 <- se[1]
  se_b  <- se[2]
  # Delta method: a2 = 2^b - 1, da2/db = 2^b * ln 2.
  se_a2 <- abs(2^b_hat * log(2)) * se_b

  # Profiled noise scale at the optimum (mode form: join on n, no site_map).
  pred <- calculate_lrmsd_n_msa(spm_pp_mode, a1_hat, a2_hat)
  obs <- observed_data %>%
    dplyr::select(n, lrmsd_n_obs)
  cmp <- obs %>%
    inner_join(pred, by = "n") %>%
    mutate(
      nlrmsd_n_obs = lrmsd_n_obs - mean(lrmsd_n_obs),
      nlrmsd_n_msa = lrmsd_n_msa - mean(lrmsd_n_msa)
    )
  residuals <- cmp$nlrmsd_n_obs - cmp$nlrmsd_n_msa
  sigma_hat <- sqrt(mean(residuals^2))

  # Goodness-of-fit primitives (flat/mean-only null; k counts sigma). D^2/AIC/BIC are
  # derived from these by gof_lrmsd_n_msa_ml(). deviance = sum(resid^2) (= glm's
  # residual deviance for a Gaussian, = nobs * sigma_hat^2); null_deviance vs a
  # constant profile.
  deviance      <- sum(residuals^2)
  null_deviance <- calculate_null_deviance(cmp$lrmsd_n_obs)

  list(
    a1            = unname(a1_hat),
    a2            = unname(a2_hat),
    logLik        = -opt$value,
    deviance      = deviance,
    null_deviance = null_deviance,
    nobs          = length(residuals),
    k             = 3L,
    sigma_hat     = sigma_hat,
    cov           = cov,
    se_a1         = unname(se_a1),
    se_a2         = unname(se_a2),
    convergence   = opt$convergence
  )
}

# Flat/mean-only null deviance: the total sum of squares of the observed profile about
# its mean, sum((y - mean(y))^2) (glm's null.deviance for a Gaussian; = (n-1)*var(y)).
# It is the deviance of the best constant prediction, so it depends only on the data, not
# on any model -- the same null for every model fit to that data, which is what makes
# D^2 = 1 - deviance/null_deviance comparable. Centres internally, so the result is the
# same whether the caller passes raw or already-centred y. Axis-agnostic, pure.
#' @noRd
calculate_null_deviance <- function(y) {
  n <- length(y)
  if (n < 1L) stop("y must be non-empty")
  nd <- sum((y - mean(y))^2)
  if (nd <= 0) stop("null deviance is zero (observations are all equal)")
  nd
}

# ---- fitting objective: profiled Gaussian log-likelihood (internal) --------------
# The criterion the fitter optimises. Internal (@noRd): users never call the raw
# likelihood; the goodness-of-fit surface (AIC etc.) is separate public functions.
# Kept a distinct function so the criterion stays swappable (a future
# stochastic-likelihood tree plugs in here without touching the fitter).

# Log-likelihood of an observed per-site divergence profile at (a1, a2): mean-centre
# model prediction and observations, compare under a Gaussian noise model whose scale
# is profiled out (sigma = sqrt(mean(r^2))). Maximised by fit_lrmsd_i_msa_ml.
# spm_pp: preprocess_spm() output; observed_data:
# a tibble {pdb_site, lrmsd_i_obs}. Returns one numeric log-likelihood.
#' @noRd
calculate_loglik_lrmsd_i_msa <- function(spm_pp, observed_data, a1, a2) {

  # Generate model predictions (keyed by the internal response-site index i)
  predictions <- calculate_lrmsd_i_msa(spm_pp, a1, a2)

  # Translate user-supplied pdb_site to the internal index i via the model's
  # site_map. pdb_site is the structure-anchored key; i is model-internal.
  site_map <- spm_pp$site_map
  observations <- observed_data %>%
    dplyr::select(pdb_site, lrmsd_i_obs)

  unknown <- setdiff(observations$pdb_site, site_map$pdb_site)
  if (length(unknown) > 0) {
    stop("observed_data has pdb_site value(s) not present in the model: ",
         paste(unknown, collapse = ", "))
  }

  observations <- observations %>%
    inner_join(site_map, by = "pdb_site") %>%
    dplyr::select(i, lrmsd_i_obs)

  # Match predictions with observations
  comparison <- observations %>%
    inner_join(predictions, by = "i") %>%
    mutate(
      nlrmsd_i_obs = lrmsd_i_obs - mean(lrmsd_i_obs),
      nlrmsd_i_msa = lrmsd_i_msa - mean(lrmsd_i_msa)
    )

  # Calculate residuals
  residuals <- comparison$nlrmsd_i_obs - comparison$nlrmsd_i_msa

  # Estimate sigma from residuals. sigma is profiled out at each (a1, a2): the
  # value below is the closed-form ML estimate for r_i ~ N(0, sigma^2), namely
  # sqrt(mean(r^2)) (divisor n). (Previously sd(residuals), divisor n-1, which is
  # not the profile MLE; the (a1, a2) argmax is unchanged but sigma-derived
  # quantities — logLik, SEs — were slightly off.)
  sigma <- sqrt(mean(residuals^2))

  # Calculate log-likelihood using residuals
  log_lik <- sum(dnorm(residuals, 0, sigma, log = TRUE))

  return(log_lik)
}

# Mode-form counterpart of calculate_loglik_lrmsd_i_msa: scores over normal modes,
# no pdb_site mapping (mode index n used directly). spm_pp_mode: preprocess_spm_mode()
# output; observed_data: a tibble {n, lrmsd_n_obs}. Maximised by fit_lrmsd_n_msa_ml.
#' @noRd
calculate_loglik_lrmsd_n_msa <- function(spm_pp_mode, observed_data, a1, a2) {

  # Generate model predictions (keyed by the response-mode index n)
  predictions <- calculate_lrmsd_n_msa(spm_pp_mode, a1, a2)

  # Modes are not structure-anchored: n is the model index directly (no site_map).
  observations <- observed_data %>%
    dplyr::select(n, lrmsd_n_obs)

  unknown <- setdiff(observations$n, predictions$n)
  if (length(unknown) > 0) {
    stop("observed_data has mode index(es) not present in the model: ",
         paste(unknown, collapse = ", "))
  }

  # Match predictions with observations
  comparison <- observations %>%
    inner_join(predictions, by = "n") %>%
    mutate(
      nlrmsd_n_obs = lrmsd_n_obs - mean(lrmsd_n_obs),
      nlrmsd_n_msa = lrmsd_n_msa - mean(lrmsd_n_msa)
    )

  # Calculate residuals
  residuals <- comparison$nlrmsd_n_obs - comparison$nlrmsd_n_msa

  # Profiled-out sigma: closed-form ML estimate sqrt(mean(r^2)) (divisor n), same
  # as the site form (see calculate_loglik_lrmsd_i_msa).
  sigma <- sqrt(mean(residuals^2))

  log_lik <- sum(dnorm(residuals, 0, sigma, log = TRUE))

  return(log_lik)
}
