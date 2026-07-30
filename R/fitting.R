# MSA model fitting — maximum-likelihood point estimation
# Maximises the profiled Gaussian log-likelihood (calculate_loglik_lrmsd_i_msa) over
# (a1, log2(a2+1)) and returns a point estimate plus an asymptotic covariance.

# ---- axis-blind fitting core ------------------------------------------------------
# The site (i, keyed to pdb_site via site_map) and mode (n, keyed directly) fits are the
# same machinery on a resolved (index, obs) pair. These cores hold that shared machinery;
# they receive ALREADY-RESOLVED, ALREADY-VALIDATED predictions and observations (both bare
# two-column tibbles keyed by the same integer index column named `idx`), so they carry no
# notion of site vs mode, no site_map, and no unknown-key check -- those stay in the
# exported wrappers (the wrapper owns key resolution + its own unknown-key error).

#' Match a resolved prediction/observation pair and mean-centre both
#'
#' Inner-joins predictions to observations on the shared integer key `idx`, then centres
#' each over the MATCHED support (its own mean over the joined rows). Axis-blind: the
#' caller has already resolved any label->index mapping and named both value columns
#' generically. The matched-set centring is the fit's centring domain (distinct from the
#' full-support centring of the calculate_nlrmsd_* forward maps).
#'
#' @param predictions Tibble `(idx, pred)`: model lrmsd keyed by the internal index.
#' @param obs Tibble `(idx, obs)`: observed lrmsd keyed by the same index.
#' @return A tibble of the matched rows with `obs`, `pred`, and their centred forms
#'   `obs_c`, `pred_c` (residuals = `obs_c - pred_c`).
#' @noRd
match_lrmsd_obs_pred <- function(predictions, obs) {
  obs %>%
    inner_join(predictions, by = "idx") %>%
    mutate(obs_c = obs - mean(obs), pred_c = pred - mean(pred))
}

#' Axis-blind profiled Gaussian log-likelihood from a resolved (pred, obs) pair
#'
#' The shared core of `calculate_loglik_lrmsd_{i,n}_msa()`: mean-centre matched model and
#' observations, compare under a Gaussian whose scale is profiled out
#' (`sigma = sqrt(mean(r^2))`). Receives resolved, validated frames; performs no key
#' resolution or unknown-key check.
#'
#' @param predictions Tibble `(idx, pred)`.
#' @param obs Tibble `(idx, obs)`.
#' @return A single numeric log-likelihood.
#' @noRd
loglik_lrmsd_msa_core <- function(predictions, obs) {
  cmp <- match_lrmsd_obs_pred(predictions, obs)
  residuals <- cmp$obs_c - cmp$pred_c
  sigma <- sqrt(mean(residuals^2))
  sum(dnorm(residuals, 0, sigma, log = TRUE))
}

#' Axis-blind goodness-of-fit primitives from a resolved (pred, obs) pair
#'
#' The shared sigma/deviance block of `fit_lrmsd_{i,n}_msa_ml()`: matched residuals, the
#' profiled noise scale, residual deviance, and the flat-null deviance. Axis-blind.
#'
#' @param predictions Tibble `(idx, pred)` at the fit's point estimate.
#' @param obs Tibble `(idx, obs)`.
#' @return A list `sigma_hat`, `deviance`, `null_deviance`, `nobs`.
#' @noRd
fit_gof_primitives <- function(predictions, obs) {
  cmp <- match_lrmsd_obs_pred(predictions, obs)
  residuals <- cmp$obs_c - cmp$pred_c
  list(
    sigma_hat     = sqrt(mean(residuals^2)),
    deviance      = sum(residuals^2),
    null_deviance = calculate_null_deviance(cmp$obs),
    nobs          = length(residuals)
  )
}

#' Axis-blind ML optimisation + asymptotic covariance
#'
#' The shared machinery of `fit_lrmsd_{i,n}_msa_ml()`: validate the box, build the
#' grid-max start, optimise the supplied negative log-likelihood in `(a1, log2(a2+1))`
#' coordinates, form the asymptotic covariance from the Hessian, and assemble the return
#' list. The axis-specific likelihood and goodness-of-fit are injected as closures.
#'
#' @param nll Negative profiled log-likelihood, a function of `theta = c(a1, log2(a2+1))`.
#' @param gof_fn A function of `(a1_hat, a2_hat)` returning the `fit_gof_primitives()` list
#'   at the optimum.
#' @param a1_range,log2_a2_plus1_range,init,grid_n The box + start controls (see the
#'   exported fits).
#' @return The fit list (a1, a2, logLik, deviance, null_deviance, nobs, k, sigma_hat, cov,
#'   se_a1, se_a2, convergence).
#' @noRd
fit_lrmsd_msa_ml_core <- function(nll, gof_fn,
                                  a1_range, log2_a2_plus1_range,
                                  init, grid_n) {
  if (length(a1_range) != 2 || a1_range[1] >= a1_range[2]) {
    stop("a1_range must be a vector of length 2 with min < max")
  }
  if (length(log2_a2_plus1_range) != 2 ||
      log2_a2_plus1_range[1] >= log2_a2_plus1_range[2]) {
    stop("log2_a2_plus1_range must be a vector of length 2 with min < max")
  }

  lower <- c(a1_range[1], log2_a2_plus1_range[1])
  upper <- c(a1_range[2], log2_a2_plus1_range[2])

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
  se_a2 <- abs(2^b_hat * log(2)) * se_b   # delta: a2 = 2^b - 1, da2/db = 2^b * ln 2

  gof <- gof_fn(a1_hat, a2_hat)

  list(
    a1            = unname(a1_hat),
    a2            = unname(a2_hat),
    logLik        = -opt$value,
    deviance      = gof$deviance,
    null_deviance = gof$null_deviance,
    nobs          = gof$nobs,
    k             = 3L,
    sigma_hat     = gof$sigma_hat,
    cov           = cov,
    se_a1         = unname(se_a1),
    se_a2         = unname(se_a2),
    convergence   = opt$convergence
  )
}

# ---- axis-specific observation resolvers (the ONLY axis-aware boundary logic) ------
# Each turns the user's observed_data into the canonical (idx, obs) frame the cores
# consume, and owns its axis's unknown-key error. SITE: pdb_site -> internal index i via
# site_map. MODE: n is the index directly (no map). These are the sole locus of the
# site/mode asymmetry in fitting.

#' Resolve site observations (pdb_site) to the internal index (site boundary)
#'
#' Translates user `pdb_site` labels to the model-internal site index `i` via
#' `spm$site_map`, erroring on any `pdb_site` not present in the model. Returns the
#' canonical `(idx, obs)` frame the fitting cores consume.
#'
#' @param spm The `spm` object (its `site_map` keys pdb_site -> i).
#' @param observed_data Tibble with `pdb_site` and `lrmsd_i_obs`.
#' @return A tibble `(idx, obs)`.
#' @noRd
resolve_site_obs <- function(spm, observed_data) {
  observations <- observed_data %>% dplyr::select(pdb_site, lrmsd_i_obs)
  unknown <- setdiff(observations$pdb_site, spm$site_map$pdb_site)
  if (length(unknown) > 0) {
    stop("observed_data has pdb_site value(s) not present in the model: ",
         paste(unknown, collapse = ", "))
  }
  observations %>%
    inner_join(spm$site_map, by = "pdb_site") %>%
    dplyr::transmute(idx = i, obs = lrmsd_i_obs)
}

#' Resolve mode observations (n) to the internal index (mode boundary)
#'
#' Modes are not structure-anchored: the mode index `n` IS the internal index. Errors on
#' any observed `n` not present in the model's predictions. Returns the canonical
#' `(idx, obs)` frame.
#'
#' @param predictions The mode predictions tibble carrying `n` (the valid index set).
#' @param observed_data Tibble with `n` and `lrmsd_n_obs`.
#' @return A tibble `(idx, obs)`.
#' @noRd
resolve_mode_obs <- function(predictions, observed_data) {
  observations <- observed_data %>% dplyr::select(n, lrmsd_n_obs)
  unknown <- setdiff(observations$n, predictions$n)
  if (length(unknown) > 0) {
    stop("observed_data has mode index(es) not present in the model: ",
         paste(unknown, collapse = ", "))
  }
  observations %>% dplyr::transmute(idx = n, obs = lrmsd_n_obs)
}

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
#' @param spm A single-point-mutation `spm` object from [generate_spm_data()] (its `site_map` keys the fit to PDB residues).
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
#' spm <- generate_spm_data(znb_wt, pdb_site_active = c(99,101,103,162,181,184,193,223), seed = 1024)
#' ml <- fit_lrmsd_i_msa_ml(spm, znb_profile)
#' c(a1 = ml$a1, a2 = ml$a2)
#' }
fit_lrmsd_i_msa_ml <- function(spm,
                       observed_data,
                       a1_range = c(0, 10),
                       log2_a2_plus1_range = c(0, 13),
                       init = NULL,
                       grid_n = 25) {
  # Site boundary: the nll and the sigma/GoF block both resolve observations pdb_site -> i
  # via site_map (with the site-axis unknown-key error). The axis-blind core owns the
  # optimisation, covariance, and return assembly.
  nll <- function(theta) {
    -calculate_loglik_lrmsd_i_msa(spm, observed_data,
                          a1 = theta[1], a2 = 2^theta[2] - 1)
  }
  gof_fn <- function(a1_hat, a2_hat) {
    obs  <- resolve_site_obs(spm, observed_data)
    pred <- calculate_lrmsd_i_msa(spm, a1_hat, a2_hat) %>%
      dplyr::transmute(idx = i, pred = lrmsd_i_msa)
    fit_gof_primitives(pred, obs)
  }
  fit_lrmsd_msa_ml_core(nll, gof_fn, a1_range, log2_a2_plus1_range, init, grid_n)
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
#' @param spm A single-point-mutation `spm` object from [generate_spm_data()] (energy_data +
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
#' spm <- generate_spm_data(znb_wt, pdb_site_active = c(99,101,103,162,181,184,193,223), seed = 1024)
#' ml <- fit_lrmsd_n_msa_ml(spm, znb_profile_n)
#' c(a1 = ml$a1, a2 = ml$a2)
#' }
fit_lrmsd_n_msa_ml <- function(spm,
                       observed_data,
                       a1_range = c(0, 10),
                       log2_a2_plus1_range = c(0, 13),
                       init = NULL,
                       grid_n = 25) {
  # Mode boundary: n is the internal index directly (no site_map). Same axis-blind core;
  # the mode unknown-key check is on n, done inside resolve_mode_obs.
  nll <- function(theta) {
    -calculate_loglik_lrmsd_n_msa(spm, observed_data,
                          a1 = theta[1], a2 = 2^theta[2] - 1)
  }
  gof_fn <- function(a1_hat, a2_hat) {
    pred <- calculate_lrmsd_n_msa(spm, a1_hat, a2_hat)
    obs  <- resolve_mode_obs(pred, observed_data)
    fit_gof_primitives(dplyr::transmute(pred, idx = n, pred = lrmsd_n_msa), obs)
  }
  fit_lrmsd_msa_ml_core(nll, gof_fn, a1_range, log2_a2_plus1_range, init, grid_n)
}

#' Flat/mean-only null deviance of an observed profile
#'
#' The total sum of squares of the observed profile about its mean,
#' `sum((y - mean(y))^2)` (glm's `null.deviance` for a Gaussian; `= (n-1)*var(y)`). It
#' is the deviance of the best constant prediction, so it depends only on the data, not
#' on any model -- the same null for every model fit to that data, which is what makes
#' `D^2 = 1 - deviance/null_deviance` comparable. Centres internally, so the result is
#' the same whether the caller passes raw or already-centred `y`. Axis-agnostic, pure.
#'
#' @param y Numeric vector of observed profile values (raw or already mean-centred).
#' @return A single numeric: the null deviance. Stops if `y` is empty or all-equal.
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

#' Profiled Gaussian log-likelihood of a per-site divergence profile
#'
#' Log-likelihood of an observed per-site divergence profile at `(a1, a2)`: mean-centre
#' the model prediction and observations, compare under a Gaussian noise model whose
#' scale is profiled out (`sigma = sqrt(mean(r^2))`). Maximised by [fit_lrmsd_i_msa_ml()].
#'
#' @param spm A single-point-mutation `spm` object from [generate_spm_data()] (its `site_map` keys the fit to PDB residues).
#' @param observed_data Tibble with columns `pdb_site` and `lrmsd_i_obs`.
#' @param a1 Stability selection strength.
#' @param a2 Activity selection strength.
#' @return A single numeric log-likelihood.
#' @noRd
calculate_loglik_lrmsd_i_msa <- function(spm, observed_data, a1, a2) {
  # Resolve user pdb_site -> internal site index i (site-only, at the boundary), with the
  # site-axis unknown-key check; then hand canonical (idx, obs)/(idx, pred) to the core.
  obs         <- resolve_site_obs(spm, observed_data)             # tibble(idx, obs)
  predictions <- calculate_lrmsd_i_msa(spm, a1, a2) %>%
    dplyr::transmute(idx = i, pred = lrmsd_i_msa)
  loglik_lrmsd_msa_core(predictions, obs)
}

#' Profiled Gaussian log-likelihood of a per-mode divergence profile
#'
#' Mode-form counterpart of `calculate_loglik_lrmsd_i_msa()`: scores over normal modes,
#' with no `pdb_site` mapping (the mode index `n` is used directly). Maximised by
#' [fit_lrmsd_n_msa_ml()].
#'
#' @param spm A single-point-mutation `spm` object from [generate_spm_data()].
#' @param observed_data Tibble with columns `n` (mode index) and `lrmsd_n_obs`.
#' @param a1 Stability selection strength.
#' @param a2 Activity selection strength.
#' @return A single numeric log-likelihood.
#' @noRd
calculate_loglik_lrmsd_n_msa <- function(spm, observed_data, a1, a2) {
  # Mode index n is the internal index directly (no site_map); resolve + unknown-key check
  # at the boundary, then hand canonical frames to the axis-blind core.
  predictions <- calculate_lrmsd_n_msa(spm, a1, a2)
  obs         <- resolve_mode_obs(predictions, observed_data)     # tibble(idx, obs)
  loglik_lrmsd_msa_core(dplyr::transmute(predictions, idx = n, pred = lrmsd_n_msa), obs)
}
