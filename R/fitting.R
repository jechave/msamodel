# MSA model fitting — maximum-likelihood point estimation
#
# The whole job: find the (a1, a2) that make the model's predicted divergence profile
# match an observed one, then report how well they match.
#
# Four functions, in dependency order — each one a reduction of the one above:
#   residuals_lrmsd_msa()  matched, mean-centred residuals at one (a1, a2)
#   loglik_lrmsd_msa()     those residuals -> profiled Gaussian log-likelihood (the criterion)
#   gof_lrmsd_msa()        those residuals -> the goodness-of-fit row
#   fit_lrmsd_msa()        optimise the criterion, then build the fit object
#
# All four are AXIS-BLIND: they take a bare [mutant x response] divergence matrix and an
# ALREADY-RESOLVED (idx, obs) frame, so they carry no notion of site vs mode, no site_map,
# and no unknown-key check. Resolving observations to the internal index is the exported
# wrappers' job (bottom of this file), done ONCE before optimising -- so the unknown-key
# error fires at the boundary rather than on every likelihood evaluation.
#
# The wrappers pass VALUES (a matrix, a resolved frame), not functions. An earlier design
# injected `nll` and `gof_fn` closures into a "core"; that hid the fact that the
# likelihood and the goodness of fit are the same residuals reduced two ways, and it
# recomputed them on a second path at the optimum.

#' Matched, mean-centred residuals at one (a1, a2)
#'
#' The single definition of "how far the model is from the data" -- everything else in the
#' fit is a reduction of this. Evaluates the forward map at `(a1, a2)`, inner-joins it to
#' the observations on the internal index, and centres BOTH over the matched support (each
#' by its own mean over the joined rows).
#'
#' The matched-support centring is the fit's centring domain, and it differs deliberately
#' from the full-support centring of `nlrmsd_msa()`: the fit compares only the responses
#' actually observed, whereas prediction centres over every response in the model.
#'
#' @param dr2mat A `[mutant x response]` divergence matrix (`dr2mat_site` or `dr2mat_mode`).
#' @param energy_data The per-mutant energy tibble (for the fixation weights).
#' @param obs Tibble `(idx, obs)`: observations already keyed by the internal index.
#' @param a1,a2 Selection strengths.
#' @return A list: `resid` (numeric vector, `obs - pred` after centring) and `obs_matched`
#'   (the raw observed values on the matched rows, for the null deviance).
#' @noRd
residuals_lrmsd_msa <- function(dr2mat, energy_data, obs, a1, a2) {
  v    <- lrmsd_msa(dr2mat, energy_data, a1, a2)
  cmp  <- inner_join(obs, tibble::tibble(idx = seq_along(v), pred = v), by = "idx")
  list(resid       = (cmp$obs - mean(cmp$obs)) - (cmp$pred - mean(cmp$pred)),
       obs_matched = cmp$obs)
}

#' Profiled Gaussian log-likelihood at one (a1, a2)
#'
#' The criterion the fitter maximises. Compares matched, centred model and observations
#' under a Gaussian whose scale is profiled out (`sigma = sqrt(mean(r^2))`), so the only
#' free parameters are `(a1, a2)`.
#'
#' @inheritParams residuals_lrmsd_msa
#' @return A single numeric log-likelihood.
#' @noRd
loglik_lrmsd_msa <- function(dr2mat, energy_data, obs, a1, a2) {
  r     <- residuals_lrmsd_msa(dr2mat, energy_data, obs, a1, a2)$resid
  sigma <- sqrt(mean(r^2))
  sum(dnorm(r, 0, sigma, log = TRUE))
}

#' Goodness-of-fit row from the residuals at the optimum
#'
#' Turns the same residuals the likelihood scores into the reportable summary. Attached by
#' the fitter as `fit$gof`; there is no accessor to call.
#'
#' `D2` is deviance-explained (`= 1 - Var(resid)/Var(obs)`, since both deviances are sums
#' of squares on the same `n`). It is at most `1` but has **no lower bound**: negative
#' means the prediction is worse than a flat/mean-only null. Returned unclamped -- a
#' negative `D2` is a real signal that the fit is poor, not an error. `AIC`/`BIC` are
#' per-fit numbers, meaningful only *compared* against another model fit at its own maximum.
#'
#' @param r The `residuals_lrmsd_msa()` list at the fit's point estimate.
#' @param logLik The profiled Gaussian log-likelihood at that point.
#' @param k Free-parameter count for AIC/BIC.
#' @return A one-row tibble: `D2`, `AIC`, `BIC`, `logLik`, `deviance`, `null_deviance`,
#'   `nobs`, `k`, `sigma_hat`.
#' @noRd
gof_lrmsd_msa <- function(r, logLik, k) {
  deviance      <- sum(r$resid^2)
  null_deviance <- calculate_null_deviance(r$obs_matched)
  nobs          <- length(r$resid)
  tibble::tibble(
    D2            = 1 - deviance / null_deviance,
    AIC           = -2 * logLik + 2 * k,
    BIC           = -2 * logLik + k * log(nobs),
    logLik        = logLik,
    deviance      = deviance,
    null_deviance = null_deviance,
    nobs          = nobs,
    k             = k,
    sigma_hat     = sqrt(mean(r$resid^2))
  )
}

#' Fit the MSA model to a resolved observation frame
#'
#' The shared fitter both exported wrappers delegate to, top to bottom: validate the box,
#' optimise, form the asymptotic covariance, build the goodness-of-fit row, return.
#'
#' Optimisation runs in the coordinates `theta = (a1, log2(a2 + 1))` -- the coordinates in
#' which the prior is uniform, and in which `a2 = 2^theta2 - 1 >= 0` holds by construction.
#'
#' @inheritParams residuals_lrmsd_msa
#' @param a1_range,log2_a2_plus1_range Length-2 `[min, max]` box bounds.
#' @param init Optional length-2 unnamed numeric start `c(a1, log2(a2+1))`; `NULL` uses a
#'   deterministic coarse grid-max over the box.
#' @param grid_n Points per axis for the default grid-max start.
#' @param call The exported wrapper's `match.call()`, stored as `$call`.
#' @return The fit list (see [fit_lrmsd_msa_site()] for the documented shape).
#' @noRd
fit_lrmsd_msa <- function(dr2mat, energy_data, obs,
                          a1_range, log2_a2_plus1_range, init, grid_n, call) {
  # --- 1. the objective, and the box it is optimised over ---------------------------
  nll <- function(theta) -loglik_lrmsd_msa(dr2mat, energy_data, obs,
                                           a1 = theta[1], a2 = 2^theta[2] - 1)

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
    grid <- expand.grid(theta1 = seq(lower[1], upper[1], length.out = grid_n),
                        theta2 = seq(lower[2], upper[2], length.out = grid_n))
    ll   <- apply(grid, 1L, function(row) -nll(c(row[["theta1"]], row[["theta2"]])))
    init <- as.numeric(grid[which.max(ll), c("theta1", "theta2")])
  } else {
    if (length(init) != 2) stop("init must be a length-2 numeric c(a1, log2(a2+1))")
    if (!is.null(names(init))) stop("init must be an unnamed positional c(a1, log2(a2+1))")
    if (init[1] < lower[1] || init[1] > upper[1] ||
        init[2] < lower[2] || init[2] > upper[2]) {
      stop("init must lie within the box [a1_range] x [log2_a2_plus1_range]")
    }
  }

  # --- 2. optimise ------------------------------------------------------------------
  opt <- optim(init, nll, method = "L-BFGS-B", lower = lower, upper = upper)

  theta_hat <- opt$par
  a1_hat    <- theta_hat[1]                 # theta1 IS a1 (identity)
  a2_hat    <- 2^theta_hat[2] - 1

  # --- 3. asymptotic covariance, and SEs on the natural scale -----------------------
  H   <- optimHess(theta_hat, nll)
  cov <- tryCatch(solve(H), error = function(e) NULL)
  if (is.null(cov) || any(!is.finite(cov)) || any(diag(cov) <= 0)) {
    stop("Hessian at the ML optimum is singular or not positive-definite; ",
         "cannot form the asymptotic covariance. The likelihood may be flat in ",
         "one direction (e.g. a parameter pinned at a box bound).")
  }
  se    <- sqrt(diag(cov))
  se_a1 <- se[1]                                     # da1/dtheta1 = 1
  se_a2 <- abs(2^theta_hat[2] * log(2)) * se[2]      # delta method: da2/dtheta2

  # --- 4. goodness of fit, from the residuals at the optimum ------------------------
  logLik <- -opt$value
  gof    <- gof_lrmsd_msa(residuals_lrmsd_msa(dr2mat, energy_data, obs, a1_hat, a2_hat),
                          logLik = logLik,
                          k      = 3L)   # a1, a2, and the profiled sigma

  # --- 5. the fit object: estimate at the top, fit quality in $gof -------------------
  list(
    a1          = a1_hat,
    a2          = a2_hat,
    logLik      = logLik,
    cov         = cov,
    se_a1       = se_a1,
    se_a2       = se_a2,
    convergence = opt$convergence,
    gof         = gof,
    call        = call
  )
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


# ---- axis-specific observation resolvers (the ONLY axis-aware boundary logic) ------
# Each turns the user's two observation vectors into the canonical (idx, obs) frame the
# cores consume, and owns its axis's unknown-key error. SITE: pdb_site -> internal index
# via site_map. MODE: the mode number IS the index (no map). These are the sole locus of
# the site/mode asymmetry in fitting.
#
# Observations arrive as a (key, value) VECTOR PAIR, not a tibble: the package never
# looks up a column name in the user's data, so the caller's own frame can name its
# columns whatever it likes. It also makes length disagreement checkable -- a tibble's
# columns are equal-length by construction, so a misaligned pairing is invisible there.

#' Validate a (key, value) observation vector pair (shared boundary check)
#'
#' Fail-loud gate both resolvers run before doing anything: equal lengths, non-empty, and
#' no `NA` in either vector. Axis-blind; `key_name` only shapes the message.
#'
#' @param key The observation key vector (pdb_site labels, or mode numbers).
#' @param value The observed lrmsd vector.
#' @param key_name Name of the key argument, for the error message.
#' @return Invisibly `NULL`; stops with a clear message on any violation.
#' @noRd
check_obs_vectors <- function(key, value, key_name) {
  if (length(key) != length(value)) {
    stop("`", key_name, "` and `lrmsd_obs` must have the same length; got ",
         length(key), " and ", length(value), ".")
  }
  if (length(key) == 0L) {
    stop("`", key_name, "` and `lrmsd_obs` must be non-empty.")
  }
  if (anyNA(key))   stop("`", key_name, "` must not contain NA.")
  if (anyNA(value)) stop("`lrmsd_obs` must not contain NA.")
  invisible(NULL)
}

#' Resolve site observations (pdb_site) to the internal index (site boundary)
#'
#' Translates user `pdb_site` labels to the model-internal site index via
#' `spm$site_map`, erroring on any `pdb_site` not present in the model. Returns the
#' canonical `(idx, obs)` frame the fitting cores consume.
#'
#' @param spm The `spm` object (its `site_map` keys pdb_site -> site).
#' @param pdb_site Integer vector of PDB residue numbers.
#' @param lrmsd_obs Numeric vector of observed log divergences, same length.
#' @return A tibble `(idx, obs)`.
#' @noRd
resolve_site_obs <- function(spm, pdb_site, lrmsd_obs) {
  check_obs_vectors(pdb_site, lrmsd_obs, "pdb_site")
  unknown <- setdiff(pdb_site, spm$site_map$pdb_site)
  if (length(unknown) > 0) {
    stop("observed data has pdb_site value(s) not present in the model: ",
         paste(unknown, collapse = ", "))
  }
  tibble::tibble(pdb_site = pdb_site, obs = lrmsd_obs) %>%
    inner_join(spm$site_map, by = "pdb_site") %>%
    dplyr::transmute(idx = site, obs = obs)
}

#' Resolve mode observations to the internal index (mode boundary)
#'
#' Modes are not structure-anchored: the mode number IS the internal index. Errors on
#' any observed mode not present in the model's predictions. Returns the canonical
#' `(idx, obs)` frame.
#'
#' @param valid_modes Integer vector of the mode indices the model carries.
#' @param mode Integer vector of observed mode indices.
#' @param lrmsd_obs Numeric vector of observed log divergences, same length.
#' @return A tibble `(idx, obs)`.
#' @noRd
resolve_mode_obs <- function(valid_modes, mode, lrmsd_obs) {
  check_obs_vectors(mode, lrmsd_obs, "mode")
  unknown <- setdiff(mode, valid_modes)
  if (length(unknown) > 0) {
    stop("observed data has mode index(es) not present in the model: ",
         paste(unknown, collapse = ", "))
  }
  tibble::tibble(idx = mode, obs = lrmsd_obs)
}


#' Maximum-likelihood point fit of the MSA model to a site-axis profile
#'
#' Maximises the profiled Gaussian log-likelihood (`loglik_lrmsd_msa()`, evaluated on
#' `spm$dr2mat_site`)
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
#' @param spm A single-point-mutation `spm` object from [generate_spm()] (its `site_map` keys the fit to PDB residues).
#' @param pdb_site Integer vector of PDB residue numbers identifying the observations.
#'   Observations are a `(pdb_site, lrmsd_obs)` vector pair rather than a data frame, so
#'   the columns of your own table can be named anything. May cover a subset of the
#'   model's sites; every value must exist in the model.
#' @param lrmsd_obs Numeric vector of observed log structural divergences (the fit
#'   target), the same length as `pdb_site`.
#' @param a1_range Length-2 `[min, max]` box bound for `a1`.
#' @param log2_a2_plus1_range Length-2 `[min, max]` box bound for `log2(a2 + 1)`.
#' @param init Optional length-2 numeric start `c(a1, log2(a2+1))`. When `NULL`
#'   (default), a deterministic coarse grid-max of the likelihood over the box is
#'   used as the start (robust against a bad local start; cheap).
#' @param grid_n Number of points per axis for the default grid-max start (ignored
#'   when `init` is supplied).
#' @return A list. The top level is the **point estimate**; goodness of fit is the
#'   `gof` tibble (computed here, so there is no accessor to call):
#'   \describe{
#'     \item{a1, a2}{Point estimate of the stability (`a1`) and activity (`a2`)
#'       selection strengths, on the natural scale (the paper's `aS`, `aA`).}
#'     \item{logLik}{Profiled Gaussian log-likelihood at the optimum.}
#'     \item{cov}{2×2 covariance matrix on the `(a1, log2(a2+1))` scale.}
#'     \item{se_a1, se_a2}{Standard errors; `se_a2` via the delta method.}
#'     \item{convergence}{`optim` convergence code (0 = success).}
#'     \item{gof}{One-row tibble of fit quality: `D2` (deviance explained vs the flat
#'       null, `1 - deviance/null_deviance`; at most `1` but **unbounded below** --
#'       negative means worse than a constant profile, a real signal, not an error),
#'       `AIC`, `BIC`, and the primitives `logLik`, `deviance`, `null_deviance`,
#'       `nobs` (matched observations scored), `k` (`3`: `a1`, `a2`, profiled `sigma`),
#'       `sigma_hat`. `AIC`/`BIC` are per-fit numbers, meaningful only compared against
#'       another model fit at its own maximum.}
#'     \item{call}{The matched call that produced this fit -- the fit itself is
#'       axis-free, so this is what records which fitter made it.}
#'   }
#' @seealso [predict_profiles()] (propagate the fit to a banded profile),
#'   [fit_lrmsd_msa_mode()] (the mode counterpart),
#'   `loglik_lrmsd_msa()` (the objective).
#' @family fitting
#' @export
#' @examples
#' \dontrun{
#' ex   <- function(f) system.file("extdata", f, package = "msamodel")
#' wt   <- set_enm(bio3d::read.pdb(ex("1znb_A.pdb")), node = "ca",
#'                       model = "ming_wall", d_max = 10.5, frustrated = FALSE)
#' act  <- readr::read_csv(ex("znb_active_site.csv"))
#' spm  <- generate_spm(wt, pdb_site_active = act$pdb_site, seed = 1024)
#'
#' obs <- readr::read_csv(ex("znb_lrmsd_obs_site.csv"))
#' ml  <- fit_lrmsd_msa_site(spm, obs$pdb_site, obs$lrmsd_obs)
#' c(a1 = ml$a1, a2 = ml$a2)
#' }
fit_lrmsd_msa_site <- function(spm,
                       pdb_site,
                       lrmsd_obs,
                       a1_range = c(0, 10),
                       log2_a2_plus1_range = c(0, 13),
                       init = NULL,
                       grid_n = 25) {
  # The only axis-aware step: map pdb_site -> the internal site index (and raise the
  # site-axis unknown-key error), ONCE, before optimising. Everything after is axis-blind.
  fit_lrmsd_msa(spm$dr2mat_site, spm$energy_data,
                obs = resolve_site_obs(spm, pdb_site, lrmsd_obs),
                a1_range, log2_a2_plus1_range, init, grid_n,
                call = match.call())
}

#' Maximum-likelihood point fit of the MSA model to a mode-axis profile
#'
#' Mode counterpart of [fit_lrmsd_msa_site()]: maximises the profiled Gaussian
#' log-likelihood over `(a1, a2)` by numerical optimisation, returning a point estimate
#' plus an asymptotic covariance from the Hessian at the optimum. Identical machinery to
#' the site fit -- the same axis-blind objective `loglik_lrmsd_msa()`, evaluated on
#' `spm$dr2mat_mode` instead of `spm$dr2mat_site`; the response index is the mode (no `site_map`
#' / `pdb_site`).
#'
#' The optimiser works in the same coordinates as the site fit: `a1` and
#' `log2(a2 + 1)` (so `a2 = 2^(log2(a2+1)) - 1 >= 0`), on the box `a1_range` ×
#' `log2_a2_plus1_range` — the coordinates in which the prior is uniform. The
#' returned covariance `cov` is on the `(a1, log2(a2+1))` scale; the standard error
#' of `a2` is obtained by the delta method (`da2/d(log2(a2+1)) = 2^(log2(a2+1)) * ln 2`).
#'
#' @param spm A single-point-mutation `spm` object from [generate_spm()] (energy_data +
#'   the `dr2mat_mode` matrix; no `site_map`).
#' @param mode Integer vector of mode indices identifying the observations. Observations
#'   are a `(mode, lrmsd_obs)` vector pair rather than a data frame, so the columns of
#'   your own table can be named anything. Every value must exist in the model.
#' @param lrmsd_obs Numeric vector of observed log structural divergences (the fit
#'   target), the same length as `mode`.
#' @param a1_range Length-2 `[min, max]` box bound for `a1`.
#' @param log2_a2_plus1_range Length-2 `[min, max]` box bound for `log2(a2 + 1)`.
#' @param init Optional length-2 numeric start `c(a1, log2(a2+1))`. When `NULL`
#'   (default), a deterministic coarse grid-max of the likelihood over the box is
#'   used as the start (robust against a bad local start; cheap).
#' @param grid_n Number of points per axis for the default grid-max start (ignored
#'   when `init` is supplied).
#' @inherit fit_lrmsd_msa_site return
#' @seealso [fit_lrmsd_msa_site()] (the site counterpart),
#'   [predict_profiles()] (propagate the fit to a banded profile),
#'   `loglik_lrmsd_msa()` (the objective).
#' @family fitting
#' @export
#' @examples
#' \dontrun{
#' ex   <- function(f) system.file("extdata", f, package = "msamodel")
#' wt   <- set_enm(bio3d::read.pdb(ex("1znb_A.pdb")), node = "ca",
#'                       model = "ming_wall", d_max = 10.5, frustrated = FALSE)
#' act  <- readr::read_csv(ex("znb_active_site.csv"))
#' spm  <- generate_spm(wt, pdb_site_active = act$pdb_site, seed = 1024)
#'
#' # Synthetic observations: no empirical per-mode profile exists (see the file name).
#' obs <- readr::read_csv(ex("znb_lrmsd_obs_mode_syn.csv"))
#' ml  <- fit_lrmsd_msa_mode(spm, obs$mode, obs$lrmsd_obs)
#' c(a1 = ml$a1, a2 = ml$a2)
#' }
fit_lrmsd_msa_mode <- function(spm,
                       mode,
                       lrmsd_obs,
                       a1_range = c(0, 10),
                       log2_a2_plus1_range = c(0, 13),
                       init = NULL,
                       grid_n = 25) {
  # The only axis-aware step: the mode number IS the internal index (no site_map), so this
  # just range-checks it, ONCE, before optimising. Everything after is axis-blind.
  fit_lrmsd_msa(spm$dr2mat_mode, spm$energy_data,
                obs = resolve_mode_obs(seq_len(ncol(spm$dr2mat_mode)), mode, lrmsd_obs),
                a1_range, log2_a2_plus1_range, init, grid_n,
                call = match.call())
}
