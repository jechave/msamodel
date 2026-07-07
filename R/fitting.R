# MSA model fitting — maximum-likelihood point estimation
# Point-estimate counterpart of the AGQ posterior fit (fit_lrmsd_i_msa_agq, below).
# Shares the same profiled Gaussian log-likelihood (calculate_loglik_lrmsd_i_msa) and
# the same (a1, log2(a2+1)) coordinates / box bounds, so the two arms are directly
# comparable.

#' Maximum-likelihood point fit of the lrmsd_i MSA model
#'
#' Point-estimation counterpart of [fit_lrmsd_i_msa_agq()]: maximises the same
#' profiled Gaussian log-likelihood (`calculate_loglik_lrmsd_i_msa()`) over `(a1, a2)`
#' by numerical
#' optimisation, returning a point estimate plus an asymptotic covariance from the
#' Hessian at the optimum. Where AGQ returns a full posterior over the selection
#' strengths, this returns a single point plus asymptotic errors; intended for large
#' proteins and path simulations. This is **not** a Bayesian fit and returns no
#' posterior sample.
#'
#' The optimiser works in the same coordinates as the AGQ fit: `a1` and
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
#'     \item{sigma_hat}{Profiled noise scale `sqrt(mean(residuals^2))` at the optimum.}
#'     \item{cov}{2×2 covariance matrix on the `(a1, log2(a2+1))` scale.}
#'     \item{se_a1, se_a2}{Standard errors; `se_a2` via the delta method.}
#'     \item{convergence}{`optim` convergence code (0 = success).}
#'   }
#' @seealso [fit_lrmsd_i_msa_agq()] (the Bayesian/posterior counterpart),
#'   `calculate_loglik_lrmsd_i_msa()` (the shared objective).
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
  # Validate box bounds (same contract as fit_lrmsd_i_msa_agq) -- fail loud.
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

  list(
    a1          = unname(a1_hat),
    a2          = unname(a2_hat),
    logLik      = -opt$value,
    sigma_hat   = sigma_hat,
    cov         = cov,
    se_a1       = unname(se_a1),
    se_a2       = unname(se_a2),
    convergence = opt$convergence
  )
}

#' Maximum-likelihood point fit of the lrmsd_n MSA model (mode form)
#'
#' Mode counterpart of [fit_lrmsd_i_msa_ml()]: maximises the mode-form profiled
#' Gaussian log-likelihood (`calculate_loglik_lrmsd_n_msa()`) over `(a1, a2)` by
#' numerical optimisation, returning a point estimate plus an asymptotic covariance
#' from the Hessian at the optimum. Identical machinery to the site fit; the response
#' index is the mode `n` (no `site_map` / `pdb_site`). There is no Bayesian mode
#' counterpart.
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
#'     \item{sigma_hat}{Profiled noise scale `sqrt(mean(residuals^2))` at the optimum.}
#'     \item{cov}{2×2 covariance matrix on the `(a1, log2(a2+1))` scale.}
#'     \item{se_a1, se_a2}{Standard errors; `se_a2` via the delta method.}
#'     \item{convergence}{`optim` convergence code (0 = success).}
#'   }
#' @seealso [fit_lrmsd_i_msa_ml()] (the site counterpart),
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

  list(
    a1          = unname(a1_hat),
    a2          = unname(a2_hat),
    logLik      = -opt$value,
    sigma_hat   = sigma_hat,
    cov         = cov,
    se_a1       = unname(se_a1),
    se_a2       = unname(se_a2),
    convergence = opt$convergence
  )
}

# Gauss-Hermite nodes and weights (physicists' convention, weight e^{-x^2}) via
# Golub-Welsch: eigen-decomposition of the symmetric Jacobi matrix. Returns nodes
# `x` and weights `w` with sum(w) = sqrt(pi). Small n; built directly.
#' @noRd
gauss_hermite <- function(n) {
  if (n < 1L) stop("n must be a positive integer")
  if (n == 1L) return(list(x = 0, w = sqrt(pi)))
  i <- 1:(n - 1)
  a <- sqrt(i / 2)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1)] <- a
  J[cbind(i + 1, i)] <- a
  e <- eigen(J, symmetric = TRUE)
  x <- e$values
  w <- sqrt(pi) * e$vectors[1, ]^2
  o <- order(x)
  list(x = x[o], w = w[o])
}


#' Posterior fit of the lrmsd_i MSA model by adaptive Gauss-Hermite quadrature
#'
#' The Bayesian arm of the MSA fit: it approximates the posterior of the selection
#' strengths `(a1, a2)` by adaptive Gauss-Hermite quadrature, under a uniform prior
#' in `(a1, log2(a2 + 1))` and the profiled Gaussian likelihood
#' (`calculate_loglik_lrmsd_i_msa()`). The quadrature is *referenced* to the Laplace
#' (Gaussian) approximation that [fit_lrmsd_i_msa_ml()] already computes: nodes are
#' placed at that Gaussian's optimal locations, so a handful of likelihood
#' evaluations reproduce the posterior. For the bundled `znb_profile` the posterior
#' is near-Gaussian in `(a1, log2(a2 + 1))`, and even a 5x5 grid (25 evaluations)
#' matches a dense reference posterior essentially exactly. The fit is fully
#' deterministic -- no seed, burn-in, or autocorrelation.
#'
#' The quadrature is built on the unconstrained scale `t = log2(a2 + 1)` (where the
#' prior is flat and the posterior is near-Gaussian) and results are reported on the
#' natural scale `a2`, transformed at the boundary -- the standard
#' constrained/unconstrained convention. Node weights are stored as **unnormalized
#' log masses** (`log_weight`); a probability mass is reparameterization-invariant,
#' so the node table is reported entirely on the natural `(a1, a2)` scale and any
#' posterior expectation is the plain weighted sum
#' `sum(exp(log_weight) * g(a1, a2)) / sum(exp(log_weight))`. The only quantity kept
#' on the `t` scale is the Laplace covariance (see `laplace$cov`).
#'
#' @param spm_pp Preprocessed data from [preprocess_spm()] (must include `site_map`).
#' @param observed_data Tibble with columns `pdb_site` and `lrmsd_i_obs` (the fit
#'   target), as documented for `calculate_loglik_lrmsd_i_msa()`.
#' @param n_nodes Number of Gauss-Hermite nodes per axis (the quadrature uses an
#'   `n_nodes` x `n_nodes` tensor grid). Default 7 (49 evaluations); the posterior
#'   moments are already accurate at 5, the extra nodes give the credible band margin.
#' @param a1_range,log2_a2_plus1_range,init,grid_n Passed to [fit_lrmsd_i_msa_ml()]
#'   to obtain the Laplace reference (box bounds, optional start, grid-max start size).
#' @return A list of class `"msa_agq"` with the posterior summary and the quadrature:
#'   \describe{
#'     \item{a1, a2}{Posterior means of the stability (`a1`) and activity (`a2`)
#'       selection strengths, natural scale (the paper's `aS`, `aA`).}
#'     \item{sd_a1, sd_a2}{Posterior standard deviations, natural scale.}
#'     \item{ci_a1, ci_a2}{Length-2 95% credible intervals, natural scale.}
#'     \item{nodes}{Tibble `a1, a2, log_weight` -- the quadrature nodes on the natural
#'       scale and their unnormalized log posterior masses. The propagation primitive.}
#'     \item{n_nodes}{Nodes per axis.}
#'     \item{laplace}{`list(mu, cov)`, the Gaussian reference used to place nodes, on
#'       the `(a1, log2(a2 + 1))` scale (`cov` carries matching `dimnames`).}
#'     \item{log_evidence}{Log marginal likelihood from the quadrature.}
#'   }
#' @seealso [fit_lrmsd_i_msa_ml()] (the Laplace reference / point-estimate arm),
#'   [predict_lrmsd_i_msa_agq()] (propagate the posterior to a banded profile).
#' @family fitting
#' @export
#' @examples
#' \dontrun{
#' pp  <- preprocess_spm(znb_spm)
#' agq <- fit_lrmsd_i_msa_agq(pp, znb_profile)
#' c(a1 = agq$a1, a2 = agq$a2)
#' }
fit_lrmsd_i_msa_agq <- function(spm_pp,
                                observed_data,
                                n_nodes = 7,
                                a1_range = c(0, 10),
                                log2_a2_plus1_range = c(0, 13),
                                init = NULL,
                                grid_n = 25) {
  if (length(n_nodes) != 1 || n_nodes < 1 || n_nodes != round(n_nodes)) {
    stop("n_nodes must be a single positive integer")
  }

  # Laplace reference: reuse the ML fitter's optimum + covariance. mu and cov are on
  # the t = (a1, log2(a2 + 1)) scale (cov already documented as such by the ML fit).
  ml <- fit_lrmsd_i_msa_ml(spm_pp, observed_data,
                           a1_range = a1_range,
                           log2_a2_plus1_range = log2_a2_plus1_range,
                           init = init, grid_n = grid_n)
  mu <- c(a1 = ml$a1, log2_a2_plus1 = log2(ml$a2 + 1))
  S  <- ml$cov
  dimnames(S) <- list(c("a1", "log2_a2_plus1"), c("a1", "log2_a2_plus1"))
  L  <- chol(S)  # S = t(L) %*% L

  # Tensor Gauss-Hermite grid, mapped to N(mu, S): y = mu + sqrt(2) * t(L) %*% z.
  gh <- gauss_hermite(n_nodes)
  Z  <- as.matrix(expand.grid(z1 = gh$x, z2 = gh$x))
  Wq <- as.vector(t(outer(gh$w, gh$w)))                 # tensor GH weights
  Y  <- t(mu + sqrt(2) * t(L) %*% t(Z))                 # nodes in (a1, log2_a2_plus1)
  a1_nodes <- Y[, 1]
  a2_nodes <- 2^Y[, 2] - 1

  # Log-likelihood at each node.
  ll <- vapply(seq_len(nrow(Y)), function(k)
    calculate_loglik_lrmsd_i_msa(spm_pp, observed_data, a1_nodes[k], a2_nodes[k]),
    numeric(1))

  # Change of measure: the GH rule integrates against e^{-|z|^2}, so divide that
  # reference out (the + rowSums(Z^2) term) to recover the flat-measure posterior
  # mass. logW are unnormalized log masses; normalize via the max-shift.
  z2   <- rowSums(Z^2)
  logW <- log(Wq) + z2 + (ll - max(ll))
  log_norm   <- max(logW) + log(sum(exp(logW - max(logW))))
  log_weight <- logW - log_norm                         # sum(exp(log_weight)) = 1
  w <- exp(log_weight)

  # Posterior moments on the natural scale (mass is reparameterization-invariant).
  m_a1 <- sum(a1_nodes * w); m_a2 <- sum(a2_nodes * w)
  sd_a1 <- sqrt(sum((a1_nodes - m_a1)^2 * w))
  sd_a2 <- sqrt(sum((a2_nodes - m_a2)^2 * w))
  ci_a1 <- weighted_quantile(a1_nodes, w, c(0.025, 0.975))
  ci_a2 <- weighted_quantile(a2_nodes, w, c(0.025, 0.975))

  # Log marginal likelihood (evidence): integral of exp(ll) over the flat measure.
  # det Jacobian of the affine map is sqrt(2)^d * det(L) = sqrt(2)^d * sqrt(det(S)).
  d <- 2
  log_det_jac  <- d * log(sqrt(2)) + 0.5 * log(det(S))
  log_evidence <- max(logW) + log(sum(exp(logW - max(logW)))) + log_det_jac

  structure(list(
    a1 = m_a1, a2 = m_a2,
    sd_a1 = sd_a1, sd_a2 = sd_a2,
    ci_a1 = unname(ci_a1), ci_a2 = unname(ci_a2),
    nodes = tibble(a1 = a1_nodes, a2 = a2_nodes, log_weight = log_weight),
    n_nodes = as.integer(n_nodes),
    laplace = list(mu = mu, cov = S),
    log_evidence = log_evidence
  ), class = "msa_agq")
}

# ---- fitting objective: profiled Gaussian log-likelihood (internal) --------------
# The criterion both fitters optimise/integrate. Internal (@noRd): users never call
# the raw likelihood; a future goodness-of-fit surface (AIC etc.) will be separate
# public functions. Kept a distinct function so the criterion stays swappable (a
# future stochastic-likelihood tree plugs in here without touching the fitters).

# Log-likelihood of an observed per-site divergence profile at (a1, a2): mean-centre
# model prediction and observations, compare under a Gaussian noise model whose scale
# is profiled out (sigma = sqrt(mean(r^2))). Maximised by fit_lrmsd_i_msa_ml,
# integrated by fit_lrmsd_i_msa_agq. spm_pp: preprocess_spm() output; observed_data:
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
  # quantities — logLik, SEs, posterior width — were slightly off.)
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
