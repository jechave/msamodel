#' MSA model fitting — maximum-likelihood point estimation
#' Point-estimate counterpart of the MCMC fit (R/msa_bayesian_analysis.R). Shares
#' the same profiled Gaussian log-likelihood (calculate_loglik_msa) and the same
#' (a1, log2(a2+1)) coordinates / box bounds, so the two arms are directly
#' comparable.

#' Maximum-likelihood point estimate of the MSA selection parameters
#'
#' Point-estimation counterpart of [run_mcmc_msa()]: maximises the same profiled
#' Gaussian log-likelihood ([calculate_loglik_msa()]) over `(a1, a2)` by numerical
#' optimisation, returning a point estimate plus an asymptotic covariance from the
#' Hessian at the optimum. Much faster than MCMC; intended for large proteins and
#' path simulations. This is **not** a Bayesian fit and returns no posterior sample.
#'
#' The optimiser works in the same coordinates as the MCMC: `a1` and
#' `b = log2(a2 + 1)` (so `a2 = 2^b - 1 >= 0`), on the box `a1_range` ×
#' `log2_a2_plus1_range`. The returned covariance `cov` is on the `(a1, b)` scale;
#' the standard error of `a2` is obtained by the delta method
#' (`da2/db = 2^b * ln 2`).
#'
#' @param spm_pp Preprocessed data from [preprocess_spm()] (must include `site_map`).
#' @param observed_data Tibble with columns `pdb_site` and `lrmsd_obs` (the fit
#'   target), as documented for [calculate_loglik_msa()].
#' @param a1_range Length-2 `[min, max]` box bound for `a1`.
#' @param log2_a2_plus1_range Length-2 `[min, max]` box bound for `log2(a2 + 1)`.
#' @param init Optional length-2 numeric start `c(a1, log2(a2+1))`. When `NULL`
#'   (default), a deterministic coarse grid-max of the likelihood over the box is
#'   used as the start (robust against a bad local start; cheap).
#' @param grid_n Number of points per axis for the default grid-max start (ignored
#'   when `init` is supplied).
#' @return A list with the point estimate and asymptotic uncertainty:
#'   \describe{
#'     \item{a1, a2}{Point estimate (`a2` on the natural scale).}
#'     \item{logLik}{Profiled Gaussian log-likelihood at the optimum.}
#'     \item{sigma_hat}{Profiled noise scale `sqrt(mean(residuals^2))` at the optimum.}
#'     \item{cov}{2×2 covariance matrix on the `(a1, log2(a2+1))` scale.}
#'     \item{se_a1, se_a2}{Standard errors; `se_a2` via the delta method.}
#'     \item{convergence}{`optim` convergence code (0 = success).}
#'     \item{par_fit}{The fitted coordinates `c(a1, b)` with `b = log2(a2+1)`.}
#'   }
#' @seealso [run_mcmc_msa()] (the Bayesian counterpart), [calculate_loglik_msa()]
#'   (the shared objective).
#' @family fitting
#' @export
#' @examples
#' \dontrun{
#' pp <- preprocess_spm(znb_spm)
#' ml <- fit_msa_ml(pp, znb_profile)
#' c(a1 = ml$a1, a2 = ml$a2)
#' }
fit_msa_ml <- function(spm_pp,
                       observed_data,
                       a1_range = c(0, 10),
                       log2_a2_plus1_range = c(0, 13),
                       init = NULL,
                       grid_n = 25) {
  # Validate box bounds (same contract as run_mcmc_msa) -- fail loud.
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
  # calculate_loglik_msa already mean-centers both profiles, so this IS the
  # MCMC's likelihood (no separate centering here).
  nll <- function(theta) {
    -calculate_loglik_msa(spm_pp, observed_data,
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
  dr2_i_msa <- calculate_dr2_i_msa(spm_pp, a1_hat, a2_hat)
  pred <- dr2_i_msa %>%
    mutate(lrmsd_i_msa = log(sqrt(dr2_i))) %>%
    dplyr::select(i, lrmsd_i_msa)
  obs <- observed_data %>%
    dplyr::select(pdb_site, lrmsd_obs) %>%
    inner_join(spm_pp$site_map, by = "pdb_site") %>%
    dplyr::select(i, lrmsd_obs)
  cmp <- obs %>%
    inner_join(pred, by = "i") %>%
    mutate(
      nlrmsd_obs   = lrmsd_obs - mean(lrmsd_obs),
      nlrmsd_i_msa = lrmsd_i_msa - mean(lrmsd_i_msa)
    )
  residuals <- cmp$nlrmsd_obs - cmp$nlrmsd_i_msa
  sigma_hat <- sqrt(mean(residuals^2))

  list(
    a1          = unname(a1_hat),
    a2          = unname(a2_hat),
    logLik      = -opt$value,
    sigma_hat   = sigma_hat,
    cov         = cov,
    se_a1       = unname(se_a1),
    se_a2       = unname(se_a2),
    convergence = opt$convergence,
    par_fit     = c(a1 = unname(a1_hat), b = unname(b_hat))
  )
}
