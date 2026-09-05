#' @details
#' The MSA model predicts how far each part of a protein structure drifts under
#' neutral evolution, given the structure itself and a list of its active-site
#' residues. Divergence is read on two axes -- per residue and per normal mode --
#' and every profile function returns both at once, as a list with a `$site` and a
#' `$mode` tibble.
#'
#' Two selection strengths govern the model: `a1` (selection on stability) and `a2`
#' (selection on activity, i.e. proximity to the active site). Both are
#' non-negative, and `0` switches that pressure off. You either choose them, or
#' estimate them from an observed divergence profile.
#'
#' @section The workflow:
#'
#' The scan is the one slow step; everything downstream of it is fast.
#'
#' \describe{
#'   \item{1. Build the elastic network}{[penm::set_enm()] -- msamodel does not do
#'     this; it belongs to the penm package, which you attach yourself.}
#'   \item{2. Run the mutation scan}{[generate_spm()] -- mutates every site many
#'     times and records how the structure responds. Returns the `spm` object every
#'     later function reads.}
#'   \item{3a. Evaluate at chosen `(a1, a2)`}{[calculate_profiles()] for the
#'     divergence profile, [calculate_decomposition()] for the four nested model
#'     variants and the mutation/stability/activity split.}
#'   \item{3b. Or fit `(a1, a2)` to your data}{[fit_lrmsd_msa_site()] on the residue
#'     axis, [fit_lrmsd_msa_mode()] on the mode axis. Both return the estimate, its
#'     standard errors, and a goodness-of-fit summary in `$gof`.}
#'   \item{4. Propagate the fit}{[predict_profiles()] and [predict_decomposition()]
#'     give the same quantities as step 3a, with delta-method standard errors
#'     carried through from the fit.}
#' }
#'
#' [pfix_msa()] exposes the model's elementary quantity, a single mutant's fixation
#' probability.
#'
#' @section Two metrics:
#'
#' Every profile function takes a `metric` argument. `"lrmsd"` is the profile as the
#' model predicts it, keeping its overall level; `"nlrmsd"` is the same profile with
#' its mean subtracted off -- its *shape*. An observed divergence profile pins down
#' shape but not level, so `"nlrmsd"` is the one to compare against data.
#'
#' @section Getting started:
#'
#' `vignette("msamodel")` is a short end-to-end tour of both axes. See also
#' `vignette("msamodel-explore")` and `vignette("msamodel-fit")`.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom dplyr across all_of any_of arrange bind_rows everything filter group_by group_split inner_join left_join mutate rename select starts_with summarise transmute where
#' @importFrom tidyr complete fill pivot_longer pivot_wider
#' @importFrom tibble as_tibble tibble
#' @importFrom purrr discard map_dfr possibly
#' @importFrom rlang sym !!
#' @importFrom stats cor dnorm dunif fitted loess median optim optimHess predict qnorm quantile residuals rnorm runif sd var
#' @importFrom magrittr %>%
## usethis namespace: end
NULL

# NOTE: R CMD check reports "no visible binding for global variable" for column
# names referenced bare inside dplyr/tidyr verbs (NSE). It is a NOTE, not an error
# or a warning -- the code is correct. A utils::globalVariables() call used to
# suppress it; that list had grown to 47 names of which only 11 corresponded to a
# real complaint, so it was asserting things about the package that were not true.
# The note is left visible rather than silenced.
#
# The names R actually cannot resolve, in 2 functions:
#   generate_spm     : m, j, ddg_dv_jm, ddg_tds_jm, ddgact_dv_jm, ddgact_tds_jm,
#                      ddg, ddgact
#   resolve_site_obs : site, obs
#
# (ddg / ddgact are assigned by mutate() and then named bare in select(); the
# _dv_ / _tds_ names are the raw scan columns the reshape reads.)
