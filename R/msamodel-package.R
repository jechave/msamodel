#' @details
#' Three things recur throughout the package.
#'
#' **The selection strengths, `a1` and `a2`.** `a1` sets the strength of selection on
#' stability, `a2` the strength of selection on activity, that is on proximity to the
#' active site. Both are non-negative, and `0` switches that pressure off. They are in
#' mol/kcal, the inverse of the unit of the energy changes they multiply (see
#' [generate_spm()]). They are either set to chosen values, in [calculate_profiles()]
#' and [calculate_decomposition()], or estimated from an observed profile by
#' [fit_lrmsd_msa_site()] and [fit_lrmsd_msa_mode()] and then passed to
#' [predict_profiles()] and [predict_decomposition()].
#'
#' **The two representations, `site` and `mode`.** A mutation displaces the structure,
#' and that displacement is written in two bases: per residue, and per normal mode of
#' the elastic network. Both hold the same divergence. Every profile function returns
#' both at once, as a list with a \code{$site} tibble, keyed by `site` and `pdb_site`,
#' and a \code{$mode} tibble, keyed by `mode`.
#'
#' **The two metrics, `lrmsd` and `nlrmsd`.** `lrmsd` is the divergence profile as the
#' model predicts it, keeping its overall level. `nlrmsd` is the same profile with its
#' mean subtracted off. An observed divergence profile fixes shape but not level, so
#' `nlrmsd` is the one to compare against data.
#'
#' @section The workflow:
#'
#' \describe{
#'   \item{1. Build the elastic network}{[penm::set_enm()], which belongs to the penm
#'     package, not to msamodel.}
#'   \item{2. Run the mutation scan}{[generate_spm()] gives every site each other
#'     letter of the sequence alphabet and records how the structure responds.
#'     Returns the `spm` object every later function reads.}
#'   \item{3a. Evaluate at chosen `(a1, a2)`}{[calculate_profiles()] for the
#'     divergence profiles, [calculate_decomposition()] for the four nested models
#'     and the mutation/stability/activity split.}
#'   \item{3b. Or estimate `(a1, a2)` from data}{[fit_lrmsd_msa_site()] fits an
#'     observed site profile, [fit_lrmsd_msa_mode()] an observed mode profile. Both
#'     return the estimate, its standard errors, and a goodness-of-fit summary in
#'     \code{$gof}.}
#'   \item{4. Propagate the fit}{[predict_profiles()] and [predict_decomposition()]
#'     give the same quantities as step 3a, with delta-method standard errors
#'     carried through from the fit.}
#' }
#'
#' [pfix_msa()] exposes the model's elementary quantity, a single mutant's fixation
#' probability.
#'
#' @section Getting started:
#'
#' `vignette("msamodel")` runs the whole workflow, scan to decomposition. See also
#' `vignette("msamodel-explore")` and `vignette("msamodel-fit")`.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom dplyr across all_of any_of arrange bind_rows everything filter group_by group_split inner_join left_join mutate rename select starts_with summarise transmute where
#' @importFrom tidyr complete fill pivot_longer pivot_wider
#' @importFrom tibble as_tibble tibble
#' @importFrom purrr discard map_dfr possibly
#' @importFrom rlang sym !! .data
#' @importFrom stats cor dnorm dunif fitted loess median optim optimHess predict qnorm quantile residuals rnorm runif sd var
#' @importFrom magrittr %>%
## usethis namespace: end
NULL

# Column names inside dplyr verbs are written `.data$col` (and as strings in
# select()), not bare. A bare name is looked up as a column first and then in the
# surrounding environment, so a missing column can silently pick up an unrelated
# variable; `.data$col` can only mean the column and errors when it is absent.
# R CMD check used to report the bare ones as "no visible binding for global
# variable" (in generate_spm and resolve_site_obs). They were rewritten rather than
# declared in utils::globalVariables(): an earlier such list had grown to 47 names,
# only 11 of them real, and it muted the check for those names package-wide.
