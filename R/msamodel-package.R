#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom dplyr across all_of any_of arrange bind_rows everything filter group_by group_split inner_join left_join mutate rename select starts_with summarise transmute where
#' @importFrom tidyr complete fill pivot_longer pivot_wider
#' @importFrom tibble as_tibble tibble
#' @importFrom purrr discard map_dfr possibly
#' @importFrom rlang sym !!
#' @importFrom stats cor dnorm dunif fitted loess median optim optimHess predict qnorm quantile residuals rnorm runif sd var
#' @importFrom penm set_enm get_site get_pdb_site get_mutant_site ddg_dv ddg_tds ddgact_dv ddgact_tds
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
#   generate_spm_data : m, j, ddg_dv_jm, ddg_tds_jm, ddgact_dv_jm, ddgact_tds_jm,
#                       ddg, ddgact
#   resolve_site_obs  : site, obs
#
# (ddg / ddgact are assigned by mutate() and then named bare in select(); the
# _dv_ / _tds_ names are the raw scan columns the reshape reads.)
