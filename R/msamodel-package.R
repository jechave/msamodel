#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom dplyr across all_of any_of arrange bind_rows everything filter group_by group_split inner_join left_join mutate rename select starts_with summarise transmute where
#' @importFrom tidyr complete fill pivot_longer pivot_wider
#' @importFrom tibble as_tibble tibble
#' @importFrom purrr discard map_dfr possibly
#' @importFrom rlang sym !!
#' @importFrom stats cor dnorm dunif fitted loess median optim optimHess predict qnorm quantile residuals rnorm runif sd var
#' @importFrom penm set_enm get_site get_pdb_site get_mutant_site get_dactive get_cn get_msf_site ddg_dv ddg_tds ddgact_dv ddgact_tds
#' @importFrom magrittr %>%
## usethis namespace: end
NULL

# Quiet R CMD check's "no visible binding for global variable" NOTEs from
# tidyverse non-standard evaluation (column names referenced bare inside
# dplyr/tidyr verbs). These are data-column / grouping names, not real globals.
utils::globalVariables(c(
  "a1", "a2", "dactive", "ddg_dv_jm", "ddg_jm", "ddg_tds_jm", "ddgact_dv_jm",
  "ddgact_jm", "ddgact_tds_jm", "dr2_i", "dr2_ijm", "dr2_n", "dr2_njm",
  "i", "j", "lower", "lrmsd", "lrmsd_i_ma", "lrmsd_i_mm", "lrmsd_i_ms", "lrmsd_i_msa",
  "lrmsd_n_ma", "lrmsd_n_mm", "lrmsd_n_ms", "lrmsd_n_msa",
  "lrmsd_i_obs", "lrmsd_n_obs", "lrmsd_n_true", "m", "mode", "model", "n",
  "lrmsd_i_msa_mean", "lrmsd_i_msa_lower", "lrmsd_i_msa_upper", "log_weight",
  "nlrmsd_i_msa", "nlrmsd_i_obs", "nlrmsd_n_msa", "nlrmsd_n_obs",
  "parameter", "pdb_site",
  "phi_act", "phi_mut", "phi_stab", "sample_id",
  "upper", "value", "variable"
))
