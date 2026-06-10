# Precompute heavy results for vignettes/msamodel-intro.Rmd
#
# The vignette demonstrates the full msamodel pipeline (SPM scan, profiles,
# parameter sweeps, Bayesian MCMC fit, decomposition). Running that at build
# time would make `R CMD check` slow and would require `penm` as a build-time
# dependency. Instead we run it once here and cache the (small) summary objects
# the vignette needs; the .Rmd loads `vignette_cache.rds` and only draws plots.
#
# Regenerate intentionally:
#   Rscript vignettes/precompute.R
#
# This is excluded from the package build via .Rbuildignore.

suppressMessages(devtools::load_all("."))
suppressMessages({
  library(dplyr)
  library(tidyr)
})

set.seed(1024)

# ---- Data ------------------------------------------------------------------
# The shipped znb_* datasets are the cached result of the "Setup from PDB"
# section (load_protein -> get_active_site -> setup_enm -> generate_spm_data).
data("znb_wt", package = "msamodel")
data("znb_spm", package = "msamodel")
data("znb_profile", package = "msamodel")
data("znb_dataset", package = "msamodel")

active <- get_active_site("1znb_A", znb_dataset)
pdb_site_active <- active$pdb_site_active

# Preprocessed SPM (energies + dr2 matrix + site_map) used everywhere below.
pp <- preprocess_spm(znb_spm)

# Per-site structural / flexibility properties (lrmsf, dactive). add_site_properties
# joins on the internal index i, so build a one-row-per-site frame keyed by i.
site_ids <- as.integer(colnames(pp$dr2mat))
site_props <- add_site_properties(tibble(i = site_ids), znb_wt, pdb_site_active) %>%
  mutate(lrmsf = log(sqrt(msf))) %>%
  left_join(pp$site_map, by = "i") %>%
  select(i, pdb_site, dactive, lrmsf, msf, cn, shell)

# Helper: raw lrmsd profile for given (a1, a2), with pdb_site + properties.
lrmsd_profile <- function(a1, a2) {
  calculate_dr2i_msa(pp, a1, a2) %>%
    mutate(lrmsd = log(sqrt(dr2_i))) %>%
    left_join(site_props, by = "i") %>%
    select(i, pdb_site, lrmsd, dactive, lrmsf)
}

# ---- 3. Profile at an illustrative (a1, a2) --------------------------------
a1_demo <- 2
a2_demo <- 500
profile_demo <- lrmsd_profile(a1_demo, a2_demo)

# ---- 4. Parameter sweeps (whole raw lrmsd) ---------------------------------
a1_values <- c(0, 1, 2, 4, 8)     # at a2 = 0
a2_values <- c(0, 31, 127, 511, 2047)  # at a1 = 0 (2^k - 1 spacing)

sweep_a1 <- purrr::map_dfr(a1_values, function(a1)
  lrmsd_profile(a1, 0) %>% mutate(a1 = a1))
sweep_a2 <- purrr::map_dfr(a2_values, function(a2)
  lrmsd_profile(0, a2) %>% mutate(a2 = a2))

# ---- 5. Bayesian fit -------------------------------------------------------
# Modest but real run (documented in the vignette as "for demonstration").
fit <- run_msa_bayesian_analysis(
  spm            = znb_spm,
  observed_data  = znb_profile,
  n_mcmc_iter    = 4000,
  n_burnin       = 1000
)

parameter_summary     <- fit$parameter_summary
prediction_summary    <- fit$prediction_summary
decomposition_summary <- fit$decomposition_summary

# Fitted-vs-observed comparison table, mean-centered (nlrmsd) — exactly what the
# likelihood compares. Posterior-mean MSA prediction per site, joined to observed
# by pdb_site, then both centered on the COMPARED subset.
fit_msa <- prediction_summary %>%
  filter(variable == "lrmsd_msa") %>%
  select(pdb_site, lrmsd_fit = mean)

comparison <- znb_profile %>%
  inner_join(fit_msa, by = "pdb_site") %>%
  left_join(select(site_props, pdb_site, dactive, lrmsf), by = "pdb_site") %>%
  mutate(
    nlrmsd_obs = lrmsd_obs - mean(lrmsd_obs),
    nlrmsd_fit = lrmsd_fit - mean(lrmsd_fit)
  ) %>%
  select(pdb_site, nlrmsd_obs, nlrmsd_fit, dactive, lrmsf)

# ---- Save cache ------------------------------------------------------------
vignette_cache <- list(
  meta = list(
    pdb_chain       = "1znb_A",
    pdb_site_active = pdb_site_active,
    a1_demo = a1_demo, a2_demo = a2_demo,
    a1_values = a1_values, a2_values = a2_values,
    n_mcmc_iter = 4000, n_burnin = 1000,
    n_sites = nrow(site_props), n_obs = nrow(znb_profile)
  ),
  site_props            = site_props,
  profile_demo          = profile_demo,
  sweep_a1              = sweep_a1,
  sweep_a2              = sweep_a2,
  parameter_summary     = parameter_summary,
  prediction_summary    = prediction_summary,
  decomposition_summary = decomposition_summary,
  comparison            = comparison
)

out <- file.path("vignettes", "vignette_cache.rds")
saveRDS(vignette_cache, out)
message("Wrote ", out)
