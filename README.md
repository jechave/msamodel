
<!-- README.md is generated from README.Rmd. Please edit that file and knit. -->

# msamodel

<!-- badges: start -->

[![R-CMD-check](https://github.com/jechave/msamodel/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jechave/msamodel/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The `msamodel` package implements the Mutation-Stability-Activity (MSA)
model of protein structural evolution. In the MSA model, mutations
perturb the structure and become fixed or are lost depending on their
effects on stability and activity, weighted by stability and activity
selection strengths, respectively. Given a single structure and its
active-site residues, the model predicts the resulting structural
divergence profiles, per residue and per normal mode. The two selection
strengths can be set to chosen values, to explore how each constraint
shapes the profiles, or estimated from observed profiles by maximum
likelihood. In addition, predicted profiles can be analysed in terms of
nested “what-if” scenarios in which one or both selection constraints
are dropped, and can be decomposed as a sum of mutation, stability, and
activity contributions.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("jechave/msamodel", dependencies = TRUE)
```

The elastic network model comes from
[penm](https://github.com/jechave/penm), installed alongside. Neither
package is on CRAN.

## Example

Mutate every site, fit the model to an observed profile, predict:

``` r
library(msamodel)
library(penm)

ex     <- function(f) system.file("extdata", f, package = "msamodel")
active <- read.csv(ex("1d6o_A_active_site.csv"))
obs    <- read.csv(ex("1d6o_A_lrmsd_obs_site.csv"))

wt  <- penm::set_enm(bio3d::read.pdb(ex("1d6o_A.pdb")), node = "ca",
                     model = "ming_wall", d_max = 10.5, frustrated = FALSE)
spm <- generate_spm(wt, n_mutations = 10,
                    pdb_site_active = active$pdb_site, ensemble = 1L)

fit  <- fit_lrmsd_msa_site(spm, obs$pdb_site, obs$lrmsd_obs)
pred <- predict_profiles(fit, spm, metric = "nlrmsd")$site
```

<img src="man/figures/README-profile-1.png" width="100%" />

Observed divergence (points), fitted model with 95% band, active-site
residues dashed. The two fitted selection strengths are 0.23 on
stability and 79 on activity, accounting for 60% of the variance in the
observed profile.

The same fit decomposes the profile into its three contributions:

``` r
dec <- predict_decomposition(fit, spm, metric = "nlrmsd")$site
```

<img src="man/figures/README-decomposition-1.png" width="100%" />

The three contributions add up to the profile above. Their balance
varies from family to family.

## Interface

`generate_spm()` mutates every site and records the structural response.
Everything else works from its output.

`calculate_profiles()` and `calculate_decomposition()` evaluate the
model at selection strengths set to chosen values, to explore how each
constraint shapes the profiles. `predict_profiles()` and
`predict_decomposition()` evaluate it at a fit obtained from data, with
error bands. Fits come from `fit_lrmsd_msa_site()` and
`fit_lrmsd_msa_mode()`.

All four return a `$site` table and a `$mode` table, holding the
profiles on the residue axis and on the normal-mode axis.

Inputs are a PDB file, active-site residue numbers, and an observed
profile to fit against. Measuring that profile, by superimposing
homologs over an alignment, is outside this package.

## Documentation

``` r
?msamodel                        # API, indexed by workflow step
vignette("msamodel")             # both axes, end to end
vignette("msamodel-explore")     # profiles at chosen parameters, both axes
vignette("msamodel-fit")         # fitting to data, with confidence bands
```

## Reference

Echave J, Carpentier M (2026). Why structural divergence varies among
residues in enzyme evolution: contributions of mutation, stability, and
activity constraints. *Molecular Biology and Evolution* 43(7), msag162.
<https://doi.org/10.1093/molbev/msag162>

## License

MIT © Julian Echave
