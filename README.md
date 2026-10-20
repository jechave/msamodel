
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

Run this at an R prompt, not in the shell:

``` r
# install.packages("remotes")
remotes::install_github("jechave/msamodel", dependencies = TRUE)
```

The elastic network model comes from
[penm](https://github.com/jechave/penm), installed alongside. Neither
package is on CRAN.

Every code block below is R as well, to be run in the same session after
`library(msamodel)`.

## Interface

The inputs are a PDB file and a vector of active-site residue numbers.
Fitting also takes an observed divergence profile, as two vectors: the
residue numbers, and the divergence at each. Measuring that profile, by
superimposing homologs over an alignment, is outside this package.

| Function | Does | Returns |
|----|----|----|
| `generate_spm()` | mutates every site and records the structural response | the scan every other function reads |
| `fit_lrmsd_msa_site()` | estimates the two selection strengths by fitting an observed site profile | the fit, with standard errors and `$gof` |
| `fit_lrmsd_msa_mode()` | estimates them by fitting an observed mode profile | the fit, with standard errors and `$gof` |
| `calculate_profiles()` | evaluates the model at selection strengths given as arguments | the profiles |
| `calculate_decomposition()` | decomposes the model at selection strengths given as arguments | the nested models and the three contributions |
| `predict_profiles()` | evaluates the model at a fit | the profiles, with standard errors |
| `predict_decomposition()` | decomposes the model at a fit | the nested models and the three contributions, with standard errors |

The `calculate_*` pair explores how each constraint shapes the profiles;
the `predict_*` pair evaluates the model where the data puts it.

Every profile comes in two representations, and the `calculate_*` and
`predict_*` functions return both: the same divergence written in two
bases. The `$site` tibble gives it residue by residue, the `$mode`
tibble gives it per normal mode. Their `metric` argument chooses the
quantity: `"lrmsd"` keeps the profile’s overall level, `"nlrmsd"`
centres it on its mean, which is the form a fit is made on.

## Documentation

`?msamodel` is the package page: the workflow in order, with a link to
each function.

Three vignettes, each a worked example on the shipped protein:

``` r
vignette("msamodel")             # the whole workflow, scan to decomposition
vignette("msamodel-explore")     # the model at chosen selection strengths
vignette("msamodel-fit")         # estimating those strengths from data
```

## Example

`penm::set_enm()` builds the elastic network from the structure, and
`generate_spm()` mutates every site and records how the structure
responds.

``` r
library(msamodel)
library(penm)

ex     <- function(f) system.file("extdata", f, package = "msamodel")
active <- read.csv(ex("1d6o_A_active_site.csv"))
obs    <- read.csv(ex("1d6o_A_lrmsd_obs_site.csv"))

wt  <- penm::set_enm(bio3d::read.pdb(ex("1d6o_A.pdb")), node = "ca",
                     model = "ming_wall", d_max = 10.5, ensemble = 1L)
spm <- generate_spm(wt, pdb_site_active = active$pdb_site)
```

`fit_lrmsd_msa_site()` estimates the two selection strengths from an
observed profile, given as two vectors: the residue numbers and the
divergence at each. The fit carries the estimates, their standard
errors, and a goodness-of-fit summary in `$gof`.

``` r
fit <- fit_lrmsd_msa_site(spm, obs$pdb_site, obs$lrmsd_obs)

c(a1 = fit$a1, a2 = fit$a2)
#>         a1         a2 
#>  0.3675052 87.8156581
fit$gof
#> # A tibble: 1 × 9
#>      D2   AIC   BIC logLik deviance null_deviance  nobs     k sigma_hat
#>   <dbl> <dbl> <dbl>  <dbl>    <dbl>         <dbl> <int> <int>     <dbl>
#> 1 0.602  114.  122.  -53.9     17.2          43.1   107     3     0.401
```

`predict_profiles()` evaluates the model at that fit. It returns a
`$site` and a `$mode` tibble, each holding the profile and its standard
error in one representation.

``` r
pred <- predict_profiles(fit, spm, metric = "nlrmsd")

pred$site
#> # A tibble: 107 × 4
#>     site pdb_site nlrmsd_msa nlrmsd_msa_se
#>    <int>    <int>      <dbl>         <dbl>
#>  1     1        1     -0.142        0.123 
#>  2     2        2     -0.473        0.109 
#>  3     3        3     -0.126        0.124 
#>  4     4        4     -0.337        0.0925
#>  5     5        5     -0.138        0.0938
#>  6     6        6      0.384        0.134 
#>  7     7        7      0.388        0.118 
#>  8     8        8      0.650        0.105 
#>  9     9        9      0.713        0.103 
#> 10    10       10      0.493        0.128 
#> # ℹ 97 more rows
pred$mode
#> # A tibble: 315 × 3
#>     mode nlrmsd_msa nlrmsd_msa_se
#>    <int>      <dbl>         <dbl>
#>  1     1       2.36        0.0586
#>  2     2       2.27        0.0568
#>  3     3       1.98        0.0536
#>  4     4       1.98        0.0611
#>  5     5       1.98        0.0437
#>  6     6       1.77        0.0611
#>  7     7       1.99        0.0431
#>  8     8       1.98        0.0422
#>  9     9       1.82        0.142 
#> 10    10       1.77        0.0574
#> # ℹ 305 more rows
```

Plotting both: the observed profile exists only per residue, so the
points appear in the left panel alone.

<img src="man/figures/README-profiles-1.png" width="100%" />

`predict_decomposition()` returns more columns in the same two
representations: the profiles of the four nested models, in `nlrmsd_mm`,
`nlrmsd_ms`, `nlrmsd_ma` and `nlrmsd_msa`; the three contributions, in
`nphi_mut`, `nphi_stab` and `nphi_act`; and a standard error for each.

``` r
dec <- predict_decomposition(fit, spm, metric = "nlrmsd")

dec$site
#> # A tibble: 107 × 16
#>     site pdb_site nlrmsd_mm nlrmsd_ms nlrmsd_ma nlrmsd_msa nphi_mut nphi_stab
#>    <int>    <int>     <dbl>     <dbl>     <dbl>      <dbl>    <dbl>     <dbl>
#>  1     1        1   -0.0304   -0.0381   -0.137      -0.142  -0.0304  -0.00770
#>  2     2        2   -0.201    -0.324    -0.362      -0.473  -0.201   -0.123  
#>  3     3        3   -0.215    -0.262    -0.0920     -0.126  -0.215   -0.0471 
#>  4     4        4   -0.310    -0.426    -0.245      -0.337  -0.310   -0.115  
#>  5     5        5   -0.202    -0.243    -0.103      -0.138  -0.202   -0.0415 
#>  6     6        6   -0.0932   -0.0473    0.327       0.384  -0.0932   0.0459 
#>  7     7        7   -0.0199    0.0135    0.355       0.388  -0.0199   0.0334 
#>  8     8        8    0.131     0.135     0.665       0.650   0.131    0.00402
#>  9     9        9    0.164     0.174     0.729       0.713   0.164    0.0101 
#> 10    10       10   -0.0554   -0.0529    0.522       0.493  -0.0554   0.00245
#> # ℹ 97 more rows
#> # ℹ 8 more variables: nphi_act <dbl>, nlrmsd_mm_se <dbl>, nlrmsd_ms_se <dbl>,
#> #   nlrmsd_ma_se <dbl>, nlrmsd_msa_se <dbl>, nphi_mut_se <dbl>,
#> #   nphi_stab_se <dbl>, nphi_act_se <dbl>
dec$mode
#> # A tibble: 315 × 15
#>     mode nlrmsd_mm nlrmsd_ms nlrmsd_ma nlrmsd_msa nphi_mut nphi_stab nphi_act
#>    <int>     <dbl>     <dbl>     <dbl>      <dbl>    <dbl>     <dbl>    <dbl>
#>  1     1      2.94      2.94      2.32       2.36     2.94   0.00167  -0.578 
#>  2     2      2.74      2.78      2.23       2.27     2.74   0.0443   -0.508 
#>  3     3      2.32      2.36      1.94       1.98     2.32   0.0442   -0.386 
#>  4     4      2.28      2.30      1.93       1.98     2.28   0.0267   -0.323 
#>  5     5      2.15      2.16      1.95       1.98     2.15   0.0129   -0.183 
#>  6     6      2.06      2.08      1.72       1.77     2.06   0.0202   -0.316 
#>  7     7      2.03      2.05      1.95       1.99     2.03   0.0162   -0.0615
#>  8     8      1.92      1.94      1.95       1.98     1.92   0.0220    0.0380
#>  9     9      1.90      1.99      1.67       1.82     1.90   0.0880   -0.165 
#> 10    10      1.86      1.93      1.71       1.77     1.86   0.0670   -0.152 
#> # ℹ 305 more rows
#> # ℹ 7 more variables: nlrmsd_mm_se <dbl>, nlrmsd_ms_se <dbl>,
#> #   nlrmsd_ma_se <dbl>, nlrmsd_msa_se <dbl>, nphi_mut_se <dbl>,
#> #   nphi_stab_se <dbl>, nphi_act_se <dbl>
```

MM, MS and MSA form a progression: MM has mutation alone, MS adds
selection on stability, and MSA adds selection on activity. Each is the
previous one with one more constraint switched on. MA, activity without
stability, lies outside this progression and is left out of the figure.

<img src="man/figures/README-nested-1.png" width="100%" />

The MSA profile splits into three contributions, one per constraint,
each the increment its step adds: mutation is MM, stability is MS minus
MM, and activity is MSA minus MS.

<img src="man/figures/README-decomposition-1.png" width="100%" />

Dashed lines mark active-site residues, and bands are 95% intervals. The
mode panels show the first 50 of 315 modes. `nlrmsd` is centred over all
315, and these are the modes carrying the largest divergence, so the
cropped panels do not appear centred on zero.

## Reference

Echave J, Carpentier M (2026). Why structural divergence varies among
residues in enzyme evolution: contributions of mutation, stability, and
activity constraints. *Molecular Biology and Evolution* 43(7), msag162.
<https://doi.org/10.1093/molbev/msag162>

## License

MIT © Julian Echave
