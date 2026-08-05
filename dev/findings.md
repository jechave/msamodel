# Durable findings

Knowledge about the model and the package that is expensive to rediscover. Facts
and reasoning only — no roadmap, no status, no history. `dev/LOG.md` holds the
history; `dev/ideas.md` holds things to maybe do later.

Function names below are checked against `R/` as of 2026-08-05.

---

## Why the package reports centred quantities (`nlrmsd`)

Fitting `lrmsd_obs = a0 + lrmsd_model + noise` with an explicit level parameter `a0`
was considered and dropped.

The reason: **it is not known how the level is shared among the nested variants
MM/MS/MA/MSA.** Any level-carrying quantity (`lrmsd_mm` and siblings) is therefore
uninterpretable, while the centred quantities (`nlrmsd_mm` …) are invariant to that
unknown *whatever the sharing structure turns out to be*.

So mean-centring is not a workaround for a missing parameter. It is the standard
move of reporting only quantities invariant to a nuisance you cannot pin down.

A related objection — that `nlrmsd` has "no canonical support" because the protein
has 228 sites but only 225 are matched — is weaker than it first appears. Needing to
*choose* a support is not the same as being ill-defined: fix the reference set once
(the 225 matched sites, which the fit already uses) and `nlrmsd(i)` is well-defined
on all 228, with the 3 unobserved sites (pdb_site 243/244/245) taking values
relative to that reference.

## The precompute-and-reweight architecture

Three stages, easy to conflate:

1. **SPM physics** — `generate_spm_core()` (`@noRd`, `R/spm.R:48`). The expensive
   part: ENM + mutation scans via `penm`, producing per-mutant `dr2` and ΔΔG.
   Computed once, independent of `a1`/`a2`.
2. **Reshaping** — `preprocess_spm()` / `preprocess_spm_mode()` (`@noRd`,
   `R/spm.R:199` / `R/spm.R:246`), composed inside the public `generate_spm_data()`.
   Not a computation: filters `m > 0`, sums energy columns into `ddg_jm`/`ddgact_jm`,
   builds the `[mutant × site]` and `[mutant × mode]` matrices and the `i ↔ pdb_site`
   `site_map`. Cheap, deterministic, `a1`/`a2`-independent.
3. **Reweighting** — the axis-blind primitive `dr2_msa(dr2_mat, energy_data, a1, a2)`
   (`R/model.R:85`), which takes either `dr2_ijm` or `dr2_njm`. The site/mode split
   lives in the dispatch table `axis_branches()` (`R/api.R:65`), not in separate
   per-axis calculators.

`generate_spm_data()` returns a classed `spm` list `{energy_data, dr2_ijm, dr2_njm,
site_map}` (stages 1+2 fused, `R/spm.R:173-180`).

**In the forward map, `a1`/`a2` enter through one path only**: `weights_jm()`
(`R/model.R:69`) → `pfix_msa()` (`R/model.R:40`), where
`pstab = pmin(exp(-a1*ddg), 1)` and `pact = pmin(exp(-a2*ddgact), 1)`. Nothing in
stages 1–2 sees them.

Two qualifications, so this is not mistaken for a global "only place":

- The **band machinery** consumes the parameters separately, as a coordinate
  transform rather than a reweighting: `as_theta(a1, a2) = c(a1, log2(a2 + 1))`
  (`R/predict_band_helpers.R:73`), used by `var_param_delta()` for the parameter arm
  of the standard errors, with the inverse in `R/fitting.R:120-133`.
- The **nested variants** dispatch the same parameters into four combinations —
  `(0,0)`, `(a1,0)`, `(0,a2)`, `(a1,a2)` — for MM/MS/MA/MSA (`R/model.R:126-129`,
  `R/predict_band_helpers.R:167-169`). Still one path, evaluated four times.

### Why this generalizes

Every divergence quantity in the model has the form

```
X = sum(pfix_jm * X_jm) / sum(pfix_jm)
```

The `X_jm` are mutant properties computed *before* selection (stage 1), so they do
not depend on `a1`/`a2`; the parameters enter only through the per-mutant `pfix_jm`.

Consequence: each new quantity is just **another value-matrix reweighted by the same
weights**. Adding an arm means adding `_jm` columns to the SPM, not new machinery.

### The limit of this property

The factorization depends on each mutant being scored independently of the others.
Any route where substitutions are *sequential* — each fixation depending on
accumulated state — breaks it, because the weights no longer factorize per-mutant
and there is no fixed matrix to reweight.

## The divergence grid

The full model spans divergence quantities × response axis {site `i`, mode `n`}.
Only the structure row is implemented.

| quantity | site form | mode form | penm function | implemented? |
|---|---|---|---|---|
| structure divergence `dr2` | `dr2_ijm` | `dr2_njm` | `delta_structure_dr2i` / `dr2n` | yes |
| motion divergence `dh` | `dh_ijm` | `dh_njm` | `delta_motion_dhi` / `dhn` | no |
| mode fluctuation `nh` | — | `nh_njm` | `delta_motion_nhn` | no |

Verified 2026-08-05: no `dh_ijm`/`dh_njm`/`nh_njm` anywhere in `R/`. The only
fluctuation-adjacent code is `msf = get_msf_site(wt)` in `R/site_properties.R:34`,
a static wild-type descriptor, not a divergence arm.

## The `znb_*` fixture is generated, not copied

`data-raw/prepare_znb_data.R` builds the embedded datasets by calling the package's
own code (`setup_enm()`, then `generate_spm_data()`), so the fixture stays consistent
with the package rather than with any external source. It is frozen — regenerate it
intentionally, not incidentally.

`tests/testthat/test-spm-generate.R` catches drift between the SPM-generation code
and the fixture, so that test and the data-prep script must stay in sync. Seed and
scan constants are deliberately duplicated in both files.
