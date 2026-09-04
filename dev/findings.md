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

1. **SPM physics** — `generate_spm_core()` (`@noRd`, `R/spm.R:54`). The expensive
   part: ENM + mutation scans via `penm`, producing per-mutant `dr2` and ΔΔG.
   Computed once, independent of `a1`/`a2`.
2. **Reshaping** — inlined in the public `generate_spm()` (`R/spm.R:171`). Not a
   computation: filters `m > 0`, sums energy columns into `ddg`/`ddgact`, builds the
   `[mutant × site]` and `[mutant × mode]` matrices and the `site ↔ pdb_site`
   `site_map`. Cheap, deterministic, `a1`/`a2`-independent. (Until 2026-08-12 this
   lived in two `@noRd` helpers, `preprocess_spm()` / `preprocess_spm_mode()`; they
   computed identical `energy_data` twice and had one caller each, so they were
   folded in.)
3. **Reweighting** — the axis-blind primitive `dr2_msa(dr2mat, energy_data, a1, a2)`
   (`R/model.R:87`), which takes either `dr2mat_site` or `dr2mat_mode`. The site/mode
   split is not in per-axis calculators: `calculate_*` and `predict_*` write the two
   axes out explicitly, passing `spm$dr2mat_site` / `spm$dr2mat_mode` directly.

`generate_spm()` returns a classed `spm` list `{energy_data, dr2mat_site,
dr2mat_mode, site_map, mode_map}` (stages 1+2 fused).

**Naming boundary.** The public object's names carry no index letters; the internal
scan's do. `generate_spm_core()` keeps `dr2_ijm` / `dr2_njm` list-columns (there the
letters state the shape exactly: one row per mutant `(j,m)`, each cell a vector over
response sites `i` / modes `n`), and `generate_spm()` is where they become
`dr2mat_site` / `dr2mat_mode`.

**In the forward map, `a1`/`a2` enter through one path only**: `weights_jm()`
(`R/model.R:69`) → `pfix_msa()` (`R/model.R:40`), where
`pstab = pmin(exp(-a1*ddg), 1)` and `pact = pmin(exp(-a2*ddgact), 1)`. Nothing in
stages 1–2 sees them.

Two qualifications, so this is not mistaken for a global "only place":

- The **se machinery** (`R/predict_se.R`) consumes the parameters separately, as a
  coordinate transform rather than a reweighting: `c(fit$a1, log2(fit$a2 + 1))`,
  written inline in each `var_param_*` function, with the inverse in
  `R/fitting.R:120-133`.
- The **nested variants** dispatch the same parameters into four combinations —
  `(0,0)`, `(a1,0)`, `(0,a2)`, `(a1,a2)` — for MM/MS/MA/MSA (`R/model.R:126-129`,
  and `var_spm_nested_nlrmsd()` in `R/predict_se.R`). Still one path, evaluated four
  times.

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
| structure divergence `dr2` | `dr2mat_site` | `dr2mat_mode` | `delta_structure_dr2i` / `dr2n` | yes |
| motion divergence `dh` | `dh_ijm` | `dh_njm` | `delta_motion_dhi` / `dhn` | no |
| mode fluctuation `nh` | — | `nh_njm` | `delta_motion_nhn` | no |

Verified 2026-08-05: no `dh_ijm`/`dh_njm`/`nh_njm` anywhere in `R/`. The only
fluctuation-adjacent code is `msf = get_msf_site(wt)` in `R/site_properties.R:34`,
a static wild-type descriptor, not a divergence arm.

## Example data: two proteins, two scripts, and no test fixtures

The package ships **no datasets** — there is no `data/` directory, and since
2026-09-04 no cached test fixtures either.

**User-facing example FILES** — `inst/extdata/`, named `<pdb>_<chain>_*` throughout:

- `1d6o_A.pdb` + `1d6o_A_{active_site,lrmsd_obs_site,lrmsd_obs_mode_syn}.csv`
  (107 residues, peptidyl-prolyl isomerase FKBP1A). Built by
  `data-raw/prepare_1d6o_data.R`. **This is the one everything uses** — vignettes,
  roxygen examples, and the test suite.

A 228-residue `1znb_A` set (metallo-beta-lactamase type 2) was shipped until 2026-09-04
and removed once the docs and tests had all moved to 1d6o_A and nothing read it. Its
files and `data-raw/prepare_znb_data.R` are recoverable from git history if a larger
example is ever wanted.

The `_syn` in the mode files is a warning label: no empirical per-mode profile exists,
so those are model output plus seeded noise. They exercise machinery, never science.

**Test fixtures** — none. `tests/testthat/helper-setup.R` builds the ENM and the scan
per run (~7.4 s on 1d6o_A) rather than loading a cache. The reason is not speed but
correctness: the suite's core is a set of drift pins on the public verbs' output, and a
cached `.rds` cannot notice a regression in `generate_spm()` — the pins would stay green
while the generator drifted underneath them. Building closes that by construction.
`generate_spm` is also pinned directly, on a cheap 2-mutation scan, so a generator change
is reported as one rather than inferred from a downstream verb failing.

This replaced `tests/testthat/fixtures/{znb_wt,znb_spm}.rds` (21 MB, shipped to every
user), `make-znb-fixtures.R`, `test-spm-generate.R`, the `MSAMODEL_FULL_TESTS` gate and
pre-commit gate 3 — all of which existed to police the cache. The built tarball went
from 24.6 MB to 3.2 MB.

`helper-setup.R` and `data-raw/prepare_1d6o_data.R` each declare the nine ENM/SPM
constants, with identical values, because `data-raw/` must not read `tests/`. Nothing
enforces the match — change both by hand.

The active-site vector lives only in the shipped CSV, read by both pipelines. It was
previously maintained in three places at once.

**Reading shipped files from tests: use `system.file()`, never a relative path.**
Under `R CMD check` the tests run against a BUILT package, where `inst/extdata/` has
been promoted to `extdata/` and `inst/` does not exist. `test_path("..", "..",
"inst", ...)` therefore passes under `devtools::test()` and ERRORs under
`devtools::check()`. `system.file()` is correct in both, because `pkgload` shims it
during `load_all()`.

## What the test suite claims, and what it does not

Rebuilt 2026-09-04: 21 `test_that` blocks in five files (~30 s), replacing 63 blocks
in eleven (~75 s). The organising rule is that a test earns its place only if it can
fail in a way you would not otherwise notice — a wrong divergence profile is visible on
sight; a swapped `phi_stab`/`phi_act` that still sums correctly is not.

- **`test-pins.R`** — whole-object snapshots of `generate_spm` and the six public verbs.
  These pin **stability, never correctness**: they freeze whatever the code produced, and
  a wrong value freezes as happily as a right one. Correctness came from the paper and
  the delta-method validation study (`dev/reports/`), which was research, done once, and
  is deliberately not re-run per suite.
- **`test-contract.R`** — the `pdb_site` → site-index join, on scrambled order and a
  strict subset. Necessary because the pins exercise one input: a positional
  implementation would return plausible, wrong numbers for any other shape of input.
  1znb_A gave this property by accident (it numbers from 20, so `pdb_site != site`);
  1d6o_A numbers 1..107, so the case is now constructed on purpose.
- **`test-errors.R`** — the fail-loud paths, encoding the policy that an unexpected
  input errors rather than returning a silent NA.
- **`test-invariants.R`** — identities that hold for any input, so they survive a
  deliberate re-freeze and still constrain the new values.

**Tolerances are set from measurement, not left at the default.** The two-basis identity
(`rowSums(dr2mat_site) == rowSums(dr2mat_mode)`) agrees to 3.0e-14 relative — 134x machine
epsilon — so it is pinned at 1e-12; the default was loose enough to pass a 1% error on a
single cell. Four other identities were measured to hold exactly and use
`expect_identical`.

**Snapshot descriptions must not be edited.** A snapshot is keyed by its `test_that`
description: renaming one orphans the recorded value, and testthat writes a NEW snapshot
with a warning rather than failing — the comparison silently stops happening.

**Known limit:** the pins catch a 1e-6 multiplicative error but not 1e-9, and an additive
constant is invisible to the centred (`n*`) columns because centring removes it — the
uncentred `lrmsd` snapshot is what catches that case.

## A CSV round-trip is not bit-exact for doubles

`write_csv()` emits ~15 significant digits, so a value read back differs from the
in-memory original by ~1e-16 (measured: 8.9e-16 max on the mode profile). 17
significant digits would round-trip losslessly, but that is unreadable in a
user-facing example file.

`data-raw/prepare_1d6o_data.R` therefore **reads the mode profile back from the CSV
and uses that** for everything downstream, so the shipped file is the source of
truth rather than a lossy copy of a value nothing else can see. The reproducibility
test (`test-fit-ml-mode.R`) uses `tolerance = 1e-12` for this reason and says so in
place — real drift moves numbers far above that.
