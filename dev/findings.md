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

## Example data and test fixtures are separate, and so are their pipelines

The package ships **no datasets** — there is no `data/` directory. Two different
kinds of artifact, in two places, built by two scripts that do not read each other:

**User-facing example FILES** — `inst/extdata/`: `1znb_A.pdb`,
`znb_active_site.csv` (`pdb`, `chain`, `pdb_site`; one residue per row),
`znb_lrmsd_obs_site.csv` (empirical), `znb_lrmsd_obs_mode_syn.csv` (synthetic —
the `_syn` is the warning label, since no empirical per-mode profile exists). Built
by `data-raw/prepare_znb_data.R`. They are files, not `.rda`, so the vignettes read
them with the same `read_csv()` a user runs on their own data.

**Test fixtures** — `tests/testthat/fixtures/`: `znb_wt.rds`, `znb_spm.rds`, built
by `make-znb-fixtures.R` beside them (the layout R Packages (2e) prescribes). No
vignette uses either; they exist so the suite does not spend ~21 s per run
regenerating the scan. `helper-fixtures.R` loads them once for all test files.

Both are generated by calling the package's own code (`setup_enm()`, then
`generate_spm()`), so they stay consistent with the package rather than with an
external source. The scan fixture is frozen — regenerate it intentionally, not
incidentally.

`make-znb-fixtures.R` **owns** the ENM/SPM constants; `test-spm-generate.R` sources
that file to read them, so there is one declaration rather than two under a "MUST
match" comment. `data-raw/` declares its own copy deliberately — it must not depend
on `tests/`, and it regenerates the scan rather than borrowing the fixture. That
duplication is checked by the drift guard, which compares a fresh scan to the
fixture.

The active-site vector lives only in `znb_active_site.csv`. It was previously
maintained in three places at once (data-raw, the test, `znb_dataset`).

**Reading shipped files from tests: use `system.file()`, never a relative path.**
Under `R CMD check` the tests run against a BUILT package, where `inst/extdata/` has
been promoted to `extdata/` and `inst/` does not exist. `test_path("..", "..",
"inst", ...)` therefore passes under `devtools::test()` and ERRORs under
`devtools::check()`. `system.file()` is correct in both, because `pkgload` shims it
during `load_all()`.

## The SPM drift guard, and what each tier can see

`tests/testthat/test-spm-generate.R` holds two checks against the fixture, and the
difference between them matters:

- **Always-on**: reproduces ONE mutant `(j=1, m=3)` from its seeded recipe. Cheap,
  and it locates its row by lookup — so it survives, and therefore cannot detect, a
  change that reorders rows or alters which are included.
- **Gated** behind `MSAMODEL_FULL_TESTS=true` (~21 s): regenerates all 2280 rows and
  compares the whole object. This is the one that sees structural change.

There is **no CI in this repo**. The gated check runs from `.githooks/pre-commit`
gate 3, automatically, on commits that stage `R/spm.R`, `R/enm_setup.R`, the fixture
recipe, or `data-raw/prepare_znb_data.R` — plus by hand at a milestone.

## A CSV round-trip is not bit-exact for doubles

`write_csv()` emits ~15 significant digits, so a value read back differs from the
in-memory original by ~1e-16 (measured: 8.9e-16 max on the mode profile). 17
significant digits would round-trip losslessly, but that is unreadable in a
user-facing example file.

`data-raw/prepare_znb_data.R` therefore **reads the mode profile back from the CSV
and uses that** for everything downstream, so the shipped file is the source of
truth rather than a lossy copy of a value nothing else can see. The reproducibility
test (`test-fit-ml-mode.R`) uses `tolerance = 1e-12` for this reason and says so in
place — real drift moves numbers far above that.
