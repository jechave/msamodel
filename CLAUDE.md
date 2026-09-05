# CLAUDE.md

Guidance for Claude Code when working in this repository.

## What this package is

`msamodel` implements the Mutation-Stability-Activity (MSA) model of neutral
structural evolution: it predicts site-specific structural divergence in protein
evolution from single-point-mutation scans (SPMs). Structure/active-site input is a
bio3d pdb object + a plain `pdb_site_active` integer vector.

The package computes structure-divergence profiles (by site and by mode) and fits
them to data: SPM generation, divergence profiles, parameter estimation (ML) with
delta-method standard errors, and a decomposition into mutation/stability/activity
contributions (the `phi_*` / `nphi_*` columns).

It depends on `penm` (an `Imports:` dependency) for the ENM machinery — `set_enm()`,
mutation scans, response matrices. **penm is a dependency, never a source**: call
`penm::fn()`; never copy or wrap its code, and never rename a penm function to fit a
local convention.

## Test discipline

A permanent test must be able to **fail for a real reason**. Before trusting any new
invariant/value test, run a negative control: feed a deliberately-wrong input and
confirm the assertion goes red. If it can't go red, it's a tautology — drop it. A
"numbers didn't change" refactor-invariance check is a one-time scratchpad artifact,
not a permanent suite test.

## House conventions

- **Naming:** `snake_case`.
- **`dr2`-family naming — index letters inside, axis words outside.** Two conventions,
  split at `generate_spm()`; both are deliberate, neither should be "unified".
  - **Internal, per-mutant quantities** use the index signature `dr2_<indices>`: one
    underscore, then exactly the free indices the object spans, in order **response
    index (`i` site / `n` mode), mutated site `j`, mutation `m`**. A reduction over an
    axis drops that axis's letter (`dr2_ijm` → `dr2_i`). Transform prefixes (`l`=log,
    `n`=normalized) and source labels (`_msa`/`_obs`) are kept, index still
    underscore-set: `lrmsd_i`, `nlrmsd_i_msa`, `dr2_i_msa`. These live in `R/spm.R`.
  - **The public `spm` object's parts** use an axis word instead: `dr2mat_site`,
    `dr2mat_mode` (`[mutant x response]` matrices), with `site_map` / `mode_map`. These
    are what the rest of the package handles — `model.R`, `fitting.R`, `predict.R`,
    `predict_se.R` — and they are correct as-is. Do **not** rename them to index
    letters.
  - `generate_spm()` is the crossing point: it stacks the per-mutant `dr2_ijm` /
    `dr2_njm` list-columns into `dr2mat_site` / `dr2mat_mode` (`R/spm.R:206-246`,
    which documents the boundary in place).

  This governs names msamodel **creates**, not what it **calls** — call penm's own names
  (`penm::delta_structure_dr2i`) directly; never wrap or rename a dependency. A future
  motion arm follows the same split (`dh_*` / `nh_*` internally).
- **Roxygen on every function.** `#' @export` for public API, `#' @noRd` for internal
  helpers. Document `@param`, `@return`/`@returns`, add `@family` tags to group
  related functions. `@examples` in `\dontrun{}` when they need real data.
- **Imports live in one place:** `R/msamodel-package.R` carries the `"_PACKAGE"` doc
  and the `@importFrom` directives. Re-export `%>%` via `#' @importFrom magrittr %>%`.
- **Pipes split by layer:** `R/` uses magrittr `%>%` (a real `Imports:`); the vignettes
  use the native `|>` throughout. Match the layer you are editing.
- **Tibbles** (not data.frames) for tabular returns; tidyverse for data manipulation.

## Development commands

```bash
# INNER LOOP (constant) — load_all + targeted tests; this is the working loop
Rscript -e "devtools::load_all()"
Rscript -e "testthat::test_file('tests/testthat/test-fit-ml.R')"    # one file
# (or devtools::test(filter='fit-ml') for a name-filtered subset)

# Document (regenerate NAMESPACE + man/) — ONLY after a roxygen/@importFrom change
Rscript -e "devtools::document()"

# Full suite — the gate before a code/data/roxygen commit
Rscript -e "devtools::test()"

# AT A MILESTONE ONLY — full check. No baseline count is recorded anywhere: the
# old "0E/1W/3N" one here went stale (its 1 WARNING was the installed-size one,
# cleared 2026-08-11 when the datasets moved out, 25.8 MB -> 368 KB, see NEWS.md).
# Read the current result off the run itself rather than trusting a number nobody
# has re-measured. Standing accepted note: the NSE "no visible binding" one, left
# visible instead of suppressed with globalVariables() — see the comment in
# R/msamodel-package.R.
Rscript -e "devtools::check()"

# Build the vignettes and render them to HTML for review -> doc/<name>.html
# (doc/ and Meta/ are git- and build-ignored). This RUNS the vignette code, so it
# is also the check that the vignettes still work against the current package.
Rscript -e "devtools::build_vignettes()"

# One vignette on its own, when iterating (from inside vignettes/)
Rscript -e "rmarkdown::render('msamodel-explore.Rmd')"
```

NAMESPACE and everything under `man/` are roxygen-generated — edit the roxygen
comments and run `document()`, never hand-edit them.

### Vignettes — HARD RULES

- **Vignettes are ordinary `.Rmd` files that run at build time**, like any R package.
  There is no precompute step, no `.Rmd.orig`, and no committed figures: `vignettes/`
  holds exactly the four `.Rmd` sources. `R CMD check` executes them, so a vignette that
  breaks against the current code **fails the check** — which is the point.
- **Never commit/push `vignettes/` without the user's explicit HTML approval.** Do NOT
  commit any change touching `vignettes/` until the **user** has said they read the
  rendered HTML and approved it. "Tests pass", or *you* glancing at a render, do not
  satisfy this — only the user's explicit OK. Build with
  `Rscript -e "devtools::build_vignettes()"`, tell the user exactly which
  `doc/<name>.html` to open, then STOP and wait. `.githooks/pre-commit` gate 1 enforces
  this mechanically (`VIGNETTE_APPROVED=1`).
- **`doc/` and `Meta/` are build output** — git-ignored and `.Rbuildignore`d. Never
  commit them; regenerate on demand.
- **Never `rm` untracked paths inside the repo.** Report what they are and let the user
  decide. On 2026-08-06 four stray `*_files/` dirs were `rm -rf`'d as "byproducts"; they
  were the figures a knit had just produced — the only evidence of what it did. Nothing
  enforces this now (the guard hook was removed 2026-09-02, along with the precompute
  that made those figures exist): **it is discipline, not a mechanism.** `git rm` on a
  tracked file is fine — git keeps the blob and the deletion is reviewable in the diff.

### Example data — one protein, one script

The package ships **no datasets**; there is no `data/`, and since 2026-09-04 there are
no cached test fixtures either.

- **User-facing example FILES** → `inst/extdata/`: `1d6o_A.pdb` +
  `1d6o_A_{active_site,lrmsd_obs_site,lrmsd_obs_mode_syn}.csv` (107 residues), built by
  `data-raw/prepare_1d6o_data.R`. Vignettes, roxygen examples and tests all read them via
  `system.file()` — the same call a user makes on their own data. This is the only
  example protein: a 228-residue `1znb_A` set was removed 2026-09-04 once the docs and
  tests had moved off it (recoverable from git history if a larger example is wanted).
- **Test fixtures** → none. `tests/testthat/helper-setup.R` BUILDS the ENM and the scan
  on every run (~7.4 s) rather than loading a cache, and owns the ENM/SPM constants for
  the test side. It also builds `spm_small` (`n_mutations = 2`, ~1.8 s), used by exactly
  one snapshot — the pin on `generate_spm()` itself, so a generator change is reported
  as a generator change rather than only as a downstream verb failing.

`data-raw/` must not read from `tests/`. If a new artifact seems to need both, that is a
signal it belongs to one side only — decide which, do not bridge them.

**Duplicated constants, still open on one side.** `data-raw/prepare_1d6o_data.R` and
`helper-setup.R` each declare the nine ENM/SPM constants, with identical values, because
`data-raw/` may not read `tests/`. Nothing enforces the match, so editing one side leaves
the suite green while `inst/extdata/` ships from different constants — change both by
hand. (The old `make-znb-fixtures.R` was a third declaration; deleting it removed one of
the three, not the tension.)

**Why the scan is built, not cached.** The suite's core is a set of drift pins on the
public verbs' output. That design only works if the scan is rebuilt per run: a cached
`.rds` cannot notice a regression in `generate_spm()`, so the pins would stay green while
the generator drifted underneath them. Building closes it by construction — which is why
`test-spm-generate.R`, `MSAMODEL_FULL_TESTS` and pre-commit gate 3 are all gone. They
existed only to police the cache.

**Reading shipped files from tests: `system.file()`, never `test_path("..","..")`.**
Under `check()` the tests run against a built package where `inst/extdata/` has
become `extdata/` — a relative path passes `test()` and ERRORs in `check()`.

## Working style

- Don't hand-wave uncertainty — flag anything unverified. Only say
  "confirmed/verified/fixed" when you actually checked.
- Own contradictions across turns instead of smoothing them over.
- Push back on bad ideas rather than silently implementing them.
- Terse is good; real uncertainty should be stated explicitly.

## Working agreement

The section above states dispositions. Dispositions do not survive contact with
producing output at volume — on 2026-08-06 every rule above was already written and
every one was violated. These four are **procedures with checkable outcomes**. Each
exists because its absence cost the maintainer hours that day.

1. **Outside the plan: say it, do not do it.**
   If something not in the plan looks like it wants changing — even trivially, even
   while already editing that file — report it and stop. Approval of a plan is not
   approval of what the plan reminded you of. *(2026-08-06: the plan said MOVE
   `calculate_*`; they were rewritten, and readable code became opaque.)*

2. **The claim may not exceed the check. Every report carries a "not verified" line.**
   Name the scope: which files, which vignette, which axis. "All four" and "the
   vignettes" are the phrases to distrust in your own drafts. Put the output on screen
   in the same message as the claim about it. A summary table of green rows implies the
   whole job is green — anything unverified gets its own row saying so. *(2026-08-06:
   "only sessionInfo reflow" was true of the `.Rmd` files and false of the figures,
   which had never been looked at.)*

3. **One question at a time. Questions are not bundled with recommendations.**
   If a decision is needed, ask it alone and wait. Do not attach a recommendation, a
   second question, or a plan for what happens after the answer.

4. **An unexplained artifact is evidence. Never tidy it away.**
   An unexpected file, directory, or output is information about what just happened.
   Report what it is and what you think created it; let the maintainer decide.
   *(2026-08-06: four unexplained directories were deleted as "byproducts" — they were
   the only record of what a knit had produced. A hook enforced this for a while; it
   was removed 2026-09-02 with the precompute that created those directories, so this
   is discipline again, not a mechanism.)*
