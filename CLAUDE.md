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

## Development context files

None of these is read at session start — the user says what the session is about.

- `dev/findings.md` — durable model/package knowledge that is expensive to
  rediscover (the precompute architecture, why the package reports centred
  quantities, the divergence grid).
- `dev/LOG.md` — append-only history, newest first. Consult when you need history.
  Entries describe the code *as it was on their date*; older function names are
  frequently obsolete.
- `dev/ideas.md` — parking lot for things to maybe do later. Not a plan.

For what is currently true of the code, read the code.

## Test discipline

A permanent test must be able to **fail for a real reason**. Before trusting any new
invariant/value test, run a negative control: feed a deliberately-wrong input and
confirm the assertion goes red. If it can't go red, it's a tautology — drop it. A
"numbers didn't change" refactor-invariance check is a one-time scratchpad artifact,
not a permanent suite test. The **`/test-review` skill** holds the full checklist.

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

# AT A MILESTONE ONLY — full check. The old "0E/1W/3N" baseline recorded here was
# stale: the 1 WARNING was the installed-size one, cleared 2026-08-11 when the
# datasets moved out (25.8 MB -> 368 KB, see NEWS.md). Rather than carry a number
# nobody has re-measured, take the current baseline from the newest check() entry
# in dev/LOG.md, and update that entry when you run one. Standing accepted note:
# the NSE "no visible binding" one, left visible instead of suppressed with
# globalVariables() — see the comment in R/msamodel-package.R.
Rscript -e "devtools::check()"

# Knit a vignette .Rmd.orig -> .Rmd. ALWAYS via this wrapper (it cds into vignettes/
# so figures land where the .Rmd points, and prints the figure timestamps it wrote).
# It ALSO runs devtools::install() itself, so there is no separate install step --
# note that install becomes the build every other project here gets from
# library(msamodel). <name> is required; no argument is an error.
# A PreToolUse hook blocks hand-rolled knitr::knit() calls.
Rscript dev/knit-vignette.R mode-analysis

# Preview a pre-rendered vignette as standalone HTML, WITHOUT touching the repo
Rscript dev/preview-vignette.R mode-analysis   # -> dev/preview/<name>.html (git-ignored)
```

NAMESPACE and everything under `man/` are roxygen-generated — edit the roxygen
comments and run `document()`, never hand-edit them.

### Vignettes — HARD RULES

- **Never commit/push `vignettes/` without the user's explicit HTML approval.** Do
  NOT commit any change touching `vignettes/` (`.Rmd`, `.Rmd.orig`, `_files/`) until
  the **user** has said they previewed the rendered HTML and approved it. Re-knitting,
  "tests pass", or *you* glancing at a preview do not satisfy this — only the user's
  explicit OK does. When vignettes change: regenerate **every** affected preview, tell
  the user exactly which `dev/preview/<name>.html` to open, then STOP and wait. The
  `.githooks/pre-commit` hook enforces this mechanically.
- **`.Rmd.orig` precompute pattern.** Edit `<name>.Rmd.orig` → **`Rscript dev/knit-vignette.R
  <name>`** → preview with `dev/preview-vignette.R` → user approves → commit. Install the
  working tree first if the vignette calls new functions (`library(msamodel)` loads the
  *installed* package). **Never render the committed `vignettes/<name>.Rmd` in place** —
  `html_vignette` renders `--self-contained`, which base64-embeds figures and then
  DELETES the `<name>_files/` dir the shipped `.Rmd` references. Use
  `dev/preview-vignette.R` (renders a temp-dir copy; leaves repo figures intact).
- **Always knit via `dev/knit-vignette.R`, never a hand-rolled `knitr::knit()` call.**
  knitr writes figures RELATIVE TO THE WORKING DIRECTORY, so knitting from the package
  root silently writes them to `<root>/<name>_files/` and leaves the real
  `vignettes/<name>_files/` STALE — the knitted `.Rmd` then points at old images and
  nothing errors (2026-08-06). The wrapper `cd`s correctly and prints the figure
  timestamps it wrote, so "the figures were regenerated" is evidence you can see.
  A `PreToolUse` hook (`.claude/hooks/guard-bash.sh`) blocks the hand-rolled form.
- **Never `rm` untracked paths inside the repo.** Report what they are and let the user
  decide. On 2026-08-06 four stray `*_files/` dirs were `rm -rf`'d as "byproducts"; they
  were the figures the knit had just produced — the only evidence of what it did.
  `guard-bash.sh` helps but is cruder than the rule: it **text-matches** `rm` and
  `git clean` and blocks them, exempting `git rm` and paths under `/tmp/`,
  `/private/tmp/`, or `scratchpad`. It never consults git, so it cannot tell tracked
  from untracked, and it will block a legitimate `rm` outside the repo too. **If it
  fires, that is a stop — report and ask; do not rephrase the command to slip past it.**
  **`git rm` on a tracked file is fine** — git keeps the blob and the deletion is
  reviewable in the diff; don't hand an approved deletion back to the user to type.

### Example data vs test fixtures — two places, two scripts

The package ships **no datasets**; there is no `data/`. Keep the two kinds apart:

- **User-facing example FILES** → `inst/extdata/` (`1znb_A.pdb`,
  `znb_active_site.csv`, `znb_lrmsd_obs_site.csv`, `znb_lrmsd_obs_mode_syn.csv`),
  built by `data-raw/prepare_znb_data.R`. Vignettes read them via `system.file()`,
  the same call a user makes on their own data.
- **Test fixtures** → `tests/testthat/fixtures/` (`znb_wt.rds`, `znb_spm.rds`),
  built by `make-znb-fixtures.R` beside them, which **owns** the ENM/SPM constants.

`data-raw/` must not read from `tests/`, and the fixture recipe must not write to
`inst/`. If a new artifact seems to need both, that is a signal it belongs to one
side only — decide which, do not bridge them.

**Open tension, not yet resolved — do not "tidy" either side.** The nine ENM/SPM
constants `make-znb-fixtures.R:27-36` owns are *also* declared, with identical values, at
`data-raw/prepare_znb_data.R:93-102`, which needs them to regenerate the scan
(`:108-116`) and may not read `tests/`. So the second declaration is forced by the
independence rule, not sloppiness — but nothing enforces the match (no test reads
`data-raw/`), so editing one side leaves the suite green while `inst/extdata/` ships from
different constants. The active site had this same shape and was fixed by moving it to
`inst/extdata/znb_active_site.csv`, which both sides read legally
(`data-raw/prepare_znb_data.R:43-46`); the constants have had no equivalent home chosen.
Until one is, **change both by hand — and do not add a comment that presents that as the
settled convention.**

The scan fixture is **frozen** — regenerate it intentionally, not incidentally
(`Rscript tests/testthat/fixtures/make-znb-fixtures.R`). `test-spm-generate.R`
catches drift; its full 2280-row check is gated behind `MSAMODEL_FULL_TESTS` and
fires automatically from `.githooks/pre-commit` gate 3 on SPM/ENM commits. The
always-on cheap check verifies ONE row of 2280 and cannot see a reordering.

**Reading shipped files from tests: `system.file()`, never `test_path("..","..")`.**
Under `check()` the tests run against a built package where `inst/extdata/` has
become `extdata/` — a relative path passes `test()` and ERRORs in `check()`.

## Working style

- Don't hand-wave uncertainty — flag anything unverified. Only say
  "confirmed/verified/fixed" when you actually checked.
- Own contradictions across turns instead of smoothing them over.
- Push back on bad ideas rather than silently implementing them.
- Terse is good; real uncertainty should be stated explicitly.
- **Tool choice:** read files with `Read` (use `offset`/`limit` for big files). Avoid
  `sed`/`head`/`cat`/`tail`/`awk`/`echo` in Bash for reading — they trigger permission
  prompts here.

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
   the only record of what a knit had produced. Now also blocked by
   `.claude/hooks/guard-bash.sh`.)*
