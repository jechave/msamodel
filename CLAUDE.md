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
- **`dr2`-family index-signature convention.** Every `dr2`-family name msamodel
  *creates* is `dr2_<indices>`: one underscore, then exactly the free indices the
  object spans, in order **response index (`i` site / `n` mode), mutated site `j`,
  mutation `m`**. A reduction over an axis drops that axis's letter (`dr2_ijm` →
  `dr2_i`). Transform prefixes (`l`=log, `n`=normalized) and source labels
  (`_msa`/`_obs`) are kept, index still underscore-set: `lrmsd_i`, `nlrmsd_i_msa`,
  `dr2_i_msa`. This governs names msamodel **creates**, not what it **calls** — call
  penm's own names (`penm::delta_structure_dr2i`) directly; never wrap or rename a
  dependency. Applies to future motion-arm names (`dh_*`, `nh_*`) too.
- **Roxygen on every function.** `#' @export` for public API, `#' @noRd` for internal
  helpers. Document `@param`, `@return`/`@returns`, add `@family` tags to group
  related functions. `@examples` in `\dontrun{}` when they need real data.
- **Imports live in one place:** `R/msamodel-package.R` carries the `"_PACKAGE"` doc
  and the `@importFrom` directives. Re-export `%>%` via `#' @importFrom magrittr %>%`.
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

# AT A MILESTONE ONLY — full check (baseline: 0 errors, 1 warning, 2 notes,
# deliberately accepted, GitHub-only not CRAN).
Rscript -e "devtools::check()"

# Install locally (only when a vignette needs the working tree's new functions)
Rscript -e "devtools::install()"

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
- **`.Rmd.orig` precompute pattern.** Edit `<name>.Rmd.orig` → knit to `<name>.Rmd`
  (`Rscript -e 'knitr::knit("vignettes/<name>.Rmd.orig", output="vignettes/<name>.Rmd")'`)
  → preview with `dev/preview-vignette.R` → user approves → commit. Install the working
  tree first if the vignette calls new functions (`library(msamodel)` loads the
  *installed* package). **Never render the committed `vignettes/<name>.Rmd` in place** —
  `html_vignette` renders `--self-contained`, which base64-embeds figures and then
  DELETES the `<name>_files/` dir the shipped `.Rmd` references. Use
  `dev/preview-vignette.R` (renders a temp-dir copy; leaves repo figures intact).

Embedded `znb_*` datasets are (re)generated by `data-raw/prepare_znb_data.R`. The
fixture is **frozen** — regenerate it intentionally, not incidentally.
`test-spm-generate.R` catches drift between the SPM-generation code and the fixture,
so the data-prep script and that test must stay in sync.

## Working style

- Don't hand-wave uncertainty — flag anything unverified. Only say
  "confirmed/verified/fixed" when you actually checked.
- Own contradictions across turns instead of smoothing them over.
- Push back on bad ideas rather than silently implementing them.
- Terse is good; real uncertainty should be stated explicitly.
- **Tool choice:** read files with `Read` (use `offset`/`limit` for big files). Avoid
  `sed`/`head`/`cat`/`tail`/`awk`/`echo` in Bash for reading — they trigger permission
  prompts here.
