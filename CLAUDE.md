# CLAUDE.md

Guidance for Claude Code when working in this repository.

## What this package is

`msamodel` implements the Mutation-Stability-Activity (MSA) model of neutral
structural evolution: it predicts site-specific structural divergence in protein
evolution from single-point-mutation scans (SPMs).

**v0.1 is shipped.** It computes structure-divergence profiles (structure × site)
and fits them to data: SPM generation, structure-divergence profiles, a1/a2 grid
exploration, Bayesian parameter estimation (MCMC), and site-level decomposition
(the `shap_*` columns — a misnomer; see the v0.2 rename in the roadmap).

The package was migrated out of a paper-specific research project, and **much
capability is still unmigrated** — most notably the model's motion/mode arm (only
the structure × site cell of a larger divergence grid is implemented). The
`dev/plan.md` roadmap tracks what each future version adds.

It depends on the `penm` package (an `Imports:` dependency) for the ENM
machinery — `set_enm()`, mutation scans, response matrices. `penm` is consumed as
an installed library; its source is not part of this repo and is not edited here.

## How planning works here (two tiers)

The v0.1 plan was a single 575-line blueprint that tried to be both the durable
spec and the step-by-step execution script — so every execution surprise forced an
edit to a huge document. We deliberately split that:

- **`dev/plan.md`** — the **coarse, durable roadmap**: what each version *is* (goal
  + rough scope) and findings that must not be re-discovered (e.g. the
  precomputation property). High-level and stable; edited only when a version's
  *goal* changes. Read it for orientation.
- **`dev/PROGRESS.md`** — the checklist for the **version currently in flight**.
  Written from that version's detailed plan when the version starts; dormant
  between versions.
- **`dev/LOG.md`** — append-only history (sessions, decisions, dead ends). Holds
  the "why" behind the current state and the v0.1 detail the roadmap no longer
  carries.

**Detailed planning is done per-version, at execution time** — not up front. When a
version starts: enter plan mode, read the specific code that version touches, write
the executable detail then (and write `dev/PROGRESS.md` from it). Deep code
exploration belongs where it's actionable.

**Record decisions before coding.** This still holds — it now means *write the
per-version detailed plan before coding*, and touch the roadmap only if the
version's goal changed. Don't fold surprises into ad-hoc code; plan them first.

When this file and `dev/plan.md` disagree, the roadmap wins — update this file or
flag the contradiction.

## The migration source is READ-ONLY

Continuing migration reads from a **local frozen snapshot** of the source project
that lives inside this package:

```
tmp_src/
```

A one-time copy of the live project (`.../lab/projects/MSA-neutral-structure-evolution`),
made so the migration cannot touch the original — the live path is intentionally
**not** reachable from this session. `tmp_src/` is `.Rbuildignore`d. It is a
snapshot, so it goes stale: re-copy deliberately if the live source changes.
**Treat `tmp_src/` as read-only** — copy *out of* it, never edit *in* it.

## Migration rules (apply whenever migration continues)

This is a migration, not a refactor — stay close to the source. When copying a file
from `tmp_src/R/` into `R/`:

1. Remove `library(...)` calls (top-level and inline). All package usage goes
   through `Imports:` + `@importFrom` (or `pkg::fn()`).
2. Keep function signatures identical. Do not rename functions or variables or
   "improve" logic as part of a copy. (Deliberate, recorded API changes — like the
   v0.2 `phi_*` rename — are separate, planned work, not folded into a copy.)
3. Namespacing: **Option A** by default — package-level `@importFrom` + bare calls
   (smallest diff). **Exception:** a Suggested package (e.g. ggplot2) cannot be
   `@importFrom`'d — qualify calls as `ggplot2::...` behind a `requireNamespace`
   guard.
4. Keep `%>%` (magrittr, re-exported).
5. Strip the non-standard `@requires` roxygen tags (every source file has them;
   roxygen2 doesn't recognise them). Keep standard tags.
6. Add `@export` to public top-level functions; helpers stay unexported
   (`#' @noRd`).
7. Fix catalogued bugs as you reach each file — record them in the version's
   detailed plan (e.g. the tree route's `p_act`→`p_ma` bug, the `akima`→`interp`
   swap; see `dev/plan.md` findings).

No file renames are folded into migration copies (renaming is a later deliberate
refactor, decided when files are actually rearranged).

## House conventions

- **Naming:** `snake_case`, following the source.
- **Roxygen on every function.** `#' @export` for public API, `#' @noRd` for
  internal helpers. Document `@param`, `@return`/`@returns`, and add `@family` tags
  to group related functions. `@examples` wrapped in `\dontrun{}` when they need
  real data.
- **Imports live in one place:** `R/msamodel-package.R` carries the `"_PACKAGE"`
  doc and the `@importFrom` directives for the whole package. Re-export `%>%` via
  `#' @importFrom magrittr %>%`.
- **Tibbles** (not data.frames) for tabular returns; tidyverse for data
  manipulation.

## Tracking progress

- **`dev/PROGRESS.md`** — checklist for the in-flight version (bare substep titles,
  `[ ]`/`[x]`). Rewritten from that version's detailed plan when it starts.
- **`dev/LOG.md`** — append-only history, newest first. Reverts, scope changes,
  decisions, dead ends live here.

**After finishing any substep:** (a) tick it in `dev/PROGRESS.md`, and (b) add a
one-line dated entry to `dev/LOG.md`. Do this as the substep finishes — tied to
*finishing*, not to ending a session.

## Development commands

```bash
# Document (regenerate NAMESPACE + man/ from roxygen) — run after any roxygen change
Rscript -e "devtools::document()"

# Run all tests / a single test file
Rscript -e "devtools::test()"
Rscript -e "testthat::test_file('tests/testthat/test-msa-mcmc.R')"

# Full check (v0.1 baseline: 0 errors, 1 warning, 2 notes — all deliberately
# accepted, GitHub-only not CRAN; see dev/LOG.md)
Rscript -e "devtools::check()"

# Install locally
Rscript -e "devtools::install()"
```

Embedded `znb_*` datasets are (re)generated by running
`data-raw/prepare_znb_data.R`. The fixture is **frozen**: regenerate it
intentionally, not incidentally. `test-spm-generate.R` catches drift between the
SPM-generation code and the embedded fixture, so the parameters in the data-prep
script and that test must stay in sync.

NAMESPACE and everything under `man/` are roxygen-generated — edit the roxygen
comments and run `document()`, never hand-edit them.

## Working style (from project memory)

- Don't hand-wave uncertainty — flag anything unverified. Only say
  "confirmed/verified/fixed" when you actually checked.
- Own contradictions across turns instead of smoothing them over.
- Push back on bad ideas rather than silently implementing them.
- Terse is good; real uncertainty should be stated explicitly.
- **Tool choice:** read files with `Read` (use `offset`/`limit` for big files).
  Avoid `sed`/`head`/`cat`/`tail`/`awk`/`echo` in Bash for reading — they trigger
  permission prompts in this environment.
