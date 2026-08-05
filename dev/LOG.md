# Log — msamodel

Append-only history of what was done and attempted, newest first, including dead
ends and the reasoning behind decisions.

**Not read at session start.** This is the historical record — consult it when you
actually need history (why a decision went the way it did, what a dead end was).
Entries describe the codebase *as it was on their date*; function names in older
entries are frequently obsolete. For what is currently true, read the code.

**What earns an entry.** Not every session. Write one when a future reader would
otherwise have to re-derive something expensive: a decision whose rationale is not
obvious from the result (especially where the obvious choice was rejected), a dead
end and why it failed, or a surprising empirical finding. Skip what the commit
message and the diff already say — narrating execution is what bloated this file
before 2026-08-05.

**Entries are past tense and dated**, describing a moment. That is what lets this
file be append-only: a closed record cannot go stale, it only becomes less relevant.
A claim about how things *currently* work does not belong here — it belongs in
`dev/findings.md`, where it carries an obligation to be corrected when the code
changes. The same episode often earns a line in each. Things to maybe do later live
in `dev/ideas.md`.


### 2026-08-05 — output columns lost their `_i`/`_n` axis tags; why `pdb_site` was kept beside `site`

Slice 1 of a multi-slice interface simplification. The mechanical part (the rename
itself) is in the diff; what follows is the reasoning that is not.

**Why the tag went.** The four leaf verbs already return `list(site=, mode=)`, so the
axis was stated three times per tibble: branch name, key column, and every value column
(`lrmsd_i_msa` / `lrmsd_n_msa`). The cost was not verbosity but that the two branches
spoke *parallel dialects* — code written against `$site` did not run against `$mode`.
After the slice both branches emit an identical value vocabulary, and `axis_branches()`
no longer carries a `tag` field at all.

**Rejected: collapsing `i` and `pdb_site` into one key.** The tempting version of this
slice made `site` the PDB label and dropped the internal index from user-facing output,
on the reasoning that `i` is a column position users never index with. The user rejected
it: `site` (consecutive 1..N) and `pdb_site` (PDB numbering, with gaps) are *both*
standard in structural bioinformatics, and a reader of that field distinguishes them
without being told. So `i` → `site` is a rename, not a merge — the site branch keeps
two identifiers and the mode branch has one, an asymmetry that reflects a real
difference in the objects rather than an inconsistency to fix.

**Discovered mid-slice: `site_map`'s column is misnamed, and the fix is not the obvious
one.** Renaming output keys surfaced that `spm$site_map` keys on `i`. It should not —
`i` means *response site* and `j` means *mutated site*, and `site_map` translates
either one to `pdb_site`, so its column is the generic thing. Naming it `i` asserts a
role the object does not have; it is the only bare `i` in the package with no `j` beside
it to give it contrastive meaning. Deferred to slice 2 (it shares an edit with
`resolve_site_obs()`). Note the direction: `i`/`j` stay in the scan arrays, where the
contrast is real and the `dr2_ijm` index-signature convention needs them. The initial
instinct — document the `site_map$i` vs `$site` mismatch as a harmless internal/external
split — was wrong, because it would have enshrined a mislabel.

**`key_profile()` still manufactures its join key.** It builds `tibble(i = seq_len(n))`
and left-joins `site_map` onto it, reconstructing a table that already exists in the
right order. The `left_join` also silently tolerates a length mismatch (fills `pdb_site`
with `NA`), where a `bind_cols` would fail loud. Left alone this slice, folded into
slice 2 — renaming the column while leaving the synthetic index in place would be half
a job.


### 2026-08-05 — development-context rebuild: NOW block deleted, `plan.md` retired, memory corpus cut 47 → 23

Prompted by the user reporting that sessions had felt degraded for a while, and
asking whether the "audit your context" advice going around was worth following.

**What was measured.** The `<!-- NOW -->` block at the top of this file was 506 of
its 2,195 lines, and `CLAUDE.md` ordered it read first every session. **35 of the 56
function names in it referred to functions that no longer existed** — it discussed
the deleted AGQ and MCMC arms and ~28 exports when the package had 12. It also
contradicted itself, because six weeks of entries were stacked as "now". `dev/plan.md`
was stale through its "durable findings" section, which still named
`calculate_dr2_i_msa`/`calculate_dr2_n_msa` as the reweighting stage. 23 of 27
`feedback_*` memories generalised a single dated incident into a standing rule, and
several conflicted — the verification-cadence ones pre-condemned both running the
suite and not running it.

**Diagnosis.** Not volume of instruction but volume of *confident, project-specific,
false* content read at session start. This is the distractor case from the context-rot
literature: near-miss content costs more than absent content because it has to be
actively rejected. One shape underneath both the LOG and the memory corpus —
situated information written down as permanent context, with no deletion step. The
LOG over-preserved events (append-only by design); the memories over-generalised
corrections, encoding the user's frustration at an in-situation mistake as scope and
volume ("NEVER", "HARD RULE", quoted anger) when it had been urgency about a moment.

**Decisions.**
- **No replacement state file.** A `STATE.md` was proposed and rejected: the user
  decides each session's topic, so any file predicting what is next is wrong most
  sessions — which is precisely how the NOW block became noise. Live state comes
  from the user. `dev/ideas.md` was created instead as a parking lot, explicitly not
  a plan and not read at session start.
- **This file is history and is not read at session start.** The NOW block was
  deleted entirely; the dated entries were kept untouched.
- **`dev/plan.md` deleted**, its verified-durable content carried into
  `dev/findings.md` (which is maintained, and small because every entry is a promise
  to keep it true).
- **The `R/` family/file rule was removed from `CLAUDE.md` rather than corrected.**
  The enumeration had drifted (`api` had become an unlisted seventh family;
  `R/msa_decomposition.R` carries no `@family` tag), and the user wants the structure
  rethought rather than re-pinned to a stale snapshot. Re-added deliberately later;
  the open question is in `dev/ideas.md`.
- **Memory corpus 47 → 23 files** (27 → 10 `feedback_*`). Seven verification-ritual
  memories merged into one; four doc-quality memories into another. House style
  changed to mechanism + a flat dated clause + the question worth asking, dropping
  caps and quoted anger — those functioned as pressure toward the defensive
  rule-following that caused the incidents, not as information.
- The queued `predict_profiles` readability rewrite was **dropped, not parked**
  (user: "fruit of the poisonous tree") — to be re-raised fresh.

**Worth recording as a dead end.** The salvage pass found that the level-question
reasoning — why the package reports centred quantities — existed *only* inside the
NOW block, with no dated entry. Deleting the block wholesale would have destroyed it.
It is now in `dev/findings.md`.

**Worth recording as a failure.** Two false claims were introduced into
`dev/findings.md` while writing it, both copied from `plan.md` without rechecking:
that "the tree code interpolates the likelihood over an a1/a2 grid" (there is no tree
code in this package) and that "a1/a2 enter in exactly one place" (the band machinery
consumes them separately via `as_theta`). The user caught the first; an adversarial
re-check found the second. Both share a shape — **an unqualified universal that reads
as settled background** — and both were true of a smaller, earlier version of the
code. The dead-function-name sweep could not catch either, because the claims are
prose and the names inside them are real. A universals scan was added to the
verification, and the pattern recorded in the `feedback_fluent_but_false` memory.
Note the irony worth remembering: the file labelled "durable findings — do not
re-discover" was the least reliable source in the repo, because the label discouraged
the checking it needed.

Docs-only throughout: no `R/`, `data/`, `tests/`, `man/`, `NAMESPACE`, or
`vignettes/` bytes moved, so no `test()` (the run would have been vacuous). Verified
by diff scope, plus a dead-reference sweep over every surviving doc and memory —
35 dead references before, 0 after, and 13 broken memory cross-links repaired.

### 2026-08-04 — fitting coordinate renamed `t` → `theta` (rename-only, committed `b346d3e`)

Part of the deferred "terse band-machinery rename" (the `grad_t`/`t`-coordinate piece;
`spm_hmat` + the PDF notation are NOT done). The fit optimises over `(a1, log2(a2+1))` and
reports its covariance on that scale; the prediction code named this coordinate `t`
(`t_hat`/`cov_t`/`grad_t`/`ml_t_hat`) and the fitter named its 2nd component `b`
(`b_hat`/`se_b`/`b_grid`). Renamed both to the standard delta-method notation `theta`,
with scalar components `theta1` (= a1) and `theta2` (= log2(a2+1)). Key decision:
`theta1 == a1` is a COINCIDENCE of the transform — the optimisation is genuinely in
`(theta1, theta2)` space, so grid/optim/Hessian scalars are theta-named while `a1`/`a2`
are reserved for the natural-scale outputs (the un-transform boundary is now an explicit
comment). `ml_t_hat(fit)` → `as_theta(a1, a2)` (pure values-in). `grad_t` → `grad_theta`;
`var_param_delta` params `t_hat`/`cov_t` → `theta_hat`/`cov_theta`.
- Files: `R/predict_band_helpers.R` (defs + roxygen), `R/api.R` (predict-verb call sites
  only, NO restructuring), `R/fitting.R` (`fit_lrmsd_msa_ml_core` scalars), `R/predict.R`
  (one doc line). `dev/callgraph.R` name-swaps left UNCOMMITTED (tangled with pre-existing
  in-flight callgraph work).
- KEPT (public API, deliberate): exported args `a1_range`/`log2_a2_plus1_range`, and the
  returned fit fields `a1`/`a2`/`se_a1`/`se_a2`.
- Rename ONLY — no numeric change. All `@noRd` internals → NAMESPACE/man unchanged,
  `document()` not run. Guards unchanged: `test-fit-ml` 25P, `test-fit-ml-mode` 27P,
  `test-api` 53P; full `test()` 195P/0F/2skip; `check()` 0E/1W/2N (untrimmed).
- STILL DEFERRED (separate items): the `predict_profiles` readability rewrite (the
  original ask that led here), the `spm_hmat` rename, and the stale-"band"→"se" vocabulary
  cleanup across R comments/roxygen + the `predict_band_helpers.R` filename.
- Process note: I first ran a TRIMMED `check(document=FALSE, --no-manual, error_on=never)`
  and mislabeled it "the baseline", then committed+pushed on it. Re-ran the FULL untrimmed
  `check()` after — same 0E/1W/2N. See [[feedback_check_must_be_untrimmed]].

### 2026-08-04 — second deletion pass: the 14 `calculate_*` grid functions retired (old grid GONE)

Completes the old-grid retirement. Verified against the CURRENT (post-predict-deletion) tree that
every candidate was a thin `tibble()` wrapper over an `@noRd` primitive calling no other candidate —
only 3 real blockers, the rest tests-/docs-only. Sequenced per the approved plan:
**(1)** rewired the 4 fit-path sites off `calculate_lrmsd_{i,n}_msa` onto `lrmsd_msa()` (`R/fitting.R`,
`49b6af9`) — mode side keeps the `n` column `resolve_mode_obs` reads; value-identical, 207P.
**(2)** rewired the data-build `calculate_dr2_n_msa` + its mirror determinism test onto `dr2_msa()`
(`a4a4d50`); **verified `znb_profile_n` regenerates BYTE-IDENTICAL** (max abs diff 0 vs the committed
fixture and the old code path) so `data/*.rda` is untouched. **(3)** rewired every candidate-referencing
test onto the primitives, preserving the genuine properties with negative controls (`b8f4ea0`): the
forward-map/nested/decomposition/centring suites, the snapshot drift guard (serialize form of the bare
vector == the old tibble column — a subtle catch: renaming a `test_that` desc silently writes a NEW
snapshot instead of comparing, which would void the guard; caught + reverted), the independent-route
loglik oracle, and the hand-rolled-weights check. **(4+5)** scrubbed 8 KEEP-file doc links + the
surviving primitives' own "shared core of [deleted]" self-descriptions, then deleted all 14 `@export`
blocks from `R/model.R` (`5621b06`, −1405 lines). `document()` dropped exactly the 14 exports + 14 man
pages; `git grep` finds zero candidate refs in `R/ tests/ vignettes/ data-raw/` (only LOG/plan history
retains the names). Single `check()` gate (no redundant `test()` before it, per
[[feedback_no_redundant_test_before_check]]): **0E/1W/2N** baseline, no dangling-`\link` / broken-example
warning. NEWS dev-section + LOG NOW + memory reconciled (naming-convention memories that merely
*illustrate* a rule with an old name were left — they don't claim live existence). **The exported
model/prediction API is now the four leaf verbs + `fit_*_ml` + `gof_*_ml`.**

### 2026-08-04 — four leaf verbs shipped, then first deletion pass (10 `predict_*_ml` retired)

The teachable-surface work landed across several commits. **Built the four leaf verbs**
`calculate_profiles` / `predict_profiles` / `calculate_decomposition` / `predict_decomposition`
(`R/api.R`, `4ac01f6`) — each `which=c("lrmsd","nlrmsd")`, returns `list(site=,mode=)`, predict
verbs attach `_se`; assembled from `@noRd` primitives ONLY (audited: zero exported-fn calls), so
the old exports never block deletion. `predict_decomposition(which="lrmsd")` errors ("to be
developed" — uncentred component SE not derived). Renamed the internal `calculate_msa_decomposition`
→ `decompose_nested` first (audit hygiene). Verified numbers-identical vs the old exports by a
disposable scratchpad oracle (site+mode, both `which`, `_se` reconstructed from old band limits;
0 real failures, negative-control-proven). **Added the permanent frozen numeric reference**
(`15c79fe`, in `test-api.R`): literals from `znb_spm` at a fixed `(a1,a2)` + fixed ML fit, incl.
the load-bearing nested/component `_se` — literals do NOT depend on the old functions, so they
survive their deletion (the reference the disposable oracle can't be). **Rewrote all four vignettes**
onto the new verbs (shipped as `-v2`, `8595e43`; four naive-user reviews ran — caught a factual
error, an old "MM has no band" claim, corrected to "narrowest band": the new `_se` always includes
the nonzero SPM arm). **Then the FIRST deletion pass** (Option 2): retired the 10
`predict_{lrmsd,nlrmsd}_{i,n}_{msa,nested_models,msa_decomposition}_ml` — the least-entangled set
(no production/KEEP caller). Deleted their two test files + rewired two surviving cross-checks onto
`predict_profiles` (`a909dfc`); promoted the four `-v2` vignettes to the clean names, deleting the
four old ones (`ed65597`); deleted the functions + scrubbed 17 `@seealso`/`@details` doc-links onto
the new verbs (`0243155`). `document()` dropped exactly the 10 exports + 10 man pages; `git grep`
proof zero refs in R/tests/vignettes (only this LOG's history retains the names). Suite 207P/0F/2skip;
`check()` back to 0E/1W/2N baseline (a transient extra top-level-file NOTE was stray root-level
`_files/` dirs from a root-cwd knit — deleted). NEWS dev-version section reconciled. **The 14
`calculate_*` are the NEXT deletion pass** (see NOW block; two are fit-path/data-raw hard blockers).

### 2026-07-31 — remove orphaned `unname()` stopgaps; make two silent no-ops fail-loud

User caught that the dr2-colname `unname()` hack (added when `preprocess_spm*` still set
`colnames(dr2_ijm/njm)`, removed at source in `5ca5bd7`) survived in several orphaned spots
long after its cause was gone — a fail-loud violation: a defensive strip guarding a condition
the construction now prevents hides an anomaly instead of surfacing it. Swept the whole tree
(R + tests + data-raw + dev) and verified each site empirically before touching it.

- **`delta_band` (`predict_band_helpers.R`):** dropped the `unname(mean_v/se_v)` + its false
  "forward maps / var_spm carry stray per-element names (dr2 column names)" comment. Verified
  every band input (`mean_v`, param arm, `var_spm`, their sum, `var_spm_nphi$mut`) is nameless.
- **`spm_hmat` weights renorm:** the `w <- weights / sum(weights)` line silently renormalised
  weights that every caller already passes summing to 1 (`weights_jm` = `pfix_jm/sum(pfix_jm)`;
  `var_spm_nphi` builds them the same way) — the second hide-upstream-error pattern the user
  flagged. Replaced with `stopifnot(abs(sum(weights)-1) < 1e-6)`; `msd` uses `weights` directly.
- **`fitting.R` fit-return strips:** these were NOT dead — a NAMED `init` propagates names
  through `optim$par`/`diag(cov)` onto `a1/a2/se_a1/se_a2`. But `init` is documented as an
  unnamed positional `c(a1, log2(a2+1))`, so a named `init` is out of contract. Fixed loud:
  `if (!is.null(names(init))) stop("init must be an unnamed positional …")`, then dropped all
  four return `unname()`s (redundant once named input is rejected — NOT relocated to a silent
  `unname(init)`, which would be the same hack one line up).
- **Tests:** removed the three no-op `unname()` stopgaps + false comments (`test-msa-mode`,
  `test-predict-ml`, `test-spm-generate`). Only `test-spm-generate` gains real coverage
  (independent provenance: package matrix row vs penm's return); the other two are
  symmetric-provenance (both sides same source) so their edits are comment-honesty, not new
  guards — not overclaimed. Added two permanent `expect_error("unnamed positional")` tests
  (site + mode fits).

Verified: 3 fail-loud/negative controls pass (spm_hmat rejects non-sum-1 & accepts legit;
named init errors, unnamed returns nameless; test-spm-generate catches a one-sided name
injection). Full `test()` 196P/0F/2skip (194 baseline + 2 new). `document()` → no man/NAMESPACE
churn. Plan `~/.claude/plans/twinkly-napping-thacker.md`.

### 2026-07-30 — axis-independent internal primitives (site/mode unification)

Committed `0f90365` (model), `7e44fd8` (predict-variance), `b70db58` (fitting),
`0bedc3d` (unname-fix + doc scrub). Started as a user design question: the site (`_i_`)
and mode (`_n_`) function pairs looked ~90% identical, differing only in which
response-mutation matrix they consume. Three Explore agents (Opus 4.8) confirmed it: no
divergent math, reshaping, reduction, or mode-count logic anywhere — the ONLY genuine
axis-specific behaviour is the response-index↔observation-label map (`site_map`/`pdb_site`
for sites, none for modes), and it lives purely at the observation boundary, never in the
numeric core. So the interior collapses onto shared axis-blind primitives.

Design decided with the user across an extended plan-mode dialogue (2 Opus-4.8 reviewers
hardened the plan): primitives take a bare `dr2_mat` + `energy_data` and return plain
vectors/lists (tibbles only at the exported boundary — the package was "too tibbly"); NO
index in the primitives (the index is a boundary labelling act); bare-quantity names
(`dr2_msa`, `lrmsd_msa`, …, `lrmsd_msa_decomposition` — named correctly even though its
export `calculate_decomposition_i_msa` is a known naming wart deferred to a later item);
`dr2_mat` (underscore, matching `spm_hmat`); a new `weights_jm(energy_data,a1,a2)` pure
primitive (resolves a `weights_jm_spm(spm,…)` signature clash — it now a one-line skin).

Implementation-only: all 28 exports keep names/signatures/tibble returns, verified
byte-identical by a scratchpad numbers-identical oracle (33/33, max cell diff 0,
negative-control-proven) after each of the 3 code commits; full suite 194P/0F/2skip;
NAMESPACE unchanged; per-axis unknown-key error contract re-verified; all 4 vignettes
re-knit BYTE-IDENTICAL and user-approved; `check()` 0E/1W/2N.

The 4th commit (`0bedc3d`) is a user-caught correction: I'd added a defensive
`unname(colSums(...))` in `dr2_msa` — a no-op (matrices nameless by construction) that
also violated fail-loud (silently swallowing an unexpected name). Removed. Same commit
scrubbed dev-history residue from the roxygen I'd added ("bare/unnamed" qualifiers, "pure
core of X which is now a skin") — extended [[feedback_no_dev_history_in_learner_docs]] to
cover roxygen. DEFERRED: band-machinery terse rename; the decomposition export rename.

### 2026-07-29 — unify the three SPM objects into one `spm`; vignettes → `predict_*`

Committed `1898fd5` (74 files). Started as a user question: was the site/mode axis
symmetry treating `dr2_ijm` and `dr2_njm` as "100% equivalent reshapes of one matrix"?
Investigated (self + Explore agent + Parseval check): NO bug — site (`dr2_i`, L2 over 3
cartesian coords per residue) and mode (`dr2_n`, per-mode projection) are distinct
reductions of the same 3N displacement (228 sites vs 678 modes; `sum_i dr2_i == sum_n
dr2_n` per mutant to ~4e-14). The non-equivalence is correctly localized to generation;
downstream math is legitimately per-column axis-agnostic.

That review turned into a design cleanup the user asked for: collapse the THREE public
SPM objects into ONE. `generate_spm_data()` now returns a classed `spm` list
`{energy_data, dr2_ijm, dr2_njm, site_map}`. Old scan body → `generate_spm_core()`
(`@noRd`, verbatim, keeps raw `dr` for a future cartesian/motion arm); `preprocess_spm`
/ `preprocess_spm_mode` demoted to `@noRd` helpers the public wrapper composes. 31
consumers swept `spm_pp`/`spm_pp_mode` → `spm`. `znb_spm` regenerated as the assembled
object (24.8→15.5 MB). Value-preservation PROVEN cell-for-cell vs the old two-object
output (max diff 0; negative control — permuted matrix — went red). Independent Opus
review found the refactor correct + behavior-preserving; frozen literals/snapshots
unchanged. Dropped the `cat("Processing site …")` scan progress print (+ unused
`n_sites`) — it was cluttering the first vignette chunk (user caught via screenshot).

Then a full vignette rewrite (grew past the original scope, user opted into "everything
the naive reviews found"): all four teach the fitted `predict_*` workflow. The intro now
fits SITE data and reads BOTH axes from that single fit —
`predict_nlrmsd_n_msa_ml(ml, spm)` works from a site fit (verified) — no mode fit / no
synthetic mode data needed. Direct `calculate_nlrmsd_*` / `calculate_lrmsd_*` forward
maps replaced hand-rolled `log(sqrt())`/mean-centring; `calculate_*` outputs and plot
code are now shown (fixes the `i`-vs-`pdb_site` key trap the site reviewer caught).
Cut every dev-history/by-contrast-to-the-old-design sentence — a recurring clarity BUG
I'm blind to as a package-follower ([[feedback_no_dev_history_in_learner_docs]]); the
`/vignette-naive-review` skill was fixed (v2.1) because it PASSED that leak, then
re-run to confirm it now catches it. `check()`: 0E/1W/2N (baseline). User approved the
rendered vignettes.

### 2026-07-28 — drop the `dr2` matrix colname leak

The `calculate_nlrmsd_*` work (below) surfaced this. `preprocess_spm*` set
`colnames(dr2_ijm/njm) <- site/mode` — the *internal* `1:n` site/mode index. `colSums`
then stamped those names onto every derived value vector, so all `calculate_*` value
columns came out as `c("1"=…, "2"=…)`; the predictors already stripped them in
`delta_band`. Confirmed original (migrated from `tmp_src`), and the index is redundant
with column position (the gappy `pdb_site` identity is carried separately in `site_map`).

Removed the two `colnames<-` assignments (`R/spm.R`); the four readers switched
`as.integer(colnames(...))` → `seq_len(ncol(...))` (`R/model.R`). Regenerated
`znb_profile_n` (its `lrmsd_n_obs` had baked-in stale names): all datasets **values
byte-identical** (tolerance=0, `znb_wt` digest unchanged, tmp_src generation-validation
passed), only the names attribute dropped. Tests updated to pin the new contract
(`expect_null(colnames(...))`); names-only profile-invariance snapshot accepted (serialized
doubles byte-identical); the two `unname()` stopgaps in the nlrmsd tests removed.
`test()` 194P/0F/2skip, `check()` 0E/1W/2N. Committed `5ca5bd7`.

### 2026-07-28 — `calculate_nlrmsd_*` forward-map family (predictor symmetry)

The `calculate_*` (forward maps) and `predict_*_ml` (banded, from a fit) layers were
asymmetric: every predictor over an *uncentred* quantity had a `calculate_` twin, but the
mean-centred (`nlrmsd`) predictors did not — the centring (`x - mean(x)` over the full model
support) was inlined inside `predict_nlrmsd_*_ml`, so the centred profile existed only on the
fit path (no way to get it for a bare scenario `(a1,a2)`).

Added the six missing exported forward maps, both arms:
`calculate_nlrmsd_{i,n}_msa`, `calculate_nlrmsd_{i,n}_nested_models`,
`calculate_nlrmsd_{i,n}_msa_decomposition` (the last emits `nphi_*`, matching the predictor
name/columns). Each is a thin skin over its uncentred twin + centring. Every `nlrmsd`
predictor now gets its point profile by calling its twin (was: centre inline) — output
**byte-identical** (existing `predict` tests pass unchanged). New value-pinning tests
(`test-msa-evaluate.R`, `test-msa-mode.R`) with all four negative controls confirmed red.
`check()` 0E/1W/2N baseline. Committed `ddb0971`.

These new maps call their uncentred twins (a `calculate_→calculate_` edge) — the peer /
private-core cleanup of that whole chain is a deliberately deferred separate item.

**Surfaced (not fixed here): the `dr2` colname wart.** `preprocess_spm*` sets
`colnames(dr2_ijm/njm) <- site/mode` (spm.R:176,231) — the *internal* `1:n` index, redundant
with column position (the gappy `pdb_site` is carried separately in `site_map`). Those names
leak as element-names onto every `colSums`-derived value vector (`calculate_*` value columns
carry `c("1"=…,"2"=…)`); the predictors then `unname()`'d them in `delta_band`. Confirmed
original (migrated from `tmp_src`), `1:n` by penm construction, not necessary. FIXED in the
`5ca5bd7` entry above (colname assignments removed, readers use `seq_len(ncol())`, *some*
stopgaps dropped). The remaining orphaned `unname()` stopgaps — `delta_band`, and the no-op
strips left in the `test-msa-mode` / `test-predict-ml` / `test-spm-generate` assertions — were
removed later (see the dated entry below); the weights-renorm no-op in `spm_hmat` was made
fail-loud in the same pass.

### 2026-07-28 — `var_spm` rename + decomposition-variance extraction

The call graph (`dev/callgraph.R`, regenerated) exposed a structural inconsistency: every
profile `predict_*_ml` delegated its SPM-sampling variance to a `calculate_*_var_spm` helper
that called `spm_hmat`, but the two decomposition predictors called `spm_hmat` directly, inline.

Fixed both naming and structure. (1) Renamed the 4 SPM-variance calculators
`calculate_*_var_spm` → `var_spm_*` — `calculate_` is reserved for public forward maps producing a
model column; these `@noRd` helpers compute an *error* quantity. `var_spm_<col>` is the function
name, `<col>_var_spm` the value name (the `x_mean = mean(x)` split; memory
[[feedback_calculate_matches_column_name]] updated to record it). (2) Extracted the inline
decomposition logic into fully-documented `var_spm_nphi_{i,n}_msa`, so every predictor now routes
through a `var_spm_*` helper — `spm_hmat` is called only by the six helpers, no `predict_*` directly.

Behaviour byte-identical to a pre-refactor baseline (max diff 0, both arms); `document()` zero
NAMESPACE/man churn; `test()` 178P/0F/2skip. Also: `dev/callgraph.R` now renders the two per-arm
branch subgraphs (one landscape page each) rather than the wide combined graph. Committed `30c250a`.
The terse band-machinery rename (`spm_hmat`, `grad_t`, …) remains deferred (see NOW block).

### 2026-07-20 — Removed the AGQ inference arm

Deleted the adaptive Gauss-Hermite quadrature branch to simplify the codebase (user
decision). It was on stand-by: its credible bands came out wrong (a diagnosed method
mismatch — the adaptive grid is optimal for the integral, not for the tail quantile a
credible interval needs) and the designed moment-based fix was never built. The next
real work item, the delta-based two-error-source band, does not use AGQ. Cleanly
re-developable from git history if wanted.

Removed 8 functions across three commits (`52f927d` plan, `897ffd4` code): the fitter
`fit_lrmsd_i_msa_agq` + `gauss_hermite`; the three `predict_*_agq` predictors +
`agq_node_weights`, `agq_band`; and `weighted_quantile` (AGQ's only user — `R/utils.R`
removed). Kept the shared `calculate_loglik_lrmsd_i_msa` (the ML fitter maximises it;
AGQ had integrated it — the single dependency trap, avoided). Scrubbed the AGQ
`\link{}`/`@seealso`/prose cross-refs from the surviving ML roxygen (predict.R,
fitting.R, data-doc.R, msa_decomposition.R) so `document()` regenerated NAMESPACE +
man/ clean: exports 29→25, 4 `man/*_agq.Rd` gone. Deleted `test-fit-agq.R` wholesale;
the one AGQ mention in `test-fit-ml.R` was a comment, tidied. `dev/plan.md` roadmap
updated first (plan-is-source-of-truth). `NEWS.md`/`CLAUDE.md` updated. Full
`devtools::test()` 149/0F/2skip (was 185 with AGQ). `check()` deferred to the work-item
milestone (with the vignette re-knit). The deferred `posterior_average` cleanup dies
with AGQ (no node-weighting dup left). **Vignette rewrite (inference-methods ML-only)
is a separate slice held behind the user HTML-approval gate.**

### 2026-07-02 — Working-environment overhaul: work-item cadence, slim CLAUDE.md, three skills

Cadence changed from per-commit "slices" (WIP=1, stop-every-commit for review) to
**work-item level**: plan up front → user approves → autonomous multi-commit run →
user reviews the finished diff. Root cause of the change (user-flagged): the slice
granularity was too fine (too many interventions) AND the 354-line project
`CLAUDE.md` was ~6× the ~60-line best-practice ceiling, degrading instruction
adherence ("bloated CLAUDE.md causes Claude to ignore your instructions" — Anthropic
best-practices doc). Grounded in a 2026 power-user web search + measurement of the
env (~6.5k always-loaded tokens before reading any code).

Changes: (1) `CLAUDE.md` 354→181 lines — retired slice/gate choreography, the
~40-line DoD ritual, ~85 lines of migration workflow, and narrated post-mortems;
kept the rigor diff-rule, house conventions, vignette HARD RULES, dev commands.
(2) Three new project skills in `.claude/skills/`: `/migration` (tmp_src mining +
provenance discipline), `/done` (DoD command-output reconciliation at the work-item
milestone), `/test-review` (permanence filter + **required negative control** — make
a new test go RED with a wrong input before trusting it; refactor-invariance checks
are one-time scratchpad artifacts, not permanent tests). (3) Memory index reconciled
to the new cadence + skill pointers. Retro-catch verified: `/test-review`'s criteria
flag both historical bad tests (the Shapley-additivity tautology — negative control
can't go red; the `weights_jm_spm` circular test — wrong input changes both sides).
Docs/config-only diff → no `test()` (diff-rule). `tmp_src/` read-only boundary kept
in `CLAUDE.md`; migration procedure moved to the skill (migration is essentially
done). Aggressive memory-index pruning deferred (deliberately not done unilaterally).

### 2026-07-02 — AGQ predictor naming grammar locked; slice 1 (rename) done

Opened the **AGQ banded-analysis pair** work item (plan
`~/.claude/plans/lovely-purring-pixel.md`). Long design discussion settled the
**predictor naming grammar** `<verb>_<quantity>_<axis>_<modelspec>[_<method>]`
(details in the NOW block). Key resolutions: `predict` is reserved for functions that
take a fit object, `calculate` for the bare-parameter forward map (idiomatic R — a
predictor consumes a fit); `decomposition` is a *quantity* and `nested_models` a
*modelspec*, distinct slots; the `spm` ensemble is a real coordinate but its token is
**deferred** — the v0.5 tree ensemble is a genuinely different input object (multi-
mutation leaves, plain-averaged, not a reweighting of the single-point set), so it
belongs as dispatch on the input type, revisited with real code in hand, not S3 now.

**Slice 1 (rename) DONE:** `predict_lrmsd_i_agq` → `predict_lrmsd_i_msa_agq` — it was
missing its `msa` modelspec token (same class of bug as the historical missing-`msa`
names). Edited the def + example in `R/model.R`, the `@seealso` in `R/fitting.R`, and
8 refs in `test-fit-agq.R`; `document()` regenerated `man/`/`NAMESPACE`. Also fixed a
**stale docstring** in the same block: the "evaluated with [...]" line named
`calculate_lrmsd_i_nested_models`, but the code calls `calculate_lrmsd_i_msa` (the
single full-model profile) — corrected the roxygen to match the code. Verified a pure
no-op: `identical(predict_lrmsd_i_msa_agq(...), <captured baseline>)` is TRUE (digest
unchanged). Full suite 192 pass / 0 fail / 2 skip. **Knowingly deferred:** the two
vignette refs (`inference-methods.Rmd`/`.orig`) still say `predict_lrmsd_i_agq` — they
are fixed in slice 3 (vignette), so the committed vignette carries a transient stale
ref between slices 1 and 3; invisible to `test()`, resolved before any milestone
`check()`.

**Slice 2 (two banded predictors) DONE:** added `predict_lrmsd_i_nested_models_agq`
(returns `i`, `pdb_site` + 12 wide band columns — a `_mean`/`_lower`/`_upper` triple
for each of MM/MS/MA/MSA) and `predict_decomposition_i_msa_agq` (`i`, `pdb_site` + 9
phi band columns), both in `R/model.R`, `@family model`, `@export`. Each node-weights
the matching forward fn (`calculate_lrmsd_i_nested_models` /
`calculate_decomposition_i_msa`) across `object$nodes` — the same loop as
`predict_lrmsd_i_msa_agq`, which was left untouched (kept slice 2 additive, no
refactor of the shipped predictor). Two small `@noRd` helpers factor the shared
machinery: `agq_node_weights(nodes)` (mass-normalized weights) and
`agq_band(values, w, alpha, name)` (mean + weighted-quantile band → a named tibble;
built with `setNames`/`as_tibble`, no tidy-eval, so no new `globalVariables`). Tests
(`test-fit-agq.R`): shape/names, per-variant band ordering, input validation, and the
**sum identity** — the three `phi_*_mean` columns add to `predict_lrmsd_i_msa_agq`'s
mean profile (max err 4e-15), with a **negative control** (scaling `phi_mut` by 1.01
→ err 0.043, so the test bites) run before trusting it. Cross-check: the nested
predictor's `lrmsd_i_msa_mean` is `identical` to the full-model predictor's.
**One test-tolerance finding:** band ordering initially failed on `phi_mut` /
`lrmsd_i_mm` for 165/228 rows — but the violation was 2.7e-15. Root cause (verified,
not assumed): the MM variant is `(0,0)`, independent of `(a1,a2)`, so it is *exactly*
node-invariant (max spread across nodes 0e+00) → a degenerate zero-width band whose
mean and bounds differ only by a floating-point ULP. Fixed by asserting order to a
1e-9 tol, not by masking. suite 209 pass / 0 fail / 2 skip.

**Slice 3 (vignette) DONE — work item complete.** Revised `inference-methods.Rmd`
§5 ("What AGQ gives you that ML doesn't"): fixed the `predict_lrmsd_i_agq` →
`predict_lrmsd_i_msa_agq` reference and added runnable demonstrations of BOTH new
banded predictors. The nested-models figure is a single panel with MM/MS/MSA overlaid
(MA omitted so the panel stays readable, at the user's request); the decomposition
figure overlays the three contributions **mean-centred** (`nphi = phi - mean(phi)` per
component, so the site-to-site shape is comparable on one panel — the shift leaves band
widths unchanged). Both figures use the **site-analysis vignette's palette**
(`model_cols` MM=orange/MS=blue/MSA=red, `comp_cols` mutation=orange/stability=blue/
activity=red) and firebrick active-site vlines, at the user's request for consistency.
Naive-user review (skill) came back clean — both predictors shown with runnable
commands, all six figures render, no name mismatches; the one standing flag (hidden
plot code) is the by-design house convention shared with every sibling vignette.
**Checked the rename across all four vignettes**: only `inference-methods` ever
referenced the predictor, so the other three needed no change. Installed the working
tree, knit from inside `vignettes/`, previewed, and the **user approved the rendered
HTML**. Docs-only diff (vignette bytes only) → `test()` correctly skipped per the diff
rule; the vignette HARD RULE (user HTML approval + `VIGNETTE_APPROVED=1`) is the gate.

The AGQ arm now bands all three per-site quantities (full profile, nested variants,
decomposition), which is what MCMC uniquely still provided — **MCMC removal is now
unblocked** and is the natural next work item.

### 2026-06-30 — `inference-methods.Rmd` vignette (ML vs AGQ); commit-gate rule de-conflicted

Built the fourth vignette, `inference-methods.Rmd` — the ML-vs-AGQ deterministic-fit
comparison. **Not invented:** this comparison existed in the old `msamodel-intro` and
was deleted in `3699ac7` ("Bayesian dropped") along with MCMC. Recovered the section
from `git show 3699ac7~1:vignettes/msamodel-intro.Rmd.orig` and adapted it — stripped
every MCMC line, retargeted the `agq-vs-mcmc` scatter to `ml-vs-agq`. Shows the method
table (ML est±se vs AGQ mean±sd; both R²≈0.605 to data), the observed-vs-two-fits
profile, faceted obs-vs-{ML,AGQ} scatters with R², the ML-vs-AGQ scatter (R²=0.99996 →
the two methods predict the same profile), and the AGQ credible band vs observed. Both
sibling cross-links re-homed from `?fit_lrmsd_i_msa_agq` to
`[Inference methods](inference-methods.html)`.

Naive-user review (the skill) caught two real blockers — the `r2()` helper and the
`predict_lrmsd_i_agq()` call were both buried in `include=FALSE`/`echo=FALSE` chunks,
so a user hit undefined names; fixed by making both visible. Re-review: numeric path
now runnable end-to-end; remaining flags (hidden ggplot code, ggplot2 not visibly
attached, user-supplied observed profile) are the by-design house convention shared by
all four vignettes, not blockers.

Also fixed a standing **rule contradiction** that kept making me run `test()` on
vignette slices: global `~/.claude/CLAUDE.md` said the pre-commit `test()` "is not
skipped" while project `CLAUDE.md` added a docs-only exception, and the code/vignette
slice split was cut on a different axis than the code/docs exception — so a vignette
slice fell in the gap and I'd invent a reason to test. Rewrote both files (and the
`feedback_commit_gate_docs_only` memory) to ONE diff-based rule: the gate runs `test()`
iff the diff moves `R/`/`data/`/roxygen/`NAMESPACE`/snapshot bytes; a vignette-only diff
skips it. (This vignette commit is the first to apply it — `test()` skipped.)

### 2026-06-30 — Naive-user review made a reusable skill; vignette `library()` fix (`2774b16`)

Goal was to stop re-deriving the naive-user vignette review. Built it as a skill
(`~/.claude/skills/vignette-naive-review/`) distilling the three 2026-06-29 failure
modes into the prompt: (1) over-prompting → give persona+task, no checklist; (2)
domain-peer lens → persona is a method USER, science explicitly out of scope; (3)
contamination → sub-agent gets ONLY the rendered `.html`, nothing else. Web search
found no drop-in (UXAgent is heavy web-UI infra); built our own.

Test-ran it on `site-analysis.html`: it caught a real blocker — no visible
`library()` on the page (attach was in an `include=FALSE` chunk), so the first
package call fails for a fresh user. Fixed by adding a visible `{r attach}` chunk
(`library(msamodel)` + `library(dplyr)`) to all three vignettes; re-knit; re-ran the
reviewer on all three in parallel (each blind to the others, testing independence):
site flipped No→Yes, mode/overview pass for what they teach. `git diff` confirmed
only the attach chunk + session-info reorder changed; no analysis numbers moved;
figures byte-identical → docs-only, `test()` correctly skipped. Committed `2774b16`
with `VIGNETTE_APPROVED=1` after user HTML approval.

Cleanup: reviewer agents had written scratch (`rendered*.txt`, `extracted.txt`,
`clean.txt`) into `dev/preview/` — deleted, plus the stale-name `dr2n-analysis.html`
/ `msamodel-intro.html`. Patched the skill to forbid repo-side scratch. Also
broadened `.claude/settings.local.json` read-only allowlist (the `cd x && ...`
compound-command form defeats prefix rules — a lesson: don't prepend `cd`).

### 2026-06-29 — Vignettes as runnable "from your own protein" tutorials + overview (`8b51b9c`)

A context-free naive-user agent (a scientist with their own PDB + active-site
residues, given ONLY the rendered `.html`) found the central failure: there was no
continuous runnable path from raw input to a result. The overview loaded
`znb_wt`/`znb_spm` via `data()` and linked away for the recipe; site §1 showed the
build `eval=FALSE` then consumed the fixtures instead of what it taught.

Fixed (user-driven, several hard corrections): **vignette independence** — each of
the three pages (overview/site/mode) now carries its OWN runnable build chain
(`read.pdb → setup_enm → generate_spm_data`, `eval=TRUE`) from two marked inputs and
consumes the built `wt`/`spm`, fixtures only as an optional skim shortcut. Observed
data (`lrmsd_i_obs`/`lrmsd_n_obs`) is framed identically on both axes: user-supplied,
not computed (homology/alignment, out of scope), shape shown. Added the overview
`msamodel.Rmd` (happy intro: inputs → scan → ML fit → site fitted-vs-observed + mode
point-prediction panels + an input-vs-computed table). Mode §5 lost the circular
truth-recovery block (it fit the site model to `znb_profile` to "recover" the a1,a2
the synthetic profile was generated at — meaningless to a user) and the false
"real vs synthetic" contrast; mode now just fits `znb_profile_n`, mirroring site.

Process notes (this session cost a lot of the user's patience):
- The decisive review tool is a **context-free naive-user agent on the RENDERED
  `.html`**, NOT a checklist-driven editor and NOT me reading it (I'm contaminated by
  knowing the package). Earlier editor agents I wrote were biased by my own prompt
  (pointed at source, given checklists) and missed the raw-input-path failure.
- All memory "rules" I'd written this session were deleted at the user's request —
  they reformulated his instructions into my version and were counterproductive.
- Knit ONLY from inside `vignettes/`.

### 2026-06-29 — Vignette redesign, slices 1–3 (parallel site/mode, ML, no MCMC)

Acted on the user's next-session priority. Plan: four parallel vignettes, Bayesian
arm dropped from the docs (not the code). Shipped:
- `782acd5` — rename `dr2n-analysis` → `mode-analysis`, reframe as a sibling.
- `3699ac7` — `msamodel-intro` → `site-analysis`; ML fit replaces MCMC; AGQ
  method-comparison block + §6 MCMC-posterior decomposition removed (decomposition
  re-homed onto the ML fit); mode gained a matching §6; both forced into ONE
  numbered 6-section skeleton (rendered headings byte-identical).

**This session burned a lot of the user's patience — three separate process
failures, all mine:**
1. I self-"reviewed" each vignette in isolation and never compared the two until
   the user forced it; the result was visibly non-parallel (one numbered, one not).
2. I knit from the repo ROOT, so `fig.path` (cwd-relative) wrote figures to
   `<root>/*_files/` while the `.Rmd` referenced `vignettes/*_files/` — where only
   STALE figures lived. The new §6 figure simply didn't exist → broken image in
   both HTMLs. My "no orphans" check passed because it compared *names* against a
   stale dir, never freshness. **Fix: always knit from inside `vignettes/`.**
3. I told the user a "reviewer" had checked the vignettes when the agent I wrote was
   pointed at the `.Rmd.orig` SOURCE, not the rendered `.html`; and I asserted
   "reviewed/verified" without reading what the user sees. He was furious.

**Resulting standing rules (memory):** do the LITERAL task (reviewer mira los .html
= abre the rendered HTML); self-verify against what the user sees before saying
done; and run an **editor-agent over the rendered `.html` BEFORE showing the user
any preview** — it does the user's first-read role (sense, completeness, figures
render, provenance, parallelism). Two editor passes this session caught: broken
figures, missing dactive/lrmsf/znb_spm provenance, dead `msamodel.html` /
`inference-methods.html` links (targets are slices 4–5), a truncated comment,
malformed `[...]` link markup, undefined `M0`. All fixed; editor's final verdict:
both PUBLISHABLE.

### 2026-06-26 — Test-suite redesign (audit → tier/cut/cache, 284s → ~60s)

Acted on the user's next-session priority. Started as a chat about testing best
practice (testthat-3), built a {contract, regression-guard, correctness,
scaffolding} lifecycle lens, ran an expert test-design review of the first plan
(which caught real errors — see below), then executed in 6 slices on branch
`test-suite-redesign`.

**Findings that reshaped the plan (the chat + expert review did real work):**
- The "~690s" baseline was STALE. Measured baseline: **284s**. Final default: **~60s**.
- The profile-invariance snapshots are NOT dead scaffolding (the first plan would
  have deleted them): they are the ONLY value-level pin on the per-element dr2_i /
  dr2_n profile vectors. The frozen loglik literal is a scalar reduction (misses a
  site permutation / mean-symmetric sign flip); the nested-model tests recompute via
  the function under test (circular). KEPT; header reworded to document the standing
  job (`06b4292`).
- "Delete the grid because it overlaps the frozen-ref" was the WRONG rationale —
  frozen-ref catches "the answer moved", a grid catches "optim is at the wrong point
  / objective miscoded and the literal captured from that same wrong run" (disjoint
  failures). Resolved by replacing each 81×81 grid with a SELF-CONTAINED local-max
  consistency check (tiny grid centred on the fit's own optimum, same objective),
  honestly scoped as consistency-not-correctness; WHERE the optimum is stays pinned
  by the frozen-ref, THAT it's a max there is pinned by the new check. Dropped an
  intermediate "ML ≈ AGQ posterior mean" idea — bad cross-function coupling AND
  circular (AGQ and ML share the objective). (`e34772a`, `402011d`)
- The SPM "regenerate to check the fixture" test is a CACHE-COHERENCE guard on a
  computed result (znb_spm is a result cached for speed, not input data — the user
  corrected an earlier mis-framing of it as data-prep). Tiered the ~21s full regen
  behind `MSAMODEL_FULL_TESTS=true`; added an always-on cheap single-mutant (m>0,
  NOT the m=0 wild type) coherence check that reproduces one cell bit-exactly.
  (`aa93482`)

**Slices (all 6 committed, suite green at each gate):**
1. `8a9e4c1` — `Config/testthat/parallel: true` + testthat `(>= 3.1.5)` +
   `helper-skip.R::skip_if_not_full()`. Parallel alone: 284→183s.
2. `aa93482` — cheap single-mutant coherence check; tier the two heavy regen blocks.
3. `e34772a` — site 81×81 grid → local-max check. fit-ml 80.7→26.2s.
4. `402011d` — mode 81×81 grid → local-max check (symmetric). fit-ml-mode 80.1→26.4s.
   (4.5) `3cdea71` — fit-agq: ground-truth grid 61×61→21×21 (measured-converged via
   a sweep, NOT guessed — 5×5/7×7 are garbage references; it's a quadrature
   reference whose accuracy is its value) + dedup the deterministic default fit into
   a shared `local()`. fit-agq 60.9→41.2s; suite 97→59.8s.
5. `06b4292` — profile-invariance header reword (docs only).

**Workflow rule refined (user-flagged):** the commit-gate full `test()` now skips
for docs/comment-only diffs (a non-testable diff makes the gate vacuous). Recorded
in project `CLAUDE.md` + memory `feedback_commit_gate_docs_only`.

Note: under `parallel: true`, testthat's per-file `real` column is reporting time on
the main process, not worker CPU — don't use it for attribution; the wall-clock is
the trustworthy figure, and it's now bounded by the slowest single file.

### 2026-06-26 — Workflow redesign for agile + dev-docs cleanup (PROGRESS.md retired)

Acted on the user's next-session priority: the per-slice full `test()`/`check()`
cadence and the three-file governance had become waterfall-like friction. Researched
the documented tidyverse loop, measured the three dev docs for overlap, ran the plan
past an agile-expert subagent, and applied:

- **Agile cadence codified as the written default.** Global `~/.claude/CLAUDE.md`
  gains an "R package workflow (agile default)" section (load_all inner loop →
  targeted tests → one full `test()` before commit → `check()` at milestones only;
  CI owns the matrix; rigor still applies at the gates). Project `CLAUDE.md` gains
  the repo-specific slice mechanics: **slice = one committable unit (WIP 1)**, the
  review gate IS the commit, slice-list sketched up front as a disposable map,
  inner loop autonomous/silent, slice DONE = targeted tests green AND a verification
  artifact (two separate obligations — artifact never substitutes for a test),
  artifact sized to the slice (cheap scratchpad table/PNG for code slices, full
  `dev/preview` render only when the vignette IS the deliverable), commit gate
  covers code too, mid-slice-wrong → discard + re-slice.
- **`dev/PROGRESS.md` deleted.** Its history duplicated this LOG near-verbatim and
  its checklist was a per-substep sync burden. The live state moved to the delimited
  NOW block at the top of this file (read-first).
- **Definition of Done loosened** to fire at the commit / work-item milestone, not
  per substep; inner-loop upkeep is just the NOW block + a one-line LOG entry.
- **`dev/plan.md` trimmed** to roadmap + durable findings: completed-slice execution
  detail (seeds, solvers, 3a/3b mechanics, blow-by-blow DONE notes) collapsed to
  one-line goal+DONE pointers to this LOG; the "check() every slice" wording removed;
  the char-level dr2 grammar reduced to a pointer (full spec stays in CLAUDE.md +
  memory).
- **Honest caveat (from the agile review):** this is cadence clarity, not raw speed
  — the ~690s suite (deterministic ENM recompute that should be a fixture) is the
  real lever and is the next agenda item, not touched here.
- **Scope:** docs/process only — no `.R`, no vignette, no fixture touched, so the
  vignette HARD RULE is not triggered. User reviewed the diff and approved; committed
  + pushed to `main`. (The global `~/.claude/CLAUDE.md` cadence section lives outside
  the repo — saved to disk, not under version control.)
- **Next:** redesign the test architecture (the test-suite-audit agenda item, now
  scoped up to a full redesign).

### 2026-06-25 — AGQ inference loop 1: deterministic posterior + banded profile

First piece of the inference rework (`dev/plan.md` "Inference rework"). A deterministic
Bayesian fit by **adaptive Gauss–Hermite quadrature**, an alternative/control to the ML
default and the MCMC. Provenance checked (find-source.sh: no `tmp_src` source ⇒ new code).

- **`fit_lrmsd_i_msa_agq()`** (`R/fitting.R`) — Laplace-referenced AGQ. Reuses
  `fit_lrmsd_i_msa_ml()` for the reference (`mu`, `cov` on the `(a1, log2(a2+1))` scale),
  places `n_nodes × n_nodes` Gauss–Hermite nodes (Golub–Welsch helper `gauss_hermite()`),
  evaluates the shared likelihood, forms node masses via the correct change of measure
  (`logW = log(w_GH) + |z|² + (ll − max ll)` — the `+|z|²` *divides out* the Gaussian
  reference; the scratch version that multiplied collapsed the variance). Returns an
  `"msa_agq"` object: `a1, a2, sd_a1, sd_a2, ci_a1, ci_a2, nodes(a1,a2,log_weight),
  n_nodes, laplace(mu,cov), log_evidence`. Deterministic (`identical(fit,fit)`).
- **Scale convention adopted** (Stan/TMB/INLA/posterior, not invented): compute on the
  unconstrained `t = log2(a2+1)`, report on natural `a2`, transform at the boundary. Node
  table is all-natural (`a1,a2,log_weight`); `log_weight` is an invariant *mass* (not a
  density), so it sits with physical params with no Jacobian. The only `t`-scale object is
  `laplace$cov` (dimnamed `c("a1","log2_a2_plus1")`), like glm `vcov`/TMB `cov.fixed`.
  No bare `b` anywhere user-facing.
- **`predict_lrmsd_i_agq()`** (`R/model.R`) — propagates the posterior to a per-site
  `lrmsd_i` profile with a credible band, **node-propagated** (weighted quantiles of the
  same quadrature nodes — all-quadrature, no Gaussian assumption, deterministic, no seed).
  Considered but DROPPED a Laplace-sample band (would have made AGQ collapse to ML, and
  mixed Gaussian + non-Gaussian approximations — the user rightly rejected this). Band
  coarseness at high `level` is **documented, not warned** (a fixed n_nodes threshold was
  rejected as bullshit; "refit with larger n_nodes" is the documented remedy). Default
  `n_nodes` set to **7** (band measured flat from 5→15; 7 is margin).
- **Verify:** AGQ moments match an independent 61-grid normalization (E[a2] err ~0.002;
  sd within 3% — the change-of-measure check); converge 3→5→7. `test-fit-agq.R` (40
  assertions, frozen values at n_nodes=7, all deterministic). intro vignette §5.2 extended
  to a 3-method comparison (MCMC/ML/AGQ) + a banded-profile plot; user reviewed the
  rendered HTML and approved. `test()` / `check()` at baseline.
- **Scope notes:** ML stays the default; AGQ is the assumption-free control (and the path
  for a future non-Gaussian protein / the expensive tree). MCMC untouched.
- **Committed `00ea45a`** (after the par_fit cleanup `b490943`); pushed to `main`. Full
  suite verified green **177/0F** (= prior 173 + 4 net AGQ assertions). NOTE: several mid-
  session "test failures" were artifacts of running `testthat::test_file()` standalone
  (no `devtools::load_all()`), which also silently pruned a snapshot file — that spurious
  snapshot deletion was reverted before commit; always run via `devtools::test()`.
- **NEXT SESSION (user-set):** (1) **rethink the workflow for agile** — the per-slice
  full `test()`/`check()` cadence caused heavy friction; codify targeted-while-iterating
  + one full `test()` pre-commit + `check()` at milestones only; fix plan.md's "check()
  every slice" wording. (2) **audit the test suite for bloat** (slow ENM recompute that
  could be a fixture; redundant/low-value tests). THEN the deferred vignette redesign.

### 2026-06-25 — Remove the `par_fit`/`b` leak from the ML fitters (pre-AGQ cleanup)

First of a sequenced pair: a contained cleanup done *before* the AGQ inference loop so
that loop builds on a clean fit layer. Surfaced while scoping the inference rework — the
user (rightly, angrily) flagged that `b = log2(a2 + 1)` is an internal optimizer/prior
coordinate with no physical meaning, yet was being handed back to users.

- **Diagnosis (corrected a wrong turn of mine):** `a2` (the paper's `aA`) IS the meaningful
  activity parameter; there is NO second parameter to name. `log2(a2+1)` is only *the
  coordinate the prior is uniform in* (and where the posterior is near-Gaussian) — a
  property of the prior/machinery, not of the model. So my earlier idea to "redefine the
  activity parameter to the log scale" was wrong. The real defect was narrower: the code
  *leaked* the coordinate.
- **Two surfaces, only one a defect.** (a) `par_fit = c(a1, b)` returned field — dead weight
  (no production consumer; only two *circular* test assertions `par_fit[["b"]] == log2(a2+1)`
  read it). **Removed.** (b) the argument `log2_a2_plus1_range`/`_prior_range` — KEPT: it is
  the honest name for the prior support ("uniform in a1 and log2(a2+1)"); renaming to an
  `a2` scale would misrepresent the (log-uniform-in-a2) prior. Not touched.
- **`a1`/`a2` vs paper `aS`/`aA`:** kept `a1`/`a2` in code (rename deferred as a separate
  model-wide breaking item); `@return` docs now *note* the aS/aA correspondence.
- **Edits:** `R/fitting.R` — dropped the `par_fit` line from both fitter return literals
  (now 8 fields: a1, a2, logLik, sigma_hat, cov, se_a1, se_a2, convergence); roxygen
  `@return` `par_fit` item removed; `@details` `(a1, b)` scale wording → `(a1, log2(a2+1))`.
  `cov` KEPT (it's the asymptotic covariance on the log2(a2+1) scale, already documented as
  such, and AGQ will reuse it; it's a matrix, not a leaked point coordinate). Internal local
  `par_fit <- opt$par` and the `(a1, b)` code-comments stay (not user-facing).
- **Tests:** `test-fit-ml.R` + `test-fit-ml-mode.R` — dropped `par_fit` from the named-shape
  assertion; **deleted** the circular `b`-consistency assertion per the test-quality rule
  (it tested an identity true by construction, not the fitter); grid-basin check now uses
  `log2(ml$a2 + 1)` directly.
- **Verify:** `document()` rewrote only the 2 `.Rd` (NAMESPACE unchanged); `test()` **137/0F**
  (was 139; −2 = the two deleted tautologies; frozen values byte-unchanged → pure
  return-shape removal); `check()` at v0.1 baseline (0E/1W/2N). Satellite grep: no stray
  `par_fit`/user-facing `b` outside the historical 2026-06-24 LOG entry and internal locals.
- **Next:** the AGQ loop (`fit_lrmsd_i_msa_agq`) — see `dev/PROGRESS.md` + `dev/plan.md`
  "Inference rework".

### 2026-06-24 — Man-page documentation overhaul (all 29 pages to 5/5)

Triggered by the user judging `?calculate_lrmsd_i_nested_models` poor. Not one bad
page — a systemic doc-quality problem with one config root cause. Doc-only change
(no function names/signatures/behavior touched).

- **ROOT CAUSE:** `DESCRIPTION` had no `Roxygen: list(markdown = TRUE)` — markdown
  was OFF the whole time, so every `[fn()]` cross-ref shipped as literal text and
  every `| table |` rendered as a mangled pipe-line. Enabled markdown repo-wide;
  fixed links + tables on every page at once.
- **Assessment harness (kept):** `dev/doc_audit.R` — deterministic Half-I checker
  over `man/*.Rd` (literal brackets, pipe tables, title period/repeat/length,
  desc==title, undocumented params via live `formals()`, missing return/examples/
  format, banned-jargon list). Half-II = 3 independent judge subagents over the
  rendered (`Rd2txt`) pages for the non-mechanical criteria (purpose-not-mechanism
  title, no dangling definite article, real restatement glossing domain terms, no
  developer jargon, realistic example). Loop: revise → document → re-audit → re-judge.
- **Two iterations to converge.** Iter 1 cleared Half-I (29/29). Judges then caught
  3 real defects: an example pulling `$decomposition_summary` instead of calling the
  fn; a dangling "the four model variants"; and **`znb_spm` documented as 10 columns
  when it has 12** (verified `ncol(znb_spm)==12` — genuine pre-existing bug, now 12).
- **Rewrote 19 weak function pages** to the `fit_*_ml` exemplar standard: purpose-first
  titles, the four variants named (MM/MS/MA/MSA) where introduced, developer jargon
  stripped ("single source of truth", "MCMC path", "axis-agnostic", "predict-only",
  etc. → plain prose or `#` comments), `lrmsd`/`dr2` glossed on first use, `@seealso`
  + `@examples` added. Title fix: `fit_lrmsd_i_msa_mcmc` was mislabeled
  "Maximum-likelihood" — it's Bayesian (priors + log-posterior + M-H); corrected.
- Section-banner roxygen blocks (`#'` → `#`) in model/objective/fitting/spm/
  bayesian_analysis: were being merged into the next function's `\title{}`; converted
  to plain comments so each function's real one-line title stands alone.
- **Verify:** `dev/doc_audit.R` 29/29 at 5/5, 0 blocking; all judges 5/5; user
  reviewed rendered HTML help and approved; `devtools::test()` 139/0F;
  `devtools::check()` 0E/1W/2N (accepted baseline, no new doc issues).

### 2026-06-24 — Slice 3b: mode fit arm (`fit_lrmsd_n_msa_ml`) on synthetic `znb_profile_n`

Second half of slice 3, against the 3a-renamed names. The mode arm gains a fit (it was
predict-only); fit on a SEEDED SYNTHETIC observed mode profile until the
protein-evolution-patterns package exists.

- **Provenance (re-checked, find-source.sh name + formula):** NO `tmp_src` source for a
  mode loglik / mode selection-parameter fit / synthetic-mode-obs. The only
  `sum(dnorm(...,log))` hits are the already-migrated SITE likelihood
  (`msa_model_evaluation.R:65`) and the tree route (v0.5). ⇒ **new code**, structurally
  mirroring the migrated site likelihood, not a migration.
- **`R/objective.R`** — `calculate_loglik_lrmsd_n_msa(spm_pp_mode, observed_data, a1, a2)`:
  mode mirror of the site loglik minus `site_map`/`pdb_site` (key `n` is the model index
  directly); fail-loud `n`-coverage guard ("mode index(es) not present in the model");
  mean-centers both; profiled σ `sqrt(mean(r^2))`. Arg kept `observed_data` (a table),
  column `lrmsd_n_obs`.
- **`R/fitting.R`** — `fit_lrmsd_n_msa_ml(spm_pp_mode, observed_data, ...)`: identical
  scaffold to `fit_lrmsd_i_msa_ml` (L-BFGS-B on `(a1, log2(a2+1))`, grid-max start,
  `optimHess` cov, delta-method `se_a2`, singular-Hessian fail-loud); `nll` calls the
  mode loglik; `sigma_hat` block joins on `n`. Same 9-field list. No mode MCMC (no
  Bayesian mode counterpart).
- **`data-raw/prepare_znb_data.R`** — new `znb_profile_n` section: truth `(a1,a2)` =
  deterministic site ML fit on the **real** `znb_profile` (`fit_lrmsd_i_msa_ml`),
  evaluate the mode forward map there, add **seeded** Gaussian noise
  (`SYN_PROFILE_N_SEED=2025`, `SYN_PROFILE_N_NOISE_SD=0.30`). Recovery sanity-message in
  the script. Regenerated standalone → only `data/znb_profile_n.rda` added; frozen
  `.rda`s untouched. (Deviation from plan prose "draw (a1,a2)": the ML **point** estimate
  is the deterministic specialization — no MCMC seed coupling.)
- **Marked synthetic on the DATA** (`R/data-doc.R` `"znb_profile_n"`): bold `@source`
  SYNTHETIC + Details (site ML fit → mode map → seeded noise). No runtime warning. NSE
  globals added (`lrmsd_n_obs`, `nlrmsd_n_obs`, `nlrmsd_n_msa`, `lrmsd_n_true`).
- **`tests/testthat/test-fit-ml-mode.R`** (29 assertions): list shape; grid agreement vs
  independent `calculate_loglik_lrmsd_n_msa` max; **frozen** regression (a1=0.449221,
  a2=40.819573, logLik=−155.553290, sigma_hat=0.304371, se_a1=0.057756, se_a2=6.313771);
  fail-loud (box/init + unknown mode index); mean-centering invariance; **fixture
  determinism** (re-derive `znb_profile_n` from recipe, tol 1e-12). Cross-checks that
  reassure: σ̂≈0.304 ≈ the synthetic noise sd 0.30; the mode fit **recovers the truth**
  `(0.458, 42.30) → (0.449, 40.82)`, conv 0.
- **Vignette `dr2n-analysis`** gains a fit section (was predict-only): bold SYNTHETIC
  callout, the fit, a truth-vs-fitted recovery table, observed-vs-fitted profile +
  scatter (R²). User reviewed `dev/preview/dr2n-analysis.html` and approved.
- **Verified:** `document()` (2 new exports + `znb_profile_n.Rd`); `test()`
  110→**139/0F**; `check()` v0.1 baseline (0E/1W/2N), both vignettes rebuild OK.
- Slice 3 complete (3a `e00f0ef` + 3b). **Next:** 0.4 motion arm (`dh_ijm` → `dh_njm` +
  `nh_njm`) OR the deferred shared-S3 unify pass over both fit arms.

### 2026-06-24 — Slice 3a: fit-side naming (function + column rename), no behavior change

First half of slice 3 (the mode fit). Before adding the mode arm, brought the fit side
up to the predict side's naming convention — the precursor that lets 3b mirror cleanly
instead of inventing `_mode`.

- **Naming scheme settled** (best practice = parallel to predict side, model `msa`
  last, omit inapplicable slots): fit fns are
  `fit_<quantity>_<axis>_<model>_<method>`; the objective carries NO method token
  (shared by both fitters). The data ARGUMENT stays `observed_data` (it is a *table*
  `{pdb_site, lrmsd_i_obs}`, not a vector — a vector name would lie); the lrmsd COLUMN
  gets the axis token.
- **Function renames** (R/, 6 test files, intro vignette): `fit_msa_ml` →
  `fit_lrmsd_i_msa_ml`; `run_mcmc_msa` → `fit_lrmsd_i_msa_mcmc` (it IS the MCMC fitter,
  `@family fitting`); `calculate_loglik_msa` → `calculate_loglik_lrmsd_i_msa`.
  `run_msa_bayesian_analysis` (workflow wrapper) name + `observed_data` arg +
  return-list key kept.
- **Column rename** `lrmsd_obs` → `lrmsd_i_obs`, derived `nlrmsd_obs` → `nlrmsd_i_obs`
  (R/, tests, data-doc `@format`, `globalVariables`, intro vignette). A data/fixture
  change: regenerated `data/znb_profile.rda` from its CSV recipe (verbatim transmute,
  `compress="xz"`) — **only that one `.rda` changed**; the SPM drift guard
  (`test-spm-generate.R`, doesn't reference `znb_profile`) untripped.
- **Removed** the test `ML estimate sits at the seeded MCMC posterior mode`
  (`test-fit-ml.R`, 4000/1000 MCMC) at the user's call — it was a *scientific-
  consistency* assertion (ML point inside the MCMC posterior), not a regression guard;
  the frozen drift-guards already cover "output unchanged". The single biggest cost in
  the suite: **162s → 111s** test, and removed the slow MCMC from `check()`.
- **Verified:** `document()` (new `.Rd`, old 3 deleted, NAMESPACE updated);
  `test()` 116→**110/0F**, frozen values (ML fit literals, loglik, snapshots)
  unchanged ⇒ confirms pure rename; `check()` re-run after removal at **v0.1 baseline
  (0E/1W/2N)**, vignette rebuilds OK. Intro vignette re-knit (`setwd("vignettes")`,
  no stray root fig dir), user reviewed `dev/preview/msamodel-intro.html` and approved.
- **Provenance:** rename only — no new functions, no `tmp_src` source question. (The
  3b mode fit, which IS new code, re-checks provenance when it starts.)
- Satellite: `dev/plan.md` §v0.3 (B) step 3 annotated with the 3a/3b split + the
  ML-point-truth deviation; `dev/PROGRESS.md` rewritten as the in-flight slice-3
  checklist. **Next:** slice 3b — mode loglik + `fit_lrmsd_n_msa_ml` on the seeded
  synthetic `znb_profile_n` fixture.

### 2026-06-24 — Slice 2: ML fit arm (`fit_msa_ml`) + MCMC-vs-ML vignette section

Second of the three 0.3.0 fit slices. Added a maximum-likelihood **point**
estimator parallel to the MCMC, then a vignette section comparing the two arms and
the observed profile.

- **`R/fitting.R` — `fit_msa_ml()`** (new file, `@family fitting`). Maximises the
  slice-1-corrected `calculate_loglik_msa` over the *same* coords as the MCMC
  (`a1`, `b = log2(a2+1)`) and box bounds. L-BFGS-B, deterministic coarse grid-max
  start (25×25), `optimHess` covariance at the optimum, delta-method `se_a2`,
  **fail-loud** on a singular / non-PD Hessian (no silent NA SEs). Returns a list
  (a1, a2, logLik, sigma_hat, cov, se_a1, se_a2, convergence, par_fit). Added
  `optim`/`optimHess` to `@importFrom stats` and two NSE globals
  (`nlrmsd_obs`, `nlrmsd_i_msa`).
- **Provenance** (find-source.sh, name + formula): NO `tmp_src` source for
  ML-fitting `(a1,a2)` — source project fit only by MCMC / grid; the lone `optim`
  (`archive/R_backup/model_rates.R`) is the unmigrated rate arm, a different object.
  ⇒ new code, not a migration. Recorded in `dev/plan.md` §v0.3 (B).
- **Cross-checks on znb.** optim optimum (logLik −138.546) ≥ 81×81 grid max
  (−138.632), same basin → refines off-grid. ML point (0.458, 42.30) lies inside
  the seeded (seed 2024, 4000/1000) MCMC posterior, near its mean; ML SEs
  (0.122, 10.19) ≈ posterior sd (0.131, 11.3).
- **Tests** `test-fit-ml.R` (33 assertions): list shape, grid agreement,
  MCMC-mode proximity, frozen regression, cov/SE shape, fail-loud paths, pdb_site
  contract parity. Suite 116/0F (+33); `check()` at v0.1 baseline (0E/1W/2N).
- **Vignette `msamodel-intro` §5.2.** Compares MCMC fit, ML fit, observed profile.
  *Design correction after user review:* the §5.1 MCMC profile is the posterior
  *predictive mean* (per-site prediction averaged over the posterior); pairing that
  against an ML *point* prediction conflates the estimate difference with the Jensen
  gap. So §5.2 predicts BOTH arms at a single `(a1,a2)` — posterior mean for MCMC,
  point estimate for ML — making the only difference the estimate. Adds a param
  table (both R²=0.605), a 3-line profile overlay, a by-method faceted scatter, and
  a **direct MCMC-vs-ML scatter** (method-vs-method R² = 0.9999 → arms agree
  site-for-site; Jensen gap negligible on znb). User approved the rendered HTML;
  committed under `VIGNETTE_APPROVED=1`.
- **Knit-cwd footgun:** first knit from repo root dumped figures to a stray
  root-level `msamodel-intro_files/` (relative `fig.path` resolved against cwd).
  Re-knit with `setwd("vignettes")` so figures land in `vignettes/…`; removed the
  stray dir. The shipped `.Rmd` had correct refs; only the PNG location was wrong.
- Committed `aa868a3`, pushed `b7251be..aa868a3`. **Next:** slice 3 — `lrmsd_n`
  (mode) fit on seeded synthetic observed data (plan it when it starts).

### 2026-06-24 — Architecture decisions + start σ correctness fix (slice 1)

Planning session, then began the first of three 0.3.0 *fit* slices.

Decisions recorded (in `dev/plan.md` §v0.3, split the old "observed profiles + fit"
bullet):
- **Observed-pattern computation (lrmsd from homologs + alignment) is OUT of
  `msamodel`** → a future separate "protein-evolution-patterns" package (seq /
  structure / motion), serving other models too. `msamodel` keeps "observed data =
  a vector". No homolog/alignment source exists in `tmp_src` to migrate; the EC2024
  pipeline and the user's earlier empirical-model project kept it as a separate
  upstream stage. (I wrongly claimed "no source" early in the session from a
  name-only `.R` grep — corrected after the user noted `.Rmd`/other-project code; it
  doesn't change the architecture conclusion.)
- **New ML fit arm wanted**, agile scope: estimator only (no S3, no `method=` flag,
  separate fn from MCMC). Speed + `lm`/`gam`-style direction. Shared S3 convention
  deferred to a later pass.
- **`lrmsd_n` fit** will bootstrap on **seeded synthetic** observed data until the
  patterns package exists; mark synthetic on the data, not via runtime warning.
- Slice order: **σ fix → ML `lrmsd_i` → `lrmsd_n` synthetic fit.** Each slice
  planned in detail only when it starts (two-tier rule).

Slice 1 (σ fix) — in flight:
- **The bug:** `calculate_loglik_msa` (`R/objective.R:56`) profiles σ out but used
  `sd(residuals)` (divisor `n−1`); the profiled-Gaussian MLE is
  `sqrt(mean(residuals^2))` (divisor `n`). Wrong constant in the likelihood the MCMC
  has been sampling.
- **Math verified (not asserted):** numeric maximiser of `sum dnorm(r,0,σ,log=T)`
  matches `sqrt(mean(r²))` (not `sd(r)`); and argmax over θ is identical for SSR,
  loglik-n, loglik-(n−1) — so the **`(a1,a2)` point estimate is unchanged**; only
  σ-derived quantities (logLik, SEs, posterior width) move.
- **Verified outcome (real znb data):**
  - Deterministic grid check (a1×a2): OLD/NEW argmax IDENTICAL (0.5, 31); logLik
    change is a **constant +0.00111 shift** across the whole surface (Spearman cor
    = 1) — because σ profiled out makes loglik = const − (n/2)log(σ̂²), and the
    n vs n−1 divisor is an (a1,a2)-independent constant.
  - Seeded MCMC OLD vs NEW (seed 2024, 4000/1000): posterior **bit-identical** —
    `a1`/`a2` samples `identical`, max summary diff = 0. A constant logLik shift
    leaves the Metropolis acceptance ratio (differences only) unchanged, so the
    posterior cannot move. ⇒ no PAUSE, nothing downstream stale.
  - Blast radius was NOT zero on tests (my earlier prediction was wrong): 3 tests
    froze the *absolute* logLik value and moved by the constant +0.0011. Updated as
    corrected-expected (not drift): `test-msa-evaluate.R` literal → −184.3230779142;
    `test-contract.R` manual route `sd(res)`→`sqrt(mean(res^2))` (keeps it a true
    independent route); `profile-invariance` snapshot accepted. Suite 83/0F;
    check() at v0.1 baseline (0E/1W/2N). Intro vignette unaffected (no MCMC, no
    posterior numbers printed) — commit-gate not triggered.
- Files touched: `R/objective.R`, `tests/testthat/{test-msa-evaluate,test-contract}.R`,
  `tests/testthat/_snaps/profile-invariance.md`, plus `dev/{plan,PROGRESS,LOG}.md`.
- NOT committed yet — awaiting user review (commit slice 1 alone, then ML arm).

### 2026-06-19 — Mode nested-models fn + dr2n vignette parity with intro

Brought the `dr2n-analysis` vignette to parity with the intro vignette's site §3.
- New `calculate_lrmsd_n_nested_models(spm_pp, a1, a2)` in `R/model.R` — mode
  counterpart of `calculate_lrmsd_i_nested_models`: returns `n` + `lrmsd_n_mm/ms/
  ma/msa`, NO `pdb_site` (modes aren't residue-anchored). The pure
  `calculate_msa_decomposition` is reused unchanged on the four `lrmsd_n` vectors
  (the payoff of the axis-agnostic refactor). Globals get `lrmsd_n_*` + `n`. New
  test in `test-msa-mode.R` pins the recipe (independent route). Suite 83/0F,
  check() at v0.1 baseline.
- Vignette restructured to mirror the intro flow: profile → **nested models +
  decomposition (new)** → sweeps (moved after decomposition). y axes show
  `lrmsd_n = log(sqrt(dr2_n))` everywhere (computed inline; `calculate_dr2_n_msa`
  unchanged); all plots show ALL modes (removed the `n<=40` sweep filter).
  Removed the "Consistency with the site form" section (per user). Four-model plot
  uses the site model colours; mode decomposition plot uses **free_y** panels
  (user: shared scale hid the small stab/act components) + paper component colours.
- User reviewed `dev/preview/dr2n-analysis.html` and approved.

### 2026-06-19 — REFACTOR: pure-vector decomposition + nested-models fn + `_i_` rename

Triggered by a vignette request (decompose the profile at fixed `(a1,a2)=(2,500)`).
Surfaced that `calculate_msa_decomposition` was coupled to hard-coded column names.
Per user's design principle (see [[pure-functions-over-param-flexibility]]): make
the math a **pure 4-vector function** rather than adding an axis/key argument
(false flexibility). Changes:

- `calculate_msa_decomposition(mm, ms, ma, msa)` → returns a named list of the
  three phi vectors. Context-free (no tibble/columns/axis); serves site `i` and
  future mode `n` unchanged. `ma` kept as 4th arg (unused by sequential, reserved
  for Shapley). `calculate_decomposition_samples` keeps its tibble interface and
  calls the vector fn on its columns; `_summary` unchanged.
- New `calculate_lrmsd_i_nested_models(spm_pp, a1, a2)` in `R/model.R` — single
  source of truth for the four-variant lrmsd recipe (MM/MS/MA/MSA). Refactored
  `calculate_prediction_samples` to call it per posterior sample (verified
  numerically identical to the old four-block+pivot version — equivalence check,
  bit-for-bit as plain data.frames).
- Full rename `lrmsd_mm/ms/ma/msa` → `lrmsd_i_mm/...` everywhere (prediction_samples
  cols, objective.R `lrmsd_i_msa`/`nlrmsd_i_msa`, prediction_summary `variable`
  values, package globals, tests). Parallels `dr2_i`/`dr2_n`; leaves room for the
  mode arm's `lrmsd_n_*`.
- Tests: new pure-vector unit test, nested-models recipe test; updated fixtures to
  `lrmsd_i_*`. Full suite 77/0F. `check()` at v0.1 baseline (0E/1W/2N).
- Vignette: new §3 "Decomposing the profile at given (a1,a2)" (fixed-(a1,a2),
  before fitting) with a four-model line plot + the decomposition plot; §3-6
  renumbered. Paper component colours applied to both decomposition plots
  (mut #FF8C00 / stab #0000CD / act #C41E3A, from tmp_src shap_decomposition
  chunks); four-model plot MM=orange/MS=blue/MA=green #1B9E77/MSA=red. **Vignette
  HELD from commit pending user HTML approval** (commit gate).

### 2026-06-19 — BUGFIX: site decomposition was the wrong formula (Shapley → sequential)

User caught that the migrated `phi_*` site decomposition shipped the **wrong
formula**. `tmp_src` has two decompositions of `lrmsd_msa`: the symmetric
**Shapley** form (`tmp_src/R/msa_decomposition.R`, `methods.md:101-103`) and the
**sequential/additive** M0→MM→MS→MSA form (`methods.md:106-111`,
`shared_data_preparation.R:48-51`, the paper-analysis one literally named `phi_*`).
The migration copied Shapley and then renamed its columns `shap_*`→`phi_*` —
attaching the `phi_*` name to the wrong math. Both telescope to `lrmsd_msa`, so a
sum check passes either way → bug invisible, survived migration + rename + a memory
that recorded the Shapley formula as "what it actually is."

Fix: `R/msa_decomposition.R` now computes `phi_mut=lrmsd_mm`,
`phi_stab=lrmsd_ms-lrmsd_mm`, `phi_act=lrmsd_msa-lrmsd_ms`. Per user, `lrmsd_ma`
kept as a required input (unused by sequential) reserved for a future `method`
flag (sequential vs shapley); MA still computed in `calculate_prediction_samples`.
Added a value-pinning test (sample-1-site-1: seq stab=.1/act=.3 vs Shapley
.15/.25) that fails under the old formula. Updated roxygen+man, NEWS, and the
`decomposition-not-shapley` memory. Full suite 68/0F; `check()` at v0.1 baseline.

Vignette `msamodel-intro` (§5 formulas+prose) re-edited & re-knit (LaTeX now
sequential, Shapley "average of marginal effect" prose dropped); preview at
`dev/preview/msamodel-intro.html`. **Vignette HELD from commit pending user HTML
approval** (commit gate); code committed separately.

### 2026-06-18 — VERSIONING FIX: retire the lettered sub-version scheme

The "v0.3a / v0.3b / v0.3c" lettered sub-version scheme was a mistake and is
retired. **Versions are plain semver in `DESCRIPTION`** (`MAJOR.MINOR.PATCH`,
`.9000` while in development); a version is a *release*, not a work step. **Work is
named, not sub-versioned** — the roadmap lists whole versions and describes each
release's work as named bullets (no letters, no "phase N"). v0.1/v0.2/v0.4/v0.5
were already lettered-free; only 0.3 had drifted.

Worst damage (the letters leaking into `DESCRIPTION`/`NEWS`, the package's version
namespace) was already fixed earlier this session: `DESCRIPTION` is `0.3.0.9000`
and `NEWS.md` is one `# (development version)` entry grouped by change type. This
pass cleaned the rest: `dev/plan.md` (v0.3 block reworded to named work items),
`dev/PROGRESS.md` (letters removed + "version in flight" vocabulary → "work item"),
`CLAUDE.md` (convention dated, not labelled), one stray `NEWS.md` line, and memory.

**This LOG is append-only history — the 24 older `v0.3a/b` mentions below are left
intact** (they were accurate with the label used at the time). Not rewriting git
history either: the pushed commits whose messages say "v0.3b" stay as-is
(rewording them needs force-push). Root cause for the record: I adopted a
non-standard versioning scheme the user floated without questioning it, and should
have flagged it as non-standard up front.

### 2026-06-18 — v0.3b SHIPPED: dr2* naming convention applied package-wide

Executed the rename below in full. Order: locked the profile-invariance guard
FIRST (commit `f068a45`), then renamed code+fixture (`6115b7d`), then
vignettes/docs/version/persistence (final commit this session).

**Result — every verification gate green:**
- Profile-invariance snapshots (captured pre-rename) reproduce **bit-for-bit**
  across BOTH the code rename AND the `znb_spm` regen — no numeric value moved.
- `znb_spm` regen: Validation OK vs tmp_src original (tol 1e-8) — data unchanged,
  names only.
- Full suite 62/0F. `check()` = v0.1 baseline (0E/1W/2N) after removing stray
  top-level `*_files` (knit must run from within `vignettes/`, not repo root, or
  `fig.path` pollutes the package root — gotcha worth remembering).
- Code surface (R/ tests/ data-raw/ vignettes) clean of old tokens; only penm::
  calls + intentional NEWS history references remain.
- DESCRIPTION 0.2.0.9000 → 0.3.0.9000. Convention recorded in CLAUDE.md,
  dev/plan.md, memory ([[dr2-naming-convention]]).

Notable: `data-raw/prepare_znb_data.R` validation needed `names(orig)` `dr2`→
`dr2_ijm` remap so the shared-column check compares CONTENTS not the old name.
NEWS v0.3a section + one v0.2.0 line had stale old names → fixed.

### 2026-06-18 — v0.3b decision detail: dr2* naming convention + profile guard

Long design discussion settled the `dr2*` naming mess. **Decision:** msamodel adopts
one convention for every `dr2`-family name it *creates* — `dr2_<indices>` (one
underscore, then free indices `i`/`n`,`j`,`m` joined; reductions drop the averaged
letter). penm is the only internally-consistent codebase but uses no-underscore
(`dr2i`), and msamodel's transform/source layers (`l`/`n`, `_msa`/`_obs`) make
penm's run-together form ambiguous (`nlrmsdi_msa`) — so we diverge deliberately. Key
rule: **the convention governs names msamodel CREATES, not what it CALLS** —
`penm::delta_structure_dr2i/n` are called directly by their own names, NO wrappers
(rejected pass-through adapters as indirection-without-behavior). Local
`delta_structure_dr2` is a verified-bit-identical penm duplicate → deleted. Exported
fns normalized too: `calculate_dr2i_msa`→`calculate_dr2_i_msa`,
`calculate_dr2n_msa`→`calculate_dr2_n_msa` (breaking; GitHub-only/solo/pre-1.0).
`_msa`/`_obs` source labels stay (not index sigs). Profile cols `dr2_i`/`dr2_n`
already correct.

Per project rule (fix-plan-first), updated `dev/plan.md` v0.3b bullet BEFORE coding.
dev/plan.md quantity table (line ~163) already used `dr2_ijm`/`dr2_njm` — consistent.

**Step 0 done (guard locked before any rename):** added
`tests/testthat/test-profile-invariance.R` snapshotting the `dr2_i`/`dr2_n` profiles
+ loglik at (a1=2,a2=5) via `expect_snapshot_value(style="serialize")` — verified
json2 is NOT lossless (~3e-15) so used serialize. Recorded
`_snaps/profile-invariance.md` from CURRENT pre-rename code; stable on re-run (PASS
3, no "Adding"). This is the numeric net: any rename that moves a profile value will
fail. Committing the guard before touching names.

### 2026-06-18 — SESSION END: reorg + roadmap committed & pushed; next = v0.3b

Committed & pushed both pieces of this session to `main`:
- `56dac7d` — the R/ reorg (rule retired, files grouped by `@family`, evaluation
  split into model+objective; see the entry below for detail).
- `992e771` — roadmap reshape in `dev/plan.md`: **v0.3 is now the complete
  site+mode structural-divergence suite** (model AND observed, predict AND fit),
  NOT motion. v0.3a SHIPPED; **v0.3b = consistency tidy** (`dr2→dr2i`; use
  `penm::delta_structure_dr2i`); v0.3c+ = TBD step-by-step (observed dr2i/dr2n from
  homologous structures + alignment; fit dr2n; joint fit). **Motion promoted to its
  own v0.4**; tree route renumbered v0.4→**v0.5**.

Tree clean, in sync with origin/main. PROGRESS.md set dormant. **Next session:
start v0.3b** — enter plan mode, read `R/spm.R` + `data-raw/prepare_znb_data.R` +
the drift test, write the detailed plan (the `dr2→dr2i` rename + `penm::` helper
swap need a deliberate `znb_spm` fixture regen). Memory `project-next-session`
points here.

### 2026-06-18 — File-rename rule RETIRED; R/ reorganized by @family; evaluation split

**Decision (with user): retire the "keep `tmp_src` filenames" migration rule.**
New rule — organize `R/` files by function `@family`/role, not by source filename;
sourceless functions go in the matching family file; `dev/find-source.sh` (content
search) provides traceability instead of filenames. Updated CLAUDE.md (migration
rules), dev/plan.md (resolved decisions), and the `feedback-no-renames` memory
(now states the *opposite* of its old self). Provenance HARD RULES and the
no-function-rename rule are untouched — only *file placement* is freed.

**Naming fix (with user): `@family evaluation` was misleading** ("evaluation" reads
as *assessing*; these are the model's *forward* computation). Split into:
- `@family model` — `calculate_dr2i_msa` / `calculate_dr2n_msa`: the pfix model's
  forward map (predicted divergence profile at given `(a1,a2)`). NOT `prediction`
  — that word is taken by the fit-level `calculate_prediction_*` (posterior/ML
  profile).
- `@family objective` — `calculate_loglik_msa`: one fitting *criterion* (future
  RMSE/robust/other-likelihoods join here). Single-member family for now, by
  design.

**Reorg (pure moves + tag edits, no logic change):**
- `R/spm.R` ← `generate_spm_data` + `delta_structure_*` (`@noRd`) + `preprocess_spm`
  + `preprocess_spm_mode` (consolidates old `generate_spm_data.R` +
  `msa_bayesian_data_preparation.R`, both deleted).
- `R/model.R` ← `calculate_dr2i_msa`, `calculate_dr2n_msa`.
- `R/objective.R` ← `calculate_loglik_msa`. (Old `msa_model_evaluation.R` deleted.)
- `msa_bayesian_workflow.R` / `msa_bayesian_analysis.R` kept (genuinely MCMC,
  `@family fitting`).
- **Boundary recorded so it isn't re-litigated:** `spm` = per-mutant *raw*
  divergence (the data); `model` = pfix-weighted *profile* (the forward map).

**Verification (all green):** baseline captured first (16 exports, 23 man pages,
59 tests pass). After: `document()` → NAMESPACE unchanged, man/ file set unchanged
(only the 6 moved pages' content changed: provenance line + the evaluation trio's
`\concept`/`\seealso` retag `model`/`objective` + first-function title-header
text). `test()` 59 pass / 0 fail (= baseline). `check()` 0 errors / 1 warning /
2 notes (= accepted v0.1 baseline, no new diagnostics). `grep @family evaluation
R/` = 0. `find-source.sh 'preprocess_spm'` and `'calculate_dr2i_msa'` still resolve
to `tmp_src/` despite the new filenames — proves the rule-change premise.
**Visual vignette regression:** installed the refactored tree, re-knit both
`.Rmd.orig`; the only diff was the `sessionInfo()` version banner
(`0.2.0` → `0.2.0.9000`, a precompute artifact) — zero numerical/figure change, so
reverted the version churn. Previews in `dev/preview/` for eyeballing.

### 2026-06-18 — SESSION END / NEXT-SESSION AGENDA

**Pick up here next session.** Open topic the user wants to discuss (NOT yet
decided — do not implement, discuss first):

- **Revisit the "don't rename files when migrating" rule** (CLAUDE.md migration
  rule + [[no file renames]] memory). The rule's original purpose was
  traceability: keep migrated files matching `tmp_src/` filenames 1:1. But this
  session exposed two cracks: (a) the migration ALREADY routinely moves functions
  into files that don't match their `tmp_src` origin, so strict file-to-file
  correspondence is already broken; (b) it produced a genuine mis-housing —
  `preprocess_spm` / `preprocess_spm_mode` (general SPM reshaping, `@family spm`)
  and `calculate_dr2*_msa` (per-point prediction, `@family evaluation`) sit in
  `R/msa_bayesian_data_preparation.R` / live under `msa_bayesian_*` names, even
  though they are NOT Bayesian-specific. The `@family` roxygen tags already
  encode the *real* grouping (spm / evaluation / fitting / decomposition /
  setup / datasets), independent of filenames.
- **The proposal to weigh:** shift the mental model from "migrate whole FILES
  (keep names)" to "migrate FUNCTIONS into well-named, function-grouped files"
  (e.g. a neutral `spm_preprocess.R`, a pipeline-neutral evaluation file, leaving
  `msa_bayesian_*` for genuinely-MCMC code). Trade-off: lose easy
  file-to-tmp_src diffing (mitigated now by `dev/find-source.sh`, which finds a
  source regardless of filename) vs. gain a layout that reflects architecture.
- Decide: keep the rule, retire it, or replace it with "group by @family /
  function role; filenames need not match tmp_src." If retired, a follow-up
  reorg task (pure moves + document() + tests) would realign existing files.

**State at session end:** working tree clean, in sync with origin/main
(`a0dba7e`). Tests 59 pass. v0.3a mode code (`calculate_dr2n_msa`,
`preprocess_spm_mode`) is STAYING — user confirmed the dr2n-analysis vignette
HTML looks correct, so no further verification/revert. The migration guardrails
(unhide archive, find-source.sh, provenance rules) are in place.

### 2026-06-18 — PROCESS FIX: unhide archive + migration-discipline guardrails

After a serious failure (see below), added guardrails so it can't recur.

- **The failure:** Claude "migrated" the v0.3a mode functions poorly and then
  asserted false provenance repeatedly — claimed `calculate_dr2n_msa` was
  "package-native, no tmp_src source" when the source
  (`calculate_diff_mode_msa.R`) was sitting in `tmp_src/.archive/`. Root causes:
  (1) `.archive/` is a HIDDEN dir, so `grep -r tmp_src/` silently skipped its 18
  `.R` files and returned 0 hits → false "no source" conclusion; (2) confused
  penm (a dependency to CALL) with a migration source — even read penm internals
  to reason about provenance; (3) cloned the user's already-refactored `dr2_i`
  function, swapped `i`→`n`, and shipped it WITHOUT verifying against the
  archive's `dr2_n` — skipping the refactor-and-check the user would have done.
- **Fix 1 — unhid the archive:** `mv tmp_src/.archive tmp_src/archive`. Now
  `grep -rln dr2_njm tmp_src/` returns 3 files (was 0). Pure local rename
  (tmp_src is git-ignored + Rbuildignored). Updated all `.archive` refs in
  plan.md/LOG.md.
- **Fix 2 — `dev/find-source.sh`:** mechanical "does this have a migration
  source?" check; greps all of tmp_src/ incl. hidden dirs. Run before claiming
  any code is new.
- **Fix 3 — CLAUDE.md migration rules + memory:** migration sources are ONLY
  tmp_src/; penm is a dependency never a source; migration = restructure + VERIFY
  (reproduce old numbers), not clone-and-rename; never state provenance unchecked.
- **penm access:** kept (rule-only guardrail, no hard block) — user's call.
- **STILL OPEN (not done here):** the v0.3a mode functions on `main` remain
  UNVERIFIED against the archive. Verifying / re-migrating them, and deciding
  whether to revert the v0.3a commit, are separate user-directed follow-ups.

### 2026-06-18 — v0.3a started: mode-form structural divergence (`dr2_njm`)

- Split the roadmap's large v0.3 into slices; **v0.3a** is the first/smallest:
  add the mode-form of `dr2` only (`dr2_njm` / `dr2_n`), introducing the mode
  response axis with no new quantity and no motion. Predict-only.
- Design forks resolved in conversation before planning (the substance of the
  session): (1) keep the **slow per-mutant loop**, not penm's fast `smrs`
  (the `C·F` matrix batch) or `amrs` (analytical) — the loop is the general
  engine that survives into motion (v0.3b) and trees (v0.4) and is the
  benchmark; `dr2n` is just one more cheap reduction of the same `dr`.
  (2) **Raw SPM stays one "all-effects" object** — `dr2n` is another SPM
  column next to energies + `dr2`; the user's framing: energy changes are
  mutation effects too, and all effects of one mutation share the `(j,m)` key.
  (3) **Site/mode separation lives at the preprocess/evaluate layer**, which is
  a reshape-for-purpose step, not in the raw object: new `preprocess_spm_mode()`
  + `calculate_dr2n_msa()`; a mode analysis never touches `dr2mat`/`site_map`.
  (4) **Predict-only** — no fit (no observed mode data), v0.1 site fit untouched.
- Vignette: new **`dr2n-analysis`** topic vignette (best practice = vignettes by
  task, not feature accumulation; intro stays the narrow on-ramp). v0.3a content
  = predict-only `a1`/`a2` sweep of the `dr2_n` profile; a fitting section is
  added there later when alignment-derived `dr2_n` + a mode fit exist.
- Verified pre-planning (penm local source at `../penm`):
  `delta_structure_dr2n(wt,mut)` returns 678 finite per-mode values on `znb_wt`
  (678 = 3·228−6), same call shape as `delta_structure_dr2`; `get_mode` is NOT
  exported (use `penm:::get_mode`, as the archive did). Memory written:
  `project-penm-mutscan-tiers`, `project-msamodel-v03-mode-axis`.

### 2026-06-18 — vignette preview tooling + .Rbuildignore fix (post-v0.3a)

- **Problem:** rendering a committed `vignettes/<name>.Rmd` to HTML in place
  DESTROYS its figures — `html_vignette` uses pandoc `--self-contained`, which
  base64-embeds the figures then deletes the `<name>_files/` dir the shipped
  `.Rmd` references. Hit during the v0.3a commit (figures vanished, nearly
  committed missing).
- **Fix:** `dev/preview-vignette.R` — renders a COPY of the `.Rmd` + `_files/`
  in a temp dir and copies only the self-contained `.html` back to
  `vignettes/<name>.html` (git-ignored). The repo's `vignettes/` is never
  written to. Verified: figures present & `git status` clean before and after.
  Usage `Rscript dev/preview-vignette.R <name>`. Documented in CLAUDE.md.
  (Output goes to `dev/preview/<name>.html`, a git-ignored scratch dir — NOT
  `vignettes/`, which should hold only files that ship.)
- **Also fixed a latent v0.3a bug:** `.Rbuildignore` had a filename-specific
  `^vignettes/msamodel-intro\.Rmd\.orig$`, so the new `dr2n-analysis.Rmd.orig`
  would have SHIPPED in the build tarball. Generalized to
  `^vignettes/.*\.Rmd\.orig$` (ignores all `.orig` sources; rendered `.Rmd`
  still ships). Verified both `.orig` ignored, both `.Rmd` not.

### 2026-06-18 — v0.3a implemented (steps 1-8 + verification)

- `generate_spm_data()` now emits `mode` + `dr2n` SPM list-columns
  (`dr2n <- penm::delta_structure_dr2n(wt, mut)` per mutant; `mode <-
  seq_along(dr2n)`). Added `preprocess_spm_mode()` (builds `dr2nmat`, no
  `pdb_site` map) and `calculate_dr2n_msa()` (predict-only; identical
  mutant-axis reweighting). `get_mode` stays `penm:::` (not exported). Exports
  14 → 16.
- Regenerated `znb_spm`. The one-time `tmp_src` validation now compares only
  shared columns (`znb_spm[names(orig)]`) — the migrated SPM legitimately has
  extra columns the 2021 original lacks; shared columns still match at 1e-8.
- Tests: added `mode == penm:::get_mode(znb_wt)` guard to test-spm-generate;
  new `test-msa-mode.R` (shape; a **Parseval** site/mode total-equality check —
  basis invariance, a real invariant not a recomputation; and an
  independent-route reweighting check). Suite 45 → 59 pass.
- New `dr2n-analysis` vignette (predict-only `a1`/`a2` sweep of the per-mode
  profile; consistency check shows site_total == mode_total). NOTE: knitting
  required `devtools::install()` first — the vignette's `library(msamodel)`
  loads the INSTALLED package, which lacked the new functions until installed.
  (Same will apply to any future vignette using new code: install, then knit.)
- Generalized `.githooks/pre-commit` to guard every `vignettes/*.Rmd.orig`
  pair, not just the intro.
- DESCRIPTION → `0.2.0.9000`; NEWS development section added.
- `check()` = 0E/1W/2N (baseline). **Finding (not acted on):** the data-size
  NOTE grew — `znb_spm` 12.5 → 23.6 MB, because `dr2n` (678 modes/row) is ~3x
  wider than `dr2` (228 sites/row). Each future v0.3 quantity adds another wide
  column to the same fixture, so fixture-size management (e.g. dropping the
  redundant `dr` raw-vector column, or shipping a smaller example) is worth a
  deliberate decision before v0.3b. Not fixed here — out of v0.3a scope.
- v0.3a COMPLETE; NOT committed (awaiting user review).

### 2026-06-11 — vignette: switch to the `.Rmd.orig` precompute pattern

- **Why:** the intro vignette had drifted and violated "fail loud". Heavy compute
  lived in a *separate* `vignettes/precompute.R` → `vignette_cache.rds`; the `.Rmd`
  showed hand-typed `eval=FALSE` snippets *resembling* that script and read the
  cache. Two hand-maintained copies of the pipeline code with nothing comparing
  them → they diverged: the `sweep-code` chunk showed only `sweep_a1` while the
  cache (and a rendered figure) also had `sweep_a2`, so a figure was produced by
  code the reader never saw. Worse, every figure chunk was gated
  `eval = have_cache && has_ggplot`, so a missing/stale cache *silently skipped*
  figures instead of erroring.
- **Fix (rOpenSci `.Rmd.orig` pattern):** `vignettes/msamodel-intro.Rmd.orig` is
  now the single source of truth — every pipeline chunk is real `eval=TRUE` and
  computes inline; the shipped `msamodel-intro.Rmd` is *pre-rendered* from it via
  `knitr::knit("vignettes/msamodel-intro.Rmd.orig", output="vignettes/msamodel-intro.Rmd")`.
  Shown code == executed code (one copy), so the drift class is gone structurally.
  Deleted `precompute.R` and `vignette_cache.rds`. The MCMC fit runs at 500 iter /
  100 burn-in (fast render). The only `eval=FALSE` chunk is the SPM-scan recipe in
  §1, explicitly labelled as the code that produced the shipped `znb_spm` (already
  drift-guarded by `test-spm-generate.R`).
- **No drift test (deliberate):** a "recompute-and-assert-equal" test (the guard
  `znb_spm` uses) is wrong for a vignette — a vignette is *meant* to change across
  versions, so such a test fails on every intentional edit. The `.Rmd.orig` pattern
  already removes the root cause (two copies). `set.seed(1024)` stays at the top of
  `.Rmd.orig` only to keep the committed render stable, not as an asserted invariant.
- **Stale-render guard = a pre-commit hook** (`.githooks/pre-commit`, enabled via
  `git config core.hooksPath .githooks`). It blocks a commit that stages
  `msamodel-intro.Rmd.orig` *without* also staging the re-rendered
  `msamodel-intro.Rmd`. It does NOT recompute/diff content (no brittleness, doesn't
  fight intentional edits) — it only enforces "re-knit was staged alongside the
  source edit." `--no-verify` bypasses for genuine no-render edits. Tested: blocks
  when only `.orig` staged; passes when both staged or `.orig` untouched. (A buried
  comment alone was no safeguard — it can't tell you that you forgot.) `.githooks`
  is `.Rbuildignore`d. The re-knit command is also in a comment atop `.Rmd.orig` for
  reference.
- **Wiring:** `.Rbuildignore` now ignores `msamodel-intro.Rmd.orig` (was
  `precompute.R`). `fig.path = "msamodel-intro_files/figure-html/"` so figures don't
  land in a bare `figure/` dir (which `R CMD check` flags as a knitr leftover —
  that was the one new NOTE, now cleared). Figures are checked in (un-ignored in
  `vignettes/.gitignore`).
- **Verify:** knit ran all chunks incl. the live MCMC (no skips); rendered `.Rmd`
  carries embedded results (`a1≈0.455`, `a2≈44.4`) and references both sweeps;
  `tools::pkgVignettes()` sees exactly one vignette; `check()` back at baseline
  (0 errors, 1 warning, 2 notes); the built `doc/msamodel-intro.html` embeds all
  10 figures as base64 (0 broken path refs). NOTE: the 2026-06-11 v0.2 entry below
  describes the now-removed `precompute.R`/cache mechanism — superseded by this
  entry, not edited (append-only).

### 2026-06-11 — v0.2 COMPLETE: all four API changes executed

- **All four v0.2 changes done; `check()` = 0 errors, 1 warning, 2 notes — same
  as the v0.1 baseline (no new issues).** Tests 55 → 45 pass (the drop removed
  the grid test file + the `get_active_site` parse test + 3 contract grid lines;
  no coverage lost — the dropped assertions were redundant with the manual i-join
  test and the prediction/decomposition pdb_site checks). Exports 18 → 14.
  DESCRIPTION bumped to 0.2.0; NEWS.md written (changes flagged breaking).
- **Step 1 (rename):** `shap_* → phi_*` in `R/msa_decomposition.R` (3 columns +
  both validation `required_cols` + `starts_with("phi_")`), `test-decomposition.R`,
  vignette prose + `labs_comp`, roxygen prose in `msa_bayesian_workflow.R`,
  `DESCRIPTION` Description field (also dropped its "grid-based exploration"
  clause), `CLAUDE.md`. NO `data/*.rda` regen. `man/*.Rd` regenerated. `tmp_src/`
  untouched.
- **Step 2 (bio3d object):** deleted `load_protein` (one-line `read.pdb`
  wrapper). Added `inherits(pdb,"pdb")` guard to `setup_enm` with a clear error
  pointing at `read.pdb`/`read.cif`. Updated callers (test-spm-generate, vignette
  precompute + Rmd) to `bio3d::read.pdb`. **mmCIF verified empirically:**
  `bio3d::get.pdb("1znb", format="cif")` → `read.cif` (class "pdb", 458 sites for
  the full multi-chain CIF) → `setup_enm` OK, `get_msf_site` works. Caveats
  recorded in NEWS: bio3d flags `read.cif` as beta; it warned helix/sheet records
  unparsed (`set_enm` doesn't need them).
- **Step 3 (active-site vector):** deleted `get_active_site` → `R/pdb_utils.R`
  emptied → `git rm`'d the file. Inlined `c(99,101,103,162,181,184,193,223)`
  (from `data-raw/raw/dataset_1znb_A.csv`) in `data-raw/prepare_znb_data.R`,
  `test-spm-generate.R` (`PDB_SITE_ACTIVE` const), `vignettes/precompute.R` +
  `.Rmd`. Kept `znb_dataset` shipped; rewrote its data-doc as *illustrative*
  (its `@seealso get_active_site` would've broken check). Moved `bio3d`
  Imports→Suggests and dropped `stringr` from Imports + the two `@importFrom`
  lines (both unused in `R/` after the deletions).
- **Step 4 (drop grid):** `git rm` `R/msa_a1a2grid_workflow.R` +
  `tests/testthat/test-a1a2grid.R`. Removed the grid lines from
  `test-contract.R` (pdb_site-on-output already covered by the prediction/
  decomposition checks; pdb_site alignment by the independent manual-join test).
  Replaced vignette §6 with a `preprocess_spm` once → `map_dfr` over
  `calculate_dr2i_msa` example (makes reshape-once explicit). Reworded
  `site_properties.R` `@param site_data`. globalVariables unchanged (`dr2_msa`
  still used in `calculate_loglik_msa`; `dr2_mm/ms/ma` were never listed).
- **Verify:** `document()` clean (14 exports, 4 dropped `.Rd` deleted, no stale
  bio3d/stringr in NAMESPACE); `test()` 45 pass / 0 fail / 0 warn; reran
  `precompute.R` (~100s MCMC) — cache now holds `phi_*` (confirmed no `shap` left
  anywhere in the `.rds`), so the built vignette matches the renamed columns;
  `check()` clean vs baseline. The MCMC is seeded (`precompute.R:20`
  `set.seed(1024)`, base-R RNG), so the rerun is deterministic: every numeric
  value in the cache is bit-for-bit identical to the committed one (verified
  parameter_summary + all.equal(prediction_summary) TRUE); only the `component`
  factor labels changed `shap_*`→`phi_*`. (No drift — an earlier claim of
  "slight drift" was wrong.)

### 2026-06-11 — v0.2 started: API decisions (plan-docs-first)

- **v0.2 entered plan mode; four API decisions made with the user, recorded in
  `dev/plan.md` BEFORE any code** (fix-plan-first). No new model capability;
  exports 18 → 14; all four are breaking; `tmp_src/` not touched.
  1. **`shap_*` → `phi_*`** rename (misnomer; φ in the paper). Load-bearing edits
     confined to `R/msa_decomposition.R` (cols + validation + `starts_with`),
     `test-decomposition.R`, vignette labels/prose, roxygen prose. No `data/*.rda`
     regen (no `shap_` columns there) — BUT verified `vignette_cache.rds` carries
     `shap_act/mut/stab` as *values* in its `component` column, so
     `precompute.R` must be rerun (~100s MCMC; demo R²/summaries will drift
     slightly — cosmetic, expected).
  2. **Structure input = bio3d pdb object, NOT a path** (changed from the
     roadmap's earlier "accept a `.pdb` path"). Rationale (user-driven): the real
     boundary is already `setup_enm(pdb, ...)`; `load_protein` is a one-line
     `read.pdb` wrapper forcing an ID→filename layout. Drop it; user calls
     `bio3d::read.pdb()`/`read.cif()`. Decouples from file I/O and gets **mmCIF
     for free** (bio3d reads mmCIF via `read.cif`; the package no longer sniffs
     extensions/tracks formats). Add `inherits(pdb,"pdb")` guard to `setup_enm`.
     CAVEAT to verify empirically (penm not in repo): does `set_enm` work on a
     `read.cif` object, or need a `.pdb`-only field? If the latter, dial back the
     mmCIF docs claim.
  3. **Active-site input = plain `pdb_site_active` integer vector.** Drop
     `get_active_site` (only thing forcing the `dataset_ec2024.csv` lookup);
     downstream (`generate_spm_data`, `add_site_properties`) already take the
     vector. Inline 1znb_A's vector `c(99,101,103,162,181,184,193,223)` in
     data-raw/tests/vignette. **Keep `znb_dataset`** shipped as *illustrative* of
     the source CSV format (no longer consumed by package code).
  4. **Drop the grid API** (`define_selection_grid` +
     `calculate_dr2i_msa_a1a2grid`). No real consumer — fit + vignette sweeps use
     `calculate_dr2i_msa(pp, ...)` directly; the grid fn also re-ran
     `preprocess_spm` internally (violates reshape-once). A sweep is a 3-line
     `map_dfr` in user code. Delete `R/msa_a1a2grid_workflow.R` +
     `test-a1a2grid.R`; re-prove the `pdb_site`-alignment property in
     `test-contract.R` via `calculate_dr2i_msa` + `site_map` (not the grid route).
- **Also tightened the durable "precomputation" finding in `dev/plan.md`** — it
  conflated three stages the user wants kept distinct: SPM *physics*
  (`generate_spm_data`, the actual precomputation, once) vs *reshaping*
  (`preprocess_spm`, cheap, not a computation) vs *reweighting*
  (`calculate_dr2i_msa`, the only a1/a2-dependent step).
- `dev/PROGRESS.md` rewritten as the v0.2 checklist (Step 0 done). No v0.2 code
  written yet.

### 2026-06-10 — planning restructured + full tmp_src survey

- **Restructured planning into two tiers (user-driven).** Lesson from v0.1: the
  575-line `dev/plan.md` was both blueprint and execution script, so every
  execution surprise forced a plan edit. Now: `dev/plan.md` = COARSE durable
  roadmap (what each version is) + durable findings; detailed planning is done
  per-version at execution time (in plan mode, against the actual code) and tracked
  in `dev/PROGRESS.md`. Rewrote `dev/plan.md` (dropped the 13 v0.1 sections — that
  job is done; history is here), repurposed `dev/PROGRESS.md` as the in-flight-
  version checklist (currently dormant), and rewrote `CLAUDE.md` for the post-v0.1
  reality + the two-tier model.
- **Full read of `tmp_src/` (.R + .Rmd) — much more unmigrated than §12 tracked.**
  The old §12 only listed unmigrated `R/*.R` files. Surveyed all four regions:
  (A) deferred R/ model code [assessment 3-file bundle, allotment, protein decomp,
  viz]; (B) `Rmd/` paper analyses [figure library, AlphaFold2 validation,
  sequence-vs-structure divergence — new feature signals]; (C) `someday_maybe/tree/`
  trajectory route; (D) **`archive/` — the model's MOTION/MODE arm**, which the
  user wants for future projects and which was buried/never analyzed.
  *(2026-06-18: this dir was renamed `.archive/`→`archive/` — unhidden so
  `grep -r` finds it; path updated here to stay valid.)*
- **KEY FINDING (durable, now in plan.md): the precomputation property generalizes
  for the SPM-mean route but NOT the tree.** The model is a 2×2 grid {structure
  dr2, motion dh} × {site i, mode n} (+ nh_n); v0.1 implements only structure×site.
  All archive quantities are `X = sum(pfix_jm * X_jm)/sum(pfix_jm)` with a1/a2 only
  in per-mutant pfix → the v0.1 `[mutant×site]` reweighting matrix generalizes to
  more value-matrices reweighted by the same weights. So motion/mode (D) ≈ "more
  SPM columns + reweighting", mechanically close to v0.1. The TREE (C) uses
  sequential substitutions → weights don't factorize → it interpolates over an
  a1/a2 grid (akima) instead; no clean precomputation. Hence roadmap order
  v0.3=motion/mode (D), v0.4=tree (C).
- **Roadmap set (user decisions):** v0.2 = API-only (phi rename, PDB path input,
  pdb_site-vector active sites, grid keep/demote/drop); v0.3 = motion/mode (D),
  read rest of `archive/` incl. model_rates.R at its start; v0.4 = tree (C);
  later/unordered = assessment+allotment+protein-decomp, viz/figure-lib, AF2,
  sequence divergence. No v0.2 code written yet.

### 2026-06-10 (earlier)
- **§9 step 10 done — `devtools::check()` as clean as v0.1 gets (GitHub-only, not
  CRAN).** Final: **0 errors, 1 warning, 2 notes**, tests pass. Fixed the two
  genuinely-fixable items: (N2) `.Rbuildignore`d top-level `CLAUDE.md`; (N4)
  silenced the package-wide tidyverse-NSE "no visible binding" NOTE via
  `utils::globalVariables(c(...))` in `R/msamodel-package.R` (28 data-column /
  grouping names — not real globals). **Accepted (deliberate, user is GitHub-only):**
  LazyData WARNING + installed-size NOTE (both the known `znb_spm`~13 MB/`znb_wt`~6 MB
  data-size story) and the spurious future-timestamp NOTE (clock/network skew).
  **Dead end:** tried removing `LazyData: true` to clear the WARNING — broke 19
  tests (`object 'znb_spm' not found`: tests use the `znb_*` datasets bare, relying
  on lazy auto-load). Reverted; not worth editing the whole test suite to chase a
  cosmetic warning. **v0.1 migration COMPLETE** — all 10 plan §9 steps done. Next:
  the v0.2 backlog in plan §12 (rename "Shapley"→phi_*, reconsider the grid API,
  the deferred assessment/allotment/protein-decomp/trajectory/viz layers).
- **Vignette revisions (user review of step 9).** Per user feedback on the rendered
  HTML: (3) added active-site vlines to the a1/a2 sweep plots; (4) moved the R²
  annotation to top-left on the obs-vs-fitted profile (was overlapping); added a
  "four model variants" table (MM/MS/MA/MSA) — they were referenced but never
  defined; (5) rewrote the decomposition section — **dropped "Shapley" and the
  wrong "baseline" gloss** (the latter was my invention, grounded in nothing — the
  `shap_mut` term is the mutation/no-selection component, which still varies among
  sites; a baseline would be residue-independent), added the explicit φ formulas
  (φ_mut, φ_stab, φ_act with the ½-averaging) + the additivity identity; (6)
  demoted the grid section to a brief mention. Flagged both v0.2 items in plan §12.
  Re-rendered + re-inspected the changed figures; cache unchanged (only plot/prose
  edited).

### 2026-06-10
- **§9 step 9 done — vignette `msamodel-intro.Rmd` written, builds clean.** Full
  plotted walkthrough (supersedes the old §8 "no plots, print() only" spec —
  rewrote §8 first, fix-plan-first). **Build mode: precompute + cache** —
  `vignettes/precompute.R` runs the heavy pipeline once (SPM preprocess, profiles,
  a1/a2 sweeps, a 4000-iter/1000-burn MCMC fit, decomposition) and saves
  `vignettes/vignette_cache.rds` (96 KB); the `.Rmd` shows the pipeline code as
  `eval=FALSE` and draws every plot off the cache. Keeps check fast, no `penm`/MCMC
  at build time. `precompute.R` is `.Rbuildignore`d; the cache ships. ggplot2 is
  Suggested, so every plot chunk is guarded by `requireNamespace("ggplot2")` +
  `eval = have_cache && has_ggplot`.
  **Sections:** (1) setup from PDB (load_protein→get_active_site→setup_enm→
  generate_spm_data, eval=FALSE, noted as reproducing the shipped znb_wt/znb_spm);
  (2) profile at a demo (a1=2,a2=500) — raw lrmsd vs residue (active sites marked),
  vs lrmsf, vs dactive (+geom_smooth); (3) a1-sweep (a2=0) and a2-sweep (a1=0),
  whole raw lrmsd; (4) Bayesian fit — parameter summary + observed-vs-fitted on
  **mean-centered nlrmsd** (what the likelihood compares) as profile/scatter/vs-
  lrmsf/vs-dactive, each annotated with **R² (≈0.605)** computed inline (v0.1 has
  no gof fns — §1); (5) Shapley decomposition (shap_mut/stab/act, active sites
  marked — activity-selection troughs visibly coincide); (6) grid API note.
  **Decisions (user):** gof = R²; centering differs by section (raw for profile/
  sweep, centered for fit-vs-obs); full PDB setup path shown. **Verify:** ran
  precompute (~99 s, MCMC dominates) → cache written; visually inspected all 10
  figures as standalone PNGs (Read tool) and iterated (added active-site vlines to
  the decomposition); `rmarkdown::render` → 1 MB self-contained HTML, all 46 chunks
  ran; **`devtools::check()` = 0 errors**, 1 warning + 4 notes, **all preexisting /
  step-10** (LazyData + installed-size = the known znb_spm/znb_wt data-size issue;
  top-level CLAUDE.md; future-timestamp; the package-wide tidyverse NSE
  "no visible binding"). Confirmed via git that this work touched only
  `.Rbuildignore`, `dev/plan.md`, `vignettes/` — DESCRIPTION/data untouched, so the
  warning/notes are not vignette-introduced. Only step 10 (`check()` clean) remains;
  the data-size warning/note is the main thing to resolve there.

### 2026-06-09 (later)
- **§9 step 8 done — test suite written (8 files, 55 tests, all pass).** Scaffolded
  `tests/testthat.R` + `tests/testthat/`. Files: spm-generate (the fixture-drift
  guard — reproduces znb_wt/znb_spm from the §5 params, ~20s, dominates runtime),
  spm-preprocess (dr2mat shape + the site_map mapping the pdb_site contract relies
  on), msa-evaluate (dr2i finiteness + a FROZEN loglik literal -184.3241923285),
  msa-mcmc (prior-range validation errors + deterministic fix_a1/fix_a2 constancy
  + shape), decomposition (validation stop()s + the pdb_site optional branch both
  ways + summary shape — additivity NOT asserted, it's tautological), a1a2grid
  (54-row grid + output cols + pdb_site alignment), workflow-endtoend (7-element
  return, real pipeline), contract (unknown-pdb_site error, subset-coverage finite
  loglik, independent-route mapping equality, pdb_site on outputs). Suite: 55 pass
  / 0 fail / 0 warn, ~29s. **Failure-sanity check:** perturbed the frozen loglik
  and the dr2mat shape → the relevant tests FAILED, then reverted → clean. Confirms
  assertions aren't vacuous.
- **Plan §7 reconciled for step 8 (tests).** Updated the test matrix to the
  pdb_site contract and to vetted assertions (the prior matrix had tautological /
  untestable rows). Key changes: dropped the Shapley-additivity assertion (it is
  algebraically the implementation — can't fail except by editing the formula);
  replaced the loglik "regression" with a FROZEN literal (captured once, not
  recomputed in-test); dropped the "acceptance rate > 0" claim (`accepts` isn't
  returned); added deterministic `fix_a1/fix_a2` and validation-error tests for
  MCMC; made decomposition test validation + the pdb_site optional-branch + shape.
  Added **`test-contract.R`** (8th file) for the pdb_site contract: unknown-pdb_site
  error, subset-coverage finite loglik, independent-route mapping equality (non-
  circular), pdb_site on outputs. v0.1 test count 7 → 8. Done before writing tests
  (fix-plan-first). See memory `feedback-test-quality`.

### 2026-06-09
- **Fixed the `observed_data` site-key contract (pdb_site, not i) + documented
  datasets (step 7).** Surfaced while prepping step 7: `observed_data` is
  *user-provided* input (public arg of `run_msa_bayesian_analysis`), but the
  contract required the user to supply `i` — the model's internal 1..N response-
  site renumbering — which a user can't reliably produce, and the fit joined on it
  with no validation (silent-misalignment risk). For 1znb_A it worked by luck
  (verified 225/225 pdb_site agreement), not by construction.
  **Decision (user): full fix overriding CLAUDE.md §4 "rename nothing"** — justified
  as a public-API correctness fix. New contract: `observed_data = {pdb_site,
  lrmsd_obs}`; package maps pdb_site→i internally via a new `site_map` returned by
  `preprocess_spm`, erroring on any unknown pdb_site. All six site-keyed objects
  now also carry `pdb_site` alongside `i` (prediction_samples/summary,
  decomposition_samples/summary, observed_data passthrough, grid output). Internal
  i/j convention (i=response site, j=perturbed site) unchanged. Files:
  `R/msa_bayesian_data_preparation.R` (site_map), `R/msa_model_evaluation.R`
  (map+validate+join), `R/msa_bayesian_analysis.R`, `R/msa_decomposition.R`
  (pdb_site optional — synthetic test input without it still works),
  `R/msa_a1a2grid_workflow.R`, `R/msa_bayesian_workflow.R` (docs).
  **`znb_profile` regenerated** as `{pdb_site, lrmsd_obs}` (dropped i; dropped the
  derivable dactive/lrmsf — verified reproducible from znb_wt to machine
  precision: dactive=get_dactive, lrmsf=log(sqrt(get_msf_site))). SPM one-time
  validation still passed on re-run.
  **Step 7:** wrote `R/data-doc.R` documenting all 5 datasets (`@family datasets`,
  accurate `@format`; znb_profile framed as the {pdb_site, lrmsd_obs} fit target).
  **Verify:** correctness gate — new contract gives loglik identical to the old
  i-join reference for 1znb_A (-184.3242, all.equal TRUE); negative test — unknown
  pdb_site errors clearly; `document()` clean, NAMESPACE still 18 exports;
  `load_all()` clean; full workflow returns 7 elems with pdb_site on every
  site-keyed output. Updated `dev/plan.md` §2/§5 (contract note + deliberate §4
  exception) BEFORE coding. v0.2 maybe: make pdb_site the canonical key throughout.

### 2026-06-06
- **§9 step 6 done — embedded `znb_*` data, SPM generated + validated.** Vendored
  the small inputs (`inst/extdata/1znb_A.pdb`; `data-raw/raw/dataset_1znb_A.csv`
  1 row; `data-raw/raw/profiles_1znb_A.csv` 225 rows), wrote
  `data-raw/prepare_znb_data.R` (generates `znb_spm` via the package's own
  `generate_spm_data`, params as named constants flagged to match the step-8
  test), ran it. **One-time validation PASSED:** `all.equal(znb_spm,
  readRDS(tmp_src/.../1znb_A_spm.rds), tol=1e-8)` is TRUE — the migrated
  `generate_spm_data` reproduces the source SPM (migration-correctness proof;
  block fenced as removable). Saved 5 `.rda` via `use_data(compress="xz")`.
  Sanity: znb_spm 2508 rows, znb_profile 225, 8 active sites, dr2mat 2280×228.
  Added `^data-raw$` to `.Rbuildignore`; wrote `data-raw/source_files.txt`. Scan
  runtime ~20 s. **Flag for step 10:** `znb_spm.rda` ≈ 13 MB and `znb_wt.rda`
  ≈ 6 MB even xz-compressed — will trigger an `R CMD check` data-size NOTE; needs
  a decision later (subsample? fewer mutations? accept the NOTE?). Not a step-6
  blocker.
- **Plan change — §5 embedded-data approach: generate the SPM, don't copy it.**
  Rejected the original §5 (copy `tmp_src/.../1znb_A_spm.rds`, point `SRC` at
  `tmp_src/`/the live project): not self-consistent — once `tmp_src/` is deleted,
  `data-raw/` would have no in-repo recipe to rebuild `data/`, and pointing at the
  live project isn't reproducibility. **New approach:** vendor only the *small*
  inputs (PDB → `inst/extdata/`, the 1znb_A rows of the two CSVs →
  `data-raw/raw/`), and **generate** `znb_spm` with the package's own
  `generate_spm_data()`. Validate the generated SPM **once** against the
  `tmp_src/` original via `stopifnot(all.equal(...))` (migration-correctness
  proof, fenced as removable) — this is NOT the permanent test; the
  regenerate-and-compare drift guard remains step 8's `test-spm-generate.R`.
  **Feasibility verified before deciding:** recovered the generation params from
  `tmp_src/scripts/03_scan_mutants_all_cases.R` (n_mutations=10, model="lfenm",
  sigma=0.3, min_sd=2, seed=1024; ENM node="ca", ming_wall, d_max=10.5,
  frustrated=FALSE); confirmed the scan is seed-deterministic — penm's
  `get_mutant_site_lfenm` does `set.seed(seed + site_mut*mutation)` per mutant
  (tested: same seed → identical mutant xyz, different seed → different;
  penm 0.2.0.9000, R 4.2). The .rds carries no provenance, so the params are
  hard-coded in the script and must match the step-8 test (option a). **Edited
  `dev/plan.md` §5 (rewrite), §2 (layout: znb_spm "generated"; data-raw adds
  `raw/`), §9 step 6, §13 decision 5.** Done as a standalone turn *before* writing
  any step-6 code, so a fresh context reads the corrected plan (per
  `CLAUDE.md` "read plan.md before doing anything"). Step-6 code execution is the
  next turn.

### 2026-06-05
- **§9 step 4 done (all sub-steps 4.1–4.7), in one batch + step 5.** Copied the 7
  v0.1 source files from `tmp_src/R/` into `R/` keeping source filenames (per the
  no-rename change below), in dependency order. Applied §4 rules: stripped all
  `@requires` tags; added `@export` to the 18 public fns + `@family` groups; the
  2 helpers (`delta_structure_dr`, `delta_structure_dr2`) marked `@noRd`. No
  inline `library()` existed to strip. **Bug fixes:** §11.1 — deleted the
  duplicate `define_selection_grid` (the `a1_range, a2_range` override at source
  lines 106–119), keeping the real defaulted one; §11.3 — fixed the
  `"lower,Quite upper"` roxygen typo in `calculate_prediction_summary`'s
  `@return`. **4.7 assessment removal:** dropped the two
  `calculate_model_comparison_samples/summary` blocks + their return-list entries
  from `run_msa_bayesian_analysis` → 7-element return; updated its `@return`.
  **Imports:** added `transmute` + `group_split` to the dplyr `@importFrom` in
  `msamodel-package.R` (used by `msa_decomposition.R`). **Verify:**
  `devtools::document()` clean (no roxygen warnings; 18 `.Rd` written, none for
  the `@noRd` helpers); NAMESPACE has exactly 18 exports, `delta_structure_*`
  absent; `devtools::load_all()` clean. Behavioral/numeric checks deferred to
  steps 6–8 (no data/tests yet).
- **Plan change — no file renames in v0.1.** Dropped the §2 rename list (10
  renames) before starting §9 step 4. Files are now copied from `tmp_src/R/`
  keeping their **source filenames**. Rationale: source names keep provenance
  obvious during the migration, and some rename rationales were v0.1-scoped (e.g.
  `msa_bayesian_workflow.R` → `msa_workflow.R` "because only one workflow in
  v0.1") — an assumption v0.2 breaks, which would force a *second* rename (churn
  the migration explicitly avoids). File naming/layout is now an explicit later
  refactor, not folded into the copy. Edited `dev/plan.md` §2 (rename block + `R/`
  layout tree) and §9 step-4 sub-steps; regenerated the §9-step-4 labels in
  `dev/PROGRESS.md` to match. Deferred-file renames mentioned elsewhere in the
  plan are left as advisory (decided when those files are actually migrated).

### 2026-06-04
- **Docs split:** `dev/PROGRESS.md` was doing two jobs (append-only log + a
  checklist meant to mirror plan.md §9); the re-scope made the checklist
  numbering diverge from the plan and lose meaning. Split into three files:
  `dev/plan.md` (normative, unchanged role; §9 step-4 re-labelled `4.1`–`4.7`),
  `dev/PROGRESS.md` (pure checklist, bare 1:1 mirror of §9), and this
  `dev/LOG.md` (history). Updated the `## Tracking progress` convention in
  `CLAUDE.md` and the cross-references.
- **v0.1 RE-SCOPED (major):** narrowed v0.1 to "compute divergence profiles + fit
  to data". Traced the real call graph from `run_msa_bayesian_analysis` + the grid
  entry points to decide scope. **Kept:** setup, SPM, profiles/eval, grid, MCMC
  fitting + the four model profiles, site-level Shapley decomposition (wired into
  the workflow). **Deferred to vNext:** the model-fit *assessment* layer
  (`model_comparison_functions.R` + its sole-consumer `loess_compare.R` + the
  `utils.R` error metrics), *allotment* (unsettled), protein-level decomposition,
  trajectory route, visualization. Decisions: keep site decomposition; defer
  assessment (NOT the profiles — those are kept); pull `utils.R` + `loess_compare.R`
  back out. **Reverted** R/utils.R + R/loess_compare.R + their man pages + 7
  exports (forward commit, no history rewrite); `document()`/`load_all()` clean,
  NAMESPACE now 0 exports, R/ = just `msamodel-package.R`. **Rewrote `dev/plan.md`
  §1/§2/§6/§7/§8/§9/§11/§12** to match (exports 39 → 18; tests 9 → 7; workflow
  edited to drop assessment calls → 7-elem return). Reconciled `CLAUDE.md` and
  the progress docs. The §9 step-4 sub-step list is the new dependency order.
  Next: resume migration at §9 step 4.1 (`pdb_utils.R`, `enm_setup.R`,
  `site_properties.R`) against the re-scoped plan.
- **§9 step 4.2 done** (later REVERTED in the re-scope above): copied
  `compare_loess_fits.R` → `R/loess_compare.R` (renamed per §2).
  `compare_loess_fits` + `compare_loess_rmse` exported with
  `@family loess comparison`; `plot_loess_comparison` kept `@noRd` (deferred to
  v0.3), inline `library(ggplot2)` stripped, every ggplot2 call qualified
  `ggplot2::` behind the existing `requireNamespace` guard (§4 rule 2/8).
  Dropped the orphaned `@param evaluation_points` (not a real arg — doc-only).
  `gridExtra` left out of Suggests and no `globalVariables()` shim added (both
  decided: minimal change, the only consumer is the deferred internal plot fn).
  `document()` wrote 2 `.Rd` (none for the `@noRd` fn) + NAMESPACE (then 7
  exports); `load_all()` clean; smoke test of both exported fns returned finite
  values with the expected result names.
- **Toolchain decision (not a §9 step):** investigated whether to update the R
  stack first. Finding: R itself is old (4.2.1; latest 4.6.0) but the packages
  are current (dplyr 1.1.4, ggplot2 4.0.1, etc.). **Decision: do NOT upgrade R
  now** — 4.2.1 satisfies the plan (`R >= 3.5`); revisit after v0.1. When done,
  upgrade R **side-by-side via `rig`** (keep 4.2 for the MSA paper project),
  reinstall packages into the new library, rebuild `penm` from source — that's
  the safe path (single shared library; all 436 pkgs built under 4.2.x). The
  only risk to penm/paper is upgrading *R*, not packages.
  **Updated roxygen2 7.2.3 → 8.0.0** (compiles under R 4.2; `Depends: R >= 4.1`):
  resolves the version-mismatch warning. 8.0.0 records its version as the new
  `Config/roxygen2/version: 8.0.0` field and removed the old `RoxygenNote` line.
  8.0.0 also enforces **one-line `@importFrom`** — flattened the wrapped dplyr/
  stats/penm directives in `msamodel-package.R`. `document()` + `load_all()`
  clean; 5 utils exports intact (those exports later reverted in the re-scope).
- **§9 step 4.1 done** (later REVERTED in the re-scope above): copied `utils.R`
  (`rmse_trend`, `rmse`, `r2`, `mrr`, `mrr_trend`). Stripped the `@requires`
  tags, added `@export` + `@family error metrics` to all 5, kept
  signatures/logic identical (no `library()` calls in source). `document()`
  wrote 5 `.Rd` + NAMESPACE (5 exports); `load_all()` clean. Bare
  `loess`/`predict`/`cor` resolve via the existing `@importFrom stats`.
- **§9 step 3 done:** added `R/msamodel-package.R` with the package-wide
  `@importFrom` directives (Option A namespacing — imports declared once;
  copied function files use bare calls). Inventory verified against the 16
  migrated source files; ggplot2 excluded (Suggests). `document()` populated
  NAMESPACE (incl. `importFrom(rlang,"!!")`, `importFrom(magrittr,"%>%")`, no
  exports yet); `load_all()` clean.
- Set up source isolation: froze a snapshot of the source project to
  `tmp_src/` (build- and git-ignored), removed the live source from the
  session's additional directories, repointed all docs at `tmp_src/`.
- Wrote `CLAUDE.md`.
- Reviewed `dev/plan.md` against the actual source and reconciled it: fixed the
  `load_protein` test, made `plot_loess_comparison` internal (ggplot/Suggests),
  added a rule to strip non-standard `@requires` roxygen tags, corrected the
  `library()` inventory, added the `R/data-doc.R` step, folded the resolved
  pre-flight questions into §13.
- Reorganized: `dev/plan.md` is now the single normative spec; deleted
  `dev/RESUME.md`; added `dev/PROGRESS.md`.
- `git init` on `main`. **Initial commit done** (`020ac2b`): skeleton +
  `dev/plan.md`, `dev/PROGRESS.md`, `CLAUDE.md`. `tmp_src/` and
  `.claude/settings.local.json` excluded.
- **Remote set up:** public repo `github.com/jechave/msamodel` (SSH origin),
  `main` pushed and tracking `origin/main`. Verified `tmp_src/` not pushed.

### 2026-04-22
- Created the package skeleton: `DESCRIPTION` (plan §3), `LICENSE` + `LICENSE.md`
  (MIT), empty roxygen-managed `NAMESPACE`, `.Rbuildignore`. `R/` left empty.
- Moved `plan.md` into `dev/`.
- Source-verification pass; answered the pre-flight questions.
