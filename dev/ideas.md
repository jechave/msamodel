# Ideas / maybe later

A parking lot for things worth doing eventually. **Not a plan and not a status
file** — nothing here is scheduled, ordered, or claimed to be current. Not read at
session start; open it when you want to pick something up.

Delete an entry when it is done, or when it stops being a good idea.

---

## Terse band-machinery names

`spm_hmat`, `grad_theta`, `var_param_delta`, `as_theta` (all in
`R/predict_band_helpers.R`, some also used in `R/api.R`) are named in a terse style
the author rejects — as are the internal `M`/`G`/`Gc`/`h`/`w` inside `spm_hmat`.

The `t` → `theta` coordinate piece of this was done (2026-08-04); the rest was not.

Coupled to the notation in `dev/reports/spm_band_validation.Rmd`, which calls
`spm_hmat` by name in 7 places under `load_all()` and mirrors the terse notation
throughout. Rename the code first so the validated behaviour is preserved, then
follow the PDF notation to match.

## "band" → "se" vocabulary

The uncertainty vocabulary is inconsistent: `R/` carries ~28 uses of
band/bands/banded/banding, while the public API returns `_se` columns. Worth
settling on one word. Touches `R/predict_band_helpers.R` (including its filename),
roxygen, and the validation report.

## Rethink the `R/` family / file structure

The `CLAUDE.md` rule fixing how `R/` files map to `@family` tags was **removed
2026-08-05** pending a rethink, deliberately not replaced with a corrected version.

Two concrete inputs for that thinking:

- `api` had become a seventh family the rule never listed (`R/api.R:101,145,198,252`)
  — and it holds the main public surface.
- `R/msa_decomposition.R` carries no `@family` tag at all (`decompose_nested`,
  `R/msa_decomposition.R:28`).

Fixing the missing tag is an `R/` + roxygen change, so it needs `document()` and the
full-suite gate; it was deliberately left out of the docs-only cleanup.

## `predict_decomposition` errors on its own default

`predict_decomposition(fit, spm)` cannot succeed as written: `match.arg` defaults
`which` to `"lrmsd"`, and that branch is an explicit `stop()` saying the uncentred
component band is not yet derived (`R/api.R:263`). Either derive the uncentred
component SE, or make `"nlrmsd"` the default so a bare call works.

## Joint site + mode fit

A single `(a1, a2)` estimated from BOTH observed profiles at once, rather than two
independent fits. Mechanically cheap — the likelihoods are additive and
`fit_lrmsd_msa()` already takes a resolved `(idx, obs)` frame — but it is a new
estimator, not repackaging, and three questions have to be answered first:

- **Sigma.** Each axis profiles out its own `sigma`. Summing the two profiled
  log-likelihoods lets each self-normalise, so their relative weight is set implicitly
  by their residual scales. Defensible, but then `k` is 4, not 3, and everything using
  `k` (AIC, BIC) is wrong if it stays hardcoded.
- **Independence.** Summing log-likelihoods asserts the two profiles are independent.
  They are not — `dr2_ijm` and `dr2_njm` are two projections of the same SPM scan — so
  the joint SEs would be optimistic by an unknown amount.
- **The GoF row.** `nobs` becomes sites + modes and `null_deviance` a sum over two
  incommensurable supports, so `D2` loses its "fraction of profile variance explained"
  reading.

Worth deciding first whether the goal is a *better single estimate* or a *consistency
check between profiles* — those pull in opposite directions (pool vs. fit separately and
compare). The interface question was settled in discussion: it would be its own verb, not
a flag on the existing fitters.

## `add_site_properties()` — join key, and whether it belongs here

Two separate doubts, both raised 2026-08-05 and deliberately deferred:

- It joins on the internal `site` index; `pdb_site` would be the unambiguous key. That
  needs `get_pdb_site(wt)` instead of `get_site(wt)`, so it changes the function's
  contract, not just the `by=`.
- It uses `left_join`, so a `site_data` row with no match silently gets `NA` for
  `dactive`/`cn`/`msf`/`shell` rather than erroring — the quiet-fallback shape
  `CLAUDE.md` asks to fail loud on.

Underneath both: the user is unsure this function belongs in `msamodel` at all.

## `dev/tmp2/` scratch

~14 MB of untracked throwaway scratch. Delete when convenient.

## v0.4 — motion arm

The motion/fluctuation quantities (`dh_ijm`, `dh_njm`, `nh_njm`) are unimplemented;
see `dev/findings.md` for the grid and why the existing reweighting machinery should
extend to them. Available in penm's slow per-mutant loop only. Not scoped.
