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

## Flatten the error-band plumbing

The fitters were flattened 2026-08-05 (`f504a8f`); the band machinery
(`R/predict_band_helpers.R` + the band code in `R/api.R`) still has the shape they had —
anonymous closures over `theta` built per column and handed to `var_param_delta()`, with
the forward map re-evaluated each time. Probably the same treatment.

Coupled to the two entries above on the terse band-helper names and the band/`_se`
vocabulary; structure first, then naming.

## Are the six shipped `znb_*` datasets right?

Six exported datasets, six man pages. Is that more surface than the package needs, or
about right? `znb_spm` (15.5 MB) and `znb_wt` (6.3 MB) are 22 of the 25.8 MB installed
size and the reason for the standing `check()` WARNING + size NOTE.

Not a rename job: dropping an exported dataset is a user-visible break, and `znb_spm`
in particular saves users a minutes-long scan.

## Joint site + mode fit

One `(a1, a2)` estimated from both observed profiles at once, rather than two independent
fits. Cheap to wire (the likelihoods are additive), but it is a new estimator, not
repackaging: each axis profiles out its own `sigma`, the two profiles are not independent
(both are projections of the same SPM scan, so joint SEs would be optimistic), and `D2`
loses its meaning over two incommensurable supports.

Undecided and worth settling first: is the goal a better single estimate, or a
consistency check between the two profiles? Those pull in opposite directions.

## `add_site_properties()`

Does this belong in `msamodel` at all? If it stays: it joins on the internal `site`
index where `pdb_site` would be unambiguous, and its `left_join` silently yields `NA`
properties for an unmatched row instead of failing loud.

## `dev/tmp2/` scratch

~14 MB of untracked throwaway scratch. Delete when convenient.

## v0.4 — motion arm

The motion/fluctuation quantities (`dh_ijm`, `dh_njm`, `nh_njm`) are unimplemented;
see `dev/findings.md` for the grid and why the existing reweighting machinery should
extend to them. Available in penm's slow per-mutant loop only. Not scoped.
