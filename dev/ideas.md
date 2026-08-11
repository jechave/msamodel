# Ideas / maybe later

A parking lot for things worth doing eventually. **Not a plan and not a status
file** — nothing here is scheduled, ordered, or claimed to be current. Not read at
session start; open it when you want to pick something up.

Delete an entry when it is done, or when it stops being a good idea.

---

## Terse names in the se machinery

`spm_hmat` and `grad_theta` (`R/predict_se.R`) are named in a terse style the author
rejects — as are the internal `M`/`G`/`Gc`/`h`/`w` inside `spm_hmat`. Both were kept
deliberately when that file was flattened: `spm_hmat` because
`dev/reports/spm_band_validation.Rmd` calls it by name in 7 places under `load_all()`,
`grad_theta` because it names real math. (`var_param_delta` and `as_theta` are gone.)

Rename the code first so the validated behaviour is preserved, then follow the PDF
notation to match.

## "band" → "se" vocabulary in the validation report

Purged from `R/` when the se machinery was flattened. `dev/reports/spm_band_validation.Rmd`
still uses band/bands throughout, including in its filename.

## Rethink the `R/` family / file structure

The `CLAUDE.md` rule fixing how `R/` files map to `@family` tags was **removed
2026-08-05** pending a rethink, deliberately not replaced with a corrected version.

Two concrete inputs for that thinking:

- `api` is a seventh family the rule never listed, and it holds the main public surface:
  the four verbs now split across `R/model.R` (`calculate_*`) and `R/predict.R`
  (`predict_*`), all still tagged `@family api`.
- `R/msa_decomposition.R` carries no `@family` tag at all (`decompose_nested`,
  `R/msa_decomposition.R:28`).

Fixing the missing tag is an `R/` + roxygen change, so it needs `document()` and the
full-suite gate; it was deliberately left out of the docs-only cleanup.

## Derive the uncentred (lrmsd) nested-model and component SEs

`se_nested_lrmsd()` and `se_components_lrmsd()` (`R/predict_se.R`) are stubs that
`stop("... to be developed")`, so `predict_decomposition(metric = "lrmsd")` errors. The
maths was never derived. The stubs are where it goes.

## Flatten `calculate_profiles` / `calculate_decomposition` too

The `predict_*` verbs were flattened and made explicit (both axes written out, no
`lapply`, no closures hiding a calculation). The `calculate_*` verbs were deliberately
moved byte-for-byte in that same work and still have the old shape — `lapply` over
`axis_branches()`. Same treatment, when there is appetite for it.

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

## v0.4 — motion arm

The motion/fluctuation quantities (`dh_ijm`, `dh_njm`, `nh_njm`) are unimplemented;
see `dev/findings.md` for the grid and why the existing reweighting machinery should
extend to them. Available in penm's slow per-mutant loop only. Not scoped.
