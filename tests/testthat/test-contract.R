# The observation site-key contract.
#
# A user supplies observations as two vectors, (pdb_site, lrmsd_obs). resolve_site_obs()
# maps pdb_site to the model's internal 1-based site index via spm$site_map.
#
# WHY THIS IS WORTH TESTING when the drift pins already cover the verbs: the pins
# exercise ONE input -- the shipped profile, in file order, covering every site. They
# cannot see behaviour that depends on the shape of the input. If the join were replaced
# by a positional assumption (taking lrmsd_obs to be in site order), the pinned snapshot
# would move and you would investigate -- but only because the shipped profile happens to
# be complete and ordered. Hand a scrambled or partial profile to a positional
# implementation and it returns plausible, wrong numbers at every site, silently.
#
# On this example protein pdb_site runs 1..107, so pdb_site == site and the two keys are
# interchangeable BY ACCIDENT. These tests construct the distinguishing cases on purpose
# -- scrambled order, strict subset -- rather than relying on a PDB whose numbering
# happens to start at 20.

test_that("a scrambled, partial observed profile resolves by key, not by position", {
  # 40 of 107 sites in shuffled order. Under a positional implementation the returned
  # obs would follow input order; under the key join each obs must land on ITS OWN site.
  set.seed(1)
  take <- sample(nrow(obs_site), 40)
  scrambled <- obs_site[take, ]

  got <- msamodel:::resolve_site_obs(spm, scrambled$pdb_site, scrambled$lrmsd_obs)

  expect_equal(nrow(got), 40L)
  # Every returned (idx, obs) pair must carry the value the caller supplied for that
  # pdb_site. Negative control: taking obs in input order (positional) breaks this,
  # because the join reorders rows.
  expected <- scrambled$lrmsd_obs[match(got$idx, scrambled$pdb_site)]
  expect_equal(got$obs, expected)
})

test_that("a partial profile still yields a finite fit", {
  # The model spans more sites than the data measures -- the ordinary case for real
  # alignment data, which rarely observes every modelled residue.
  subset_obs <- obs_site[1:60, ]
  expect_lt(nrow(subset_obs), nrow(spm$site_map))   # genuinely a subset

  ml <- fit_lrmsd_msa_site(spm, subset_obs$pdb_site, subset_obs$lrmsd_obs)
  expect_true(is.finite(ml$a1) && is.finite(ml$a2) && is.finite(ml$logLik))
  # nobs is the MATCHED count, not the model's full support. With a strict subset these
  # two differ, so this cannot pass by coincidence.
  expect_equal(ml$gof$nobs, 60L)
  expect_lt(ml$gof$nobs, nrow(spm$site_map))
})

test_that("an unknown pdb_site is an error, not a silent drop", {
  bad <- obs_site$pdb_site
  bad[1] <- 999999L
  expect_error(fit_lrmsd_msa_site(spm, bad, obs_site$lrmsd_obs),
               "not present in the model")
  # Control: the same call with valid keys does not error, so the guard is not blanket.
  expect_no_error(fit_lrmsd_msa_site(spm, obs_site$pdb_site, obs_site$lrmsd_obs))
})

test_that("an unknown mode is an error on the mode axis too", {
  bad <- obs_mode$mode
  bad[1] <- 999999L
  expect_error(fit_lrmsd_msa_mode(spm, bad, obs_mode$lrmsd_obs),
               "not present in the model")
})
