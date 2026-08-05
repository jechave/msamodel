# Frozen regression baseline for the parameter-uncertainty band, used by
# test-predict-uncertainty.R. These are the PRE-REFACTOR predict_*_ml band values (the
# code before the SPM-arm / `uncertainty` toggle work item, when the band was
# parameter-only), captured once from git HEAD at the start of that item. They are
# hard-coded LITERALS -- never recomputed by the code under test -- so that
# `uncertainty = "parameter"` reproducing them is a real regression lock: a change in the
# parameter arm makes the test go red. First six residues of the znb site fit
# (fit_lrmsd_msa_site(znb_spm, znb_profile$pdb_site, znb_profile$lrmsd_i_obs)).
PREREF <- list(
  lrmsd_lower  = c(-2.319570163, -3.291483151, -3.224439846, -3.541502590,
                   -3.646696482, -3.548340622),
  lrmsd_upper  = c(-2.280454884, -3.237787517, -3.154508664, -3.494742268,
                   -3.584651334, -3.474623011),
  nlrmsd_lower = c(1.725007407, 0.7753382911, 0.8556994599, 0.5045577431,
                   0.4256395535, 0.5306804704),
  nlrmsd_upper = c(1.908672292, 0.9290957884, 0.9990567765, 0.6929021462,
                   0.5767173766, 0.6800606434)
)
