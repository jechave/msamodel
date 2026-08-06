#!/usr/bin/env Rscript
#
# ############################################################################
# #  THIS SCRIPT INSTALLS THE PACKAGE. It runs devtools::install(), which     #
# #  OVERWRITES the msamodel in your R library with the current working       #
# #  tree. Every other R session and project that calls library(msamodel)     #
# #  will get this build. If you had a released version installed on          #
# #  purpose, it is gone.                                                     #
# #                                                                           #
# #  It also DELETES vignettes/<name>_files/ before knitting.                  #
# ############################################################################
#
# Knit vignettes/<name>.Rmd.orig -> vignettes/<name>.Rmd, guaranteeing that
# EVERYTHING you then look at was produced by THIS run: the printed numbers, the
# figures, and the code that computed them.
#
# Three guarantees, and the reason each is needed:
#
#   1. INSTALL FIRST. Vignette chunks call library(msamodel), which loads the
#      INSTALLED package, not your working tree. Knitting without installing
#      silently recomputes everything with the OLD code and stamps fresh
#      timestamps on it -- an inconsistency that looks perfectly healthy.
#
#   2. WIPE THE FIGURE DIR. Figures live as PNG files that the shipped .Rmd only
#      REFERENCES. Nothing ties a PNG to the knit that wrote the .Rmd text, so a
#      leftover figure is silently reused: new numbers, old plots, no error.
#      Deleting the directory first makes a stale figure impossible rather than
#      merely unlikely.
#
#   3. RUN FROM THE PACKAGE ROOT, or stop. knitr writes figures relative to the
#      WORKING DIRECTORY. On 2026-08-06 a knit run from the package root wrote
#      every figure to <root>/<name>_files/ and left vignettes/<name>_files/ a
#      week stale; the previews embedded old images and nothing errored. This
#      script therefore refuses to run anywhere but the root, and sets the cwd to
#      vignettes/ itself for the knit. Deriving the root from the script's own
#      location was tried and rejected: it is more machinery, and it has sharp
#      edges (R encodes spaces in --file= as "~+~", which this package's path
#      hits). Failing loudly on the wrong cwd is more robust than being clever
#      about finding the right one.
#
# Verification is on the ARTIFACT, not on this script's bookkeeping: after
# knitting, the figure paths are read OUT OF the knitted .Rmd, and each one must
# exist and postdate the start of this run. Anything missing or stale is an error.
#
# Cost, accepted deliberately: a full install plus a full recompute on every run.
# Nothing is reused. That is the point.
#
# Usage (MUST be run from the package root):
#   Rscript dev/knit-vignette.R mode-analysis
#   Rscript dev/knit-vignette.R                 # every vignettes/*.Rmd.orig
#
# Does NOT render HTML. Preview after this:
#   Rscript dev/preview-vignette.R <name>

args <- commandArgs(trailingOnly = TRUE)

# ---- 1. refuse to run anywhere but the package root ------------------------
# Deliberately the least clever option available: no path derivation, just a
# check that fails loudly. Getting the cwd wrong is the bug this script exists to
# prevent, so it must not be silently survivable.
if (!file.exists("DESCRIPTION") || !dir.exists("vignettes")) {
  stop("Run this from the package root -- the directory holding DESCRIPTION and ",
       "vignettes/.\n  Currently in: ", getwd(),
       "\n  Usage: Rscript dev/knit-vignette.R <name>", call. = FALSE)
}
root    <- normalizePath(getwd(), mustWork = TRUE)
vig_dir <- file.path(root, "vignettes")

message("package root : ", root)

# ---- 2. which vignettes ----------------------------------------------------
if (length(args) >= 1) {
  names <- sub("\\.Rmd(\\.orig)?$", "", basename(args))
} else {
  names <- sub("\\.Rmd\\.orig$", "", list.files(vig_dir, pattern = "\\.Rmd\\.orig$"))
}
if (length(names) == 0) stop("No vignettes/*.Rmd.orig found.")
for (name in names) {
  if (!file.exists(file.path(vig_dir, paste0(name, ".Rmd.orig")))) {
    stop("No such source: vignettes/", name, ".Rmd.orig")
  }
}
message("vignettes    : ", paste(names, collapse = ", "))

# ---- 3. install the working tree so library(msamodel) sees current code -----
message("\n>>> installing the working tree (this OVERWRITES your installed msamodel) ...")
devtools::install(root, quiet = TRUE, upgrade = "never")
message(">>> installed: msamodel ", as.character(utils::packageVersion("msamodel")))

# Everything below must postdate this instant to count as produced by this run.
run_started <- Sys.time()

# ---- 4. wipe figures, then knit with the cwd inside vignettes/ --------------
for (name in names) {
  fig_dir <- file.path(vig_dir, paste0(name, "_files"))
  if (dir.exists(fig_dir)) {
    unlink(fig_dir, recursive = TRUE)
    message("wiped: vignettes/", name, "_files/")
  }
}

knit_one <- function(name) {
  old <- setwd(vig_dir)                 # cwd inside vignettes/: fig.path resolves here
  on.exit(setwd(old), add = TRUE)
  knitr::knit(paste0(name, ".Rmd.orig"), output = paste0(name, ".Rmd"), quiet = TRUE)
  message("knitted: vignettes/", name, ".Rmd")
}
invisible(lapply(names, knit_one))

# ---- 5. verify the ARTIFACT: every figure the .Rmd references is from this run
# Read the references out of the knitted .Rmd rather than listing a directory we
# assume is the right one -- that assumption is what failed on 2026-08-06.
cat("\n================ FIGURES REFERENCED BY THE KNITTED .Rmd ================\n")
problems <- character()

for (name in names) {
  rmd   <- file.path(vig_dir, paste0(name, ".Rmd"))
  lines <- readLines(rmd, warn = FALSE)
  refs  <- unique(unlist(regmatches(
    lines, gregexpr("[A-Za-z0-9._-]+_files/figure-html/[A-Za-z0-9._-]+\\.(png|jpg|jpeg|svg|pdf)",
                    lines))))

  cat("\n", name, ":\n", sep = "")
  if (length(refs) == 0) {
    cat("   (references no figures)\n")
    next
  }
  for (ref in refs) {
    f <- file.path(vig_dir, ref)
    if (!file.exists(f)) {
      cat(sprintf("   %-52s *** MISSING ***\n", ref))
      problems <- c(problems, paste0(name, ": missing ", ref))
    } else {
      m <- file.info(f)$mtime
      fresh <- m >= run_started
      cat(sprintf("   %-52s %s  %s\n", ref, format(m),
                  if (fresh) "from this run" else "*** STALE ***"))
      if (!fresh) problems <- c(problems, paste0(name, ": stale ", ref))
    }
  }
}

cat("\n=======================================================================\n")
if (length(problems)) {
  cat("\n")
  stop("Figure verification FAILED:\n  ", paste(problems, collapse = "\n  "),
       "\n\nThe knitted .Rmd references figures that are missing or predate this run. ",
       "Do NOT preview or commit this state.", call. = FALSE)
}
cat("All referenced figures were produced by this run, against a freshly\n")
cat("installed working tree. The .Rmd and its figures are consistent.\n")
cat("\nNext: Rscript dev/preview-vignette.R <name>   (then the USER approves the HTML)\n")
