#!/usr/bin/env Rscript
# Preview a pre-rendered vignette as standalone HTML WITHOUT touching the repo.
#
# Why this exists: `html_vignette` renders with pandoc `--self-contained`, which
# base64-embeds the figures into the HTML and then deletes the intermediate
# `<name>_files/` directory. Rendering a committed vignette *in place* therefore
# DESTROYS the figures the shipped `<name>.Rmd` references. (We hit this once and
# had to re-knit to recover them.)
#
# This script sidesteps that: it copies the already-knit `<name>.Rmd` and its
# `<name>_files/` figures into a throwaway temp dir, renders there, and copies
# only the resulting self-contained `.html` to `dev/preview/` (a scratch dir,
# outside the package source). The real `vignettes/<name>.Rmd` and
# `<name>_files/` are never written to, so the figures the package ships are
# safe, and `vignettes/` holds only files that actually ship.
#
# It does NOT knit the .Rmd.orig -> .Rmd (that is the separate, deliberate
# precompute step). Run knit first if you changed the .orig; then preview.
#
# Usage (from the package root):
#   Rscript dev/preview-vignette.R mode-analysis
#   Rscript dev/preview-vignette.R              # previews every vignettes/*.Rmd
#
# Output: dev/preview/<name>.html  (git-ignored; for local viewing only).

args <- commandArgs(trailingOnly = TRUE)

vig_dir <- "vignettes"
if (!dir.exists(vig_dir)) {
  stop("Run this from the package root (no ", vig_dir, "/ directory here).")
}

# Scratch output dir for previews (git-ignored; .Rbuildignore already excludes
# all of dev/). Kept out of vignettes/ so the source dir holds only files that
# ship.
preview_dir <- file.path("dev", "preview")
dir.create(preview_dir, showWarnings = FALSE, recursive = TRUE)

# Resolve which vignettes to preview.
if (length(args) >= 1) {
  names <- sub("\\.Rmd$", "", basename(args))
} else {
  rmds <- list.files(vig_dir, pattern = "\\.Rmd$", full.names = FALSE)
  names <- sub("\\.Rmd$", "", rmds)
}
if (length(names) == 0) stop("No vignettes found to preview.")

preview_one <- function(name) {
  rmd <- file.path(vig_dir, paste0(name, ".Rmd"))
  if (!file.exists(rmd)) stop("No such rendered vignette: ", rmd,
                              " (knit the .Rmd.orig first).")
  figdir <- file.path(vig_dir, paste0(name, "_files"))

  # Scratch dir: render the COPY here so the real vignettes/ is untouched.
  scratch <- file.path(tempdir(), paste0("preview-", name))
  unlink(scratch, recursive = TRUE)
  dir.create(scratch, recursive = TRUE)
  file.copy(rmd, file.path(scratch, basename(rmd)))
  if (dir.exists(figdir)) {
    file.copy(figdir, scratch, recursive = TRUE)
  }

  out_html <- file.path(scratch, paste0(name, ".html"))
  rmarkdown::render(
    input = file.path(scratch, basename(rmd)),
    output_format = "rmarkdown::html_vignette",
    output_file = basename(out_html),
    quiet = TRUE
  )

  # Self-contained HTML -> scratch dir under dev/ (NOT vignettes/, which holds
  # only files that ship).
  dest <- file.path(preview_dir, paste0(name, ".html"))
  file.copy(out_html, dest, overwrite = TRUE)
  message("Preview written: ", dest,
          "  (repo vignettes/", name, "_files/ untouched)")
}

invisible(lapply(names, preview_one))
