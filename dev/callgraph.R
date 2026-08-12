#!/usr/bin/env Rscript
# dev/callgraph.R — regenerate the package's function call graph.
#
#   Rscript dev/callgraph.R
#
# Produces (git-ignored, in dev/preview/):
#   callgraph.pdf   — a 3-page vector document, one PIPELINE STAGE per page, each sized
#                     to print legibly on a single Letter sheet (landscape):
#                       p1  build & fit   (setup -> spm -> fit -> gof)
#                       p2  calculate     (the bare-(a1,a2) forward-map verbs)
#                       p3  predict       (the fit-driven banded verbs + band machinery)
#
# HOW IT WORKS (no hand-drawn edges): mvbutils::foodweb does static analysis of the
# loaded package's function bodies to get the real caller -> callee edges; those are
# written to Graphviz DOT and rendered with `dot`. Every arrow means "A calls B". Each
# page is an induced subgraph: a stage's entry-point functions plus their transitive
# callees, so the shared primitive layer each stage reaches is shown (and repeats, grey,
# across pages). Re-run after any change to R/ to refresh the picture.
#
# STRUCTURE THIS SHOWS (post old-grid retirement): the exported surface is the four
# axis-agnostic leaf verbs (calculate_profiles / predict_profiles / calculate_decomposition
# / predict_decomposition) + fit_*_ml + gof_*_ml + the setup/spm helpers. They fan down
# into the shared @noRd forward-map primitives (dr2_msa, lrmsd_msa, nlrmsd_msa, the
# nested_models / decomposition primitives, weights_jm) and the se machinery in
# R/predict_se.R (se_profile_* / se_nested_* / se_components_*, over var_param_* and
# var_spm_*, on grad_theta + spm_hmat). Nodes are coloured by ROLE
# (setup / spm / model-primitive / fitting / prediction / api), NOT by site/mode axis —
# the axis-specific grid was deleted, so there is no longer a site arm vs a mode arm.
#
# PRINT-LEGIBLE by design: high-contrast (black text, black arrows, saturated fills),
# each page scaled to fill a Letter sheet at a readable font. dashed border = @noRd
# internal; solid = exported.
#
# Requires: mvbutils + qpdf (R pkgs) and graphviz's `dot` on PATH
# (`brew install graphviz`).

suppressMessages({
  library(mvbutils)
})

# --- locate the package root (works whether run from root or dev/) -----------------
here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
root <- if (!is.null(here)) normalizePath(file.path(here, "..")) else getwd()
if (!file.exists(file.path(root, "DESCRIPTION"))) root <- getwd()

# --- fail loud if a dependency is missing (don't silently produce nothing) ----------
if (nzchar(Sys.which("dot")) == FALSE) {
  stop("graphviz `dot` not found on PATH. Install it first: brew install graphviz")
}
if (!requireNamespace("qpdf", quietly = TRUE)) {
  stop("R package `qpdf` not installed (needed to merge the pages). install.packages('qpdf')")
}

# --- load the package under development ---------------------------------------------
suppressMessages(devtools::load_all(root, quiet = TRUE))
ns   <- getNamespace("msamodel")
objs <- ls(ns)
fns  <- objs[vapply(objs, function(o) is.function(get(o, envir = ns)), logical(1))]

# --- static call analysis -> caller x callee matrix (m[a, b] == 1 means a calls b) --
# NB: mvbutils::foodweb opens a graphics device even with plotting = FALSE, which
# dumps a stray Rplots.pdf in the working dir. Route it to a temp file and clean up.
.tmpdev <- tempfile(fileext = ".pdf")
grDevices::pdf(.tmpdev)
fw <- foodweb(where = ns, funs = fns, plotting = FALSE)
grDevices::dev.off()
unlink(.tmpdev)
m  <- fw$funmat
nm <- rownames(m)

# --- classify each node by ROLE (the layer it belongs to) ---------------------------
# The role decides the fill colour. Exported-vs-internal decides the border (solid vs
# dashed). Roles are matched by name pattern against the current function inventory;
# anything unmatched falls to "other" (magenta) so a new function is loud, not hidden.
exported <- getNamespaceExports("msamodel")

role_of <- function(x) {
  if (x == "set_enm") return("setup")
  if (x %in% c("generate_spm", "generate_spm_core")) return("spm")
  if (grepl("^predict_", x) || grepl("^se_", x) ||
      grepl("^var_param_", x) || grepl("^var_spm_", x) ||
      x %in% c("grad_theta", "spm_hmat",
               "key_profile", "validate_ml_fit")) return("prediction")
  if (grepl("^fit_", x) || grepl("loglik", x) ||
      x %in% c("residuals_lrmsd_msa", "gof_lrmsd_msa", "check_obs_vectors",
               "calculate_null_deviance",
               "resolve_site_obs", "resolve_mode_obs")) return("fitting")
  if (x %in% c("calculate_profiles", "calculate_decomposition")) return("api")
  # everything else = the model / forward-map primitive layer
  if (x %in% c("pfix_msa", "weights_jm", "weights_jm_spm", "dr2_msa", "lrmsd_msa",
               "nlrmsd_msa", "lrmsd_nested_models", "nlrmsd_nested_models",
               "lrmsd_msa_decomposition", "nlrmsd_msa_decomposition",
               "decompose_nested", "axis_branches", "delta_structure_dr")) return("model")
  "other"
}
role <- vapply(nm, role_of, character(1))

# --- high-contrast, print-first palette, keyed by ROLE ------------------------------
# Saturated fills, BLACK node text and BLACK arrows so it survives a b/w laser and holds
# contrast in colour. dashed border = @noRd internal; solid = exported.
ROLE_COL <- list(
  setup      = list(fill = "#cfe8cf", border = "#245c24"),  # green
  spm        = list(fill = "#cfe0f2", border = "#1f4e79"),  # blue
  api        = list(fill = "#f6cfe0", border = "#8a1f52"),  # pink  (calculate_* leaf verbs)
  prediction = list(fill = "#f2dca6", border = "#7a4e00"),  # gold  (predict_* + band machinery)
  fitting    = list(fill = "#e3d4f2", border = "#4a237a"),  # purple
  model      = list(fill = "#d9d7d2", border = "#3a3f4b"),  # grey  (shared primitive spine)
  other      = list(fill = "#ff66cc", border = "#000000")   # loud: an unclassified function
)

# One node-declaration line, coloured by role.
node_line <- function(x, indent = "  ") {
  is_exp <- x %in% exported
  col    <- ROLE_COL[[role[[x]]]]
  sprintf(paste0("%s\"%s\" [shape=box, style=\"filled,%s\", color=\"%s\", penwidth=1.6, ",
                 "fillcolor=\"%s\", fontcolor=\"black\"];"),
          indent, x, if (is_exp) "solid" else "dashed", col$border, col$fill)
}

# Transitive closure: all functions reachable from `roots` via the call matrix.
reachable <- function(roots) {
  seen <- character(0)
  frontier <- intersect(roots, nm)
  while (length(frontier)) {
    seen <- union(seen, frontier)
    nxt  <- character(0)
    for (a in frontier) nxt <- union(nxt, colnames(m)[m[a, ] > 0])
    frontier <- setdiff(nxt, seen)
  }
  seen
}

# Induced-subgraph DOT for one pipeline stage, sized to fill a Letter landscape sheet.
stage_dot <- function(roots, title) {
  keep <- reachable(roots)
  L <- c(
    "digraph stage {",
    # ratio=fill + explicit size FORCES the drawing to the box's proportions so nodes and
    # text are large on paper (these graphs are wide-and-shallow; left alone dot renders a
    # thin strip). nodesep/ranksep give the fill room to work.
    "  size=\"10,7.5\"; ratio=fill; margin=0;",
    "  labelloc=\"t\";",
    sprintf("  label=\"%s\";", title),
    "  fontname=\"Helvetica-Bold\"; fontsize=16;",
    "  rankdir=TB;",
    "  graph [splines=true, nodesep=0.5, ranksep=1.0];",
    "  node  [fontname=\"Courier-Bold\", fontsize=15, margin=\"0.18,0.10\"];",
    "  edge  [color=\"black\", penwidth=1.2, arrowsize=0.9];"
  )
  for (x in keep) L <- c(L, node_line(x))
  for (a in keep) for (b in intersect(colnames(m)[m[a, ] > 0], keep)) {
    L <- c(L, sprintf("  \"%s\" -> \"%s\";", a, b))
  }
  c(L, "}")
}

# --- the three pipeline stages (entry points; callees follow transitively) ----------
build_fit_roots <- c("set_enm", "generate_spm",
                     "fit_lrmsd_msa_site", "fit_lrmsd_msa_mode")
calculate_roots <- c("calculate_profiles", "calculate_decomposition")
predict_roots   <- c("predict_profiles", "predict_decomposition")

# --- render each page to its own temp PDF, then merge --------------------------------
outdir <- file.path(root, "dev", "preview")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

render_pdf <- function(dot_lines) {
  dotfile <- tempfile(fileext = ".dot")
  pdffile <- tempfile(fileext = ".pdf")
  writeLines(dot_lines, dotfile)
  ok <- system2("dot", c("-Tpdf", shQuote(dotfile), "-o", shQuote(pdffile)))
  unlink(dotfile)
  if (ok != 0) stop("dot failed to render a page")
  pdffile
}

legend <- "roles: green=setup  blue=spm  pink=calculate(api)  gold=predict  purple=fitting  grey=model-primitive  |  solid=exported  dashed=@noRd  |  arrow: A calls B"
pages <- c(
  render_pdf(stage_dot(build_fit_roots,
                       sprintf("msamodel — build & fit stage (setup -> spm -> fit -> gof)\\n%s", legend))),
  render_pdf(stage_dot(calculate_roots,
                       sprintf("msamodel — calculate stage (bare-(a1,a2) forward-map verbs)\\n%s", legend))),
  render_pdf(stage_dot(predict_roots,
                       sprintf("msamodel — predict stage (fit-driven banded verbs + band machinery)\\n%s", legend)))
)
on.exit(unlink(pages), add = TRUE)

pdf_out <- file.path(outdir, "callgraph.pdf")
qpdf::pdf_combine(pages, output = pdf_out)

cat(sprintf("call graph: %d functions, %d edges — 3-page print-legible PDF (by pipeline stage)\n  %s\n",
            length(nm), sum(m > 0), pdf_out))

# Loud check: any function that fell through role classification (would render magenta).
unclassified <- nm[role == "other"]
if (length(unclassified)) {
  cat("NOTE: unclassified functions (render magenta — assign a role in role_of()):\n  ",
      paste(unclassified, collapse = ", "), "\n")
}
