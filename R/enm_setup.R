#' Set up Elastic Network Model (ENM) for a protein structure
#'
#' @param pdb A bio3d pdb object, e.g. from \code{bio3d::read.pdb()} (legacy
#'   \code{.pdb}) or \code{bio3d::read.cif()} (mmCIF). The package takes the
#'   structure object directly and does no file I/O, so any source that yields a
#'   bio3d pdb object works (PDB, AlphaFold-DB, etc.).
#' @param node Node type for ENM ("ca" or "sc")
#' @param model ENM model type
#' @param d_max Maximum distance cutoff
#' @param frustrated Whether to use frustrated model
#' @return PDB structure with ENM setup
#' @family setup
#' @export
setup_enm <- function(pdb, node = "sc", model = "ming_wall",
                      d_max = 10.5, frustrated = FALSE) {
  if (!inherits(pdb, "pdb")) {
    stop("`pdb` must be a bio3d pdb object (class \"pdb\"), e.g. from ",
         "bio3d::read.pdb() or bio3d::read.cif(); got an object of class ",
         paste(class(pdb), collapse = "/"), ".")
  }
  pdb %>% set_enm(node = node, model = model,
                  d_max = d_max, frustrated = frustrated)
}
