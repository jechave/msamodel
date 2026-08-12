#' Set up Elastic Network Model (ENM) for a protein structure
#'
#' Builds the elastic network model on a bio3d pdb object via [penm::set_enm()],
#' returning the structure object augmented with the ENM (network, Hessian, normal
#' modes) that the SPM and divergence routines need. The package entry point for
#' turning a structure into a modellable object; does no file I/O.
#'
#' @param pdb A bio3d pdb object, e.g. from \code{bio3d::read.pdb()} (legacy
#'   \code{.pdb}) or \code{bio3d::read.cif()} (mmCIF). The package takes the
#'   structure object directly and does no file I/O, so any source that yields a
#'   bio3d pdb object works (PDB, AlphaFold-DB, etc.).
#' @param node Node representation for the network: `"ca"` (C-alpha atoms) or
#'   `"sc"` (side-chain centroids).
#' @param model Name of the elastic network model to build.
#' @param d_max Distance cutoff: pairs of nodes closer than this are connected.
#' @param frustrated Whether to use the frustrated form of the model.
#' @return The input structure augmented with the elastic network model (network,
#'   Hessian, and normal modes), ready for mutation scans and site properties.
#' @seealso [generate_spm_data()], which consumes the object returned here.
#' @family setup
#' @examples
#' \dontrun{
#' pdb <- bio3d::read.pdb(system.file("extdata", "1znb_A.pdb", package = "msamodel"))
#' wt  <- setup_enm(pdb, node = "ca", model = "ming_wall", d_max = 10.5)
#' }
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
