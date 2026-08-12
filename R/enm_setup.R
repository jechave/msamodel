#' Build a protein's elastic network model
#'
#' Represents the protein as beads joined by springs: one node per residue, a spring
#' between every pair of nodes closer than `d_max`. Diagonalising the resulting
#' spring-constant matrix gives the protein's normal modes.
#'
#' The model itself is penm's; this function validates the input and forwards to
#' [penm::set_enm()], which is the authority on what each parameter means.
#'
#' @param pdb A bio3d pdb object, e.g. from \code{bio3d::read.pdb()} (legacy
#'   \code{.pdb}) or \code{bio3d::read.cif()} (mmCIF).
#' @param node What plays the part of a residue: `"ca"` (alpha carbons), `"sc"`
#'   (side-chain centroids, the default) or `"cb"` (beta carbons). The aliases
#'   `"calpha"`, `"side_chain"` and `"beta"` also work; anything else is an error.
#' @param model Which spring rules to use: `"anm"`, `"ming_wall"` (the default),
#'   `"hnm"`, `"hnm0"`, `"pfanm"` or `"reach"`. They differ in how a contact's spring
#'   constant is set — see [penm::set_enm()].
#' @param d_max Contact cutoff in Ångström: node pairs closer than this get a spring.
#'   A workable value depends on `node` — penm's own examples pair `"ca"` with 10.5,
#'   `"cb"` with 12.0 and `"sc"` with 12.5.
#' @param frustrated Whether frustration is included when the spring-constant matrix
#'   is calculated. Defaults to `FALSE`; see [penm::set_enm()].
#' @return The protein, carrying its elastic network model: the nodes, the contact
#'   network and its spring-constant matrix, and the normal modes.
#' @seealso [penm::set_enm()] (does the work, and defines these parameters);
#'   [generate_spm()], which consumes the object returned here.
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
