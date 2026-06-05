#' Set up Elastic Network Model (ENM) for a protein structure
#'
#' @param pdb PDB structure
#' @param node Node type for ENM ("ca" or "sc")
#' @param model ENM model type
#' @param d_max Maximum distance cutoff
#' @param frustrated Whether to use frustrated model
#' @return PDB structure with ENM setup
#' @family setup
#' @export
setup_enm <- function(pdb, node = "sc", model = "ming_wall",
                      d_max = 10.5, frustrated = FALSE) {
  pdb %>% set_enm(node = node, model = model,
                  d_max = d_max, frustrated = frustrated)
}
