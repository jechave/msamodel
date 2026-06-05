#' Load a protein structure file
#'
#' @param pdb_chain PDB chain identifier (e.g., "1btl_A")
#' @param data_dir Directory containing PDB files
#' @return Loaded PDB structure
#' @family setup
#' @export
load_protein <- function(pdb_chain, data_dir = "data/pdb_chain") {
  pdb_fn <- file.path(data_dir, paste0(pdb_chain, ".pdb"))
  read.pdb(pdb_fn)
}

#' Get active site residues for a protein from dataset_ec2024.csv
#'
#' @param pdb_chain PDB chain identifier
#' @param dataset Dataset containing active site information
#' @return List with pdb_site_active and site_active vectors
#' @family setup
#' @export
get_active_site <- function(pdb_chain, dataset) {
  # Filter for the specific pdb_chain and parse comma-separated strings
  asite <- dataset %>%
    filter(pdb_chain == !!pdb_chain)

  # Check if the pdb_chain exists
  if (nrow(asite) == 0) {
    stop("No active site data found for pdb_chain: ", pdb_chain)
  }
  if (nrow(asite) > 1) {
    warning("Multiple entries found for pdb_chain: ", pdb_chain, ". Using the first row.")
    asite <- asite[1, ]
  }

  # Extract vectors directly from the single row
  list(
    pdb_site_active = as.integer(str_split(asite$pdb_site_active, ",")[[1]]),
    site_active = as.integer(str_split(asite$site_active, ",")[[1]])
  )
}
