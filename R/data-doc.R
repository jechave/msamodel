#' Embedded datasets for the example protein 1znb_A
#'
#' The `znb_*` datasets provide a complete, self-contained example (zinc
#' beta-lactamase, PDB chain 1znb_A) for the MSA-model workflow. They are built by
#' \code{data-raw/prepare_znb_data.R}.
#'
#' @name msamodel-data
#' @keywords internal
NULL

#' Raw PDB structure for 1znb_A
#'
#' The bio3d structure object for the example protein, as read from
#' \code{inst/extdata/1znb_A.pdb}.
#'
#' @format A bio3d \code{pdb} object (a list with components \code{atom},
#'   \code{xyz}, \code{helix}, \code{sheet}, \code{calpha}, \code{call}).
#' @source RCSB PDB entry 1ZNB, chain A. See \code{data-raw/prepare_znb_data.R}.
#' @family datasets
"znb_pdb"

#' Elastic network model (wild type) for 1znb_A
#'
#' The penm ENM wild-type object for the example protein, used as the reference
#' structure for mutation scans and as the source of site properties.
#'
#' @format A penm ENM object (class \code{prot}): a list with components
#'   \code{param}, \code{nodes}, \code{graph}, \code{eij}, \code{kmat},
#'   \code{nma}.
#' @details Built with \code{setup_enm(znb_pdb, node = "ca", model = "ming_wall",
#'   d_max = 10.5, frustrated = FALSE)}.
#' @source Derived from \code{znb_pdb}. See \code{data-raw/prepare_znb_data.R}.
#' @seealso \code{\link{setup_enm}}
#' @family datasets
"znb_wt"

#' Single-point-mutation scan (SPM) for 1znb_A
#'
#' Per-mutation energetic and structural-response data from a single-point
#' mutation scan of every site of the example protein. This is the primary input
#' to the divergence-profile and fitting functions (via \code{preprocess_spm}).
#'
#' @format A tibble with 2508 rows (228 sites times 11 states, mutation index
#'   m = 0..10) and 10 columns:
#' \describe{
#'   \item{j}{Integer index of the mutated (perturbed) site.}
#'   \item{m}{Integer mutation replicate (0 = wild type).}
#'   \item{ddg_dv_jm}{Stability change, energy term.}
#'   \item{ddg_tds_jm}{Stability change, entropy (-T*dS) term.}
#'   \item{ddgact_dv_jm}{Activity change, energy term.}
#'   \item{ddgact_tds_jm}{Activity change, entropy term.}
#'   \item{site}{List column: vector of internal site indices.}
#'   \item{pdb_site}{List column: vector of PDB residue numbers.}
#'   \item{dr}{List column: Cartesian displacement vector (length 3 * n_sites).}
#'   \item{dr2}{List column: squared displacement per site (length n_sites).}
#' }
#' @details Generated (not copied) by \code{generate_spm_data(znb_wt,
#'   n_mutations = 10, model = "lfenm", sigma = 0.3, min_sd = 2, seed = 1024)} and
#'   validated once against the source project's precomputed scan. See
#'   \code{data-raw/prepare_znb_data.R}.
#' @seealso \code{\link{generate_spm_data}}, \code{\link{preprocess_spm}}
#' @family datasets
"znb_spm"

#' Observed structural-divergence profile for 1znb_A
#'
#' The per-site observed divergence used as the fitting target for the example
#' protein. This is the template for the \code{observed_data} argument of
#' \code{run_msa_bayesian_analysis}.
#'
#' @format A tibble with 225 rows and 2 columns:
#' \describe{
#'   \item{pdb_site}{Integer PDB residue number (the structure-anchored site key).}
#'   \item{lrmsd_obs}{Observed log structural divergence (log RMSD).}
#' }
#' @details The fit is keyed by \code{pdb_site}: the package maps it internally to
#'   its response-site index. Other per-site descriptors (e.g. distance to the
#'   active site, log RMSF) are intentionally not stored here, as they are
#'   recomputable from \code{znb_wt} (see \code{add_site_properties}; log RMSF is
#'   \code{log(sqrt(get_msf_site(wt)))}).
#' @source Echave & Carpentier (2024) dataset, filtered to pdb_chain == "1znb_A".
#'   See \code{data-raw/prepare_znb_data.R}.
#' @seealso \code{\link{run_msa_bayesian_analysis}}, \code{\link{add_site_properties}}
#' @family datasets
"znb_profile"

#' Active-site metadata for 1znb_A
#'
#' One-row dataset giving the catalytic/active-site residues of the example
#' protein, used by \code{get_active_site}.
#'
#' @format A tibble with 1 row and 4 columns:
#' \describe{
#'   \item{mcsa_id}{M-CSA entry identifier (character).}
#'   \item{pdb_chain}{PDB chain identifier, "1znb_A" (character).}
#'   \item{pdb_site_active}{Comma-separated PDB residue numbers of active sites.}
#'   \item{site_active}{Comma-separated internal site indices of active sites.}
#' }
#' @source Echave & Carpentier (2024) dataset, filtered to pdb_chain == "1znb_A".
#'   See \code{data-raw/prepare_znb_data.R}.
#' @seealso \code{\link{get_active_site}}
#' @family datasets
"znb_dataset"
