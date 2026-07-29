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
#' Model-ready single-point-mutation scan of every site of the example protein: the
#' \code{spm} object every divergence-profile and fitting function consumes directly.
#' It carries, for the 2280 mutants (228 sites times 10 replicates, wild-type states
#' dropped), the per-mutant energy changes and both response-axis divergence matrices.
#'
#' @format An \code{spm} object (a list of class \code{"spm"}) with four elements:
#' \describe{
#'   \item{energy_data}{Tibble, one row per mutant, columns \code{j} (mutated site),
#'     \code{m} (replicate), \code{ddg_jm} (stability change), \code{ddgact_jm}
#'     (activity change).}
#'   \item{dr2_ijm}{Numeric \code{[mutant x site]} matrix of per-site squared
#'     displacements; each site value is the squared displacement of a residue
#'     (L2 over its three Cartesian coordinates). 2280 rows times 228 site columns.}
#'   \item{dr2_njm}{Numeric \code{[mutant x mode]} matrix of per-mode squared
#'     displacements; each value is the squared projection of the displacement onto a
#'     normal mode. 2280 rows times 678 mode columns.}
#'   \item{site_map}{Tibble mapping the site index \code{i} to its PDB residue number
#'     \code{pdb_site} (228 rows). Modes are not residue-anchored, so there is no map
#'     for \code{dr2_njm}.}
#' }
#' @details Generated (not copied) by \code{generate_spm_data(znb_wt,
#'   n_mutations = 10, model = "lfenm", sigma = 0.3, min_sd = 2,
#'   pdb_site_active = ..., seed = 1024)} and validated once against the source
#'   project's precomputed scan. See \code{data-raw/prepare_znb_data.R}.
#' @seealso \code{\link{generate_spm_data}}
#' @family datasets
"znb_spm"

#' Observed structural-divergence profile for 1znb_A
#'
#' The per-site observed divergence used as the fitting target for the example
#' protein. This is the template for the \code{observed_data} argument of the
#' fitter \code{\link{fit_lrmsd_i_msa_ml}}.
#'
#' @format A tibble with 225 rows and 2 columns:
#' \describe{
#'   \item{pdb_site}{Integer PDB residue number (the structure-anchored site key).}
#'   \item{lrmsd_i_obs}{Observed log structural divergence (log RMSD).}
#' }
#' @details The fit is keyed by \code{pdb_site}: the package maps it internally to
#'   its response-site index. Other per-site descriptors (e.g. distance to the
#'   active site, log RMSF) are intentionally not stored here, as they are
#'   recomputable from \code{znb_wt} (see \code{add_site_properties}; log RMSF is
#'   \code{log(sqrt(get_msf_site(wt)))}).
#' @source Echave & Carpentier (2024) dataset, filtered to pdb_chain == "1znb_A".
#'   See \code{data-raw/prepare_znb_data.R}.
#' @seealso \code{\link{fit_lrmsd_i_msa_ml}}, \code{\link{add_site_properties}}
#' @family datasets
"znb_profile"

#' SYNTHETIC observed structural-divergence profile (mode form) for 1znb_A
#'
#' Per-mode observed divergence used as the fitting target for the mode arm
#' (\code{fit_lrmsd_n_msa_ml}). This is the template for the \code{observed_data}
#' argument of the mode fit, keyed by the mode index \code{n}.
#'
#' \strong{This data is SYNTHETIC, not empirical.} No empirical mode profile exists
#' yet (computing observed \code{dr2_n} from homologous structures + an alignment is
#' out of scope for \pkg{msamodel} -- it belongs to a future
#' protein-evolution-patterns package). It is a stand-in so the mode fit can be
#' exercised end-to-end until that package exists.
#'
#' @format A tibble with 678 rows (one per normal mode) and 2 columns:
#' \describe{
#'   \item{n}{Integer mode index (the response-mode key; modes are not
#'     structure-anchored, so there is no \code{pdb_site}).}
#'   \item{lrmsd_n_obs}{Synthetic observed log structural divergence per mode.}
#' }
#' @details Generated by fitting the SITE model to the real \code{znb_profile} by
#'   maximum likelihood (\code{fit_lrmsd_i_msa_ml}) for a realistic \code{(a1, a2)},
#'   evaluating the mode forward map (\code{calculate_dr2_n_msa}) at that point, and
#'   adding seeded Gaussian noise (seed 2025, sd 0.30). See
#'   \code{data-raw/prepare_znb_data.R}.
#' @source SYNTHETIC -- not empirical. See Details and
#'   \code{data-raw/prepare_znb_data.R}.
#' @seealso \code{\link{fit_lrmsd_n_msa_ml}}, \code{\link{calculate_dr2_n_msa}}
#' @family datasets
"znb_profile_n"

#' Active-site metadata for 1znb_A
#'
#' One-row dataset giving the catalytic/active-site residues of the example
#' protein. Illustrative only: it shows the source active-site table format, but
#' no package function consumes it. The model takes active sites as a plain
#' integer vector of PDB residue numbers (the \code{pdb_site_active} argument of
#' \code{generate_spm_data} / \code{add_site_properties}); for 1znb_A that vector
#' is \code{c(99, 101, 103, 162, 181, 184, 193, 223)} (the \code{pdb_site_active}
#' column below, parsed).
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
#' @seealso \code{\link{generate_spm_data}}, \code{\link{add_site_properties}}
#' @family datasets
"znb_dataset"
