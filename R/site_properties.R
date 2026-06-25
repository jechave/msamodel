#' Join structural site properties onto a per-site table
#'
#' Adds structural descriptors for each site to an existing per-site table, joining
#' by the site index. The descriptors are the distance to the nearest active site,
#' the contact number, the mean square fluctuation, and a discrete shell (a banding
#' of the active-site distance). Useful for relating predicted or observed
#' divergence to a site's structural environment.
#'
#' @param site_data A per-site table keyed by the site index `i` (for example the
#'   output of [calculate_dr2_i_msa()]).
#' @param wt Wild-type protein structure with an elastic network model, as
#'   returned by [setup_enm()].
#' @param pdb_site_active Integer vector of active-site residue numbers (PDB
#'   numbering).
#' @return The input table with four columns added: `dactive` (distance to the
#'   nearest active site), `cn` (contact number), `msf` (mean square fluctuation),
#'   and `shell` (a factor banding `dactive` into levels 0-7).
#' @seealso [setup_enm()] (builds the `wt` input); [calculate_dr2_i_msa()] (a
#'   typical source of `site_data`).
#' @family setup
#' @examples
#' \dontrun{
#' pp <- preprocess_spm(znb_spm)
#' prof <- calculate_dr2_i_msa(pp, a1 = 1, a2 = 1)
#' add_site_properties(prof, znb_wt, pdb_site_active = c(99, 101, 103, 162))
#' }
#' @export
add_site_properties <- function(site_data, wt, pdb_site_active) {
  # Create tibble with site properties
  site_props <- tibble(
    i = get_site(wt),
    dactive = get_dactive(wt, pdb_site_active),
    cn = get_cn(wt),
    msf = get_msf_site(wt),
    shell = cut(dactive, breaks = c(-.1, 2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 100), drop = FALSE)
  ) %>%
    mutate(shell = factor(shell, levels = levels(shell), labels = seq(from = 0, to = 7, by = 1)))

  # Join with the site_data
  site_data %>%
    left_join(site_props, by = "i")
}
