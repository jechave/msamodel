#' Add distance and other site-specific properties to results
#'
#' Enhances site-specific data with additional properties such as distance to
#' active site, contact number, mean square fluctuation, and shell classification.
#'
#' @param site_data Site-specific data keyed by the internal site index `i`
#'   (e.g. the output of `calculate_dr2_i_msa()`)
#' @param wt Wild type protein structure with ENM
#' @param pdb_site_active Active site residue numbers (integer vector, PDB
#'   numbering)
#' @return Enhanced tibble with additional site properties:
#'         - dactive: Distance to active site
#'         - cn: Contact number
#'         - msf: Mean square fluctuation
#'         - shell: Shell classification (0-7)
#' @family setup
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
