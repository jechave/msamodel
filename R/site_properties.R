#' Join structural site properties onto a per-site table
#'
#' Adds structural descriptors for each site to an existing per-site table, joining
#' by the site index. The descriptors are the distance to the nearest active site,
#' the contact number, the mean square fluctuation, and a discrete shell (a banding
#' of the active-site distance). Useful for relating predicted or observed
#' divergence to a site's structural environment.
#'
#' @param site_data A per-site table keyed by the site index `site` (for example the
#'   `$site` tibble from [calculate_profiles()]).
#' @param wt Wild-type protein structure with an elastic network model, as
#'   returned by [setup_enm()].
#' @param pdb_site_active Integer vector of active-site residue numbers (PDB
#'   numbering).
#' @return The input table with four columns added: `dactive` (distance to the
#'   nearest active site), `cn` (contact number), `msf` (mean square fluctuation),
#'   and `shell` (a factor banding `dactive` into levels 0-7).
#' @seealso [setup_enm()] (builds the `wt` input); [calculate_profiles()] (a
#'   typical source of `site_data`).
#' @family setup
#' @examples
#' \dontrun{
#' ex   <- function(f) system.file("extdata", f, package = "msamodel")
#' wt   <- setup_enm(bio3d::read.pdb(ex("1znb_A.pdb")), node = "ca", d_max = 10.5)
#' act  <- readr::read_csv(ex("znb_active_site.csv"))
#' spm  <- generate_spm_data(wt, pdb_site_active = act$pdb_site, seed = 1024)
#' prof <- calculate_profiles(spm, a1 = 1, a2 = 1)$site
#' add_site_properties(prof, wt, pdb_site_active = act$pdb_site)
#' }
#' @export
add_site_properties <- function(site_data, wt, pdb_site_active) {
  # Create tibble with site properties
  site_props <- tibble(
    site = get_site(wt),
    dactive = get_dactive(wt, pdb_site_active),
    cn = get_cn(wt),
    msf = get_msf_site(wt),
    shell = cut(dactive, breaks = c(-.1, 2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 100), drop = FALSE)
  ) %>%
    mutate(shell = factor(shell, levels = levels(shell), labels = seq(from = 0, to = 7, by = 1)))

  # Join with the site_data
  site_data %>%
    left_join(site_props, by = "site")
}
