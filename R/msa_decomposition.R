#' Split a divergence profile into mutation, stability, and activity contributions
#'
#' Decomposes a structural-divergence profile into the part attributable to
#' mutation alone and the parts added by each selection pressure. The split is
#' sequential along the model progression no-selection -> stability -> activity:
#' `phi_mut` is the mutation-only divergence (the MM variant), `phi_stab` is the
#' extra divergence explained when stability selection is added (MS minus MM), and
#' `phi_act` is the further change when activity selection is added (MSA minus MS).
#' The three contributions sum exactly to the full-model profile.
#'
#' The inputs are four numeric vectors of equal length -- the divergence under each
#' of the four model variants. The function works on the values alone and does not
#' care whether they are indexed by site or by mode, so the same call serves both;
#' the caller supplies whichever four columns it holds and attaches the result.
#'
#' The public entry points that own the four-variant evaluation and attach the
#' result are [calculate_decomposition_i_msa()] / [calculate_decomposition_n_msa()]
#' (forward, at one `(a1, a2)`) and [predict_decomposition_i_msa_agq()] (posterior
#' bands, propagated across quadrature nodes); this kernel is the pure math they
#' share and is internal.
#'
#' @param mm,ms,ma,msa Numeric vectors of equal length giving the divergence under
#'   the four model variants: mutation-only (`mm`), plus-stability (`ms`),
#'   plus-activity (`ma`), and the full model (`msa`). `ma` is accepted but not used
#'   by the current sequential formula; it is kept so an alternative decomposition
#'   that needs it can be added later without changing the signature.
#' @return A named list of three numeric vectors -- `phi_mut`, `phi_stab`,
#'   `phi_act` -- one value per input element.
#' @noRd
calculate_msa_decomposition <- function(mm, ms, ma, msa) {
  # Sequential M0 -> MM -> MS -> MSA decomposition. `ma` is unused here but kept
  # in the signature so the future Shapley `method` is a non-breaking addition.
  list(
    phi_mut  = mm,
    phi_stab = ms - mm,
    phi_act  = msa - ms
  )
}
