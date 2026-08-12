# The elastic network model is penm's, and msamodel does not wrap it.
#
# There used to be a setup_enm() wrapper here. It took a file path, and when that was
# refactored to take a bio3d pdb object it became a pure forward to
# penm::set_enm() -- a copy of penm's signature (so penm's defaults could drift
# out of sync unnoticed) plus one class check. The check belongs in penm, which
# owns the function; validating penm's input from here fixes it for msamodel
# users only.
#
# So set_enm() is re-exported instead: `library(msamodel)` is enough to call it,
# it carries penm's own documentation and defaults, and there is nothing here to
# keep in sync.

#' @importFrom penm set_enm
#' @export
penm::set_enm
