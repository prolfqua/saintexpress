#' Run the pure-R SAINTexpress engine
#'
#' Dispatches to the spectral-count or intensity scorer.
#'
#' @param si A named list with elements `inter`, `prey`, `bait` (data frames).
#' @param mode `"spc"` for spectral counts, `"int"` for intensities.
#' @param optimizer `"base"` (default; uses [stats::optim()]) or `"nloptr"`
#'   (uses NLopt COBYLA via the optional **nloptr** package).
#' @param ... Additional arguments forwarded to the internal scorers
#'   `saint_spc_r()` or `saint_int_r()` (e.g. `R`, `L`, `frequency`,
#'   `p2p_mapping`).
#' @return A data frame in the SAINTexpress `list.txt` shape.
#' @examples
#' si <- list(
#'   inter = data.frame(
#'     ipId = c("IP1", "IP1", "IP2", "IP2", "IP3", "IP4"),
#'     baitId = c("BaitA", "BaitA", "BaitA", "BaitA", "Ctrl", "Ctrl"),
#'     preyId = c("Prey1", "Prey2", "Prey1", "Prey2", "Prey1", "Prey2"),
#'     quant = c(20, 1, 18, 1, 1, 1)
#'   ),
#'   prey = data.frame(
#'     preyId = c("Prey1", "Prey2"),
#'     preyLength = c(500, 500),
#'     preyGeneId = c("Gene1", "Gene2")
#'   ),
#'   bait = data.frame(
#'     ipId = c("IP1", "IP2", "IP3", "IP4"),
#'     baitId = c("BaitA", "BaitA", "Ctrl", "Ctrl"),
#'     CorT = c("T", "T", "C", "C")
#'   )
#' )
#' run_saint(si, mode = "spc")
#' @export
run_saint <- function(si,
                      mode = c("spc", "int"),
                      optimizer = c("base", "nloptr"),
                      ...) {
  mode <- match.arg(mode)
  optimizer <- match.arg(optimizer)
  validate_saint_input(si)
  if (mode == "spc") {
    saint_spc_r(si, optimizer = optimizer, ...)
  } else {
    saint_int_r(si, optimizer = optimizer, ...)
  }
}
