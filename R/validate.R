#' Validate a SAINT input list
#'
#' Checks that `si` is a list with exactly the elements `inter`, `prey`, `bait`,
#' and that each is a data frame.
#'
#' @param si A named list with elements `inter`, `prey`, `bait`.
#' @return `TRUE` invisibly on success; otherwise stops with a message.
#' @examples
#' si <- list(
#'   inter = data.frame(ipId = "IP1", baitId = "BaitA", preyId = "Prey1"),
#'   prey = data.frame(preyId = "Prey1"),
#'   bait = data.frame(ipId = "IP1", baitId = "BaitA", CorT = "T")
#' )
#' validate_saint_input(si)
#' @export
validate_saint_input <- function(si) {
  if (!is.list(si) || is.null(names(si)) ||
      !identical(sort(names(si)), c("bait", "inter", "prey"))) {
    stop("si must be a named list with elements 'inter', 'prey', 'bait'.")
  }
  for (nm in c("inter", "prey", "bait")) {
    if (!is.data.frame(si[[nm]])) {
      stop("si$", nm, " must be a data.frame.")
    }
  }
  invisible(TRUE)
}
