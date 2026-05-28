#' Validate a SAINT input list
#'
#' Checks that `si` is a list with exactly the elements `inter`, `prey`, `bait`,
#' and that each is a data frame.
#'
#' @param si A named list with elements `inter`, `prey`, `bait`.
#' @return `TRUE` invisibly on success; otherwise stops with a message.
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
