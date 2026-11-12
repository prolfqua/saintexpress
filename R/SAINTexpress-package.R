#' SAINTexpress: Pure R Implementation of SAINTexpress
#'
#' Scores bait-prey interactions in affinity purification-mass spectrometry
#' (AP-MS) experiments with the SAINTexpress model (Teo et al. 2014). For each
#' prey detected in the purification of a bait, SAINT compares the prey's
#' spectral counts or intensities in the bait purifications with those in
#' negative-control purifications and returns the probability that the
#' interaction is genuine.
#'
#' The package reimplements the spectral-count (`SAINTexpress-spc`) and
#' intensity (`SAINTexpress-int`) engines of the native C++ SAINTexpress 3.6.3
#' in R. Its output has the columns of the native `list.txt` file.
#'
#' Main functions:
#' * [run_saint()] scores a SAINT input list.
#' * [validate_saint_input()] checks a SAINT input list before scoring.
#'
#' A subset of the real AP-MS dataset TIP49 (the example shipped with
#' SAINTexpress) and the native SAINTexpress 3.6.3 results for it are in
#' `inst/extdata/tip49`; see the vignette.
#'
#' @references Teo G, Liu G, Zhang J, Nesvizhskii AI, Gingras AC, Choi H
#'   (2014). SAINTexpress: Improvements and additional features in
#'   Significance Analysis of INTeractome software. Journal of Proteomics,
#'   100, 37-43. \doi{10.1016/j.jprot.2013.10.023}
#' @importFrom stats coef complete.cases dnorm lm median optim var
#' @keywords internal
"_PACKAGE"
