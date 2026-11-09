#' Score AP-MS interactions with SAINTexpress
#'
#' Runs the pure-R SAINTexpress engine: the spectral-count model
#' (`SAINTexpress-spc`) for `mode = "spc"`, the intensity model
#' (`SAINTexpress-int`) for `mode = "int"`. Both agree with native
#' SAINTexpress 3.6.3 up to small numerical differences; the package tests
#' compare them on the TIP49 data in `inst/extdata/tip49`.
#'
#' @param si A named list with the data frames `inter`, `prey` and `bait`, in
#'   any order. See the section "Input tables".
#' @param mode `"spc"` for spectral counts, `"int"` for intensities.
#' @param optimizer Optimizer for the model parameters (`beta1`, `gamma`):
#'   `"base"` (default) uses [stats::optim()] with method `"L-BFGS-B"`;
#'   `"nloptr"` uses NLopt COBYLA, as native SAINTexpress does, and needs the
#'   **nloptr** package. Both give nearly identical scores.
#' @param R Integer >= 1, default 100. The most replicates averaged into
#'   `AvgP`: for a bait with more than `R` replicates, `AvgP` is the mean of
#'   the `R` highest replicate probabilities. Native option `-R`.
#' @param L Integer >= 1, default 100. Spectral counts only: the number of
#'   control IPs, taken in decreasing order of the prey's count, used to
#'   estimate the prey's background (false-interaction) distribution. Native
#'   option `-L`. Checked but not used for `mode = "int"`, as in native
#'   SAINTexpress-int.
#' @param p2p_mapping `NULL` (default) or a list of known prey-prey
#'   interactions, used as a prior that boosts the scores of preys whose
#'   partners are also detected with the bait. The list has one element per
#'   row of `si$prey`, in the same order; element `i` is an integer vector with
#'   the row numbers of the preys linked to prey `i` (`integer(0)` for none).
#'   It plays the role of the optional fourth input file of native
#'   SAINTexpress (for example GO complexes). With `NULL`, `TopoAvgP` and
#'   `TopoMaxP` equal `AvgP` and `MaxP`.
#'
#' @section Input tables:
#' Columns are matched by name when all required names are present; otherwise
#' by position, in the order listed (as in the headerless native input files).
#' A column that already carries a required name must then be at that name's
#' position. Further columns are ignored.
#'
#' `inter`, one row per prey detected in an IP (purification):
#' | Column | Type | Content |
#' |---|---|---|
#' | `ipId` | character | IP identifier, as in `bait$ipId` |
#' | `baitId` | character | Bait of the IP, as in `bait$baitId` |
#' | `preyId` | character | Prey identifier, as in `prey$preyId` |
#' | `quant` | numeric | Spectral count, whole number >= 0 (`"spc"`), or intensity, > 0 (`"int"`) |
#'
#' `prey`, one row per prey:
#' | Column | Type | Content |
#' |---|---|---|
#' | `preyId` | character | Unique prey identifier |
#' | `preyLength` | numeric | Protein length (optional, not used in scoring) |
#' | `preyGeneId` | character | Gene name, reported as `PreyGene` (optional, defaults to `preyId`) |
#'
#' The `prey` table is matched by name when it has both `preyId` and
#' `preyGeneId`. Otherwise the first column is `preyId`; with two columns the
#' second is `preyLength` if numeric, else `preyGeneId`; with three or more,
#' the second is `preyLength` and the third `preyGeneId`.
#'
#' `bait`, one row per IP:
#' | Column | Type | Content |
#' |---|---|---|
#' | `ipId` | character | Unique IP identifier |
#' | `baitId` | character | Bait name; replicate IPs share it |
#' | `CorT` | character | `"T"` for a test IP, `"C"` for a negative control |
#'
#' Every `inter$ipId` must be in `bait`, every `inter$preyId` in `prey`, each
#' prey at most once per IP, and `inter$baitId` must agree with `bait`. At
#' least one test and one control IP are needed.
#'
#' @return A data frame with one row per bait-prey pair, with the columns of
#' the native SAINTexpress `list.txt`:
#' | Column | Content |
#' |---|---|
#' | `Bait`, `Prey`, `PreyGene` | Bait, prey and prey gene name |
#' | `Spec` / `Intensity` | Prey value per replicate IP, joined by vertical bars (int: normalised, `.` = absent) |
#' | `SpecSum` / `IntensitySum` | Sum over the replicates |
#' | `AvgSpec` / `AvgIntensity` | Mean over the replicates |
#' | `NumReplicates` | Number of replicate IPs of the bait |
#' | `ctrlCounts` / `ctrlIntensity` | Prey value per control IP, joined by vertical bars |
#' | `AvgP` | Probability of a true interaction, averaged over the top `R` replicates |
#' | `MaxP` | Highest replicate probability |
#' | `TopoAvgP`, `TopoMaxP` | `AvgP` and `MaxP` with the `p2p_mapping` prior |
#' | `SaintScore` | `max(AvgP, TopoAvgP)` |
#' | `logOddsScore` / `OddsScore` | Lowest replicate log odds of true vs false interaction |
#' | `FoldChange` | Mean over the replicates divided by the control (background) mean |
#' | `BFDR` | Bayesian false discovery rate of the pair, from `AvgP` |
#' | `boosted_by` | Linked preys from `p2p_mapping` that support the pair, joined by vertical bars |
#'
#' Names before the slash are for `mode = "spc"`, after it for `"int"`.
#' @references Teo G, Liu G, Zhang J, Nesvizhskii AI, Gingras AC, Choi H
#'   (2014). SAINTexpress: Improvements and additional features in
#'   Significance Analysis of INTeractome software. Journal of Proteomics,
#'   100, 37-43. \doi{10.1016/j.jprot.2013.10.023}
#' @seealso [validate_saint_input()]
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
#'
#' # Prey1 and Prey2 are known to interact
#' run_saint(si, mode = "spc", p2p_mapping = list(2L, 1L))
#' @export
run_saint <- function(
    si,
    mode = c("spc", "int"),
    optimizer = c("base", "nloptr"),
    R = 100L,
    L = 100L,
    p2p_mapping = NULL
) {
    mode <- match.arg(mode)
    optimizer <- match.arg(optimizer)
    validate_saint_input(si, mode = mode)
    si <- si[c("inter", "prey", "bait")]
    .saint_check_count(R, "R")
    .saint_check_count(L, "L")
    .saint_check_p2p_mapping(p2p_mapping, nrow(si$prey))
    if (optimizer == "nloptr" && !requireNamespace("nloptr", quietly = TRUE)) {
        stop("optimizer = \"nloptr\" needs the nloptr package.")
    }
    ops <- if (mode == "spc") .saint_spc_ops(L) else .saint_int_ops()
    .saint_run(si, ops, optimizer = optimizer, R = R, p2p_mapping = p2p_mapping)
}
