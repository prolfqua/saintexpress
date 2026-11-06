#' Build SAINT input from a SummarizedExperiment
#'
#' Converts a `SummarizedExperiment` with preys in rows and affinity
#' purifications (IPs) in columns into the SAINT input list used by
#' [run_saint()]. For a `QFeatures` object, pass one of its assays, for
#' example `qf[["proteins"]]`.
#'
#' Entries that are `NA` or not positive count as not observed and are left
#' out of `inter`.
#'
#' @param se A `SummarizedExperiment` with row and column names. The row names
#'   become `preyId`, the column names `ipId`.
#' @param assay Name or index of the assay with spectral counts or intensities.
#' @param baitId Name of the `colData` column with the bait of each IP.
#' @param CorT Name of the `colData` column marking test (`"T"`) and control
#'   (`"C"`) IPs.
#' @param preyGeneId Name of the `rowData` column with the gene name of each
#'   prey.
#' @return A named list with the data frames `inter`, `prey` and `bait`.
#' @seealso [saint_scores_to_se()]
#' @examples
#' counts <- matrix(
#'     c(20, 1, 18, 1, 1, 1, 0, 1),
#'     nrow = 2,
#'     dimnames = list(c("Prey1", "Prey2"), c("IP1", "IP2", "IP3", "IP4"))
#' )
#' se <- SummarizedExperiment::SummarizedExperiment(
#'     assays = list(counts = counts),
#'     rowData = data.frame(preyGeneId = c("Gene1", "Gene2")),
#'     colData = data.frame(
#'         baitId = c("BaitA", "BaitA", "Ctrl", "Ctrl"),
#'         CorT = c("T", "T", "C", "C")
#'     )
#' )
#' si <- se_to_saint_input(se)
#' run_saint(si, mode = "spc")
#' @export
se_to_saint_input <- function(se, assay = 1L, baitId = "baitId", CorT = "CorT", preyGeneId = "preyGeneId") {
    col_data <- SummarizedExperiment::colData(se)
    row_data <- SummarizedExperiment::rowData(se)
    missing_cols <- c(setdiff(c(baitId, CorT), names(col_data)), setdiff(preyGeneId, names(row_data)))
    if (length(missing_cols)) {
        stop("se lacks the colData or rowData column(s): ", paste(missing_cols, collapse = ", "))
    }

    bait <- data.frame(
        ipId = colnames(se),
        baitId = as.character(col_data[[baitId]]),
        CorT = as.character(col_data[[CorT]]),
        stringsAsFactors = FALSE
    )
    quant <- as.matrix(SummarizedExperiment::assay(se, assay))
    observed <- which(!is.na(quant) & quant > 0, arr.ind = TRUE)
    observed_prey <- observed[, 1L]
    observed_ip <- observed[, 2L]
    inter <- data.frame(
        ipId = colnames(se)[observed_ip],
        baitId = bait$baitId[observed_ip],
        preyId = rownames(se)[observed_prey],
        quant = quant[observed],
        stringsAsFactors = FALSE
    )
    prey <- data.frame(
        preyId = rownames(se),
        preyGeneId = as.character(row_data[[preyGeneId]]),
        stringsAsFactors = FALSE
    )
    list(inter = inter, prey = prey, bait = bait)
}

#' Store SAINT scores in a SummarizedExperiment
#'
#' Reshapes the result of [run_saint()] into a `SummarizedExperiment` with
#' preys in rows and baits in columns. Every numeric column of the result, for
#' example `AvgP`, `SaintScore`, `FoldChange` and `BFDR`, becomes an assay;
#' prey-bait pairs without a score are `NA`. `PreyGene` is stored in the
#' `rowData`. Character columns, such as the `|`-separated replicate values,
#' are not kept.
#'
#' @param scores A data frame returned by [run_saint()].
#' @return A `SummarizedExperiment` with one assay per numeric score column.
#' @seealso [se_to_saint_input()]
#' @examples
#' si <- list(
#'     inter = data.frame(
#'         ipId = c("IP1", "IP1", "IP2", "IP2", "IP3", "IP4"),
#'         baitId = c("BaitA", "BaitA", "BaitA", "BaitA", "Ctrl", "Ctrl"),
#'         preyId = c("Prey1", "Prey2", "Prey1", "Prey2", "Prey1", "Prey2"),
#'         quant = c(20, 1, 18, 1, 1, 1)
#'     ),
#'     prey = data.frame(preyId = c("Prey1", "Prey2"), preyGeneId = c("Gene1", "Gene2")),
#'     bait = data.frame(
#'         ipId = c("IP1", "IP2", "IP3", "IP4"),
#'         baitId = c("BaitA", "BaitA", "Ctrl", "Ctrl"),
#'         CorT = c("T", "T", "C", "C")
#'     )
#' )
#' scores_se <- saint_scores_to_se(run_saint(si, mode = "spc"))
#' SummarizedExperiment::assay(scores_se, "BFDR")
#' @export
saint_scores_to_se <- function(scores) {
    preys <- sort(unique(scores$Prey))
    baits <- sort(unique(scores$Bait))
    cells <- cbind(match(scores$Prey, preys), match(scores$Bait, baits))
    numeric_cols <- names(scores)[vapply(scores, is.numeric, logical(1))]
    assays <- lapply(numeric_cols, function(col) {
        values <- matrix(NA_real_, length(preys), length(baits), dimnames = list(preys, baits))
        values[cells] <- scores[[col]]
        values
    })
    names(assays) <- numeric_cols
    SummarizedExperiment::SummarizedExperiment(
        assays = assays,
        rowData = data.frame(PreyGene = scores$PreyGene[match(preys, scores$Prey)], row.names = preys),
        colData = data.frame(Bait = baits, row.names = baits)
    )
}
