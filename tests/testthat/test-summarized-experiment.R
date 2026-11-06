make_se <- function() {
    counts <- matrix(
        c(20, 1, 0, 18, NA, 1, 1, 1, 2, 0, 1, 1),
        nrow = 3,
        dimnames = list(prey = c("Prey1", "Prey2", "Prey3"), ip = c("IP1", "IP2", "IP3", "IP4"))
    )
    SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts),
        rowData = data.frame(gene = c("Gene1", "Gene2", "Gene3")),
        colData = data.frame(bait = c("BaitA", "BaitA", "Ctrl", "Ctrl"), type = c("T", "T", "C", "C"))
    )
}

test_that("se_to_saint_input builds the three SAINT tables", {
    si <- se_to_saint_input(make_se(), baitId = "bait", CorT = "type", preyGeneId = "gene")
    expect_identical(names(si), c("inter", "prey", "bait"))
    expect_true(validate_saint_input(si))
    expect_identical(si$bait$ipId, c("IP1", "IP2", "IP3", "IP4"))
    expect_identical(si$bait$CorT, c("T", "T", "C", "C"))
    expect_identical(si$prey$preyGeneId, c("Gene1", "Gene2", "Gene3"))
    # the two zeros and the NA entry are not observed
    expect_identical(nrow(si$inter), 9L)
    expect_false(any(si$inter$ipId == "IP1" & si$inter$preyId == "Prey3"))
    expect_false(any(si$inter$ipId == "IP2" & si$inter$preyId == "Prey2"))
    expect_identical(unique(si$inter$baitId[si$inter$ipId == "IP3"]), "Ctrl")
})

test_that("se_to_saint_input names missing columns", {
    expect_error(se_to_saint_input(make_se(), baitId = "bait", CorT = "CorT", preyGeneId = "gene"), "CorT")
    expect_error(se_to_saint_input(make_se(), baitId = "bait", CorT = "type"), "preyGeneId")
})

test_that("scores from the SummarizedExperiment match scores from the tables", {
    si <- se_to_saint_input(make_se(), baitId = "bait", CorT = "type", preyGeneId = "gene")
    si_table <- si
    si_table$inter <- si$inter[order(si$inter$preyId, si$inter$ipId), ]
    expect_equal(run_saint(si, mode = "spc"), run_saint(si_table, mode = "spc"), ignore_attr = TRUE)
})

test_that("saint_scores_to_se reshapes scores to prey x bait assays", {
    si <- se_to_saint_input(make_se(), baitId = "bait", CorT = "type", preyGeneId = "gene")
    scores <- run_saint(si, mode = "spc")
    scores_se <- saint_scores_to_se(scores)
    expect_s4_class(scores_se, "SummarizedExperiment")
    expect_identical(colnames(scores_se), "BaitA")
    expect_setequal(rownames(scores_se), scores$Prey)
    expect_true(all(c("AvgP", "SaintScore", "FoldChange", "BFDR") %in% SummarizedExperiment::assayNames(scores_se)))
    expect_false("Spec" %in% SummarizedExperiment::assayNames(scores_se))
    bfdr <- SummarizedExperiment::assay(scores_se, "BFDR")
    expect_equal(unname(bfdr[cbind(scores$Prey, scores$Bait)]), scores$BFDR)
    expect_identical(
        SummarizedExperiment::rowData(scores_se)$PreyGene,
        scores$PreyGene[match(rownames(scores_se), scores$Prey)]
    )
})
