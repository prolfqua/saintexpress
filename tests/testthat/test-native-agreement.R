# Compares run_saint() on the TIP49 subset with the list.txt written by the native
# SAINTexpress 3.6.3 binaries (inst/extdata/tip49, see inst/scripts/tip49.md).
# Both outputs are rounded (spc: 2 decimals, int: 3 decimals), so the
# tolerances are at least one rounding step.

expect_native_agreement <- function(mode, optimizer, tolerance) {
    native <- read_tip49_native(mode)
    r_out <- run_tip49(mode, optimizer = optimizer)

    expect_named(r_out, names(native))
    expect_equal(nrow(r_out), nrow(native))
    comparison <- merge(
        native,
        r_out,
        by = c("Bait", "Prey", "PreyGene"),
        suffixes = c("_native", "_r")
    )
    expect_equal(nrow(comparison), nrow(native))

    for (column in names(tolerance)) {
        native_value <- comparison[[paste0(column, "_native")]]
        r_value <- comparison[[paste0(column, "_r")]]
        delta <- abs(native_value - r_value)
        if (column == "FoldChange") {
            delta <- delta / pmax(1, abs(native_value))
        }
        expect_lte(max(delta), tolerance[[column]], label = paste(mode, optimizer, column))
    }
}

spc_tolerance <- c(
    AvgP = 0.011,
    MaxP = 0.011,
    TopoAvgP = 0.011,
    TopoMaxP = 0.011,
    SaintScore = 0.011,
    BFDR = 0.011,
    logOddsScore = 0.011,
    FoldChange = 1e-8
)

# Measured maxima on the subset: probabilities 0.001, BFDR 0.016,
# OddsScore 0.124, relative FoldChange 0.045 (full TIP49: 0.004, 0.02, 0.086,
# 0.067).
int_tolerance <- c(
    AvgP = 0.005,
    MaxP = 0.005,
    TopoAvgP = 0.005,
    TopoMaxP = 0.005,
    SaintScore = 0.005,
    BFDR = 0.025,
    OddsScore = 0.15,
    FoldChange = 0.07
)

test_that("spc scores agree with native SAINTexpress-spc on TIP49", {
    expect_native_agreement("spc", "base", spc_tolerance)
})

test_that("int scores agree with native SAINTexpress-int on TIP49", {
    expect_native_agreement("int", "base", int_tolerance)
})

test_that("nloptr optimizer agrees with native SAINTexpress on TIP49", {
    skip_if_not_installed("nloptr")
    expect_native_agreement("spc", "nloptr", spc_tolerance)
    expect_native_agreement("int", "nloptr", int_tolerance)
})

test_that("TIP49 output keeps the native self-interaction and control formats", {
    spc <- run_tip49("spc")
    self_spc <- spc[spc$Bait == "ACTR8" & spc$Prey == "ACTR8", ]
    expect_equal(self_spc$Spec, "0|0")
    expect_true(all(grepl("|", spc$ctrlCounts, fixed = TRUE)))
    expect_equal(spc$BFDR[[which.max(spc$AvgP)]], 0)

    int <- run_tip49("int")
    self_int <- int[int$Bait == "ACTR8" & int$Prey == "ACTR8", ]
    expect_equal(self_int$Intensity, ".|.")
    expect_true(any(grepl(".", int$ctrlIntensity, fixed = TRUE)))
    expect_equal(int$BFDR[[which.max(int$AvgP)]], 0)
})
