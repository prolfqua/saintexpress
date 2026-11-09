test_that("run_saint(mode = 'spc') returns a non-empty data frame", {
    si <- .make_synthetic_si()
    out <- run_saint(si, mode = "spc", optimizer = "base")
    expect_s3_class(out, "data.frame")
    expect_gt(nrow(out), 0)
    expect_true(all(c("Bait", "Prey", "AvgP", "BFDR") %in% names(out)))
    expect_true(all(out$BFDR >= 0 & out$BFDR <= 1))
})
