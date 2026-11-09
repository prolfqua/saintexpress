test_that("run_saint(mode = 'int') returns a non-empty data frame", {
    si <- .make_synthetic_si_int()
    out <- run_saint(si, mode = "int", optimizer = "base")
    expect_s3_class(out, "data.frame")
    expect_gt(nrow(out), 0)
    expect_true(all(c("Bait", "Prey", "AvgP", "BFDR") %in% names(out)))
})

test_that("constant complete controls use the native median-SD fallback", {
    si <- .make_synthetic_si_int()
    constant_control <- si$inter$preyId == "Prey1" &
        si$inter$baitId %in% c("Ctrl1", "Ctrl2")
    si$inter$quant[constant_control] <- 1e4

    expect_warning(
        out <- run_saint(si, mode = "int", optimizer = "base"),
        regexp = paste0(
            "1 complete control profile has zero variance.*Prey1.*",
            "median control SD"
        )
    )

    expect_true(any(out$Prey == "Prey1"))
    numeric_output <- vapply(out, is.numeric, logical(1))
    expect_true(all(vapply(out[numeric_output], function(x) all(is.finite(x)), logical(1))))
})

test_that("all constant complete controls give an informative error", {
    si <- .make_synthetic_si_int()
    all_controls <- si$inter$baitId %in% c("Ctrl1", "Ctrl2")
    si$inter$quant[all_controls] <- 1e4

    expect_error(
        run_saint(si, mode = "int", optimizer = "base"),
        regexp = paste0(
            "all 6 complete control profiles have zero variance.*",
            "Prey1.*Prey6"
        )
    )
})
