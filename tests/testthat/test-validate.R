test_that("validate_saint_input accepts valid input for both modes", {
    expect_true(validate_saint_input(.make_synthetic_si(), mode = "spc"))
    expect_true(validate_saint_input(.make_synthetic_si_int(), mode = "int"))
    expect_true(validate_saint_input(read_tip49(), mode = "spc"))
})

test_that("list elements may come in any order", {
    si <- .make_synthetic_si()
    shuffled <- si[c("bait", "inter", "prey")]
    expect_true(validate_saint_input(shuffled))
    expect_identical(run_saint(shuffled), run_saint(si))
})

test_that("columns are matched by name, whatever their order", {
    si <- .make_synthetic_si()
    reordered <- si
    reordered$inter <- si$inter[, c("preyId", "ipId", "quant", "baitId")]
    reordered$bait <- si$bait[, c("CorT", "baitId", "ipId")]
    expect_identical(run_saint(reordered), run_saint(si))
})

test_that("unnamed columns are matched by position", {
    si <- .make_synthetic_si()
    positional <- lapply(si, function(df) stats::setNames(df, paste0("V", seq_along(df))))
    expect_identical(run_saint(positional), run_saint(si))
})

test_that("a required column name at the wrong position is an error", {
    si <- .make_synthetic_si()
    si$inter <- stats::setNames(
        si$inter[, c("preyId", "ipId", "baitId", "quant")],
        c("preyId", "ipId", "baitId", "intensity")
    )
    expect_error(validate_saint_input(si), "inter: column 'ipId' is not at position 1")

    si <- .make_synthetic_si()
    si$prey <- si$prey[, c("preyLength", "preyId")]
    expect_error(validate_saint_input(si), "prey: column 'preyId' is not at position 1")

    si <- .make_synthetic_si()
    si$prey <- data.frame(Protein = si$prey$preyId, preyGeneId = "G", preyLength = 500)
    expect_error(validate_saint_input(si), "prey: column 'preyLength' is not at position 2")
})

test_that("factor id columns are scored like character ids", {
    si <- .make_synthetic_si()
    si$prey$preyGeneId[[1]] <- "BaitA"
    as_factors <- si
    as_factors$inter[c("ipId", "baitId", "preyId")] <- lapply(si$inter[c("ipId", "baitId", "preyId")], factor)
    as_factors$bait[] <- lapply(si$bait, factor)
    as_factors$prey$preyGeneId <- factor(si$prey$preyGeneId)
    out <- run_saint(as_factors)
    expect_identical(out, run_saint(si))
    self <- out[out$Bait == "BaitA" & out$Prey == "Prey1", ]
    expect_equal(self$Spec, "0|0")
})

test_that("validate_saint_input rejects malformed lists", {
    expect_error(validate_saint_input(list(inter = data.frame(), prey = data.frame())), "named list")
    expect_error(validate_saint_input(list(a = 1, b = 2, c = 3)), "named list")
    si <- .make_synthetic_si()
    si$bait <- "not a df"
    expect_error(validate_saint_input(si), "si\\$bait must be a data.frame")
    si <- .make_synthetic_si()
    si$prey <- si$prey[0, ]
    expect_error(validate_saint_input(si), "si\\$prey has no rows")
})

test_that("validate_saint_input checks the columns of inter", {
    si <- .make_synthetic_si()
    si$inter <- si$inter[, c("ipId", "baitId", "preyId")]
    expect_error(validate_saint_input(si), "inter: needs columns ipId, baitId, preyId, quant")

    si <- .make_synthetic_si()
    si$inter$quant <- as.character(si$inter$quant)
    expect_error(validate_saint_input(si), "'quant' must be numeric")

    si <- .make_synthetic_si()
    si$inter$quant[[1]] <- NA
    expect_error(validate_saint_input(si), "'quant' must be numeric and finite, without missing values")

    si <- .make_synthetic_si_int()
    si$inter$quant[[1]] <- Inf
    expect_error(validate_saint_input(si, mode = "int"), "'quant' must be numeric and finite")

    si <- .make_synthetic_si()
    si$inter$quant[[1]] <- -1
    expect_error(validate_saint_input(si, mode = "spc"), "spectral counts in 'quant' must be whole numbers >= 0")

    si <- .make_synthetic_si()
    si$inter$quant[[1]] <- 1.9
    expect_error(validate_saint_input(si, mode = "spc"), "spectral counts in 'quant' must be whole numbers >= 0")
    expect_true(validate_saint_input(si, mode = "int"))

    si <- .make_synthetic_si()
    si$inter <- rbind(si$inter, si$inter[1, ])
    expect_error(validate_saint_input(si), "inter: prey Prey1 appears more than once in IP IP1")

    si <- .make_synthetic_si_int()
    si$inter$quant[[1]] <- 0
    expect_error(validate_saint_input(si, mode = "int"), "intensities in 'quant' must be > 0")

    si <- .make_synthetic_si()
    si$inter$preyId[[1]] <- NA
    expect_error(validate_saint_input(si), "inter: 'preyId' must be an id column")
})

test_that("validate_saint_input checks bait and prey tables", {
    si <- .make_synthetic_si()
    si$bait$CorT[[1]] <- "X"
    expect_error(validate_saint_input(si), "'CorT' must be 'T' \\(test\\) or 'C' \\(control\\)")

    si <- .make_synthetic_si()
    si$bait$CorT <- "T"
    expect_error(validate_saint_input(si), "at least one test \\('T'\\) and one control \\('C'\\)")

    si <- .make_synthetic_si()
    si$prey$preyId[[2]] <- si$prey$preyId[[1]]
    expect_error(validate_saint_input(si), "prey: duplicated preyId")

    si <- .make_synthetic_si()
    si$bait$ipId[[2]] <- si$bait$ipId[[1]]
    expect_error(validate_saint_input(si), "bait: duplicated ipId")

    si <- .make_synthetic_si()
    si$prey$preyGeneId[[1]] <- NA
    expect_error(validate_saint_input(si), "prey: 'preyGeneId' must be an id column")
})

test_that("validate_saint_input checks that the tables fit together", {
    si <- .make_synthetic_si()
    si$inter$ipId[[1]] <- "IP99"
    expect_error(validate_saint_input(si), "inter: ipId IP99 is not in the bait table")

    si <- .make_synthetic_si()
    si$inter$preyId[[1]] <- "Prey99"
    expect_error(validate_saint_input(si), "inter: preyId Prey99 is not in the prey table")

    si <- .make_synthetic_si()
    si$inter$baitId[[1]] <- "BaitZ"
    expect_error(validate_saint_input(si), "inter: baitId differs from the bait table")
})
