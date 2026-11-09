test_that("R limits the replicates averaged into AvgP", {
    for (mode in c("spc", "int")) {
        default <- run_tip49(mode)
        top1 <- run_tip49(mode, R = 1L)
        expect_true(any(default$AvgP != default$MaxP), label = mode)
        expect_equal(top1$AvgP, top1$MaxP, label = mode)
    }
})

test_that("L sets the controls used in spc and is ignored in int", {
    expect_false(identical(run_tip49("spc", L = 2L)$AvgP, run_tip49("spc")$AvgP))
    expect_identical(run_tip49("int", L = 2L), run_tip49("int"))
})

test_that("p2p_mapping boosts preys whose linked partners are detected", {
    si <- read_tip49()
    p2p <- tip49_p2p_mapping(si)
    linked_preys <- si$prey[[1]][lengths(p2p) > 0]

    for (mode in c("spc", "int")) {
        plain <- run_tip49(mode)
        expect_equal(plain$TopoAvgP, plain$AvgP, label = mode)
        expect_true(all(plain$boosted_by == ""), label = mode)

        topo <- run_tip49(mode, p2p = TRUE)
        expect_equal(topo$AvgP, plain$AvgP, label = mode)
        expect_true(any(topo$TopoAvgP != topo$AvgP), label = mode)
        expect_equal(topo$SaintScore, pmax(topo$TopoAvgP, topo$AvgP), label = mode)

        boosted <- topo$boosted_by != ""
        expect_true(any(boosted), label = mode)
        expect_true(all(topo$Prey[boosted] %in% linked_preys), label = mode)
        boosting <- unlist(strsplit(topo$boosted_by[boosted], "|", fixed = TRUE))
        expect_true(all(boosting %in% linked_preys), label = mode)
    }
})

test_that("run_saint rejects invalid R, L and p2p_mapping", {
    si <- .make_synthetic_si()
    expect_error(run_saint(si, R = 0L), "R must be a single positive integer")
    expect_error(run_saint(si, R = 1.5), "R must be a single positive integer")
    expect_error(run_saint(si, R = c(1L, 2L)), "R must be a single positive integer")
    expect_error(run_saint(si, L = 0L), "L must be a single positive integer")
    expect_error(run_saint(si, L = NA), "L must be a single positive integer")

    nprey <- nrow(si$prey)
    expect_error(run_saint(si, p2p_mapping = list(2L)), "p2p_mapping must be NULL")
    expect_error(
        run_saint(si, p2p_mapping = c(list(nprey + 1L), rep(list(integer()), nprey - 1L))),
        "p2p_mapping must be NULL"
    )
    expect_error(run_saint(si, p2p_mapping = seq_len(nprey)), "p2p_mapping must be NULL")
})
