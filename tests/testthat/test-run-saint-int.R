.make_synthetic_si_int <- function(seed = 1) {
  set.seed(seed)
  preys <- paste0("Prey", 1:6)
  baits <- c("BaitA", "BaitA", "BaitB", "BaitB", "Ctrl1", "Ctrl2")
  ips <- paste0("IP", seq_along(baits))
  cort <- c("T", "T", "T", "T", "C", "C")

  inter_rows <- list()
  for (i in seq_along(baits)) {
    for (p in preys) {
      if (baits[i] == "Ctrl1" || baits[i] == "Ctrl2") {
        x <- rexp(1, 1 / 1e4)
      } else if (baits[i] == "BaitA" && p %in% c("Prey1", "Prey2")) {
        x <- rexp(1, 1 / 1e7)
      } else if (baits[i] == "BaitB" && p %in% c("Prey3", "Prey4")) {
        x <- rexp(1, 1 / 5e6)
      } else {
        x <- rexp(1, 1 / 1e4)
      }
      if (x > 0) {
        inter_rows[[length(inter_rows) + 1]] <-
          data.frame(ipId = ips[i], baitId = baits[i], preyId = p, quant = x, stringsAsFactors = FALSE)
      }
    }
  }
  inter <- do.call(rbind, inter_rows)
  prey <- data.frame(preyId = preys, preyLength = 500L, preyGeneId = preys, stringsAsFactors = FALSE)
  bait <- data.frame(ipId = ips, baitId = baits, CorT = cort, stringsAsFactors = FALSE)
  list(inter = inter, prey = prey, bait = bait)
}

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
