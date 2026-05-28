.make_synthetic_si <- function(seed = 1) {
  set.seed(seed)
  preys <- paste0("Prey", 1:6)
  baits <- c("BaitA", "BaitA", "BaitB", "BaitB", "Ctrl1", "Ctrl2")
  ips <- paste0("IP", seq_along(baits))
  cort <- c("T", "T", "T", "T", "C", "C")

  inter_rows <- list()
  for (i in seq_along(baits)) {
    for (p in preys) {
      if (baits[i] == "Ctrl1" || baits[i] == "Ctrl2") {
        count <- rpois(1, 0.5)
      } else if (baits[i] == "BaitA" && p %in% c("Prey1", "Prey2")) {
        count <- rpois(1, 20)
      } else if (baits[i] == "BaitB" && p %in% c("Prey3", "Prey4")) {
        count <- rpois(1, 15)
      } else {
        count <- rpois(1, 0.5)
      }
      if (count > 0) {
        inter_rows[[length(inter_rows) + 1]] <-
          data.frame(ipId = ips[i], baitId = baits[i], preyId = p,
                     quant = count, stringsAsFactors = FALSE)
      }
    }
  }
  inter <- do.call(rbind, inter_rows)
  prey <- data.frame(preyId = preys, preyLength = 500L, preyGeneId = preys,
                     stringsAsFactors = FALSE)
  bait <- data.frame(ipId = ips, baitId = baits, CorT = cort,
                     stringsAsFactors = FALSE)
  list(inter = inter, prey = prey, bait = bait)
}

test_that("run_saint(mode = 'spc') returns a non-empty data frame", {
  si <- .make_synthetic_si()
  out <- run_saint(si, mode = "spc", optimizer = "base")
  expect_s3_class(out, "data.frame")
  expect_gt(nrow(out), 0)
  expect_true(all(c("Bait", "Prey", "AvgP", "BFDR") %in% names(out)))
  expect_true(all(out$BFDR >= 0 & out$BFDR <= 1))
})
