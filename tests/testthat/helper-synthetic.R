# Small synthetic SAINT inputs: 6 preys, BaitA/BaitB with 2 replicate IPs each, 2 control IPs.

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
                    data.frame(ipId = ips[i], baitId = baits[i], preyId = p, quant = count, stringsAsFactors = FALSE)
            }
        }
    }
    inter <- do.call(rbind, inter_rows)
    prey <- data.frame(preyId = preys, preyLength = 500L, preyGeneId = preys, stringsAsFactors = FALSE)
    bait <- data.frame(ipId = ips, baitId = baits, CorT = cort, stringsAsFactors = FALSE)
    list(inter = inter, prey = prey, bait = bait)
}

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
