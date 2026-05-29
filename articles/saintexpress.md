# Scoring simulated SAINT input with the pure-R engine

This vignette walks through scoring a simulated affinity-purification
experiment with
[`saintexpress::run_saint()`](https://prolfqua.github.io/saintexpress/reference/run_saint.md).
No native binary is involved — everything runs in pure R.

## Simulate a small AP-MS experiment

We define two real baits (`BaitA`, `BaitB`) with two replicates each,
plus two controls. Six preys exist; `BaitA` enriches for `Prey1`/`Prey2`
and `BaitB` for `Prey3`/`Prey4`. The remaining preys are background.

``` r
library(saintexpress)

simulate_si <- function(seed = 42, mode = c("spc", "int")) {
  mode <- match.arg(mode)
  set.seed(seed)
  preys  <- paste0("Prey", 1:6)
  baits  <- c("BaitA", "BaitA", "BaitB", "BaitB", "Ctrl1", "Ctrl2")
  ips    <- paste0("IP", seq_along(baits))
  cort   <- c("T", "T", "T", "T", "C", "C")

  draw <- function(bait, prey) {
    if (bait %in% c("Ctrl1", "Ctrl2")) {
      if (mode == "spc") stats::rpois(1, 0.5) else stats::rexp(1, 1 / 1e4)
    } else if (bait == "BaitA" && prey %in% c("Prey1", "Prey2")) {
      if (mode == "spc") stats::rpois(1, 20) else stats::rexp(1, 1 / 1e7)
    } else if (bait == "BaitB" && prey %in% c("Prey3", "Prey4")) {
      if (mode == "spc") stats::rpois(1, 15) else stats::rexp(1, 1 / 5e6)
    } else {
      if (mode == "spc") stats::rpois(1, 0.5) else stats::rexp(1, 1 / 1e4)
    }
  }

  rows <- list()
  for (i in seq_along(baits)) {
    for (p in preys) {
      q <- draw(baits[i], p)
      if (q > 0) {
        rows[[length(rows) + 1]] <- data.frame(
          ipId = ips[i], baitId = baits[i], preyId = p, quant = q,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  list(
    inter = do.call(rbind, rows),
    prey  = data.frame(preyId = preys, preyLength = 500L, preyGeneId = preys,
                       stringsAsFactors = FALSE),
    bait  = data.frame(ipId = ips, baitId = baits, CorT = cort,
                       stringsAsFactors = FALSE)
  )
}

si_spc <- simulate_si(mode = "spc")
str(si_spc, max.level = 2)
#> List of 3
#>  $ inter:'data.frame':   19 obs. of  4 variables:
#>   ..$ ipId  : chr [1:19] "IP1" "IP1" "IP1" "IP1" ...
#>   ..$ baitId: chr [1:19] "BaitA" "BaitA" "BaitA" "BaitA" ...
#>   ..$ preyId: chr [1:19] "Prey1" "Prey2" "Prey4" "Prey6" ...
#>   ..$ quant : int [1:19] 26 17 1 1 22 22 2 2 15 10 ...
#>  $ prey :'data.frame':   6 obs. of  3 variables:
#>   ..$ preyId    : chr [1:6] "Prey1" "Prey2" "Prey3" "Prey4" ...
#>   ..$ preyLength: int [1:6] 500 500 500 500 500 500
#>   ..$ preyGeneId: chr [1:6] "Prey1" "Prey2" "Prey3" "Prey4" ...
#>  $ bait :'data.frame':   6 obs. of  3 variables:
#>   ..$ ipId  : chr [1:6] "IP1" "IP2" "IP3" "IP4" ...
#>   ..$ baitId: chr [1:6] "BaitA" "BaitA" "BaitB" "BaitB" ...
#>   ..$ CorT  : chr [1:6] "T" "T" "T" "T" ...
```

## Validate the input tables

``` r
validate_saint_input(si_spc)
```

## Spectral-count scoring

``` r
scores_spc <- run_saint(si_spc, mode = "spc", optimizer = "base")
scores_spc[, c("Bait", "Prey", "AvgP", "BFDR", "SaintScore")]
#>    Bait  Prey AvgP BFDR SaintScore
#> 1 BaitA Prey1 1.00 0.00       1.00
#> 2 BaitA Prey2 1.00 0.00       1.00
#> 3 BaitA Prey4 0.00 0.26       0.00
#> 4 BaitA Prey6 0.05 0.12       0.05
#> 5 BaitA Prey5 0.38 0.00       0.38
#> 6 BaitB Prey3 1.00 0.00       1.00
#> 7 BaitB Prey4 1.00 0.00       1.00
#> 8 BaitB Prey2 0.00 0.26       0.00
#> 9 BaitB Prey6 0.00 0.26       0.00
```

The true interactors (`Prey1`/`Prey2` for `BaitA`, `Prey3`/`Prey4` for
`BaitB`) should land at the top of each bait by `AvgP`:

``` r
top_per_bait <- function(df) {
  df <- df[order(df$Bait, -df$AvgP), ]
  do.call(rbind, by(df, df$Bait, head, 2))
}
top_per_bait(scores_spc)[, c("Bait", "Prey", "AvgP", "BFDR")]
#>          Bait  Prey AvgP BFDR
#> BaitA.1 BaitA Prey1    1    0
#> BaitA.2 BaitA Prey2    1    0
#> BaitB.6 BaitB Prey3    1    0
#> BaitB.7 BaitB Prey4    1    0
```

## Intensity scoring

The same simulator, with `mode = "int"`, produces continuous abundances.

``` r
si_int <- simulate_si(mode = "int")
scores_int <- run_saint(si_int, mode = "int", optimizer = "base")
top_per_bait(scores_int)[, c("Bait", "Prey", "AvgP", "BFDR")]
#>           Bait  Prey AvgP BFDR
#> BaitA.1  BaitA Prey1    1    0
#> BaitA.2  BaitA Prey2    1    0
#> BaitB.9  BaitB Prey3    1    0
#> BaitB.10 BaitB Prey4    1    0
```

## See also

For running the original native SAINTexpress C++ binary on the same
input, see the companion package
[`saintexpressbin`](https://github.com/prolfqua/saintexpressbin).
