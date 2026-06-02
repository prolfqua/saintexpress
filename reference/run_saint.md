# Run the pure-R SAINTexpress engine

Dispatches to the spectral-count or intensity scorer.

## Usage

``` r
run_saint(si, mode = c("spc", "int"), optimizer = c("base", "nloptr"), ...)
```

## Arguments

- si:

  A named list with elements `inter`, `prey`, `bait` (data frames).

- mode:

  `"spc"` for spectral counts, `"int"` for intensities.

- optimizer:

  `"base"` (default; uses
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html)) or `"nloptr"`
  (uses NLopt COBYLA via the optional **nloptr** package).

- ...:

  Additional arguments forwarded to the internal scorers `saint_spc_r()`
  or `saint_int_r()` (e.g. `R`, `L`, `frequency`, `p2p_mapping`).

## Value

A data frame in the SAINTexpress `list.txt` shape.

## Examples

``` r
si <- list(
  inter = data.frame(
    ipId = c("IP1", "IP1", "IP2", "IP2", "IP3", "IP4"),
    baitId = c("BaitA", "BaitA", "BaitA", "BaitA", "Ctrl", "Ctrl"),
    preyId = c("Prey1", "Prey2", "Prey1", "Prey2", "Prey1", "Prey2"),
    quant = c(20, 1, 18, 1, 1, 1)
  ),
  prey = data.frame(
    preyId = c("Prey1", "Prey2"),
    preyLength = c(500, 500),
    preyGeneId = c("Gene1", "Gene2")
  ),
  bait = data.frame(
    ipId = c("IP1", "IP2", "IP3", "IP4"),
    baitId = c("BaitA", "BaitA", "Ctrl", "Ctrl"),
    CorT = c("T", "T", "C", "C")
  )
)
run_saint(si, mode = "spc")
#>    Bait  Prey PreyGene  Spec SpecSum AvgSpec NumReplicates ctrlCounts AvgP MaxP
#> 1 BaitA Prey1    Gene1 20|18      38      19             2        1|0    1    1
#> 2 BaitA Prey2    Gene2   1|1       2       1             2        0|1    0    0
#>   TopoAvgP TopoMaxP SaintScore logOddsScore FoldChange BFDR boosted_by
#> 1        1        1          1        24.22         38    0           
#> 2        0        0          0        -0.19          2    0           
```
