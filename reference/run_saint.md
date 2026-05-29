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
