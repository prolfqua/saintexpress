# saintexpress

`saintexpress` provides a pure-R implementation of the SAINTexpress
scoring engines for spectral counts and intensities.

Documentation: <https://prolfqua.github.io/saintexpress/>

This package is intentionally focused on the R implementation. It does
not ship or run native SAINTexpress binaries. Native execution is
handled by the separate
[`saintexpressbin`](https://github.com/prolfqua/saintexpressbin)
package, and prolfqua integration is handled by
[`prolfquasaint`](https://github.com/prolfqua/prolfquasaint).

## Installation

``` r
# install.packages("remotes")
remotes::install_github("prolfqua/saintexpress")
```

## Usage

Prepare a SAINT input list with three data frames named `inter`, `prey`,
and `bait`:

``` r
si <- list(
  inter = inter,
  prey = prey,
  bait = bait
)
```

Validate the input shape:

``` r
saintexpress::validate_saint_input(si)
```

Run the spectral-count engine:

``` r
result <- saintexpress::run_saint(si, mode = "spc")
```

Run the intensity engine:

``` r
result <- saintexpress::run_saint(si, mode = "int")
```

`result` is a data frame in the SAINTexpress `list.txt` output shape.

## Optimizers

The default optimizer is base R:

``` r
saintexpress::run_saint(si, mode = "spc", optimizer = "base")
```

The optional `nloptr` optimizer can be used when the `nloptr` package is
installed:

``` r
saintexpress::run_saint(si, mode = "spc", optimizer = "nloptr")
```

## Vignette

A worked example with simulated AP-MS data is included as a package
vignette:

``` r
vignette("saintexpress", package = "saintexpress")
```

It walks through
[`validate_saint_input()`](https://prolfqua.github.io/saintexpress/reference/validate_saint_input.md)
and
[`run_saint()`](https://prolfqua.github.io/saintexpress/reference/run_saint.md)
for both `spc` and `int` modes on a 6-prey/4-bait synthetic experiment
with known true interactors, and mirrors the structure of the companion
[`saintexpressbin`](https://github.com/prolfqua/saintexpressbin)
vignette so the two engines can be compared side by side.

## Testing Scope

This package uses small synthetic inputs for unit tests. TIP49 reference
fixtures and native-vs-R comparison tests live in
[`prolfquasaint`](https://github.com/prolfqua/prolfquasaint), where both
`saintexpress` and `saintexpressbin` are integrated.

## References

SAINTexpress and the SAINT model are described in the original
publications:

- Teo G, Liu G, Zhang J, Nesvizhskii AI, Gingras AC, Choi H (2014).
  SAINTexpress: improvements and additional features in Significance
  Analysis of INTeractome software. *Journal of Proteomics* 100:37-43.
  <https://doi.org/10.1016/j.jprot.2013.10.023>
- Choi H, Glatter T, Gstaiger M, Nesvizhskii AI (2012). SAINT-MS1:
  protein-protein interaction scoring using label-free intensity data in
  affinity purification-mass spectrometry experiments. *Journal of
  Proteome Research* 11:2619-2624. <https://doi.org/10.1021/pr201185r>
- Choi H, Larsen B, Lin ZY, Breitkreutz A, Mellacheruvu D, Fermin D, Qin
  ZS, Tyers M, Gingras AC, Nesvizhskii AI (2011). SAINT: probabilistic
  scoring of affinity purification-mass spectrometry data. *Nature
  Methods* 8:70-73. <https://doi.org/10.1038/nmeth.1541>
