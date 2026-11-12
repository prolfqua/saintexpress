# saintexpress

`saintexpress` provides a pure-R implementation of the SAINTexpress scoring
engines for spectral counts and intensities.

Documentation: <https://prolfqua.github.io/saintexpress/>

This package is intentionally focused on the R implementation. It does not ship
or run native SAINTexpress binaries. Native execution is handled by the separate
[`saintexpressbin`](https://github.com/prolfqua/saintexpressbin) package, and
prolfqua integration is handled by
[`prolfquasaint`](https://github.com/prolfqua/prolfquasaint).

## Installation

From CRAN:

```r
install.packages("saintexpress")
```

Development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("prolfqua/saintexpress")
```

## Usage

Prepare a SAINT input list with three data frames named `inter`, `prey`, and
`bait`:

```r
si <- list(
  inter = inter,
  prey = prey,
  bait = bait
)
```

Validate the input shape:

```r
saintexpress::validate_saint_input(si, mode = "spc")
```

Run the spectral-count engine:

```r
result <- saintexpress::run_saint(si, mode = "spc")
```

Run the intensity engine:

```r
result <- saintexpress::run_saint(si, mode = "int")
```

`result` is a data frame in the SAINTexpress `list.txt` output shape.

## Optimizers

The default optimizer is base R:

```r
saintexpress::run_saint(si, mode = "spc", optimizer = "base")
```

The optional `nloptr` optimizer can be used when the `nloptr` package is
installed:

```r
saintexpress::run_saint(si, mode = "spc", optimizer = "nloptr")
```

## Vignette

```r
vignette("saintexpress", package = "saintexpress")
```

The vignette scores a subset of the TIP49 example data of SAINTexpress with both models, compares the scores with native SAINTexpress, explains the optimizer and the tuning arguments `R`, `L` and `p2p_mapping`, and shows how to build the three SAINT tables from a quantity matrix with sample and feature tables, or from one long table.

## Testing Scope

The tests run on the TIP49 subset in `inst/extdata/tip49`, together with the `list.txt` files the native SAINTexpress 3.6.3 programs write for it, and on small synthetic inputs. They check that the R scores agree with the native output for both models and both optimizers.

## Origin, Attribution and AI Assistance

`saintexpress` is not an independent invention of the SAINT model. It is a
port of the SAINTexpress reference implementation to pure R.

- The original SAINTexpress C++ sources (`SAINT-MRF-spc` and `SAINT-MRF-int`)
  are Copyright (C) 2011 Hyungwon Choi and Damian Fermin, and are distributed
  under the GNU General Public License, version 3 or later. The upstream
  distribution used as the reference is `SAINTexpress_v3.6.3__2018-03-09.tar.gz`
  from <https://sourceforge.net/projects/saint-apms/files/>.
- The scoring logic in this package was derived from those C++ sources so that
  the numerical output matches the native `list.txt` results. It is a
  derivative work of the upstream implementation, not a clean-room
  reimplementation from the publications alone.
- The translation from C++ to R was carried out with the help of an AI coding
  assistant. Every ported routine was subsequently reviewed by the maintainer
  and checked against the output of the native SAINTexpress binaries; the
  comparison tests live in `tests/`.

Because this package is a derivative work of GPL-3 licensed sources,
`saintexpress` is released under the same terms: GNU General Public License,
version 3 or later. It cannot be redistributed under a more permissive licence.

Credit for the SAINT model and its reference implementation belongs to the
original authors listed in the References section below.

## References

SAINTexpress and the SAINT model are described in the original publications:

- Teo G, Liu G, Zhang J, Nesvizhskii AI, Gingras AC, Choi H (2014).
  SAINTexpress: improvements and additional features in Significance Analysis
  of INTeractome software. *Journal of Proteomics* 100:37-43.
  <https://doi.org/10.1016/j.jprot.2013.10.023>
- Choi H, Glatter T, Gstaiger M, Nesvizhskii AI (2012). SAINT-MS1:
  protein-protein interaction scoring using label-free intensity data in
  affinity purification-mass spectrometry experiments. *Journal of Proteome
  Research* 11:2619-2624. <https://doi.org/10.1021/pr201185r>
- Choi H, Larsen B, Lin ZY, Breitkreutz A, Mellacheruvu D, Fermin D, Qin ZS,
  Tyers M, Gingras AC, Nesvizhskii AI (2011). SAINT: probabilistic scoring of
  affinity purification-mass spectrometry data. *Nature Methods* 8:70-73.
  <https://doi.org/10.1038/nmeth.1541>
