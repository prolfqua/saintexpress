# Changelog

## SAINTexpress 1.0.0

- First CRAN release.

- The package is renamed from saintexpress to SAINTexpress, the name of
  the program it implements: load it with
  [`library(SAINTexpress)`](https://prolfqua.github.io/SAINTexpress/).

- The GitHub repository and documentation site moved to the new name:
  <https://github.com/prolfqua/SAINTexpress> and
  <https://prolfqua.github.io/SAINTexpress/>.

- SAINTexpress suggests no Bioconductor package any more: Biobase and
  BiocStyle are gone from Suggests.

- The vignette no longer scores an `ExpressionSet`. It builds the three
  SAINT tables from a quantity matrix with a sample table (`col_data`)
  and a feature table (`row_data`), and from one long table with one row
  per protein and sample, and shows that both give the same scores.

## SAINTexpress 0.99.4

- Input columns are matched by name. Before, the first four columns of
  `inter` and the first three of `bait` were renamed by position, so
  correctly named columns in another order gave wrong scores without an
  error. Tables without the required names, such as the headerless
  native input files, are still matched by position; a required name at
  the wrong position is now an error.

- `validate_saint_input(si, mode)` checks columns, types, finite `quant`
  (whole numbers for spectral counts), `CorT` values, duplicated ids,
  each prey at most once per IP, and that the tables fit together. If it
  passes,
  [`run_saint()`](https://prolfqua.github.io/SAINTexpress/reference/run_saint.md)
  accepts the input.
  [`run_saint()`](https://prolfqua.github.io/SAINTexpress/reference/run_saint.md)
  accepts the elements of `si` in any order.

- Id columns given as factors are treated like character ids. Before, a
  factor `baitId` stopped the bait’s own prey from being masked, as
  native SAINTexpress does.

- [`run_saint()`](https://prolfqua.github.io/SAINTexpress/reference/run_saint.md)
  has the documented arguments `R`, `L` and `p2p_mapping` instead of
  `...`, and checks them. `frequency` is removed: it had no effect, as
  in native SAINTexpress 3.6.3. `L` applies to spectral counts only.

- `optimizer = "nloptr"` stops with an error when nloptr is not
  installed, instead of using
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html) without notice.

- A subset of the TIP49 example data of SAINTexpress (4 baits, 9
  controls, 78 preys) and the native SAINTexpress 3.6.3 results for it
  are in `inst/extdata/tip49`; the tests and the vignette use it.

- New package help page `?saintexpress`;
  [`?run_saint`](https://prolfqua.github.io/SAINTexpress/reference/run_saint.md)
  documents the input and output columns.

- The vignette introduces AP-MS interaction scoring, runs on the TIP49
  data, compares the scores with native SAINTexpress, explains the
  tuning arguments and shows a wrapper that scores an `ExpressionSet`
  and returns the scores as one (Biobase in Suggests).

## SAINTexpress 0.99.2

- The vignette shows how to score data held in `SummarizedExperiment`,
  `QFeatures`, `ExpressionSet` or `anndataR`: the quantity matrix,
  sample table and feature table become the SAINT input in a few lines
  of base R, and the scores go back into a prey x bait matrix.
  saintexpress depends on none of these packages.

## SAINTexpress 0.99.0

- Licence changed from MIT to GPL (\>= 3). The package is a derivative
  work of the GPL-3-or-later SAINTexpress C++ sources, so a permissive
  licence was not applicable. Hyungwon Choi and Damian Fermin are now
  credited as authors and copyright holders.

- Attribution added: the README now documents that this package is an
  AI-assisted port of the GPL-3 licensed SAINTexpress C++ sources by
  Hyungwon Choi and Damian Fermin.

- BiocStyle moved from Imports to Suggests; it is only used to build the
  vignette.

- Author credits updated: Alexey Nesvizhskii is credited as an author,
  Witold Wolski is now the maintainer.

- Intensity scoring now handles constant complete control profiles with
  the median-control-SD fallback used by native SAINTexpress, while
  reporting the affected prey identifiers.

- Initial pure-R implementation of SAINTexpress spectral-count and
  intensity scoring with SAINT input validation.
