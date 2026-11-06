# Build SAINT input from a SummarizedExperiment

Converts a `SummarizedExperiment` with preys in rows and affinity
purifications (IPs) in columns into the SAINT input list used by
[`run_saint()`](https://prolfqua.github.io/saintexpress/reference/run_saint.md).
For a `QFeatures` object, pass one of its assays, for example
`qf[["proteins"]]`.

## Usage

``` r
se_to_saint_input(
  se,
  assay = 1L,
  baitId = "baitId",
  CorT = "CorT",
  preyGeneId = "preyGeneId"
)
```

## Arguments

- se:

  A `SummarizedExperiment` with row and column names. The row names
  become `preyId`, the column names `ipId`.

- assay:

  Name or index of the assay with spectral counts or intensities.

- baitId:

  Name of the `colData` column with the bait of each IP.

- CorT:

  Name of the `colData` column marking test (`"T"`) and control (`"C"`)
  IPs.

- preyGeneId:

  Name of the `rowData` column with the gene name of each prey.

## Value

A named list with the data frames `inter`, `prey` and `bait`.

## Details

Entries that are `NA` or not positive count as not observed and are left
out of `inter`.

## See also

[`saint_scores_to_se()`](https://prolfqua.github.io/saintexpress/reference/saint_scores_to_se.md)

## Examples

``` r
counts <- matrix(
    c(20, 1, 18, 1, 1, 1, 0, 1),
    nrow = 2,
    dimnames = list(c("Prey1", "Prey2"), c("IP1", "IP2", "IP3", "IP4"))
)
se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    rowData = data.frame(preyGeneId = c("Gene1", "Gene2")),
    colData = data.frame(
        baitId = c("BaitA", "BaitA", "Ctrl", "Ctrl"),
        CorT = c("T", "T", "C", "C")
    )
)
si <- se_to_saint_input(se)
run_saint(si, mode = "spc")
#>    Bait  Prey PreyGene  Spec SpecSum AvgSpec NumReplicates ctrlCounts AvgP MaxP
#> 1 BaitA Prey1    Gene1 20|18      38      19             2        1|0    1    1
#> 2 BaitA Prey2    Gene2   1|1       2       1             2        1|1    0    0
#>   TopoAvgP TopoMaxP SaintScore logOddsScore FoldChange BFDR boosted_by
#> 1        1        1          1        24.22         38    0           
#> 2        0        0          0        -1.99          1    0           
```
