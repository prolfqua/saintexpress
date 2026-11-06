# Store SAINT scores in a SummarizedExperiment

Reshapes the result of
[`run_saint()`](https://prolfqua.github.io/saintexpress/reference/run_saint.md)
into a `SummarizedExperiment` with preys in rows and baits in columns.
Every numeric column of the result, for example `AvgP`, `SaintScore`,
`FoldChange` and `BFDR`, becomes an assay; prey-bait pairs without a
score are `NA`. `PreyGene` is stored in the `rowData`. Character
columns, such as the `|`-separated replicate values, are not kept.

## Usage

``` r
saint_scores_to_se(scores)
```

## Arguments

- scores:

  A data frame returned by
  [`run_saint()`](https://prolfqua.github.io/saintexpress/reference/run_saint.md).

## Value

A `SummarizedExperiment` with one assay per numeric score column.

## See also

[`se_to_saint_input()`](https://prolfqua.github.io/saintexpress/reference/se_to_saint_input.md)

## Examples

``` r
si <- list(
    inter = data.frame(
        ipId = c("IP1", "IP1", "IP2", "IP2", "IP3", "IP4"),
        baitId = c("BaitA", "BaitA", "BaitA", "BaitA", "Ctrl", "Ctrl"),
        preyId = c("Prey1", "Prey2", "Prey1", "Prey2", "Prey1", "Prey2"),
        quant = c(20, 1, 18, 1, 1, 1)
    ),
    prey = data.frame(preyId = c("Prey1", "Prey2"), preyGeneId = c("Gene1", "Gene2")),
    bait = data.frame(
        ipId = c("IP1", "IP2", "IP3", "IP4"),
        baitId = c("BaitA", "BaitA", "Ctrl", "Ctrl"),
        CorT = c("T", "T", "C", "C")
    )
)
scores_se <- saint_scores_to_se(run_saint(si, mode = "spc"))
SummarizedExperiment::assay(scores_se, "BFDR")
#>       BaitA
#> Prey1     0
#> Prey2     0
```
