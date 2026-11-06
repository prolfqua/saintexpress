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
simulate_si <- function(seed = 42, mode = c("spc", "int")) {
  mode <- match.arg(mode)
  set.seed(seed)
  preys  <- paste0("Prey", 1:6)
  baits  <- c("BaitA", "BaitA", "BaitB", "BaitB", "Ctrl1", "Ctrl2")
  ips    <- paste0("IP", seq_along(baits))
  cort   <- c("T", "T", "T", "T", "C", "C")

  # One row per IP x prey; recover the bait for each IP.
  grid <- expand.grid(ipId = ips, preyId = preys, stringsAsFactors = FALSE)
  grid$baitId <- baits[match(grid$ipId, ips)]

  # A prey is enriched ("signal") only for its true bait; everything else is background.
  signal <- (grid$baitId == "BaitA" & grid$preyId %in% c("Prey1", "Prey2")) |
            (grid$baitId == "BaitB" & grid$preyId %in% c("Prey3", "Prey4"))

  # Draw every quantity in one vectorized call.
  if (mode == "spc") {
    lambda <- ifelse(signal, ifelse(grid$baitId == "BaitA", 20, 15), 0.5)
    grid$quant <- stats::rpois(nrow(grid), lambda)
  } else {
    mu <- ifelse(signal, ifelse(grid$baitId == "BaitA", 1e7, 5e6), 1e4)
    grid$quant <- stats::rexp(nrow(grid), 1 / mu)
  }

  inter <- grid[grid$quant > 0, c("ipId", "baitId", "preyId", "quant")]
  rownames(inter) <- NULL
  list(
    inter = inter,
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
#>   ..$ ipId  : chr [1:19] "IP1" "IP2" "IP4" "IP6" ...
#>   ..$ baitId: chr [1:19] "BaitA" "BaitA" "BaitB" "Ctrl2" ...
#>   ..$ preyId: chr [1:19] "Prey1" "Prey1" "Prey1" "Prey1" ...
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
saintexpress::validate_saint_input(si_spc)
```

## Spectral-count scoring

``` r
scores_spc <- saintexpress::run_saint(si_spc, mode = "spc", optimizer = "base")
scores_spc[, c("Bait", "Prey", "AvgP", "BFDR", "SaintScore")]
#>    Bait  Prey AvgP BFDR SaintScore
#> 1 BaitA Prey1 1.00 0.00       1.00
#> 2 BaitB Prey1 0.00 0.16       0.00
#> 3 BaitA Prey2 1.00 0.00       1.00
#> 4 BaitB Prey3 1.00 0.00       1.00
#> 5 BaitA Prey4 0.00 0.16       0.00
#> 6 BaitB Prey4 1.00 0.00       1.00
#> 7 BaitA Prey5 0.00 0.16       0.00
#> 8 BaitB Prey6 0.22 0.00       0.22
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
#> BaitA.3 BaitA Prey2    1    0
#> BaitB.4 BaitB Prey3    1    0
#> BaitB.6 BaitB Prey4    1    0
```

## Volcano plot

The `FoldChange` column summarizes enrichment over the estimated
background. We plot it as an effect-size coordinate (`EFC`, here
`log2(FoldChange)`) against the `BFDR` evidence coordinate.

``` r
plot_volcano <- function(scores, title) {
  volcano_data <- scores
  volcano_data$EFC <- log2(volcano_data$FoldChange)
  volcano_data$neg_log10_bfdr <- -log10(pmax(volcano_data$BFDR, 1e-16))
  volcano_data$is_hit <- volcano_data$BFDR <= 0.05 & volcano_data$EFC >= 1

  plot(
    volcano_data$EFC,
    volcano_data$neg_log10_bfdr,
    pch = 19,
    col = ifelse(volcano_data$is_hit, "#B23A48", "#3D6F99"),
    xlab = "EFC (log2 fold change)",
    ylab = "-log10(BFDR)",
    main = title
  )
  abline(h = -log10(0.05), lty = 2, col = "grey40")
  abline(v = 1, lty = 2, col = "grey40")
  text(
    volcano_data$EFC,
    volcano_data$neg_log10_bfdr,
    labels = paste(volcano_data$Bait, volcano_data$Prey, sep = ":"),
    pos = 3,
    cex = 0.7
  )
}

plot_volcano(scores_spc, "SAINTexpress spectral-count volcano plot")
```

![](saintexpress_files/figure-html/volcano-spc-1.png)

## Intensity scoring

The same simulator, with `mode = "int"`, produces continuous abundances.

``` r
si_int <- simulate_si(mode = "int")
scores_int <- saintexpress::run_saint(si_int, mode = "int", optimizer = "base")
top_per_bait(scores_int)[, c("Bait", "Prey", "AvgP", "BFDR")]
#>          Bait  Prey AvgP BFDR
#> BaitA.1 BaitA Prey1    1    0
#> BaitA.3 BaitA Prey2    1    0
#> BaitB.6 BaitB Prey3    1    0
#> BaitB.8 BaitB Prey4    1    0
```

``` r
plot_volcano(scores_int, "SAINTexpress intensity volcano plot")
```

![](saintexpress_files/figure-html/volcano-int-1.png)

## SummarizedExperiment and QFeatures

Quantitative proteomics data in Bioconductor are often held in a
`SummarizedExperiment`, or in a `QFeatures` object that collects several
of them.
[`se_to_saint_input()`](https://prolfqua.github.io/saintexpress/reference/se_to_saint_input.md)
builds the SAINT input from a `SummarizedExperiment` with preys in rows
and IPs in columns: the `colData` gives the bait and the test/control
flag of each IP, the `rowData` the gene name of each prey.

Here we put the simulated spectral counts into a `QFeatures` object:

``` r
counts <- tapply(si_spc$inter$quant, si_spc$inter[c("preyId", "ipId")], sum)
prey <- si_spc$prey[match(rownames(counts), si_spc$prey$preyId), ]
bait <- si_spc$bait[match(colnames(counts), si_spc$bait$ipId), ]
se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = counts),
  rowData = data.frame(preyGeneId = prey$preyGeneId),
  colData = data.frame(baitId = bait$baitId, CorT = bait$CorT)
)
qf <- QFeatures::QFeatures(list(proteins = se))
qf
#> An instance of class QFeatures (type: bulk) with 1 set:
#> 
#>  [1] proteins: SummarizedExperiment with 6 rows and 6 columns
```

The `proteins` assay converts to SAINT input and scores as before:

``` r
si_qf <- saintexpress::se_to_saint_input(qf[["proteins"]], assay = "counts")
scores_qf <- saintexpress::run_saint(si_qf, mode = "spc")
top_per_bait(scores_qf)[, c("Bait", "Prey", "AvgP", "BFDR")]
#>          Bait  Prey AvgP BFDR
#> BaitA.1 BaitA Prey1    1    0
#> BaitA.2 BaitA Prey2    1    0
#> BaitB.5 BaitB Prey3    1    0
#> BaitB.6 BaitB Prey4    1    0
```

[`saint_scores_to_se()`](https://prolfqua.github.io/saintexpress/reference/saint_scores_to_se.md)
returns the scores as a `SummarizedExperiment` with preys in rows and
baits in columns, one assay per score. Unscored prey-bait pairs are
`NA`.

``` r
scores_se <- saintexpress::saint_scores_to_se(scores_qf)
scores_se
#> class: SummarizedExperiment 
#> dim: 6 2 
#> metadata(0):
#> assays(11): SpecSum AvgSpec ... FoldChange BFDR
#> rownames(6): Prey1 Prey2 ... Prey5 Prey6
#> rowData names(1): PreyGene
#> colnames(2): BaitA BaitB
#> colData names(1): Bait
SummarizedExperiment::assay(scores_se, "BFDR")
#>       BaitA BaitB
#> Prey1  0.00  0.16
#> Prey2  0.00    NA
#> Prey3    NA  0.00
#> Prey4  0.16  0.00
#> Prey5  0.16    NA
#> Prey6    NA  0.00
```

## See also

For running the original native SAINTexpress C++ binary on the same
input, see the companion package
[`saintexpressbin`](https://github.com/prolfqua/saintexpressbin).

## Session information

``` r
sessionInfo()
#> R version 4.6.2 (2026-10-29)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.5 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] tidyr_1.3.2                 sass_0.4.10                
#>  [3] MsCoreUtils_1.24.0          generics_0.1.4             
#>  [5] SparseArray_1.12.3          stringi_1.8.9              
#>  [7] lattice_0.23-1              magrittr_2.0.5             
#>  [9] digest_0.6.39               saintexpress_0.99.1        
#> [11] evaluate_1.0.5              grid_4.6.2                 
#> [13] fastmap_1.2.0               plyr_1.8.9                 
#> [15] jsonlite_2.0.0              Matrix_1.7-6               
#> [17] ProtGenerics_1.44.0         purrr_1.2.2                
#> [19] lazyeval_0.2.3              textshaping_1.0.5          
#> [21] jquerylib_0.1.4             abind_1.4-8                
#> [23] cli_3.6.6                   rlang_1.3.0                
#> [25] XVector_0.52.0              Biobase_2.72.0             
#> [27] cachem_1.1.0                DelayedArray_0.38.2        
#> [29] yaml_2.3.12                 otel_0.2.0                 
#> [31] QFeatures_1.22.0            S4Arrays_1.12.1            
#> [33] tools_4.6.2                 reshape2_1.4.5             
#> [35] dplyr_1.2.1                 SummarizedExperiment_1.42.0
#> [37] BiocGenerics_0.58.1         MultiAssayExperiment_1.38.0
#> [39] vctrs_0.7.3                 R6_2.6.1                   
#> [41] matrixStats_1.5.0           stats4_4.6.2               
#> [43] lifecycle_1.0.5             stringr_1.6.0              
#> [45] Seqinfo_1.2.0               S4Vectors_0.50.3           
#> [47] fs_2.1.0                    htmlwidgets_1.6.4          
#> [49] IRanges_2.46.0              clue_0.3-68                
#> [51] MASS_7.3-66                 ragg_1.5.2                 
#> [53] cluster_2.1.8.3             pkgconfig_2.0.3            
#> [55] desc_1.4.3                  pillar_1.11.1              
#> [57] pkgdown_2.2.1               bslib_0.12.0               
#> [59] Rcpp_1.1.2                  glue_1.8.1                 
#> [61] systemfonts_1.3.2           tidyselect_1.2.1           
#> [63] tibble_3.3.1                xfun_0.62                  
#> [65] GenomicRanges_1.64.0        MatrixGenerics_1.24.0      
#> [67] knitr_1.52                  AnnotationFilter_1.36.0    
#> [69] igraph_2.3.4                htmltools_0.5.9            
#> [71] rmarkdown_2.32              compiler_4.6.2
```
