# saintexpress 0.99.1

* `se_to_saint_input()` builds the SAINT input from a `SummarizedExperiment` with preys in rows and IPs in columns, for example an assay of a `QFeatures` object. `saint_scores_to_se()` returns the `run_saint()` scores as a `SummarizedExperiment` with preys in rows, baits in columns and one assay per score. The vignette shows both with `QFeatures`.

# saintexpress 0.99.0

* Licence changed from MIT to GPL (>= 3). The package is a derivative work of the GPL-3-or-later SAINTexpress C++ sources, so a permissive licence was not applicable. Hyungwon Choi and Damian Fermin are now credited as authors and copyright holders.

* Attribution added: the README now documents that this package is an AI-assisted port of the GPL-3 licensed SAINTexpress C++ sources by Hyungwon Choi and Damian Fermin.

* BiocStyle moved from Imports to Suggests; it is only used to build the vignette.

* Author credits updated: Alexey Nesvizhskii is credited as an author, Witold Wolski is now the maintainer.

* Intensity scoring now handles constant complete control profiles with the median-control-SD fallback used by native SAINTexpress, while reporting the affected prey identifiers.

* Initial pure-R implementation of SAINTexpress spectral-count and intensity scoring with SAINT input validation.
