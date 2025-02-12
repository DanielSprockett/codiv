
<!-- README.md is generated from README.Rmd. Please edit that file -->

# codiv

Performs Host-Microbe Codiversification Scans

<!-- badges: start -->
<!-- badges: end -->

------------------------------------------------------------------------

The `codiv` R package performs codiversification scans between pairs of
host and symbiont phylogenetic trees. The main functions accepts a host
tree, a symbiont tree, and a data.frame that links which symbionts were
isolated from which hosts. It outputs a correlation coefficient for each
node of the symbiont tree denoted how well it’s topology is correlated
with the host tree. Nodes with a high degree of correlation are
consistent with codiversification between hosts and those symbionts.
Other helper functions are also included.

[![DOI](https://zenodo.org/badge/931169650.svg)](https://doi.org/10.5281/zenodo.14859771)

## Installation

You can install the development version of codiv from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("DanielSprockett/codiv")
```

Citation:

> Daniel D. Sprockett, Brian A. Dillard, Abigail A. Landers, Jon G.
> Sanders, Andrew H. Moeller. (2025) **Recent genetic drift in the
> co-diversified gut bacterial symbionts of laboratory mice.** *Nature
> Communications* Accepted.
