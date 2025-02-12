
<!-- README.md is generated from README.Rmd. Please edit that file -->

# codiv

Performs Host-Microbe Codiversification Scans

<!-- badges: start -->
<!-- badges: end -->

------------------------------------------------------------------------

Package Description: The `codiv` R package performs codiversification
scans between pairs of host and symbiont phylogenetic trees. The main
functions accepts a host tree, a symbiont tree, and a data.frame that
links which symbionts were isolated from which hosts. It outputs a
correlation coefficient for each node of the symbiont tree denoted how
well it’s topology is correlated with the host tree. Nodes with a high
degree of correlation are consistent with codiversification between
hosts and those symbionts. Other helper functions are also included.

## Installation

You can install the development version of codiv from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("DanielSprockett/codiv")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(codiv)
## basic example code
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

This project borrows heavily from previous work. For the time being,
please cite:

> Sanders, Jon G., Daniel D. Sprockett, Yingying Li, Deus Mjungu,
> Elizabeth V. Lonsdorf, Jean-Bosco Ndjango, Alexander V. Georgiev, John
> H. Hart, Crickette Sanz, David Morgan, Martine Peeters, Beatrice H.
> Hahn, Andrew H. Moeller.(2023)**Widespread extinctions of
> co-diversified gut bacterial symbionts from humans** *Nature
> Microbiology*: Accepted in principle.
