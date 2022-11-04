
<!-- README.md is generated from README.Rmd. Please edit that file -->

## singleCellHaystack

<!-- badges: start -->

[![R-CMD-check](https://github.com/alexisvdb/singleCellHaystack/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/alexisvdb/singleCellHaystack/actions/workflows/R-CMD-check.yaml)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/singleCellHaystack)](https://cran.r-project.org/package=singleCellHaystack)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/singleCellHaystack)](https://cran.r-project.org/package=singleCellHaystack)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/grand-total/singleCellHaystack)](https://cran.r-project.org/package=singleCellHaystack)
<!-- badges: end -->

:warning: We are updating `singleCellHaystack`. For the version
described [here](https://doi.org/10.1038/s41467-020-17900-3), please use
branch “binary”. The version on CRAN is also the same binary version.
The master branch on GitHub is now the updated version. :warning:

`singleCellHaystack` is a package for predicting differentially active
features (e.g. genes) in single-cell and spatial transcriptomics and
genomics data. While `singleCellHaystack` originally focused on the
prediction of differentially expressed genes (DEGs; see
[here](https://doi.org/10.1038/s41467-020-17900-3)), we have updated the
method and made it more generally applicable (coming soon on bioRxiv).
It can now also be used for finding differentially accessible genomic
regions in scATAC-seq, DEGs along a trajectory, spatial DEGs, or any
other features with non-random levels of activity inside any input space
(1D, 2D, or \>2D). It does so without relying on clustering of samples
into arbitrary clusters. `singleCellHaystack` uses Kullback-Leibler
Divergence to find features that have patterns of activity in subsets of
samples that are non-randomly positioned inside any input space.

## Citations

-   Our manuscript describing the updated, more generally applicable
    version of `singleCellHaystack` will be available on bioRxiv soon.

-   Our manuscript describing the original implementation of
    `singleCellHaystack` ([version
    0.3.4](https://github.com/alexisvdb/singleCellHaystack/tree/binary))has
    been published in [Nature
    Communications](https://doi.org/10.1038/s41467-020-17900-3).

If you use `singleCellHaystack` in your research please cite our work
using:

Vandenbon A, Diez D (2020). “A clustering-independent method for finding
differentially expressed genes in single-cell transcriptome data.”
*Nature Communications*, *11*(1), 4318. doi: 10.1038/s41467-020-17900-3
(URL: <https://doi.org/10.1038/s41467-020-17900-3>).

## Documentation and Demo

:warning: We are in the process of updating this documentation :warning:

Our [documentation](https://alexisvdb.github.io/singleCellHaystack/)
includes a few example applications showing how to use our package:

-   [Toy
    example](https://alexisvdb.github.io/singleCellHaystack/articles/a01_toy_example.html)
-   [Single-cell
    RNA-seq](https://alexisvdb.github.io/singleCellHaystack/articles/examples/a02_example_scRNAseq.html)
-   [Spatial transcriptomics using
    Visium](https://alexisvdb.github.io/singleCellHaystack/articles/examples/a03_example_spatial_visium.html)
-   [Spatial transcriptomics using Slide-seq
    V2](https://alexisvdb.github.io/singleCellHaystack/articles/examples/a04_example_spatial_slideseqV2.html)
-   [MOCA 100k
    cells](https://alexisvdb.github.io/singleCellHaystack/articles/examples/a05_moca_100k.html)
    Needs update…
-   [Predicting DEGs along a
    trajectory](https://alexisvdb.github.io/singleCellHaystack/articles/examples/a06_trajectory.html)
    Coming soon…
-   [Analysis of gene set
    activities](https://alexisvdb.github.io/singleCellHaystack/articles/examples/a07_gene_sets.html)
    Coming soon…
-   Anything else?

## Installation

You can install the updated version of `singleCellHaystack` from the
GitHub repository as shown below. Typical installation times should be
less than 1 minute.

``` r
require(remotes)
remotes::install_github("alexisvdb/singleCellHaystack")
```

For the binary version of `singleCellHaystack` as described
[here](https://doi.org/10.1038/s41467-020-17900-3), you can use one of
these options:

You can install the released version of `singleCellHaystack` from
[CRAN](https://CRAN.R-project.org/package=singleCellHaystack) with:

``` r
install.packages("singleCellHaystack")
```

Or, install from the binary branch on GitHub:

``` r
require(remotes)
remotes::install_github("alexisvdb/singleCellHaystack@binary")
```

## System Requirements

### Hardware Requirements

`singleCellHaystack` requires only a standard computer with sufficient
RAM to support running R or RStudio. Memory requirements depend on the
size of the input dataset.

### Software Requirements

This package has been tested on Windows (Windows 10), macOS (Mojave
10.14.1 and Catalina 10.15.1), and Linux (CentOS 6.9 and Ubuntu 19.10).

`singleCellHaystack` depends on the following packages: splines (3.6.0),
ggplot2 (3.2.0), reshape2 (1.4.3).
