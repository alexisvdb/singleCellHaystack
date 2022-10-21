
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

`singleCellHaystack` is a package for predicting differentially active
features (e.g. genes) in single-cell and spatial transcriptomics and
genomics data. While `singleCellHaystack` originally focused on the
prediction of differentially expressed genes (DEGs; see
[here](https://doi.org/10.1038/s41467-020-17900-3)), we have updated the
method and made it more generally applicable ([LINK ON
BIORXIV](https:placeholder)). It can now also be used for finding
differentially accessible genomic regions in scATAC-seq, DEGs along a
trajectory, spatial DEGs, or any other features with non-random levels
of activity inside any input space (1D, 2D, or \>2D). It does so without
relying on clustering of samples into arbitrary clusters.
`singleCellHaystack` uses Kullback-Leibler Divergence to find features
that have patterns of activity in subsets of samples that are
non-randomly positioned inside any input space.

## Citations

- Our manuscript describing the updated, more generally applicable
  version of `singleCellHaystack` will be available on bioRxiv soon.

- Our manuscript describing the original implementation of
  `singleCellHaystack` ([version
  0.3.4](https://github.com/alexisvdb/singleCellHaystack/tree/binary))has
  been published in [Nature
  Communications](https://doi.org/10.1038/s41467-020-17900-3).

If you use `singleCellHaystack` in your research please cite our work
using:

Vandenbon A, Diez D (2020). “A clustering-independent method for finding
differentially expressed genes in single-cell transcriptome data.”
*Nature Communications*, *11*(1), 4318. <doi:10.1038/s41467-020-17900-3>
<https://doi.org/10.1038/s41467-020-17900-3>.

## Documentation and Demo

Our [documentation](https://alexisvdb.github.io/singleCellHaystack/)
includes a few example applications showing how to use our package:

- [Toy
  example](https://alexisvdb.github.io/singleCellHaystack/articles/a01_toy_example.html)
- [Multi-dimensional
  coordinates](https://alexisvdb.github.io/singleCellHaystack/articles/examples/a02_example_highD_default.html)
- [Advanced mode on multi-dimensional
  coordinates](https://alexisvdb.github.io/singleCellHaystack/articles/examples/a03_example_highD_advanced.html)
- [Spatial
  transcriptomics](https://alexisvdb.github.io/singleCellHaystack/articles/examples/a04_example_spatial_transcriptomics.html)
- [MOCA 100k
  cells](https://alexisvdb.github.io/singleCellHaystack/articles/examples/a05_moca_100k.html)
- [2D t-SNE
  coordinates](https://alexisvdb.github.io/singleCellHaystack/articles/examples/a06_example_tsne2D_default.html)
- [Advanced mode on 2D t-SNE
  coordinates](https://alexisvdb.github.io/singleCellHaystack/articles/examples/a07_example_tsne2D_advanced.html)

## Installation

You can install the released version of `singleCellHaystack` from
[CRAN](https://CRAN.R-project.org/package=singleCellHaystack) with:

``` r
install.packages("singleCellHaystack")
```

You can also install `singleCellHaystack` from the GitHub repository as
shown below. Typical installation times should be less than 1 minute.

``` r
require(remotes)
remotes::install_github("alexisvdb/singleCellHaystack")
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
