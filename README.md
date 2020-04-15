
<!-- README.md is generated from README.Rmd. Please edit that file -->
singleCellHaystack
------------------

`singleCellHaystack` is a package for finding surprising needles (=genes) in haystacks (=single cell transcriptome data). Single-cell RNA-seq (scRNA-seq) data is often converted to fewer dimensions using Principal Component Analysis (PCA) and represented in 2-dimentional plots (e.g. t-SNE or UMAP plots). `singleCellHaystack` can be used for finding genes that are expressed in subsets of cells that are non-randomly distributed in these multi-dimensional (first principal components) spaces or 2D (t-SNE, UMAP) representations.

Our manuscript about `singleCellHaystack` is now availabe on [bioRxiv](https://www.biorxiv.org/content/10.1101/557967v3).

Documentation and Demo
----------------------

Our [documentation](https://alexisvdb.github.io/singleCellHaystack/) includes a few example applications showing how to use our package:

-   [Application on toy example](articles/a01_toy_example.html)
-   [Application on multi-dimensional coordinates](articles/a02_example_highD_default.html)
-   [Application of the advanced mode on multi-dimensional coordinates](articles/a03_example_highD_advanced.html)
-   [Application on 2D t-SNE coordinates](articles/a04_example_tsne2D_default.html)
-   [Application of the advanced mode on 2D t-SNE coordinates](articles/a05_example_tsne2D_advanced.html)

System Requirements
-------------------

### Hardware Requirements

`singleCellHaystack` requires only a standard computer with sufficient RAM to support running R or RStudio. Memory requirements depend on the size of the input dataset.

### Software Requirements

This package has been tested on Windows (Windows 10), macOS (Mojave 10.14.1 and Catalina 10.15.1), and Linux (CentOS 6.9 and Ubuntu 19.10).

Installation
------------

<!-- You can install the released version of singleCellHaystack from [CRAN](https://CRAN.R-project.org) with: -->
<!-- ``` r -->
<!-- install.packages("singleCellHaystack") -->
<!-- ``` -->
You can install the `singleCellHaystack` from the GitHub repository as shown below. Typical installation times should be less than 1 minute.

``` r
require(remotes)
remotes::install_github("alexisvdb/singleCellHaystack")
```

`singleCellHaystack` depends on the following packages: splines (3.6.0), ggplot2 (3.2.0), reshape2 (1.4.3), and data.table (1.12.2).
