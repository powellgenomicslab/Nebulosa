# Nebulosa

<img src="man/figure/logo.png" align="right" height="280"/>

[![Build Status](https://travis-ci.org/powellgenomicslab/Nebulosa.svg?branch=master)](https://travis-ci.org/powellgenomicslab/Nebulosa)
![https://www.tidyverse.org/lifecycle/#maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)

## Motivation

Due to the sparsity observed in single-cell data (e.g. RNA-seq, ATAC-seq), the 
visualization of cell features (e.g. gene, peak) is frequently affected and 
unclear, especially when it is overlaid with clustering to annotate cell types. 
Nebulosa is an R package to visualize data from single cells based on kernel 
density estimation. It aims to recover the signal from dropped-out features by 
incorporating the similarity between cells allowing a “convolution” of the cell 
features.

## Installation 

`Nebulosa` is available on `Bioconductor` and can be 
installed as follows:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("Nebulosa")
```

See [Nebulosa](https://bioconductor.org/packages/devel/bioc/html/Nebulosa.html) 
for more details.

You can install the developing version of `Nebulosa` from github via `devtools`:

```R
devtools::install_github("powellgenomicslab/Nebulosa")
```

## Vignettes

Nebulosa can use `Seurat` and `SingleCellExperiment` objects. See the 
corresponding vignette:

- [Seurat](https://bioconductor.org/packages/devel/bioc/vignettes/Nebulosa/inst/doc/nebulosa_seurat.html)
- [OSCA-Bioconductor](https://bioconductor.org/packages/devel/bioc/vignettes/Nebulosa/inst/doc/introduction.html)
