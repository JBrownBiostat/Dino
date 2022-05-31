# Dino
## Dino is an R package for the normalization of single-cell RNA-seq data using a flexible mixture of negative binomials model of expression

This package was developed by *Jared Brown* in *Christina Kendziorski's* lab at the University of Wisconsin-Madison.

Normalization to remove technical or experimental artifacts is critical in the analysis of single-cell RNA-sequencing experiments, even those for which unique molecular identifiers (UMIs) are available. The majority of methods for normalizing single-cell RNA-sequencing data adjust average expression in sequencing depth, but allow the variance and other properties of the gene-specific expression distribution to be non-constant in depth, which often results in reduced power and increased false discoveries in downstream analyses. This problem is exacerbated by the high proportion of zeros present in most datasets.

To address this, *Dino* constructs a flexible negative-binomial mixture model of gene expression. The data are then normalized by sampling from the posterior distribution of expected expression conditional on observed sequencing depth.

## Installation

`Dino` is now available on `BioConductor` and can be easily installed from that repository by running:

```
# Install Bioconductor if not present, skip otherwise
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Dino package
BiocManager::install("Dino")

# View vignette from R
browseVignettes("Dino")
```

`Dino` is also available from Github, and bug fixes, patches, and updates are available there first. To install `Dino` from Github, run

```
devtools::install_github('JBrownBiostat/Dino')
```

*Note:* building the vignette can take a few minutes. If you do not require the vignette, consider running with `build_vignettes = FALSE` to save time.

## Vignette

In addition to the option to view the package vignette from `R` (see above), a compiled vinette is also available from the `Dino` page on BioConductor: [http://www.bioconductor.org/packages/release/bioc/html/Dino.html](http://www.bioconductor.org/packages/release/bioc/html/Dino.html)

The vignette includes a fuller description of use cases including code examles as well as the underlying methematics of the method.

## Quick start

Following installation, the single funtion, `Dino` can be used to return a normalized matrix of gene expression data:

```
normMat <- Dino(rawMat)
```

For further details on implementation, options, and variations, consult the vignette available by running:

```
vignette("Dino")
```

## Code repository

In addition to [BioConductor](http://www.bioconductor.org/packages/release/bioc/html/Dino.html) and GitHub, *Dino* is further freely available for download from Zenodo:
[https://zenodo.org/record/4897558#.YLjjnW5Okko](https://zenodo.org/record/4897558#.YLjjnW5Okko)

## Citation

If you use `Dino` in your analysis, please cite our paper:

Brown, J., Ni, Z., Mohanty, C., Bacher, R., and Kendziorski, C. (2021). Normalization by distributional resampling of high throughput single-cell RNA-sequencing data. Bioinformatics, 37, 4123-4128. [https://academic.oup.com/bioinformatics/article/37/22/4123/6306403](https://academic.oup.com/bioinformatics/article/37/22/4123/6306403)

## Contact

With questions, comments, or concerns regarding the `Dino` package, please consider opening an issue on Github. You can also contact us directly:

Jared Brown: ![](/vignettes/JBrownEmail.jpg)

Christina Kendziorski: ![](/vignettes/CKendzEmail.jpg)
