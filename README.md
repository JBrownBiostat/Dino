# Dino
## Dino is an R package for the normalization of single-cell RNA-seq data using a flexible mixture of negative binomials model of expression

This package was developed by *Jared Brown* in *Christina Kendziorski's* lab at the University of Wisconsin-Madison.

Normalization to remove technical or experimental artifacts is critical in the analysis of single-cell RNA-sequencing experiments, even those for which unique molecular identifiers (UMIs) are available. The majority of methods for normalizing single-cell RNA-sequencing data adjust average expression in sequencing depth, but allow the variance and other properties of the gene-specific expression distribution to be non-constant in depth, which often results in reduced power and increased false discoveries in downstream analyses. This problem is exacerbated by the high proportion of zeros present in most datasets.

To address this, *Dino* constructs a flexible negative-binomial mixture model of gene expression. The data are then normalized by sampling from the posterior distribution of expected expression conditional on observed sequencing depth.

## Installation

*Dino* can be installed from GitHub:

```
devtools::install_github("JBrownBiostat/Dino", build_vignettes = TRUE)
```

*Note:* building the vignette can take a few minutes. If you do not require the vignette, consider running with *build_vignettes = FALSE* to save time.

## Quick start

Following installation, the single funtion, **Dino** can be used to return a normalized matrix of gene expression data:

```
normMat <- Dino(rawMat)
```

For further details on implementation, options, and variations, consult the vignette available by running:

```
vignette("Dino")
```

## Code repository

In addition to GitHub, *Dino* is further freely available for download from Zenodo:
[https://zenodo.org/record/4897558#.YLjjnW5Okko](https://zenodo.org/record/4897558#.YLjjnW5Okko)

## Citation

If you use *Dino* in your analysis, please cite our paper:

Brown, J., Ni, Z., Mohanty, C., Bacher, R., and Kendziorski, C. (2020). Normalization by distributional resampling of high throughput single-cell RNA-sequencing data. bioRxiv. [https://doi.org/10.1101/2020.10.28.359901](https://doi.org/10.1101/2020.10.28.359901).

and consider citing our release on Zenodo:

doi: [https://doi.org/10.5281/zenodo.4897558](https://doi.org/10.5281/zenodo.4897558)

## Contact

Jared Brown: brownj AT biostat DOT wisc DOT edu

Christina Kendziorski: kendzior AT biostat DOT wisc DOT edu
