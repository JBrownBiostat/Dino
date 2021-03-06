% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Dino.R
\name{SeuratFromDino}
\alias{SeuratFromDino}
\title{Create Seurat object from Dino normalized data}
\usage{
SeuratFromDino(counts, doNorm = TRUE, doLog = TRUE, ...)
}
\arguments{
\item{counts}{A numeric matrix of count data, either raw (eg. UMIs) or
normalized expression.}

\item{doNorm}{A logical indicating whether to normalize the input
\emph{counts} data before exporting results to a \emph{Seurat} object.
By default, it is assumed that the contents of \emph{counts} raw
expression which should be normalized.}

\item{doLog}{A logical indicating whether normalized counts should be log
transformed with a psuedocount of 1 prior to export.}

\item{...}{Further arguments to pass to \emph{Dino}}
}
\value{
\code{SeuratFromDino} returns a Seurat object using Dino normalized
    and log transformed expression (default) for downstream analysis in the
    Seurat pipeline.

    If \emph{returnMeta = T} is passed to \emph{Dino}, then \emph{depth} and
    \emph{slope} results are stored in the \emph{Misc} slot under the
    names \emph{depth} and \emph{slope} respectively.
}
\description{
\code{SeuratFromDino} is a wrapper simplifying the export of
    \code{Dino} normalized counts to a \emph{Seurat} object for secondary
    analysis.
}
\examples{
# raw data
data("pbmcSmall")
str(pbmcSmall)

# run Dino on raw expression matrix, output Seurat object
pbmcSmall_Seurat <- SeuratFromDino(pbmcSmall)
str(pbmcSmall_Seurat)

}
\references{
Brown, J., Ni, Z., Mohanty, C., Bacher, R. and Kendziorski, C.
    (2020). "Normalization by distributional resampling of high throughput
    single-cell RNA-sequencing data." bioRxiv.
    \href{https://doi.org/10.1101/2020.10.28.359901}{https://doi.org/10.1101/2020.10.28.359901}

    Satija, R., Farrell, J.A., Gennert, D., Schier, A.F. and Regev, A.
    (2015). "Spatial reconstruction of single-cell gene expression data."
    Nat. Biotechnol., 33, 495–502.
    \href{https://doi.org/10.1038/nbt.3192}{https://doi.org/10.1038/nbt.3192}
}
\author{
Jared Brown
}
