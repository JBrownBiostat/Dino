#' Normalize scRNAseq data
#'
#' @description \code{Dino} removes cell-to-cell variation in observed
#' counts due to the effects of sequencing depth from single-cell mRNA
#' sequencing experiments. \code{Dino} was particularly designed with UMI
#' based protocols in mind, but is applicable to non-UMI based chemistries
#' in the library preparation stage of sequencing.
#'
#' @usage Dino(counts, nCores = 2, prec = 3, minNZ = 10,
#'     nSubGene = 1e4, nSubCell = 1e4, depth = NULL, slope = NULL,
#'     minSlope = 1/2, maxSlope = 2, clusterSlope = TRUE,
#'     returnMeta = FALSE, doRQS = FALSE,
#'     emPar = list(maxIter = 100, tol = 0.1, conPar = 15, maxK = 100), ...)
#'
#' @param counts A numeric matrix object of expression counts - usually in
#'     dgCMatrix format for memory efficiency. Column names denote cells
#'     (samples or droplets) and row names denote genes.
#' @param nCores A non-negative integer scalar denoting the number of cores
#'     which should be used. Setting nCores to 0 uses all cores as determined by
#'     running \code{parallel::detectCores()}
#' @param prec A positive integer denoting the number of decimals to which to
#'     round depth (if estimated internally via \code{depth = NULL}) and
#'     normalized counts for computational efficiency.
#' @param minNZ A positive integer denoting the minimum number of non-zero
#'     counts for a gene to be normalized by the Dino algorithm. It is
#'     recommended to pre-filter the \emph{counts} matrix such that all genes
#'     meet this threshold. Otherwise, genes with fewer than \emph{minNZ}
#'     non-zeros will be scaled by depth for normalization.
#' @param nSubGene A positive integer denoting the number of genes to subset
#'     for calculation of \emph{slope}.
#' @param nSubCell A positive integer denoting the number of samples to subset
#'     for calculation of \emph{slope} and the EM algorithm.
#' @param depth A numeric vector of length equal to the columns of counts.
#'     \emph{depth} denotes a median-centered, log-scale measure of cell-wise
#'     sequencing depth. \code{Dino} defaults to defining depth as the
#'     (within-cell) sum of counts across genes, followed by a log and
#'     median-centering transformation.
#' @param slope A numeric scalar denoting the count-depth relationship on
#'     the log-log scale. Typical values are close to 1 (implying a unit
#'     increase in depth corresponds to a unit increase in expected counts on
#'     the log-log scale), but may be higher, particularly in the case of
#'     non-UMI protocols. \code{Dino} defaults to estimating \emph{slope}
#'     internally.
#' @param minSlope A numeric scalar denoting the minimum slope. Fitted slopes
#'     below this value will return a warning and be set to 1
#' @param maxSlope A numeric scalar denoting the maximum slope. Fitted slopes
#'     above this value will return a warning and be set to 1
#' @param clusterSlope A logical indicating whether cells should be
#'     pre-clustered prior to calculation of slope. Under the default where
#'     cells are pre-clustered, cluster is used as a factor in the regression.
#' @param returnMeta A logical indicating whether metadata (sequencing depth
#'     and slope) should be returned.
#' @param doRQS A logical indicating how normalization resampling is to be done.
#'     By default (F), normalization is done by resampling from the full
#'     posterior distribution. Alternately, restricted quantile sampling (RQS)
#'     can be performed to enforce stronger preservation of expression ranks
#'     in normalized data. Currently RQS is considered experimental.
#' @param emPar A list of parameters to send to the EM algorithm.
#'     \emph{maxIter} denotes the maximum number of model updates. \emph{tol}
#'     denotes the cutoff threshold for reductions in the log likelihood
#'     function. \emph{conPar} denotes the concentration parameter for the
#'     resampling. \emph{conPar = 1} implies full resampling from the fitted
#'     distribution. As \emph{conPar} increases, the normalized expression
#'     converges to the scale-factor normalized values. \emph{maxK} denotes the
#'     maximum number of mixture components in the mixture model.
#' @param ... Additional parameters to pass to \code{Scran::quickCluster}.
#'
#' @return \code{Dino} by default returns a matrix of normalized expression
#'     with identical dimensions as \emph{counts}. If \emph{returnMeta = TRUE},
#'     then \code{Dino} returns a list of normalized expression, sequencing
#'     depth, and slope.
#'
#' @references Brown, J., Ni, Z., Mohanty, C., Bacher, R. and Kendziorski, C.
#'     (2020) "Normalization by distributional resampling of high throughput
#'     single-cell RNA-sequencing data." bioRxiv.
#'     \href{https://doi.org/10.1101/2020.10.28.359901}{https://doi.org/10.1101/2020.10.28.359901}
#'
#' @export
#'
#' @examples
#' # raw data
#' data("pbmcSmall")
#' str(pbmcSmall)
#'
#' # run Dino on raw expression matrix
#' pbmcSmall_Norm <- Dino(pbmcSmall)
#' str(pbmcSmall_Norm)
#'
#' @author Jared Brown
#'
#' @importFrom Matrix Matrix
#' @importFrom Matrix colSums
#' @importFrom Matrix rowSums
#' @importFrom Matrix t
Dino <- function(counts, nCores = 2, prec = 3, minNZ = 10,
            nSubGene = 1e4, nSubCell = 1e4, depth = NULL, slope = NULL,
            minSlope = 1/2, maxSlope = 2, clusterSlope = TRUE,
            returnMeta = FALSE, doRQS = FALSE,
            emPar = list(maxIter = 100, tol = 0.1, conPar = 15, maxK = 100),
            ...) {
    ## Perform argument checks
    checkOut <- check_DinoIn(counts, nCores, prec, minNZ,
        nSubGene, nSubCell, depth, slope, returnMeta)
    counts <- t(checkOut$counts)
    nCores <- checkOut$nCores
    prec <- checkOut$prec
    minNZ <- checkOut$minNZ
    nSubGene <- checkOut$nSubGene
    nSubCell <- checkOut$nSubCell
    depth <- checkOut$depth
    slope <- checkOut$slope
    returnMeta <- checkOut$returnMeta


    ## Set depth
    depth <- estDepth_func(depth, counts, prec)
    depthRep <- calcDepthRep(depth)


    ## Calculate slope
    slope <- estSlope_func(
        slope, counts, depth, nSubGene, nSubCell, clusterSlope, nCores, minSlope,
        maxSlope, ...
    )


    ## Resample counts
    normCounts <- dinoResamp_func(
        counts, nSubCell, depth, minNZ, slope, prec, doRQS, emPar, nCores,
        depthRep
    )

    colnames(normCounts) <- colnames(counts)
    rownames(normCounts) <- rownames(counts)
    normCounts <- t(normCounts)

    message("Done")
    if(returnMeta == FALSE) {
        return(normCounts)
    } else {
        return(list(
            normMat = normCounts,
            depth = depth,
            slope = slope
        ))
    }
}


# dinoResamp_func is a helper function to resample expression
dinoResamp_func <- function(
    counts, nSubCell, depth, minNZ, slope, prec, doRQS, emPar, nCores,
    depthRep
) {
    if(nrow(counts) > nSubCell) {
        subInd <- sample.int(nrow(counts), nSubCell)
        depthRep <- calcDepthRep(depth[subInd])
    } else {
        subInd <- seq_len(nrow(counts))
    }
    dinoGenes <- which(colSums(counts[subInd, ] > 0) >= minNZ)
    libGenes <- which(!(seq_len(ncol(counts[subInd, ])) %in% dinoGenes))
    if(length(libGenes) > 0){
        warning(cat("Some genes have expression non-zero expression below ",
                    "'minNZ' when subsampled\nand will be normalized via
            scale-factor"))
        libCounts <- counts[, libGenes, drop = FALSE]
        countList <- splitGenes(counts[, dinoGenes, drop = FALSE])
    } else {
        countList <- splitGenes(counts[, dinoGenes, drop = FALSE])
    }
    prll <- setPar(nCores, countList)

    message("Normalizing counts by resampling")
    if(length(libGenes) > 0) {
        normCounts_lib <- libCounts / exp(depth)
        colnames(normCounts_lib) <- colnames(counts)[libGenes]
        normCounts_Dino <- parResampCounts(countList, depth, depthRep, slope,
                                           prec, subInd, prll, doRQS, emPar)
        colnames(normCounts_Dino) <- colnames(counts)[dinoGenes]
        normCounts <- cbind(normCounts_Dino, normCounts_lib)
        normCounts <- normCounts[, colnames(counts)]
    } else {
        normCounts <- parResampCounts(countList, depth, depthRep, slope,
                                      prec, subInd, prll, doRQS, emPar)
    }

    return(normCounts)
}


# estSlope_func is a helper function to estimate regression slope
estSlope_func <- function(
    slope, counts, depth, nSubGene, nSubCell, clusterSlope, nCores, minSlope,
    maxSlope, ...
) {
    if(is.null(slope)) {
        message("Calculating regression slope")
        slope <- calcSlope(counts, depth, nSubGene, nSubCell, clusterSlope,
                           nCores, ...)
        if(slope < minSlope) {
            warning(cat("Fitted slope (",
                        round(slope, 3),
                        ") below minSlope (",
                        round(minSlope, 3),
                        ") -- setting slope = 1"))
            slope = 1
        }
        if(slope > maxSlope) {
            warning(cat("Fitted slope (",
                        round(slope, 3),
                        ") above maxSlope (",
                        round(maxSlope, 3),
                        ") -- setting slope = 1"))
            slope = 1
        }
    }

    return(slope)
}


# estDepth_func is a helper function to estimate sequencing depth
estDepth_func <- function(depth, counts, prec) {
    if(is.null(depth)) {
        message("Computing sequencing depth")
        depth <- log(rowSums(counts))
        depth <- depth - median(depth)
        depth <- round(depth, prec)
    }

    return(depth)
}


#' Create Seurat object from Dino normalized data
#'
#' @description \code{SeuratFromDino} is a wrapper simplifying the export of
#'     \code{Dino} normalized counts to a \emph{Seurat} object for secondary
#'     analysis.
#'
#' @usage SeuratFromDino(counts, doNorm = TRUE, doLog = TRUE, ...)
#'
#' @param counts A numeric matrix of count data, either raw (eg. UMIs) or
#'     normalized expression.
#' @param doNorm A logical indicating whether to normalize the input
#'     \emph{counts} data before exporting results to a \emph{Seurat} object.
#'     By default, it is assumed that the contents of \emph{counts} raw
#'     expression which should be normalized.
#' @param doLog A logical indicating whether normalized counts should be log
#'     transformed with a psuedocount of 1 prior to export.
#' @param ... Further arguments to pass to \emph{Dino}
#'
#' @return \code{SeuratFromDino} returns a Seurat object using Dino normalized
#'     and log transformed expression (default) for downstream analysis in the
#'     Seurat pipeline.
#'
#'     If \emph{returnMeta = T} is passed to \emph{Dino}, then \emph{depth} and
#'     \emph{slope} results are stored in the \emph{Misc} slot under the
#'     names \emph{depth} and \emph{slope} respectively.
#'
#' @references Brown, J., Ni, Z., Mohanty, C., Bacher, R. and Kendziorski, C.
#'     (2020). "Normalization by distributional resampling of high throughput
#'     single-cell RNA-sequencing data." bioRxiv.
#'     \href{https://doi.org/10.1101/2020.10.28.359901}{https://doi.org/10.1101/2020.10.28.359901}
#'
#'     Satija, R., Farrell, J.A., Gennert, D., Schier, A.F. and Regev, A.
#'     (2015). "Spatial reconstruction of single-cell gene expression data."
#'     Nat. Biotechnol., 33, 495–502.
#'     \href{https://doi.org/10.1038/nbt.3192}{https://doi.org/10.1038/nbt.3192}
#'
#' @export
#'
#' @examples
#' # raw data
#' data("pbmcSmall")
#' str(pbmcSmall)
#'
#' # run Dino on raw expression matrix, output Seurat object
#' pbmcSmall_Seurat <- SeuratFromDino(pbmcSmall)
#' str(pbmcSmall_Seurat)
#'
#' @author Jared Brown
#'
#' @importFrom Seurat CreateAssayObject
#' @importFrom Seurat CreateSeuratObject
#' @importFrom Seurat Misc
SeuratFromDino <- function(counts, doNorm = TRUE, doLog = TRUE, ...) {
    if(doNorm) {
        normOut <- Dino(counts, ...)
        argList <- list(...)
        if("returnMeta" %in% names(argList)) {
            if(argList$returnMeta) {
                if(doLog) {
                    normAssay <-
                        CreateAssayObject(data = log1p(normOut$normMat))
                } else {
                    normAssay <-
                        CreateAssayObject(data = normOut$normMat)
                }
                normAssay <- CreateSeuratObject(normAssay)
                Seurat::Misc(normAssay, "depth") <- normOut$depth
                Seurat::Misc(normAssay, "slope") <- normOut$slope
            } else {
                if(doLog) {
                    normAssay <- CreateAssayObject(data = log1p(normOut))
                } else {
                    normAssay <- CreateAssayObject(data = normOut)
                }
                normAssay <- CreateSeuratObject(normAssay)
            }
        } else {
            if(doLog) {
                normAssay <- CreateAssayObject(data = log1p(normOut))
            } else {
                normAssay <- CreateAssayObject(data = normOut)
            }
            normAssay <- CreateSeuratObject(normAssay)
        }
    } else {
        if(doLog) {
            normAssay <- CreateAssayObject(data = log1p(counts))
        } else {
            normAssay <- CreateAssayObject(data = counts)
        }
        normAssay <- CreateSeuratObject(normAssay)
    }
    return(normAssay)
}


#' Run Dino normalization on a SingleCellExperiment dataset
#'
#' @description \code{Dino_SCE} is a wrapper simplifying the application of the
#'     \emph{Dino} method to data formatted as a \emph{SingleCellExperiment}
#'
#' @usage Dino_SCE(SCE, ...)
#'
#' @param SCE A \emph{SingleCellExperiment} object with unnormalized count data
#'     (eg. raw  UMIs) in the \emph{assays} slot under the name \emph{counts}.
#' @param ... Further arguments to pass to \emph{Dino}
#'
#' @return \code{Dino_SCE} returns a \emph{SingleCellExperiment} object using
#'     Dino normalized expression in the \emph{assays} slot under the
#'     \emph{normcounts} name for downstream analysis.
#'
#'     If \emph{returnMeta = T} is passed to \emph{Dino}, then \emph{depth} and
#'     \emph{slope} results are stored in the \emph{metadata} slot under the
#'     names \emph{depth} and \emph{slope} respectively.
#'
#' @references Brown, J., Ni, Z., Mohanty, C., Bacher, R. and Kendziorski, C.
#'     (2020). "Normalization by distributional resampling of high throughput
#'     single-cell RNA-sequencing data." bioRxiv.
#'     \href{https://doi.org/10.1101/2020.10.28.359901}{https://doi.org/10.1101/2020.10.28.359901}
#'
#'     Amezquita, R.A., Lun, A.T.L., Becht, E., Carey, V.J., Carpp, L.N.,
#'     Geistlinger, L., Marini, F., Rue-Albrecht, K., Risso, D., Soneson, C.,
#'     et al. (2020). "Orchestrating single-cell analysis with Bioconductor."
#'     Nat. Methods, 17, 137–145.
#'     \href{https://doi.org/10.1038/s41592-019-0654-x}{https://doi.org/10.1038/s41592-019-0654-x}
#'
#' @export
#'
#' @examples
#' # raw data
#' data("pbmcSmall")
#' str(pbmcSmall)
#'
#' # format as SingleCellExperiment
#' library(SingleCellExperiment)
#' pbmc_SCE <- SingleCellExperiment(assays = list("counts" = pbmcSmall))
#'
#' # Run Dino
#' pbmc_SCE <- Dino_SCE(pbmc_SCE)
#' str(pbmc_SCE)
#' str(normcounts(pbmc_SCE))
#'
#' @author Jared Brown
#'
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment assayNames
#' @importFrom methods is
#' @importFrom S4Vectors metadata
Dino_SCE <- function(SCE, ...) {
    if(!is(SCE, "SingleCellExperiment")) {
        stop("'SCE' is not a SingleCellExperiment")
    }
    if(is.null(assayNames(SCE))) {
        stop("no 'counts' assay present")
    }
    if(!("counts" %in% assayNames(SCE))) {
        stop("no 'counts' assay present")
    }

    argList <- list(...)
    if("returnMeta" %in% names(argList)) {
        if(argList$returnMeta) {
            countMat <- counts(SCE)
            normOut <- Dino(countMat, ...)
            SingleCellExperiment::normcounts(SCE) <- normOut$normMat
            S4Vectors::metadata(SCE)[["depth"]] <- normOut$depth
            S4Vectors::metadata(SCE)[["slope"]] <- normOut$slope
        } else {
            countMat <- counts(SCE)
            SingleCellExperiment::normcounts(SCE) <- Dino(countMat, ...)
        }
    } else {
        countMat <- counts(SCE)
        SingleCellExperiment::normcounts(SCE) <- Dino(countMat, ...)
    }

    return(SCE)
}
