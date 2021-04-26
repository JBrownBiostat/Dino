data("pbmcSmall")

test_that("Test Dino 'counts' input", {
    # Commented out to save time for Bioconductor submission
    # expect_that(Dino(pbmcSmall, nCores = 2), is_a("dgCMatrix"))

    expect_error(Dino(seq_len(10)),
                 "'counts' not a matrix; dim returns null")

    expect_error(Dino(pbmcSmall - 1),
                 "'counts' contains negative values")

    expect_error(Dino(cbind(pbmcSmall, 0)),
                 "'counts' contains samples with all-zero genes")
})

test_that("Test Dino 'minNZ' input", {
    expect_error(Dino(pbmcSmall, minNZ = seq_len(10)),
                 "'minNZ' doesn't have length 1")

    expect_error(Dino(pbmcSmall, minNZ = 1),
                 "'minNZ' less than 2, choose a larger value")

    expect_error(Dino(pbmcSmall[, seq_len(5)]),
                 "all genes have fewer than 'minNZ' non-zero values")
})

test_that("Test Dino 'depth' input", {
    expect_error(Dino(pbmcSmall[, seq_len(20)], depth = log(colSums(pbmcSmall))),
                 "'depth' has length not equal to the number of samples")
})

