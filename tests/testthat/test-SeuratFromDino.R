data("pbmcSmall")

test_that("Test SeuratFromDino", {
    expect_that(SeuratFromDino(pbmcSmall, nCores = 2), is_a("Seurat"))
})
