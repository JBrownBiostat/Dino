data("pbmcSmall")
library(SingleCellExperiment)

test_that("Test Dino_SCE", {
    pbmc_SCE <- SingleCellExperiment(assays = list("counts" = pbmcSmall))

    expect_that(Dino_SCE(pbmc_SCE, nCores = 2), is_a("SingleCellExperiment"))
})
