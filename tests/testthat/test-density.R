test_that("calculate_density returns a vector of size n", {
    em <- Seurat::Embeddings(SeuratObject::pbmc_small)[, 1:2]
    w <- Seurat::GetAssayData(SeuratObject::pbmc_small)["CD8A", ]
    res <- Nebulosa:::calculate_density(w, em, method = "ks")

    testthat::expect_vector(res, size = nrow(em))
})



test_that("calculate-density returns a list with raw estimates", {
    em <- Seurat::Embeddings(SeuratObject::pbmc_small)[, 1:2]
    w <- Seurat::GetAssayData(SeuratObject::pbmc_small)["CD8A", ]
    res <- Nebulosa:::calculate_density(w, em, method = "ks", map = FALSE)

    testthat::expect_type(res, "list")
})
