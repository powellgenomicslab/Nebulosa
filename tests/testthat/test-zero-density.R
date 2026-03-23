test_that("calculate_density returns zeros for all-zero weights", {
    em <- Seurat::Embeddings(SeuratObject::pbmc_small)[, 1:2]
    w <- rep(0, nrow(em))

    res_ks <- Nebulosa:::calculate_density(w, em, method = "ks")
    res_wkde <- Nebulosa:::calculate_density(w, em, method = "wkde")

    expect_equal(res_ks, rep(0, nrow(em)))
    expect_equal(res_wkde, rep(0, nrow(em)))
})

test_that("plot_density handles an all-zero metadata feature", {
    data <- SeuratObject::pbmc_small
    data$all_zero_meta <- 0

    expect_is(plot_density(data, "all_zero_meta"), "ggplot")
    expect_is(plot_density(data, "all_zero_meta", raster = FALSE), "ggplot")
})
