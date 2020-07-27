test_that("plot_density returns a ggplot object", {
    expect_is(plot_density(Seurat::pbmc_small, "CD8A", reduction = "pca"), "ggplot")
})

test_that("plot_density returns error if no feature is provided", {
    expect_error(plot_density(Seurat::pbmc_small))
})

test_that("plot_density returns error if feature is not present", {
    expect_error(plot_density(Seurat::pbmc_small, "test"))
})

test_that("plot_density returns error if non-existent reduction is provided", {
    expect_error(plot_density(Seurat::pbmc_small, "CD8A", reduction = "test_function"))
})

test_that("plot_density returns error if more than 2 dimensions are provided", {
    expect_error(plot_density(Seurat::pbmc_small, dims = 1:3))
})

test_that("plot_density returns error if only one dimension provided", {
    expect_error(plot_density(Seurat::pbmc_small, dims = 1))
})

test_that("plot_density returns error if data slot is not present", {
    expect_error(plot_density(Seurat::pbmc_small, "CD8A", slot = "test_function"))
})

test_that("plot_density returns error if unknown method is provided", {
    expect_error(plot_density(Seurat::pbmc_small, "CD8A", method = "test_method"))
})
