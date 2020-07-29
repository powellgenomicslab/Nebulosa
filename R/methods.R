#' @title Plot gene-weighted 2D kernel density
#' @author Jose Alquicira-Hernandez
#' @param object Seurat or SingleCellExperiment object
#' @param features Features (e.g. genes) to visualize
#' @param slot Type of data: \code{counts} or\code{data} for Seurat objects and
#'  \code{counts},
#' \code{logcounts}, or \code{normcounts} for SingleCellExperiment objects
#' @param joint Return joint density plot? By default \code{FALSE}
#' @param reduction Name of the reduction to visualize. If not provided, last
#'  computed reduction is visualized
#' @param dims Vector of length 2 specifying the dimensions to be plotted. By
#'  default, the first two dimensions are considered.
#' @param method Kernel density estimation method:
#' \itemize{
#' \item \code{ks}: Computes density using the \code{kda} function from the
#'  \code{ks} package.
#' \item \code{wkde}: Computes density using a modified version of the
#' \code{kde2d} function from the \code{MASS}
#' package to allow weights. Bandwidth selection from the \code{ks} package
#'  is used instead.
#' \item \code{sm}: Computes density using the \code{sm.density} function from
#'  the \code{sm} package
#' }
#' @param adjust Numeric value to adjust to bandwidth. Default: 1. Not available
#'  for \code{ks} method
#' @param size Size of the geom to be plotted (e.g. point size)
#' @param shape Shape of the geom to be plotted
#' @param combine Create a single plot? If \code{FALSE}, a list with ggplot
#'  objects is returned
#' @param pal String specifying the viridis color palette to use.
#' Options:
#' \itemize{
#' \item \code{viridis}
#' \item \code{magma}
#' \item \code{cividis}
#' \item \code{inferno}
#' \item \code{plasma}
#' }
#' @return A scatterplot from a given reduction showing the gene-weighted
#'  density
#' @export
#' @examples
#'
#' data <- Seurat::pbmc_small
#' plot_density(data, "CD3E")
setGeneric("plot_density", function(object, features, slot = NULL,
                                    joint = FALSE, reduction = NULL,
                                    dims = c(1, 2),
                                    method = c("ks", "wkde", "sm"),
                                    adjust = 1, size = 1, shape = 16,
                                    combine = TRUE, pal = "viridis")
    standardGeneric("plot_density"))

#' @importFrom Seurat GetAssayData Reductions Embeddings FetchData
#' @export
#' @describeIn plot_density Plot gene-weighted 2D kernel density
setMethod("plot_density", signature("Seurat"),
          function(object,
                   features, slot = NULL, joint = FALSE, reduction = NULL,
                   dims = c(1, 2), method = c("ks", "wkde", "sm"), adjust = 1,
                   size = 1, shape = 16, combine = TRUE, pal = "viridis") {

              # Validate dimensions -----
              .validate_dimensions(dims)
              # Validate existence of reduction
              if (!is.null(reduction)) {
                  if (!(reduction %in% names(slot(object, "reductions")))) {
                      stop("No reduction named '", reduction,
                           "' found in object")
                  }
              }

              # Set up default reduction -----
              reduction_list <- Reductions(object)
              if (!length(reduction_list)){
                  stop("No reduction has been computed!")
              }
              if (is.null(reduction)) {
                  reduction <- reduction_list[length(reduction_list)]
              }
              cell_embeddings <- as.data.frame(Embeddings(
                  slot(object, "reductions")[[reduction]])
              )

              # Search for dimensions -----
              cell_embeddings <- .search_dimensions(dims, cell_embeddings,
                                                    reduction)
              # Set up default assay -----
              if (is.null(slot)) slot <- "data"

              # Match kde method
              method <- match.arg(method)

              # Extract metadata
              metadata <- slot(object, "meta.data")

              # Determine type of feature and extract
              if (all(features %in% colnames(metadata))) {
                  vars <- metadata[, features, drop = FALSE]
              }else{
                  exp_data <- FetchData(object, vars = features, slot = slot)
                  vars <- .extract_feature_data(exp_data, features)
              }
              .plot_final_density(vars, cell_embeddings, features, joint,
                                  method, adjust, shape, size, pal, combine)
          })



#' @importFrom SingleCellExperiment reducedDims reducedDim colData
#' @importFrom Matrix Matrix t
#' @export
#' @describeIn plot_density Plot gene-weighted 2D kernel density
setMethod("plot_density", signature("SingleCellExperiment"),
          function(object,
                   features, slot = NULL, joint = FALSE, reduction = NULL,
                   dims = c(1, 2), method = c("ks", "wkde", "sm"), adjust = 1,
                   size = 1, shape = 16, combine = TRUE, pal = "viridis") {

              # Validate dimensions -----
              .validate_dimensions(dims)
              # Validate existence of reduction
              if (!is.null(reduction)) {
                  reductions <- names(reducedDims(object))
                  if (!(reduction %in% reductions)) {
                      stop("No reduction named '", reduction,
                           "' found in object")
                  }
              }

              # Set up default reduction -----
              reduction_list <- names(reducedDims(object))
              if (!length(reduction_list)){
                  stop("No reduction has been computed!")
              }

              if (is.null(reduction)) {
                  reduction <- reduction_list[length(reduction_list)]
              }
              cell_embeddings <- reducedDim(object, reduction)

              # Search for dimensions -----
              cell_embeddings <- .search_dimensions(dims, cell_embeddings,
                                                    reduction)
              # Set up default assay -----
              if (is.null(slot)) slot <- "logcounts"

              # Match kde method
              method <- match.arg(method)

              # Extract metadata
              metadata <- as.data.frame(colData(object))

              # Determine type of feature and extract
              if (all(features %in% colnames(metadata))) {
                  vars <- metadata[, features, drop = FALSE]
              }else{
                  if (slot == "data") slot <- "logcounts"
                  assays <- names(object@assays)
                  if (!slot %in% assays) stop(slot, " assay not found")
                  exp_data <- Matrix::t(object@assays@data[[slot]])
                  vars <- .extract_feature_data(exp_data, features)
              }
              .plot_final_density(vars, cell_embeddings, features, joint,
                                  method, adjust, shape, size, pal, combine)
          })
