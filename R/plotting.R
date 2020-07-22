#' @title Plot gene-weighted 2D kernel density
#' @author Jose Alquicira-Hernandez
#' @param object Seurat or SingleCellExperiment object
#' @param features Features (e.g. genes) to visualize
#' @param thr Numeric threshold to weight the cell density. All cells above this threshold are considered to have the
#' same weight. By default, the gene expression values are used as weights.
#' @param slot Type of data: \code{counts} or\code{data} for Seurat objects and \code{counts},
#' \code{logcounts}, or \code{normcounts} for SingleCellExperiment objects
#' @param reduction Name of the reduction to visualize. If not provided, last computed reduction is visualized
#' @param dims Vector of length 2 specifying the dimensions to be plotted. By default, the first two dimensions are considered.
#' @param method Kernel density estimation method:
#' \itemize{
#' \item \code{ks}: Computes density using the \code{kda} function from the \code{ks} package.
#' Ideal for small-medium datasets (up to 70,000 cells in a desktop computer). Use other methods for larger datasets.
#' \item \code{wkde}: Computes density using a modified version of the \code{kde2d} function from the \code{MASS}
#' package to allow weights. Bandwidth selection from the \code{ks} package is used instead.
#' \item \code{sm}: Computes density using the \code{sm.density} function from the \code{sm} package
#' }
#' @param adjust Numeric value to adjust to bandwidth. Default: 1. Not available for \code{ks} method
#' @param size Size of the geom to be plotted (e.g. point size)
#' @param shape Shape of the geom to be plotted
#' @param combine Create a single plot? If \code{FALSE}, a list with ggplot objects is returned
#' @return A scatterplot from a given reduction showing the gene-weighted density
#' @importFrom Matrix Matrix t
#' @importFrom SingleCellExperiment reducedDims reducedDim colData
#' @importFrom Seurat GetAssayData Reductions Embeddings FetchData
#' @importFrom patchwork wrap_plots
#' @export

plot_density <- function(object, features, thr = NULL, slot = NULL, reduction = NULL, dims = c(1,2),
                         method = c("ks", "wkde", "sm"), adjust = 1, size = 1, shape = 16, combine = TRUE){


  # Validate object ---------------------------------------------------------

  if(!(is(object, "Seurat") | is(object, "SingleCellExperiment"))){
    stop("Input object must be of 'Seurat' class")
  }

  # Validate dimensions -----------------------------------------------------

  if(length(dims) != 2) stop("Only two dimensions can be plotted")


  # Test existence of reduction ---------------------------------------------

  if(!is.null(reduction)){

    if(is(object, "Seurat")){

      if(!(reduction %in% names(slot(object, "reductions")))){
        stop("No reduction named '", reduction , "' found in object")
      }

    }else if(is(object, "SingleCellExperiment")){

      reductions <- names(reducedDims(object))

      if(!(reduction %in% reductions)){
        stop("No reduction named '", reduction , "' found in object")
      }

    }

  }

  # Set default reductions --------------------------------------------------

  if(is(object, "Seurat")){

    reduction_list <- Reductions(object)
    if(!length(reduction_list)) stop("No reduction has been computed!")

    if(is.null(reduction)){
      reduction <- reduction_list[length(reduction_list)]
    }
      cell_embeddings <- as.data.frame(Embeddings(slot(object, "reductions")[[reduction]]))


  }else if(is(object, "SingleCellExperiment")){

    reduction_list <- names(reducedDims(object))
    if(!length(reduction_list)) stop("No reduction has been computed!")

    if(is.null(reduction)){
      reduction <- reduction_list[length(reduction_list)]
    }
    cell_embeddings <- reducedDim(object, reduction)

  }


  # Validate dimensions
  i <- dims %in% seq_len(ncol(cell_embeddings))
  if(!all(i)){
    missing_dims <- dims[which(!i)]
    stop(paste("Dimension(s) ", missing_dims, " not present in", reduction, "\n  "))
  }


  cell_embeddings <- cell_embeddings[,dims]


  # Set default expression data ---------------------------------------------

  if(is.null(slot)){

    if(is(object, "Seurat")){
      slot <- "data"
    }else if(is(object, "SingleCellExperiment")){
      slot <- "logcounts"
    }

  }


  # Match method argument ---------------------------------------------------

  method <- match.arg(method)

  # Extract metadata --------------------------------------------------------

  if(is(object, "Seurat")){

    metadata <- slot(object, "meta.data")

  }else if(is(object, "SingleCellExperiment")){

    metadata <- as.data.frame(colData(object))
  }


  # Extract features --------------------------------------------------------

  # Test existence of feature in metadata
  if(all(features %in% colnames(metadata))){

    vars <- metadata[,features, drop = FALSE]
    feature_type <- "meta"

    # Extract whole expression data
  }else{

    if(is(object, "Seurat")){

      exp_data <- FetchData(object, vars = features, slot = slot)

    }else if(is(object, "SingleCellExperiment")){

      if(slot == "data") slot <- "logcounts"
      assays <- names(object@assays)

      if(!slot %in% assays) stop(slot, " assay not found")

      exp_data <- Matrix::t(object@assays@data[[slot]])

    }

    # Extract data for input features
    i <- colnames(exp_data) %in% features
    feature_type <- "feature"
    # Test existence of feature in gene expression data
    j <- !features %in% colnames(exp_data)
    if(any(j)) stop("'", paste(features[j], collapse = ", "), "' feature(s) not present in meta.data or expression data")
    vars <- exp_data[,i, drop = FALSE]
    vars <- vars[, features, drop = FALSE]

  }




  # Plot data
  dim_names <- colnames(cell_embeddings)

  if(!is.null(thr)){
    i <- vars > thr
    if(all(i)) stop("All values are above thr = ", thr)

    vars[i] <- 1
    vars[!i] <- 0
  }

  if(ncol(vars) > 1){
    res <- apply(vars, 2, calculate_density, cell_embeddings, method, adjust)
    p <- mapply(plot_density_, as.list(as.data.frame(res)), colnames(res),
                MoreArgs = list(cell_embeddings, dim_names, shape, size, "Density"), SIMPLIFY = FALSE)

    z <- apply(res, 1, prod)
    z_max <- apply(res, 2, max)
    z_max <- Reduce(`*`, z_max)
    joint_label <- paste0(paste(features, "+", sep = ""), collapse = " ")
    pz <- plot_density_(z, joint_label, cell_embeddings, dim_names, shape, size, "Joint density", joint = z_max)


    if(combine){
      p <- wrap_plots(p) + pz
    }else{
      p <- c(p, list(pz))
      names(p) <- c(features, joint_label)
    }

  }else{

    z <- calculate_density(vars[,1], cell_embeddings, method, adjust)
    p <- plot_density_(z, features, cell_embeddings, dim_names, shape, size, "Density")

  }


  p

}


#' @title Weighted 2D kernel density estimation
#' @author Jose Alquicira-Hernandez
#' @param x Dimension 1
#' @param y Dimension 2
#' @param w Weight variable
#' @param h vector of bandwidths for x and y directions.
#' Defaults to normal reference bandwidth (ks::hpi).
#' A scalar value will be taken to apply to both directions.
#' @param adjust Bandwidth adjustment
#' @param n Number of grid points in each direction. Can be scalar or a length-2 integer vector.
#' @param lims The limits of the rectangle covered by the grid as c(xl, xu, yl, yu).
#' @return A list of three components.
#' \itemize{
#' \item \code{x, y}	The x and y coordinates of the grid points, vectors of length n.
#' \item \code{z} An n[1] by n[2] matrix of the weighted estimated density: rows correspond to the value of x, columns to the value of y.
#' }
#' @importFrom Matrix Matrix
#' @importFrom stats dnorm
#' @importFrom methods is
#' @examples
#'
#' set.seed(1)
#' x <- rnorm(100)
#'
#' set.seed(2)
#' y <- rnorm(100)
#'
#' set.seed(3)
#' w <- sample(c(0,1), 100, replace = TRUE)
#'
#' dens <- Nebulosa:::wkde2d(x, y, w)
#'




wkde2d <- function(x, y, w, h, adjust = 1, n = 100, lims = c(range(x), range(y))){

  # Validate values and dimensions
  nx <- length(x)
  if(!all(all(nx == length(y)), all(nx == length(w))))
    stop("data vectors must be the same length")
  if(any(!is.finite(x)) || any(!is.finite(y)))
    stop("missing or infinite values in the data are not allowed")
  if (any(!is.finite(lims)))
    stop("only finite values are allowed in 'lims'")

  h <- c(ks::hpi(x),
         ks::hpi(y))
  h <- h * adjust

  # Get grid
  gx <- seq.int(lims[1L], lims[2L], length.out = n)
  gy <- seq.int(lims[3L], lims[4L], length.out = n)

  # weight
  ax <- outer(gx, x, "-") / h[1L]
  ay <- outer(gy, y, "-") / h[2L]

  w <- Matrix::Matrix(rep(w, n), nrow = n, ncol = nx, byrow = TRUE)

  z <- Matrix::tcrossprod(dnorm(ax) * w, dnorm(ay) * w) /
    (sum(w) * h[1L] * h[2L])

  dens <- list(x = gx, y = gy, z = z)
  dens
}


get_dens <- function(data, dens, method){

  if(method == "ks"){

    ix <- findInterval(data[, 1], dens$eval.points[[1]])
    iy <- findInterval(data[, 2], dens$eval.points[[2]])
    ii <- cbind(ix, iy)
    z <- dens$estimate[ii]

  }else if(method == "wkde"){

    ix <- findInterval(data[, 1], dens$x)
    iy <- findInterval(data[, 2], dens$y)
    ii <- cbind(ix, iy)
    z <- dens$z[ii]

  }else{

    ix <- findInterval(data[, 1], dens$eval.points[, 1])
    iy <- findInterval(data[, 2], dens$eval.points[, 2])
    ii <- cbind(ix, iy)
    z <- dens$estimate[ii]

  }

  z

}


#' @title Estimate weighted kernel density
#' @author Jose Alquicira-Hernandez
#' @param w Vector with weights for each observation
#' @param x Matrix with dimensions where to calculate the density from. Only the first two dimensions will be used
#' @param method Kernel density estimation method:
#' \itemize{
#' \item \code{ks}: Computes density using the \code{kda} function from the \code{ks} package.
#' Ideal for small-medium datasets (up to 70,000 cells in a desktop computer). Use other methods for larger datasets.
#' \item \code{wkde}: Computes density using a modified version of the \code{kde2d} function from the \code{MASS}
#' package to allow weights. Bandwidth selection from the \code{ks} package is used instead.
#' \item \code{sm}: Computes density using the \code{sm.density} function from the \code{sm} package
#' }
#' @param adjust Numeric value to adjust to bandwidth. Default: 1. Not available for \code{ks} method
#' @param map Whether to map densities to individual observations
#' @return If \code{map} is \code{TRUE}, a vector with corresponding densities for each observation is returned. Otherwise,
#' a war a list with the density estimates from the selected method is returned.
#' @importFrom ks kde hpi
#' @importFrom sm sm.density
#' @examples
#'
#' dens <- Nebulosa:::calculate_density(iris[,3],iris[,1:2], method = "wkde")
#'

calculate_density <- function(w, x, method, adjust = 1, map = TRUE){

  if(method == "ks"){

    dens <- kde(x[,1:2],
                w = w / sum(w) * length(w))

  }else if(method == "wkde"){

    dens <- wkde2d(x = x[, 1],
                   y = x[, 2],
                   w =  w / sum(w) * length(w),
                   adjust = adjust)


  }else{
    h1 <- hpi(x[, 1]) * adjust
    h2 <- hpi(x[, 2]) * adjust
    h <-  c(h1, h2)

    dens <- quiet(sm.density(x[, 1:2],
                             weights = w,
                             h = c(h1, h2),
                             structure.2d = "separate",
                             nbins = 0,
                             display = "none"))

  }

  if(map){
    get_dens(x, dens, method)
  }else{
    dens
  }


}

#' @title Plot density estimates
#' @author Jose Alquicira-Hernandez
#' @param z Vector with density values for each cells
#' @param feature Name of the feature being plotted
#' @param cell_embeddings Matrix with cell embeddings
#' @param dim_names Names of the dimensions from the cell embeddings
#' @param shape Geom shape
#' @param size Geom size
#' @param legend_title String used as legend title
#' @param joint Maximum product of joint densities. For visualization purposes.
#' @return A ggplot object
#' @importFrom ggplot2 ggplot aes_string geom_point xlab ylab ggtitle labs guide_legend theme element_text element_line element_rect element_blank scale_color_viridis_c scale_color_gradientn
#' @importFrom viridis viridis

plot_density_ <- function(z, feature, cell_embeddings, dim_names, shape, size, legend_title, joint = 0){

  ggplot(data.frame(cell_embeddings, feature = z)) +
    aes_string(dim_names[1], dim_names[2], color = "feature") +
    geom_point(shape = shape, size = size) +
    xlab(gsub("_", " ", dim_names[1])) +
    ylab(gsub("_", " ", dim_names[2])) +
    ggtitle(feature) +
    labs(color = guide_legend(legend_title)) +
    theme(text = element_text(size = 14),
          panel.background = element_blank(),
          axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(color = "black"),
          axis.line = element_line(size = 0.25),
          strip.background = element_rect(color = "black", fill =  "#ffe5cc")) -> p

  if(joint){
    p <- p + scale_color_gradientn(colors = viridis(10),  limits = c(0, joint))
  }else{
    p <- p + scale_color_viridis_c()
  }

}


quiet <- function(x){
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}
