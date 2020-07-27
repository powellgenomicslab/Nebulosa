#' @title Weighted 2D kernel density estimation
#' @author Jose Alquicira-Hernandez
#' @param x Dimension 1
#' @param y Dimension 2
#' @param w Weight variable
#' @param h vector of bandwidths for x and y directions.
#' Defaults to normal reference bandwidth (ks::hpi).
#' A scalar value will be taken to apply to both directions.
#' @param adjust Bandwidth adjustment
#' @param n Number of grid points in each direction. Can be scalar or a
#' length-2 integer vector.
#' @param lims The limits of the rectangle covered by the grid as
#' c(xl, xu, yl, yu).
#' @return A list of three components.
#' \itemize{
#' \item \code{x, y} The x and y coordinates of the grid points, vectors of
#' length n.
#' \item \code{z} An n[1] by n[2] matrix of the weighted estimated density:
#' rows correspond to the value of x, columns to the value of y.
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
#' w <- sample(c(0, 1), 100, replace = TRUE)
#'
#' dens <- Nebulosa:::wkde2d(x, y, w)
wkde2d <- function(x, y, w, h, adjust = 1, n = 100,
                   lims = c(range(x), range(y))) {

    # Validate values and dimensions
    nx <- length(x)
    if (!all(all(nx == length(y)), all(nx == length(w)))) {
        stop("data vectors must be the same length")
    }
    if (any(!is.finite(x)) || any(!is.finite(y))) {
        stop("missing or infinite values in the data are not allowed")
    }
    if (any(!is.finite(lims))) {
        stop("only finite values are allowed in 'lims'")
    }

    h <- c(
        ks::hpi(x),
        ks::hpi(y)
    )
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


get_dens <- function(data, dens, method) {
    if (method == "ks") {
        ix <- findInterval(data[, 1], dens$eval.points[[1]])
        iy <- findInterval(data[, 2], dens$eval.points[[2]])
        ii <- cbind(ix, iy)
        z <- dens$estimate[ii]
    } else if (method == "wkde") {
        ix <- findInterval(data[, 1], dens$x)
        iy <- findInterval(data[, 2], dens$y)
        ii <- cbind(ix, iy)
        z <- dens$z[ii]
    } else {
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
#' @param x Matrix with dimensions where to calculate the density from. Only
#' the first two dimensions will be used
#' @param method Kernel density estimation method:
#' \itemize{
#' \item \code{ks}: Computes density using the \code{kda} function from the
#'  \code{ks} package.
#' Ideal for small-medium datasets (up to 70,000 cells in a desktop computer).
#'  Use other methods for larger datasets.
#' \item \code{wkde}: Computes density using a modified version of the
#'  \code{kde2d} function from the \code{MASS}
#' package to allow weights. Bandwidth selection from the \code{ks} package
#'  is used instead.
#' \item \code{sm}: Computes density using the \code{sm.density} function from
#'  the \code{sm} package
#' }
#' @param adjust Numeric value to adjust to bandwidth. Default: 1. Not available
#'  for \code{ks} method
#' @param map Whether to map densities to individual observations
#' @return If \code{map} is \code{TRUE}, a vector with corresponding densities
#'  for each observation is returned. Otherwise,
#' a war a list with the density estimates from the selected method is returned.
#' @importFrom ks kde hpi
#' @importFrom sm sm.density
#' @examples
#'
#' dens <- Nebulosa:::calculate_density(iris[, 3], iris[, 1:2], method = "wkde")
calculate_density <- function(w, x, method, adjust = 1, map = TRUE) {
    if (method == "ks") {
        dens <- kde(x[, c(1, 2)],
                    w = w / sum(w) * length(w)
        )
    } else if (method == "wkde") {
        dens <- wkde2d(
            x = x[, 1],
            y = x[, 2],
            w = w / sum(w) * length(w),
            adjust = adjust
        )
    } else {
        h1 <- hpi(x[, 1]) * adjust
        h2 <- hpi(x[, 2]) * adjust
        h <- c(h1, h2)

        dens <- quiet(sm.density(x[, c(1, 2)],
                                 weights = w,
                                 h = c(h1, h2),
                                 structure.2d = "separate",
                                 nbins = 0,
                                 display = "none"
        ))
    }

    if (map) {
        get_dens(x, dens, method)
    } else {
        dens
    }
}
