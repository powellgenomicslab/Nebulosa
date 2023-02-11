.validate_dimensions <- function(dims){
    if (length(dims) != 2) stop("Only two dimensions can be plotted")
}


.search_dimensions <- function(dims, cell_embeddings, reduction) {
    # Validate dimensions
    i <- dims %in% seq_len(ncol(cell_embeddings))
    if (!all(i)) {
        missing_dims <- dims[which(!i)]
        stop("Dimension(s) ", missing_dims, " not present in", reduction, "\n")
    }

    cell_embeddings <- cell_embeddings[, dims]
    cell_embeddings
}



.extract_feature_data <- function(exp_data, features) {
    # Extract data for input features
    i <- colnames(exp_data) %in% features

    # Test existence of feature in gene expression data
    j <- !features %in% colnames(exp_data)
    if (any(j)) {
        stop(
            "'", paste(features[j], collapse = ", "),
            "' feature(s) not present in meta.data or expression data"
        )
    }
    vars <- exp_data[, i, drop = FALSE]
    vars <- vars[, features, drop = FALSE]
    vars
}


#' @importFrom patchwork wrap_plots
.plot_final_density <- function(vars, cell_embeddings, features, joint, method,
                                adjust, shape, size, pal, combine, raster, ...) {
    dim_names <- colnames(cell_embeddings)
    if (ncol(vars) > 1) {
        res <- apply(vars, 2, calculate_density, 
                     cell_embeddings, method, adjust)
        p <- mapply(plot_density_, as.list(as.data.frame(res)), colnames(res),
                    MoreArgs = list(cell_embeddings, dim_names, shape, size,
                                    "Density", pal = pal, raster, ...), 
                    SIMPLIFY = FALSE)

        if(joint){

            z <- apply(res, 1, prod)
            joint_label <- paste0(paste(features, "+", sep = ""), 
                                  collapse = " ")
            pz <- plot_density_(z, joint_label, cell_embeddings,
                                dim_names,
                                shape,
                                size,
                                "Joint density",
                                pal = pal,
                                raster,
                                ...
            )


            if (combine) {
                p <- wrap_plots(p) + pz
            } else {
                p <- c(p, list(pz))
                names(p) <- c(features, joint_label)
            }
        }else{
            if (combine) {
                p <- wrap_plots(p)
            } else {
                names(p) <- features
            }
        }

    } else {
        z <- calculate_density(vars[, 1], cell_embeddings, method, adjust)
        p <- plot_density_(z,
                           features,
                           cell_embeddings,
                           dim_names,
                           shape,
                           size,
                           "Density",
                           pal = pal,
                           raster,
                           ...
        )
    }
    p
}
