#' @title Plot density estimates
#' @author Jose Alquicira-Hernandez
#' @param z Vector with density values for each cells
#' @param feature Name of the feature being plotted
#' @param cell_embeddings Matrix with cell embeddings
#' @param dim_names Names of the dimensions from the cell embeddings
#' @param shape Geom shape
#' @param size Geom size
#' @param legend_title String used as legend title
#' @param pal String specifying the viridis color palette to use
#' @param ... Further scale arguments passed to scale_color_viridis_c
#' @return A ggplot object
#' @importFrom ggplot2 ggplot aes_string geom_point xlab ylab ggtitle labs
#' guide_legend theme element_text element_line element_rect element_blank
#' scale_color_viridis_c scale_color_gradientn
plot_density_ <- function(z, feature, cell_embeddings, dim_names, shape, size,
                          legend_title,
                          pal = c(
                              "viridis", "magma", "cividis",
                              "inferno", "plasma", "mako",
                              "rocket","turbo"
                          ), ...) {
    p <- ggplot(data.frame(cell_embeddings, feature = z)) +
        aes_string(dim_names[1], dim_names[2], color = "feature") +
        geom_point(shape = shape, size = size) +
        xlab(gsub("_", " ", dim_names[1])) +
        ylab(gsub("_", " ", dim_names[2])) +
        ggtitle(feature) +
        labs(color = guide_legend(legend_title)) +
        theme(
            text = element_text(size = 14),
            panel.background = element_blank(),
            axis.text.x = element_text(color = "black"),
            axis.text.y = element_text(color = "black"),
            axis.line = element_line(size = 0.25),
            strip.background = element_rect(color = "black", fill = "#ffe5cc")
        )

    pal <- match.arg(pal)
    p + scale_color_viridis_c(option = pal, ...)

}
