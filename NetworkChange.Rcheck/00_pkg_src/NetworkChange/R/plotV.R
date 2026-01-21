## Avoid R CMD check notes for ggplot2 NSE
utils::globalVariables(c("Time", "Value", "Dimension"))

#' Plot of layer-specific network generation rules.
#'
#' Plot layer-specific network generation rules using ggplot2.
#' Uses colorblind-friendly viridis palette for publication-quality output.
#'
#' @param OUT Output of networkchange objects.
#' @param main The title of plot
#' @param point_size Point size (default: 3)
#' @param line_size Line width (default: 1)
#' @return A ggplot2 plot object
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline labs
#' @importFrom ggplot2 scale_color_viridis_d
#' @importFrom tidyr pivot_longer
#'
#' @export
#'
#' @examples
#'
#'    \dontrun{set.seed(1973)
#'    \## generate an array with two constant blocks
#'    Y <- MakeBlockNetworkChange(n=10, shape=10, T=40, type ="constant")
#'    out0 <- NetworkStatic(Y, R=2, mcmc=10, burnin=10,
#'    verbose=10, UL.Normal = "Orthonormal")
#'    \## latent node positions
#'    plotV(out0)
#'    }



plotV <- function(OUT, main = "", point_size = 3, line_size = 1) {
    Vmat <- attr(OUT, "Vmat")
    R <- attr(OUT, "R")
    Y <- attr(OUT, "Z")
    Time_points <- dim(Vmat)[2] / R

    ## Compute posterior means
    vmean <- apply(Vmat, 2, mean)
    V <- matrix(vmean, Time_points, R)

    ## Create data frame for ggplot
    df <- data.frame(
        Time = rep(1:Time_points, R),
        Value = as.vector(V),
        Dimension = factor(rep(paste0("Dimension ", 1:R), each = Time_points))
    )

    ## Create ggplot
    p <- ggplot(df, aes(x = Time, y = Value, color = Dimension, group = Dimension)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        geom_line(linewidth = line_size, alpha = 0.8) +
        geom_point(size = point_size, alpha = 0.9) +
        scale_color_viridis_d(option = "D", end = 0.9) +
        labs(
            title = if(nchar(main) > 0) main else "Layer-Specific Network Generation Rules",
            x = "Time",
            y = expression(V),
            color = "Latent\nDimension"
        ) +
        theme_networkchange() +
        theme(legend.position = "bottom")

    return(p)
}
