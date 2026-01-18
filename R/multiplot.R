#' Printing multiple ggplots in one file (DEPRECATED)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function has been deprecated in favor of [patchwork::wrap_plots()].
#' Please use `patchwork::wrap_plots(plotlist, ncol = cols)` instead.
#'
#' @param ...  A list of ggplot objects separated by commas.
#' @param plotlist A list of ggplot objects
#' @param cols The number of columns.
#' @param layout A matrix specifying the layout. If present, 'cols' is ignored.
#'
#' @return A plot object
#'
#' @importFrom patchwork wrap_plots
#'
#' @export
#'
multiplot <- function(..., plotlist = NULL, cols = 1, layout = NULL) {
    .Deprecated("patchwork::wrap_plots",
                msg = "multiplot() is deprecated. Please use patchwork::wrap_plots() instead.")

    plots <- c(list(...), plotlist)
    numPlots <- length(plots)

    if (numPlots == 1) {
        return(plots[[1]])
    }

    ## Use patchwork for modern layout
    patchwork::wrap_plots(plots, ncol = cols)
}
