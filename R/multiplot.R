#' Printing multiple ggplots in oone file
#' 
#' Print multiple ggplots in one file. Slightly modified for packaging from the original version in the web. 
#'
#' @param ...  A list of ggplot objects separated by commas.
#' @param plotlist A list of ggplot objects
#' @param cols The number of columns.
#' @param layout A matrix specifying the layout. If present, 'cols' is ignored.
#'
#' @author \url{http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/}
#' 
#' @return A plot object
#'
#' @export
#'
#'
multiplot <- function (..., plotlist = NULL, cols = 1, layout = NULL) 
{
    plots <- c(list(...), plotlist)
    numPlots = length(plots)
    if (is.null(layout)) {
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)), 
            ncol = cols, nrow = ceiling(numPlots/cols), byrow=TRUE)
    }
    if (numPlots == 1) {
        print(plots[[1]])
    }
    else {
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), 
            ncol(layout))))
        for (i in 1:numPlots) {
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row, 
                layout.pos.col = matchidx$col))
        }
    }
}
