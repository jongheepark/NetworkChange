#' Plot of network array data
#'
#' Plot network array data using ggplot2 and patchwork for layout.
#'
#' @param Y network array data
#' @param n.graph number of subgraphs. Default is 4.
#' @param node.size node size.  Default is 2.
#' @param node.color node color.  Default is "brown."
#' @param edge.alpha transparency of edge.  Default is 0.5.
#' @param edge.size edge size.  Default is 0.2.
#' @param edge.color edge color.  Default is "grey."
#'
#' @references  Jong Hee Park and Yunkyun Sohn. 2020. "Detecting Structural Change
#' in Longitudinal Network Data." \emph{Bayesian Analysis}. Vol.15, No.1, pp.133-157.
#'
#' @return A patchwork plot object
#'
#' @importFrom patchwork wrap_plots
#' @importFrom network network
#' @importFrom GGally ggnet2
#'
#' @export
#'
#' @examples
#'
#'    \dontrun{
#'    set.seed(1973)
#'    ## generate an array with two constant blocks
#'    Y <- MakeBlockNetworkChange(n=10, shape=1, T=20, type ="split")
#'    plotnetarray(Y)
#'    }


plotnetarray <- function(Y, n.graph = 4, node.size = 2,
                          node.color = "brown", edge.alpha = 0.5,
                          edge.size = 0.2, edge.color = "grey"){
    K <- dim(Y)
    multigraph <- vector("list", n.graph)
    dist <- round(seq(1, K[3], length = n.graph))

    for(g in 1:n.graph){
        net <- network(Y[,,dist[g]], directed = FALSE)
        multigraph[[g]] <- ggnet2(net, node.size = node.size,
                                  node.color = node.color,
                                  edge.size = edge.size, edge.color = edge.color) +
            ggtitle(paste("t =", dist[g])) +
            theme_networkchange() +
            theme(aspect.ratio = 1)
    }
    ## Use patchwork for layout instead of gridExtra
    patchwork::wrap_plots(multigraph, nrow = round(sqrt(n.graph)))
}



