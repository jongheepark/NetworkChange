#' Contour plot of latent node positions
#'
#' Draw a contour plot of latent node positions
#'
#' @param out Output of networkchange objects.
#'
#' @param n.cluster The number of latent block in network for k-means clustering
#'
#' @param vertex.names The vertex names
#' 
#' @param main The title of plot

#' @return A plot object
#'
#' @export
#'
#'
#' @examples
#'
#'    \dontrun{set.seed(1973)
#'     ## generate an array with three blocks
#' set.seed(11173)
#'
#' ## The number of node in each block is
#' N <- 5
#'
#' ## Generate block-splitting network time series data
#' Yarr <- MakeBlockNetworkChange(n=N, break.point = .5, base.prob=.05, block.prob=.5, shape=1, T=20, type ="split")
#'
#' ## Pick is a homophily network with 3 blocks
#' Y1 <- Yarr[,,16]
#'
#' ## Y2 is a mirror image of Y1, which is heterophilous
#' Y2 <- 1-Yarr[,,20]
#' diag(Y2) <- 0
#'
#' ## Combine them into a multilayer network array
#' Y <- abind(Y1, Y2, along=3)
#'
#' G <- 100 ## Small mcmc scans to save time
#' out0 <- NetworkStatic(Y, R=2,  mcmc = G, burnin = G, constant=FALSE, verbose= G, degree.normal="eigen")
#' plot3 <- Kmeans(out0, n.cluster=3, main="(C) Recovered Three Blocks")
#'
#' require (sna)
#' g1 <- network(Y[,,1], directed = FALSE)
#' g2 <- network(Y[,,2], directed = FALSE)
#'
#' require(ggnet)
#' plot1 <- ggnet2(g1, node.size = 4, label.size=3, label = plot3$data$names, color = rep(c("red", "green", "blue"), each=N),
#' label.color = "white", alpha = 0.5) + ggtitle("(A) Network Layer 1")
#' plot2 <- ggnet2(g2, node.size = 4, label.size=3, label = plot3$data$names, color = rep(c("red", "green", "blue"), each=N),
#' label.color = "white", alpha = 0.5) + ggtitle("(B) Network Layer 2")
#' multiplot(plot1, plot2, plot3, cols=3)
#'
#' }

Kmeans <- function(out, n.cluster=2, vertex.names=NULL, main="",...){
    MU <- out
    R <- attr(out, "R")
    K <- dim(MU)
    n <- K[1]
    V <- matrix(apply(attr(out, "Vmat"), 2, mean), K[3], 2)
    tmp <- eigen(MU[,,1])
    
    L <- diag(tmp$val[order(-abs(tmp$val))[seq(1, R, length = R)]]/n, nrow = R)
    U <- tmp$vec[, order(-abs(tmp$val))[seq(1, R, length = R)], 
                 drop = FALSE] * sqrt(n)
    cls <- kmeans(U, n.cluster, nstart = 20)
    
    ggdata <- data.frame(U1 = U[,1], U2 = U[,2], V1 = V[1], V2 = V[2])
    ggdata$cluster <- as.factor(cls$cluster)
    if(is.null(vertex.names)){
        vertex.names=1:n
    }
    ggdata$names <- vertex.names
        
    g.out <- ggplot(ggdata, aes(U1, U2, color = ggdata$cluster, label= ggdata$names)) + geom_point(size=5) +
        geom_text(size = 3, label=ggdata$names, color = "white") + theme(legend.position="none") +
        ggtitle(main)
    
    return(g.out)    
}
