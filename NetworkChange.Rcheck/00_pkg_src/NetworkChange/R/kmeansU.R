#' K-mean clustering of latent node positions
#'
#' K-mean clustering of latent node positions
#' 
#' @param out Output of networkchange objects.
#' @param R Number of latent space dimensions
#' @param n.cluster Number of latent cluster
#' @param layer Layer id for the cluster analysis
#' @param main Title
#'
#' @return A plot object
#' 
#' @export
#'
#' @examples
#'
#'    \dontrun{set.seed(1973)
#'    ## generate an array with two constant blocks
#'    Y <- MakeBlockNetworkChange(n=10, shape=10, T=10, type ="constant")
#'    out0 <- NetworkStatic(Y, R=2, mcmc=10, burnin=10,
#'    verbose=10, UL.Normal = "Orthonormal")
#'    ## latent node positions
#'    kmeansU(out0)
#'    }



kmeansU <- function(out, R = 2, n.cluster=3, layer = 1,
                    main=""){
    MU <- out
    Z <- attr(out, "Z")
    K <- dim(Z)
    n <- dim(Z)[1]
    V <- matrix(apply(attr(out, "Vmat"), 2, mean), K[3], 2)
    tmp <- eigen(MU[,,layer])
    
    L <- diag(tmp$val[order(-abs(tmp$val))[seq(1, R, length = R)]]/n, nrow = R)
    U <- tmp$vec[, order(-abs(tmp$val))[seq(1, R, length = R)], 
                 drop = FALSE] * sqrt(n)
    cls <- kmeans(U, n.cluster, nstart = 20)
    
    ggdata <- data.frame(U1 = U[,1], U2 = U[,2], V1 = V[1], V2 = V[2])
    ggdata$cluster <- as.factor(cls$cluster)
    return(ggdata)    
}
