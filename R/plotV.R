#' Plot of layer-specific network generation rules.
#'
#' Plot layer-specific network generation rules.
#' 
#' @param OUT Output of networkchange objects.
#' @param main The title of plot
#' @param cex point size
#' @return A plot object
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



plotV <- function (OUT, main = "", cex = 2) {
    ##   par(mar = c(5, 4, 2.4, 2.2))
    
    Vmat <- attr(OUT, "Vmat")
    R <- attr(OUT, "R")
    Y <- attr(OUT, "Z")
    T <- dim(Vmat)[2]/R
    my.cols = rainbow(R)
    vmean <- apply(Vmat, 2, mean)
    V <- matrix(vmean, T, R)
    ## Vmat1 <- Vmat[, 1:T]
    ## Vmat2 <- Vmat[, (T + 1):(2 * T)]
    ## V1 <- apply(Vmat1, 2, mean)
    ## V2 <- apply(Vmat2, 2, mean)
    
    plot(1:T, V[,1], type = "n", main = "", ylim=range(V), 
         ylab = expression(V), xlab = "Time", xaxt = "n", yaxt = "n")
    axis(1); axis(2); grid( col="grey40")
    abline(h=0, lty=3)
    for(i in 1:R){
        lines(1:T, V[,i], lwd=2, col = addTrans(my.cols[i],100))
        points(1:T,V[,i], cex = cex, pch=19, col = addTrans(my.cols[i],150))
    }
    legend("bottomright", legend=paste0("Dimension ", 1:R), pch = c(19), 
               lwd=1, lty=c(1), 
               col = my.cols,
               inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="n", cex=1)
}
