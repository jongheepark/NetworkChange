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
#'    ## generate an array with two constant blocks
#'    Y <- MakeBlockNetworkChange(n=10, shape=10, T=40, type ="constant")
#'    out0 <- NetworkStatic(Y, R=2, mcmc=10, burnin=10,
#'    verbose=10, UL.Normal = "Orthonormal")
#'    ## latent node positions
#'    plotV(out0)
#'    }



plotV <- function (OUT, main = "", cex = 2) 
{
  Vmat <- attr(OUT, "Vmat")
  Y <- attr(OUT, "Z")
  T <- dim(Vmat)[2]/2
  my.cols = rainbow(T)
  Vmat1 <- Vmat[, 1:T]
  Vmat2 <- Vmat[, (T + 1):(2 * T)]
  V1 <- apply(Vmat1, 2, mean)
  V2 <- apply(Vmat2, 2, mean)
  par(mar = c(5, 4, 2.4, 2.2))
  plot(1:T, V1, type = "n", main = "", ylim=range(c(V1, V2)), 
       ylab = expression(V), xlab = "Time", xaxt = "n", yaxt = "n")
  axis(1); axis(2); grid( col="grey40")
  abline(h=0, lty=3)
  lines(1:T, V1, lwd=2, col = addTrans("red",100))
  lines(1:T, V2, lwd=2, col = addTrans("blue",100))
  points(1:T, V1, cex = cex, pch=19, col = addTrans("red",150))
  points(1:T, V2, cex = cex, pch=21, col = addTrans("blue",150))
  add_legend("top", legend=c("1st", "2nd"), pch = c(19, 21), 
             lwd=1, lty=c(1, 1), bty="n",
             col = c(addTrans("red",100), addTrans("blue",100)),
             horiz=TRUE, cex=1)
}
