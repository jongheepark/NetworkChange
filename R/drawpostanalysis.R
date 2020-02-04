#' Plot of latent node cluster
#'
#' Plot latent node cluster
#' 
#' @param mcmcout NetworkChange output
#' @param Y Input raw data
#' @param point.cex node point size.  Default is 3. 
#' @param text.cex node label size.  Default is 3.
#' @param segment.size segment size.  Default is 0.1.
#' @param n.cluster number of cluster. Default is 3. 
#'
#' @references   Jong Hee Park and Yunkyun Sohn. 2020. "Detecting Structural Change
#' in Longitudinal Network Data." \emph{Bayesian Analysis}. Forthcoming.
#' 
#' @return A plot object
#'
#' @importFrom reshape melt
#' @importFrom ggrepel geom_text_repel
#' @importFrom rlang .data
#' 
#' @export
#'
#' @examples
#'
#'    \dontrun{
#'    set.seed(1973)
#'    ## generate an array with two constant blocks
#'    data(MajorAlly)
#'    Y <- MajorAlly
#'    fit <- NetworkChange(newY, R=2, m=2, mcmc=G, initial.s = initial.s,
#'           burnin=G, verbose=0, v0=v0, v1=v1)
#'    drawPostAnalysis(fit, Y, n.cluster=c(4, 4, 3))
#'    }

drawPostAnalysis <- function(mcmcout, Y, point.cex=3,  text.cex=3,
                              segment.size = 0.1, n.cluster = NULL){
  m <- attr(mcmcout, "m")
  mcmc <- attr(mcmcout, "mcmc")
  U <- attr(mcmcout, "U")
  V <- attr(mcmcout, "V")
  R <- attr(mcmcout, "R")
  K <- dim(Y);
  T <- K[3]
  if(is.null(dimnames(Y))){
     dimnames(Y)[[1]]<-c(1:dim(Y)[1])
     dimnames(Y)[[2]]<-c(1:dim(Y)[2])
     dimnames(Y)[[3]]<-c(1:dim(Y)[3])
  }
  median.s <- ceiling(apply(attr(mcmcout, "Smat"), 2, median))
  Year <- as.numeric(dimnames(Y)[[3]])
  y <- ts(Year, start=Year[1])
  ns <- m + 1
  time <- K[[3]]
  Year <- as.numeric(dimnames(Y)[[3]])

  if(is.null(n.cluster)){
    n.cluster <- rep(3, ns)
  }
  ## regime specific country names
  names <- as.list(rep(NA, ns))
  for(t in 1:ns){
    names[[t]] <- dimnames(Y)[[1]]
  }

  ## load network data
  end <- c(which(diff(median.s) == 1), time)
  start <- c(1, which(diff(median.s) == 1)+1)
  N <- K[1]

  net <- array(NA, dim=c(N, N, ns))
  D <- t(sapply(1:time, function(t){V[t, order(V[t,], decreasing=TRUE)]}))
  Dmat <- as.data.frame(cbind(Year, D))
  colnames(Dmat) <- c("Year", "1st", "2nd")
  D.long <- reshape::melt(Dmat, id.vars="Year")
  D.regime <- matrix(NA, ns, R)

  for(t in 1:ns){
    net[,,t] <- apply(Y[,,start[t]:end[t]], 1:2, sum)
    D.regime[t, ] <- apply(Dmat[start[t]:end[t], 2:3], 2, mean)
  }

  ## multiplot object
  U.list <- .df.list <- title.list <- time.period <- as.list(rep(NA, ns))
  p.list <- as.list(rep(NA, ns))
  median.s <- ceiling(apply(attr(mcmcout, "Smat"), 2, median))
  First <- Second <- Size <- Names <- Cluster <- NULL
  for(i in 1:ns){
    time.period[[i]] <- paste0(range(Year[median.s == i])[1], "-", range(Year[median.s == i])[2])
    U.list[[i]] <- U[[i]]
    cls <- kmeans(U.list[[i]], n.cluster[i], nstart = 20)
    ## U.list[[i]] <- matrix(apply(out[[i]], 2, mean), dim(Y)[1], R)
    .df.list[[i]] <- data.frame(First = U.list[[i]][,1], Second = U.list[[i]][,2],
                               Size = sqrt((U.list[[i]][,1])^2 + (U.list[[i]][,2])^2),
                               Names = names[[i]],
                               Cluster = factor(cls$cluster))
    title.list[[i]] <- paste0("Regime ", i, " (", time.period[[i]], ")")
    p.list[[i]] <- ggplot(.df.list[[i]], aes(x=First, y = Second, label=Names, color=Cluster)) +
      geom_point(aes(colour = Cluster, alpha=1/2), size = point.cex, show.legend=F) +
      scale_size_continuous(guide = FALSE) +
      ggtitle(title.list[[i]]) +
      labs(x = paste("v[1] = ", round(D.regime[i, 1], 2))) +
      labs(y = paste("v[2] = ", round(D.regime[i, 2], 2)))  +
      ggrepel::geom_text_repel(size = text.cex, segment.size = segment.size, segment.color = alpha("navy", 1/6),
                      colour = "navy", aes(label = Names)) +
      theme(axis.title=element_text(size=8), plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5))

  }
  attr(p.list, "title.list") <- title.list
  attr(p.list, "D.long") <- D.long
  attr(p.list, "median.s") <- median.s
  return(p.list)
}
