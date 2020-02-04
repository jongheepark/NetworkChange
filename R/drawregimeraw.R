#' Plot of network by hidden regime
#'
#' Plot latent node cluster
#' 
#' @param mcmcout NetworkChange output
#' @param Y Input raw data
#'
#' @references   Jong Hee Park and Yunkyun Sohn. 2019. "Detecting Structural Change
#' in Longitudinal Network Data." \emph{Bayesian Analysis}. Forthcoming.
#'
#' @importFrom igraph graph.adjacency get.data.frame
#' @importFrom qgraph qgraph
#' 
#' @return A plot object
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
#'    drawRegimeRaw(fit, newY)
#'    }

drawRegimeRaw <- function(mcmcout, Y){    
    m <- attr(mcmcout, "m")
    mcmc <- attr(mcmcout, "mcmc")
    U <- attr(mcmcout, "U")
    V <- attr(mcmcout, "V")
    R <- attr(mcmcout, "R")
    K <- dim(Y);
    T <- K[3]

    ## median estimate of hidden states
    median.s <- ceiling(apply(attr(mcmcout, "Smat"), 2, median))

    ## dim names for time layer
    if(is.null(dimnames(Y)[[3]])){
        Year <- 1:T
    }else{
        Year <- as.numeric(dimnames(Y)[[3]])
    }
    y <- ts(Year, start=Year[1])
    ns <- m + 1
    time <- K[[3]]
    
    ## regime specific country names
    names <- as.list(rep(NA, ns))
    if(is.null(dimnames(Y)[[3]])){
        for(t in 1:ns){
            names[[t]] <- 1:K[1]
        }
    }else{
        for(t in 1:ns){
            names[[t]] <- dimnames(Y)[[1]]
        }
    }
 
    end <- c(which(diff(median.s) == 1), length(Year))
    start <- c(1, which(diff(median.s) == 1)+1)

    net <- array(NA, dim=c(dim(Y)[1], dim(Y)[2], ns))
    for(t in 1:ns){
        net[,,t] <- apply(Y[,,start[t]:end[t]], 1:2, sum)
    }

    ## title.list
    title.list <- time.period <- as.list(rep(NA, ns))
    p.list <- as.list(rep(NA, ns))
    median.s <- ceiling(apply(attr(mcmcout, "Smat"), 2, median))
    start = 1; 
    for(i in 1:ns){
        time.period[[i]] <- paste0(range(Year[median.s == i])[1], "-", range(Year[median.s == i])[2])
        title.list[[i]] <- paste0("Regime ", i, " (", time.period[[i]], ")")
    }
    
    par(mfrow=c(1,ns)); par (mar=c(1,1,2,1), mgp=c(2,.7,0), tck=.02)
    ver.size <- 5
    for(t in 1:ns){
        g <- graph.adjacency(net[,,t], weighted=TRUE, mode="undirected")
        df <- get.data.frame(g)
        qgraph(df, gray=TRUE, labels= names[[t]],
               label.cex = 1, arrows=FALSE, title=title.list[[t]], normalize=FALSE)
    }
       
}
