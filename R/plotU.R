#' Plot of latent node positions
#'
#' Plot latent node positions.
#' 
#' @param OUT Output of networkchange objects.
#' @param Year Starting of the time period. If NULL, 1.
#' @param main The title of plot
#' @param names Node names. If NULL, use natural numbers.
#'
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
#'    plotU(out0)
#'    }


plotU <- function(OUT, Year=NULL, names=NULL, main=""){
    m <- attr(OUT, "m")
    mcmc <- attr(OUT, "mcmc")
    Z <- attr(OUT, "Z")
    K <- dim(Z)
    R <- attr(OUT, "R")
    if(is.null(Year)){
        y <- Year <- ts(1:K[3])
    } else{
        y <- ts(Year)
    }
    if(is.null(names)){
        names <- 1:K[1]
    } 
    ns <- m + 1
    First <- Second <- Size <- Names <- NA
    if(m == 0){
        x <- OUT
        x.mean <-  apply(x, 1:2, mean)
        tmp <- eigen(x.mean)
        U <- tmp$vec[, order(-abs(tmp$val))[seq(1, R, length = R)], 
                     drop = FALSE] * sqrt(K[1])
        df <- data.frame(First = U[,1], Second = U[,2],
                         Size = sqrt((U[,1])^2 + (U[,2])^2),
                         Names = names)
        title <- paste0("Latent Space of No Break Model")
        ggplot(df, aes(x=First, y = Second, label=Names)) + geom_point(size = df$Size+1, colour = alpha("red", 1/5)) +
            ggtitle(title) + geom_text(size = df$Size, colour = "navy") + 
                theme(plot.title = element_text(lineheight=.8, face="bold"))
    }
    else{      
        ## plot
        U.list <- df.list <- p.list <- title.list <- time.period <- as.list(rep(NA, ns))
        median.s <- apply(attr(OUT, "Smat"), 2, median)
        
        for(i in 1:ns){
            time.period[[i]] <- paste0(range(Year[median.s == i])[1], "-", range(Year[median.s == i])[2])
            x <- OUT[[i]][,,median.s == i]
            x.mean <-  apply(x, 1:2, mean)
            tmp <- eigen(x.mean)
            U.list[[i]] <- tmp$vec[, order(-abs(tmp$val))[seq(1, R, length = R)], 
                                   drop = FALSE] * sqrt(K[1])
            ## U.list[[i]] <- matrix(apply(out[[i]], 2, mean), dim(Y)[1], R)
            df.list[[i]] <- data.frame(First = U.list[[i]][,1], Second = U.list[[i]][,2],
                                       Size = sqrt((U.list[[i]][,1])^2 + (U.list[[i]][,2])^2),
                                       Names = names)
            title.list[[i]] <- paste0("Latent Space of Break ", i, " (", time.period[[i]], ")")
            p.list[[i]] <- ggplot(df.list[[i]], aes(x=First, y = Second, label=Names)) +
                geom_point(size = df.list[[i]]$Size+1, colour = alpha("red", 1/5)) +
                    ggtitle(title.list[[i]]) +
                        ## geom_text(colour = "navy", aes(label = Names)) + 
                        geom_text(size = df.list[[i]]$Size, colour = "navy", aes(label = Names)) + 
                            theme(plot.title = element_text(lineheight=.8, face="bold"))
        }
        if(ns < 5){
             multiplot(plotlist = p.list, cols=ns)
        }
        else{
            multiplot(plotlist = p.list, cols=ceiling(ns/2))
        }
    }
}
