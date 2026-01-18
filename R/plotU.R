#' Plot of latent node positions
#'
#' Plot latent node positions using ggplot2 with patchwork for multi-panel layouts.
#' Uses colorblind-friendly palettes and cross-platform compatible fonts.
#'
#' @param OUT Output of networkchange objects.
#' @param Time Starting of the time period. If NULL, 1.
#' @param main The title of plot
#' @param names Node names. If NULL, use natural numbers.
#' @param label.prob Label print threshold. 0.9 is the default.
#' @return A plot object
#'
#' @importFrom ggrepel geom_label_repel
#' @importFrom patchwork wrap_plots
#'
#' @export
#'
#' @examples
#'
#'    \dontrun{
#'    set.seed(1973)
#'    \## generate an array with two constant blocks
#'    Y <- MakeBlockNetworkChange(n=10, shape=10, T=40, type ="constant")
#'    out0 <- NetworkStatic(Y, R=2, mcmc=10, burnin=10,
#'    verbose=10, UL.Normal = "Orthonormal")
#'    \## latent node positions
#'    plotU(out0)
#'    }


plotU <- function(OUT, Time=NULL, names=NULL, main=NULL, label.prob=0.9){
    m <- attr(OUT, "m")
    mcmc <- attr(OUT, "mcmc")
    Z <- attr(OUT, "Z")
    K <- dim(Z)
    R <- attr(OUT, "R")
    if(is.null(Time)){
        y <- Time <- ts(1:K[3])
    } else{
        y <- ts(Time)
    }
    if(is.null(names)){
        names <- 1:K[1]
    }
    ns <- m + 1
    First <- Second <- Size <- Names <- Name <- NA

    if(m == 0){
        x <- OUT
        x.mean <- apply(x, 1:2, mean)
        tmp <- eigen(x.mean)
        U <- tmp$vec[, order(-abs(tmp$val))[seq(1, R, length = R)],
                     drop = FALSE] * sqrt(K[1])
        df <- data.frame(First = U[,1], Second = U[,2],
                         Size = sqrt((U[,1])^2 + (U[,2])^2),
                         Names = names)
        if(is.null(main)){
            title <- "Latent Space of No Break Model"
        } else {
            title <- main
        }
        df$Name <- ifelse(df$Size > quantile(df$Size, probs=label.prob), df$Names, NA)
        p <- ggplot(df, aes(x=First, y=Second, label=Name)) +
            geom_point(size = df$Size+1, colour = alpha("red", 0.2)) +
            ggtitle(title) +
            geom_label_repel(aes(label=Name),
                             color="navy",
                             box.padding=0.35, point.padding=0.5, segment.color="grey70",
                             fontface="bold", size=2) +
            theme_networkchange()
        return(p)
    } else {
        ## Multi-regime plots
        U.list <- df.list <- p.list <- title.list <- time.period <- vector("list", ns)
        median.s <- apply(attr(OUT, "Smat"), 2, median)

        for(i in 1:ns){
            time.period[[i]] <- paste0(range(Time[median.s == i])[1], "-", range(Time[median.s == i])[2])
            x <- OUT[,,median.s == i]
            x.mean <- apply(x, 1:2, mean)
            tmp <- eigen(x.mean)
            U.list[[i]] <- tmp$vec[, order(-abs(tmp$val))[seq(1, R, length = R)],
                                   drop = FALSE] * sqrt(K[1])
            df.list[[i]] <- data.frame(First = U.list[[i]][,1], Second = U.list[[i]][,2],
                                       Size = sqrt((U.list[[i]][,1])^2 + (U.list[[i]][,2])^2),
                                       Names = names)
            if(is.null(main)){
                title.list[[i]] <- paste0("Latent Space of Regime ", i, " (", time.period[[i]], ")")
            } else {
                title.list[[i]] <- main
            }

            df.list[[i]]$Name <- ifelse(df.list[[i]]$Size >
                                        quantile(df.list[[i]]$Size, probs=label.prob),
                                        df.list[[i]]$Names, NA)

            p.list[[i]] <- ggplot(df.list[[i]], aes(x=First, y=Second, label=Name)) +
                geom_point(size = df.list[[i]]$Size+1, colour = alpha("red", 0.2)) +
                ggtitle(title.list[[i]]) +
                geom_label_repel(aes(label=Name),
                                 color="navy",
                                 box.padding=0.35, point.padding=0.5, segment.color="grey70",
                                 fontface="bold", size=2) +
                theme_networkchange()
        }

        ## Use patchwork for layout instead of multiplot
        ncols <- if(ns < 5) ns else ceiling(ns/2)
        combined_plot <- patchwork::wrap_plots(p.list, ncol = ncols)
        return(combined_plot)
    }
}
