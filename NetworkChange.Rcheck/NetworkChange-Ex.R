pkgname <- "NetworkChange"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('NetworkChange')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("BreakDiagnostic")
### * BreakDiagnostic

flush(stderr()); flush(stdout())

### Name: BreakDiagnostic
### Title: Detect a break number using different metrics
### Aliases: BreakDiagnostic

### ** Examples

   ## Not run: 
##D    set.seed(19333)
##D    ## Generate an array (15 by 15 by 20) with a block merging transition
##D    Y <- MakeBlockNetworkChange(n=5, T=20, type ="merge")
##D 
##D    ## Fit 3 models (no break, one break, and two break) for break number detection 
##D    detect <- BreakDiagnostic(Y, R=2, break.upper = 2)
##D    
##D    ## Look at the graph
##D    detect[[1]]; print(detect[[2]])
##D 
## End(Not run)






cleanEx()
nameEx("BreakPointLoss")
### * BreakPointLoss

flush(stderr()); flush(stdout())

### Name: BreakPointLoss
### Title: Compute the Average Loss of Hidden State Changes from Expected
###   Break Points
### Aliases: BreakPointLoss

### ** Examples

   ## Not run: 
##D    set.seed(1973)
##D    ## Generate an array (30 by 30 by 40) with block transitions
##D    from 2 blocks to 3 blocks
##D    Y <- MakeBlockNetworkChange(n=10, T=40, type ="split")
##D    G <- 100 ## Small mcmc scans to save time
##D 
##D    ## Fit multiple models for break number detection using Bayesian model comparison
##D    out0 <- NetworkStatic(Y, R=2, mcmc=G, burnin=G, verbose=G, Waic=TRUE)
##D    out1 <- NetworkChange(Y, R=2, m=1, mcmc=G, burnin=G, verbose=G, Waic=TRUE)
##D    out2 <- NetworkChange(Y, R=2, m=2, mcmc=G, burnin=G, verbose=G, Waic=TRUE)
##D    out3 <- NetworkChange(Y, R=2, m=3, mcmc=G, burnin=G, verbose=G, Waic=TRUE)
##D 
##D    ## The most probable model given break number 0 to 3 and data is out1 according to WAIC 
##D    out <- BreakPointLoss(out0, out1, out2, out3, waic=TRUE)
##D 
##D    print(out[["ave.loss"]])
## End(Not run)





cleanEx()
nameEx("NetworkChange")
### * NetworkChange

flush(stderr()); flush(stdout())

### Name: NetworkChange
### Title: Changepoint analysis of a degree-corrected multilinear tensor
###   model
### Aliases: NetworkChange

### ** Examples


   ## Not run: 
##D    set.seed(1973)
##D    \## Generate an array (30 by 30 by 40) with block transitions
##D    from 2 blocks to 3 blocks
##D    Y <- MakeBlockNetworkChange(n=10, T=40, type ="split")
##D    G <- 100 ## Small mcmc scans to save time
##D 
##D    \## Fit multiple models for break number detection using Bayesian model comparison
##D    out0 <- NetworkStatic(Y, R=2, mcmc=G, burnin=G, verbose=G, Waic=TRUE)
##D    out1 <- NetworkChange(Y, R=2, m=1, mcmc=G, burnin=G, verbose=G, Waic=TRUE)
##D    out2 <- NetworkChange(Y, R=2, m=2, mcmc=G, burnin=G, verbose=G, Waic=TRUE)
##D    out3 <- NetworkChange(Y, R=2, m=3, mcmc=G, burnin=G, verbose=G, Waic=TRUE)
##D    outlist <- list(out0, out1, out2, out3)
##D 
##D    \## The most probable model given break number 0 to 3 and data is out1 according to WAIC 
##D    WaicCompare(outlist)
##D 
##D    plotU(out1)
##D  
##D    plotV(out1)
##D    
## End(Not run)



cleanEx()
nameEx("NetworkChangeRobust")
### * NetworkChangeRobust

flush(stderr()); flush(stdout())

### Name: NetworkChangeRobust
### Title: Changepoint analysis of a degree-corrected multilinear tensor
###   model with t-distributed error
### Aliases: NetworkChangeRobust

### ** Examples


   ## Not run: 
##D    set.seed(1973)
##D    ## Generate an array (30 by 30 by 40) with block transitions
##D    from 2 blocks to 3 blocks
##D    Y <- MakeBlockNetworkChange(n=10, T=40, type ="split")
##D    G <- 100 ## only 100 mcmc scans to save time
##D    ## Fit models
##D    out1 <- NetworkChangeRobust(Y, R=2, m=1, mcmc=G, burnin=G, verbose=G)
##D    ## plot latent node positions
##D    plotU(out1)
##D    ## plot layer-specific network generation rules
##D    plotV(out1)
##D    
## End(Not run)



cleanEx()
nameEx("NetworkStatic")
### * NetworkStatic

flush(stderr()); flush(stdout())

### Name: NetworkStatic
### Title: Degree-corrected multilinear tensor model
### Aliases: NetworkStatic

### ** Examples


   ## Not run: 
##D    set.seed(1973)
##D 
##D    \## generate an array with three constant blocks
##D    Y <- MakeBlockNetworkChange(n=10, shape=10, T=10, type ="constant")
##D    G <- 100 ## Small mcmc scans to save time
##D    out0 <- NetworkStatic(Y, R=2, mcmc=G, burnin=G, verbose=G)
##D 
##D    \## recovered latent blocks
##D    Kmeans(out0, n.cluster=3, main="Recovered Blocks")
##D 
##D    \## contour plot of latent node positions
##D    plotContour(out0)
##D 
##D    \## plot latent node positions
##D    plotU(out0)
##D 
##D    \## plot layer-specific network connection rules
##D    plotV(out0)
##D    
## End(Not run)




cleanEx()
nameEx("drawPostAnalysis")
### * drawPostAnalysis

flush(stderr()); flush(stdout())

### Name: drawPostAnalysis
### Title: Plot of latent node cluster
### Aliases: drawPostAnalysis

### ** Examples


   ## Not run: 
##D    set.seed(1973)
##D    ## generate an array with two constant blocks
##D    data(MajorAlly)
##D    Y <- MajorAlly
##D    fit <- NetworkChange(newY, R=2, m=2, mcmc=G, initial.s = initial.s,
##D           burnin=G, verbose=0, v0=v0, v1=v1)
##D    drawPostAnalysis(fit, Y, n.cluster=c(4, 4, 3))
##D    
## End(Not run)



cleanEx()
nameEx("drawRegimeRaw")
### * drawRegimeRaw

flush(stderr()); flush(stdout())

### Name: drawRegimeRaw
### Title: Plot of network by hidden regime
### Aliases: drawRegimeRaw

### ** Examples


   ## Not run: 
##D    set.seed(1973)
##D    ## generate an array with two constant blocks
##D    data(MajorAlly)
##D    Y <- MajorAlly
##D    fit <- NetworkChange(newY, R=2, m=2, mcmc=G, initial.s = initial.s,
##D           burnin=G, verbose=0, v0=v0, v1=v1)
##D    drawRegimeRaw(fit, newY)
##D    
## End(Not run)



cleanEx()
nameEx("kmeansU")
### * kmeansU

flush(stderr()); flush(stdout())

### Name: kmeansU
### Title: K-mean clustering of latent node positions
### Aliases: kmeansU

### ** Examples


   ## Not run: 
##D set.seed(1973)
##D    ## generate an array with two constant blocks
##D    Y <- MakeBlockNetworkChange(n=10, shape=10, T=10, type ="constant")
##D    out0 <- NetworkStatic(Y, R=2, mcmc=10, burnin=10,
##D    verbose=10, UL.Normal = "Orthonormal")
##D    ## latent node positions
##D    kmeansU(out0)
##D    
## End(Not run)



cleanEx()
nameEx("plotContour")
### * plotContour

flush(stderr()); flush(stdout())

### Name: plotContour
### Title: Contour plot of latent node positions
### Aliases: plotContour

### ** Examples


   ## Not run: 
##D set.seed(1973)
##D    \## generate an array with two constant blocks
##D    Y <- MakeBlockNetworkChange(n=10, shape=10, T=40, type ="constant")
##D    out0 <- NetworkStatic(Y, R=2, mcmc=10, burnin=10,
##D    verbose=10, UL.Normal = "Orthonormal")
##D    \## contour plot of latent node positions
##D    plotContour(out0)
##D    
## End(Not run)



cleanEx()
nameEx("plotU")
### * plotU

flush(stderr()); flush(stdout())

### Name: plotU
### Title: Plot of latent node positions
### Aliases: plotU

### ** Examples


   ## Not run: 
##D    set.seed(1973)
##D    \## generate an array with two constant blocks
##D    Y <- MakeBlockNetworkChange(n=10, shape=10, T=40, type ="constant")
##D    out0 <- NetworkStatic(Y, R=2, mcmc=10, burnin=10,
##D    verbose=10, UL.Normal = "Orthonormal")
##D    \## latent node positions
##D    plotU(out0)
##D    
## End(Not run)



cleanEx()
nameEx("plotV")
### * plotV

flush(stderr()); flush(stdout())

### Name: plotV
### Title: Plot of layer-specific network generation rules.
### Aliases: plotV

### ** Examples


   ## Not run: 
##D set.seed(1973)
##D    \## generate an array with two constant blocks
##D    Y <- MakeBlockNetworkChange(n=10, shape=10, T=40, type ="constant")
##D    out0 <- NetworkStatic(Y, R=2, mcmc=10, burnin=10,
##D    verbose=10, UL.Normal = "Orthonormal")
##D    \## latent node positions
##D    plotV(out0)
##D    
## End(Not run)



cleanEx()
nameEx("plotnetarray")
### * plotnetarray

flush(stderr()); flush(stdout())

### Name: plotnetarray
### Title: Plot of network array data
### Aliases: plotnetarray

### ** Examples


   ## Not run: 
##D    set.seed(1973)
##D    ## generate an array with two constant blocks
##D    Y <- MakeBlockNetworkChange(n=10, shape=1, T=20, type ="split")
##D    plotnetarray(Y)
##D    
## End(Not run)



cleanEx()
nameEx("theme_networkchange")
### * theme_networkchange

flush(stderr()); flush(stdout())

### Name: theme_networkchange
### Title: NetworkChange ggplot2 Theme
### Aliases: theme_networkchange

### ** Examples

## Not run: 
##D library(ggplot2)
##D ggplot(mtcars, aes(mpg, wt)) + geom_point() + theme_networkchange()
## End(Not run)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
