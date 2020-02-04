#################################################
## functions.R
#################################################


## Gram-Schmidt orthogonalization
GramSchmidt <- function(U){
    R <- dim(U)[2]
    N <- dim(U)[1]
    q <- list()
    A <- matrix(U[,1], N, 1)
    q[[1]] <- (1/sqrt(sum(A^2)))*A
    for(r in 2:R){
        B <- U[,r] - (q[[1]]%*%U[,r])%*%q[[1]]
        q[[r]] <- (1/sqrt(sum(B^2)))*B 
    }
    out <- matrix(unlist(q), N, R)
    return(out)
}
               

marglike <- function(obj){
    ## obj <- list(out1, out2, out3)
    k <- length(obj)
    out <- matrix(NA, k, 3)
    for(i in 1:k){
        out[i, 1] <- attr(obj[[i]], "loglike")
        out[i, 2] <- attr(obj[[i]], "logmarglike.upper")
        out[i, 3] <- attr(obj[[i]], "logmarglike")
    }
    colnames(out) <- c("loglike", "logmarglike.upper", "logmarglike")
    rownames(out) <- names(obj)
    return(out)
}

normalize <- function(x){
  if(length(which(is.na(x)))==0) {
    out <- (x-mean(x))/sd(x)
  }
  else{
    out <- (x-mean(x,na.rm=T))/sd(x,na.rm=T)
  }
  return(out)
}


## normalization of U and V
Unormal <- function(U1){
    if(is.list(U1)){
        R <- dim(U1[[1]])[2]
    }
    else{
        R <- dim(U1)[2]
    }
    if(R ==1){
        su <- sqrt(apply(U1^2, 2, sum))
        ## su <- apply(U1, 2, sum)
        U2 <- U1*(1/su)

    }
    else{
        ## Euclidean norm
        su <- sqrt(apply(U1^2, 2, sum))
        ## su <- apply(U1, 2, sum)
        U2 <- U1%*%diag(1/su)
    }
  return(U2)
}

## centering U
Ucenter <- function(U1){
    if(is.list(U1)){
        R <- dim(U1[[1]])[2]
    }
    else{
        R <- dim(U1)[2]
    }
    if(R ==1){
        U2 <- scale(U1, scale = FALSE)

    }
    else{
        ## centering
        U2 <- scale(U1, scale = FALSE)
    }
  return(U2)
}


Vnormal <- function(V1, U1){
    if(is.list(U1)){
        R <- dim(U1[[1]])[2]
    }
    else{
        R <- dim(U1)[2]
    }
    if(R ==1){
        su <- sqrt(apply(U1^2, 2, sum))
        V2 <- V1*su^2
        
    }
    else{
        su <- sqrt(apply(U1^2, 2, sum))
        V2 <- V1%*%diag(su^2)
    }
    return(V2)
}


switchg <-  function(s1){
  ## for one step ahead case only
  ## example: s1<-c(1,1,1,2,2,2,3,4,4,4,4,5,5)
  ## switchg(s1)
  s<-max(s1)
  out<-matrix(0,s,s)
  ## all P(i,i+1) s is 1
  for (i in 1:(s-1)){
    out[i,i+1]<-1}
  ## diagonal elements is (# of occurrence-1)
  diag(out)<-table(s1)-1
  return(out)
}


"rdirichlet.cp" <- function(n, alpha){
  ## "rdirichlet.cp" picks n random deviates from the Dirichlet function
  ## with shape parameters alpha
  ## Note that alpha can contain zero to deal with transition matrix rowwise.
  ## It returns a matrix.
  
  col <-  length(alpha)
  out <-  matrix(NA, n, col)
  out[,which(alpha==0)]<- 0
  a   <-  alpha[alpha!=0]
  l   <-  length(a);
  x   <-  matrix(rgamma(l*n,a), ncol=l, byrow=TRUE);
  sm  <-  x%*%rep(1,l);
  dir.out <-  x/as.vector(sm);
  ## combine zero and nonzero parts prob
  out[,which(alpha!=0)] <- dir.out
  return(out)
}

trans.mat.prior <- function(m, n, a=NULL, b=NULL){
  if (!is.null(a)|!is.null(b)){
    a <- a
    b <- b
  }
  else {
    expected.duration <- round(n/(m+1))
    b <- 0.1
    a <- b*expected.duration
  }
  
  ## make a transition matrix
  trans<-diag(1, m+1)
  ## put a as diagonal elements except the last row
  diag(trans)[1:m]<-rep(a, m)
  ## put b in trans[i, i+1]
  for (i in 1:m){trans[i, i+1]<-b}
  return(trans)
}

###########################
## a bunch of plot functions
###########################

## for arrange_ggplot2.R	 
vp.layout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)


addTrans <- function(color,trans){
    ## thanks to Sacha Epskamp
    ## originally posted in
    ## http://stackoverflow.com/questions/12995683/any-way-to-make-plot-points-in-scatterplot-more-transparent-in-r
    ## This function adds transparancy to a color.
    ## Define transparancy with an integer between 0 and 255
    ## 0 being fully transparant and 255 being fully visable
    ## Works with either color and trans a vector of equal length,
    ## or one of the two of length 1.
    
    if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
    if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
    if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
    
    num2hex <- function(x)
        {
            hex <- unlist(strsplit("0123456789ABCDEF",split=""))
            return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
        }
    rgb <- rbind(col2rgb(color),trans)
    res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
    return(res)
}



## from http://stackoverflow.com/questions/3932038/plot-a-legend-outside-of-the-plotting-area-in-base-graphics
add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}


## code by Gelman and Vehtari (2014)
colVars <- function(a) {
    n <- dim(a)[[1]]; c <- dim(a)[[2]];
    return(.colMeans(((a - matrix(.colMeans(a, n, c), nrow = n, ncol =
                                      c, byrow = TRUE)) ^ 2), n, c) * n / (n - 1))}

waic <- function(log_lik){
    ## log_lik <- extract (stanfit, "log_lik")$log_lik
    dim(log_lik) <- if (length(dim(log_lik))==1) c(length(log_lik),1) else
    c(dim(log_lik)[1], prod(dim(log_lik)[2:length(dim(log_lik))]))
    S <- nrow(log_lik)
    n <- ncol(log_lik)
    lpd <- log(colMeans(exp(log_lik)))
    p_waic <- colVars(log_lik)
    p_waic1 <- 2*(log(colMeans(exp(log_lik))) - colMeans(log_lik))
    elpd_waic <- lpd - p_waic
    waic <- -2*elpd_waic
    waic2 <- -2*(lpd - p_waic1)
    loo_weights_raw <- 1/exp(log_lik-max(log_lik))
    loo_weights_normalized <- loo_weights_raw/
        matrix(colMeans(loo_weights_raw),nrow=S,ncol=n,byrow=TRUE)
    loo_weights_regularized <- pmin (loo_weights_normalized, sqrt(S))
    elpd_loo <- log(colMeans(exp(log_lik)*loo_weights_regularized)/
                        colMeans(loo_weights_regularized))
    p_loo <- lpd - elpd_loo
    pointwise <- cbind(waic,waic2, lpd,p_waic,p_waic1, elpd_waic,p_loo,elpd_loo)
    total <- colSums(pointwise)
    ## this is strange? SE = s/sqrt(n)
    se <- sqrt(n*colVars(pointwise))
    ## se <- sqrt(colVars(pointwise)/n)
    waic.ci <- c(total[1] - 1.96*se[1], total[1] + 1.96*se[1])
    return(list(waic=total["waic"], waic2=total["waic2"],  elpd_waic=total["elpd_waic"],
                p_waic=total["p_waic"], p_waic1 = total["p_waic1"], elpd_loo=total["elpd_loo"], p_loo=total["p_loo"],
                pointwise=pointwise, total=total, se=se, waic.ci=waic.ci))
}
#

############################################
## The following codes are from eigenmodel
## https://cran.r-project.org/web/packages/eigenmodel/
############################################
rsmn <- function (m, n)
{
  matrix(rnorm(m * n), m, n)
}

##### Create array out of factors
M.U <- function(U)
{
  r <- dim(U[[1]])[2]
  m <- length(U)
  n <- sapply(U,length)/r
  
  M <- array(0,dim=n)
  
  for(k in 1:r)
  {
    tmp <- U[[1]][,k]
    for(j in (2:m))
    {
      tmp <- outer(tmp,U[[j]][,k])
    }
    M <- M+tmp
  }
  M
}
#####
rMVNorm <- function(n,mu,Sigma)
{
  E <- matrix(rnorm(n*length(mu)),n,length(mu))
  t(  t(E%*%chol(Sigma)) +c(mu))
}
##

##
rsmn <- function (m, n)
{
  matrix(rnorm(m * n), m, n)
}
rmn <- function (M = 0, Srow, Scol)
{
  m  <-  dim(Srow)[1]
  n  <-  dim(Scol)[1]
  tmp  <-  eigen(Srow)
  Srow.h  <-  tmp$vec %*% diag(sqrt(tmp$val),nrow=m) %*% t(tmp$vec)
  tmp  <-  eigen(Scol)
  Scol.h  <-  tmp$vec %*% diag(sqrt(tmp$val),nrow=n) %*% t(tmp$vec)
  Z  <-  rsmn(m, n)
  Srow.h %*% Z %*% Scol.h + M
}
rwish <- function(S0, nu)
{
  # sample from a Wishart distribution
  sS0 <- chol(S0)
  Z <- matrix(rnorm(nu*dim(S0)[1]),nu, dim(S0)[1])%*%sS0
  t(Z)%*%Z
  # expectation is S0*nu
  # S0 appears as S0^{-1} in the exponent of the Wishart density
}

## modified UDS.als
UV.lsq <-function(Y,R, U, V, tol=1e-5, iter.limit = 50)
{
  m <- dim(Y)[1] ; n <- dim(Y)[3]
  M0 <- M.U(list(U,U,V))
  rdif <- 1  

  Y0 <- Y
  for(k in 1:n) {
      Y0[,,k][!upper.tri(Y0[,,k])]<- 0
  }
  Y0 <- aperm(Y0,c(3,1,2))
  
  iter <- 1
  while(rdif>tol & iter < iter.limit) {
      
      for(i in sample(m)){
          Ui<-U ; Ui[i,]<-0
          VU<- aperm(array(apply(Ui,1,"*",t(V)),dim=c(R,n,m)),c(3,2,1))
          zi<-Y[i,,]
          L<- apply(VU*array(rep(zi,R),dim=c(m,n,R)),3,sum)
          Q<- (t(Ui)%*%Ui ) * ( t(V)%*%V )
          U[i,]<-solve(Q)%*%L
      }
      
      if(R == 1){
          Q<-((t(U)%*%U)^2-
                         matrix(sum(apply(U^2,1,function(x){x%*%t(x)})),R,R))/2
      }else{
          Q<-((t(U)%*%U)^2-
                         matrix(apply(apply(U^2,1,function(x){x%*%t(x)}),1,sum),R,R))/2
      }
      
      UU<-aperm(array( apply(U,1,"*",t(U)) ,dim=c(R,m,m) ),c(2,3,1))
      ZUU<-array(apply(UU,3,function(x){apply(Y0,1,"*",x)}),
                 dim=c(m,m,n,R))
      L<-apply(ZUU,c(3,4),sum)
      V<-L%*%solve(Q)
      
      M1<-M.U(list(U,U,V))  
      rdif<- mean( (M1-M0)^2 )/mean(M0^2) 
      M0<-M1
      if(iter%%10 == 0){
          cat("Initializing using the Alternating Least Squares Method: # of iteration = ", iter, "\n")
      }
      iter <- iter + 1
  }
  list(U=U,V=V, iter=iter)
}

plotU.separate <- function(OUT, Year=NULL, names=NULL, main=""){
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
        p.list <- NA
        x <- OUT
        x.mean <-  apply(x, 1:2, mean)
        tmp <- eigen(x.mean)
        U <- tmp$vec[, order(-abs(tmp$val))[seq(1, R, length = R)], 
                     drop = FALSE] * sqrt(K[1])
        df <- data.frame(First = U[,1], Second = U[,2],
                         Size = sqrt((U[,1])^2 + (U[,2])^2),
                         Names = names)
        title <- paste0("Latent Space of No Break Model")
        p.list <- ggplot(df, aes(x=First, y = Second, label=Names)) + geom_point(size = df$Size+1, colour = alpha("red", 1/5)) +
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
        ## if(ns < 5){
        ##      multiplot(plotlist = p.list, cols=ns)
        ##   }else{
        ##       multiplot(plotlist = p.list, cols=ceiling(ns/2))
        ##       }
    }
    return(p.list)
}

