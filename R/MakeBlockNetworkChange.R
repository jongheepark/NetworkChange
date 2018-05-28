#' Build a synthetic block-structured temporal data with breaks
#'
#' MakeBlockNetworkChange generates a block-structured temporal data with breaks.
#'
#'
#' @param n The number of nodes within a block. The total number of nodes is n*block.number.
#'
#' @param break.point The point of break. 0 indicates the beginning, 0.5 indicates the middle,
#' and 1 indicates the end. 
#' @param break.point1 The point of the first break in "merge-split" or "split-merge". Any number between 0 and 0.5 can be chosen.
#' For example, 0 indicates #' the beginning, 0.25 indicates the 1/4th point, and 0.5 indicates the half point. 
#' @param break.point2 The point of the second breakin "merge-split" or "split-merge". Any number between 0.5 and 1 can be chosen.
#' For example, 0.5 indicates the beginning, 0.75 indicates the 3/4th point, and 1 indicates the end point. 
#' @param base.prob The probability of link among non-block members.
#' @param block.prob The probability of link among within-block members.
#' @param shape The speed of breaks.
#' The larger shape is, the faster the transition is. shape > 0 and shape < 8. 
#' @param T The length of time.
#' @param type The type of network changes. Options are "constant", "merge", "split",
#' "merge-split", "split-merge." If "constant" is chosen, the number of breaks is zero.
#' If "merge" or "split" is chosen, the number of breaks is one.
#' If either "merge-split" or "split-merge" is chosen, the number of breaks is two.
#' 
#' @return output  An output of \code{MakeBlockNetworkChange} contains a symmetric block-structured temporal network data set with breaks.
#'
#' @export
#'
#'

MakeBlockNetworkChange <- function(n=10, break.point = 0.5, 
                                   base.prob=.05, block.prob=.5,
                                   shape=1, T= 40, break.point1 = 0.25, break.point2 = 0.75,
                                   type ="merge" 
                                   ){
    call <- match.call()
    block.number = 3
    N <- block.number*n
    cpn = block.number -1

    b0 <- matrix(1, n, n)
    x = seq(-20, 20, length=T)
    intercept <- as.numeric(quantile(x, probs=break.point))
    sv = (1/(1+exp(intercept*shape - shape*x))-1/(1+exp(intercept*shape - shape*min(x))))*n
        ## (1/(1+exp(intercept - shape*max(x)))-1/(1+exp(intercept - shape*min(x))))
    ## plot(sv)
    
    if(type == "split-merge" |type ==  "merge-split"){
        ## break.point1 <- .50
        end <- ceiling(T/2)
        x1 = seq(-20,20,length=end)
        intercept1 <- as.numeric(quantile(x1, probs=break.point1/0.5))
        sv1 = (1/(1+exp(intercept1*shape - shape*x1))-1/(1+exp(intercept1*shape - shape*min(x1))))*n
        ## (1/(1+exp(intercept1*shape - shape*max(x1)))-1/(1+exp(intercept1*shape - shape*min(x1))))*n
        sv.split <- round(sort(sv1))
        
        ## break.point2<- .50
        x2 = seq(-20,20, length=(T - end))
        intercept2 <- as.numeric(quantile(x2, probs=(break.point2 - 0.5)/0.5))
        sv2 = (1/(1+exp(intercept2*shape - shape*x2))-1/(1+exp(intercept2*shape - shape*min(x2))))*n
            ## (1/(1+exp(intercept2*shape - shape*max(x2)))-1/(1+exp(intercept2*shape - shape*min(x2))))*n
        sv.merge <- round(sort(sv2, decreasing=TRUE))
    }

    sv = round(sv)
    data.array <- array(NA, dim = c(N,N,T))
    if(type == "merge"){
        sv <- sort(sv, decreasing=TRUE)
        for (t in 1:T){
            ## heterohpilic block till t<tp; homophilic blcok from t=tp         
            b1=matrix(1, sv[t], sv[t])
            b2=matrix(1, 2*n-sv[t], 2*n-sv[t])
            
            ## initial block structure with three homophilic groups
            B=rbind(cbind(matrix(1,sv[t],sv[t]), matrix(0,sv[t],2*n-sv[t]),matrix(0,sv[t],n)),
                cbind(matrix(0,2*n-sv[t],sv[t]),b2,matrix(0,2*n-sv[t],n)),cbind(1-b0,1-b0,b0))  
            A=(matrix(runif(N*N),N) < (B*(block.prob + base.prob) + (1-B)*base.prob))*1
            A[ row(A) <= col(A) ] <- 0
            A=A+t(A)
            data.array[,,t] <- A
        }
    }else if(type == "split-merge"){
        ## split
        for (t in 1:end){
            ## heterohpilic block till t<tp; homophilic blcok from t=tp         
            b1=matrix(1,sv.split[t],sv.split[t])
            b2=matrix(1,2*n-sv.split[t],2*n-sv.split[t]) 
            ## initial block structure with two homophilic groups
            B = rbind(
                cbind(matrix(1, sv.split[t], sv.split[t]),
                      matrix(0, sv.split[t], 2*n-sv.split[t]),
                      matrix(0, sv.split[t], n) ),
                cbind(matrix(0, 2*n-sv.split[t], sv.split[t]), b2, matrix(0,2*n-sv.split[t],n)),
                cbind(1-b0,1-b0,b0))  
            A=(matrix(runif(N*N),N)<(B*(block.prob + base.prob)+(1-B)*base.prob))*1
            A[ row(A) <= col(A) ] <- 0
            A=A+t(A)
            data.array[,,t] <- A           
        }
        ## merge
        for (s in (end+1):T){
            t <- s - end 
            ## heterohpilic block till t<tp; homophilic blcok from t=tp         
            b1=matrix(1,sv.merge[t],sv.merge[t])
            b2=matrix(1,2*n-sv.merge[t],2*n-sv.merge[t]) 
            ## initial block structure with two homophilic groups
            B=rbind(cbind(matrix(1,sv.merge[t],sv.merge[t]),matrix(0,sv.merge[t],2*n-sv.merge[t]),matrix(0,sv.merge[t],n)),
                cbind(matrix(0,2*n-sv.merge[t],sv.merge[t]),b2,matrix(0,2*n-sv.merge[t],n)),cbind(1-b0,1-b0,b0))  
            A=(matrix(runif(N*N),N)<(B*(block.prob + base.prob)+(1-B)*base.prob))*1
            A[ row(A) <= col(A) ] <- 0
            A=A+t(A)
            data.array[,,s] <- A           
        }
    }else if(type == "merge-split"){
       for (t in 1:end){
        ## merge
            ## heterohpilic block till t<tp; homophilic blcok from t=tp         
            b1=matrix(1,sv.merge[t],sv.merge[t])
            b2=matrix(1,2*n-sv.merge[t],2*n-sv.merge[t]) 
            ## initial block structure with two homophilic groups
            B=rbind(cbind(matrix(1,sv.merge[t],sv.merge[t]),matrix(0,sv.merge[t],2*n-sv.merge[t]),matrix(0,sv.merge[t],n)),
                cbind(matrix(0,2*n-sv.merge[t],sv.merge[t]),b2,matrix(0,2*n-sv.merge[t],n)),cbind(1-b0,1-b0,b0))  
            A=(matrix(runif(N*N),N)<(B*(block.prob + base.prob)+(1-B)*base.prob))*1
            A[ row(A) <= col(A) ] <- 0
            A=A+t(A)
            data.array[,,t] <- A           
        }

       ## split
        for (s in (end+1):T){
            t <- s - end 
            ## heterohpilic block till t<tp; homophilic blcok from t=tp         
            b1=matrix(1,sv.split[t],sv.split[t])
            b2=matrix(1,2*n-sv.split[t],2*n-sv.split[t]) 
            ## initial block structure with two homophilic groups
            B = rbind(
                cbind(matrix(1, sv.split[t], sv.split[t]),
                      matrix(0, sv.split[t], 2*n-sv.split[t]),
                      matrix(0, sv.split[t], n) ),
                cbind(matrix(0, 2*n-sv.split[t], sv.split[t]), b2, matrix(0,2*n-sv.split[t],n)),
                cbind(1-b0,1-b0,b0))  
            A=(matrix(runif(N*N),N)<(B*(block.prob + base.prob)+(1-B)*base.prob))*1
            A[ row(A) <= col(A) ] <- 0
            A=A+t(A)
            data.array[,,s] <- A           
        }
           
    }else if(type == "split"){
        sv <- sort(sv)
        for (t in 1:T){
            ## heterohpilic block till t<tp; homophilic blcok from t=tp         
            b1=matrix(1,sv[t],sv[t])
            b2=matrix(1,2*n-sv[t],2*n-sv[t]) 
            ## initial block structure with two homophilic groups
            B=rbind(cbind(matrix(1,sv[t],sv[t]),matrix(0,sv[t],2*n-sv[t]),matrix(0,sv[t],n)),
                cbind(matrix(0,2*n-sv[t],sv[t]),b2,matrix(0,2*n-sv[t],n)),cbind(1-b0,1-b0,b0))  
            A=(matrix(runif(N*N),N)<(B*(block.prob + base.prob)+(1-B)*base.prob))*1
            A[ row(A) <= col(A) ] <- 0
            A=A+t(A)
            data.array[,,t] <- A
        }
  
    } else if(type == "constant"){
        sv <- rep(max(sv), length(sv))
        for (t in 1:T){
            ## heterohpilic block till t<tp; homophilic blcok from t=tp         
            b1=matrix(1,sv[t],sv[t])
            b2=matrix(1,2*n-sv[t],2*n-sv[t]) 
            ## initial block structure with two homophilic groups
            B=rbind(cbind(matrix(1,sv[t],sv[t]),matrix(0,sv[t],2*n-sv[t]),matrix(0,sv[t],n)),
                cbind(matrix(0,2*n-sv[t],sv[t]),b2,matrix(0,2*n-sv[t],n)),cbind(1-b0,1-b0,b0))  
            A=(matrix(runif(N*N),N)<(B*(block.prob + base.prob)+(1-B)*base.prob))*1
            A[ row(A) <= col(A) ] <- 0
            A=A+t(A)
            data.array[,,t] <- A
        }
        
    } else{
        stop("Nothing was called.")
    }
    if(type == "split-merge"){
        true.ps1 <- sv.split/(max(sv.split))
        true.s1 <- sort(ifelse(true.ps1 > .5, 1, 2))
        true.ps2 <- sv.merge/(max(sv.merge))
        true.s2 <- ifelse(true.ps2 > .5, 2, 3)
        true.s <- c(true.s1, true.s2)
        ## true.ps <- c(true.ps1, 2- true.ps2)/2 ##
        true.ps <- c(true.ps1, true.ps2)
    } else if(type ==  "merge-split"){
        true.ps1 <- sv.merge/(max(sv.merge))
        true.s1 <- ifelse(true.ps1 > .5, 1, 2)
        true.ps2 <- sv.split/(max(sv.split))
        true.s2 <- sort(ifelse(true.ps2 > .5, 2, 3))
        true.s <- c(true.s1, true.s2)
        ## true.ps <- c(1 - true.ps1, true.ps2 + 1)/2
        true.ps <- c(true.ps1, true.ps2)
    } else{
        true.ps <- sv/(max(sv))
        true.s <- ifelse(true.ps > .5, 1, 2)
    }
    out <- data.array
    attr(out, "call") <- call
    attr(out, "type") <- type
    attr(out, "true.ps") <- true.ps
    attr(out, "true.s") <- sort(true.s)
    return(out)
}
