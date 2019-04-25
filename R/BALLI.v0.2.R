###########################################
################## BALLI ##################
###########################################
########## Made by Kyungtaek Park #########
######### Created on 26 Dec 2017 ##########
######## final edited on 24 Apr 2019 #######
###########################################

### Set Class ###
#	Created by Kyungtaek Park on 26 Dec 2017
# Last modified 4 Apr 2019

#' Class TecVarList
#' Class \code{TecVarList} holds technical variance
#' @name TecVarList-class
#' @rdname TecVarList-class
#' @exportClass TecVarList

setClass("TecVarList", representation("list"))

#' Class Balli
#' Class \code{Balli} holds results from BALLI
#' @name Balli-class
#' @rdname Balli-class
#' @exportClass Balli

setClass("Balli", representation("list"))

#' Class LargeDataObject
#' Class \code{LargeDataObject} holds large data such as technical variance and results from BALLI fit
#' @name LargeDataObject-class
#' @rdname LargeDataObject-class
#' @exportClass LargeDataObject

setClass("LargeDataObject")

setIs("TecVarList","LargeDataObject")
setIs("Balli","LargeDataObject")

#  Print and show method large data objects derived from code written in limma package
setMethod("show","LargeDataObject",
          function(object)
          {
            cat("An object of class \"",class(object),"\"\n",sep="")
            for (what in names(object)) {
              x <- object[[what]]
              cat("$",what,"\n",sep="")
              printHead(x)
              cat("\n")
            }
            for (what in setdiff(slotNames(object),".Data")) {
              x <- slot(object,what)
              if(length(x) > 0) {
                cat("@",what,"\n",sep="")
                printHead(x)
                cat("\n")
              }
            }
          })

#' Technical Variance Estimation
#' @description Estimate technical variance by using voom-trend. The code is derived from voom function in limma package
#' @param counts a DGEList object
#' @param design design matrix with samples in row and coefficient(s) to be estimated in column
#' @param lib.size numeric vector containing total library sizes for each sample
#' @param span width of the lowess smoothing window as a proportion
#' @param ...	other arguments are passed to lmFit.
#' @return an TecVarList object with the following components:
#' \item{targets}{matrix containing covariables, library sizes and normalization foctors of each sample}
#' \item{design}{design matrix with samples in row and covariable(s) to be estimated in column}
#' \item{logcpm}{logcpm values of each gene and each sample}
#' \item{tecVar}{estimated techical variance of each gene and each sample}
#' @examples
#' expr <- data.frame(t(sapply(1:1000,function(x)rnbinom(20,mu=500,size=50))))
#' group <- c(rep("A",10),rep("B",10))
#' design <- model.matrix(~group, data = expr)
#' dge <- DGEList(counts=expr, group=group)
#' dge <- calcNormFactors(dge)
#' tecVarEstim(dge,design)
#' @export

# Derived from code written in limma package
#	Created by Kyungtaek Park on 26 Dec 2017
# Last modified 4 Apr 2019

tecVarEstim <- function(counts,design=NULL,lib.size=NULL,span=0.5,...)
{
  out <- list()

  # Check counts
  if(is(counts,"DGEList")) {
    logcpm <- cpm(counts,log=T)
    out$targets <- counts$samples
    if(is.null(design) && diff(range(as.numeric(counts$sample$group)))>0) design <- model.matrix(~group,data=counts$samples)
    if(is.null(lib.size)) lib.size <- with(counts$samples,lib.size*norm.factors)
    counts <- counts$counts
  } else {
    stop("counts must be DEGList")
  }

  n <- nrow(counts)
  if(n < 2L) stop("Need at least two genes to fit a mean-variance trend")

  # Fit linear model to log2-counts-per-million
  y <- t(log2(t(counts+0.5)/(lib.size+1)*1e6))
  fit <- lmFit(y,design,...)

  # Fit lowess trend to sqrt-standard-deviations by log-count-size
  sx <- fit$Amean+mean(log2(lib.size+1))-log2(1e6)
  sy <- sqrt(fit$sigma)

  allzero <- rowSums(counts)==0
  if(any(allzero)) {
    sx <- sx[!allzero]
    sy <- sy[!allzero]
  }

  l <- lowess(sx,sy,f=span)
  f <- approxfun(l, rule=2)

  # Find individual quarter-root fitted counts
  if(fit$rank < ncol(design)) {
    j <- fit$pivot[1:fit$rank]
    fitted.values <- fit$coef[,j,drop=FALSE] %*% t(fit$design[,j,drop=FALSE])
  } else {
    fitted.values <- fit$coef %*% t(fit$design)
  }
  fitted.cpm <- 2^fitted.values
  fitted.count <- 1e-6 * t(t(fitted.cpm)*(lib.size+1))
  fitted.logcount <- log2(fitted.count)

  # Apply trend to individual observations
  ssquare <- f(fitted.logcount)^4
  dim(ssquare) <- dim(fitted.logcount)

  inverse.mu.hat <- 1/(2^fitted.logcount*(1+1/2*log(2)^2*ssquare))
  colnames(inverse.mu.hat) <- colnames(y)

  out$targets$lib.size <- lib.size
  out$design <- design
  out$logcpm <- logcpm
  out$tecVar <- inverse.mu.hat

  new("TecVarList",out)

}

#' balliFit
#' @description Estimates likelihood and Bartlett correction factor using BALLI algorithm of each gene
#' @param y_mat numeric vector containing log-cpm values of each gene and each sample
#' @param x_mat design matrix with samples in row and covariable(s) to be estimated in column
#' @param tecVar numeric vector containing estimated technical variance of a gene of each sample
#' @param intVar numeric vector designating interest variable(s) which is(are) column number(s) of x_mat
#' @param full logical value designating full model (TRUE) or reduced model (FALSE).
#' @param cfault initial value of index showing whether converged (0) or not (1).
#' @param miter maximum number of iteration to converge.
#' @param conv threshold for convergence
#' @return following components are estimated
#' \item{ll}{log-likelihoods}
#' \item{beta}{coefficients of interested variable(s)}
#' \item{alpha}{coefficients of nuisance variable(s)}
#' \item{BCF}{Bartlett's correction factor}
#' \item{cfault}{index whether converged or not}
#' @examples
#' expr <- data.frame(t(sapply(1:1000,function(x)rnbinom(20,mu=500,size=50))))
#' group <- c(rep("A",10),rep("B",10))
#' design <- model.matrix(~group, data = expr)
#' dge <- DGEList(counts=expr, group=group)
#' dge <- calcNormFactors(dge)
#' tV <- tecVarEstim(dge,design)
#' gtv <- tV$tecVar[1,]
#' gdat <- data.frame(logcpm=tV$logcpm[1,],design,tecVar=gtv)
#' gy <- matrix(unlist(gdat[,1]),ncol=1)
#' gx <- matrix(unlist(gdat[,2:(ncol(gdat)-1)]),ncol=ncol(gdat)-2)
#' balliFit(y_mat=gy,x_mat=gx,tecVar=gtv,intVar=2,full=TRUE,cfault=0,miter=200,conv=1e-6)
#' @export

#	Created by Kyungtaek Park on 26 Dec 2017
# Last modified 4 Apr 2019

balliFit <- function(y_mat,x_mat,tecVar,intVar=2,full=T,cfault=0,miter=200,conv=1e-6) {

  TT <- length(y_mat)
  xp_tilda <- x_mat[,intVar,drop=F]
  xnp_tilda <- x_mat[,-intVar,drop=F]

  # Initial Value
  sigma2_m1 <- var(y_mat)
  varCov <- diag(c(tecVar)+c(sigma2_m1))
  vi <- ginv(varCov)
  Anpti <- solve(t(xnp_tilda) %*% vi %*% xnp_tilda)
  xpp_tilda <- (diag(TT)-xnp_tilda %*% Anpti %*% t(xnp_tilda) %*% vi) %*% xp_tilda
  Appti <- solve(t(xpp_tilda) %*% vi %*% xpp_tilda)

  beta <- solve(t(x_mat) %*% vi %*% x_mat) %*% t(x_mat) %*% vi %*% y_mat
  psi <- beta[intVar,]
  xi <- Anpti %*% t(xnp_tilda) %*% vi %*% (y_mat - xpp_tilda %*% psi)
  z <- y_mat - xpp_tilda %*% psi - xnp_tilda %*% xi

  d1vs <- diag(TT)
  d1vis <- -vi %*% d1vs %*% vi
  d1xppts <- -xnp_tilda %*% Anpti %*% t(xnp_tilda) %*% d1vis %*% xpp_tilda
  d1ls <- -1/2*sum(diag(vi %*% d1vs))-1/2*t(z) %*% d1vis %*% z + t(psi) %*% t(d1xppts) %*% vi %*% z

  # Fisher Scoring Method
  ii <- 1
  diffs <- 1

  while(diffs > conv) {

    if(full) {
      psi <- Appti %*% t(xpp_tilda) %*% vi %*% (y_mat - xnp_tilda %*% xi)
    } else {
      psi <- matrix(rep(0,length(intVar)),ncol=1)
    }

    xi <- Anpti %*% t(xnp_tilda) %*% vi %*% (y_mat - xpp_tilda %*% psi)
    z <- y_mat - xpp_tilda %*% psi - xnp_tilda %*% xi

    sigma2_i_mat <- - 1/2 * sum(diag(d1vis %*% d1vs)) + t(psi) %*% t(d1xppts) %*% vi %*% d1xppts %*% psi
    sigma2_m2 <- sigma2_m1 + 1/sigma2_i_mat * d1ls

    varCov <- diag(c(tecVar)+c(sigma2_m2))
    vi <- ginv(varCov)
    Anpti <- solve(t(xnp_tilda) %*% vi %*% xnp_tilda)
    xpp_tilda <- (diag(TT)-xnp_tilda %*% Anpti %*% t(xnp_tilda) %*% vi) %*% xp_tilda
    Appti <- solve(t(xpp_tilda) %*% vi %*% xpp_tilda)
    d1vis <- -vi %*% d1vs %*% vi
    d1xppts <- -xnp_tilda %*% Anpti %*% t(xnp_tilda) %*% d1vis %*% xpp_tilda
    d1ls <- -1/2*sum(diag(vi %*% d1vs))-1/2*t(z) %*% d1vis %*% z + t(psi) %*% t(d1xppts) %*% vi %*% z

    diffs <- abs(sigma2_m2-sigma2_m1)

    if( ii == miter ) {

      iii <- 1
      diffs <- 1
      sigma2_m2 <- sigma2_m1

      while(diffs > conv) {

        if(full) {
          psi <- Appti %*% t(xpp_tilda) %*% vi %*% (y_mat - xnp_tilda %*% xi)
        } else {
          psi <- matrix(rep(0,length(intVar)),ncol=1)
        }

        xi <- Anpti %*% t(xnp_tilda) %*% vi %*% (y_mat - xpp_tilda %*% psi)
        z <- y_mat - xpp_tilda %*% psi - xnp_tilda %*% xi

        sigma.ll <- function(ss) {1/2*sum(log(tecVar+ss))+1/2*t(z) %*% solve(diag(c(tecVar+ss))) %*% z}
        # An integer code. 0 indicates successful completion (which is always the case for "SANN" and "Brent").
        sigma2_m2 <- optim(sigma2_m1,sigma.ll,method="Brent",lower=0,upper=1e+10)$par

        varCov <- diag(c(tecVar)+c(sigma2_m2))
        vi <- ginv(varCov)
        Anpti <- solve(t(xnp_tilda) %*% vi %*% xnp_tilda)
        xpp_tilda <- (diag(TT)-xnp_tilda %*% Anpti %*% t(xnp_tilda) %*% vi) %*% xp_tilda
        Appti <- solve(t(xpp_tilda) %*% vi %*% xpp_tilda)
        d1vis <- -vi %*% d1vs %*% vi
        d1xppts <- -xnp_tilda %*% Anpti %*% t(xnp_tilda) %*% d1vis %*% xpp_tilda

        diffs <- abs(sigma2_m2-sigma2_m1)

        if (sigma2_m2 < 1e-8) {
          sigma2_m2 <- 0
          break
        }

        if( iii == miter ) {
          cfault <- 1
          break
        } else {
          iii <- iii + 1
          sigma2_m1 <- sigma2_m2
        }

      }

      break

    } else {
      ii <- ii + 1
      sigma2_m1 <- sigma2_m2
    }
  }

  if(sigma2_m2 <= 0) {

    sigma2_m2 <- 0
    varCov <- diag(c(tecVar))
    vi <- ginv(varCov)
    Anpti <- solve(t(xnp_tilda) %*% vi %*% xnp_tilda)
    xpp_tilda <- (diag(TT)-xnp_tilda %*% Anpti %*% t(xnp_tilda) %*% vi) %*% xp_tilda
    Appti <- solve(t(xpp_tilda) %*% vi %*% xpp_tilda)
    d1vs <- matrix(0,nrow=nrow(vi),ncol=ncol(vi))
    d1vis <- d1vs
    d1xppts <- -xnp_tilda %*% Anpti %*% t(xnp_tilda) %*% d1vis %*% xpp_tilda

    psi <- beta[intVar,]
    xi <- Anpti %*% t(xnp_tilda) %*% vi %*% (y_mat - xpp_tilda %*% psi)

    if(full) {
      psi <- Appti %*% t(xpp_tilda) %*% vi %*% (y_mat - xnp_tilda %*% xi)
    } else {
      psi <- matrix(rep(0,length(intVar)),ncol=1)
    }

    xi <- Anpti %*% t(xnp_tilda) %*% vi %*% (y_mat - xpp_tilda %*% psi)
    z <- y_mat - xpp_tilda %*% psi - xnp_tilda %*% xi

  }

  ll <- -TT/2*log(2*pi) - 1/2*sum(log(diag(varCov)))-1/2*t(z) %*% vi %*% z

  if(full) {
    return(list(ll=ll,beta=psi,alpha=xi,cfault=cfault))
  } else {
    if(sigma2_m2 != 0) {
      d2viss <- -2 * d1vis %*% d1vs %*% vi

      D <- matrix(1/2*sum(diag(d1vis %*% d1vs)))
      M <- matrix(sum(diag(Appti %*% (t(xpp_tilda) %*% d2viss %*% xpp_tilda + 2 * t(d1xppts) %*% d1vis %*% xpp_tilda))))
      P <- matrix(sum(diag(t(xpp_tilda) %*% d1vis %*% xpp_tilda %*% Appti %*% t(xpp_tilda) %*% d1vis %*% xpp_tilda %*% Appti)))
      tau <- sum(diag(Appti %*% t(xpp_tilda) %*% d1vis %*% xpp_tilda))
      gamma <- 0
      nu <- sum(diag(Anpti %*% t(xnp_tilda) %*% d1vis %*% xnp_tilda))

      C <- sum(diag(solve(D) %*% (-1/2*M + 1/4*P - 1/2*(gamma+nu) %*% t(tau))))

    } else {
      C <- 0
    }
    return(list(ll=ll,BCF=C,cfault=cfault))
  }

}


#' BALLI
#' @description DEG analysis using BALLI algorithm
#' @param object a TecVarList object
#' @param intV numeric vector designating interest variable(s) which is(are) column number(s) of design matrix
#' @param logcpm logcpm values for each gene and each sample
#' @param tecVar estimated technical variance values for each gene and each sample
#' @param design design matrix with samples in row and covariable(s) to be estimated in column
#' @param numCores number of cores to be used for multithreding. If NULL, a single core is used
#' @param threshold threshold for convergence
#' @param maxiter maximum number of iteration to converge of estimated biological variance. If not, biological variance is estimated by using Brent method
#' @return an Balli object including Result and topGenes list. Following components are shown by Result (same order of genes with input data) and topGenes (ordered by pBALLI in Result) :
#' \item{log2FC}{log2 fold changes of interest variable(s)}
#' \item{lLLI}{log-likelihoods estimated by LLI}
#' \item{lBALLI}{log-likelihoods estimated by BALLI}
#' \item{pLLI}{p-values estimated by LLI}
#' \item{pBALLI}{p-values estimated by BALLI}
#' \item{BCF}{Bartlett's correction factor}
#' expr <- data.frame(t(sapply(1:1000,function(x)rnbinom(20,mu=500,size=50))))
#' group <- c(rep("A",10),rep("B",10))
#' design <- model.matrix(~group, data = expr)
#' dge <- DGEList(counts=expr, group=group)
#' dge <- calcNormFactors(dge)
#' tV <- tecVarEstim(dge,design)
#' balli(tV,intV=2)
#' @export

#	Created by Kyungtaek Park on 26 Dec 2017
# Last modified 4 Apr 2019

balli <- function(object,intV=2,logcpm=NULL,tecVar=NULL,design=NULL,numCores=NULL,threshold=1e-6,maxiter=200) {

  out <- list()

  if(is(object,"TecVarList")) {
    # Check logcpm
    if(is.null(logcpm)) logcpm <- object$logcpm
    # Check technical variance
    if(is.null(tecVar)) tecVar <- object$tecVar
    # Check design
    if(is.null(design)) design <- object$design
  } else {
    if (is.null(logcpm)) stop("logcpm must be designated or object must be TecVarList")
    if (is.null(tecVar)) stop("tecVar must be designated or object must be TecVarList")
    if (is.null(design)) stop("design must be designated or object must be TecVarList")
  }

  design <- as.matrix(design)

  balliFitAllGenes <- function(gidx){
    gtv <- tecVar[gidx,]
    gdat <- data.frame(logcpm=logcpm[gidx,],design,tecVar=gtv)
    gy <- matrix(unlist(gdat[,1]),ncol=1)
    gx <- matrix(unlist(gdat[,2:(ncol(gdat)-1)]),ncol=ncol(gdat)-2)

    LL_h1 <- balliFit(y_mat=gy,x_mat=gx,tecVar=gtv,intVar=intV,miter=maxiter,conv=threshold,full=T)
    if (LL_h1$cfault == 1) {
      message("number of iteration increases")
      stop("check")
      LL_h1 <- balliFit(y_mat=gy,x_mat=gx,tecVar=gtv,intVar=intV,conv=threshold,full=T,miter=500)
    }
    LL_h0 <- balliFit(y_mat=gy,x_mat=gx,tecVar=gtv,intVar=intV,miter=maxiter,conv=threshold,full=F)
    if (LL_h0$cfault == 1) {
      message("number of iteration increases")
      stop("check")
      LL_h0 <- balliFit(y_mat=gy,x_mat=gx,tecVar=gtv,intVar=intV,conv=threshold,full=F,miter=500)
    }

    LL_full <- LL_h1$ll
    beta <- LL_h1$beta
    LL_red <- LL_h0$ll
    BCF <- LL_h0$BCF

    LR <- -2*(LL_red-LL_full)
    BCLR <- LR/(1+BCF/length(intV))
    pLR <- pchisq(LR,df=length(intV),lower.tail=F)
    pBCLR <- pchisq(BCLR,df=length(intV),lower.tail=F)

    return(c(beta,LR,BCLR,pLR,pBCLR,BCF))
  }

  if(is.null(numCores)) {
    res <- t(sapply(1:nrow(logcpm),balliFitAllGenes))
  } else {
    res <- do.call(rbind,mclapply(1:nrow(logcpm),balliFitAllGenes,mc.cores=numCores))
  }

  res <- as.data.frame(res)
  rownames(res) <- rownames(logcpm)
  colnames(res) <- c(paste0("log2FC_",colnames(design)[intV]),"lLLI","lBALLI","pLLI","pBALLI","BCF")

  adjpLLI <- p.adjust(res$pLLI,method="fdr")
  adjpBALLI <- p.adjust(res$pBALLI,method="fdr")

  top_res <- data.frame(res[,c(paste0("log2FC_",colnames(design)[intV]),"pLLI","pBALLI")],adjpLLI=adjpLLI,adjpBALLI=adjpBALLI)
  top_res <- top_res[order(top_res$pBALLI),]

  out$Result <- res
  out$topGenes <- top_res

  new("Balli",out)

}

#'@importFrom edgeR DGEList calcNormFactors cpm
#'@importFrom limma lmFit printHead
#'@importFrom MASS ginv
#'@importFrom parallel mclapply
#'@importFrom stats approxfun lowess model.matrix optim p.adjust pchisq var
#'@importFrom methods is new slot slotNames
NULL




