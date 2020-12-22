#' Target Generating Function
#'
#' Generates different types of targets. The following types can be used:
#' "Identity" The identity matrix
#' "vI" A scaled identity matrix, where v is 1/mean(diag(S))
#' "Regression" A regression based target matrix. Use parameter: safety.scaling, nfolds.small, nfolds.number to modify. Need Y for this type of target.
#' "Eigenvalue" An Eigenvalue based approach. Warning: In some cases the values might get to large. Use parameter: fraction, const to modify.
#' "MSC" Maximal Single Correlation approach.
#' @param Y Data matrix: n x p matrix, where p is the number of predictors and n the number of observations
#' @param S Sample covariance matrix: symmetric p x p matrix (zero centered data)
#' @param type Type of the output Target. Possible values: "Identity","vI","Regression","Eigenvalue","MSC". The details to each type are in the description.
#' @param cor Should correlation matrices be used instead of covariance matrices. Default is set to FALSE.
#' @param fraction Needed for Eigenvalue type. The cutoff for which Eigenvalues are treated to be 0. Default is set to 1e-4
#' @param const Needed for Eigenvalue type. Default is set to 1
#' @param safety.scaling Needed for Regression type. Default is set to 1
#' @param nfolds.small Needed for Regression type. If true 10 fold is used. Default is set to TRUE.
#' @param nfolds.number Needed for Regression type. Number of folds used.
#' @importFrom stats cor cov cov2cor
#' @importFrom glmnet cv.glmnet
#' @return Returns the target matrix corresponding to the input type.
#' @references
#' Carel F.W. Peeters,  Anders Ellern Bilgrau and Wessel N. van Wieringen (2019).
#' rags2ridges: Ridge Estimation of Precision Matrices from High-Dimensional Data.
#' CRAN. https://CRAN.R-project.org/package=rags2ridges
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 10
#' p <- 5
#' Y <- array(rnorm(n*p),dim=c(n,p)) # data-matrix
#' target <- target(Y=Y,type="Identity")
#' }
#' @export

target = function(Y=NULL,S=NULL,type,cor=FALSE,fraction=1e-4,const=1,safety.scaling=1,nfolds.small=TRUE,nfolds.number=NA){
  if(is.null(S)==FALSE){p=nrow(S)}
  else if(is.null(Y)==FALSE){ p=ncol(Y)}
  else{ stop("Input either Y or S") }

  if(cor==FALSE){f <- function(x){cov(x)}}else{ f <- function(x){cor(x)}}

  if(type=="Identity"){ Target=diag(p)}
  else if(type=="vI"){ if(is.null(S)){ S=f(Y)}; Target=1/mean(diag(S))*diag(p)}
  else if(type=="Regression"){ if(is.null(Y)){stop("Need Y for this type")}
    safety.scaling = safety.scaling; nfolds.small = nfolds.small; nfolds.number = nfolds.number
    if(nfolds.small==T&!is.na(nfolds.number)){warning("nfolds.number set to 10 as nfolds.small = T was given as input")}
    if(safety.scaling<0){return("safety.scaling should be non-negative")}

    tar <- rep(NA, ncol(Y))

    for(i in 1:ncol(Y)){
      if(nfolds.small==T){LASSOfit <- cv.glmnet(x=Y[,-i], y=Y[,i])}
      if(nfolds.small==F){LASSOfit <- cv.glmnet(x=Y[,-i], y=Y[,i], nfolds = nfolds.number)}
      tar[i] <- min(LASSOfit$cvm + safety.scaling * LASSOfit$cvsd)
    }

    Target=diag(1/tar)

  }
  else if(type=="Eigenvalue"){ if(is.null(S)){ S=f(Y)};
    Eigs   <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
    const  <- mean(1/(Eigs[Eigs >= Eigs[1]*fraction]))
    Target <- const * diag(p)
  }
  else if(type=="MSC"){  if(is.null(S)){ S=f(Y)};  Target=diag(1/(diag(S)*(1-apply(cov2cor(S), 2, function(x){sort(abs(x))[p-1]})^2)))}
  else{ stop("Please input a valid type for the Target") }

  return(Target=Target)
}
