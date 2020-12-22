#' Crossvalidation Function for Graphical Elastic Net
#'
#' Calculates with nfold crossvalidation for a given set of scalars in the lambda vector the optimal penalty parameter to a given data set Y.
#' This parameter is then applied to the sample covariance S of Y and the resulting Theta and W is returned.
#' If cor=TRUE the procedure is done with correlation matrices instead of covariance matrices.
#' @param nfold Number of crossvalidation folds: integer value
#' @param Y Data matrix: n x p matrix, where p is the number of predictors and n the number of observations
#' @param lambda Penalization parameter: p-length vector of values to be minimized over
#' @param alpha Tuning paremeter: scalar value between 0 and 1. 0 corresponds to ridge (L2) type penalty and 1 corresponds to lasso (L1) type penalty. Values inbetween lead to elastic net type penalties.
#' @param ind Index order: n-length permutation vector. The order of observations put into the crossvalidation sets.
#' @param type Type of Target matrix. Possible values: NULL,"Identity","vI","Regression","Eigenvalue","MSC". Default is set to NULL.
#' @param Target Fixed Target matrix: p x p matrix. Is ignored if type is not NULL. Default is set to NULL.
#' @param outer.maxit Maximal number of iterations for outer loop. Default is set to 1000.
#' @param outer.thr Threshold for convergence of the outer loop. Default is set to 1e-5.  Iterations stop when average absolute parameter change is less than outer.thr*average(abs(offdiag(S)))
#' @param inner.maxit Maximal number of iterations for inner loop. Default is set to 1000.
#' @param inner.thr Threshold for inner loop. Default is equal to outer.thr.
#' @param penalize.diagonal Should the diagonal of Theta be penalized. Default is set to TRUE.
#' @param cor Should correlation matrices be used instead of covariance matrices. Default is set to FALSE.
#' @param rope Should the rope algorithm be used. If TRUE, alpha is automatically set to 0. Default is set to FALSE.
#' @return A list with the following
#'   \item{optimal}{Optimal value from lambda vector decided by crossvalidation}
#'   \item{lambda}{Sorted lambda vector from the input}
#'   \item{CV}{Average crossvalidation values of the sorted lambda vector}
#'   \item{Theta}{Estimated precision matrix for the optimal lambda value}
#'   \item{W}{Estimated covariance matrix for the optimal lambda value. Omitted if rope=TRUE.}
#'   \item{niter}{Number of iterations of outer.loop used. Omitted if rope=TRUE.}
#'   \item{del}{Change in parameter value at convergence. Omitted if rope=TRUE.}
#'   \item{conv}{Returns if the algorithm converged before reaching the maximal number of outer iterations. Omitted if rope=TRUE.}
#' @references
#' Jerome Friedman, Trevor Hastie and Rob Tibshirani (2019).
#' Graphical Lasso: Estimation of Gaussian Graphical Models.
#' CRAN. http://www-stat.stanford.edu/~tibs/glasso
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 10
#' p <- 5
#' Y <- array(rnorm(n*p),dim=c(n,p)) # data-matrix
#' object <- crossvalidation(nfold=5,Y=Y,ind=sample(1:n),lambda=seq(0.1,0.4,length.out=20))
#' }
#' @importFrom stats var cor
#' @export
#' @useDynLib GLassoElnetFast


crossvalidation <-function(nfold=5L,Y,lambda,alpha,ind=NULL,type=NULL,Target=NULL,outer.maxit=1000,outer.thr=1e-5,inner.maxit=1000,inner.thr=outer.thr/10,penalize.diagonal=TRUE,cor=FALSE,rope=FALSE){

  if(is.null(Y)){ print("Y is required as input"); return(); }
  if(is.vector(lambda) == FALSE){ print("lambda should be a vector"); return(); }else{lambda <- sort(lambda)}

  l = function(Theta,Sigma){
    res <- determinant(Theta, logarithm=TRUE)$modulus - sum(diag(Sigma%*%Theta))
    return(res)
  }
  n = nrow(Y)
  p = ncol(Y)

  if(cor==TRUE){ f = function(x){return(cor(x))} }else{ f = function(x){return(var(x))} }
  if(is.null(type)){ if(is.null(Target)){Target <- matrix(0,p,p)}; fTarget <- function(Y){Target}  }else{ fTarget <- function(Y){target(Y,type=type,cor=cor)} }


  if(is.null(ind)){ ind=c(1:n) }

  indnr <- rep(0L, nfold+1)

  indnr[1]=0
  rest = (n - (n %% nfold))/nfold
  for(i in 1:(nfold-1)){
    indnr[i+1]=i*rest
  }
  indnr[nfold+1]=n
  Targets <- array(NA,c(p,p,nfold))
  for(i in 1:nfold){
    indtest <- ind[(indnr[i]+1):indnr[i+1]]
    Targets[,,i] <- fTarget(Y[-indtest,])
  }

  CV <- rep(0, length(lambda))
  for(j in 1:length(lambda)){
    for(i in 1:nfold){
      indtest <- ind[(indnr[i]+1):indnr[i+1]]
      Strain <- f(Y[-indtest,])
      if(rope==FALSE){
        Thetatrain <- gelnet(S=Strain,lambda=lambda[j],alpha=alpha,outer.thr=outer.thr,penalize.diagonal=penalize.diagonal,Target=Targets[,,i])$Theta
      }else{
        Thetatrain <- rope(S=Strain,lambda=lambda[j],Target=Targets[,,i])
      }
      CV[j]=CV[j]+l(Thetatrain, f(Y[indtest,]))
    }
  }
  k=which.max(CV)
  if( (k==1) || (k==length(lambda)) ){ warning("The optimal lambda value is a boundary value. The range of possible values should be increased") }
  optimal=lambda[k]
  if(rope==FALSE){
    object <- gelnet(S=f(Y),lambda=optimal,alpha=alpha,outer.thr=outer.thr,penalize.diagonal=penalize.diagonal,Target=fTarget(Y))
    return(list(optimal=optimal,lambda=lambda,CV=CV/nfold,Theta=object$Theta,W=object$W,niter=object$niter,del=object$del,conv=object$conv))
  }else{
    Theta <- rope(S=f(Y),lambda=optimal,Target=fTarget(Y))
    return(list(optimal=optimal,lambda=lambda,CV=CV/nfold,Theta=Theta))
  }
}
