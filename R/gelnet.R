#' Gelnet Function
#'
#' Estimates a sparse precision matrix Theta and the corresponding covariance matrix W using an elastic net type penalty.
#' @param S Sample covariance matrix: symmetric p x p matrix
#' @param lambda  Penalization parameter: non-negative scalar, symmetric p x p matrix or p-length vector. In the latter case, the penalty matrix has ij-th element sqrt(lambda[i]*lambda[j])
#' @param alpha Tuning paremeter: scalar value between 0 and 1. 0 corresponds to ridge (L2) type penalty and 1 corresponds to lasso (L1) type penalty. Values inbetween lead to elastic net type penalties.
#' @param zero (Optional) Indices of entries of precision matrix to be constrained to be zero. The input should be a matrix with two columns, each row indicating the indices of elements to be constrained to be zero. The solution must be symmetric, so you need only specify one of (i,j) and (j,i). An entry in the zero matrix overrides any entry in the rho matrix for a given element.
#' @param Theta (Optional) Precision matrix warm start: symmetric positive definit p x p matrix used as a warm start. If Theta is too far from the correct solution, the algorithm may no converge.
#' @param W (Optional) Covariance matrix warm start: symmetric positive definit p x p matrix used as a warm start. If W is too far from the correct solution, the algorithm may no converge.
#' @param Target (Optional) Target matrix: non-negative diagonal p x p matrix.
#' @param outer.maxit Maximal number of iterations for outer loop. Default is set to 1000.
#' @param outer.thr Threshold for convergence of the outer loop. Default is set to 1e-5.  Iterations stop when average absolute parameter change is less than outer.thr*average(abs(offdiag(S)))
#' @param inner.maxit Maximal number of iterations for inner loop. Default is set to 1000.
#' @param inner.thr Threshold for inner loop. Default is equal to outer.thr.
#' @param penalize.diagonal Should the diagonal of Theta be penalized. Default is set to TRUE.
#' @param active Should an active set of variables in the inner loop be used to increase computational speed. Default is set to FALSE.
#' @return
#' A list with the following
#' \item{Theta}{Estimated precision matrix}
#' \item{W}{Estimated covariance matrix}
#' \item{niter}{Number of iterations of outer.loop used.}
#' \item{del}{Change in parameter value at convergence.}
#' \item{conv}{Returns if the algorithm converged before reaching the maximal number of outer iterations}
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
#' S <- cov(Y) # sample covariance matrix
#' object <- gelnet(S=S,lambda=0.5,alpha=0.5)
#' }
#' @export
#' @useDynLib GLassoElnetFast


gelnet <-function(S,lambda,alpha,zero=NULL,Theta=NULL,W=NULL,Target=NULL,outer.maxit=1000,outer.thr=1e-5,inner.maxit=1000,inner.thr=outer.thr,penalize.diagonal=TRUE,active=FALSE){

  BIG <- 10e9

  p <- nrow(S)
  #if (is.null(S)) { print("S is required as input"); return(); }else { p=nrow(S) }
  #check1<-(nrow(S)==ncol(S))&&(S==t(S));
  #if (check1==0) { print("S should be square and symmetric"); return(); }

  if(!is.matrix(lambda) & length(lambda) != 1 & length(lambda) != nrow(S))
  {stop("Wrong number of elements in lambda")}

  if(is.vector(lambda) && length(lambda) > 1){ lambda = matrix(sqrt(lambda))%*%sqrt(lambda)}
  if(length(lambda) == 1){lambda = matrix(lambda,ncol=p,nrow=p)}

  if(is.vector(Target) && length(Target) > 1){Target=diag(Target)}



  if(!is.null(zero)){
    if(!is.matrix(zero)){ zero=matrix(zero,nrow=TRUE)}
    for(k in 1:nrow(zero)){
      i=zero[k,1]
      j=zero[k,2]
      lambda[i,j]=BIG
      lambda[j,i]=BIG
    }
  }


  thr <- outer.thr; thr2 <- inner.thr
  lambda1 <- diag(lambda)*alpha; lambda2 <- diag(lambda)*(1-alpha)
  conv <- TRUE



  if(is.null(Target)) {Target <- matrix(0,p,p)}
  U.vec <- rep(0,p)
  if(penalize.diagonal==TRUE){
    if(is.null(Theta)){ Theta <- matrix(0,p,p)
      for(i in 1:p){
        if(S[i,i]+lambda1[i] > 0){
          if(Target[i,i] < 1/(S[i,i]+lambda1[i])){
            U.vec[i] <- lambda1[i]
            if(alpha==1){ Theta[i,i] <- 1/(S[i,i]+lambda1[i])}
            else{Theta[i,i] <- (-S[i,i]-lambda1[i]+lambda2[i]*Target[i,i] + sqrt((S[i,i]+lambda1[i]-lambda2[i]*Target[i,i])^2 +4*lambda2[i]) ) / (2*lambda2[i])}
          }
          else{
            U.vec[i] <- 0
            Theta[i,i] <- Target[i,i]
          }
        }
        else{
          U.vec[i] <- 0
          Theta[i,i] <- Target[i,i]
        }
      }
    }
    if(is.null(W)){ W <- S; diag(W) <- diag(S) + U.vec  + lambda2*diag(Theta-Target) }
  }
  else{
    if(is.null(Theta)){ Theta <- diag( 1/diag(S) ) }
    if(is.null(W)){ W <- S }else{ diag(W) <- diag(S) }
  }

  #p1 <- nrow(Theta); p11 <- ncol(Theta); p2 <- nrow(W); p22 <- ncol(W);
  #if (p1!=p11) { print("check dimensions of Theta"); return(); }
  #if (p2!=p22) { print("check dimensions of W"); return(); }

  #check2 <- (p1==p2)&&(p1==p);
  #if (check2==0) { print("check dimensions of Theta,W,S"); return(); }
  niter <- 0L
  dlz <- 0

  mode(p)="integer"
  mode(S)="double"
  mode(lambda)="double"
  mode(alpha)="double"
  mode(Theta)="double"
  mode(W)="double"
  mode(Target)="double"
  mode(outer.maxit)="integer"
  mode(thr)="double"
  mode(inner.maxit)="integer"
  mode(thr2)="double"
  mode(niter)="integer"
  mode(penalize.diagonal)="logical"
  mode(active)="logical"
  mode(dlz)="double"



  loop<-.Fortran("gelnet",p,S,lambda,alpha,Theta,W,Target,outer.maxit,thr,inner.maxit,thr2,niter,penalize.diagonal,active,dlz,PACKAGE="GLassoElnetFast")

  # Theta <- loop$TTh; W <- loop$Wm
  # del <- loop$dlz

  #if(anyNA(Theta) | anyNA(W)){ print("error"); return(); }

  if(loop[[12]] == outer.maxit){ conv <- FALSE }
  #return(list(Theta=Theta,W=W,niter=1,del=1,conv=1))
  return(list(Theta=loop[[5]],W=loop[[6]],niter=loop[[12]],del=loop[[15]],conv=conv))
}
