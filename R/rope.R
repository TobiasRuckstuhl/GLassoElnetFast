#' Rope Function
#'
#' Calculates explicitly a sparse precision matrix Theta using a ridge type penalty.
#' @param S Sample covariance matrix: symmetric p x p matrix
#' @param lambda Penalization parameter: non-negative scalar
#' @param Target (Optional) Target matrix: A symmetric and positive definite p x p matrix. The estimate approaches the target matrix as the penalty parameter increases.
#' @return
#' \item{Theta}{Estimated precision matrix}
#' @references
#' M. O. Kuismin, J. T. Kemppainen & M. J. Sillanpaa (2017).
#' Precision Matrix Estimation With ROPE.
#' Journal of Computational and Graphical Statistics. https://doi.org/10.1080/10618600.2016.1278002
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 10
#' p <- 5
#' Y <- array(rnorm(n*p),dim=c(n,p)) # data-matrix
#' S <- cov(Y) # sample covariance matrix
#' Theta <- rope(S=S,lambda=0.5)
#' }
#' @export

rope = function(S,lambda,Target=NULL){

  p = ncol(S)
  if(is.null(Target)) Target = matrix(0,p,p)
  S = S - lambda*Target
  LM = eigen(S)
  L = LM$values
  M = LM$vectors
  Lam = 2/(L + sqrt(L^2 + 4*lambda))
  Lam = diag(sort(Lam,decreasing=T),p)
  hatTheta = M[,p:1]%*%Lam%*%t(M[,p:1])

  return(hatTheta)

}
