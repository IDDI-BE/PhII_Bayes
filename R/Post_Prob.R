

#' @title Bayesian posterior probability P(Theta>cutoff|data;prior)
#'
#' @description This function calculates the Bayesian posterior probability P(Theta>cutoff|data;prior)
#'
#' @param n number of patients for which results are available
#' @param x number of successes
#' @param beta_par two shape parameters c(alpha,beta) for beta distribution
#' @param p Proportion corresponding with some hypothesis, e.g. null hypothesis
#' @examples
#'\dontrun{
#' Post_Prob(n=68,x=52,beta_par=c(1,1),p=0.66)
#' }
#'
#' @export

Post_Prob<-function(n,x,beta_par,p){
  1-pbeta(p,x+beta_par[1],n-x+beta_par[2])
}