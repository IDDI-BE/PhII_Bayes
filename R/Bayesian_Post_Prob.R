

#' @title Bayesian posterior probability P(Theta>cutoff|data;prior)
#'
#' @description This function calculates the Bayesian posterior probability P(Theta>cutoff|data;prior)
#'
#' @param n number of patients for which results are available
#' @param succ number of successes
#' @param beta_par two shape parameters c(alpha,beta) for beta distribution
#' @param cutoff Proportion corresponding with some hypothesis, e.g. null hypothesis
#' @examples
#'\dontrun{
#' Bayesian_post_prob(n=68,succ=52,beta_par=c(1,1),cutoff=0.66)
#' }
#'
#' @export
#' @importFrom stats pbeta pbinom rbinom
#' @importFrom utils write.csv

Bayesian_Post_Prob<-function(n,succ,beta_par,cutoff){
  1-pbeta(cutoff,succ+beta_par[1],n-succ+beta_par[2])
}