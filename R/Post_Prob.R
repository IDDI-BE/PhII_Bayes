
#' @title Bayesian posterior probability, given 1 observed proportion
#'
#' @description This function calculates the Bayesian posterior probability 
#' P(Delta>Dcut|n,x) ~beta(beta_par[1]+x,beta_par[2]+n-x) with prior beta distribution beta_par
#' @param n number of patients
#' @param x number of successes
#' @param beta_par two shape parameters c(alpha,beta) for beta distribution
#' @param Dcut Proportion corresponding with some hypothesis
#' @examples
#'\dontrun{
#' Post_Prob.f(n=40,x=16,beta_par=c(0.6,0.4),Dcut=0.6) # Example Lee 2008, p97, Table 1 
#' (first line of Panel B: Y=0)
#' }
#'
#' @export

Post_Prob<-function(n,x,beta_par,Dcut){
  1-pbeta(Dcut,x+beta_par[1],n-x+beta_par[2])
}