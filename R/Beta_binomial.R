
#' @title Beta-binomial distribution of second-stage binary outcomes
#'
#' @description This function calculates the probability density for second-stage binary outcomes
#' @param N total number of patients (stage1 + stage2)
#' @param n1 total number of patients in stage1
#' @param x1 number of successes in stage1
#' @param beta_par two shape parameters c(alpha,beta) for beta distribution
#' @return a dataframe with two variables: 
#' \itemize{
#' \item x2: all possible outcomes (number of successes) for stage2
#' \item prediction: density for each outcome
#' }
#' @examples
#'\dontrun{
#' Beta_binomial(N=40,n1=23,x1=16,beta_par=c(1,1))
#' }
#'@importFrom VGAM dbetabinom.ab
#' @export


Beta_binomial<-function(N,n1,x1,beta_par){
  n2<-N-n1
  predictions<-VGAM::dbetabinom.ab(x=seq(0,n2), size=n2, shape1=x1+beta_par[1],shape2=n1-x1+beta_par[2])
  return(data.frame(x2=seq(0,n2),predictions=predictions))
}