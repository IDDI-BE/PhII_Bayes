

#' @title Posterior predictive probability of a successful trial
#'
#' @description Posterior predictive probability of a successful trial, given number of successes
#'     at the end of stage1
#' @param N total number of patients (stage1 + stage2)
#' @param n1 total number of patients in stage1
#' @param x1 number of successes in stage1
#' @param beta_par two shape parameters c(alpha,beta) for beta distribution
#' @param p Proportion corresponding with some hypothesis, e.g. null hypothesis 
#' @param Delta_T Threshold for succesful trial in terms of Bayesian posterior probability
#'     P(theta>p|x1+x2)>Delta_T
#' @return Probability of a succesful trial
#' @examples
#'\dontrun{
#' Pred_Prob(N=68,n1=38,x1=27,beta_par=c(1  ,1)  ,p=0.664,Delta_T=0.95)
#' Pred_Prob(N=40,n1=23,x1=16,beta_par=c(0.6,0.4),p=0.6  ,Delta_T=0.9 )
#' }
#'
#'@references Lee JJ, Liu DD.A predictive probability design for phase II cancer clinical trialsClinical Trials 2008; 5: 93–106
#'@importFrom VGAM dbetabinom.ab
#'@export

Pred_Prob<-function(N,n1,x1,beta_par,p,Delta_T){
  n2<-N-n1
  i<-seq(0,n2)
  PP_i<- (VGAM::dbetabinom.ab(x=i,size=n2,shape1=beta_par[1]+x1,shape2=beta_par[2]+n1-x1))*
    as.numeric((1-pbeta(p,beta_par[1]+x1+i,beta_par[2]+N-x1-i))>Delta_T);
  PP<-sum(PP_i)
  return(PP)
}

