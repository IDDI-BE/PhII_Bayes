

#' @title Posterior predictive probability, given 1 observed proportion
#'
#' @description Posterior predictive probability for posterior probability outcome
#' at the end of the trial, given currently observed number of patients and successes
#' Definition: Sum(P(X2=i|x1)I[P(Delta>Dcut|x1,X2=i)>PostProb]) for all X2=i,...,n2,
#' with P(Delta>Dcut|x1,X2=i) ~beta(beta_par[1]+x1+i,beta_par[2]+N-x1-i) with prior distribution beta_par
#' with X2~beta-binomial(n2,beta_par[1]+x1,beta_par[2]+n1-x1)
#' @param N total number of patients at the end of the study
#' @param n1 number of patients currently observed
#' @param x1 number of successes currently observed
#' @param beta_par two shape parameters c(alpha,beta) for prior beta distribution
#' @param Dcut Proportion corresponding with some hypothesis 
#' @param PostProb Threshold for outcome at the end of the study, in terms of Bayesian posterior probability
#'     P(theta>Dcut|x1+x2)>PostProb
#' @return Predictive Probability for outcome at the end of the study
#' @examples
#'\dontrun{
#' Pred_Prob(N=68,n1=38,x1=27,beta_par=c(1  ,1)  ,Dcut=0.664,PostProb=0.95)
#' Pred_Prob(N=40,n1=23,x1=16,beta_par=c(0.6,0.4),Dcut=0.6  ,PostProb=0.9 ) # Example Lee 2008, 
#' p97, Table 1 (PP in text)
#' }
#'
#'@references Lee JJ, Liu DD.A predictive probability design for phase II cancer clinical trialsClinical Trials 2008; 5: 93–106
#'@importFrom VGAM dbetabinom.ab
#'@export

Pred_Prob<-function(N,n1,x1,beta_par,Dcut,PostProb){
  n2<-N-n1
  i<-seq(0,n2)
  PP_i<- (VGAM::dbetabinom.ab(x=i,size=n2,shape1=beta_par[1]+x1,shape2=beta_par[2]+n1-x1))*
    as.numeric((1-pbeta(Dcut,beta_par[1]+x1+i,beta_par[2]+N-x1-i))>PostProb);
  PP<-sum(PP_i)
  return(PP)
}

