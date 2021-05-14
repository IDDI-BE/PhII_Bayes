
#' @title Bayesian posterior probability, given an observed difference of 2 proportions
#'
#' @description This function calculates the Bayesian posterior probability 
#' P(Delta>Dcut|n,x) with Delta=Delta_exp-Delta_ctrl, with Delta_exp~beta(beta_par_exp[1]+x,beta_par_exp[2]+n_exp-x_exp) 
#' with prior beta distribution beta_par, and Delta_ctrl~beta(beta_par_ctrl[1]+x,beta_par_ctrl[2]+n_ctrl-x_ctrl) 
#' @param n_exp number of patients in experimental arm (scalar)
#' @param n_ctrl number of patients in control arm (scalar)
#' @param x_exp number of successes in experimental arm (scalar)
#' @param x_ctrl number of successes in control arm (scalar)
#' @param beta_par_exp two shape parameters c(alpha,beta) for prior beta distribution experimental arm (scalar)
#' @param beta_par_ctrl two shape parameters c(alpha,beta) for prior beta distribution control arm (scalar)
#' @param p_exp binomial parameter experimental arm, corresponding with some hypothesis (scalar)
#' @param p_ctrl binomial parameter control arm, corresponding with some hypothesis (scalar)
#' @param Dcut Difference between two proportions corresponding with some hypothesis (can be a vector)
#' @param distrisize Size of sampled distributions (the larger, the better)
#' @examples
#' Post_Prob_R(p_exp=0.475,p_ctrl=0.475-5/50,n_exp=50,n_ctrl=50,distrisize=10^3,
#' D_cut=c(0,0.05,0.075,0.1),beta_par_exp=c(1,1),beta_par_ctrl=c(1,1))
#' @export

Post_Prob_R<-function(p_exp,p_ctrl,n_exp,n_ctrl,distrisize,D_cut,beta_par_exp,beta_par_ctrl){
  
  df<-expand.grid(p_exp,p_ctrl,paste0("n_exp:",n_exp,";n_ctrl:",n_ctrl),D_cut)
  names(df)<-c("p_exp","p_ctrl","N","D_cut")
  df$prob<-NA
  counter<-0
  
  distr_exp  <- rbeta(n=distrisize,shape1= p_exp *n_exp  + beta_par_exp [1], shape2=n_exp  - p_exp *n_exp  + beta_par_exp [2])
  distr_ctrl <- rbeta(n=distrisize,shape1= p_ctrl*n_ctrl + beta_par_ctrl[1], shape2=n_ctrl - p_ctrl*n_ctrl + beta_par_ctrl[2])
  distr_diff<-distr_exp-distr_ctrl
        
  for(i in 1:length(D_cut)){
    cut<-D_cut[i]
    prob<-sum(distr_diff>cut)/distrisize
    df[i,"prob"]<-prob
  }
  
  return(df)
}