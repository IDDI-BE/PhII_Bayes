
#' @title Bayesian predictive probability, given an observed difference of 2 proportions
#'
#' @description This function calculates the Bayesian predictive probability 
#'     P(Delta>Dcut|N_exp,N_ctrl,n1_exp,n1_ctrl,x1) with Delta=Delta_exp-Delta_ctrl
#' @param n1_exp number of patients in experimental arm (scalar) at interim
#' @param n1_ctrl number of patients in control arm (scalar) at interim
#' @param N_exp number of patients in experimental arm (scalar) at end of study
#' @param N_ctrl number of patients in control arm (scalar) at end of study
#' @param p_exp observed proportion experimental arm
#' @param p_ctrl observed proportion control arm
#' @param distrisize Size of sampled distributions (the larger, the better)
#' @param nsim Number of simulation to sample from predictive distribution
#' @param Dcut True difference between two proportions (can be a vector)
#' @param PostProb Threshold for outcome at the end of the study, in terms of Bayesian posterior probability
#'     P(theta>Dcut|x1+x2)>PostProb, with x1 the difference in proportions at interim and x2 at the final
#' @param beta_par_exp two shape parameters c(alpha,beta) for prior beta distribution experimental arm
#' @param beta_par_ctrl two shape parameters c(alpha,beta) for prior beta distribution control arm
#' @return Predictive Probability for outcome at the end of the study
#' @examples
#' Pred_Prob_R(p_exp=0.475,p_ctrl=0.475-0.08,N_exp=100,N_ctrl=100,n1_exp=50,n1_ctrl=50,
#' distrisize=10^3,nsim=10^3,PostProb=0.83,Dcut=0,beta_par_exp=c(1,1),beta_par_ctrl=c(1,1))
#' @importFrom stats rbeta
#' @importFrom utils txtProgressBar 
#' @importFrom utils setTxtProgressBar
#' @export

Pred_Prob_R <- function(p_exp,p_ctrl,N_exp,N_ctrl,n1_exp,n1_ctrl,distrisize=1000,nsim=1000,Dcut,PostProb,beta_par_exp,beta_par_ctrl){
  
  postprob_sim <- rep(NA,nsim) # Empty vector for filling in 'for' loop
  pb <-  txtProgressBar(min = 0, max = nsim, style = 3)
  
  for (i in 1:nsim){
    
    setTxtProgressBar(pb, i)
    
    x1_exp  <- p_exp *n1_exp
    x1_ctrl <- p_ctrl*n1_ctrl
    
    # Sample 1 value from Posterior beta distribution at interim and generate results for remainder of study, for 1 simulated study
    sim_exp        <- rbinom(n=N_exp -n1_exp ,size=1,prob=rbeta(n=1,shape1= x1_exp  + beta_par_exp [1], shape2=n1_exp  - x1_exp  + beta_par_exp [2]))
    sim_ctrl       <- rbinom(n=N_ctrl-n1_ctrl,size=1,prob=rbeta(n=1,shape1= x1_ctrl + beta_par_ctrl[1], shape2=n1_ctrl - x1_ctrl + beta_par_ctrl[2]))
    
    # Based on combination (observed data + generated): calculate posterior predictive distribution
    distr_exp_sim  <- rbeta(n=distrisize,shape1= (x1_exp  + sum(sim_exp )) + beta_par_exp [1], shape2=N_exp  - (x1_exp  + sum(sim_exp )) + beta_par_exp [2]);
    distr_ctrl_sim <- rbeta(n=distrisize,shape1= (x1_ctrl + sum(sim_ctrl)) + beta_par_ctrl[1], shape2=N_ctrl - (x1_ctrl + sum(sim_ctrl)) + beta_par_ctrl[2]);
    
    # Distribution of difference (not parametric, so sampled)
    distr_diff_sim <- distr_exp_sim-distr_ctrl_sim
    
    # Calculate posterior probability(true difference>0|observed data at futility combined with simulated data)
    postprob_sim[i]<- sum(distr_diff_sim>Dcut)/distrisize
  }
  
  # Proportion of simulated studies where "posterior probability(true difference>0|observed data at futility combined with simulated data)">=0.8
  predpow<-sum(postprob_sim>=PostProb)/nsim
  return(predpow)
}