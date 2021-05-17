

#' @title Bayesian efficacy monitoring for the difference of two proportions
#' @description
#' Comparative Bayesian Efficacy Monitoring Via Predictive Probability (BEMPR) or Bayesian Efficacy Monitoring Via Posterior Probability (BEMPO) \cr
#' The frequentist properties can be investigated for combinations of \cr
#' 1) Different sample sizes
#' 2) Different true response rates
#' #
#' BEMPO interim futility stopping rule is P(Delta<=Delta_fut)>P_fut
#' BEMPO interim efficacy stopping rule is P(Delta>Delta_eff)>=P_eff
#' BEMPO final efficacy stopping rule is   P(Delta>Delta_fin)>=P_fin
#'
#' BEMPR interim futility stopping rule is PP<P_fut
#' BEMPR interim efficacy stopping rule is PP>=P_eff
#' BEMPR final efficacy stopping rule is P(Delta>Delta_fin)>=P_fin
#' @param N sample size to investigate operating characteristics: can be a scalar or a vector
#' @param ar allocation ratio experimental:control arm
#' @param p_exp response rate  experimental arm to investigate operating characteristics: can be a scalar or a vector
#' @param p_ctrl response rate  control arm to investigate operating characteristics: can be a scalar or a vector
#' @param design_fut choice between "futPR" (predictive probability) and "futPO" (posterior probability)
#' @param design_eff choice between "effPR" (predictive probability) and "effPO" (posterior probability)
#' @param interim Information fraction at interim, e.g. c(0.5,0.75)
#' @param Delta_fut This may be a scalar or a vector. Must only be specified for BEMPO \cr
#' @param P_fut This may be a scalar or a vector \cr
#' @param Delta_eff This may be a scalar or a vector. Must only be specified for BEMPO \cr
#' @param P_eff This may be a scalar or a vector \cr
#' @param Delta_fin This must be a scalar\cr
#' @param P_fin This must be a scalar\cr
#' @param beta_par_exp two shape parameters c(alpha,beta) for prior beta distribution experimental arm
#' @param beta_par_ctrl two shape parameters c(alpha,beta) for prior beta distribution control arm
#' @param nsim number of simulations
#' @param distrisize Size of sampled distributions (the larger, the better)
#' @param PP_nsim Number of simulation to sample from predictive distribution
#' @return a list of 3 data.frames: first with design parameters + operating characteristics ($param_simul) , \cr
#' one with futility decision rules ($Fut_rules), and one with efficacy decision rules ($Eff_rules)
#' \itemize{
#' \item param_simul: dataframe with all input parameters + operating characteristics
#' \item param_simul$scenario: combination of vector N and p
#' \item param_simul$N: sample size for which frequentist properties simulated
#' \item param_simul$p: binomial parameter (true proportion) for which frequentist properties simulated
#' \item param_simul$param: string with all input parameters
#' \item param_simul$N_interim: number of interim analyses (final analysis not included)
#' \item param_simul$interim: string with number of patients per interim analysis
#' \item param_simul$pow: proportion of simulations where \eqn{H0} was rejected
#' \item param_simul$eff_fin: proportion of simulations where \eqn{H0} was rejected only at the final analysis
#' \item param_simul$eff_stop: proportion of simulations where \eqn{H0} was rejected at an interim analysis
#' \item param_simul$fut: proportion of simulations where \eqn{H0} was not rejected
#' \item param_simul$eff_fin: proportion of simulations where decision of not rejecting \eqn{H0} was at the final analysis
#' \item param_simul$fut_stop: proportion of simulations where trial stopped for futility at an interim analysis
#' \item param_simul$N_avg: average number of patients for simulated datasets
#' \item param_simul$int_nrx: number of patients at x'th interim analysis
#' \item param_simul$effx: proportion of simulations where trial stopped for efficacy at xth interim analysis
#' \item param_simul$futx: proportion of simulations where trial stopped for futility at xth interim analysis
#' \item param_simul$decision rules: for efficacy: if >=x: stop for efficacy; for futility: if <=x: stop for futility
#' \item Fut_rules: dataframe with futility rules for each interim
#' \item Fut_rules$n: number of patients at interim
#' \item Fut_rules$x: number of successes at interim (Futility if observed successes <=x)
#' \item Fut_rules$P_Bayes/Fut_rules$PP: actual posterior probability or PP at cutoff
#' \item Eff_rules: dataframe with efficacy rules for each interim
#' \item Eff_rules$n: number of patients at interim
#' \item Eff_rules$x: number of successes at interim (Efficacy if observed successes >=x)
#' \item Eff_rules$P_Bayes/Fut_rules$PP: actual posterior probability or PP at cutoff
#' }
#'@references Lee JJ, Liu DD.A predictive probability design for phase II cancer clinical trialsClinical Trials 2008; 5: 93–106
#'@importFrom VGAM dbetabinom.ab
#'@importFrom utils txtProgressBar 
#'@importFrom utils setTxtProgressBar
#'@export
#'
#' @examples #Check versus https://biostatistics.mdanderson.org/shinyapps/BEMPO/
#' test1<-BEMPP(N=15,p=0.5, design="BEMPO", interim_type="fix",interim=c(5,10),Delta_fut=0.3,
#'  P_fut=0.7,Delta_eff=0.3,P_eff=0.9,Delta_fin=0.3,P_fin=0.8,Beta_dis=c(0.5,0.5),nsim=10)
#' test2<-BEMPP(N=15,p=0.5, design="BEMPO", interim_type="fix",interim=c(5,10),Delta_fut=0.3,
#'  P_fut=1  ,Delta_eff=0.3,P_eff=1  ,Delta_fin=0.3,P_fin=0.8,Beta_dis=c(0.5,0.5),nsim=10)
#' test3<-BEMPP(N=15,p=0.5, design="BEMPR", interim_type="fix",interim=c(5,10),P_fut=0.3,
#'  P_eff=0.9,Delta_fin=0.3,P_fin=0.7,Beta_dis=c(0.5,0.5),nsim=10)
#' test4<-BEMPP(N=15,p=0.5, design="BEMPR", interim_type="fix",interim=c(5,10),P_fut=0,P_eff=1,
#'  Delta_fin=0.3,P_fin=0.8,Beta_dis=c(0.5,0.5),nsim=10)


#-----------------------------------------------------------------------------------------------------------------------#
# Function to calculate                                                                                                 #
# 1) Decision rules                                                                                                     #
# 2) operating characteristics for Bayesian Efficacy Monitoring Via Predictive (BEMPR) or Posterior Probability (BEMPO) #
#-----------------------------------------------------------------------------------------------------------------------#

# N=c(100,100);p_exp=c(0.5,0.55);p_ctrl=c(0.4,0.45); interim=c(0.25,0.5,0.75); ar=1; Delta_fut=0.05;P_fut=0.3;Delta_eff=0.2;P_eff=0.8;Delta_fin=0.3;P_fin=0.95; 
# distrisize=10^3; PP_nsim=10^3; design_fut="futPO"; design_eff="effPO";nsim=100

BM_R<-function(N,ar,p_exp,p_ctrl,design_fut,design_eff,interim,Delta_fut=NULL,P_fut,Delta_eff=NULL,P_eff,Delta_fin,P_fin,beta_par_exp,beta_par_ctrl,nsim,distrisize=10^3,PP_nsim=10^3){

  # Checks
  
  if ((length(p_exp) != length(p_ctrl)) | (length(p_exp) != length(N)) | (length(p_ctrl) != length(N))){stop("N, p_exp and p_ctrl input vectors must be of equal length")}
  
  # Create dataframe with empty vectors for operating characteristics

  param<-cbind(N,p_exp,p_ctrl,pow=0,eff_fin=0,eff_stop=0,fut=0,fut_fin=0,fut_stop=0,N_avg=NA,scenario=1:length(N)) # pow=proportion(simulations) where decision efficacy at interim or final analysis
                                                                                                  # eff_fin=proportion(simulations) where decision efficacy at final analysis
                                                                                                  # eff_stop=proportion(simulations) where decision efficacy at interim analysis
                                                                                                  # fut=proportion(simulations) where decision of futility at interim or final analysis
                                                                                                  # fut_fin=proportion(simulations) where decision futility at final analysis
                                                                                                  # fut_stop=proportion(simulations) where decision futility at interim analysis
                                                                                                  # N_avg=average number of patients over all simulations
                                                                                                  # scenario= scenario (combination of p and N)
  
  # Create scenario number
  scenario<-0

  for (a in 1:length(N)){

    p_exp.l<-p_exp[a]
    p_ctrl.l<-p_ctrl[a]
    N.l<-N[a]
    
    N_ctrl.l<-floor(N.l*(1/(ar+1)))
    N_exp.l <-N.l-N_ctrl.l
    
    scenario<-scenario+1
    print(paste0("scenario ",scenario,":p_exp=",p_exp.l,";p_ctrl=",p_ctrl.l,";N=",N.l))

    # Get number of patients for each interim analysis, possibly specific by N parameter (in case of 'perc' and 'cont')
    #------------------------------------------------------------------------------------------------------------------
    
    interim.f<-floor(interim*N.l)
    param[,"interim"]<-paste(interim.f,collapse =',')

    param[,"N_interim"]<-length(interim.f)  # Number of interims (without final) for in table

    if (length(Delta_fut)==1){Delta_fut<-rep(Delta_fut,length(interim.f))} # If P_fut is a scalar, change it into a vector
    if (length(Delta_eff)==1){Delta_eff<-rep(Delta_eff,length(interim.f))} # If P_eff is a scalar, change it into a vector

    if (length(P_fut)==1)    {P_fut    <-rep(P_fut    ,length(interim.f))} # If P_fut is a scalar, change it into a vector
    if (length(P_eff)==1)    {P_eff    <-rep(P_eff    ,length(interim.f))} # If P_eff is a scalar, change it into a vector


    # Start simulations
    sim<-data.frame(sim=1:nsim,fut_int=0,fut_int_nr=0,eff_int=0,eff_int_nr=0,fin=0,N_actual=NA)
    pb <-  txtProgressBar(min = 0, max = nsim, style = 3) # set progress bar

    for (b in 1:nsim){
      
      setTxtProgressBar(pb, b)
      
      exp  <- rbinom(n=N_exp.l ,size=1,prob=p_exp.l )
      ctrl <- rbinom(n=N_ctrl.l,size=1,prob=p_ctrl.l)

      for (c in 1:length(interim.f)){

        n1_ctrl.l<-floor(interim.f[c]*(1/(ar+1)))
        n1_exp.l <-interim.f[c]-n1_ctrl.l
        
        if (sim[b,]$fut_int==0 & sim[b,]$eff_int==0){ # if no decision of futility or efficacy yet

          sim[b,]$N_actual<-interim.f[c]
          succ_exp.l <-sum(exp [1:n1_exp.l ]);ph_exp.l =succ_exp.l /n1_exp.l
          succ_ctrl.l<-sum(ctrl[1:n1_ctrl.l]);ph_ctrl.l=succ_ctrl.l/n1_ctrl.l

          if (design_fut=="futPR" & !P_fut[c]==0) { # Conditions for no futility interim
            
            PP<-Pred_Prob_R(p_exp=ph_exp.l,p_ctrl=ph_ctrl.l,N_exp=N_exp.l,N_ctrl=N_ctrl.l,n1_exp=n1_exp.l,n1_ctrl=n1_ctrl.l,distrisize=distrisize,nsim=PP_nsim,PostProb=P_fin,Dcut=Delta_fin,beta_par_exp=c(1,1),beta_par_ctrl=c(1,1))
            if (PP<P_fut[c]){# first check interim futility
              sim[c,]$fut_int_nr <- c   # futility indicator for specific interim analysis
              sim[c,]$fut_int    <- 1}  # futility indicator at interim
          }
          
          if (design_fut=="futPO" & !P_fut[c]==1) { # Conditions for no futility interim
            
            PO<-Post_Prob_R(p_exp=ph_exp.l,p_ctrl=ph_ctrl.l,n_exp=n1_exp.l,n_ctrl=n1_ctrl.l,distrisize=distrisize,Dcut=Delta_fut[c],beta_par_exp=c(1,1),beta_par_ctrl=c(1,1))
            if (PO$prob<P_fut[c]){# first check interim futility
              sim[c,]$fut_int_nr <- c   # futility indicator for specific interim analysis
              sim[c,]$fut_int    <- 1}  # futility indicator at interim
          }


          if ( (design_eff=="effPR" & !P_eff[c]==1)) { # Conditions for no efficacy interim
            
            PP<-Pred_Prob_R(p_exp=ph_exp.l,p_ctrl=ph_ctrl.l,N_exp=N_exp.l,N_ctrl=N_ctrl.l,n1_exp=n1_exp.l,n1_ctrl=n1_ctrl.l,distrisize=distrisize,nsim=PP_nsim,PostProb=P_fin,Dcut=Delta_fin,beta_par_exp=c(1,1),beta_par_ctrl=c(1,1))
            if (PP>=P_eff[c]){# then check interim efficacy
              sim[c,]$eff_int_nr <- c   # futility indicator for specific interim analysis
              sim[c,]$eff_int    <- 1}  # futility indicator at interim
          }
          
          if ( (design_eff=="effPO" & !P_eff[c]==1)) { # Conditions for no efficacy interim
            
            PO<-Post_Prob_R(p_exp=ph_exp.l,p_ctrl=ph_ctrl.l,n_exp=n1_exp.l,n_ctrl=n1_ctrl.l,distrisize=distrisize,Dcut=Delta_eff[c],beta_par_exp=c(1,1),beta_par_ctrl=c(1,1))
            if (PO$prob>=P_eff[c]){# first check interim futility
              sim[c,]$eff_int_nr <- c   # futility indicator for specific interim analysis
              sim[c,]$eff_int    <- 1}  # futility indicator at interim
          }

        } # 'if' accolade (no decision of futility or efficacy yet)

      } #c-loop [number of interim analyses]

      if (sim[c,]$fut_int==0 & sim[c,]$eff_int==0){ # Only final analysis if no stop at interim analyses for futility or efficacy
        
        succ_exp.l <-sum(exp [1:N_exp.l ]);ph_exp.l  =succ_exp.l /N_exp.l
        succ_ctrl.l<-sum(ctrl[1:N_ctrl.l]);ph_ctrl.l =succ_ctrl.l/N_ctrl.l
        
        PO<-Post_Prob_R(p_exp=p_exp.l,p_ctrl=p_ctrl.l,n_exp=n1_exp.l,n_ctrl=n1_ctrl.l,distrisize=distrisize,Dcut=Delta_fin,beta_par_exp=c(1,1),beta_par_ctrl=c(1,1))
        
        sim[c,]$fin<-as.numeric(PO$prob>=P_fin)
        sim[c,]$N_actual<-N.l}

    } #c-loop [number of simulations]

    param[scenario,]$pow      <- (sum(sim$fin)+sum(sim$eff_int))/nsim  # Frequency efficacy decision (whether at final or at interim)
    param[scenario,]$eff_fin  <- sum(sim$fin)/nsim                     # Frequency efficacy decision at final
    param[scenario,]$eff_stop <- sum(sim$eff_int)/nsim                 # Frequency efficacy decision at interim
    param[scenario,]$fut      <- (nsim-sum(sim$fin)-sum(sim$eff_int))/nsim
    param[scenario,]$fut_fin  <- (nsim-sum(sim$fin)-sum(sim$eff_int)-sum(sim$fut_int))/nsim
    param[scenario,]$fut_stop <- sum(sim$fut_int)/nsim                 # Frequency futility decision at interim
    param[scenario,]$N_avg    <- mean(sim$N_actual)                    # Average number of patients

    for(i in 1:length(interim.f)){
      param[scenario,paste0("int_nr",i)] <- interim.f[i]
      param[scenario,paste0("eff"   ,i)] <- sum(sim$eff_int_nr==i)/nsim
      param[scenario,paste0("fut"   ,i)] <- sum(sim$fut_int_nr==i)/nsim
    }

  } #a-loop [N]

  # Order dataset, and return results
  firstcols<-c("scenario","N","p_exp","p_ctrl","N_interim","interim","pow","eff_fin","eff_stop","fut","fut_fin","fut_stop","N_avg")
  lastcols<-colnames(param)[(! colnames(param) %in% firstcols)]

  return(param[,c(firstcols,lastcols)])
} # end function

#res<-BM_R(N=c(100,100),ar=1, p_exp=c(0.5,0.55),p_ctrl=c(0.4,0.45), interim=c(0.25,0.5,0.75), Delta_fut=0.05,P_fut=0.3,Delta_eff=0.2,P_eff=0.8,Delta_fin=0.3,P_fin=0.95, 
#distrisize=10^3, PP_nsim=10^3, design_fut="futPO", design_eff="effPO",nsim=1000)