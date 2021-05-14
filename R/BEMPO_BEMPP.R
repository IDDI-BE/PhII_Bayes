
#--------------------------------------------------------------------------------------------#
# Functions to calculate decision rules for posterior probability, for efficacy and futility #
#--------------------------------------------------------------------------------------------#

# BEMPO interim efficacy stopping rule is P(Delta>Dcut)>=PostProb
# PostProb=1 means never stopping for efficacy
# PostProb=0 means always stopping for efficacy
#-------------------------------------------------------------------

POST_dec_eff<-function(n_,PostProb,Dcut,beta_par){ # For efficacy
  x<-P_Bayes<-0
  
  while(P_Bayes<PostProb & x<=n_){
    P_Bayes<-1-pbeta(Dcut,shape1=x+beta_par[1],shape2=n_-x+beta_par[2])
    x<-x+1
  }
  
  if (PostProb==0 | PostProb==1){return(list(P_Bayes=NA,x=NA))}
  else {return(list(P_Bayes=P_Bayes,x=x-1))}
}

# Check: POST_dec_eff(n_=5,PostProb=0.9,Dcut=0.3,beta_par=c(0.5,0.5))  # default example on https://biostatistics.mdanderson.org/shinyapps/BEMPR/
# Check: POST_dec_eff(n_=10,PostProb=0.9,Dcut=0.3,beta_par=c(0.5,0.5))
# Check: POST_dec_eff(n_=15,PostProb=0.8,Dcut=0.3,beta_par=c(0.5,0.5))
# Check: POST_dec_eff(n_=5,PostProb=1,Dcut=0.3,beta_par=c(0.5,0.5)) #NA
# Check: POST_dec_eff(n_=5,PostProb=0,Dcut=0.3,beta_par=c(0.5,0.5)) #NA

# BEMPO interim futility stopping rule is P(Delta<=D_cut)>PostProb
# PostProb=1 means never stopping for futility
# PostProb=0 means always stopping for futility
#---------------------------------------------------------------------

POST_dec_fut<-function(n_,PostProb,Dcut,beta_par){ # For futility
  x<-n_
  P_Bayes<-0
  
  while(P_Bayes<=PostProb & x>=0){
    P_Bayes<- pbeta(Dcut,shape1=x+beta_par[1],shape2=n_-x+beta_par[2])
    x<-x-1
  }
  
  if (PostProb==0 | PostProb==1){return(list(P_Bayes=NA,x=NA))}
  else {return(list(P_Bayes=P_Bayes,x=x+1))}
}

# Check: POST_dec_fut(n_=5,PostProb=0.7,Dcut=0.3,beta_par=c(0.5,0.5)) # default example on https://biostatistics.mdanderson.org/shinyapps/BEMPR/
# Check: POST_dec_fut(n_=10,PostProb=0.7,Dcut=0.3,beta_par=c(0.5,0.5))
# Check: POST_dec_fut(n_=5,PostProb=0,Dcut=0.3,beta_par=c(0.5,0.5)) #NA
# Check: POST_dec_fut(n_=5,PostProb=1,Dcut=0.3,beta_par=c(0.5,0.5)) #NA


#-------------------------------------------------------------------------------------------------------#
# Functions to calculate decision rules for posterior predictive probability, for efficacy and futility #
#-------------------------------------------------------------------------------------------------------#

# BEMPR interim efficacy stopping rule is PP>=PredProb
#
# PredProb=1 means never stopping for efficacy
# PredProb=0 means always stopping for efficacy
#-------------------------------------------

PP_dec_eff<-function(PredProb,...){ # For efficacy
  args<-list(...)
  x<-PP<-0
  
  while(PP<PredProb & x<=args$n1){
    PP<-Pred_Prob(N=args$N,n1=args$n1,x1=x,beta_par=args$beta_par,Dcut=args$Dcut,PostProb=args$PostProb)
    x<-x+1
  }
  
  if (PredProb==0 | PredProb==1){return(list(PP=NA,x=NA))}
  else {return(list(PP=PP,x=x-1))}
}

# Check: PP_dec_eff(PredProb=0.9,N=15,n1=5,beta_par=c(0.5,0.5),Dcut=0.3,PostProb=0.7) # default example on https://biostatistics.mdanderson.org/shinyapps/BEMPR/
# Check: PP_dec_eff(PredProb=0.9,N=15,n1=10,beta_par=c(0.5,0.5),Dcut=0.3,PostProb=0.7)
# Check: POST_dec_eff(n_=15,PostProb=0.7,Dcut=0.3,beta_par=c(0.5,0.5))
# Check: PP_dec_eff(PredProb=0,N=15,n1=10,beta_par=c(0.5,0.5),Dcut=0.3,PostProb=0.7) #NA
# Check: PP_dec_eff(PredProb=1,N=15,n1=10,beta_par=c(0.5,0.5),Dcut=0.3,PostProb=0.7) #NA

# BEMPR interim futility stopping rule is PP<PredProb
#
# PredProb=0 means never stopping for futility
# PredProb=1 means always stopping for futility
#---------------------------------------------------

PP_dec_fut<-function(PredProb,...){ # For futility
  args<-list(...)
  x<-args$n1
  PP<-1
  
  while(PP>=PredProb & x>=0){
    PP<-Pred_Prob(N=args$N,n1=args$n1,x1=x,beta_par=args$beta_par,Dcut=args$Dcut,PostProb=args$PostProb)
    x<-x-1
  }
  if (PredProb==0 | PredProb==1){return(list(PP=NA,x=NA))}
  else {return(list(PP=PP,x=x+1))}
  
}


# Check: PP_dec_fut(PredProb=0.3,N=15,n1=5 ,beta_par=c(0.5,0.5),Dcut=0.3,PostProb=0.7) # default example on https://biostatistics.mdanderson.org/shinyapps/BEMPR/
# Check: PP_dec_fut(PredProb=0.3,N=15,n1=10,beta_par=c(0.5,0.5),Dcut=0.3,PostProb=0.7)
# Check: PP_dec_fut(PredProb=0  ,N=15,n1=10,beta_par=c(0.5,0.5),Dcut=0.3,PostProb=0.7) #NA
# Check: PP_dec_fut(PredProb=1  ,N=15,n1=10,beta_par=c(0.5,0.5),Dcut=0.3,PostProb=0.7) #NA

#' @title Bayesian efficacy monitoring
#' @description
#' Single arm Bayesian Efficacy Monitoring Via Predictive Probability (BEMPR) or Bayesian Efficacy Monitoring Via Posterior Probability (BEMPO) \cr
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
#' @param p response rate to investigate operating characteristics: can be a scalar or a vector
#' @param design choice between "BEMPR" (predictive probability) and "BEMPO" (posterior probability)
#' @param interim_type choice between 3 character types: "perc", "cont" or "fix". \cr
#' Only one sequence at once can be investigated: \cr
#' 1) "perc": percentage(s) of information, e.g. c(0.5,0.75), specified as a vector via the 'interim' parameter \cr
#' 2) "fix": fixed number of patients, e.g. c(10,20), specified as a vector via the 'interim' parameter \cr
#' 3) "cont": continuous monitoring. Run-in without monitoring (e.g. first 10) is specified in the "interimstart" \cr
#' parameter. Run-out without monitoring (e.g. last 5) is specified in the "interimstop" parameter.
#' @param interim This must be a vector, only to be specified when the parameter 'interim type' is "perc" or "fix" \cr
#' e.g. c(0.5,0.75) when 'interim type' is "perc", or e.g. c(10,20) when 'interim type' is "fix"
#' @param cohortsize This must be a scalar, indicating the frequency of continuous monitoring. \cr
#' Only needs to be defined in case the parameter 'interim_type'="cont"
#' @param interimstart This must be a scalar, specifying run-in without monitoring. \cr
#' Only needs to be defined in case the parameter 'interim_type'="cont". Default is 10
#' @param interimstop This must be a scalar, specifying run-out without monitoring. \cr
#' Only needs to be defined in case the parameter 'interim_type'="cont". Default is 5.
#' @param Delta_fut This may be a scalar or a vector. Must only be specified for BEMPO \cr
#' @param P_fut This may be a scalar or a vector \cr
#' @param Delta_eff This may be a scalar or a vector. Must only be specified for BEMPO \cr
#' @param P_eff This may be a scalar or a vector \cr
#' @param Delta_fin This must be a scalar\cr
#' @param P_fin This must be a scalar\cr
#' @param Beta_dis two parameters Beta(alpha,beta) of the prior Beta distribution
#' @param nsim number of simulations
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

# N=15;p=0.5; design="BEMPO"; interim_type="fix";interim=c(5,10);Delta_fut=0.3;P_fut=0.7; Delta_eff=0.3;P_eff=0.9;Delta_fin=0.3;P_fin=0.8;Beta_dis=c(0.5,0.5);nsim=10
# N=15;p=0.5; design="BEMPR"; interim_type="fix";interim=c(5,10);              P_fut=0.7;               P_eff=0.9;Delta_fin=0.3;P_fin=0.8;Beta_dis=c(0.5,0.5);nsim=10

BEMPP<-function(N,p,design,interim_type,interim,cohortsize=NULL,interimstart=NULL,interimstop=NULL,Delta_fut=NULL,P_fut,Delta_eff=NULL,P_eff,Delta_fin,P_fin,Beta_dis,nsim){

  # Create dataframe with all input parameters + empty vectors for operating characteristics

  if (design=="BEMPR"){
    param<-cbind(paste0("type=",design,
                        "/interim_type=",paste0(interim_type,":(",paste(interim,collapse =',')),")",
                        "/P_fut=(",paste(P_fut,collapse =','),")",
                        "/P_eff=(",paste(P_eff,collapse =','),")",
                        "/D_fin=(",Delta_fin,")",
                        "/P_fin=(",P_fin,")",
                        "/Beta=(",paste(Beta_dis,collapse =','),")",
                        "/nsim=",nsim))
  }

  if (design=="BEMPO"){
    param<-cbind(paste0("type=",design,
                        "/interim_type=",interim_type,
                        "/D_fut=(",paste(Delta_fut,collapse =','),")",
                        "/P_fut=(",paste(P_fut,collapse =','),")",
                        "/D_eff=(",paste(Delta_eff,collapse =','),")",
                        "/P_eff=(",paste(P_eff,collapse =','),")",
                        "/D_fin=(",Delta_fin,")",
                        "/P_fin=(",P_fin,")",
                        "/Beta=(",paste(Beta_dis,collapse =','),")",
                        "/nsim=",nsim))
  }

  param<-data.frame(data.table::CJ(p,N,param))
  param<-cbind(param,pow=0,eff_fin=0,eff_stop=0,fut=0,fut_fin=0,fut_stop=0,N_avg=NA,scenario=1:dim(param)[1]) # pow=proportion(simulations) where decision efficacy at interim or final analysis
                                                                                                  # eff_fin=proportion(simulations) where decision efficacy at final analysis
                                                                                                  # eff_stop=proportion(simulations) where decision efficacy at interim analysis
                                                                                                  # fut=proportion(simulations) where decision of futility at interim or final analysis
                                                                                                  # fut_fin=proportion(simulations) where decision futility at final analysis
                                                                                                  # fut_stop=proportion(simulations) where decision futility at interim analysis
                                                                                                  # N_avg=average number of patients over all simulations
                                                                                                  # scenario= scenario (combination of p and N)
  # Create scenario number
  scenario<-0

  for (a in 1:length(p)){
    for (b in 1:length(N)){
      p.l<-p[a]
      N.l<-N[b]
      scenario<-scenario+1

      # Get number of patients for each interim analysis, possibly specific by N parameter (in case of 'perc' and 'cont')
      #------------------------------------------------------------------------------------------------------------------

      if (interim_type=="perc"){
        interim.f<-floor(interim*N.l)
        param[,"interim"]<-paste(interim.f,collapse =',')
      }

      if (interim_type=="cont"){
        interim.f<-seq(interimstart,N.l-interimstop,cohortsize)
        param[,"interim"]<-paste0(interimstart,"-",N.l-interimstop," by",cohortsize)
      }

      if (interim_type=="fix"){
        interim.f<-interim
        param[,"interim"]<-paste(interim.f,collapse =',')
      }

      param[,"N_interim"]<-length(interim.f)  # Number of interims (without final) for in table

      if (length(Delta_fut)==1){Delta_fut<-rep(Delta_fut,length(interim.f))} # If P_fut is a scalar, change it into a vector
      if (length(Delta_eff)==1){Delta_eff<-rep(Delta_eff,length(interim.f))} # If P_eff is a scalar, change it into a vector

      if (length(P_fut)==1)    {P_fut    <-rep(P_fut    ,length(interim.f))} # If P_fut is a scalar, change it into a vector
      if (length(P_eff)==1)    {P_eff    <-rep(P_eff    ,length(interim.f))} # If P_eff is a scalar, change it into a vector

      # Get decision rules for futility at each interim
      #------------------------------------------------

      if (design=="BEMPR"){

        F_decision<- data.frame(n=interim.f,x=NA,PP=NA,scenario=scenario)
        for (i in 1:length(interim.f)){
          result<- PP_dec_fut(PredProb=P_fut[i],N=N.l,n1=interim.f[i],beta_par=Beta_dis,Dcut=Delta_fin,PostProb=P_fin)
          F_decision[i,"x"]  <- result$x
          F_decision[i,"PP"] <- result$PP
        }

        E_decision<- data.frame(n=c(interim.f,N.l),x=NA,P_Bayes_PP=NA,scenario=scenario)
        for (i in 1:length(interim.f)){
          result<- PP_dec_eff(PredProb=P_eff[i],N=N.l,n1=interim.f[i],beta_par=Beta_dis,Dcut=Delta_fin,PostProb=P_fin)
          E_decision[i,"x"]          <- result$x
          E_decision[i,"P_Bayes_PP"] <- result$PP
        }
        result<- POST_dec_eff(n_=N.l,PostProb=P_fin,Dcut=Delta_fin,beta_par=Beta_dis)
        E_decision[length(interim.f)+1,"x"]          <- result$x
        E_decision[length(interim.f)+1,"P_Bayes_PP"] <- result$P_Bayes
      }

      if (design=="BEMPO"){

        F_decision<- data.frame(n=interim.f,x=NA,P_Bayes=NA,scenario=scenario)
        for (i in 1:length(interim.f)){
          result<- POST_dec_fut(n_=F_decision[i,"n"],PostProb=P_fut[i],Dcut=Delta_fut[i],beta_par=Beta_dis)
          F_decision[i,"x"]      <-result$x
          F_decision[i,"P_Bayes"]<-result$P_Bayes
        }


        E_decision<- data.frame(n=c(interim.f,N.l),x=NA,P_Bayes=NA,scenario=scenario)
        for (i in 1:length(interim.f)){
          result<- POST_dec_eff(n_=E_decision[i,"n"],PostProb=P_eff[i],Dcut=Delta_eff[i],beta_par=Beta_dis)
          E_decision[i,"x"]      <-result$x
          E_decision[i,"P_Bayes"]<-result$P_Bayes
        }
        result<- POST_dec_eff(n_=N.l,PostProb=P_fin,Dcut=Delta_fin,beta_par=Beta_dis)
        E_decision[length(interim.f)+1,"x"]      <-result$x
        E_decision[length(interim.f)+1,"P_Bayes"]<-result$P_Bayes
      }

      Cut_Fut<-F_decision$x  # For simulations
      Cut_Eff<-E_decision$x

      if (scenario==1){    # For output
        Futility<- F_decision
        Efficacy<- E_decision
      }

      if (scenario>1){
        Futility<- rbind(Futility,F_decision)
        Efficacy<- rbind(Efficacy,E_decision)
      }

      # Start simulations
      sim<-data.frame(sim=1:nsim,fut_int=0,fut_int_nr=0,eff_int=0,eff_int_nr=0,fin=0,N_actual=NA)

      if (nsim>0){
        for (c in 1:nsim){

          data<-rbinom(n=N.l,size=1,prob=p.l)

          for (d in 1:length(interim.f)){

            if (sim[c,]$fut_int==0 & sim[c,]$eff_int==0){ # if no decision of futility or efficacy yet

              sim[c,]$N_actual<-interim.f[d]
              succ<-sum(data[1:interim.f[d]])

              if ( (design=="BEMPR" & !P_fut[d]==0) | (design=="BEMPO" & !P_fut[d]==1) ) { # Conditions for no futility interim
                if (succ<=Cut_Fut[d]){ # first check interim futility
                  sim[c,]$fut_int_nr <- d   # futility indicator for specific interim analysis
                  sim[c,]$fut_int    <- 1}  # futility indicator at interim
              }


              if ( (design=="BEMPR" & !P_eff[d]==1) | (design=="BEMPO" & !P_eff[d]==1) ) { # Conditions for no efficacy interim
                if (succ>=Cut_Eff[d]){ # then check interim efficacy
                  sim[c,]$eff_int_nr <- d   # futility indicator for specific interim analysis
                  sim[c,]$eff_int    <- 1}  # futility indicator at interim
              }

            } # 'if' accolade (no decision of futility or efficacy yet)

          } #d-loop [number of interim analyses]

          if (sim[c,]$fut_int==0 & sim[c,]$eff_int==0){ # Only final analysis if no stop at interim analyses for futility or efficacy
            succ<-sum(data)
            sim[c,]$fin<-as.numeric(succ>=Cut_Eff[length(Cut_Eff)])
            sim[c,]$N_actual<-N.l}

        } #c-loop [number of simulations]
      } # 'if' accolade (nsim>0)

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

    } #b-loop [N]

  } #a-loop [p]

  # Order dataset, and return results
  firstcols<-c("scenario","N","p","param","N_interim","interim","pow","eff_fin","eff_stop","fut","fut_fin","fut_stop","N_avg")
  lastcols<-colnames(param)[(! colnames(param) %in% firstcols)]

  return(list(param_simul=param[,c(firstcols,lastcols)], Fut_rules=Futility, Eff_rules=Efficacy))
} # end function