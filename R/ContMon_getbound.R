
#' @title Get boundaries for continuous trial monitoring
#'
#' @description This function calculates boundaries for every sample size
#'
#' @details Two parameters should be specified: either 'target' and 'maxtox'
#'     either 'target' and 'P_target'. If the Bayesian posterior probability, given the data
#'     is higher than 'P_target', i.e. P(toxicity rate>target|c+1,n)>=P_target, the boundary is crossed.
#'     For every 'n', the boundary 'c' is given, which corresponds with the maximum allowed
#'     toxicity. So the boundary is crossed at 'c+1'
#' @details P_Bayes= Bayesian posterior probability P(toxicity rate>target|c,n)
#' @details P_obs= c/n or 'allowable proportion': if number of events <=c then continuation, if c+1: boundary crossed
#' @details 'Bayes_result' gives the Bayesian posterior probability P(toxicity rate<=target|c+1,n), so given that
#'     this boundary 'c+1' is exceeded, for assumptions, specified by the 'tox_ass' parameter
#' @details Prob(X>n|n,tox_true=...) gives the Prob(X>=c+1|different assumptions ('toxass'))
#' @details Cum Prob(...) gives the cumulative probability of at least 1x/2x/3x and 4x crossing the boundary
#'     by n patients
#'
#' @param SS total sample size, e.g. calculated for primary efficacy endpoint
#' @param target target toxicity rate
#' @param Beta_dis two parameters Beta(alpha,beta) of the prior Beta distribution
#' @param P_target P(toxicity rate>target|data)>=P_target as decision threshold
#' @param maxtox maximum allowable observed toxicity rate, where boundary not crossed
#' @param prt1 value for number of patients in starting cohort where alternative safety rule is followed, e.g. 3+3
#' @param tox_ass vector of rates for calculating Bayesian probabilities P(toxicity rate<=tox_ass[]|boundary crossed) and frequentist
#'      operating characteristics
#' @param sim binary indicator: simulations for frequentist operating characteristics (0=NO; 1=YES)
#' @param nsim number of simulated datasets
#' @param out path for storing csv output files
#' @examples
#'\dontrun{
#' ContMon_getbound(SS=19,target=0.2,Beta_dis=c(1,1),P_target=NULL,maxtox=0.33,
#'     tox_ass=seq(0.10,0.5,0.05),sim=1,nsim=100,prt1=6,out=NULL)
#' }
#'
#' @export
#' @importFrom stats pbeta pbinom rbinom
#' @importFrom utils write.csv


# Determine cutoffs for safety review
#------------------------------------

# SS=19;target=0.2;P_target=NULL;maxtox=0.33;tox_ass=seq(0.10,0.5,0.05);sim=1;nsim=100; prt1=6;out=NULL

ContMon_getbound<-function(SS,target,Beta_dis,P_target=NULL,maxtox=NULL,prt1=0,tox_ass=seq(0.10,0.5,0.05),sim=0,nsim,out){

  n_ass <-length(tox_ass)

  # Define empty dataframes for results
  #------------------------------------

  result       <- data.frame(target=target,n=1:SS,c=rep(NA,SS),P_target=rep(NA,SS),P_Bayes=rep(NA,SS),P_obs=rep(1,SS))

  Bayes_result <- data.frame(B=matrix(rep(NA,length(tox_ass)*SS),ncol=length(tox_ass)))

  sim_result   <- data.frame(S1=matrix(rep(NA,length(tox_ass)*SS),ncol=length(tox_ass)),  # For frequentist Prob(X>c|n,tox_true=tox_ass[i])
                             S2=matrix(rep(NA,length(tox_ass)*SS),ncol=length(tox_ass)),
                             S3=matrix(rep(NA,length(tox_ass)*SS),ncol=length(tox_ass)),
                             S4=matrix(rep(NA,length(tox_ass)*SS),ncol=length(tox_ass)),
                             S5=matrix(rep(NA,length(tox_ass)*SS),ncol=length(tox_ass)))

  # Determine decision rule
  #------------------------

  if (is.null(P_target)){

    P_target <- 1 # Start value for P(tox>target|data)>=P_target

    while (max(result[(prt1+1):SS,"P_obs"])>maxtox){

      P_target<-P_target-0.01

      for (i in 1:SS){

        c       <- -1
        P_Bayes <- 0

        while (P_Bayes<P_target){

          result[i,"c"]       <- c
          result[i,"P_Bayes"] <- round(P_Bayes,3)

          c                   <- c+1
          P_Bayes<-1-pbeta(target,shape1=c+Beta_dis[1],shape2=i-c+Beta_dis[2])  # P(toxicity>target|c,n), with Beta(1,1) prior

        }
      }

      result[,"P_obs"]    <- round(result$c/result$n,2)
      result[,"P_target"] <- P_target

    }
  }


  if (is.null(maxtox)){

    for (i in 1:SS){

      c       <- -1
      P_Bayes <- 0

      while (P_Bayes<P_target){

        result[i,"c"]       <- c
        result[i,"P_Bayes"] <- round(P_Bayes,3)

        c                   <- c+1
        P_Bayes<-1-pbeta(target,shape1=c+Beta_dis[1],shape2=i-c+Beta_dis[2])  # P(toxicity>target|c,n), with Beta(1,1) prior

      }
    }

    result[,"P_obs"]    <- round(result$c/result$n,2)
    result[,"P_target"] <- P_target

  }


  # Calculate Bayesian posterior probability P(toxicity<=given rate|c+1,n), so at decision to stop, for set of assumptions
  #----------------------------------------------------------------------------------------------------------------------

  c<-result$c
  n<-result$n

  for (i in (prt1+1):SS){
    for (j in 1:length(tox_ass)){
      Bayes_result[i,j]=round(pbeta(tox_ass[j],shape1=(c[i]+1)+Beta_dis[1],shape2=i-(c[i]+1)+Beta_dis[2]),2)
    }
  }

  # Calculate probability of trigger at least once, twice, three times and 4 times in a frequentist setting
  #--------------------------------------------------------------------------------------------------------

  if(sim==1){
    for (i in 1:length(tox_ass)){

      sim_result[,i]<-c(rep(NA,prt1),round(1-pbinom(c[(prt1+1):SS],n[(prt1+1):SS],tox_ass[i]),2))
      names(sim_result)[i]<-paste0("Prob(X>c|n,tox_true=",tox_ass[i],")")

      simdata_1 <-data.frame(matrix(rep(NA,SS*nsim),ncol=nsim))  # Trigger at least once
      simdata_2 <-data.frame(matrix(rep(NA,SS*nsim),ncol=nsim))  # Trigger at least twice
      simdata_3 <-data.frame(matrix(rep(NA,SS*nsim),ncol=nsim))  # Trigger at least 3 times
      simdata_4 <-data.frame(matrix(rep(NA,SS*nsim),ncol=nsim))  # Trigger at least 4 times

      for (j in 1:nsim){
        DRT     <- rbinom(n=SS,size=1,prob=tox_ass[i])
        DRT_cum <- cumsum(DRT)
        exceed  <- DRT_cum>c

        if(prt1!=0){
          simdata_1[,j]<-c(rep(NA,prt1),as.numeric(cumsum(exceed[(prt1+1):SS])>=1))
          simdata_2[,j]<-c(rep(NA,prt1),as.numeric(cumsum(exceed[(prt1+1):SS])>=2))
          simdata_3[,j]<-c(rep(NA,prt1),as.numeric(cumsum(exceed[(prt1+1):SS])>=3))
          simdata_4[,j]<-c(rep(NA,prt1),as.numeric(cumsum(exceed[(prt1+1):SS])>=4))
        }

        if(prt1==0){
          simdata_1[,j]<-c(             as.numeric(cumsum(exceed)             >=1))
          simdata_2[,j]<-c(             as.numeric(cumsum(exceed)             >=2))
          simdata_3[,j]<-c(             as.numeric(cumsum(exceed)             >=3))
          simdata_4[,j]<-c(             as.numeric(cumsum(exceed)             >=4))
        }

        print(paste0(i,"_",j)) # simulations counter
      }

      sim_result[,  n_ass+i] <- round(apply(simdata_1, 1, sum)/nsim,2)
      sim_result[,2*n_ass+i] <- round(apply(simdata_2, 1, sum)/nsim,2)
      sim_result[,3*n_ass+i] <- round(apply(simdata_3, 1, sum)/nsim,2)
      sim_result[,4*n_ass+i] <- round(apply(simdata_4, 1, sum)/nsim,2)

      names(sim_result)[  n_ass+i]<-paste0("Cum Prob(1x X>c|n;n>",prt1,";tox_true=",tox_ass[i],")")
      names(sim_result)[2*n_ass+i]<-paste0("Cum Prob(2x X>c|n;n>",prt1,";tox_true=",tox_ass[i],")")
      names(sim_result)[3*n_ass+i]<-paste0("Cum Prob(3x X>c|n;n>",prt1,";tox_true=",tox_ass[i],")")
      names(sim_result)[4*n_ass+i]<-paste0("Cum Prob(4x X>c|n;n>",prt1,";tox_true=",tox_ass[i],")")
    }

  }

  write.csv(result      ,paste0(out,"result_"      ,Sys.Date(),".csv"))
  write.csv(Bayes_result,paste0(out,"Bayes_result_",Sys.Date(),".csv"))
  if(sim==1){
    write.csv(sim_result  ,paste0(out,"sim_result_"  ,Sys.Date(),".csv"))
  }

  if (sim==0){return(list(results=cbind(result,Bayes_result)           ,assumptions=tox_ass))}
  if (sim==1){return(list(results=cbind(result,Bayes_result,sim_result),assumptions=tox_ass))}
}


