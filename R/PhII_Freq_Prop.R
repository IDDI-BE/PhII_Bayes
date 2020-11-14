
# Function to calculate predictive probability

Pred_Prob.f<-function(N,n1,x1,beta_par,p,Delta_T_){
  n2<-N-n1
  i<-seq(0,n2)
  PP_i<- (VGAM::dbetabinom.ab(x=i,size=n2,shape1=beta_par[1]+x1,shape2=beta_par[2]+n1-x1))*
    as.numeric((1-pbeta(p,beta_par[1]+x1+i,beta_par[2]+N-x1-i))>Delta_T_);
  PP<-sum(PP_i)
  return(PP)
}

# Function to calculate posterior probability
Post_Prob.f<-function(n,x,beta_par,p){
  1-pbeta(p,x+beta_par[1],n-x+beta_par[2])
}


#' @title PhII Bayesian design with predictive probability
#' @description
#' A flexible PhII, signle-arm study design, based on Bayesian predictive probability \cr
#' The frequentist properties can be investigated of combinations of \cr
#' 1) Different sample sizes
#' 2) Different Delta_T
#' 3) Different boundary (efficacy boundary Delta_U and futility boundary Delta_L are linked \cr
#' and must be specified as a pair of boundaries. So the parameters Delta_U and Delta_L \cr
#' should have exactly the same structure. 
#' @param N this must be a vector: of different total sample sizes. Different N's can be investigated
#' @param Delta_T This must be vector: of different threshold values, defining the posterior \cr
#' probability P(Delta>=p0|data at the final analysis), for a succesful trial. 
#' @param interim_type choice between 3 character types: "Perc", "cont" or "fix". \cr
#' Only one sequence at once can be investigated: \cr
#' 1) "Perc": percentage(s) of information, e.g. c(0.5,0.75), specified as a vector via the 'interim' parameter \cr
#' 2) "fix": fixed number of patients, e.g. c(10,20), specified as a vector via the 'interim' parameter \cr
#' 3) "cont": continuous monitoring. Run-in without monitoring (e.g. first 10) is specified in the "interimstart" \cr
#' parameter. Run-out without monitoring (e.g. last 5) is specified in the "interimstop" parameter.
#' @param interim This must be a vector, only to be specified when the parameter 'interim type' is "Perc" or "fix" \cr
#' e.g. c(0.5,0.75) when 'interim type' is "Perc", or e.g. c(10,20) when 'interim type' is "fix"
#' @param cohortsize This must be a scalar, indicating the frequency of continuous monitoring. \cr
#' Only needs to be defined in case the parameter 'interim_type'="cont"
#' @param interimstart This must be a scalar, specifying run-in without monitoring. \cr
#' Only needs to be defined in case the parameter 'interim_type'="cont". Default is 10
#' @param interimstop This must be a scalar, specifying run-out without monitoring. \cr
#' Only needs to be defined in case the parameter 'interim_type'="cont". Default is 5.
#' @param Delta_L Vector for futility ("Lower") boundary. This may be a scalar, if the parameter "interim"='cont' \cr
#' or if there is only one interim analysis. This corresponds with a threshold value for the predictive probability \cr
#' of a succesful result at the final analysis (Posterior probability >=Delta_T). Default value=0 (no interim for futility).
#' @param Delta_U Vector for efficacy ("Upper") boundary. This may be a scalar, if the parameter "interim"='cont' \cr
#' or if there is only one interim analysis. This corresponds with a threshold value for the predictive probability \cr
#' of a succesful result at the final analysis (Posterior probability >=Delta_T). Default value=1 (no interim for efficacy).
#' @param p0 probability of the uninteresting response (null hypothesis \eqn{H0})
#' @param pa probability of the interesting response (alternative hypothesis Ha)
#' @param Beta_dis two parameters Beta(alpha,beta) of the prior Beta distribution
#' @param Freqprop "pow" for power (simulation under Ha) or "typeI" for simulation under \eqn{H0}
#' @param nsim number of simulations
#' @return a data.frame with elements
#' \itemize{
#' \item N: Different sample sizes for which frequentist properties simulated
#' \item Delta_T: Different thresholds for which frequentist properties simulated
#' \item Boundaries: Different boundaries (combination of Delta_L and Delta_U) for which frequentist properties are simulated
#' \item N_interim: number of interim analyses (final analysis not included)
#' \item pow: proportion of simulations where \eqn{H0} was rejected
#' \item eff_fin: proportion of simulations where \eqn{H0} was rejected only at the final analysis
#' \item eff_stop: proportion of simulations where \eqn{H0} was rejected at an interim analysis
#' \item fut_stop: proportion of simulations where trial stopped for futility at an interim analysis
#' \item int_nrx: number of patients at x'th interim analysis
#' \item effx: proportion of simulations where trial stopped for efficacy at xth interim analysis
#' \item futx: proportion of simulations where trial stopped for futility at xth interim analysis
#' }
#'@references Lee JJ, Liu DD.A predictive probability design for phase II cancer clinical trialsClinical Trials 2008; 5: 93–106
#'@importFrom VGAM dbetabinom.ab
#'@export
#'
#' @examples
#' \donttest{
#' ## Example 1
#' result<-getPPdesign(N=seq(80,100,5),Delta_T=c(0.95),Delta_L=list(c(0.2,0.2)),
#' Delta_U=list(c(0.8,0.8)),interim_type="Perc",interim=c(0.5,0.7),
#' p0=0.4,pa=0.6,Beta_dis=c(0.4,0.6),Freqprop="pow",nsim=100)
#'                     
#' ## Example 2 
#' result<-getPPdesign(N=seq(80,100,5),Delta_T=c(0.95),Delta_L=0.3,Delta_U=1,
#' interim_type="cont",cohortsize=3,
#' p0=0.4,pa=0.6,Beta_dis=c(0.4,0.6),Freqprop="pow",nsim=100)
#' 
#' ## Example 3 
#' result<-getPPdesign(N=c(68),Delta_T=c(0.95),interim_type="fix",
#' interim=c(30,40),Delta_L=list(c(0.1,0.05),c(0.1,0.1)),Delta_U=list(c(1,1),c(1,1)),
#' p0=0.66,pa=0.8,Beta_dis=c(0.66,0.33),Freqprop="pow",nsim=10)
#' 
#' ## Example 4 
#' result<-getPPdesign(N=c(68),Delta_T=c(0.95),interim_type="fix",
#' interim=c(30,40),Delta_L=list(c(0.1,0.05),c(0.1,0.1)),Delta_U=list(c(1,1),c(1,1)),
#' p0=0.66,pa=0.8,Beta_dis=c(0.66,0.33),Freqprop="typeI",nsim=10)
#' 
#' ## Example 5 
#' result<-getPPdesign(N=c(68),Delta_T=c(0.95),interim_type="fix",
#' interim=c(30),Delta_L=list(c(0.1)),Delta_U=list(c(1,1)),
#' p0=0.66,pa=0.8,Beta_dis=c(0.66,0.33),Freqprop="pow",nsim=10)
#' 
#' result<-getPPdesign(N=c(68),Delta_T=c(0.95),interim_type="fix",
#' interim=c(30),Delta_L=list(c(0.1)),Delta_U=list(c(1,1)),
#' p0=0.66,pa=0.8,Beta_dis=c(0.66,0.33),Freqprop="typeI",nsim=10) 
#' }

getPPdesign<-function(N,Delta_T,interim_type,interim,cohortsize=NULL,interimstart=10,interimstop=5,Delta_L=0,Delta_U=1,p0,pa,Beta_dis,Freqprop,nsim){
   
  if (Freqprop=="pow")   {p_<-pa}
  if (Freqprop=="typeI") {p_<-p0}
  
  # Create dataframe with all parameters + columns for 1) power, 
  #                                                    2) overall stop futility, 
  #                                                    3) stop futility by interim
  
  Boundaries<-cbind(paste0("D_L=",Delta_L,"/D_U=",Delta_U)) 
  
  param<-data.frame(data.table::CJ(N,Delta_T,Boundaries)) 
  param<-cbind(param,pow=0,eff_fin=0,eff_stop=0,fut_stop=0,index=1:dim(param)[1]) 

  # Create index number
  index<-0
  
  for (a in 1:length(N)){

    N_loop<-N[a]
    # Get number of patients for each interim analysis, possibly specific by N parameter (in case of 'Perc' and 'cont')
    if (interim_type=="Perc"){
      interim.f<-floor(interim*N_loop)}
    
    if (interim_type=="cont"){
      interim.f<-seq(interimstart,N_loop-interimstop,cohortsize)}
    
    if (interim_type=="fix"){
      interim.f<-interim}
    
    for (b in 1:length(Delta_T)){
      D_T<-Delta_T[b]
      
      for (c in 1:length(Delta_L)){ 
        D_L<-Delta_L[[c]]
        D_U<-Delta_U[[c]]
        if (length(D_L)==1){D_L<-rep(D_L,length(interim.f))} # If Delta_L is a scalar, change it into a vector
        if (length(D_U)==1){D_U<-rep(D_U,length(interim.f))} # If Delta_U is a scalar, change it into a vector
        
        index<-index+1
        sim<-data.frame(sim=1:nsim,fut=0,fut_int=0,eff=0,eff_int=0,fin=0)
        param[index,"N_interim"]<-length(interim.f)  # Number of interims (without final) for in table
        
        for (d in 1:nsim){
          data<-rbinom(n=N_loop,size=1,prob=p_)
          
          for (e in 1:length(interim.f)){
            
            if (sim[d,]$fut==0 & sim[d,]$eff==0){
              data_int<-data[1:interim.f[e]]
              PP<-Pred_Prob.f(N=N_loop,n1=length(data_int),x1=sum(data_int),beta_par=Beta_dis,p=p0,Delta_T_=D_T)
              
              if (PP<=D_L[e]){ # first check interim futility 
                sim[d,]$fut_int <- e   # futility indicator for specific interim analysis
                sim[d,]$fut     <- 1}  # overall futility indicator
              
              if (PP>=D_U[e]){ # then check interim efficacy
                sim[d,]$eff_int <- e   # futility indicator for specific interim analysis
                sim[d,]$eff     <- 1}  # overall futility indicator
              
            } # 'if' accolade
            
          } #e-loop [number of interim analyses]
          
          if (sim[d,]$fut==0 & sim[d,]$eff==0){ # Only final analysis if no stop at interim analyses
            sim[d,]$fin<-as.numeric(Post_Prob.f(n=N_loop,x=sum(data),beta_par<-Beta_dis,p=p0)>D_T)}
          
        } #d-loop [number of simulations]
          
        param[index,]$pow      <-(sum(sim$fin)+sum(sim$eff))/nsim
        param[index,]$eff_fin  <-sum(sim$fin)/nsim
        param[index,]$fut_stop <-sum(sim$fut)/nsim
        param[index,]$eff_stop <-sum(sim$eff)/nsim
        
        for(i in 1:length(interim.f)){
          param[index,paste0("int_nr",i)] <- interim.f[i]
          param[index,paste0("eff"   ,i)] <- sum(sim$eff_int==i)/nsim
          param[index,paste0("fut"   ,i)] <- sum(sim$fut_int==i)/nsim}
        
      } # c-loop[D_L and D_U: pair of boundaries]    
    } #b-loop [D_T]       
  } #a-loop [N]
  
  # Order dataset, and return results
  firstcols<-c("N","Delta_T","Boundaries","N_interim","pow","eff_fin","eff_stop","fut_stop")
  lastcols<-colnames(param)[(! colnames(param) %in% firstcols) & (! colnames(param) %in% c("index"))]
  
  return(param[,c(firstcols,lastcols)])
} # end function


