% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BEMPO_BEMPP.R
\name{BEMPP}
\alias{BEMPP}
\title{Bayesian efficacy monitoring of 1 binomial parameter}
\usage{
BEMPP(
  N,
  p,
  design,
  interim_type,
  interim,
  cohortsize = NULL,
  interimstart = NULL,
  interimstop = NULL,
  Delta_fut = NULL,
  P_fut,
  Delta_eff = NULL,
  P_eff,
  Delta_fin,
  P_fin,
  Beta_dis,
  nsim
)
}
\arguments{
\item{N}{sample size to investigate operating characteristics: can be a scalar or a vector}

\item{p}{response rate to investigate operating characteristics: can be a scalar or a vector}

\item{design}{choice between "BEMPR" (predictive probability) and "BEMPO" (posterior probability)}

\item{interim_type}{choice between 3 character types: "perc", "cont" or "fix". \cr
Only one sequence at once can be investigated: \cr
1) "perc": percentage(s) of information, e.g. c(0.5,0.75), specified as a vector via the 'interim' parameter \cr
2) "fix": fixed number of patients, e.g. c(10,20), specified as a vector via the 'interim' parameter \cr
3) "cont": continuous monitoring. Run-in without monitoring (e.g. first 10) is specified in the "interimstart" \cr
parameter. Run-out without monitoring (e.g. last 5) is specified in the "interimstop" parameter.}

\item{interim}{This must be a vector, only to be specified when the parameter 'interim type' is "perc" or "fix" \cr
e.g. c(0.5,0.75) when 'interim type' is "perc", or e.g. c(10,20) when 'interim type' is "fix"}

\item{cohortsize}{This must be a scalar, indicating the frequency of continuous monitoring. \cr
Only needs to be defined in case the parameter 'interim_type'="cont"}

\item{interimstart}{This must be a scalar, specifying run-in without monitoring. \cr
Only needs to be defined in case the parameter 'interim_type'="cont". Default is 10}

\item{interimstop}{This must be a scalar, specifying run-out without monitoring. \cr
Only needs to be defined in case the parameter 'interim_type'="cont". Default is 5.}

\item{Delta_fut}{This may be a scalar or a vector. Must only be specified for BEMPO \cr}

\item{P_fut}{This may be a scalar or a vector \cr}

\item{Delta_eff}{This may be a scalar or a vector. Must only be specified for BEMPO \cr}

\item{P_eff}{This may be a scalar or a vector \cr}

\item{Delta_fin}{This must be a scalar\cr}

\item{P_fin}{This must be a scalar\cr}

\item{Beta_dis}{two parameters Beta(alpha,beta) of the prior Beta distribution}

\item{nsim}{number of simulations}
}
\value{
a list of 3 data.frames: first with design parameters + operating characteristics ($param_simul) , \cr
one with futility decision rules ($Fut_rules), and one with efficacy decision rules ($Eff_rules)
\itemize{
\item param_simul: dataframe with all input parameters + operating characteristics
\item param_simul$scenario: combination of vector N and p
\item param_simul$N: sample size for which frequentist properties simulated
\item param_simul$p: binomial parameter (true proportion) for which frequentist properties simulated
\item param_simul$param: string with all input parameters
\item param_simul$N_interim: number of interim analyses (final analysis not included)
\item param_simul$interim: string with number of patients per interim analysis
\item param_simul$pow: proportion of simulations where \eqn{H0} was rejected
\item param_simul$eff_fin: proportion of simulations where \eqn{H0} was rejected only at the final analysis
\item param_simul$eff_stop: proportion of simulations where \eqn{H0} was rejected at an interim analysis
\item param_simul$fut: proportion of simulations where \eqn{H0} was not rejected
\item param_simul$eff_fin: proportion of simulations where decision of not rejecting \eqn{H0} was at the final analysis
\item param_simul$fut_stop: proportion of simulations where trial stopped for futility at an interim analysis
\item param_simul$N_avg: average number of patients for simulated datasets
\item param_simul$int_nrx: number of patients at x'th interim analysis
\item param_simul$effx: proportion of simulations where trial stopped for efficacy at xth interim analysis
\item param_simul$futx: proportion of simulations where trial stopped for futility at xth interim analysis
\item param_simul$decision rules: for efficacy: if >=x: stop for efficacy; for futility: if <=x: stop for futility
\item Fut_rules: dataframe with futility rules for each interim
\item Fut_rules$n: number of patients at interim
\item Fut_rules$x: number of successes at interim (Futility if observed successes <=x)
\item Fut_rules$P_Bayes/Fut_rules$PP: actual posterior probability or PP at cutoff
\item Eff_rules: dataframe with efficacy rules for each interim
\item Eff_rules$n: number of patients at interim
\item Eff_rules$x: number of successes at interim (Efficacy if observed successes >=x)
\item Eff_rules$P_Bayes/Fut_rules$PP: actual posterior probability or PP at cutoff
}
}
\description{
Single arm Bayesian Efficacy Monitoring Via Predictive Probability (BEMPR) or Bayesian Efficacy Monitoring Via Posterior Probability (BEMPO) \cr
The frequentist properties can be investigated for combinations of \cr
1) Different sample sizes
2) Different true response rates
#
BEMPO interim futility stopping rule is P(Delta<=Delta_fut)>P_fut
BEMPO interim efficacy stopping rule is P(Delta>Delta_eff)>=P_eff
BEMPO final efficacy stopping rule is   P(Delta>Delta_fin)>=P_fin

BEMPR interim futility stopping rule is PP<P_fut
BEMPR interim efficacy stopping rule is PP>=P_eff
BEMPR final efficacy stopping rule is P(Delta>Delta_fin)>=P_fin
}
\examples{
#Check versus https://biostatistics.mdanderson.org/shinyapps/BEMPO/
test1<-BEMPP(N=15,p=0.5, design="BEMPO", interim_type="fix",interim=c(5,10),Delta_fut=0.3,
 P_fut=0.7,Delta_eff=0.3,P_eff=0.9,Delta_fin=0.3,P_fin=0.8,Beta_dis=c(0.5,0.5),nsim=10)
test2<-BEMPP(N=15,p=0.5, design="BEMPO", interim_type="fix",interim=c(5,10),Delta_fut=0.3,
 P_fut=1  ,Delta_eff=0.3,P_eff=1  ,Delta_fin=0.3,P_fin=0.8,Beta_dis=c(0.5,0.5),nsim=10)
test3<-BEMPP(N=15,p=0.5, design="BEMPR", interim_type="fix",interim=c(5,10),P_fut=0.3,
 P_eff=0.9,Delta_fin=0.3,P_fin=0.7,Beta_dis=c(0.5,0.5),nsim=10)
test4<-BEMPP(N=15,p=0.5, design="BEMPR", interim_type="fix",interim=c(5,10),P_fut=0,P_eff=1,
 Delta_fin=0.3,P_fin=0.8,Beta_dis=c(0.5,0.5),nsim=10)
}
\references{
Lee JJ, Liu DD.A predictive probability design for phase II cancer clinical trialsClinical Trials 2008; 5: 93–106
}
