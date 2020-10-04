
#' @title Plot properties of continuous safety monitoring boundaries
#'
#' @description This function plots the properties of continuous safety monitoring boundaries
#'
#' @details Freq_prob_date gives the P(X>c+1|n,tox_true) for different n and for different assumptions,
#'     specified by the 'tox_ass' parameter
#' @details Bayes_probgives the Bayesian posterior probability P(toxicity rate<=target|c+1,n), so given that
#'     this boundary 'c+1' is exceeded, for different assumptions, specified by the 'tox_ass' parameter
#' @details Cum_Freq_prob_1 gives the cumulative probability of at least 1x crossing the boundary
#' @details Cum_Freq_prob_2 gives the cumulative probability of at least 1x/2x/3x and 4x crossing the boundary
#'
#' @param data_in resulting list (containing data with results, and vector with assumptions for calculating
#'     properties) from the 'getbound' function
#' @param out path for storing png output files
#'
#' @examples
#'\dontrun{
#' result<-getbound(SS=39,target=0.2,P_target=0.9,maxtox=NULL,tox_ass=seq(0.10,0.5,0.05),
#'     sim=1,nsim=100,prt1=6,out=NULL)
#' plotbound(data_in=result,out="C:/Users/kdhollander/Downloads/")
#' }
#'
#' @export
#' @importFrom grDevices dev.off png
#' @importFrom graphics grid legend lines par plot

plotbound<-function(data_in,out){

  data    <- data_in[[1]]
  tox_ass <- data_in[[2]]

  SS      <- dim(data)[1]
  n_ass   <- length(tox_ass)

  if (n_ass>=3) {cols <- RColorBrewer::brewer.pal(n_ass,"Spectral")}
  if (n_ass==2) {cols <- c("springgreen4","steelblue4")}
  if (n_ass==1) {cols <- "springgreen4"}


  # Plot 1: P(X>c|n,true toxicity rate), so the probability of crossing the boundary at each n
  #-------------------------------------------------------------------------------------------

  png(filename=paste0(out,"Freq_prob_",Sys.Date(),".png"),width=600,height=400)

  loc<-dim(data)[2]-n_ass*5+1
  plot(1:SS,data[,loc], type="l",xlab="Number of patients treated",
       ylab="Frequentist P[X>=c+1|n,true toxicity rate]",ylim=c(0,1),col=cols[1],lwd=2)
  if (length(tox_ass)>1){
  for (i in 1:(n_ass-1)){
    lines(1:SS,data[,loc+i],col=cols[i+1],lwd=2)
  }}
  grid()
  legend(x="bottomleft",lwd=2,bty="n",paste0("True toxicity rate = ",tox_ass),col=cols)
  dev.off()


  # Plot 2: Bayesian posterior probability P(toxicity<=given rate|c+1,n), so at decision to stop, for set of assumptions
  #--------------------------------------------------------------------------------------------------------------------

  png(filename=paste0(out,"Bayes_prob_",Sys.Date(),".png"),width=600,height=400)

  loc<-7
  plot(1:SS,data[,loc], type="l",xlab="Number of patients treated",
       ylab="Bayesian posterior P(rate<=true toxicity rate|c+1,n)",ylim=c(0,1),col=cols[1],lwd=2)
  if (length(tox_ass)>1){
  for (i in 1:(n_ass-1)){
    lines(1:SS,data[,loc+i],col=cols[i+1],lwd=2)
  }}
  grid()
  legend(x="bottomleft",lwd=2,bty="n",paste0("True toxicity rate = ",tox_ass),col=cols)
  dev.off()

  # Plot 3: P(X>c at least once|n), so the probability of crossing the boundary at least once by n patients
  #--------------------------------------------------------------------------------------------------------

  png(filename=paste0(out,"Cum_Freq_prob_1_",Sys.Date(),".png"),width=600,height=400)

  loc<-dim(data)[2]-n_ass*4+1
  plot(1:SS,data[,loc], type="l",xlab="Number of patients treated",
       ylab="Frequentist Cum P[once X>=c+1|n,true toxicity rate]",ylim=c(0,1),col=cols[1],lwd=2)
  if (length(tox_ass)>1){
    for (i in 1:(n_ass-1)){
      lines(1:SS,data[,loc+i],col=cols[i+1],lwd=2)
    }}
  grid()
  legend(x="bottomleft",lwd=2,bty="n",paste0("True toxicity rate = ",tox_ass),col=cols)
  dev.off()


  # Plot 4: P(X>c at least once|n), so the probability of crossing the boundary at least once, 2x, 3x and 4x by n patients
  #-----------------------------------------------------------------------------------------------------------------------

  png(filename=paste0(out,"Cum_Freq_prob_2_",Sys.Date(),".png"),width=1200,height=800)

  par(mfrow=c(2,2))

  loc<-dim(data)[2]-n_ass*4+1
  plot(1:SS,data[,loc], type="l",xlab="Number of patients treated",
       ylab="Frequentist Cum P[once X>=c+1|n,true toxicity rate]",ylim=c(0,1),col=cols[1],lwd=2)
  if (length(tox_ass)>1){
    for (i in 1:(n_ass-1)){
      lines(1:SS,data[,loc+i],col=cols[i+1],lwd=2)
    }}
  grid()
  legend(x="bottomleft",lwd=2,bty="n",paste0("True toxicity rate = ",tox_ass),col=cols)

  loc<-dim(data)[2]-n_ass*3+1
  plot(1:SS,data[,loc], type="l",xlab="Number of patients treated",
       ylab="Frequentist Cum P[twice X>=c+1|n,true toxicity rate]",ylim=c(0,1),col=cols[1],lwd=2)
  if (length(tox_ass)>1){
    for (i in 1:(n_ass-1)){
      lines(1:SS,data[,loc+i],col=cols[i+1],lwd=2)
    }}
  grid()
  legend(x="bottomleft",lwd=2,bty="n",paste0("True toxicity rate = ",tox_ass),col=cols)

  loc<-dim(data)[2]-n_ass*2+1
  plot(1:SS,data[,loc], type="l",xlab="Number of patients treated",
       ylab="Frequentist Cum P[3 times X>=c+1|n,true toxicity rate]",ylim=c(0,1),col=cols[1],lwd=2)
  if (length(tox_ass)>1){
    for (i in 1:(n_ass-1)){
      lines(1:SS,data[,loc+i],col=cols[i+1],lwd=2)
    }}
  grid()
  legend(x="bottomleft",lwd=2,bty="n",paste0("True toxicity rate = ",tox_ass),col=cols)

  loc<-dim(data)[2]-n_ass+1
  plot(1:SS,data[,loc], type="l",xlab="Number of patients treated",
       ylab="Frequentist Cum P[4 times X>=c+1|n,true toxicity rate]",ylim=c(0,1),col=cols[1],lwd=2)
  if (length(tox_ass)>1){
    for (i in 1:(n_ass-1)){
      lines(1:SS,data[,loc+i],col=cols[i+1],lwd=2)
    }}
  grid()
  legend(x="bottomleft",lwd=2,bty="n",paste0("True toxicity rate = ",tox_ass),col=cols)
  dev.off()
}


