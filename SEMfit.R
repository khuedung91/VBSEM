############ R script: SEMfit ############

# For fitting a simple structural equation model
# via mean field variational Bayes (MFVB) and
# Markov chain Monte Carlo (MCMC).
# The script contains an option for performing 
# bootstrap and improving the MFVB results.
# The model is fitted to the Holzinger and 
# Swineford (1939) data with a single latent
# factor (visual ability).

# Last changed: 26 NOV 2021.


# Set tasks:

doMFVB <- TRUE
doMCMC <- TRUE
doBootstrap <- TRUE
makePlots <- TRUE # set this to TRUE only if the previous three flags are TRUE

# Load MFVB function:

source("MFVBforSEM.R")

# Load required packages:

library(lavaan)
library(invgamma)
library(rstan)
library(coda)
#rstan_options(auto_write=TRUE)
#options(mc.cores=parallel::detectCores())

# Extract data:

HS1939 <- HolzingerSwineford1939
#head(HS1939)
SEMdata <- HS1939[,c('x1','x2','x3')]
y <- as.matrix(SEMdata)
m <- ncol(y)   ;   n <- nrow(y)

# Set hyperparameters:

mu.lambda <- 0   ;   sigsq.lambda <- 1
sigsq.nu <- 1e2
delta.psi <- rep(0.01,m)
delta.sigsq <- 1
kappa.psi <- rep(1,m)   ;   kappa.sigsq <- 1

##########################################
################## MFVB ##################
##########################################

if (doMFVB)
{
  # Set MFVB convergence tolerance and maximum number of iterations:
  
  tolVal <- 1e-2  
  maxIter <- 1000
  
  # Run MFVB function:
  
  MFVBoutput <- MFVBforSEM(y,n,m,mu.lambda,sigsq.lambda,sigsq.nu,delta.psi,
                           delta.sigsq,tolVal,maxIter,convChk=TRUE)
  
  # Extract MFVB results:
  
  mu.q.nu.MFVB <- MFVBoutput$"mu.q.nu.MFVB"
  sigsq.q.nu.MFVB <- MFVBoutput$"sigsq.q.nu.MFVB"
  mu.q.lambda.MFVB <- MFVBoutput$"mu.q.lambda.MFVB"
  sigsq.q.lambda.MFVB <- MFVBoutput$"sigsq.q.lambda.MFVB"
  kappa.q.sigsq.MFVB <-   MFVBoutput$"kappa.q.sigsq.MFVB"
  delta.q.sigsq.MFVB <- MFVBoutput$"delta.q.sigsq.MFVB"
  kappa.q.psi.MFVB <- MFVBoutput$"kappa.q.psi.MFVB"
  delta.q.psi.MFVB <- MFVBoutput$"delta.q.psi.MFVB"
  MFVBtime <- MFVBoutput$"MFVBtime" 

  # Print MFVB computational time:
  
  cat("MFVB converged in",MFVBtime,"seconds.","\n")
  
  # Obtain MFVB credible intervals:
  
  MFVBestVec <- c(mu.q.nu.MFVB,mu.q.lambda.MFVB[2:m],
                  delta.q.sigsq.MFVB/(kappa.q.sigsq.MFVB-2),
                  delta.q.psi.MFVB/(kappa.q.psi.MFVB-2))
  
  MFVB_CI <- array(NA,dim=c(3,9))
  MFVB_CI[2,] <- MFVBestVec
  MFVB_CI[1,1:3] <- mu.q.nu.MFVB-qnorm(0.975)*sqrt(sigsq.q.nu.MFVB)
  MFVB_CI[3,1:3] <- mu.q.nu.MFVB+qnorm(0.975)*sqrt(sigsq.q.nu.MFVB)
  MFVB_CI[1,4:5] <- mu.q.lambda.MFVB[2:m]-qnorm(0.975)*sqrt(sigsq.q.lambda.MFVB[2:m])
  MFVB_CI[3,4:5] <- mu.q.lambda.MFVB[2:m]+qnorm(0.975)*sqrt(sigsq.q.lambda.MFVB[2:m])
  
  # N.B.: the 'rate' parameter of the 'invgamma' package is actually a scale parameter. 
  MFVB_CI[c(1,3),6] <- qinvgamma(c(0.025,0.975),shape=0.5*kappa.q.sigsq.MFVB,
                                 rate=0.5*delta.q.sigsq.MFVB)
  MFVB_CI[c(1,3),7] <- qinvgamma(c(0.025,0.975),shape=0.5*kappa.q.psi.MFVB[1],
                                 rate=0.5*delta.q.psi.MFVB[1])
  MFVB_CI[c(1,3),8] <- qinvgamma(c(0.025,0.975),shape=0.5*kappa.q.psi.MFVB[2],
                                 rate=0.5*delta.q.psi.MFVB[2])
  MFVB_CI[c(1,3),9] <- qinvgamma(c(0.025,0.975),shape=0.5*kappa.q.psi.MFVB[3],
                                 rate=0.5*delta.q.psi.MFVB[3])
}

##########################################
################## MCMC ##################
##########################################

if(doMCMC)
{
  # Perform MCMC via rstan code:
  
  MCMCdata <- list(N=n,K=m,y=y,mu_lambda=mu.lambda,sig2_lambda=sigsq.lambda,
                      delta_psi=delta.psi,delta_sig2=delta.sigsq,sigma_nu=sqrt(sigsq.nu)) 
  
  HMCresults <- stan(file='MCMCforSEM.stan',data=MCMCdata,iter=15000,warmup=7500,
                     chains=1,algorithm='NUTS',control=list(adapt_delta=0.85))
  
  # Extract MCMC results:
  
  nuMCMCmat <- as.matrix(As.mcmc.list(HMCresults,pars="nu")[[1]])
  lambdaMCMCmat <- as.matrix(As.mcmc.list(HMCresults,pars="lambday")[[1]])
  sigsqMCMCmat <- as.matrix(As.mcmc.list(HMCresults,pars="sigma2")[[1]])
  psiMCMCmat <- as.matrix(As.mcmc.list(HMCresults,pars="psidiag")[[1]])
  
  # Obtain MCMC credible intervals:
  
  HMC_CI <- apply(cbind(nuMCMCmat,lambdaMCMCmat,sigsqMCMCmat,psiMCMCmat),
                  2,quantile,prob = c(0.025,0.5,0.975))
}

##########################################
############ MFVB + bootstrap ############
##########################################

if(doBootstrap)
{
  # Select number of bootstrap iterations:
  
  nBootstrap <- 1000
  
  bootMat <- array(NA,dim=c(nBootstrap,20))
  SEMdataMat <- as.matrix(SEMdata)
  for(bs in 1:nBootstrap)
  {
    # Create new bootstrap dataset by sampling with replacement:
    
    newidx <- sample(n,n,replace=T)
    yBoot <- SEMdataMat[newidx,]
    
    # Run MFVB on b-th bootstrap dataset:
    
    MFVBoutput <- MFVBforSEM(yBoot,n,m,mu.lambda,sigsq.lambda,sigsq.nu,delta.psi,
                             delta.sigsq,tolVal,maxIter,convChk=FALSE)
    
    # Extract MFVB results:
    
    mu.q.nu.MFVB.boot <- MFVBoutput$"mu.q.nu.MFVB"
    sigsq.q.nu.MFVB.boot <- MFVBoutput$"sigsq.q.nu.MFVB"
    mu.q.lambda.MFVB.boot <- MFVBoutput$"mu.q.lambda.MFVB"
    sigsq.q.lambda.MFVB.boot <- MFVBoutput$"sigsq.q.lambda.MFVB"
    kappa.q.sigsq.MFVB.boot <-   MFVBoutput$"kappa.q.sigsq.MFVB"
    delta.q.sigsq.MFVB.boot <- MFVBoutput$"delta.q.sigsq.MFVB"
    kappa.q.psi.MFVB.boot <- MFVBoutput$"kappa.q.psi.MFVB"
    delta.q.psi.MFVB.boot <- MFVBoutput$"delta.q.psi.MFVB"
    MFVBtime.boot <- MFVBoutput$"MFVBtime" 
    
    bootMat[bs,] <- c(mu.q.nu.MFVB.boot,sigsq.q.nu.MFVB.boot,mu.q.lambda.MFVB.boot,
                      sigsq.q.lambda.MFVB.boot,kappa.q.psi.MFVB.boot,delta.q.psi.MFVB.boot,
                      kappa.q.sigsq.MFVB.boot,delta.q.sigsq.MFVB.boot)
    
    # Print MFVB computational time and bootstrap iteration number:
    
    cat("MFVB bootstrap iteration",bs,"converged in",MFVBtime.boot,"seconds.","\n")
  }
  
  # Extract bootstrap results:
  
  nu.MFVB.boot <- bootMat[,1:3]
  lambda.MFVB.boot <- bootMat[,8:9] # not saving lambda1
  sigsq.MFVB.boot <- bootMat[,20]/(bootMat[,19]-2)
  psi.MFVB.boot <- bootMat[,16:18]/(bootMat[,13:15]-2)
  
  # Obtain MFVB percentile and pivotal credible intervals:
  
  bootMatPointEst <- cbind(nu.MFVB.boot,lambda.MFVB.boot,sigsq.MFVB.boot,psi.MFVB.boot)
  
  deltaMat <- sweep(bootMatPointEst,2,MFVBestVec)
  quantileMat <- apply(deltaMat,2,quantile,prob=c(0.025,0.975))
  percentileCI <- rbind(MFVBestVec+ quantileMat[1,],MFVBestVec,MFVBestVec+ quantileMat[2,])
  var.q.sigsq.MFVB <- 2*delta.q.sigsq.MFVB^2/((kappa.q.sigsq.MFVB-2)^2*(kappa.q.sigsq.MFVB-4))
  var.q.psi.MFVB <- 2*delta.q.psi.MFVB^2/((kappa.q.psi.MFVB-2)^2*(kappa.q.psi.MFVB-4))
  
  MFVBvarEst <- c(sigsq.q.nu.MFVB, sigsq.q.lambda.MFVB[2:m],var.q.sigsq.MFVB,var.q.psi.MFVB)
  
  var.q.sigsq.boot <- 2*bootMat[,20]^2/((bootMat[,19]-2)^2*(bootMat[,19]-4))
  var.q.psi.boot <- 2*bootMat[,16:18]^2/((bootMat[,13:15]-2)^2*(bootMat[,13:15]-4))
  
  MFVBvarEstBoot <- cbind(bootMat[,4:6],bootMat[,11:12],var.q.sigsq.boot,var.q.psi.boot)
  
  tMat <- deltaMat/sqrt(MFVBvarEstBoot)
  quantileMat_pivot <- apply(abs(tMat),2,quantile,prob=c(0.025,0.975),na.rm=T)
  
  pivotalCI <- rbind(MFVBestVec-sqrt(MFVBvarEst)*quantileMat_pivot[2,],MFVBestVec,
                     MFVBestVec+sqrt(MFVBvarEst)*quantileMat_pivot[2,])
}

##########################################
################# Plots ##################
##########################################

if(makePlots){
 
  # Set graphical parameters:
  
  MFVBcol <- "black"
  MFVBlwd <- 3
  MFVBpch <- 16
  MFVBcex <- 1.5
  
  MCMCcol <- "grey65"
  MCMClwd <- 7
  MCMCpch <- 16
  MCMCcex <- 2.5
  
  # Define plots main expressions:
  
  main_expressions_nu <- 
    sapply(1:3, function(i) {
      as.expression(substitute(A[B],list(A=as.name("nu"),as.name,B=(i))))
    })
  main_expressions_lambda <- 
    sapply(2:3, function(i) {
      as.expression(substitute(A[B],list(A=as.name("lambda"),as.name,B=(i))))
    })
  main_expressions_psi <- 
    sapply(1:3, function(i) {
      as.expression(substitute(A[B],list(A=as.name("psi"),as.name,B=(i))))
    })
  main_expressions <- c(main_expressions_nu,main_expressions_lambda,
                        expression(sigma^2),main_expressions_psi)
  
  # Draw plots:
  
  par(mfrow=c(3,3),mar=c(2.5,2,2.5,1))
  for(parnum in 1:9){
    xlimValInf <- min(c(percentileCI[1,parnum],pivotalCI[1,parnum],HMC_CI[1,parnum]))
    xlimValSup <- max(c(percentileCI[3,parnum],pivotalCI[3,parnum],HMC_CI[3,parnum]))
    xlimVal <- c(xlimValInf-0.1*(xlimValSup-xlimValInf),xlimValSup)
    
    xlimRange <- xlimValSup - xlimValInf
    ylimValSup <- 1
    ylimValInf <- ylimValSup/3
    ylimVal <- c(ylimValInf,ylimValSup)
    
    if(parnum %in% c(1,4,7)){
      plot(0,type="l",bty="n",ylim=ylimVal,
           xlim=c(xlimVal[1]-0.3*(xlimVal[2]-xlimVal[1]),xlimVal[2]),
           xlab="",yaxt="n",ylab='',main=main_expressions[parnum])
    }else{
      plot(0,type="l",bty="n",ylim=ylimVal,xlim=xlimVal,
           xlab="",yaxt="n",ylab='',main=main_expressions[parnum])
    }
    yVal <- ylimValInf
    lines(c(HMC_CI[1,parnum],HMC_CI[3,parnum]),rep(yVal,2),col=MCMCcol,lwd=MCMClwd)
    points(HMC_CI[2,parnum],yVal,col=MCMCcol,pch=MCMCpch,cex=MCMCcex)
    lines(c(pivotalCI[1,parnum],pivotalCI[3,parnum]),rep(yVal,2),col=MFVBcol,lwd=MFVBlwd)
    points(pivotalCI[2,parnum],yVal,col=MFVBcol,pch=MFVBpch,cex=MFVBcex)
    if(parnum %in% c(1,4,7)){
      text(x=xlimVal[1]-0.15*(xlimVal[2]-xlimVal[1]),y=mean(yVal),
           labels="MFVB piv.boot",cex=0.9)
    }
    
    yVal <- ylimValInf*2
    lines(c(HMC_CI[1,parnum],HMC_CI[3,parnum]),rep(yVal,2),col=MCMCcol,lwd=MCMClwd)
    points(HMC_CI[2,parnum],yVal,col=MCMCcol,pch=MCMCpch,cex=MCMCcex)
    lines(c(percentileCI[1,parnum],percentileCI[3,parnum]),rep(yVal,2),col=MFVBcol,lwd=MFVBlwd)
    points(percentileCI[2,parnum],yVal,col=MFVBcol,pch=MFVBpch,cex=MFVBcex)
    if(parnum %in% c(1,4,7)){
      text(x=xlimVal[1]-0.15*(xlimVal[2]-xlimVal[1]),y=mean(yVal),
           labels="MFVB per.boot",cex=0.9)
    }
    
    yVal <- ylimValInf*3
    lines(c(HMC_CI[1,parnum],HMC_CI[3,parnum]),rep(yVal,2),col=MCMCcol,lwd=MCMClwd)
    points(HMC_CI[2,parnum],yVal,col=MCMCcol,pch=MCMCpch,cex=MCMCcex)
    lines(c(MFVB_CI[1,parnum],MFVB_CI[3,parnum]),rep(yVal,2),col=MFVBcol,lwd=MFVBlwd)
    points(MFVB_CI[2,parnum],yVal,col=MFVBcol,pch=MFVBpch,cex=MFVBcex)
    if(parnum %in% c(1,4,7)){
      text(x=xlimVal[1]-0.15*(xlimVal[2]-xlimVal[1]),y=mean(yVal),labels="MFVB",cex=0.9)
    }
  }
}

############ End of SEMfit ############