########## R function: MFVBforSEM ###########

# For fitting a simple structural equation model
# via mean field variational Bayes.

# Last changed: 25 NOV 2022.


MFVBforSEM <- function(y,n,m,mu.lambda=0,sigsq.lambda=1,sigsq.nu=1e2,
                       delta.psi=rep(0.01,m),kappa.psi=rep(1,m),delta.sigsq=1,
                       kappa.sigsq=1,tolVal=1e-2,maxIter=1000,convChk=TRUE)
{
  # Initialize MFVB paramters to optimize:
  
  mu.q.nu <- rep(0,m)   ;   mu.q.nu2 <- rep(0,m)   ;   sigsq.q.nu <- rep(1,m)
  mu.q.lambda <- mu.q.lambda2 <- sigsq.q.lambda <- rep(1,m)
  sigsq.q.lambda[1] <- 1
  mu.q.recip.psi <- rep(1,m)
  mu.q.eta <- mu.q.eta2 <- rep(1,n)
  mu.q.recip.sigsq <- 1
  
  # Set MFVB fixed parameters:
  
  delta.q.psi <- kappa.q.psi <- rep(1,m)
  kappa.q.sigsq <- n + kappa.sigsq
  kappa.q.psi <- n + kappa.psi + 1 
  
  # Perform MFVB iterations:
  
  converged <- FALSE  
  itNum <- 0
  MFVBstartTime <- Sys.time()
  while (!converged)
  {
    itNum <- itNum + 1 
    
    # Do updates for q*(psi), q*(nu) and q*(lambda):
    
    for (j in 1:m) 
    {
      ySubSet = y[,j]
      
      delta.q.psi[j] <- (sum(ySubSet*ySubSet + mu.q.nu2[j] 
                             + mu.q.lambda2[j]*mu.q.eta2 - 2*ySubSet*mu.q.nu[j] 
                             - 2*ySubSet*mu.q.lambda[j]*mu.q.eta 
                             + 2*mu.q.nu[j]*mu.q.lambda[j]*mu.q.eta)
                         + (mu.q.lambda2[j] - 2*mu.lambda*mu.q.lambda[j] 
                            + mu.lambda^2)/sigsq.lambda + delta.psi[j])
      mu.q.recip.psi[j] <- kappa.q.psi[j]/delta.q.psi[j] 
      
      sigsq.q.nu[j] <- 1/(n*mu.q.recip.psi[j] + 1/sigsq.nu)
      mu.q.nu[j] <- sigsq.q.nu[j]*mu.q.recip.psi[j]*sum(ySubSet - mu.q.lambda[j]*mu.q.eta)
      mu.q.nu2[j] <- sigsq.q.nu[j] + mu.q.nu[j]^2
      
      if (j > 1)
      {
        sigsq.q.lambda[j] <- 1/(mu.q.recip.psi[j]*(sum(mu.q.eta2) + 1/sigsq.lambda))
        mu.q.lambda[j] <- sigsq.q.lambda[j]*(sum(mu.q.eta*(ySubSet - mu.q.nu[j])) 
                                             + mu.lambda/sigsq.lambda )*mu.q.recip.psi[j]
        mu.q.lambda2[j] <- sigsq.q.lambda[j] + mu.q.lambda[j]^2
      }
    }
    
    # Do updates for q*(eta):
    
    sigsq.q.eta <- rep((1/(sum(mu.q.recip.psi*mu.q.lambda2) + mu.q.recip.sigsq)),n)
    for (i in 1:n) mu.q.eta[i] <- sigsq.q.eta[i]*sum(mu.q.recip.psi*mu.q.lambda*(y[i,] - mu.q.nu))
    mu.q.eta2 <- sigsq.q.eta + mu.q.eta^2
    
    # Do updates for q*(sigma^2):
    
    delta.q.sigsq <- sum(mu.q.eta2) + delta.sigsq
    mu.q.recip.sigsq <- kappa.q.sigsq/delta.q.sigsq
    
    # Check for convergence:
    
    parVecHat <- c(mu.q.nu,mu.q.recip.psi,mu.q.lambda,delta.psi,
                   mu.q.eta,mu.q.recip.sigsq,delta.sigsq)
    
    if (itNum>1)
    {
      relErrTemp <- abs((parVecHat - parVecHatPrev)/parVecHatPrev)
      relErrTemp[is.nan(relErrTemp)] <- 0
      relErr <- max(relErrTemp)
      if(convChk==TRUE) cat("Iteration number:",itNum,";","MFVB relative error:",signif(relErr,5),"\n")
      if (relErr<tolVal)
        converged <- TRUE
    }
    if (itNum>=maxIter)
    {
      converged <- TRUE
      warning("Maximum number of MFVB iterations exceeded.\n")
    }
    
    parVecHatPrev <- parVecHat
  }
  
  MFVBendTime <- Sys.time()
  MFVBtime <- MFVBendTime - MFVBstartTime
  
  # Return results:

  MFVBoutput <- list()
  MFVBoutput$"mu.q.nu.MFVB" <- mu.q.nu
  MFVBoutput$"sigsq.q.nu.MFVB" <- sigsq.q.nu
  MFVBoutput$"mu.q.lambda.MFVB" <- mu.q.lambda
  MFVBoutput$"sigsq.q.lambda.MFVB" <- sigsq.q.lambda
  MFVBoutput$"kappa.q.sigsq.MFVB" <- kappa.q.sigsq
  MFVBoutput$"delta.q.sigsq.MFVB" <- delta.q.sigsq
  MFVBoutput$"kappa.q.psi.MFVB" <- kappa.q.psi
  MFVBoutput$"delta.q.psi.MFVB" <- delta.q.psi
  MFVBoutput$"MFVBtime" <- MFVBtime
  
  return(MFVBoutput)
}

############ End of MFVBforSEM ############
