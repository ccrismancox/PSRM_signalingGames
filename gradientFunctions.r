#### Helper functions ####
dbivnorm <- function(x1, x2, rho = 0){
  denom <- 1 - rho^2
  pdf <- 1/(2*pi * sqrt(denom)) * exp( -1/(2*denom)  * (x1^2 + x2^2 - 2 *rho * x1 *x2))
  return(pdf)
}

delPbivnorm <- function(x1, x2, rho=0){
  ## The partial derivative is from Wickens (1992) for 
  ## $\pnorm_2(x, y, rho)$ we have
  ## $Dx = \pnorm(x) * \pnorm((y-x*rho)/sqrt(1-rho^2))$
  ## This function always just put the one that is w.r.t.
  ## as the first argument.
  return(dnorm(x1)*pnorm((x2-x1*rho)/sqrt(1-rho^2)))
}


#### gradient for the PL method ####
eval_gr_qll <- function(x,PRhat,PFhat,Y,regr){
  # here x = (theta,p)
  ##EQ[,1]= pR
  ##EQ[,2]= pC
  ##EQ[,3]= pF
  M <- dim(Y)[2]
  
  
  param <-vec2U.regr(x,regr)
  
  
  c <- cStar.jo(PRhat,param)
  pR <- f.jo(PFhat, param)
  pR[pR<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps) 
  pC <- g.jo(c, param)
  pC[pC<=(.Machine$double.eps)^(.25)] <- (.Machine$double.eps)^(.25) #edited 7 March
  pC[pC>=1-(.Machine$double.eps)^(.25)] <- 1-(.Machine$double.eps)^(.25) #edited 7 March
  pF <- h.jo(c, param)/pC
  v1 <- (c-param$barWA)/param$sig
  v2 <- (c-param$bara)/param$sig
  d1 <- (param$barWA - param$bara)/(param$sig*sqrt(2))
  d2 <- (param$barWA - c)/(param$sig)
  
  r <- 1/sqrt(2)
  VApr <- (param$VA/pR - param$VA * (pR-1)/(pR^2))
  P1 <- pnorm(v1)
  P1[P1<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  P1[P1>=1-sqrt(.Machine$double.eps)] <- 1-sqrt(.Machine$double.eps)
  P2 <- pnorm(v2)
  P2[P2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  P2[P2>=1-sqrt(.Machine$double.eps)] <- 1-sqrt(.Machine$double.eps)
  D1 <- dnorm(v1)
  D2 <- dnorm(v2)
  P1P2 <-  P1*P2 #more stable than pbiv w/rho=0
  P1P2[P1P2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)

  Del_d1v1 <- delPbivnorm(d1, -v1, r)
  Del_v1d1 <- delPbivnorm(-v1, d1, r)
  PBdv <- pbivnorm(d1, -v1, rho=r)
  PBdv[PBdv<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps) 
  PBdv[PBdv>= (1-sqrt(.Machine$double.eps))] <- 1-sqrt(.Machine$double.eps) 
  WBinner <- as.numeric((-param$CB + param$VB*(-PFhat + 1) + param$barWB*PFhat)/PFhat)
  
  
  dWA4.denom <- (pC*(1 - PBdv/pC)*pnorm(WBinner))  
  dWA4.denom[dWA4.denom<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps) 
  dWA <- rbind( -D1/P1,
                P2*D1/pC,
                (sqrt(2)*Del_d1v1/2 +Del_v1d1)/PBdv,
                (pC*((-sqrt(2)*Del_d1v1/2 - Del_v1d1)/pC + P2*PBdv*D1/pC**2)*pnorm(WBinner) + (1 - PBdv/pC)*P2*pnorm(WBinner)*D1)/dWA4.denom 
  )
  
  dWA  <- apply(regr$barWA, 2, function(x){t(as.numeric(x)*t(dWA))*as.numeric(Y)})
  
  dWB.denom <- (-pnorm(WBinner) + 1)
  dWB.denom[dWB.denom<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  dWB.denom2 <- pnorm(WBinner)
  dWB.denom2[dWB.denom2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  dWB <- rbind(0,
               (-dnorm(WBinner)/dWB.denom),
               ( dnorm(WBinner)/dWB.denom2),
               ( dnorm(WBinner)/dWB.denom2)
  ) 
  
  dWB <- apply(regr$barWB, 2, function(x){t(as.numeric(x)*t(dWB))*as.numeric(Y)})
  
  
  dbara.denom <- (pC*(-pF + 1)*pnorm(WBinner))
  dbara.denom[dbara.denom<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  dbara.stuff <- (-2*P1*P2 + 2)
  dbara.stuff[dbara.stuff<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  dbara <- rbind( -D2/P2,
                  P1*D2/pC,
                  - sqrt(2)*Del_d1v1/(2*PBdv),
                  ((pC*(sqrt(2)*Del_d1v1/dbara.stuff + pF*P1*D2/pC)*pnorm(WBinner) + (-pF + 1)*P1*pnorm(WBinner)*D2)/(dbara.denom))
                  
  )
  dbara <- apply(regr$bara, 2, function(x){t(as.numeric(x)*t(dbara))*as.numeric(Y)})
  
  dVA.denom <- (pC*(-pF + 1)*pnorm(WBinner))
  dVA.denom[dVA.denom<=.Machine$double.eps] <- .Machine$double.eps
  dVA <- rbind( ((PRhat - 1)*P1*D2/PRhat + (PRhat - 1)*P2*D1/PRhat)/(P1*P2),
                (-(PRhat - 1)*P1*D2/PRhat - (PRhat - 1)*P2*D1/PRhat)/pC,
                - (PRhat - 1)*Del_v1d1/(PRhat*PBdv),
                ((pC*(-pF*((PRhat - 1)*P1*D2/PRhat + (PRhat - 1)*P2*D1/PRhat)/pC + (PRhat - 1)*Del_v1d1/(PRhat*pC))*pnorm(WBinner) + (-pF + 1)*(-(PRhat - 1)*P1*D2/PRhat - (PRhat - 1)*P2*D1/PRhat)*pnorm(WBinner))/(dVA.denom))
                
  )
  dVA  <- apply(regr$VA, 2, function(x){t(as.numeric(x)*t(dVA))*as.numeric(Y)})
  
  dVB <- rbind(0,
               -(-PFhat + 1)*dnorm(WBinner)/(PFhat*(-pnorm(WBinner) + 1)),
               (-PFhat + 1)*dnorm(WBinner)/(PFhat*pnorm(WBinner)),
               (-PFhat + 1)*dnorm(WBinner)/(PFhat*pnorm(WBinner))
  )  
  dVB  <- apply(regr$VB, 2, function(x){t(as.numeric(x)*t(dVB))*as.numeric(Y)})
  
  dSA.denom <- (pC*(1 - PBdv/pC)*pnorm(WBinner))
  dSA.denom[dSA.denom<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  dSA <- rbind( ((P1*D2 + P2*D1)/PRhat)/(P1P2),
               ((-P1*D2- P2*D1)/PRhat)/pC ,
               -Del_v1d1 /(PRhat*PBdv),
               (P1*D2 + P2*D1 - Del_v1d1)/(PRhat*(-pC+PBdv))
  )
  dSA  <- apply(regr$SA, 2, function(x){t(as.numeric(x)*t(dSA))*as.numeric(Y)})
  
  dCB.denom <- (PFhat*(pnorm(WBinner, lower=F))) 
  dCB.denom[dCB.denom<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  dCB.denom2 <- (PFhat*pnorm(WBinner))
  dCB.denom2[dCB.denom2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  dCB <- rbind(0,
               dnorm(WBinner)/dCB.denom,
               -dnorm(WBinner)/dCB.denom2,
               - dnorm(WBinner)/dCB.denom2
  )
  dCB  <- apply(regr$CB, 2, function(x){t(as.numeric(x)*t(dCB))*as.numeric(Y)})
  
  
  out <-  -cbind(dSA, dVA, dCB, dWA, dWB, dbara, dVB)
  return(colSums(out)) 
}




#### Jacobian of the PL method wrt to theta####
eval_gr_qll.i <- function(x,PRhat,PFhat,Y,regr){
  # here x = (theta,p)
  ##EQ[,1]= pR
  ##EQ[,2]= pC
  ##EQ[,3]= pF
  M <- dim(Y)[2]
  
  
  param <-vec2U.regr(x,regr)
  
  
  c <- cStar.jo(PRhat,param)
  pR <- f.jo(PFhat, param)
  pR[pR<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps) 
  pC <- g.jo(c, param)
  pC[pC<=(.Machine$double.eps)^(.25)] <- (.Machine$double.eps)^(.25) #edited 7 March
  pC[pC>=1-(.Machine$double.eps)^(.25)] <- 1-(.Machine$double.eps)^(.25) #edited 7 March
  pF <- h.jo(c, param)/pC
  v1 <- (c-param$barWA)/param$sig
  v2 <- (c-param$bara)/param$sig
  d1 <- (param$barWA - param$bara)/(param$sig*sqrt(2))
  d2 <- (param$barWA - c)/(param$sig)
  
  r <- 1/sqrt(2)
  VApr <- (param$VA/pR - param$VA * (pR-1)/(pR^2))
  P1 <- pnorm(v1)
  P1[P1<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  P1[P1>=1-sqrt(.Machine$double.eps)] <- 1-sqrt(.Machine$double.eps)
  P2 <- pnorm(v2)
  P2[P2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  P2[P2>=1-sqrt(.Machine$double.eps)] <- 1-sqrt(.Machine$double.eps)
  D1 <- dnorm(v1)
  D2 <- dnorm(v2)
  P1P2 <-  P1*P2 #more stable than pbiv w/rho=0
  P1P2[P1P2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  
  Del_d1v1 <- delPbivnorm(d1, -v1, r)
  Del_v1d1 <- delPbivnorm(-v1, d1, r)
  PBdv <- pbivnorm(d1, -v1, rho=r)
  PBdv[PBdv<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps) 
  PBdv[PBdv>= (1-sqrt(.Machine$double.eps))] <- 1-sqrt(.Machine$double.eps)
  WBinner <- as.numeric((-param$CB + param$VB*(-PFhat + 1) + param$barWB*PFhat)/PFhat)
  
  
  dWA4.denom <- (pC*(1 - PBdv/pC)*pnorm(WBinner))  
  dWA4.denom[dWA4.denom<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps) 
  dWA <- rbind( -D1/P1,
                P2*D1/pC,
                (sqrt(2)*Del_d1v1/2 +Del_v1d1)/PBdv,
                (pC*((-sqrt(2)*Del_d1v1/2 - Del_v1d1)/pC + P2*PBdv*D1/pC**2)*pnorm(WBinner) + (1 - PBdv/pC)*P2*pnorm(WBinner)*D1)/dWA4.denom  
  )
  
  dWA  <- apply(regr$barWA, 2, function(x){colSums(t(as.numeric(x)*t(dWA))*as.numeric(Y))})
  
  dWB.denom <- (-pnorm(WBinner) + 1)
  dWB.denom[dWB.denom<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  dWB.denom2 <- pnorm(WBinner)
  dWB.denom2[dWB.denom2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  dWB <- rbind(0,
               (-dnorm(WBinner)/dWB.denom),
               ( dnorm(WBinner)/dWB.denom2),
               ( dnorm(WBinner)/dWB.denom2)
  ) 
  
  dWB <- apply(regr$barWB, 2, function(x){colSums(t(as.numeric(x)*t(dWB))*as.numeric(Y))})
  
  
  dbara.denom <- (pC*(-pF + 1)*pnorm(WBinner))
  dbara.denom[dbara.denom<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  dbara.stuff <- (-2*P1*P2 + 2)
  dbara.stuff[dbara.stuff<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  dbara <- rbind( -D2/P2,
                  P1*D2/pC,
                  - sqrt(2)*Del_d1v1/(2*PBdv),
                  ((pC*(sqrt(2)*Del_d1v1/dbara.stuff + pF*P1*D2/pC)*pnorm(WBinner) + (-pF + 1)*P1*pnorm(WBinner)*D2)/(dbara.denom))
                  
  )
  dbara <- apply(regr$bara, 2, function(x){colSums(t(as.numeric(x)*t(dbara))*as.numeric(Y))})
  
  dVA.denom <- (pC*(-pF + 1)*pnorm(WBinner))
  dVA.denom[dVA.denom<=.Machine$double.eps] <- .Machine$double.eps
  dVA <- rbind( ((PRhat - 1)*P1*D2/PRhat + (PRhat - 1)*P2*D1/PRhat)/(P1*P2),
                (-(PRhat - 1)*P1*D2/PRhat - (PRhat - 1)*P2*D1/PRhat)/pC,
                - (PRhat - 1)*Del_v1d1/(PRhat*PBdv),
                ((pC*(-pF*((PRhat - 1)*P1*D2/PRhat + (PRhat - 1)*P2*D1/PRhat)/pC + (PRhat - 1)*Del_v1d1/(PRhat*pC))*pnorm(WBinner) + (-pF + 1)*(-(PRhat - 1)*P1*D2/PRhat - (PRhat - 1)*P2*D1/PRhat)*pnorm(WBinner))/(dVA.denom))
                
  )
  dVA  <- apply(regr$VA, 2, function(x){colSums(t(as.numeric(x)*t(dVA))*as.numeric(Y))})
  
  dVB <- rbind(0,
               -(-PFhat + 1)*dnorm(WBinner)/(PFhat*(-pnorm(WBinner) + 1)),
               (-PFhat + 1)*dnorm(WBinner)/(PFhat*pnorm(WBinner)),
               (-PFhat + 1)*dnorm(WBinner)/(PFhat*pnorm(WBinner))
  )  
  dVB  <- apply(regr$VB, 2, function(x){colSums(t(as.numeric(x)*t(dVB))*as.numeric(Y))})
  
  dSA.denom <- (pC*(1 - PBdv/pC)*pnorm(WBinner))
  dSA.denom[dSA.denom<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  dSA <- rbind( ((P1*D2 + P2*D1)/PRhat)/(P1P2),
                ((-P1*D2- P2*D1)/PRhat)/pC ,
                -Del_v1d1 /(PRhat*PBdv),
                (P1*D2 + P2*D1 - Del_v1d1)/(PRhat*(-pC+PBdv))
  )
  dSA  <- apply(regr$SA, 2, function(x){colSums(t(as.numeric(x)*t(dSA))*as.numeric(Y))})
  
  dCB.denom <- (PFhat*(pnorm(WBinner, lower=F))) 
  dCB.denom[dCB.denom<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  dCB.denom2 <- (PFhat*pnorm(WBinner))
  dCB.denom2[dCB.denom2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  dCB <- rbind(0,
               dnorm(WBinner)/dCB.denom,
               -dnorm(WBinner)/dCB.denom2,
               - dnorm(WBinner)/dCB.denom2
  )
  dCB  <- apply(regr$CB, 2, function(x){colSums(t(as.numeric(x)*t(dCB))*as.numeric(Y))})
  
  
  out <-  -cbind(dSA, dVA, dCB, dWA, dWB, dbara, dVB)
  return(out)
}






##### Jacobian of the PL wrt to Phat ####
eval_gr_qll.ip <- function(PRhat, PFhat, x, Y,regr){
  # here x = (theta,p)
  ##EQ[,1]= pR
  ##EQ[,2]= pC
  ##EQ[,3]= pF
  M <- dim(Y)[2]

  param <-vec2U.regr(x,regr)
  
  
  c <- cStar.jo(PRhat,param)
  pR <- f.jo(PFhat, param)
  pR[pR<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps) 
  pC <- g.jo(c, param)
  pC[pC<=(.Machine$double.eps)^(.25)] <- (.Machine$double.eps)^(.25) #edited 7 March
  pC[pC>=1-(.Machine$double.eps)^(.25)] <- 1-(.Machine$double.eps)^(.25) #edited 7 March
  pF <- h.jo(c, param)/pC
  v1 <- (c-param$barWA)/param$sig
  v2 <- (c-param$bara)/param$sig
  d1 <- (param$barWA - param$bara)/(param$sig*sqrt(2))
  d2 <- (param$barWA - c)/(param$sig)
  
  r <- 1/sqrt(2)
  VApr <- (param$VA/pR - param$VA * (pR-1)/(pR^2))
  P1 <- pnorm(v1)
  P1[P1<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  P1[P1>=1-sqrt(.Machine$double.eps)] <- 1-sqrt(.Machine$double.eps)
  P2 <- pnorm(v2)
  P2[P2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  P2[P2>=1-sqrt(.Machine$double.eps)] <- 1-sqrt(.Machine$double.eps)
  D1 <- dnorm(v1)
  D2 <- dnorm(v2)
  P1P2 <-  P1*P2 #more stable than pbiv w/rho=0
  P1P2[P1P2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  
  Del_d1v1 <- delPbivnorm(d1, -v1, r)
  Del_v1d1 <- delPbivnorm(-v1, d1, r)
  PBdv <- pbivnorm(d1, -v1, rho=r)
  PBdv[PBdv<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  PBdv[PBdv>= (1-sqrt(.Machine$double.eps))] <- 1-sqrt(.Machine$double.eps) 
  WBinner <- as.numeric((-param$CB + param$VB*(-PFhat + 1) + param$barWB*PFhat)/PFhat)
  pWBinner <- pnorm(WBinner)
  invpWBinner <- pnorm(WBinner, lower.tail = FALSE)
  pWBinner[pWBinner<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  pWBinner[pWBinner>=1-sqrt(.Machine$double.eps)] <- 1-sqrt(.Machine$double.eps)
  invpWBinner[invpWBinner<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  invpWBinner[invpWBinner>=1-sqrt(.Machine$double.eps)] <- 1-sqrt(.Machine$double.eps)
  
  dPRhat  <- rbind(
    (P1*(param$VA/PRhat + (-param$SA + param$VA*(-PRhat + 1))/PRhat**2)*D2 + P2*(param$VA/PRhat + (-param$SA + param$VA*(-PRhat + 1))/PRhat**2)*D1)/(P1P2) ,
    (-P1*(param$VA/PRhat + (-param$SA + param$VA*(-PRhat + 1))/PRhat**2)*D2 - P2*(param$VA/PRhat + (-param$SA + param$VA*(-PRhat + 1))/PRhat**2)*D1)/pC,
    (-param$VA/PRhat + (param$SA - param$VA*(-PRhat + 1))/PRhat**2)*Del_v1d1/PBdv,
    (pC*(-(-param$VA/PRhat + (param$SA - param$VA*(-PRhat + 1))/PRhat**2)*Del_v1d1/pC - (P1*(param$VA/PRhat + (-param$SA + param$VA*(-PRhat + 1))/PRhat**2)*D2+ P2*(param$VA/PRhat + (-param$SA + param$VA*(-PRhat + 1))/PRhat**2)*D1)*PBdv/pC**2)*pnorm(WBinner) + (1 - PBdv/pC)*(-P1*(param$VA/PRhat + (-param$SA + param$VA*(-PRhat + 1))/PRhat**2)*D2- P2*(param$VA/PRhat + (-param$SA + param$VA*(-PRhat + 1))/PRhat**2)*D1)*pnorm(WBinner))/(pC*(1 - PBdv/pC)*pnorm(WBinner)) 
  )
  dPRhat <- diag(colSums(dPRhat*Y)) 
  
  dPFhat  <- rbind(
    0,
    -((-param$VB +param$barWB)/PFhat + (param$CB - param$VB*(-PFhat + 1) -param$barWB*PFhat)/PFhat**2)*dnorm(WBinner)/(invpWBinner),
    ((-param$VB + param$barWB)/PFhat + (param$CB - param$VB*(-PFhat + 1) - param$barWB*PFhat)/PFhat**2)*dnorm(WBinner)/pWBinner,
    ((-param$VB + param$barWB)/PFhat + (param$CB - param$VB*(-PFhat + 1) - param$barWB*PFhat)/PFhat**2)*dnorm(WBinner)/pWBinner
  )
  dPFhat <- diag(colSums(dPFhat*Y))
  
  
  return(-cbind(dPRhat, dPFhat))
}

#### Jacobian of the Best Response Function wrt to Phat ####
dPsiDp <- function(PRhat, PFhat, x, Y,regr){
  # here x = (theta,p)
  ##EQ[,1]= pR
  ##EQ[,2]= pC
  ##EQ[,3]= pF
  M <- dim(Y)[2]
  
  param <-vec2U.regr(x,regr)
  
  
  c <- cStar.jo(PRhat,param)
  pR <- f.jo(PFhat, param)
  pR[pR<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps) 
  pC <- g.jo(c, param)
  pC[pC<=(.Machine$double.eps)^(.25)] <- (.Machine$double.eps)^(.25) #edited 7 March
  pC[pC>=1-(.Machine$double.eps)^(.25)] <- 1-(.Machine$double.eps)^(.25) #edited 7 March
  pF <- h.jo(c, param)/pC
  v1 <- (c-param$barWA)/param$sig
  v2 <- (c-param$bara)/param$sig
  d1 <- (param$barWA - param$bara)/(param$sig*sqrt(2))
  d2 <- (param$barWA - c)/(param$sig)
  
  r <- 1/sqrt(2)
  VApr <- (param$VA/pR - param$VA * (pR-1)/(pR^2))
  P1 <- pnorm(v1)
  P1[P1<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  P1[P1>=1-sqrt(.Machine$double.eps)] <- 1-sqrt(.Machine$double.eps)
  P2 <- pnorm(v2)
  P2[P2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  P2[P2>=1-sqrt(.Machine$double.eps)] <- 1-sqrt(.Machine$double.eps)
  D1 <- dnorm(v1)
  D2 <- dnorm(v2)
  P1P2 <-  P1*P2 #more stable than pbiv w/rho=0
  P1P2[P1P2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  
  Del_d1v1 <- delPbivnorm(d1, -v1, r)
  Del_v1d1 <- delPbivnorm(-v1, d1, r)
  PBdv <- pbivnorm(d1, -v1, rho=r)
  PBdv[PBdv<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps) 
  PBdv[PBdv>= (1-sqrt(.Machine$double.eps))] <- 1-sqrt(.Machine$double.eps) 
  WBinner <- as.numeric((-param$CB + param$VB*(-PFhat + 1) + param$barWB*PFhat)/PFhat)
  
  
  dPRhat <- ((-param$VB + param$barWB)/PFhat + (param$CB - param$VB*(-PFhat + 1) - param$barWB*PFhat)/PFhat**2)*dnorm(WBinner)
  dPRhat <- diag(dPRhat) 
  
  dPFhat  <- (-param$VA/PRhat + (param$SA - param$VA*(-PRhat + 1))/PRhat**2)*Del_v1d1/pC + (P1*(param$VA/PRhat + (-param$SA + param$VA*(-PRhat + 1))/PRhat**2)*D2 + P2*(param$VA/PRhat + (-param$SA + param$VA*(-PRhat + 1))/PRhat**2)*D1)*PBdv/pC**2
  dPFhat <- diag(dPFhat) 
  ZERO <- matrix(0, nrow=M, ncol=M)
  
  return( rbind(cbind(ZERO, dPRhat), cbind(dPFhat, ZERO)) ) 
}



#### Jacobian of the Best Response Function wrt to theta ####
dPsi.dTheta <- function(x, PRhat, PFhat, Y, regr){
  # here x = (theta,p)
  ##EQ[,1]= pR
  ##EQ[,2]= pC
  ##EQ[,3]= pF
  M <- dim(Y)[2]
  
  param <-vec2U.regr(x,regr)
  
  
  c <- cStar.jo(PRhat,param)
  pR <- f.jo(PFhat, param)
  pR[pR<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps) 
  pC <- g.jo(c, param)
  pC[pC<=(.Machine$double.eps)^(.25)] <- (.Machine$double.eps)^(.25) #edited 7 March
  pC[pC>=1-(.Machine$double.eps)^(.25)] <- 1-(.Machine$double.eps)^(.25) #edited 7 March
  pF <- h.jo(c, param)/pC
  v1 <- (c-param$barWA)/param$sig
  v2 <- (c-param$bara)/param$sig
  d1 <- (param$barWA - param$bara)/(param$sig*sqrt(2))
  d2 <- (param$barWA - c)/(param$sig)
  
  r <- 1/sqrt(2)
  VApr <- (param$VA/pR - param$VA * (pR-1)/(pR^2))
  P1 <- pnorm(v1)
  P1[P1<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  P1[P1>=1-sqrt(.Machine$double.eps)] <- 1-sqrt(.Machine$double.eps)
  P2 <- pnorm(v2)
  P2[P2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  P2[P2>=1-sqrt(.Machine$double.eps)] <- 1-sqrt(.Machine$double.eps)
  D1 <- dnorm(v1)
  D2 <- dnorm(v2)
  P1P2 <-  P1*P2 #more stable than pbiv w/rho=0
  P1P2[P1P2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  
  Del_d1v1 <- delPbivnorm(d1, -v1, r)
  Del_v1d1 <- delPbivnorm(-v1, d1, r)
  PBdv <- pbivnorm(d1, -v1, rho=r)
  PBdv[PBdv<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps) 
  PBdv[PBdv>= (1-sqrt(.Machine$double.eps))] <- 1-sqrt(.Machine$double.eps) 
  WBinner <- as.numeric((-param$CB + param$VB*(-PFhat + 1) + param$barWB*PFhat)/PFhat)
  
  
  dPRdSA <- matrix(0, nrow=M, ncol=ncol(regr$SA))
  dPRdVA <- matrix(0, nrow=M, ncol=ncol(regr$VA))
  dPRda <- matrix(0, nrow=M, ncol=ncol(regr$bara))
  dPRdWA <- matrix(0, nrow=M, ncol=ncol(regr$barWA))
  
  dPRdCB <- -dnorm(WBinner)/PFhat
  dPRdCB  <- apply(regr$CB, 2, function(x){t(as.numeric(x)*t(dPRdCB))})
  dPRdWB <- dnorm(WBinner)
  dPRdWB  <- apply(regr$barWB, 2, function(x){t(as.numeric(x)*t(dPRdWB))})
  dPRdVB <- (-PFhat + 1)*dnorm(WBinner)/PFhat
  dPRdVB  <- apply(regr$VB, 2, function(x){t(as.numeric(x)*t(dPRdVB))})
  
  dPFdCB <- matrix(0, nrow=M, ncol=ncol(regr$CB))
  dPFdWB <- matrix(0, nrow=M, ncol=ncol(regr$barWB))
  dPFdVB <- matrix(0, nrow=M, ncol=ncol(regr$VB))
  
  dPFdSA <- (P1*D2/PRhat + P2*D1/PRhat)*PBdv/pC**2 - Del_v1d1/(PRhat*pC)
  dPFdSA  <- apply(regr$SA, 2, function(x){t(as.numeric(x)*t(dPFdSA))})
  dPFdVA <- (P1*(PRhat - 1)*D2/PRhat + P2*(PRhat - 1)*D1/PRhat)*PBdv/pC**2 - (PRhat - 1)*Del_v1d1/(PRhat*pC)
  dPFdVA  <- apply(regr$VA, 2, function(x){t(as.numeric(x)*t(dPFdVA))})
  dPFdWA <- -P2*PBdv*D1/pC**2 + (sqrt(2)*Del_d1v1/2 + Del_v1d1)/pC
  dPFdWA  <- apply(regr$barWA, 2, function(x){t(as.numeric(x)*t(dPFdWA))})
  dPFda <- -P1*PBdv*D2/pC**2 - sqrt(2)*Del_d1v1/(-2*P1P2 + 2)
  dPFda  <- apply(regr$bara, 2, function(x){t(as.numeric(x)*t(dPFda))})
  
  return( rbind(cbind(dPRdSA, dPRdVA, dPRdCB, dPRdWA, dPRdWB, dPRda, dPRdVB),
                cbind(dPFdSA, dPFdVA, dPFdCB, dPFdWA, dPFdWB, dPFda, dPFdVB))
          )
}


#### Jacobian of the equilibrium constraint wrt to p_R for tML implementation #### 
eval_gr_fh <- function(pr, U){
  
  r = 1/sqrt(2)
  cstar <- (U$SA - (1-pr)*U$VA)/pr
  frac1 <- cstar - U$barWA
  frac2 <- cstar-U$bara
  frac3 <-  (U$barWA - U$bara)/sqrt(2)
  
  MVT <- pbivnorm(frac3, -frac1, rho=r)
  dMVT <- delPbivnorm(-frac1, frac3, rho=r)
  dN1 <- dnorm(frac1)
  dN2 <- dnorm(frac2)
  N1 <- pnorm(frac1)
  N2 <- pnorm(frac2)
  
  denom <- (pr^2)* (MVT^2)
  inner <- MVT*(N1*dN2 + N2*dN1) - (1 - N1*N2)*dMVT
  outer <- (U$SA - U$VA) * (U$CB - U$VB) * dnorm(U$barWB-U$VB + (U$CB - U$VB)*(-1 + N1*N2)/MVT)
  
  return(-outer*inner/denom)
} 




##### CMLE functions for standard errors ##### 

eval_gr_obj <- function(x, Y, regr){
  # here x = (theta,p)
  ##EQ[,1]= pR
  ##EQ[,2]= pC
  ##EQ[,3]= pF
  M <- dim(Y)[2]
  prStar <- x[(length(x)-M+1):length(x)]
  exPR = exp(-prStar)
  exPR[exPR >= (.Machine$double.xmax)^(1/4)] <- (.Machine$double.xmax)^(1/4)
  
  xP <- plogis(prStar)
  xP[xP <= sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  xP[xP>=1]<-1-sqrt(.Machine$double.eps) 
  xT <- x[1:(length(x)-M)]
  
  param <- vec2U.regr(xT,regr)
  c <- cStar.jo(xP,param)
  pR <- xP
  pC <- g.jo(c, param)
  pC[pC <= sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  pC[pC>=1]<-1-sqrt(.Machine$double.eps)  
  pF <- h.jo(c, param)/pC
  pF[pF <= sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  pF[pF>=1]<-1-sqrt(.Machine$double.eps)  
  v1 <- (c-param$barWA)/param$sig
  v2 <- (c-param$bara)/param$sig
  d1 <- (param$barWA - param$bara)/(param$sig*sqrt(2))
  d2 <- (param$barWA - c)/(param$sig)
  
  r <- 1/sqrt(2)
  P1 <- pnorm(v1)
  P1[P1<=sqrt(.Machine$double.eps)]<-sqrt(.Machine$double.eps)
  P1[P1>=1]<-1-sqrt(.Machine$double.eps)
  P2 <- pnorm(v2)
  P2[P2<=sqrt(.Machine$double.eps)]<-sqrt(.Machine$double.eps)
  P2[P2>=1]<-1-sqrt(.Machine$double.eps)
  D1 <- dnorm(v1)
  D2 <- dnorm(v2)
  P1P2 <- P1*P2
  P1P2[P1P2<=sqrt(.Machine$double.eps)]<-sqrt(.Machine$double.eps)
  P1P2[P1P2>=1]<-1-sqrt(.Machine$double.eps)
  pRratio <- (pR-1)/pR
  Del_d1v1 <- delPbivnorm(d1, -v1, r)
  Del_v1d1 <- delPbivnorm(-v1, d1, r)
  PBdv <- pbivnorm(d1, -v1, r)
  PBdv[PBdv<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  
  
  dSA = rbind( (((1 + exPR)*P1*D2 + (1 + exPR)*P2*D1)/(P1P2)),
               (((1 + exPR)*P1*D2 - (1 + exPR)*P2*D1)/pC),
               ((-1 - exPR)*Del_v1d1/PBdv),
               ((1 + exPR)*(pR*pC*(-pF*((1 + exPR)*P1*D2 + (1 + exPR)*P2*D1)/pC - (-1 - exPR)*Del_v1d1/pC) + pR*(-pF + 1)*(-(1 + exPR)*P1*D2 - (1 + exPR)*P2*D1))/(pC*(-pF + 1)))
  )
  
  dSA  <- apply(regr$SA, 2, function(x){t(as.numeric(x)*t(dSA))*as.numeric(Y)})
  
  
  dVA = rbind(  (((1 + exPR)*(pR - 1)*P1*D2 + (1 + exPR)*(pR - 1)*P2*D1)/(P1P2)),
                ((-(1 +exPR)*(pR - 1)*P1*D2 - (1 +exPR)*(pR - 1)*P2*D1)/pC),
                -((1 + exPR)*(pR - 1)*Del_v1d1/PBdv),
                ((1 + exPR)*(pR*pC*(-pF*((1 + exPR)*(pR - 1)*P1*D2 + (1 + exPR)*(pR - 1)*P2*D1)/pC + (1 + exPR)*(pR - 1)*Del_v1d1/pC) + pR*(-pF + 1)*(-(1 + exPR)*(pR - 1)*P1*D2 - (1 + exPR)*(pR - 1)*P2*D1))/(pC*(-pF + 1)))
  )
  dVA  <- apply(regr$VA, 2, function(x){t(as.numeric(x)*t(dVA))*as.numeric(Y)})
  
  dCB <- matrix(0, M*4, ncol(regr$CB))
  
  
  dWA = rbind( (-D1/P1), #y=0
               (P2*D1/pC), #y=1,
               ((sqrt(2)/2.0 * Del_d1v1 + Del_v1d1)/PBdv), #y=2
               ((1 + exPR)*(pR*pC*(pF*P2*D1/pC + (-sqrt(2)*Del_d1v1/2 - Del_v1d1)/pC) + pR*(-pF + 1)*P2*D1)/(pC*(-pF + 1)))
               
  )
  dWA  <- apply(regr$barWA, 2, function(x){t(as.numeric(x)*t(dWA))*as.numeric(Y)})
  
  
  dWB <- matrix(0, M*4, ncol(regr$barWB))
  
  
  dbara = rbind( (-D2/P2),
                 (P1*D2/pC),
                 -(sqrt(2)*Del_d1v1/(2*PBdv)) ,
                 ((1 + exPR)*(pR*pC*(sqrt(2)*Del_d1v1/(-2*P1P2 + 2) + pF*P1*D2/pC) + pR*(-pF + 1)*P1*D2)/(pC*(-pF + 1)))
  )
  dbara <- apply(regr$bara, 2, function(x){t(as.numeric(x)*t(dbara))*as.numeric(Y)})
  
  dVB <-  matrix(0, M*4, ncol(regr$VB))
  
  
  dpR <- rbind((((param$VA*pR*exPR + (-param$SA + param$VA*(-pR + 1))*exPR)*P1*D2 + (param$VA*pR*exPR + (-param$SA + param$VA*(-pR + 1))*exPR)*P2*D1)/(P1P2)),
               exPR * (1 + ((D2*P1 + D1*P2)*(param$SA - param$VA))/pC + 1/(-1 + pR) + pR),
               (exPR *(1 + exPR)* pR* (PBdv* pR + (param$SA - param$VA) *Del_v1d1))/PBdv,
               ((1 +exPR)*(pR**2*pC*(-pF + 1)*exPR + pR*pC*(-pF*((param$VA*pR*exPR + (-param$SA + param$VA*(-pR + 1))*exPR)*P1*D2 + (param$VA*pR*exPR + (-param$SA + param$VA*(-pR + 1))*exPR)*P2*D1)/pC - (-param$VA*pR*exPR + (param$SA - param$VA*(-pR + 1))*exPR)*Del_v1d1/pC) + pR*(-pF + 1)*(-(param$VA*pR*exPR + (-param$SA + param$VA*(-pR + 1))*exPR)*P1*D2 - (param$VA*pR*exPR + (-param$SA + param$VA*(-pR + 1))*exPR)*P2*D1))/(pC*(-pF + 1)))
  )
  dpR <-colSums(dpR*Y)
  out <-  - c(colSums(cbind(dSA, dVA, dCB, dWA, dWB, dbara, dVB)), dpR)
  if(anyNA(out) | any(is.infinite(out))){
    grNA <<-x
    stop("NA in gradient")
  }
  return(out) 
}


Jconst <- function(x, Y, regr){
  M <- dim(Y)[2]
  prStar <- x[(length(x)-M+1):length(x)]
  xT <- x[1:(length(x)-M)]
  xP <- plogis(prStar)
  exPR <- exp(-prStar)
  exPR[is.infinite(exPR)] <- .Machine$double.xmax #edit added March 9
  param <- vec2U.regr(xT,regr)
  c <- cStar.jo(xP,param)
  pR <- xP
  pC <- g.jo(c, param)
  pC[pC <= sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  pF <- h.jo(c, param)/pC
  v1 <- (c-param$barWA)/param$sig
  v2 <- (c-param$bara)/param$sig
  d1 <- (param$barWA - param$bara)/(param$sig*sqrt(2))
  d2 <- (param$barWA - c)/(param$sig)
  
  r <- 1/sqrt(2)
  
  
  
  P1 <- pnorm(v1)
  P1[P1<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  P2 <- pnorm(v2)
  P2[P2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  D1 <- dnorm(v1)
  D2 <- dnorm(v2)
  P1P2 <- P1*P2
  P1P2[P1P2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps) #added march 8
  P1P2[P1P2>=1-sqrt(.Machine$double.eps)] <- 1-sqrt(.Machine$double.eps) #added march 8
  Del_d1v1 <- delPbivnorm(d1, -v1, r)
  Del_v1d1 <- delPbivnorm(-v1, d1, r)
  PBdv <- pbivnorm(d1, -v1, r)
  PBdv[PBdv<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  dBig <- dnorm((pC*(-param$CB+param$VB*(-pF + 1) + param$barWB*pF)/PBdv))
  WBinner <- as.numeric((-param$CB + param$VB*(-pF+ 1) + param$barWB*pF)/pF)
  
  
  dSA = -(-pC*(-1 - exPR)*(-param$CB + param$VB*(-pF + 1) + param$barWB*pF)*Del_v1d1/PBdv**2 + pC*(param$VB*(-pF*((1 + exPR)*P1*D2 + (1 + exPR)*P2*D1)/pC - (-1 - exPR)*Del_v1d1/pC) + param$barWB*pF*((1 + exPR)*P1*D2 + (1 + exPR)*P2*D1)/pC + param$barWB*(-1 - exPR)*Del_v1d1/pC)/PBdv + (-(1 + exPR)*P1*D2 - (1 + exPR)*P2*D1)*(-param$CB + param$VB*(-pF + 1) + param$barWB*pF)/PBdv)*dBig
  dSA <- dSA * regr$SA
  
  dVA = -(pC*(1 +exPR)*(pR - 1)*(-param$CB + param$VB*(-pF + 1) + param$barWB*pF)*Del_v1d1/PBdv**2 + pC*(param$VB*(-pF*((1 +exPR)*(pR - 1)*P1*D2 + (1 +exPR)*(pR - 1)*P2*D1)/pC + (1 +exPR)*(pR - 1)*Del_v1d1/pC) + param$barWB*pF*((1 +exPR)*(pR - 1)*P1*D2 + (1 +exPR)*(pR - 1)*P2*D1)/pC - param$barWB*(1 +exPR)*(pR - 1)*Del_v1d1/pC)/PBdv + (-(1 +exPR)*(pR - 1)*P1*D2 - (1 +exPR)*(pR - 1)*P2*D1)*(-param$CB + param$VB*(-pF + 1) + param$barWB*pF)/PBdv)*dBig
  dVA <- dVA * regr$VA
  
  dCB = pC*dBig/PBdv
  dCB <- dCB * regr$CB
  
  dWA = -(pC*(-sqrt(2)*Del_d1v1/2 - Del_v1d1)*(-param$CB + param$VB*(-pF + 1) + param$barWB*pF)/PBdv**2 + pC*(param$VB*(pF*P2*D1/pC + (-sqrt(2)*Del_d1v1/2 - Del_v1d1)/pC) - param$barWB*pF*P2*D1/pC + param$barWB*(sqrt(2)*Del_d1v1/2 + Del_v1d1)/pC)/PBdv + (-param$CB + param$VB*(-pF + 1) + param$barWB*pF)*P2*D1/PBdv)*dBig
  dWA <- dWA * regr$barWA
  
  dWB <- -dBig
  dWB <- dWB * regr$barWB
  
  dbara = -(sqrt(2)*pC*(-param$CB + param$VB*(-pF + 1) + param$barWB*pF)*Del_d1v1/(2*PBdv**2) +
              pC*(param$VB*(sqrt(2)*Del_d1v1/(-2*P1P2 + 2) + pF*P1*D2/pC) - sqrt(2)*param$barWB*Del_d1v1/(-2*P1P2 + 2) - param$barWB*pF*P1*D2/pC)/PBdv +
              (-param$CB + param$VB*(-pF + 1) + param$barWB*pF)*P1*D2/PBdv)*dBig
  dbara <- dbara * regr$bara
  
  dVB = -pC*(-pF + 1)*dBig/PBdv
  dVB <- dVB * regr$VB
  
  dpR = pR**2*exPR - (-pC*(-param$VA*pR*exPR + (param$SA - param$VA*(-pR + 1))*exPR)*(-param$CB + param$VB*(-pF + 1) + param$barWB*pF)*Del_v1d1/PBdv**2 + pC*(param$VB*(-pF*((param$VA*pR*exPR + (-param$SA + param$VA*(-pR + 1))*exPR)*P1*D2 + (param$VA*pR*exPR + (-param$SA + param$VA*(-pR + 1))*exPR)*P2*D1)/pC - (-param$VA*pR*exPR + (param$SA - param$VA*(-pR + 1))*exPR)*Del_v1d1/pC) + param$barWB*pF*((param$VA*pR*exPR + (-param$SA + param$VA*(-pR + 1))*exPR)*P1*D2 + (param$VA*pR*exPR + (-param$SA + param$VA*(-pR + 1))*exPR)*P2*D1)/pC + param$barWB*(-param$VA*pR*exPR + (param$SA - param$VA*(-pR + 1))*exPR)*Del_v1d1/pC)/PBdv + (-(param$VA*pR*exPR + (-param$SA + param$VA*(-pR + 1))*exPR)*P1*D2 - (param$VA*pR*exPR + (-param$SA + param$VA*(-pR + 1))*exPR)*P2*D1)*(-param$CB + param$VB*(-pF + 1) + param$barWB*pF)/PBdv)*dBig
  
  
  Mat <- cbind(dSA, dVA, dCB, dWA, dWB, dbara, dVB, diag(as.numeric(dpR))) 
  return(Mat)
  
}
