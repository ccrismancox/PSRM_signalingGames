
LL.jo <- function(x, Y,regr){
  # here x = (theta,p)
  M <- dim(Y)[2]
  xP <- plogis(x[(length(x)-M+1):length(x)])
  xT <- x[1:(length(x)-M)]
  
  U <- vec2U.regr(xT,regr)
  EQ <- eqProbs(xP,U,T)
  OUT <- cbind(1-EQ[,2], EQ[,2]*(1-EQ[,1]), 	
               EQ[,2]*EQ[,1]*EQ[,3], EQ[,2]*EQ[,1]*(1-EQ[,3]))
  OUT[OUT <= sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  LL <- sum(log(t(OUT))*Y)
  return(-LL)
}


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
  
  
  Mat <- cbind(dSA, dVA, dCB, dWA, dWB, dbara, dVB, diag(as.numeric(dpR))) #for dense
  return(Mat)
  
}


