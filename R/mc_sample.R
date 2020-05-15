#' Sample
#' 
#' @param nrep Antal bootstraps
#' @param input Parameters
#' 
#' @importFrom deSolve ode
#' 
#' @export
mc_sample <- function(nrep, input, pop, 
                      reference_date,
                      first_date,
                      end_prediction_date,
                      contact_reduction_first_day,
                      rnd_func = crunif) {
  
  start_day <- as.numeric(first_date) - as.numeric(reference_date)
  end_day <- as.numeric(end_prediction_date) - as.numeric(reference_date)
  times <- seq(start_day, end_day, 1)
  
  time_contact_reduction_first_day <- as.numeric(contact_reduction_first_day) - 
    as.numeric(reference_date)
  
  Ntot <- sum(pop)
  
  ####################################################################

  
  SSEIRD <- function(t,x,p)
  {
    ## Unpack state by hand 
    n   <- length(x) / 19
    S   <- x[1:n]
    E1   <- x[n + (1:n)]
    I1R <- x[2*n + (1:n)] # infected at home - recovering
    I1M <- x[3*n + (1:n)] # infected at home - will go to hospital
    HR <- x[4*n + (1:n)] # at hospital - recovering
    HM <- x[5*n + (1:n)] # at hospital - will go to ICU
    CR <- x[6*n + (1:n)] # infected intensive care - recovering
    CD <- x[7*n + (1:n)] # infected intensive care - dying
    R   <- x[8*n + (1:n)]
    D   <- x[9*n + (1:n)]
    U   <- x[10*n + (1:n)]
    HC <- x[11*n + (1:n)] # cumulative hospital
    CC <- x[12*n + (1:n)] # cumulative hospital
    E2 <- x[13*n + (1:n)]
    E3 <- x[14*n + (1:n)]
    I2R <- x[15*n + (1:n)]
    I3R <- x[16*n + (1:n)]
    I2M <- x[17*n + (1:n)]
    I3M <- x[18*n + (1:n)]
    
    
    dS <- numeric(n)
    dE1 <- numeric(n)
    dI1R <- numeric(n)
    dI1M <- numeric(n)
    dHR <- numeric(n)
    dHM <- numeric(n)
    dCR <- numeric(n)
    dCD <- numeric(n)
    dR <- numeric(n)
    dD <- numeric(n)
    dU <- numeric(n)
    dHC <- numeric(n)
    dCC <- numeric(n)
    dE2 <- numeric(n)
    dE3 <- numeric(n)
    dI2R <-numeric(n)
    dI3R <-numeric(n)
    dI2M <-numeric(n)
    dI3M <-numeric(n)
    
    
    
    with(as.list(p),
         {
           # contact reduction
           #if (U[1]<Res2DK) {beta <- betaNow} else {beta <- betaNew}
           if (U[1]<time_contact_reduction_first_day) {beta <- betaNow} else {beta <- betaNew}
           
           # interactive adjust
           RR <- ER*0.01 
           
           for(i in 1:n)
           {
             dS[i]   <- - S[i] * beta[i,] %*% (I1R+I1M+I2R+I2M+I3R+I3M) * RR
             dE1[i]  <-   S[i] * beta[i,] %*% (I1R+I1M+I2R+I2M+I3R+I3M) * RR - 3*gammaEI1[i] * E1[i]
             dI1R[i] <- 3*gammaEI1[i] * E3[i]*pI1R[i]       - 3*recI1[i] * I1R[i]
             dI1M[i] <- 3*gammaEI1[i] * E3[i]*(1-pI1R[i])   - 3*gammaI12[i] * I1M[i]
             dHR[i]  <- 3*gammaI12[i] * I3M[i]*pI2R[i]      + recI3[i]*CR[i] - recI2[i]* HR[i]
             dHM[i]  <- 3*gammaI12[i] * I3M[i]*(1-pI2R[i])  - gammaI23[i] * HM[i]
             dCR[i]  <- gammaI23[i] * HM[i]*pI3R[i]         - recI3[i] * CR[i]
             dCD[i]  <- gammaI23[i] * HM[i]*(1-pI3R[i])     - muI3[i] * CD[i]
             dR[i]   <- 3*recI1[i]*I3R[i] + recI2[i]*HR[i] 
             dD[i]   <- muI3[i]*CD[i] 
             dU[i]   <- 1
             dHC[i]  <- 3*gammaI12[i] * I3M[i]
             dCC[i]  <- gammaI23[i] * HM[i]
             dE2[i]  <- 3*gammaEI1[i]*E1[i] - 3*gammaEI1[i]*E2[i]
             dE3[i]  <- 3*gammaEI1[i]*E2[i] - 3*gammaEI1[i]*E3[i]
             dI2R[i] <- 3*recI1[i]*I1R[i] - 3*recI1[i]*I2R[i]
             dI3R[i] <- 3*recI1[i]*I2R[i] - 3*recI1[i]*I3R[i]
             dI2M[i] <- 3*gammaI12[i] * I1M[i] - 3*gammaI12[i] * I2M[i]
             dI3M[i] <- 3*gammaI12[i] * I2M[i] - 3*gammaI12[i] * I3M[i]
           }
           dX <- c(dS,dE1,dI1R,dI1M,dHR,dHM,dCR,dCD,dR,dD,dU,dHC,dCC,
                   dE2,dE3,dI2R,dI3R,dI2M,dI3M)
           #if(abs(sum(dX)-2)>1e-8) stop("Non-conservation!")
           return(list(dX))
         }
    )
    
  }
  
  #withProgress(message = 'Simulerer -', value = 0, {
    
    #nrep <- as.numeric(input$nrep)
    
    for (i in 1:nrep)
    {
      ## Two daily interacting groups 
      #n <- 2
      betaNow <- array(c(rnd_func(input$r.betaWIlt),rnd_func(input$r.betaB),
                         rnd_func(input$r.betaB),rnd_func(input$r.betaWIge)),c(2,2))
      beta1 <- array(c(rnd_func(input$r.betaWIlt1),rnd_func(input$r.betaB1),
                       rnd_func(input$r.betaB1),rnd_func(input$r.betaWIge1)),c(2,2))
      # beta4 <- array(c(crunif(input$r.betaWIlt4),crunif(input$r.betaB4),
      #                  crunif(input$r.betaB4),crunif(input$r.betaWIge4)),c(2,2))
      # beta2 <- array(c(crunif(input$r.betaWIlt2),crunif(input$r.betaB2),
      #                  crunif(input$r.betaB2),crunif(input$r.betaWIge2)),c(2,2))
      # beta3 <- array(c(crunif(input$r.betaWIlt3),crunif(input$r.betaB3),
      #                  crunif(input$r.betaB3),crunif(input$r.betaWIge3)),c(2,2))
      #betaNew <- switch(as.numeric(input$scenarie),betaNow,beta1,beta2,beta3,beta4)
      betaNew <- beta1
      
      # efficiency of restrictions
      ER <- rnd_func(input$r.ER)
      
      # mortality rates for different stages
      muI1 <- c(0,0)
      muI2 <- c(0,0)
      muI3 <- c(1/rnd_func(input$r.mu31),1/rnd_func(input$r.mu32)) #c(0.14,0.14) 1/Dage til død
      
      # recovery rates for different stages
      recI1 <- c(1/rnd_func(input$r.recI11),1/rnd_func(input$r.recI12)) # 1/dage før rask udenfor hospital
      recI2 <- c(1/rnd_func(input$r.recI21),1/rnd_func(input$r.recI22)) # 1/dage på hospital før udskrivelse
      recI3 <- c(1/rnd_func(input$r.recI31),1/rnd_func(input$r.recI32)) # 1/dage på ICU før udskrivelse
      
      # rates going from one state to the next
      gammaEI1 <- c(1/rnd_func(input$r.kE1),1/rnd_func(input$r.kE2)) # rate going from E to I1
      gammaI12 <- c(1/rnd_func(input$r.k1),1/rnd_func(input$r.k2)) # rate going from I1 to I2
      gammaI23 <- c(1/rnd_func(input$r.gI231),1/rnd_func(input$r.gI232)) # rate going from I2 to I3 - move to ICU
      
      # percentage Recover or Moving on
      pI1R <- 1-c(rnd_func(input$r.pI1R1),rnd_func(input$r.pI1R2))*0.01
      pI2R <- c(rnd_func(input$r.pI2R1),rnd_func(input$r.pI2R2)) #1-c(.05,.5)
      pI3R <- c(rnd_func(input$r.pI3R1),rnd_func(input$r.pI3R2))
      
      ## Initial states
      Inf0 <- c(rnd_func(input$r.InfInit1),rnd_func(input$r.InfInit2))
      
      
      S0 <- (pop-Inf0-c(11,8))/Ntot
      HR0 <- c(11,8)/Ntot
      HM0 <- 0*S0
      CR0 <- 0*S0
      CD0 <- 0*S0
      R0 <- 0*S0
      D0 <- 0*S0
      #U0 <- c(ResDK,ResDK)
      U0 <- c(start_day, start_day)
      HC0 <- c(29,30)/Ntot
      CC0 <- c(0,3)/Ntot
      
      E10  <-c(10,19)/24*Inf0/Ntot
      E20  <-c(9,1)/24*Inf0/Ntot
      E30  <-c(0,0)/24*Inf0/Ntot
      
      I1R0 <-c(1,2)/24*Inf0/Ntot*pI1R
      I1M0 <-c(1,2)/24*Inf0/Ntot*(1-pI1R)
      
      I2R0 <-c(2,1)/24*Inf0/Ntot*pI1R
      I2M0 <-c(2,1)/24*Inf0/Ntot*(1-pI1R)
      
      I3R0 <-c(2,1)/24*Inf0/Ntot*pI1R
      I3M0 <-c(2,1)/24*Inf0/Ntot*(1-pI1R)
      
      init0 <- c(S=S0,E1=E10,I1R=I1R0,I1M=I1M0,HR=HR0,HM=HM0,
                 CR=CR0,CD=CD0,R=R0,D=D0,U=U0,HC=HC0,CC=CC0,
                 E2=E20, 
                 E3=E30 ,
                 I2R=I2R0,
                 I3R=I3R0,
                 I2M=I2M0,
                 I3M=I3M0)
      
      p <- c(betaNow=betaNow,betaNew=betaNew,gammaEI1=gammaEI1,gammaI12=gammaI12,gammaI23=gammaI23, 
             muI1=muI1, muI2=muI2, muI3=muI3,
             recI1=recI1, recI2=recI2, recI3=recI3,
             pI1R=pI1R,pI2R=pI2R,pI3R=pI3R,ER=ER) # 15 variables
      
      sol <- ode(init0,times,SSEIRD,p)
      
      if (i==1)
      {
        InfTot <- matrix(NA,nrow=NROW(sol),ncol = nrep)
        InfHos <- matrix(NA,nrow=NROW(sol),ncol = nrep)
        InfInt <- matrix(NA,nrow=NROW(sol),ncol = nrep)
      }
      
      InfTot[,i] <- rowSums(sol[,c(4:17,28:39)])*Ntot
      InfHos[,i] <- rowSums(cbind(sol[,seq(10,12,2)]*Ntot,sol[,seq(11,13,2)]*Ntot))
      InfInt[,i] <- rowSums(cbind(sol[,seq(14,16,2)]*Ntot,sol[,seq(15,17,2)]*Ntot))
      
      if (i==1)
      {
        Sens <- c(p,Inf01=Inf0[1],Inf02=Inf0[2],maxInfInt=max(InfInt[,i]),wmaxInfInt=which.max(InfInt[,i]))
      } else {
        Sens <- rbind(Sens,c(p,Inf01=Inf0[1],Inf02=Inf0[2],maxInfInt=max(InfInt[,i]),wmaxInfInt=which.max(InfInt[,i])))
      }
      
      #incProgress(1/nrep, detail = "Vent venligst")
    }
    
  #})
  
  ####################################################################
  
  parameters <- Sens
  
  colnames(InfInt) <- paste0("Rep", seq_len(nrep))
  colnames(InfHos) <- paste0("Rep", seq_len(nrep))
  colnames(InfTot) <- paste0("Rep", seq_len(nrep))
  
  sim_indlagt_int <- as_tibble(InfInt) %>% 
    mutate(Dato = reference_date + times) %>% 
    pivot_longer(-Dato, names_to = "MC", values_to = "Antal") %>% 
    mutate(MC = factor(MC, levels = colnames(InfInt)))
  
  sim_indlagt_hosp <- as_tibble(InfHos) %>% 
    mutate(Dato = reference_date + times) %>% 
    pivot_longer(-Dato, names_to = "MC", values_to = "Antal") %>% 
    mutate(MC = factor(MC, levels = colnames(InfInt)))
  
  sim_indlagt_tot <- as_tibble(InfTot) %>% 
    mutate(Dato = reference_date + times) %>% 
    pivot_longer(-Dato, names_to = "MC", values_to = "Antal") %>% 
    mutate(MC = factor(MC, levels = colnames(InfInt)))
  
  ans <- list(parameters = parameters,
              sim_indlagt_int = sim_indlagt_int,
              sim_indlagt_hosp = sim_indlagt_hosp,
              sim_indlagt_tot = sim_indlagt_tot)

  return(ans)
}