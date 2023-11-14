#' 2022-02-10
#' Script of the analysis employed to produce the results shown in figure 1 of 
#' https://arxiv.org/abs/2111.09803
#' author: Jacques Hermes <jacques.hermes@fdm.uni-freiburg.de>

# Sets current working directory to the location of this script
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))
library(numDeriv)
library(ggplot2)
library(parallel)
library(ggthemes)
library(magrittr)
library(tidyr)
library(stringr)
library(dMod)
library(ggpubr)
library(lmf)
library(pracma)
library(fBasics)
library(cowplot)

#' Large Function constituting the Augemented Kalman Smoother method 
#' appliable to different simulated datasets 
#'
#'
#' @param sim string specifying the name of the selected simulated datasets
#' @return list of "smooth" (data.frame, output of the smoothed estimates as estimated
#' using the AKS method), "std" (data.frame, the roots of the traces of the smoothed covariance matrices,
#' e.g. the standard deviations of the state-space vector components as 
#' estimated using the AKS method.), "Meas" (data.frame, analyzed data upon 
#' which the estimations of smooth and std are based), "LLold" (numeric, best value of the objective function).

SIRdatasim_para_R <- function(sim){

  # SIRD model with fluxes
  
  In <-function(I,dI,dR,dD){
    I+dI - dR-dD
  }
  Rn <- function(R,dR){
    R+dR
  }
  Dn <- function(D,dD){
    D+dD
  }
  dIn <- function(I,gamma,theta,Rt){
    I * (gamma+theta)*Rt
  }
  dRn <- function(I,gamma){
    gamma * I 
  }
  dDn <- function(I,theta){
    theta * I 
  }
  
  # Data pre-processing
  
  # Selection of the indicated simualted data set from the provide selection of data sets
  if (sim == "sim1"){
    data <- mydata$sim1
  }
  
  if (sim == "sim2"){
    data <- mydata$sim2
  }
  
  if (sim == "sim3"){
    data <- mydata$sim3
  }
  
  if (sim=="sim4"){
    data <- mydata$sim4
  }
  
  if(sim == "sim5"){
    data <- mydata$sim5
  }
  
  if (sim == "sim6"){
    data<- mydata$sim6
  }
  
  if (sim =="sim7"){
    load("SIR_flux_sim_data_Jacques_all_parameters_gamma=002.RData")
    data <- mydata$sim1
  }
  
  # Extraction of the simulated noisy incidence data from the data set
  Meas <- data.frame("v[I]" =rep(NA,length(subset(data, name=="dI")$time)),"v[R]" =rep(NA,length(subset(data, name=="dI")$time)),"v[D]" =rep(NA,length(subset(data, name=="dI")$time)),check.names=FALSE)
  Meas$"v[I]" <- subset(data, name=="dI")$value
  Meas$"v[R]" <- subset(data, name=="dR")$value
  Meas$"v[D]" <- subset(data, name=="dD")$value
  
  time <- 1:length(Meas$"v[I]")
  
  # Number of observations
  NumObs <- length(Meas$"v[I]")
  
  # Observation function
  G <- t(cbind(c(0,0,0,1,0,0,0,0,0),
               c(0,0,0,0,1,0,0,0,0),
               c(0,0,0,0,0,1,0,0,0)))

  # Original guess for state-space vector
  x0 <- c(100,0.1,0.1,1,0.1,0.1,log(0.1),log(0.01),log(6))
  
  # Original guess for state-space covariance matrix
  u0 <- c(mean(Meas$"v[I]"),mean(Meas$"v[R]"),mean(Meas$"v[D]"), mean(Meas$"v[I]"),mean(Meas$"v[R]"),mean(Meas$"v[D]"),0.01+exp(x0[7]),0.01+exp(x0[8]),1)
  P0 <- 1*diag(u0)
  
  # Evolution error guess
  Q <- 0.1*P0
  
  # Observation error guess
  R <-diag(c(mean(Meas$"v[I]"),mean(Meas$"v[R]"),mean(Meas$"v[D]")))
  
  count <- 0
  
  new <- 1
  pnew <- sum(Q)+sum(R)+sum(x0)+sum(P0)
  pold <- 1e20
  epsilon <- 100
  
  system.time(
    while (epsilon>1e-3) {
      
      # Different Lists
      
      # Predicted state-space vectors and covariance matrices
      Pred <- list()
      P_Pred <- list()
      
      # Filtered state-space vectors and covariance matrices
      Est <- data.frame( I=rep(NA,NumObs), R=rep(NA,NumObs),D=rep(NA,NumObs),"v[I]"=rep(NA,NumObs),"v[R]"=rep(NA,NumObs),"v[D]"=rep(NA,NumObs), gamma=rep(NA,NumObs),theta=rep(NA,NumObs),"R[t]"=rep(NA,NumObs),check.names=FALSE)
      Est_list <- list()
      P_Est <- list()
      
      # Smoothed state-space vectors and covariance matrices
      smooth <- data.frame( I=rep(NA,NumObs), R=rep(NA,NumObs),D=rep(NA,NumObs),"v[I]"=rep(NA,NumObs),"v[R]"=rep(NA,NumObs),"v[D]"=rep(NA,NumObs), gamma=rep(NA,NumObs),theta=rep(NA,NumObs),"R[t]"=rep(NA,NumObs),check.names=FALSE)
      std <- data.frame( I=rep(NA,NumObs), R=rep(NA,NumObs),D=rep(NA,NumObs),"v[I]"=rep(NA,NumObs),"v[R]"=rep(NA,NumObs),"v[D]"=rep(NA,NumObs), gamma=rep(NA,NumObs),theta=rep(NA,NumObs),"R[t]"=rep(NA,NumObs),check.names=FALSE)
      smooth_list <- list()
      P_smooth <- list()
      
      # List of the Kalman gains for every time step
      Gain<- list()
      
      # List for fixed lag filter covariance matrices needed for the M-step of the algorithm
      P_lag <- list()
      
      # List the jacobians for every time step
      Jacob <- list()
      
      #' Prediction function giving the state-space vector at time k as 
      #' predicted by the SIRD model
      #'
      #' @param x vector, smoothed state-space vector at time point k-1
      #' @return vector, predicted state-space vector for time point k
      
      b <- function(x){
        c(In(x[1],x[4],x[5],x[6]),
          Rn(x[2],x[5]),
          Dn(x[3],x[6]),
          dIn(x[1],exp(x[7]),exp(x[8]),exp(x[9])),
          dRn(x[1],exp(x[7])),
          dDn(x[1],exp(x[8])),
          x[7],
          x[8],
          x[9])
      }
      
      #' Jacobian function, giving the jacobian of the prediction function for 
      #' the smoothed estimate of the state-space vector at time point k-1.
      #'This matrix was computed by hand in order to facilitate
      #' the numerical calculations
      #'
      #' @param x vector, smoothed state-space vector at time point k-1
      #' @return matrix, jacobian of the prediction function for 
      #' the smoothed estimate of the state-space vector at time point k-1.
      
      jac <- function(x){
        rbind(
          c(1,0,0,1,-1,-1,0,0,0),
          c(0,1,0,0,1,0,0,0,0),
          c(0,0,1,0,0,1,0,0,0),
          c((exp(x[8])+exp(x[7]))*(exp(x[9])),0,0,0,0,0,exp(x[7])*x[1]*(exp(x[9])),exp(x[8])*x[1]*(exp(x[9])),x[1]*(exp(x[8])+exp(x[7]))*exp(x[9])),
          c(exp(x[7]),0,0,0,0,0,exp(x[7])*x[1],0,0),
          c(exp(x[8]),0,0,0,0,0,0,exp(x[8])*x[1],0),
          c(0,0,0,0,0,0,1,0,0),
          c(0,0,0,0,0,0,0,1,0),
          c(0,0,0,0,0,0,0,0,1)
        )
      }
      
      # Augmented Kalman Filter
      P<-P0
      x<-x0
      
      for(i in 1:NumObs){
        # Prediction Step
        A <- jac(x )  
        Jacob[[i]] <- A
        x_k <- b(x )
        
        P_k <-  A %*% P%*% t(A) + Q
        
        # Cholesky decomposition ensures positive definiteness of cov. matx
        P_k <- t(chol(P_k))%*%chol(P_k)
        
        Pred[[i]] <- x_k
        P_Pred[[i]] <- P_k
        
        # Observation
        y <- c(Meas$"v[I]"[i],Meas$"v[R]"[i],Meas$"v[D]"[i])
        
        # Measurement Update Step  
        K_k <- P_k %*% t(G) %*% chol2inv(chol(G %*% P_k %*% t(G) + R))
        Gain[[i]] <- K_k
        
        x_k1 <- x_k +   K_k %*% (y - G %*% x_k)
        
        # Joseph-form covariance update
        P_k1 <- (diag(1,ncol(G))-K_k%*%G)%*%P_k%*%t((diag(1,ncol(G))-K_k%*%G)) + K_k%*%R%*%t(K_k)
        
        P_k1 <- t(chol(P_k1))%*%chol(P_k1)
        
        # Next iteration step uses filtered estimates as input
        P <- P_k1
        x <- x_k1
        Est_list[[i]] <- x
        Est[i,] <- x
        P_Est[[i]]<-P
        
      }
      
      # Kalman Smoother is initialized with the last filtered estimate for the 
      # state-space vector and covariance matrix
      smooth[1,] <- tail(Est,n=1)
      smooth_list[[1]] <- tail(Est_list,n=1)[[1]]
      P_smooth[[1]] <- tail(P_Est,n=1)[[1]]
      std[1,] <- sqrt(diag(tail(P_Est,n=1)[[1]]))
      B_list <- list()
      
      # Augmented Kalman Smoother
      for(i in 2:NumObs){
        
        B <- tail(P_Est,n=i)[[1]] %*% t(tail(Jacob, n=i-1)[[1]]) %*% chol2inv(chol(tail(P_Pred, n=(i-1))[[1]]))
        
        B_list[[i-1]] <- B
        x_s <- tail(Est_list, n=i)[[1]] + B %*% (smooth_list[[i-1]]-tail(Pred, n=(i-1))[[1]])
        
        P_s <- tail(P_Est, n=i)[[1]] + B %*% (P_smooth[[i-1]]-tail(P_Pred, n=(i-1))[[1]]) %*% t(B)
        P_s <- t(chol(P_s))%*%chol(P_s)
        
        smooth_list[[i]] <- x_s
        
        smooth[i,] <- x_s
        std[i,] <- sqrt(diag(P_s))
        P_smooth[[i]] <- P_s
        
      }
      
      # Fixed lag smoother for M-step
      P_lag[[1]] <- (diag(1,ncol(G))-tail(Gain,n=1)[[1]]%*%G)%*% tail(Jacob,n=1)[[1]] %*% tail(P_Est, n=2)[[1]]
      for (i in 2:(NumObs-1)){
        P_lag[[i]] <- tail(P_Est, n=i)[[1]] %*% t(B_list[[i]]) + B_list[[i-1]] %*%(P_lag[[i-1]]- tail(Jacob, n=i-1)[[1]]%*%tail(P_Est, n=i)[[1]]) %*%t(B_list[[i]])
      }
      
      # Change order of smoother lists 
      smooth <- smooth[nrow(smooth):1,]
      std <- std[nrow(std):1,]
      smooth_list <- smooth_list[length(smooth_list):1]
      
      P_smooth <- P_smooth[length(P_smooth):1]
      P_lag <- P_lag[length(P_lag):1]
      
      time <- c(1:NumObs)
      
      # Prediction of smoothed x by b
      Pred_smooth <- list()
      Pred_smooth <- lapply(smooth_list, b)
      
      #Calculations needed for M-step
      summ <- 0
      for (i in time[-1]){ summ<- summ+(c(Meas$"v[I]"[i],Meas$"v[R]"[i],Meas$"v[D]"[i]) - G%*%smooth_list[[i]])%*%
        t(c(Meas$"v[I]"[i],Meas$"v[R]"[i],Meas$"v[D]"[i]) - G%*%smooth_list[[i]]) + G%*%P_smooth[[i]]%*%t(G)}
      
      sig <- 0 
      for (i in time[-1]){sig <- sig +P_smooth[[i]]+((smooth_list[[i]]-Pred_smooth[[i-1]])%*%t((smooth_list[[i]]-Pred_smooth[[i-1]])))+
        Jacob[[i-1]]%*%P_smooth[[i-1]]%*%t(Jacob[[i-1]])-P_lag[[i-1]]%*%t(Jacob[[i-1]])-Jacob[[i-1]]%*%t(P_lag[[i-1]])}
      
      # M-step
      Rold<-R
      R <- 1/NumObs *summ
      R <- t(chol(R)) %*% chol(R)
      
      Qold<- Q
      Q<- 1/NumObs * sig
      Q <- t(chol(Q)) %*% chol(Q)
      
      # Set initial values of next EM iteration to smoothed estimates of x and
      # P at time point 0
      x0 <- smooth_list[[1]]
      P0 <- P_smooth[[1]]
      
      count<-count+1
      
      Meas$time <- time
      
      smooth$time <- Meas$time
      std$time <- Meas$time
      
      # Stopping criterion
      pold <- pnew
      pnew <- sum(Q)+sum(R)+sum(x0)+sum(P0)
      epsilon<- abs(pold-pnew)/abs(pold)
    })
 
  return(list(smooth,std,Meas))
  
}
result <- list()
error <- list()
meas <- list()
Likelihood <- list()
c <- 1

# Function call for different data sets
for (i in c("sim1","sim2","sim3", "sim4", "sim5",
            "sim6","sim7"
            )){
  
  # Loading of data generated using an SIRD Model with different parameter time
  # courses
  load("SIR_flux_sim_data_Jacques_all_gamma=002.RData")
  dummy <- SIRdatasim_para_R(i)
  result[[c]] <- dummy[[1]]
  error[[c]] <- dummy[[2]]
  meas[[c]] <- dummy[[3]]
  
  c <- c+1
  print(paste(i, "Done"))
}

save(result, file = "Result_simulations_SIRD.RData")
save(error, file = "Error_simulations_SIRD.RData")
save(meas, file = "Meas_simulations_SIRD.RData")