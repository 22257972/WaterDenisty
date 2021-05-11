
library(dplyr)
library(reshape2)

# load("Crux_KP150_Phs1.RData")
# load("Crux_KP150_Phs2.RData")

fulldata <- merge(Crux_KP150_Phs1$full_dat, select(Crux_KP150_Phs2$full_dat, -"-11.7"), by = intersect(names(Crux_KP150_Phs1$full_dat), names(Crux_KP150_Phs2$full_dat)), all = TRUE)
n_depths <- NCOL(fulldata) - 1
n_times <- NROW(fulldata)

Burnin <- 1000
Samples <- 10000
Sampled_Range <- (Burnin+1):(Burnin+Samples)
Total_Samples <- Burnin + Samples
n_para <- 6

Beta_All <- array(data = NA, dim = c(n_para, n_times, Burnin + Samples))
Sigma_Squared_All <- rep(NA, Burnin + Samples)
Mu_All <- matrix(data = NA, nrow = n_para, ncol = Burnin + Samples)
Tau_Squared_All <- matrix(data = NA, nrow = n_para, ncol = Burnin + Samples)

X <- matrix(data = 1, nrow = n_depths, ncol = n_para)
depths <- as.numeric(names(fulldata[,-1]))
X[,2] <- (depths-mean(depths))/sd(depths)

for (i in 2:(n_para-1)) {
  X[,i + 1] <- X[,2]^i
}

Y <- matrix(data = data.matrix(fulldata[1:n_times,-1]), nrow = n_depths, ncol = n_times, byrow = TRUE)

Beta <- matrix(data = NA, nrow = n_para, ncol = n_times)

a <- 0.001
b <- 0.001
Sigma_Squared <- runif(1, 1, 100)

Mu_Means <- rep(0, n_para)
VCOV_Matrix_Mu <- diag(x = 10^4, nrow = n_para, ncol = n_para)
Mu <- t(t(rnorm(n_para, sd = 10)))

c <- 0.001
d <- 0.001
Tau_Squared <- runif(n_para, 1, 10^4)

X_T <- t(X)
X_T_Times_X <- X_T %*% X
X_T_Times_Y <- X_T %*% Y

# Simulations to improve loop efficiency:

Standard_Normal_Simulation_For_Betas <- array(data = rnorm(n = n_times * (Burnin + Samples) * n_para), dim = c(n_para, n_times, Burnin + Samples))
Gamma_Simulation_For_Sigma_Squared <- rgamma(n = Burnin + Samples, shape = (n_times * n_depths)/2 + a, rate = 1)
Standard_Normal_Simulation_For_Mus <- matrix(data = rnorm(n = n_para * (Burnin + Samples)), nrow = n_para, ncol = Burnin + Samples)
Gamma_Simulation_For_Tau_Squared <- matrix(data = rgamma(n = (Burnin + Samples) * n_para, shape = n_times/2 + c, rate = 1), nrow = n_para, ncol = Burnin + Samples)

# MCMC Loop:

timer <- Sys.time()

for (i in 1:(Burnin + Samples)) {
  
  # Beta's
  VCOV_Matrix_Inverse <- diag(1/Tau_Squared)
  V <- solve(X_T_Times_X/Sigma_Squared + VCOV_Matrix_Inverse)
  V_chol <- t(chol(V))
  SSE <- 0
  for (t in 1:n_times) {
    M <- V %*% (X_T_Times_Y[,t]/Sigma_Squared + VCOV_Matrix_Inverse %*% Mu)
    Beta[,t] <- M + V_chol %*% Standard_Normal_Simulation_For_Betas[,t,i]
    SSE <- SSE + as.numeric(t(Y[,t] - (X %*% Beta[,t])) %*% (Y[,t] - (X %*% Beta[,t])))
  }
  
  # Sigma Squared
  Sigma_Squared <- (SSE/2 + b)/Gamma_Simulation_For_Sigma_Squared[i]
  
  # Mu's
  for (j in 1:n_para) {
    Beta_Sum <- sum(Beta[j,])
    Mu_Mean_Temp <- (Tau_Squared[j] * Mu_Means[j]+VCOV_Matrix_Mu[j,j]*Beta_Sum)/(Tau_Squared[j]+ n_times * VCOV_Matrix_Mu[j,j])
    Mu_sd_Temp <- sqrt((Tau_Squared[j] * VCOV_Matrix_Mu[j,j])/(Tau_Squared[j] + n_times * VCOV_Matrix_Mu[j,j]))
    Mu[j] <- Mu_Mean_Temp + Mu_sd_Temp * Standard_Normal_Simulation_For_Mus[j,i]
  }
  
  # Tau Squared's
  for (j in 1:n_para) {
    SSE_Mu <- 0
    for (t in 1:n_times) {
      SSE_Mu <- SSE_Mu + (Beta[j,t] - Mu[j])^2
    }
    Tau_Squared[j] <- (SSE_Mu/2 + d)/Gamma_Simulation_For_Tau_Squared[j,i]
  }
  
  # Store Data
  Beta_All[,,i] <- Beta
  Sigma_Squared_All[i] <- Sigma_Squared
  Mu_All[,i] <- Mu
  Tau_Squared_All[,i] <- Tau_Squared
  
  if(i%%1000 == 0){
    cat("Iteration:", i, "/", Total_Samples, "Time Taken:", difftime(Sys.time(), timer, units='mins'), "Approx Time Remaining:", (difftime(Sys.time(), timer, units='mins'))/i * (Total_Samples - i),"Approx Total Time:", difftime(Sys.time(), timer, units='mins') + (difftime(Sys.time(), timer, units='mins'))/i * (Total_Samples - i),"\n")
  }
  
} 

# Plot chains

Time_Chosen <- 1
par(mfrow=c(4,n_para))
for (i in 1:n_para) {
  plot(Beta_All[i, Time_Chosen, Sampled_Range], type = "l", col = "red", main = paste0("beta_", i-1, " for time = ", Time_Chosen))
}
for (i in 1:n_para) {
  plot(Mu_All[i, Sampled_Range], type = "l", col = "blue", main = paste0("mu_", i-1))
}
for (i in 1:n_para) {
  plot(Tau_Squared_All[i,Sampled_Range], type = "l", col = "purple", main = paste0("tau_squared_", i-1))
}
plot(Sigma_Squared_All[Sampled_Range], type = "l", main = "sigma_squared")
