
load("crux_kp150_phs1.Rdata")
load("crux_kp150_phs2.Rdata")

library(dplyr)
library(reshape2)

fulldata <- merge(Crux_KP150_Phs1$full_dat, select(Crux_KP150_Phs2$full_dat, -"-11.7"), by = intersect(names(Crux_KP150_Phs1$full_dat), names(Crux_KP150_Phs2$full_dat)), all = TRUE)
n_depths <- NCOL(fulldata) - 1
n_times <- NROW(fulldata)
n_chains <- 4

Burnin <- 50000
Samples <- 50000
Thinning <- 5
Total_Sampled <- Burnin + Samples
Sampled_Post_Thin <- Samples/Thinning
Sampled_Range <- (Burnin+1):(Burnin+Samples)
Sampled_Range_Thin <- seq(from = 1, length.out = Samples %/% Thinning)
n_para <- 6

X <- matrix(as.numeric(names(fulldata[,-1])), nrow = n_depths, ncol = 1)
X_T <- t(X)
Y <- matrix(data = data.matrix(fulldata[1:n_times,-1]), nrow = n_depths, ncol = n_times, byrow = TRUE)

DHT <- function(X, Betas){
  tanh((X + Betas[1])/Betas[2]) + tanh((X + Betas[3])/Betas[4])
}

beta_ll_optim <- function(beta, y, alpha_curr, sigma2_curr) {
  f <- cbind(1, DHT(X, beta))
  -(sum(dnorm(y, mean = f %*% alpha_curr, sd = sqrt(sigma2_curr), log = TRUE)) +
      sum(dnorm(beta, c(75, 80, 150, 80), rep(15, 4), log = TRUE)))
}

beta_mode <- matrix(NA, nrow = n_times, ncol = 4)
chol_prop_mat <- array(NA, dim = c(4, 4, n_times))

for (i in 1:n_times) {
  beta_mode_single <- optim(
  c(75, 80, 150, 80),
  beta_ll_optim,
  y = Y[,i],
  alpha_curr = c(1025, -5),
  sigma2_curr = 0.01,
  hessian = TRUE
  )
  beta_mode[i, ] <- beta_mode_single$par
  
  tryCatch(
  chol_prop_mat[, , i] <- t(chol(solve(beta_mode_single$hessian))),
  error = function(e) {
    print(e)
  }
)
if (any(is.na(chol_prop_mat[, , i]))) {
chol_prop_mat[, , i] <- chol_prop_mat[, , i - 1]
}
}

Beta_Log_Posterior <- function(Y_Subtract_DHT, Beta, Mu_Beta, VCOV_Beta_Inverse, Sigma_Squared){
as.numeric(-1/(2 * Sigma_Squared) * (t(Y_Subtract_DHT) %*% Y_Subtract_DHT) - 0.5 * t(Beta - Mu_Beta) %*% VCOV_Beta_Inverse %*% (Beta - Mu_Beta))
}

timer <- Sys.time()

#Set up storage
Alpha_All <- array(data = NA, dim = c(2, n_times, Sampled_Post_Thin))
Beta_All <- array(data = NA, dim = c(4, n_times, Sampled_Post_Thin))

Sigma_Squared_All <- rep(NA, Sampled_Post_Thin)

Mu_Alpha_All <- matrix(data = NA, nrow = 2, ncol = Sampled_Post_Thin)
Mu_Beta_All <- matrix(data = NA, nrow = 4, ncol = Sampled_Post_Thin)

Tau_Squared_Alpha_All <- matrix(data = NA, nrow = 2, ncol = Sampled_Post_Thin)
Tau_Squared_Beta_All <- matrix(data = NA, nrow = 4, ncol = Sampled_Post_Thin)

# Initialise model parameters and starting values
Sigma_Squared <- runif(1,0.1,5)
a <- 0.001
b <- 0.001

Mu_Alpha_Means <- c(1025, -5)
VCOV_Mu_Alpha <- diag(x = c(10^2, 2^2), nrow = 2, ncol = 2)
Mu_Alpha <- as.numeric(mvtnorm::rmvnorm(1, Mu_Alpha_Means, VCOV_Mu_Alpha))

Mu_Beta_Means <- c(75, 80, 150, 80)
VCOV_Mu_Beta <- diag(x = c(15^2, 15^2, 15^2, 15^2), nrow = 4, ncol = 4)
Mu_Beta <- as.numeric(mvtnorm::rmvnorm(1, Mu_Beta_Means, VCOV_Mu_Beta))

c <- 1
d <- 1
Tau_Squared_Alpha <- runif(2, 1, 4)
Tau_Squared_Beta <- runif(4, 5, 20)

Alpha <- matrix(mvtnorm::rmvnorm(n_times, Mu_Alpha, diag(Tau_Squared_Alpha)), nrow = 2, ncol = n_times, byrow = TRUE)
Beta <- matrix(mvtnorm::rmvnorm(n_times, Mu_Beta, diag(Tau_Squared_Beta)), nrow = 4, ncol = n_times, byrow = TRUE)
Beta_Propose <- Beta

# Generate random numbers for loop
Gamma_Simulation_For_Sigma_Squared <- rgamma(n = Total_Sampled, shape = (n_times * n_depths)/2 + a, rate = 1)

Standard_Normal_Simulation_For_Mu_Alpha <- matrix(data = rnorm(n = 2 * Total_Sampled), nrow = 2, ncol = Total_Sampled)
Standard_Normal_Simulation_For_Mu_Beta <- matrix(data = rnorm(n = 4 * Total_Sampled), nrow = 4, ncol = Total_Sampled)

Gamma_Simulation_For_Tau_Squared_Alpha <- matrix(data = rgamma(n = Total_Sampled * 2, shape = n_times/2 + c, rate = 1), nrow = 2, ncol = Total_Sampled)
Gamma_Simulation_For_Tau_Squared_Beta <- matrix(data = rgamma(n = Total_Sampled * 4, shape = n_times/2 + c, rate = 1), nrow = 4, ncol = Total_Sampled)

timer <- Sys.time()

for (i in 1:Total_Sampled) {
  
  # Alpha's & Beta's
  VCOV_Alpha_Inverse <- diag(1/Tau_Squared_Alpha)
  VCOV_Beta_Inverse <- diag(1/Tau_Squared_Beta)
  SSE <- 0
  
  for (t in 1:n_times) {
    # Alpha's
    DHT_Matrix_Current <- matrix(data = c(rep(1, n_depths), DHT(X, Beta[,t])), nrow = n_depths, ncol = 2)
    V <- solve(t(DHT_Matrix_Current) %*% DHT_Matrix_Current / Sigma_Squared + VCOV_Alpha_Inverse)
    V_chol <- t(chol(V))
    M <- V %*% (t(DHT_Matrix_Current) %*% Y[,t] / Sigma_Squared + VCOV_Alpha_Inverse %*% Mu_Alpha)
    Alpha[,t] <- M + V_chol %*% rnorm(2)
    
    # Beta's
    Beta_Propose[,t] <- Beta[,t] + chol_prop_mat[,,t] %*% rnorm(4, sd = 1)/4
    
    Y_Subtract_DHT_Current <- Y[,t] - (DHT_Matrix_Current %*% Alpha[,t,drop=FALSE])
    Beta_Density_Current <- Beta_Log_Posterior(Y_Subtract_DHT_Current, Beta[,t], Mu_Beta, VCOV_Beta_Inverse, Sigma_Squared)
    
    DHT_Matrix_Propose <- matrix(data = c(rep(1, n_depths), DHT(X, Beta_Propose[,t])), nrow = n_depths, ncol = 2)
    Y_Subtract_DHT_Propose <- Y[,t] - (DHT_Matrix_Propose %*% Alpha[,t,drop=FALSE])
    Beta_Density_Propose <- Beta_Log_Posterior(Y_Subtract_DHT_Propose, Beta_Propose[,t], Mu_Beta, VCOV_Beta_Inverse, Sigma_Squared)
    
    if(log(runif(1)) < min(0, Beta_Density_Propose - Beta_Density_Current)){
      Beta[,t] <- Beta_Propose[,t]
      Y_Subtract_DHT <- Y_Subtract_DHT_Propose
    } else {
      Y_Subtract_DHT <- Y_Subtract_DHT_Current
    }
    
    if(Beta[1,t] > Beta[3,t]){
      Temp <- Beta[1,t]
      Beta[1,t] <- Beta[3,t]
      Beta[3,t] <- Temp
    }
    
    SSE <- SSE + as.numeric(t(Y_Subtract_DHT) %*% Y_Subtract_DHT)
  }
  
  # Sigma Squared
  Sigma_Squared <- (SSE/2 + b)/Gamma_Simulation_For_Sigma_Squared[i]
  
  # Mu's
  for (j in 1:2) {
    Alpha_Sum <- sum(Alpha[j,])
    Mu_Mean_Temp <- (Tau_Squared_Alpha[j] * Mu_Alpha_Means[j]+VCOV_Mu_Alpha[j,j]*Alpha_Sum)/(Tau_Squared_Alpha[j] + n_times * VCOV_Mu_Alpha[j,j])
    Mu_sd_Temp <- sqrt((Tau_Squared_Alpha[j] * VCOV_Mu_Alpha[j,j])/(Tau_Squared_Alpha[j] + n_times * VCOV_Mu_Alpha[j,j]))
    Mu_Alpha[j] <- Mu_Mean_Temp + Mu_sd_Temp * Standard_Normal_Simulation_For_Mu_Alpha[j,i]
  }
  
  for (j in 1:4) {
    Beta_Sum <- sum(Beta[j,])
    Mu_Mean_Temp <- (Tau_Squared_Beta[j] * Mu_Beta_Means[j]+VCOV_Mu_Beta[j,j]*Beta_Sum)/(Tau_Squared_Beta[j] + n_times * VCOV_Mu_Beta[j,j])
    Mu_sd_Temp <- sqrt((Tau_Squared_Beta[j] * VCOV_Mu_Beta[j,j])/(Tau_Squared_Beta[j] + n_times * VCOV_Mu_Beta[j,j]))
    Mu_Beta[j] <- Mu_Mean_Temp + Mu_sd_Temp * Standard_Normal_Simulation_For_Mu_Beta[j,i]
  }
  
  # Tau Squared's
  for (j in 1:2) {
   SSE_Mu <- 0
  for (t in 1:n_times) {
    SSE_Mu <- SSE_Mu + (Alpha[j,t] - Mu_Alpha[j])^2
  }
    Tau_Squared_Alpha[j] <- (SSE_Mu/2 + d)/Gamma_Simulation_For_Tau_Squared_Alpha[j,i]
  }
  
  for (j in 1:4) {
    SSE_Mu <- 0
  for (t in 1:n_times) {
    SSE_Mu <- SSE_Mu + (Beta[j,t] - Mu_Beta[j])^2
  }
    Tau_Squared_Beta[j] <- (SSE_Mu/2 + d)/Gamma_Simulation_For_Tau_Squared_Beta[j,i]
  }
  
  # Store Data
  output <- list()
  if(i > Burnin && i%%Thinning == 0){
    
    index <- (i - Burnin)/Thinning
    
    Alpha_All[,,index] <- Alpha
    Beta_All[,,index] <- Beta
    
    Sigma_Squared_All[index] <- Sigma_Squared
    
    Mu_Alpha_All[,index] <- Mu_Alpha
    Mu_Beta_All[,index] <- Mu_Beta
    
    Tau_Squared_Alpha_All[,index] <- Tau_Squared_Alpha
    Tau_Squared_Beta_All[,index] <- Tau_Squared_Beta
  }
  
  if(i%%5000 == 0){
  cat("Iteration:", i, "/", Total_Sampled, "Time:", difftime(Sys.time(), timer, units='mins'), "Approx Time Remaining:", (difftime(Sys.time(), timer, units='mins'))/i * (Total_Sampled - i),"Approx Total Time:", difftime(Sys.time(), timer, units='mins') + (difftime(Sys.time(), timer, units='mins'))/i * (Total_Sampled - i),"\n")
  }
} 

par(mfrow = c(4,6))
Time_Chosen <- 1

for (i in 1:2) {
  matplot(Alpha_All[i,Time_Chosen,], type = "l", 
          main = paste0("alpha_", i-1, " for time = ", Time_Chosen), lty = 1)
}

for (i in 1:4) {
  matplot(Beta_All[i,Time_Chosen,], type = "l", main = paste0("beta_", i-1, " for time = ", Time_Chosen), lty = 1)
}

for (i in 1:2) {
  matplot(Mu_Alpha_All[i,], type = "l", main = paste0("mu_alpha_", i-1), lty = 1)
}

for (i in 1:4) {
  matplot(Mu_Beta_All[i,], type = "l", main = paste0("mu_beta_", i-1), lty = 1)
}

for (i in 1:2) {
  matplot(Tau_Squared_Alpha_All[i,], type = "l", main = paste0("tau_squared_alpha_", i-1), lty = 1)
}

for (i in 1:4) {
  matplot(Tau_Squared_Beta_All[i,], type = "l", main = paste0("tau_squared_beta_", i-1), lty = 1)
}

matplot(Sigma_Squared_All, type = "l", main = paste0("sigma_squared_", 0), lty = 1)

