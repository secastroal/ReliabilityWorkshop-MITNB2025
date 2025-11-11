# Simulate data based on the 2RDM.

# In this script, we simulate data based on the 2RDM model when there is only 
# one common factor at the between- and within-level. For example, when all the
# the items are supposed to measure a construct such as positive affect.

# Set up ----
N  <- 100 # Number of participants
nT <- 90  # Number of time points
I  <- 6   # Number of items

set.seed(1234)

# Between-level parameters ----

# Intercepts
mu <- rep(0, I) 
# Variance of the between-level latent factor
var_xi <- 2 
# Variance of the between-level measurement error
var_theta <- rep(0.5, I) 
# Between-level factor loadings
bloadings <- rep(1, I) 

# True between-person Reliability
true_brel <- ((sum(bloadings)^2) * var_xi)/ 
  ((sum(bloadings)^2) * var_xi + sum(var_theta))

# Within-level parameters ----

# Mean of the within-level random factor loadings
mean_wload <- rep(c(1, 1.2, 0.8), length = I) 
# Variance of the within-level random factor loadings
var_wload  <- c(0, rlnorm(I - 1, -3, 0.4)) 
# Mean of the within-level measurement error
mean_werr <- rlnorm(I, 0.15, 0.2)
# Variance of the within-level measurement error
var_werr <- rlnorm(I, -1, 0.4)
# Mean of the autoregressive effect
mean_phi <- 0.5
# Variance of the autoregressive effect
var_phi <- 0.05
# Mean of the innovations
mean_inno <- 0.8
# Variance of the innovations
var_inno <- 0.3

# Generate random effects ----

# Random loadings
wload <- matrix(NA, nrow = N, ncol = I)

for (i in 1:I) {
  wload[, i] <- rnorm(N, mean_wload[i], sqrt(var_wload[i]))
}
rm(i)

# Random measurement error variances
werr <- matrix(NA, nrow = N, ncol = I)

for (i in 1:I) {
  location <- log(mean_werr[i]^2 / sqrt(var_werr[i] + mean_werr[i]^2))
  shape    <- sqrt(log(1 + (var_werr[i] / mean_werr[i]^2)))
  werr[, i] <- rlnorm(N, location, shape)
}
rm(i, location, shape)

# Random Autoregressive effect
phi <- rnorm(N, mean_phi, var_phi)

# Random innovation variances
location <- log(mean_inno^2 / sqrt(var_inno + mean_inno^2))
shape    <- sqrt(log(1 + (var_inno / mean_inno^2)))
inno     <- rlnorm(N, location, shape)
rm(location, shape)

# Compute True within-level reliability coefficients ----

# Observed Variance 
ovar   <- inno/(1 - phi^2)
# Squared sum of loadings
sqload <- rowSums(wload)^2
# sqload <- rowSums(wload^2)
# Sum of measurement error
ersum  <- rowSums(werr)

# Compute within person reliability
true_wrel <- (sqload * ovar)/(sqload * ovar + ersum)

# Generate Factor Scores ----

# Between-person factor scores
xi <- rnorm(N, 0, var_xi)

# Within-person factor scores
zeta <- c()

for (p in 1:N) {
  tmp <- c()
  tmp[1] <- rnorm(1, 0, sqrt(inno[p]/(1 - phi[p]^2)))
  for (i in 2:(2*nT)) {
    tmp[i] <- phi[p] * tmp[i-1] + rnorm(1, 0, sqrt(inno[p]))
  }
  tmp <- tail(tmp, nT)
  zeta <- c(zeta, tmp)
  rm(tmp, i)
}
rm(p)

# Compute observed scores ----

# Between measurement error
b_error <- matrix(NA, N, I)

for (i in 1:I) {
  b_error[, i] <- rnorm(N, 0, var_theta[i])
}
rm(i)

# Between scores
bscore <- matrix(mu, nrow = N, ncol = I, byrow = TRUE) + 
  xi %*% t(bloadings) + b_error

# Between scores expanded

bscore <- bscore[rep(1:N, each = nT), ]

# Within measurement error

w_error <- c()

for (p in 1:N) {
  
  tmp <- matrix(NA, nT, I)
  
  for (i in 1:I) {
    tmp[, i] <- rnorm(nT, 0, werr[p, i])
  }
  
  w_error <- rbind(w_error, tmp)
  
  rm(tmp)
}
rm(p, i)

# Within "true scores"

w_truescores <- c()

for (p in 1:N) {
  
  tmp <- zeta[(((p - 1) * nT) + 1):(p * nT)] %*% t(wload[p,])
  
  w_truescores <- rbind(w_truescores, tmp)
  
  rm(tmp)
}
rm(p)

# Within scores

wscore <- w_truescores + w_error

# Observed scores

data2rdm <- bscore + wscore

# Turn data to data.frame, add id and time variables ----

# Turn to data.frame
data2rdm <- data.frame(data2rdm)
names(data2rdm) <- paste0("I", 1:I)

# ID variable
id <- rep(1:N, each = nT)

# Time variable
time <- rep(1:nT, times = N)

# Create complete data
data2rdm <- data.frame(id, time, data2rdm)

# Save simulated data and true reliabilities ----

truepars <- list(
  true_brel = true_brel, 
  true_wrel = true_wrel,
  true_pars = list(
    between_loadings = bloadings,
    mean_wload = mean_wload,
    mean_phi = mean_phi,
    mean_werr = mean_werr,
    mean_inno = mean_inno,
    intercepts = mu,
    var_bf = var_xi,
    var_wload = var_wload,
    var_phi = var_phi,
    var_werr = var_werr,
    var_inno = var_inno,
    var_theta = var_theta
  )) 

write.table(data2rdm, file = "Data/data2rdm.dat")
saveRDS(truepars, file = "Data/truepars.rds")



