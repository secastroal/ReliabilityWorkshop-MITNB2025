# Computing the internal consistency reliability

# This file shows how to estimate the internal consistency reliability of 
# the simulated dataset "data2rdm" based on six/seven different approaches. 

# Set up ----

# Install required packages
packages <- c("psych", "misty", "lavaan", "dynr", "MplusAutomation")
install.packages(setdiff(packages, rownames(installed.packages())))
rm(packages)

# Load custom functions
custom_functions <- list.files("R/", full.names = TRUE)
sapply(custom_functions, source)
rm(custom_functions)

# Load data ----

# Load data
data2rdm <- read.table("Data/data2rdm.dat", header = TRUE)
# Load true parameters
truepars <- readRDS("Data/truepars.rds")

# Define environment variables needed later.
I  <- 6   # Number of items
N  <- 100 # Number of persons
nT <- 90  # Number of timepoints

# Store items' names to easily select them later.
items <- paste0("I", 1:I)

# GT and MLM ----

# The reliability based on the GT and MLM can be easily computed with the 
# function mlr or multilevel.reliability of the psych package.

t0 <- proc.time()
mlr_data2rdm <- psych::mlr(
  data2rdm, 
  grp   = "id",
  Time  = "time",
  items = items,
  lmer  = TRUE,
  lme   = FALSE,
  alpha = FALSE,
  aov   = FALSE
  )
t1 <- proc.time() - t0

t1[3]

# GT Reliability estimates
mlr_data2rdm$RkR # between-person
mlr_data2rdm$Rc  # within-person

# MLM Reliability estimates
mlr_data2rdm$RkRn # between-person
mlr_data2rdm$Rcn  # within-person

# Multilevel omega ----

# To estimate the reliability based on the multilevel CFA, one can easily use 
# the function multilevel.omega from the misty package. 

mlcfa_data2rdm <- misty::multilevel.omega(
  data2rdm[, items],
  cluster = data2rdm$id
)

# Between and within omega
mlcfa_data2rdm$result$omega

# P-technique factor analysis ----

# The P-technique implies fitting a factor model to the responses of each
# person. To do this, we wrote a custom function that automatically loops
# through each person to fit the model in lavaan. 
# The arguments of this function are the data, the label of the "id" variable,
# the names of the items, and the how the model is identified. See source for
# further details.

#!# Make sure your data is a data.frame!

pfa_data2rdm <- esm_pfa_omega(data  = data2rdm, 
                              idlab = "id", 
                              items = items, 
                              constraint = "variance")

# Summary and histogram of estimated reliabilities
summary(pfa_data2rdm$pfa_omega$pfa_omega)
hist(pfa_data2rdm$pfa_omega$pfa_omega, xlim = c(0, 1))

# Dynamic factor analysis in dynr ----

# Similarly as with the P-technique, dynamic factor analysis implies fitting a 
# factor model to the responses of each person, but in this case, the model can 
# include autoregressive effects.
# To do this, we wrote a custom function that automatically loops
# through each person to fit a lag-1 DFA in dynr. 
# The arguments of this function are the data, the label of the "id" variable,
# the label of the "time" variable, and the names of the items. See source for
# further details.

#!# Make sure your data is a data.frame!

dfa_data2rdm <- esm_dfa_omega(data    = data2rdm,
                              idlab   = "id",
                              timelab = "time",
                              items   = items)

# Summary and histogram of estimated reliabilities 
summary(dfa_data2rdm$dfa_omega$dfa_omega)
hist(dfa_data2rdm$dfa_omega$dfa_omega, xlim = c(0, 1))

# Export data to Mplus with Mplus Automation ----
# For the remaining two approaches, we use Mplus. 
# Therefore, we export the data to Mplus.

#!# This code should be run only once. It takes your data object
#!# and exports it to the path defined in filename. With this function, 
#!# missing values are represented as dots(.) in the exported data file.
#!# To avoid overwriting the file over and over, the if function first check
#!# whether the data have been exported before. If the exported data exists 
#!# already, the code is ignored. 

#!# If you make changes to your data (e.g., removing participants), delete 
#!# the existing exported data from your folder, and run the following code 
#!# again to export the updated data.  

if (!file.exists("Mplus/data2rdm.dat")) {
  MplusAutomation::prepareMplusData(
    data2rdm,
    filename = "Mplus/data2rdm.dat",
    inpfile = FALSE
  )
}

# 2RDM ----
# The 2RDM is fitted in Mplus. This can take long depending on the amount of 
# data, the number of iterations, and other factors. This analysis took 
# approximately 20 minutes.

#!# Here, we check if an output file from running the model already exists.
#!# If it does not, the model if run in Mplus from within R. However,
#!# it might be better to simply run the model directly in Mplus, because 
#!# when running in Mplus, the program shows the progress of the algorithm
#!# every 100 iterations indicating how much time takes. This helps to 
#!# adjust the expectations on how long it will take to run the model.

rdm_outfile <- "Mplus/data2rdm_2rdm.out"

if (!file.exists(rdm_outfile)) {
  MplusAutomation::runModels("Mplus/data2rdm_2rdm.inp")
}

# To compute the reliabilities in the 2RDM, the parameter estimates and factor 
# scores from the MCMC estimation are needed. We wrote a custom function that
# automatically reads into R the Mplus output of the model and computes all
# the coefficients of interest. The arguments of this function are the path to 
# the output file and the names of the items. 
rel_rdm_data2rdm <- reliability_2rdm(rdm_outfile, items = items)

# Between-person reliability
rel_rdm_data2rdm$rdm_between_reliability

# Within-person reliabilities
summary(rel_rdm_data2rdm$rdm_within_reliability$within_rel)
hist(rel_rdm_data2rdm$rdm_within_reliability$within_rel, xlim = c(0, 1))

# ME-TSO ----
# The ME-TSO is fitted in Mplus. This can take long depending on the amount of 
# data, the number of iterations, and other factors. This analysis took 
# approximately 2 minutes.

#!# Here, we check if an output file from running the model already exists.
#!# If it does not, the model if run in Mplus from within R. However,
#!# it might be better to simply run the model directly in Mplus, because 
#!# when running in Mplus, the program shows the progress of the algorithm
#!# every 100 iterations indicating how much time takes. This helps to 
#!# adjust the expectations on how long it will take to run the model.

metso_outfile <- "Mplus/data2rdm_metso.out"

if (!file.exists(metso_outfile)) {
  MplusAutomation::runModels("Mplus/data2rdm_metso.inp")
}

# To compute the reliabilities in the ME-TSO, the parameter estimates and factor 
# scores from the MCMC estimation are needed. We wrote a custom function that
# automatically reads into R the Mplus output of the model and computes all
# the coefficients of interest. The arguments of this function are the path to 
# the output file and the names of the items.
rel_metso_data2rdm <- reliability_metso(metso_outfile, items = items)

# Median reliability across persons per item
apply(rel_metso_data2rdm$metso_reliability[, items], 2, median)

# Compare person-specific reliabilities between PFA, DFA, and 2RDM ----

# Merge estimated coefficients into one dataset. 
# It is important to merge based on the participants' ID because Mplus output
# reorders the participants based on the number of observations (compliance).
rel_est <- merge(pfa_data2rdm$pfa_omega, dfa_data2rdm$dfa_omega, by = "id")
rel_est <- merge(rel_est, rel_rdm_data2rdm$rdm_within_reliability, by = "id")
rel_est <- merge(rel_est, rel_metso_data2rdm$metso_reliability, by = "id")

# Compute correlations
cor(rel_est[, c(2, 4, 6)],
    use = "pairwise")

# Histogram PFA vs. DFA
plot(rel_est$pfa_omega, 
     rel_est$dfa_omega, 
     xlim = c(0, 1),
     ylim = c(0, 1))
abline(a = 0, b = 1, col = "red")

# Histogram PFA vs. 2RDM
plot(rel_est$pfa_omega, 
     rel_est$within_rel, 
     xlim = c(0, 1),
     ylim = c(0, 1))
abline(a = 0, b = 1, col = "red")

# Histogram DFA vs. 2RDM
plot(rel_est$dfa_omega,
     rel_est$within_rel,
     xlim = c(0, 1),
     ylim = c(0, 1))
abline(a = 0, b = 1, col = "red")

#### END #### 