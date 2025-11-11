# Computing the internal consistency reliability

# This file shows how to estimate the internal consistency reliability of 
# the dataset "datelisa" based on six/seven different approaches.
# The data comes from dr. Leonie Vogelsmeier and collaborators. These data
# were collected during a ESM study that lasted 14 days. Participants received
# a notification in two-hours block to fill in the questionnaire up to six times
# a day. The questionnaire included ten items that measure momentary life 
# satisfaction. However, not all items were administered every time. On each
# occasion, participants report their responses on seven items. One item
# was an anchor item and always present. The other six were selected randomly
# out of the remaining nine.

# Set up ----

# Install required packages
packages <- c("remotes", "psych", "misty", "lavaan", "dynr", "MplusAutomation")
install.packages(setdiff(packages, rownames(installed.packages())))
rm(packages)

# Install esmpack package
if (!"esmpack" %in% rownames(installed.packages())) {
  remotes::install_github("wviechtb/esmpack")
}

# Load custom functions
custom_functions <- list.files("R/", full.names = TRUE)
sapply(custom_functions, source)
rm(custom_functions)

# Data loading and cleaning ----

# Load data
datelisa <- readRDS("Data/ELISA_CR.rds")
datelisa <- as.data.frame(datelisa)

# Store items' names.
items <- c('consider_happy', 'satisfied_life', 'change_many',
           'life_ideal', 'makes_happy', 'satisfied_situation',
           'change_nothing', 'leaves_desired', 'important_things', 'line_life')

# Keep variables of interest
# id: participants ID
# n_obs: number of the valid observation. this ignores complete missingness and nights
# Items: 'consider_happy', 'satisfied_life', 'change_many', 'life_ideal', 
#        'makes_happy', 'satisfied_situation', 'change_nothing', 'leaves_desired', 
#        'important_things', 'line_life'

# Select and reorder the variables of interest for these analyses.
datelisa <- datelisa[, c("id", "n_obs", items)]

# Rename items for future use in Mplus.
# Variables names in Mplus cannot be longer than 8 characters. They should be
# preferably shorter. To make new names, we take the first three letters of 
# each word in the item's name and paste it together. 
# For example,  "line_life" turns into "linlif".
items <- sapply(strsplit(items, split = "_"), function(x) {
  paste0(strtrim(x, 3), collapse = "")
})

names(datelisa) <- c("id", "n_obs", items)

# Reorder items so the first item is line_life, which is the anchor item.
items <- items[c(10, 1:9)]
datelisa <- datelisa[, c("id", "n_obs", items)]

# Get the number of items. Needed for later.
I <- 10

# Reshape data to have one row for the response of a person to a given item 
# at a certain timepoint.
datelisa_long <- reshape(
  datelisa,
  direction = "long",
  idvar = c("id", "n_obs"),
  timevar = "item",
  times = c(items),
  varying = list(c(items)),
  v.names = "resp"
)

rownames(datelisa_long) <- NULL

datelisa_long <- na.omit(datelisa_long)

# Compute compliance per person
compliance <- esmpack::calc.nomiss(id, id, datelisa)

# Create subset with participants who completed at least 50 observations. 
# These is needed for the approaches that are more idiosyncratic. 
# Get IDs of the participants to exclude
id_exclude <- names(which(compliance < 50))

datelisa_50 <- datelisa[!(datelisa$id %in% id_exclude), ]

# GT and MLM ----

# The reliability based on the GT and MLM can be easily computed with the 
# function mlr or multilevel.reliability of the psych package.

# We do it twice. Once with the data in its original structure, meaning that
# one row has the response of a person to all the items at a certain timepoint,
# and another time with the restructured data with one row for the response 
# of a person to a given item at a certain timepoint. 
t0 <- proc.time()
mlr_datelisa <- psych::mlr(
  datelisa, 
  grp   = "id",
  Time  = "n_obs",
  items = items,
  lmer  = TRUE,
  lme   = FALSE,
  alpha = FALSE,
  aov   = FALSE
)
t1 <- proc.time() - t0

t1[3]

t0 <- proc.time()
mlr_datelisa2 <- psych::mlr(
  datelisa_long, 
  grp   = "id",
  Time  = "n_obs",
  items = "item",
  values = "resp",
  long = TRUE,
  lmer  = TRUE,
  lme   = FALSE,
  alpha = FALSE,
  aov   = FALSE
)
t1 <- proc.time() - t0

t1[3]

# GT Reliability estimates
mlr_datelisa$RkR # between-person
mlr_datelisa$Rc  # within-person
mlr_datelisa2$RkR # between-person
mlr_datelisa2$Rc  # within-person

# MLM Reliability estimates
mlr_datelisa$RkRn # between-person
mlr_datelisa$Rcn  # within-person
mlr_datelisa2$RkRn # between-person
mlr_datelisa2$Rcn  # within-person

# Multilevel omega ----

# To estimate the reliability based on the multilevel CFA, one can easily use 
# the function multilevel.omega from the misty package. 

mlcfa_datelisa <- misty::multilevel.omega(
  datelisa[, items],
  cluster = datelisa$id,
  missing = "fiml"
)

# Between and within omega
mlcfa_datelisa$result$omega

# P-technique factor analysis ----

# The P-technique implies fitting a factor model to the responses of each
# person. To do this, we wrote a custom function that automatically loops
# through each person to fit the model in lavaan. 
# The arguments of this function are the data, the label of the "id" variable,
# the names of the items, and the how the model is identified. See source for
# further details.

#!# Make sure your data is a data.frame!

pfa_datelisa <- esm_pfa_omega(data  = datelisa_50, 
                              idlab = "id", 
                              items = items, 
                              constraint = "variance",
                              missing = "fiml",
                              em.h1.iter.max = 1000)

# Summary and histogram of estimated reliabilities
# The unidimensional P-technique factor model only converged to a proper solution
# for 7 out of 54 participants.
summary(pfa_datelisa$pfa_omega$pfa_omega)
hist(pfa_datelisa$pfa_omega$pfa_omega, xlim = c(0, 1))

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

dfa_datelisa <- esm_dfa_omega(data    = datelisa_50,
                              idlab   = "id",
                              timelab = "n_obs",
                              items   = items,
                              missing = "fiml")

# Summary and histogram of estimated reliabilities
summary(dfa_datelisa$dfa_omega$dfa_omega)
hist(dfa_datelisa$dfa_omega$dfa_omega, xlim = c(0, 1))

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

if (!file.exists("Mplus/datelisa.dat")) {
  MplusAutomation::prepareMplusData(
    datelisa_50,
    filename = "Mplus/datelisa.dat",
    inpfile = FALSE
  )
}

# 2RDM ----
# The 2RDM is fitted in Mplus. This can take long depending on the amount of 
# data, the number of iterations, and other factors. This analysis took 
# approximately 30 minutes. 

#!# Here, we check if an output file from running the model already exists.
#!# If it does not, the model if run in Mplus from within R. However,
#!# it might be better to simply run the model directly in Mplus, because 
#!# when running in Mplus, the program shows the progress of the algorithm
#!# every 100 iterations indicating how much time takes. This helps to 
#!# adjust the expectations on how long it will take to run the model.

rdm_outfile <- "Mplus/datelisa_2rdm.out"

if (!file.exists(rdm_outfile)) {
  MplusAutomation::runModels("Mplus/datelisa_2rdm.inp")
}

# To compute the reliabilities in the 2RDM, the parameter estimates and factor 
# scores from the MCMC estimation are needed. We wrote a custom function that
# automatically reads into R the Mplus output of the model and computes all
# the coefficients of interest. The arguments of this function are the path to 
# the output file and the names of the items.
rel_rdm_datelisa <- reliability_2rdm(rdm_outfile, items = items)

# Between-person reliability
rel_rdm_datelisa$rdm_between_reliability

# Within-person reliabilities
summary(rel_rdm_datelisa$rdm_within_reliability$within_rel)
hist(rel_rdm_datelisa$rdm_within_reliability$within_rel, xlim = c(0, 1))

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

metso_outfile <- "Mplus/datelisa_metso.out"

if (!file.exists(metso_outfile)) {
  MplusAutomation::runModels("Mplus/datelisa_metso.inp")
}

# To compute the reliabilities in the ME-TSO, the parameter estimates and factor 
# scores from the MCMC estimation are needed. We wrote a custom function that
# automatically reads into R the Mplus output of the model and computes all
# the coefficients of interest. The arguments of this function are the path to 
# the output file and the names of the items.
rel_metso_datelisa <- reliability_metso(metso_outfile, items = items)

# Median reliability across persons per item
apply(rel_metso_datelisa$metso_reliability[, items], 2, median)

# Compare person-specific reliabilities between PFA, DFA, and 2RDM ----

# Merge estimated coefficients into one dataset. 
# It is important to merge based on the participants' ID because Mplus output
# reorders the participants based on the number of observations (compliance).
rel_est <- merge(pfa_datelisa$pfa_omega, dfa_datelisa$dfa_omega, by = "id")
rel_est <- merge(rel_est, rel_rdm_datelisa$rdm_within_reliability, by = "id")
rel_est <- merge(rel_est, rel_metso_datelisa$metso_reliability, by = "id")

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