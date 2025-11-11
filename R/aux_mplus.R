# Auxiliary functions to read Mplus output and compute the reliability of the 
# ME-TSO and 2RDM

# readFScores ----
# Function to read the factor scores of Mplus DSEM bayesian models.

# This function is needed because the function readModel from the package
# MplusAutomation sometimes fails to read the factor scores correctly.

# outfile: output file of Mplus
# clus: name of the id or cluster variable in uppercase. Optional.
# ranef.vars: name of the random effects of interest. Optional.

readFScores <- function(outfile, clus = NULL, ranef.vars = NULL) {
  # Import output file
  mplus_out <- readLines(outfile)
  
  # Extract total number of observations
  nobs <- grep("Number of observations", 
               mplus_out, 
               value = TRUE)
  nobs <- sub("Number of observations\\s+(\\d+)", "\\1", nobs)
  nobs <- as.numeric(nobs)
  # Extract number of participants
  nid <- grep("Number of clusters", 
              mplus_out, 
              value = TRUE)
  nid <- sub("Number of clusters\\s+(\\d+)", "\\1", nid)
  nid <- as.numeric(nid)
  
  # Import factor scores file
  # Get FScores file name from outfile
  fsfile <- grep("FILE is .*\\.dat", mplus_out, value = TRUE)
  fsfile <- gsub("^\\s+FILE is\\s", "", fsfile)
  fsfile <- gsub("\\s*;$", "", fsfile)
  # FScores file path in wd.
  fsfile <- grep(fsfile, list.files(recursive = TRUE), value = TRUE)
  
  # Read FScores file into R
  fscores_raw <- scan(fsfile, what = "c", na.strings = "*")
  fscores_raw <- as.numeric(fscores_raw)
  fscores_raw <- matrix(fscores_raw,
                        nrow  = nobs,
                        ncol  = length(fscores_raw)/nobs,
                        byrow = TRUE)
  
  # Get fscores variable names
  
  # Get lines where the names of the factor scores are
  line_start <- grep("Order and format of variables", mplus_out)
  line_end   <- grep(
    "\\+ variables that have a value for each of the \\d+ imputations",
    mplus_out)
  # Extract number of iterations used to compute the factor scores
  n_fs_samples <- sub(
    ".*\\+ variables that have a value for each of the (\\d+) imputations.*", 
    "\\1", mplus_out[line_end])
  n_fs_samples <- as.numeric(n_fs_samples)
  # Obtain the name of the factor scores
  fscoresnames <- mplus_out[(line_start + 2):(line_end - 2)]
  # Identify names that need to be repeated for the number of iterations
  which_samples <- grep("^\\s+\\+", fscoresnames)
  # Clean the factor scores names by removing unnecessary spaces
  fscoresnames <- gsub("^\\s+|^\\s+\\+|^\\s+_", "", fscoresnames)
  fscoresnames <- gsub("\\s{2,}..+", "", fscoresnames)
  
  fscoresnames <- 
    lapply(fscoresnames, function(x) {
      if (x %in% fscoresnames[which_samples]) {
        return(paste0(x, sprintf("_i%03d", 1:n_fs_samples)))
      } else {
        return(x)
      }
    })
  
  fscoresnames <- unlist(fscoresnames)
  
  # Rename variables in the factor scores
  colnames(fscores_raw) <- fscoresnames
  
  # Reduce fscores data.
  if (!is.null(clus)) {
    # Select one row per participant/cluster.
    # This is useful when one wants to extract the estimated values of random
    # effects. As the factor scores of these variables are equal for each 
    # subject.
    # Get the first row of each participant in the data
    firstrow_id <- match(unique(fscores_raw[, clus]), fscores_raw[, clus])
    
    # Save only one row per ID
    fscores <- fscores_raw[firstrow_id, ]
    # Save ids/clusters to name the rows of output matrix when
    # ranef.vars are provided.
    id_names <- fscores[, clus]
    
    if (!is.null(ranef.vars)) {
      tmp_vars <- grep(paste0("^",  ranef.vars, ".*_i", collapse = "|"), 
                       fscoresnames, value = TRUE)
      fscores <- fscores[, tmp_vars]
      
      # Get unique variables
      tmp_vars <- unique(sub("_i.*$", "", tmp_vars))
      
      if (length(tmp_vars) == 1L) {
        row.names(fscores) <- id_names
      } else {
        fscores_array <- array(
          fscores,
          dim = c(nid, n_fs_samples, length(tmp_vars)),
          dimnames = list(ids     = id_names,
                          iter    = paste0("it", 1:n_fs_samples),
                          varname = tmp_vars)
        )
        
        fscores <- fscores_array
      }
    }
  } else {
    fscores <- fscores_raw
  }
  
  return(fscores)
}

# metso.var.coeff ----

# Function to compute the variance coefficients  of the me-tso. 
# Retrieved from "https://github.com/secastroal/ME-TSO"

metso.var.coeff <- function(within.parameters, between.parameters, ind.ar, id){
  
  bp <- between.parameters
  wp <- within.parameters
  
  # Coefficients across fixed situations
  # Fixed situation variances
  fxi.var <- bp$xi.var + 2 * bp$xi.var * bp$beta1 + bp$xi.var * (bp$beta1 ^ 2) + bp$omega.var 
  # Reference X fixed situations covariances
  fxi.cov <- bp$xi.var + bp$xi.var * bp$beta1
  # Reference X fixed situations correlations
  fxi.cor <- fxi.cov / sqrt(bp$xi.var * fxi.var)
  
  # Consistency of traits
  con.trait <- data.frame(round((fxi.cor) ^ 2, 2))
  # Situation specificity of traits
  spe.trait <- data.frame(round(1 - con.trait, 2))
  
  names(con.trait) <- names(spe.trait) <- paste0("Xi_j", 1:(dim(con.trait)[2]))
  
  # Total variance of Gamma
  gamma.var <- bp$xi.var * (bp$beta1 ^ 2) + bp$omega.var
  # Person-situation interaction coefficient
  int.coeff  <- data.frame(round((bp$xi.var * (bp$beta1 ^ 2)) / gamma.var, 2))
  unique.eff <- data.frame(round(bp$omega.var / gamma.var, 2))
  
  names(int.coeff) <- names(unique.eff) <- names(con.trait)
  
  coeff.fixed <- list(trait.consistency = con.trait,
                      trait.specificity = spe.trait,
                      interaction.coeff = int.coeff,
                      unique.effect     = unique.eff)
  
  # Coefficients across random situations within fixed situations
  #create list to store all the coefficients across random situations
  coeff.random <- list()
  
  # Reference situation
  # Total variance
  trait.var       <- matrix(bp$xi.var, length(id), length(bp$xi.var), byrow = TRUE)
  state.loadings  <- matrix(wp$lambda.state, length(id), length(bp$xi.var), byrow = TRUE)
  state.resid     <- matrix(wp$epsilon.var, length(id), length(bp$xi.var), byrow = TRUE)
  autoreg.effects <- as.vector(ind.ar ^ 2 / (1 - ind.ar ^ 2))
  
  tot.var <- trait.var + autoreg.effects * wp$zeta.var * state.loadings ^ 2 +
    wp$zeta.var * state.loadings ^ 2 + state.resid
  
  rel   <- data.frame(round((trait.var + autoreg.effects * wp$zeta.var * state.loadings ^ 2 +
                               wp$zeta.var * state.loadings ^ 2) / tot.var, 2)) # reliability
  con   <- data.frame(round((trait.var + autoreg.effects * wp$zeta.var * state.loadings ^ 2) / tot.var, 2)) # consistency
  pred  <- data.frame(round((trait.var) / tot.var, 2)) # predictability by trait
  upred <- data.frame(round((autoreg.effects * wp$zeta.var * state.loadings ^ 2) / tot.var, 2)) # unpredictability by trait
  spe   <- data.frame(round((wp$zeta.var * state.loadings ^ 2) / tot.var, 2)) # occasion specificity
  
  names(rel) <- names(con) <- names(pred) <-
    names(upred) <- names(spe) <- paste0("Y_", 1:length(bp$xi.var), "0")
  row.names(rel) <- row.names(con) <- row.names(pred) <-
    row.names(upred) <- row.names(spe) <- id
  
  coeff.random[[1]] <- list(reliability = rel,
                            consistency = con,
                            specificity = spe,
                            predictability   = pred,
                            unpredictability = upred)
  
  for (f in 1:(dim(fxi.var)[2])) {
    trait.var       <- matrix(fxi.var[, f], length(id), length(bp$xi.var), byrow = TRUE)
    state.loadings  <- matrix(wp$lambda.state, length(id), length(bp$xi.var), byrow = TRUE)
    state.resid     <- matrix(wp$epsilon.var, length(id), length(bp$xi.var), byrow = TRUE)
    autoreg.effects <- as.vector(ind.ar ^ 2 / (1 - ind.ar ^ 2))
    
    tot.var <- trait.var + autoreg.effects * wp$zeta.var * state.loadings ^ 2 +
      wp$zeta.var * state.loadings ^ 2 + state.resid
    
    rel   <- data.frame(round((trait.var + autoreg.effects * wp$zeta.var * state.loadings ^ 2 +
                                 wp$zeta.var * state.loadings ^ 2) / tot.var, 2)) # reliability
    con   <- data.frame(round((trait.var + autoreg.effects * wp$zeta.var * state.loadings ^ 2) / tot.var, 2)) # consistency
    pred  <- data.frame(round((trait.var) / tot.var, 2)) # predictability by trait
    upred <- data.frame(round((autoreg.effects * wp$zeta.var * state.loadings ^ 2) / tot.var, 2)) # unpredictability by trait
    spe   <- data.frame(round((wp$zeta.var * state.loadings ^ 2) / tot.var, 2)) # occasion specificity
    
    names(rel) <- names(con) <- names(pred) <-
      names(upred) <- names(spe) <- paste0("Y_", 1:length(bp$xi.var), f)
    row.names(rel) <- row.names(con) <- row.names(pred) <-
      row.names(upred) <- row.names(spe) <- id
    
    coeff.random[[f + 1]] <- list(reliability = rel,
                                  consistency = con,
                                  specificity = spe,
                                  predictability   = pred,
                                  unpredictability = upred)
  }
  
  names(coeff.random) <- paste0("Xi_j", 0:(dim(fxi.var)[2]))
  
  out <- list(fixed.situations  = coeff.fixed,
              random.situations = coeff.random)
  
  return(out)
}

# extract.metso.bpar ----

# Function to extract estimated parameters of the ME-TSO for each iteration of the 
# Bayesian samples in Mplus.
# Retrieved from https://osf.io/k56vp/

extract.metso.bpar <- function(bpar, I) {
  # Turn bayesian samples of the parameters for a given iteration
  # into a vector.
  bpar <- bpar[, grep("Parameter", names(bpar))]
  bpar <- t(bpar)
  
  # Set number of situation covariates to 0.
  D <- 0
  
  # Number of elements in the between-level var-cov matrix
  NC <- (I * (I + 1)) / 2
  
  # Retrieve within and between parameters from bpar
  within_par  <- list(lambda.state = c(1, bpar[1:(I - 1)]), # Factor Loadings
                      epsilon.var  = bpar[(I):(2 * I - 1)], # Measurement Error 
                      zeta.var     = bpar[2 * I]) # Latent State Variance
  between_par <- list(beta1 = matrix(0, I, 1), # Person-Situation Effect
                      beta0 = matrix(0, I, 1), # Person-Situation Effect 
                      xi.means = bpar[(2 * I + 1):(3 * I)], # Latent Traits Mean
                      ar.mean  = bpar[3 * I + 1], # Mean Autoregressive Effect
                      xi.var   = bpar[cumsum(1:I) + 3 * I + 1], # Latent Traits Var
                      ar.var   = bpar[4 * I + NC + 2], # AR Effect Variance
                      omega.var = 0) # Person-Situation Residual
  
  out <- list(within = within_par, between = between_par)
  
  return(out)
}

# getFScores ----

# Function to extract the factor scores from an Mplus outfile

# outfile: output file of Mplus
# clus: name of the id or cluster variable in uppercase. Optional.
# ranef.vars: name of the random effects of interest. Optional.

getFScores <- function(outfile, clus = NULL, ranef.vars = NULL) {
  # Read Mplus output
  mplus_out <- MplusAutomation::readModels(outfile)
  
  # Extract factor scores imputations of the fixed effects
  fscores_fixed <- mplus_out$bparameters$factor.score.imputations
  
  # Read factor scores of random effects
  
  # If factor scores were read correctly use MplusAutomation else use readFScores
  if (length(mplus_out$savedata) > 0L) {
    fscores_raw <- mplus_out$savedata
    
    # Reduce fscores data
    if (!is.null(clus)) {
      # Select one row per participant/cluster.
      # This is useful when one wants to extract the estimated values of random
      # effects. As the factor scores of these variables are equal for each 
      # subject.
      # Get the first row of each participant in the data
      firstrow_id <- match(unique(fscores_raw[, clus]), fscores_raw[, clus])
      
      # Save only one row per ID
      fscores_rand <- fscores_raw[firstrow_id, ]
      # Save ids/clusters to name the rows of output matrix when
      # ranef.vars are provided.
      id_names <- fscores_rand[, clus]
      
      if (!is.null(ranef.vars)) {
        tmp_vars <- grep(paste0("^",  ranef.vars, ".*_", collapse = "|"), 
                         names(fscores_raw), value = TRUE)
        fscores_rand <- fscores_rand[, tmp_vars]
        fscores_rand <- as.matrix(fscores_rand)
        
        # Get unique variables
        tmp_vars <- unique(sub(".I_.*$", "", tmp_vars))
        
        if (length(tmp_vars) == 1L) {
          row.names(fscores_rand) <- id_names
        } else {
          nid <- length(id_names)
          n_fs_samples <- ncol(fscores_rand)/length(tmp_vars)
            
          fscores_array <- array(
            fscores_rand,
            dim = c(nid, n_fs_samples, length(tmp_vars)),
            dimnames = list(ids     = id_names,
                            iter    = paste0("it", 1:n_fs_samples),
                            varname = tmp_vars)
          )
          
          fscores_rand <- fscores_array
        }
      }
    } else {
      fscores_rand <- fscores_raw
    }
  } else {
    fscores_rand <- readFScores(outfile, clus = clus, ranef.vars = ranef.vars)
  }
  
  # Return factor scores of fixed and random effects
  
  out <- list(
    fscores_fixed = fscores_fixed,
    fscores_rand  = fscores_rand
  )
}

# reliability_metso ----

# Function to compute the reliability based on the ME-TSO given the Mplus 
# output file

# outfile: output file of Mplus
# items: character with the name of the item variables 
# clus: name of the id or cluster variable in uppercase. Default "ID".
# ranef.vars: label of the random autoregressive effect. Default "AR".

reliability_metso <- function(outfile, items, clus = "ID", ranef.vars = "AR") {
  I <- length(items)
  
  # Get factor scores of the fixed and random effects
  fscores <- getFScores(outfile, clus = clus, ranef.vars = ranef.vars)
  # fixed factor scores
  fscores_fixed <- fscores$fscores_fixed
  fscores_rand  <- fscores$fscores_rand
  
  # Extract MCMC sample of parameters of interest
  metso_var_samples <- list()
  
  for (j in 1:nrow(fscores_fixed)) {
    # Extract the j-th samples of the parameter estimates
    bpar_tmp <- extract.metso.bpar(fscores_fixed[j, ], I)
    # Extract the j-th samples of the individual AR effects
    ar_tmp <- fscores_rand[, j]
    # Compute ME-TSO variance coefficients for the j-th samples
    coeff_tmp <- metso.var.coeff(
      within.parameters  = bpar_tmp$within,
      between.parameters = bpar_tmp$between,
      ind.ar = ar_tmp,
      id     = names(ar_tmp)
    )
    metso_var_samples[[j]] <- coeff_tmp
  }
  rm(bpar_tmp, ar_tmp, coeff_tmp, j)
  
  metso_rel_samples <- list()
  
  # Extract reliability samples
  metso_rel_samples <- lapply(
    metso_var_samples, function(x) {
      x$random.situations$Xi_j0$reliability
    }
  )
  # Restructure reliability samples to array
  metso_rel_samples <- array(
    data = unlist(metso_rel_samples),
    dim = c(nrow(metso_rel_samples[[1]]),
            ncol(metso_rel_samples[[1]]),
            length(metso_rel_samples)),
    dimnames = list(id = rownames(metso_rel_samples[[1]]),
                    items = items,
                    iter  = paste0("it_", 1:length(metso_rel_samples))
    )
  )
  
  # Compute reliability estimates, the median of the samples
  
  tmp <- apply(metso_rel_samples, c(1, 2), median, na.rm = TRUE)
  id <- dimnames(metso_rel_samples)[[1]]
  
  metso_rel_sum <- data.frame(id, tmp) 
  
  out <- list(
    metso_reliability = metso_rel_sum,
    metso_rel_samples = metso_rel_samples
  )
  
  return(out)
}

# reliability_2rdm ----

# Function to compute the reliability based on the unidimensional 2RDM given 
# the Mplus output file

# outfile: output file of Mplus
# items: character with the name of the item variables 
# clus: name of the id or cluster variable in uppercase. Default "ID".
# ranef.vars: labels of the random effects in the 2RDM. Not all labels are 
#             needed, only what is common to each type of random effect. 
#             It consist of a character vector with four elements: The label of 
#             the random factor loadings, the label of the autoregressive effect,
#             the label of the random measurement error variances, and the label 
#             of the random innovation variance. Default c("LAM","AR","WERR","FWVAR").
# brel: label of the between-person reliability. Default "BREL".

reliability_2rdm <- function(outfile, 
                             items, 
                             clus       = "ID", 
                             ranef.vars = c("LAM","AR","WERR","FWVAR"),
                             brel       = "BREL") {
  I <- length(items)
  
  # Get factor scores of the fixed and random effects
  fscores <- getFScores(outfile, clus = clus, ranef.vars = ranef.vars)
  # Fixed effects factor scores
  fscores_fixed <- fscores$fscores_fixed
  # Random effects factor scores
  fscores_rand  <- fscores$fscores_rand
  
  # Get random variable labels
  load_lab  <- dimnames(fscores_rand)[[3]][1:(I - 1)]
  ar_lab    <- dimnames(fscores_rand)[[3]][I]
  werr_lab  <- dimnames(fscores_rand)[[3]][(I + 1):(2 * I)]
  fwvar_lab <- dimnames(fscores_rand)[[3]][2 * I + 1]
  
  # Transform variances in log scale back to normal
  fscores_rand[, , c(werr_lab, fwvar_lab)] <- 
    exp(fscores_rand[, , c(werr_lab, fwvar_lab)])
  
  # Get number of participants and number of MCMC iterations
  nid <- dim(fscores_rand)[1]
  n_fs_samples <- dim(fscores_rand)[2]
  
  # Within-person reliability
  rdm_wrel_samples <- matrix(NA, nid, n_fs_samples)
  for (i in 1:n_fs_samples) {
    tmp <- fscores_rand[, i, ]
    # Observed Variance
    ovar <- tmp[, fwvar_lab] / (1 - tmp[, ar_lab]^2)
    # Squared sum of loadings
    sqload <- rowSums(cbind(1, tmp[, load_lab]))^2
    # Sum of measurement error
    ersum <- rowSums(tmp[, werr_lab])
    # Within-person reliability 
    rdm_wrel_samples[, i] <- (sqload * ovar) / (sqload * ovar + ersum)
  }
  rm(i, tmp)
  
  # Turn wrel samples into data.frame
  rdm_wrel_samples <- as.data.frame(rdm_wrel_samples)
  names(rdm_wrel_samples) <- paste0("it_", 1:n_fs_samples)
  row.names(rdm_wrel_samples) <- dimnames(fscores_rand)[[1]]
  # Remove within-person reliability estimates lower than 0
  rdm_wrel_samples[rdm_wrel_samples < 0] <- NA 
  
  # Compute the median across iterations 
  within_rel <- apply(rdm_wrel_samples, 1, median, na.rm = TRUE)
  id         <- dimnames(fscores_rand)[[1]]
  rdm_wrel_sum <- data.frame(id, within_rel)
  
  # Between-person reliability
  
  # Get label of the between-person reliability parameter
  brel_lab <- grep(brel, names(fscores_fixed), value = TRUE)
  
  rdm_brel_samples <- fscores_fixed[, brel_lab]
  
  rdm_brel_sum <- c(
    estimate    = median(rdm_brel_samples, na.rm = TRUE),
    posteriorSD = sd(rdm_brel_samples, na.rm = TRUE),
    lowerCI     = quantile(rdm_brel_samples, na.rm = TRUE, probs = 0.025),
    upperCI     = quantile(rdm_brel_samples, na.rm = TRUE, probs = 0.975)
  )
  
  rdm_brel_sum <- round(rdm_brel_sum, 3)
  
  out <- list(
    rdm_within_reliability = rdm_wrel_sum,
    rdm_between_reliability = rdm_brel_sum,
    rdm_wrel_samples = rdm_wrel_samples,
    rdm_brel_samples = rdm_brel_samples
  )
  
  return(out)
}




