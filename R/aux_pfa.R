# Auxiliary functions to compute the reliability based on the P-Technique
# and the Dynamic Factor Analysis.

# pfa_syntax ----
# Function to write unidimensional cfa lavaan syntax

# items: character with the name of the item variables
# constraint: either "variance", "loading" or "tau", defines how the cfa model is 
#             identified. If "variance", the variance of the latent factor is 
#             fixed at 1. If "loading", the factor loading of the first
#             item is fixed at 1. If "tau", all factor loadings are fixed at 1,
#             which is useful when there are 3 or less items. Default is "variance"

pfa_syntax <- function(items, constraint = "variance") {
  # Get number of items
  I <- length(items)
  
  # Create loadings labels
  
  load_lab <- paste0("lw", 1:I)
  
  if (constraint == "loading") {load_lab[1] <- 1}
  if (constraint == "tau") {load_lab <- rep(1, I)}
  
  # Create factor variance label
  
  fvar_lab <- "fvar"
  
  if (constraint == "variance") {fvar_lab <- 1}
  
  # Create measurement error labels
  
  meas_lab <- paste0("ew", 1:I)
  
  # Create factor model syntax 
  
  load_syntax <- paste0(load_lab, "*", items, collapse = " + ")
  
  if (constraint == "variance") {
    load_syntax <- paste(paste0("NA*", items[1]), "+", load_syntax)
  }
  
  factor_syntax <- paste("wf =~", load_syntax) 
  
  # Variances and measurement error syntaxes
  
  fvar_syntax <- paste("wf", "~~", paste0(fvar_lab, "*wf"))
  meas_syntax <- paste(items, "~~", meas_lab, "*", items, collapse = "\n")
  
  # Person specific omega syntax
  
  # Squared of Sum of loadings 
  sqload <- paste0("(", paste(load_lab, collapse = " + "), ")^2")
  # Sum of measurement error variances
  ersum <- paste(meas_lab, collapse = " + ")
  # Variance of the true score
  truevar <- paste(sqload, "*", fvar_lab)
  # Total variance
  totalvar <- paste(truevar, "+", "(", ersum, ")")
  # Omega
  omg_syntax <- paste0("omega := (", truevar, ") / (", totalvar, ")")
  
  # Model syntax
  model_syntax <- paste(
    "# P-Technique unidimensional factor model",
    factor_syntax,
    fvar_syntax,
    meas_syntax,
    "\n# Person specific omega",
    omg_syntax,
    sep = "\n"
  )
  
  return(model_syntax)
}

# esm2list ----

# Function to turn ESM data into a list. Where each element in the list is the 
# data from a specific participant.

# data: Data frame in long format including an id variable and a time variable.
# id: vector of the ID variable
# include: character or numeric vector indicating which variables to include.
# exclude: character or numeric vector indicating which variables to exclude.

esm2list <- function(data, id, include = NULL, exclude = NULL) {
  if (is.null(data)) 
    stop("Argument 'data' must be specified.")
  
  if (!is.data.frame(data)) {
    data <- data.frame(data)
  }
  
  mf <- match.call()
  mf.id <- mf[[match("id", names(mf))]]
  id <- eval(mf.id, data, enclos = sys.frame(sys.parent()))
  if (is.null(id)) 
    stop("Argument 'id' must be specified.")
  if (any(is.na(id))) 
    stop("Argument 'id' should not contain any NAs.")
  
  varnames <- names(data)
  nvars    <- length(varnames)
  
  if (!is.null(include) & !is.null(exclude)) {
    exclude <- NULL
    if (is.numeric(include)) {
      include <- varnames[include]
    }
    warning("Both arguments 'include' and 'exclude' were specified. ",
            "The following variables were included in the individual datasets: ",
            paste0(include, c(rep(", ", length(include) - 1), "."), collapse = ""))
  }
  
  if (is.null(include) & is.null(exclude)) {
    vars <- 1:length(varnames)
  }
  
  if (!is.null(include)) {
    if (!(is.character(include) | is.numeric(include))) 
      stop("Argument 'include' must either be a character or a numeric vector.")
    if (is.character(include)) {
      vars.pos <- lapply(include, function(x) {
        pos <- charmatch(x, varnames)
        if (is.na(pos)) 
          stop("Variable '", x, "' not found in the data frame.", 
               call. = FALSE)
        if (pos == 0L) 
          stop("Multiple matches for variable '", x, 
               "' in the data frame.", call. = FALSE)
        return(pos)
      })
      vars <- unique(unlist(vars.pos))
    } else {
      vars.pos <- unique(round(include))
      if (min(vars.pos) < 1 | max(vars.pos) > nvars) { 
        stop("Variables positions must be between 1 and ", 
             nvars, ".")
      }
      vars <- vars.pos
    }
  }
  
  if (!is.null(exclude)) {
    if (!(is.character(exclude) | is.numeric(exclude))) 
      stop("Argument 'exclude' must either be a character or a numeric vector.")
    if (is.character(exclude)) {
      vars.pos <- lapply(exclude, function(x) {
        pos <- charmatch(x, varnames)
        if (is.na(pos)) 
          stop("Variable '", x, "' not found in the data frame.", 
               call. = FALSE)
        if (pos == 0L) 
          stop("Multiple matches for variable '", x, 
               "' in the data frame.", call. = FALSE)
        return(pos)
      })
      vars <- (1:length(varnames))[-unique(unlist(vars.pos))]
    } else {
      vars.pos <- unique(round(exclude))
      if (min(vars.pos) < 1 | max(vars.pos) > nvars) { 
        stop("Variables positions must be between 1 and ", 
             nvars, ".")
      }
      vars <- (1:length(varnames))[-vars.pos]
    }
  }
  
  temp <- list()
  ids <- sort(unique(id))
  
  for (i in 1:length(ids)) {
    temp[[i]] <- data[id == ids[i], vars]
  }
  rm(i)
  
  if (!is.character(ids)) {
    ids <- as.character(ids)
  }
  
  names(temp) <- ids
  
  return(temp)
}

# esm_pfa_omega ----

# Function to compute the omega reliability based on the P-technique factor 
# analysis.

# data: data frame with the ESM data including the id variable and the time variable
# idlab: character indicating the name of the id variable in the data.
# items: character with the name of the item variables
# constraint: either "variance", "loading" or "tau", defines how the cfa model is 
#             identified. If "variance", the variance of the latent factor is 
#             fixed at 1. If "loading", the factor loading of the first
#             item is fixed at 1. If "tau", all factor loadings are fixed at 1,
#             which is useful when there are 3 or less items. Default is "variance"
# missing: options to handle missing values available in lavaan. Default is "listwise"
# ...: other arguments to input in lavaan.

esm_pfa_omega <- function(data, idlab, items, 
                          constraint = "variance", 
                          missing = "listwise",
                          ...) {
  
  # Turn data to list
  data_list <- esm2list(data, id = data[, idlab])
  
  # Create lavaan pfa model
  pfa_model <- pfa_syntax(items, constraint = constraint)
  
  # P-technique loop
  # This loop fit the p-technique to each person saving the lavaan output,
  # the estimated omega, and whether there is a reasonable goodness of fit
  # to each person according to the CFI, TLI, and RMSEA.
  
  pfa <- pfa_omg <- pfa_fit <- list()
  
  for (p in 1:length(data_list)) {
    cat(sprintf("This is %d-th participant in the data list.\n", p))
    tmp_dat <- data_list[[p]][, items]
    
    if (!any(eigen(cov(tmp_dat, use = "pairwise"))$values <= 0) &
        all(diag(var(tmp_dat, use = "pairwise")) > 0)) {
      pfa[[p]] <- lavaan::cfa(pfa_model, tmp_dat, missing = missing, ...)
      
      if (lavaan::lavInspect(pfa[[p]], what = "converged") &
          lavaan::lavInspect(pfa[[p]], what = "post.check")) {
        tmp_est <- lavaan::parameterEstimates(pfa[[p]])
        pfa_omg[[p]] <- tmp_est$est[tmp_est$label == "omega"]
        pfa_fit[[p]] <- all(c(
          lavaan::fitMeasures(pfa[[p]])[c("cfi", "tli")] > 0.9,
          lavaan::fitMeasures(pfa[[p]])["rmsea"] < 0.08))
      } else {
        pfa_omg[[p]] <- NA
        pfa_fit[[p]] <- NA
      }
    } else {
      pfa[[p]]  <- NULL
      pfa_omg[[p]] <- NA
      pfa_fit[[p]] <- NA
    }
  }
  rm(tmp_dat, tmp_est, p)
  
  person_output <- data.frame(names(data_list), 
                              unlist(pfa_omg), 
                              unlist(pfa_fit))
  names(person_output) <- c(idlab, "pfa_omega", "pfa_fit")
  
  out <- list(pfa_omega   = person_output, 
              data_list   = data_list,
              lavaan_list = pfa)
  
  return(out)
}

# esm_dfa_omega ----

# Function to fit and unidimensional dynamic factor model per person in dynr and 
# to compute person-specific omega.

# data: data frame with the ESM data including the id variable and the time variable
# idlab: character indicating the name of the id variable in the data.
# timelab: character indicating the name of the time variable in the data.
# items: character with the name of the item variables
# ...: other arguments to input in lavaan.

esm_dfa_omega <- function(data, idlab, timelab, items, ...){
  
  require(dynr)
  
  I <- length(items)
  
  # Fit p-technique per person to get starting values
  
  cat("Getting starting values from P-technique factor analysis\n\n")

  pfa_fit <- suppressWarnings(
    esm_pfa_omega(data          = data,
                  idlab         = idlab,
                  items         = items,
                  constraint    = "loading",
                  meanstructure = TRUE,
                  ...)
  )
    
  # Get data list 
  
  data_list <- pfa_fit$data_list
  
  # Get starting values
  
  pfa_start <- list()
  
  for (p in 1:length(data_list)) {
    if (!is.na(pfa_fit$pfa_omega$pfa_omega[p])) {
      pfa_start[[p]] <- 
        lavaan::parameterEstimates(pfa_fit$lavaan_list[[p]])[1:(3 * I + 1), "est"]
    } else {
      pfa_start[[p]] <- rep(c(1, 0.5, 0), times = c(I + 1, I, I))
    }
  }
  rm(p)
  
  # Define recipe components to input into dynr
  
  dynamics <- meas <- ecov <- initial <- list()
  
  for (p in 1:length(data_list)) {
    # Dynamic component. Defining Autoregressive effect
    dynamics[[p]] <- dynr::prep.matrixDynamics(
      values.dyn = matrix(c(0), nrow = 1), # AR effect initial value set to 0
      params.dyn = matrix(c("phi"), nrow = 1), # Autoregressive effect label
      isContinuousTime = FALSE)
    # Measurement component. Defining factor model.
    meas[[p]] <- dynr::prep.measurement(
      values.load = matrix(pfa_start[[p]][1:I], I, 1), # Loadings initial values
      params.load = matrix(c('fixed', paste0("lambda_", 2:I)), I, 1), # labels
      values.int = matrix(pfa_start[[p]][(2 * I + 2):(3 * I + 1)], I, 1), # Intercepts
      params.int = matrix(paste0("int", 1:I), I, 1), #labels
      state.names = c("Theta"), # Latent variable name
      obs.names   = items) # Items
    # Noise component. Defining measurement error and innovation error terms.
    ecov[[p]] <- dynr::prep.noise(
      values.latent = dynr::diag(pfa_start[[p]][I + 1], 1), # Innovation Variance initial value
      params.latent = dynr::diag(c("innovar"), 1), # Innovation variance label
      values.observed = dynr::diag(pfa_start[[p]][(I + 2):(2 * I + 1)], I), # Measurement error variance initial value
      params.observed = dynr::diag(paste0("merror_", 1:I), I)) # Measurement error labels
    # Initial states component. Defining number of states and initial conditions.
    initial[[p]] <- dynr::prep.initial(
      values.inistate = c(0), # Initial mean of the latent variable
      params.inistate = c('fixed'), # Fixed mean of latent variable
      values.inicov   = dynr::diag(1, 1), # Initial variance
      params.inicov   = dynr::diag('fixed', 1)) # Fixed variance
  }
  rm(p)
  
  # Fit unidimensional DFA per person with dynr
  
  dfa_fit  <- dfa_omega <- dfa_ar <- list()
  
  cat("\n\nFitting the unidimensional dynamic factor model in 'dynr'.\n\n")
   
  for (p in 1:length(data_list)) {
    cat(sprintf("\nThis is %d-th participant in the data list.\n", p))
    # Get data for participant p
    tmp <- data_list[[p]]
    # Structure data with dynr.data
    data_dynr <- dynr::dynr.data(
      tmp,
      id = idlab,
      time = timelab,
      observed = items)
    # Specify dfa model in dynr
    model_dynr_dfa <- dynr::dynr.model(
      dynamics    = dynamics[[p]],
      measurement = meas[[p]],
      noise       = ecov[[p]],
      initial     = initial[[p]],
      data        = data_dynr)
    # Estimate the dynamic factor model and compute omega
    cook_dfa <- dynr::dynr.cook(model_dynr_dfa, verbose = FALSE)
    # Save results if model "converged" (exitflag 1, 2, or 3) and if
    # the estimated autoregressive effect is lower than 1.
    if (cook_dfa$exitflag %in% c(1:3) && coef(cook_dfa)["phi"] < 1) {
      # Save dynr object
      dfa_fit[[p]] <- cook_dfa
      # Save estimated omega
      dfa_omega[[p]] <-
        (sum(c(1, coef(cook_dfa)[paste0("lambda_", 2:I)]))^2 *
           (coef(cook_dfa)["innovar"]/(1 - coef(cook_dfa)["phi"]^2))) /
        (sum(c(1, coef(cook_dfa)[paste0("lambda_", 2:I)]))^2 *
           (coef(cook_dfa)["innovar"]/(1 - coef(cook_dfa)["phi"]^2)) +
           sum(coef(cook_dfa)[paste0("merror_", 1:I)]))
      # Save autoregressive effect
      dfa_ar[[p]] <- coef(cook_dfa)["phi"]
    } else {
      dfa_fit[[p]]   <- NULL
      dfa_omega[[p]] <- NA
      dfa_ar[[p]]    <- NA
    }
    # Clear temporal objects to avoid overwriting
    rm(tmp, data_dynr, model_dynr_dfa,
       cook_dfa)
  }
  rm(p)
  
  person_output <- data.frame(names(data_list), 
                              unlist(dfa_omega), 
                              unlist(dfa_ar))
  names(person_output) <- c(idlab, "dfa_omega", "dfa_ar")
  
  out <- list(dfa_omega = person_output, 
              data_list = data_list,
              dynr_list = dfa_fit)
  
  return(out)
}





















