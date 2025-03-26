## Helper function: calculate CI with t-distribution
cal_confint_t <- function(pe, df, se){
  cihigh <- pe + qt(0.975, df)*se
  cilow <- pe - qt(0.975, df)*se
  return(cbind(cilow, cihigh))
}


## format RE estimates for presentation
format_reest <- function(reest) {
  random_effects <- as.data.frame(reest)
  
  # Process the random effects table
  random_effects_table <- random_effects %>%
    dplyr::select(grp = grp, name = var1, std_dev = sdcor, variance = vcov) %>%
    dplyr::mutate(
      grp = as.character(grp),
      grp = dplyr::recode(grp, "ij" = "Cluster-period"),  # Rename groups
      name = ifelse(is.na(name), "(Intercept)", name),
      variance_str = ifelse(is.na(variance), "", format(round(variance, 3), nsmall = 3))
    )
  
  random_effect_strings <- apply(random_effects_table, 1, function(row) {
    paste0(
      row["grp"], ": ", row["name"], "\n",
      "  Std.Dev = ", format(as.numeric(row["std_dev"]), digits = 3),
      if (row["variance_str"] != "") paste0(", var = ", row["variance_str"]) else ""
    )
  })
  
  # Collapse all group blocks into a single plain text string with newlines
  random_effect_summary <- paste(random_effect_strings, collapse = "\n\n")
  
  return(random_effect_summary)
}



## Function to process it models
get_coef <- function(model, rse_type, ss_correct){
  model_summary <- summary(model$model)
  model_vcov <- vcov(model$model)
  model_reest <-model_summary$varcor
  model_est <- model_summary$coefficients["treatment",1]
  model_se <- model_summary$coefficients["treatment",2]
  
  rse.matrix <- model_rse <- vcovCR.glmerMod(model$model, type = rse_type)
  model_rse <- sqrt(model_rse["treatment", "treatment"])
  
  ### Construct CI depending on SS correction
  if (ss_correct == T){
    dof <- model_summary$ngrps[["cluster_id"]] - 2
    modelci <- cal_confint_t(pe = model_est, df = dof, se = model_se)
    modelci_rse <- cal_confint_t(pe = model_est, df = dof, se = model_rse)
  } else {
    modelci <- model_est + c(-1.96,1.96) * model_se
    modelci_rse <- model_est + c(-1.96,1.96) * model_rse
  }
  list <- list(
    model = model$model,
    reest = format_reest(model_reest),
    est = model_est,
    se = model_se,
    rse = model_rse,
    ci = modelci,
    rseci = modelci_rse,
    lte = model_est,
    lte_se = model_se,
    lte_rse = model_rse,
    lteci = modelci,
    lterseci = modelci_rse,
    vcov_rv = rse.matrix,
    converged = performance::check_convergence(model$model)[1],
    lme4_converged = ifelse(model$model@optinfo$conv$opt==0, TRUE, FALSE),
    messages = model$model@optinfo$conv$lme4$messages
  )
  return(list)
}

## Function to process eti models
get_eticoef <- function(model, rse_type, ss_correct){
  model_summary <- summary(model$model)
  model_reest <- model_summary$varcor
  indices <- grep("exp_", rownames(model_summary$coefficients))
  index_max <- length(indices)
  lte <- model_summary$coefficients[indices,][index_max,1]
  lte_se <- model_summary$coefficients[indices,][index_max,2]
  
  
  coeffs <- model_summary$coefficients[,1][indices] # column 1 contains the estimates
  sigma.matrix <- vcov(model_summary)[indices,indices]
  
  A <- matrix(rep(1/index_max), index_max, nrow=1)
  model_est <- (A %*% coeffs)[1]
  model_se <- (sqrt(A %*% sigma.matrix %*% t(A)))[1,1]
  rse.matrix <- model_rse <- vcovCR.glmerMod(model$model, type = rse_type)
  lte_rse <- model_rse[indices,indices][index_max, index_max]
  sigmarse.matrix <- model_rse[indices,indices]
  model_rse <- (sqrt(A %*% sigmarse.matrix %*% t(A)))[1,1]
  
  
  ### obtain CI depending on SS correction
  if (ss_correct == T){
    dof <- model_summary$ngrps[["cluster_id"]] - 2
    modelci <- cal_confint_t(pe = model_est, df = dof, se = model_se)
    modelci_rse <- cal_confint_t(pe = model_est, df = dof, se = model_rse)
    lteci <- cal_confint_t(pe = lte, df = dof, se = lte_se)
    lteci_rse <- cal_confint_t(pe = lte, df = dof, se = lte_rse)
  } else {
    modelci <- model_est + c(-1.96,1.96) * model_se
    modelci_rse <- model_est + c(-1.96,1.96) * model_rse
    lteci <- lte + c(-1.96,1.96) * lte_se
    lteci_rse <- lte + c(-1.96,1.96) * lte_rse
  }
  list <- list(
    model = model$model,
    reest = format_reest(model_reest),
    est = model_est,
    se = model_se,
    rse = model_rse,
    ci = modelci,
    rseci = modelci_rse,
    lte = lte,
    lte_se = lte_se,
    lte_rse = lte_rse,
    lteci = lteci,
    lterseci = lteci_rse,
    vcov_rv = rse.matrix,
    converged = performance::check_convergence(model$model)[1],
    lme4_converged = ifelse(model$model@optinfo$conv$opt==0, TRUE, FALSE),
    messages = model$model@optinfo$conv$lme4$messages
  )
  return(list)
}


## Function to process teh models
get_tehcoef <- function(model, ss_correct){
  model_summary <- summary(model$model)
  model_reest <- model_summary$varcor
  model_vcov <- vcov(model$model)
  model_est <- model_summary$coefficients["treatment",1]
  model_se <- model_summary$coefficients["treatment",2]
  
  ### Obtain rse for estimated trt effects
  #modelrse <- vcovCR(model, type= rse_type)
  #modelrse_trt <- sqrt(modelrse["Treatment", "Treatment"])
  if (ss_correct == T){
    dof <- model_summary$ngrps[["cluster_id"]] - 2
    modelci <- cal_confint_t(pe = model_est, df = dof, se = model_se)
  } else {
    modelci <- model_est + c(-1.96,1.96) * model_se
  }
  list <- list(
    model = model$model,
    lte = "NA",
    reest = format_reest(model_reest),
    est = model_est,
    se = model_se,
    rse = "NA",
    ci = modelci,
    ci_rse = "NA",
    lte = "NA",
    lte_se = "NA",
    lte_rse = "NA",
    lteci = "NA",
    lterseci = "NA",
    converged = performance::check_convergence(model$model)[1],
    lme4_converged = ifelse(model$model@optinfo$conv$opt==0, TRUE, FALSE),
    messages = model$model@optinfo$conv$lme4$messages
  )
  return(list)
}


## Function to process ncs models
get_ncscoef <- function(model, data, rse_type, ss_correct){
  
  model_summary <- summary(model$model)
  model_reest <- model_summary$varcor

  model_est <- model$te_est
  model_se <- model$te_se
  
  # Still working on implementing the robust standard error
  S <- max(data$exposure_time)
  knots_exp <- seq(0, S, length.out=n_knots)
  ns_basis <- splines::ns(
    x = data$exposure_time,
    knots = knots_exp[2:(n_knots-1)],
    intercept = TRUE,
    Boundary.knots = knots_exp[c(1,n_knots)]
  )
  
  # Specify the indices corresponding to the spline terms
  indices <- grep("^b[0-9]+$", rownames(model_summary$coefficients))
  index_max <- length(indices)
  
  # Extract coefficient estimates and covariance matrix corresponding to spline
  # terms
  coeffs_spl <-  model_summary$coefficients[,1][indices]
  
  # Get number of unique (non-zero) exposure times
  exp_timepoints <- unique(data$exposure_time[data$exposure_time != 0])
  num_exp_timepoints <- length(exp_timepoints)
  max_exp_timepoint <- max(exp_timepoints)
  
  # Transform the spline terms into effect curve estimates (+ covariance matrix)
  B <- as.matrix(splines::ns(
    x = c(1:S),
    knots = knots_exp[2:(n_knots-1)],
    intercept = TRUE,
    Boundary.knots = knots_exp[c(1,n_knots)]
  ))
  
  class(B) <- "matrix"
  
  rse.matrix <- model_rse <- vcovCR.glmerMod(model$model, type = rse_type)
  sigmarse.matrix <- model_rse[indices,indices]
  sigmarse.matrix <- B %*% sigmarse.matrix %*% t(B)
  rse_ncs <- sqrt(diag(matrix(sigmarse.matrix, nrow = nrow(sigmarse.matrix))))
  
  
  num_estimand_timepoints <- max(exp_timepoints)
  M <- matrix(rep(1/num_exp_timepoints, num_exp_timepoints), nrow=1)
  model_rse <- (sqrt(M %*% sigmarse.matrix %*% t(M)))[1,1]
  

  
  #se_ncs <- sqrt(diag(matrix(cov_mtx, nrow = nrow(cov_mtx))))
  lte_ind <- which.max(model$effect_curve$exp_time)
  lte <- model$effect_curve$est[lte_ind]
  lte_se <-  model$effect_curve$se[lte_ind]
  lte_rse <- rse_ncs[lte_ind]
  ### obtain CI depending on SS correction
  if (ss_correct == T){
    dof <- model_summary$ngrps[["cluster_id"]] - 2
    modelci <- cal_confint_t(pe = model_est, df = dof, se = model_se)
    modelci_rse <- cal_confint_t(pe = model_est, df = dof, se = model_rse)
    lteci <- cal_confint_t(pe = lte, df = dof, se = lte_se)
    lteci_rse <- cal_confint_t(pe = lte, df = dof, se = lte_rse)
  } else {
    modelci <- model_est + c(-1.96,1.96) * model_se
    modelci_rse <- model_est + c(-1.96,1.96) * model_rse
    lteci <- lte + c(-1.96,1.96) * lte_se
    lteci_rse <- lte + c(-1.96,1.96) * lte_rse
  }
  
  list <- list(
    model = model$model,
    reest = format_reest(model_reest),
    est = model_est,
    se = model_se,
    rse = model_rse,
    ci = modelci,
    rseci = modelci_rse,
    lte = lte,
    lte_se = lte_se,
    lte_rse = lte_rse,
    lteci = lteci,
    lterseci = lteci_rse,
    vcov_rv = rse.matrix,
    converged = performance::check_convergence(model$model)[1],
    lme4_converged = ifelse(model$model@optinfo$conv$opt==0, TRUE, FALSE),
    messages = model$model@optinfo$conv$lme4$messages
  )
  return(list)
}



