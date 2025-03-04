## Helper function: calculate CI with t-distribution
cal_confint_t <- function(pe, df, se){
  cihigh <- pe + qt(0.975, df)*se
  cilow <- pe - qt(0.975, df)*se
  return(cbind(cilow, cihigh))
}


## format RE estimates for presentation
format_reest <- function(reest){
  random_effects <- as.data.frame(reest)
  
  # Process the random effects table
  random_effects_table <- random_effects %>%
    dplyr::select(grp = grp, name = var1, std_dev = sdcor, var = vcov) %>%
    dplyr::mutate(
      grp = as.character(grp),
      name = ifelse(is.na(name), "(Intercept)", name),  # Handle missing variable names
      var = ifelse(is.na(var), "", format(round(var, 3), nsmall = 3)) # Handle missing correlations
    )
  
  random_effect_strings <- apply(random_effects_table, 1, function(row) {
    paste0(
      row["grp"], ": ", row["name"], "<br>",
      "&nbsp;&nbsp;Std.Dev = ", format(as.numeric(row["std_dev"]), digits = 3),
      if (row["var"] != "") paste0(", var = ", row["var"]) else "",
      "<br>"
    )
  })
  random_effect_summary  <- random_effect_strings
  return(random_effect_summary)
}




F## Function to process it models
get_coef <- function(model, rse_type, ss_correct){
  model_summary <- summary(model)
  model_vcov <- vcov(model)
  model_reest <-model_summary$varcor
  model_est <- model_summary$coefficients["Treatment",1]
  model_se <- model_summary$coefficients["Treatment",2]
  
  model_rse <- vcovCR(model, type = rse_type)
  model_rse <- sqrt(model_rse["Treatment", "Treatment"])
  
  ### Construct CI depending on SS correction
  if (ss_correct == T){
    dof <- model_summary$ngrps[["Cluster"]] - 2
    modelci <- cal_confint_t(pe = model_est, df = dof, se = model_se)
    modelci_rse <- cal_confint_t(pe = model_est, df = dof, se = model_rse)
  } else {
    modelci <- model_est + c(-1.96,1.96) * model_se
    modelci_rse <- model_est + c(-1.96,1.96) * model_rse
  }
  list <- list(
    model = model,
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
    converged = performance::check_convergence(model)[1],
    lme4_converged = ifelse(model@optinfo$conv$opt==0, TRUE, FALSE),
    messages = model@optinfo$conv$lme4$messages
  )
  return(list)
}

## Function to process eti models
get_eticoef <- function(model, rse_type, ss_correct){
  model_summary <- summary(model)
  model_reest <-model_summary$varcor
  indices <- grep("Exposure", rownames(model_summary$coefficients))
  index_max <- length(indices)
  lte <- model_summary$coefficients[indices,][index_max,1]
  lte_se <- model_summary$coefficients[indices,][index_max,2]
  
  
  coeffs <- model_summary$coefficients[,1][indices] # column 1 contains the estimates
  sigma.matrix <- vcov(model_summary)[indices,indices]
  
  A <- matrix(rep(1/index_max), index_max, nrow=1)
  model_est <- (A %*% coeffs)[1]
  model_se <- (sqrt(A %*% sigma.matrix %*% t(A)))[1,1]
  model_rse <- vcovCR(model, type = rse_type)
  lte_rse <- model_rse[indices,indices][index_max, index_max]
  sigmarse.matrix <- model_rse[indices,indices]
  model_rse <- (sqrt(A %*% sigmarse.matrix %*% t(A)))[1,1]
  
  
  ### obtain CI depending on SS correction
  if (ss_correct == T){
    dof <- model_summary$ngrps[["Cluster"]] - 2
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
    model = model,
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
    converged = performance::check_convergence(model)[1],
    lme4_converged = ifelse(model@optinfo$conv$opt==0, TRUE, FALSE),
    messages = model@optinfo$conv$lme4$messages
  )
  return(list)
}


## Function to process teh models
get_tehcoef <- function(model, ss_correct){
  model_summary <- summary(model)
  model_reest <- model_summary$varcor
  model_vcov <- vcov(model)
  model_est <- model_summary$coefficients["Treatment",1]
  model_se <- model_summary$coefficients["Treatment",2]
  
  ### Obtain rse for estimated trt effects
  #modelrse <- vcovCR(model, type= rse_type)
  #modelrse_trt <- sqrt(modelrse["Treatment", "Treatment"])
  if (ss_correct == T){
    dof <- model_summary$ngrps[["Cluster"]] - 2
    modelci <- cal_confint_t(pe = model_est, df = dof, se = model_se)
  } else {
    modelci <- model_est + c(-1.96,1.96) * model_se
  }
  list <- list(
    model = model,
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
    converged = performance::check_convergence(model)[1],
    lme4_converged = ifelse(model@optinfo$conv$opt==0, TRUE, FALSE),
    messages = model@optinfo$conv$lme4$messages
  )
  return(list)
}


## Function to process ncs models
get_ncscoef <- function(model, data, rse_type, ss_correct, ns_basis, J, nnode){
  model_summary <- summary(model)
  model_reest <- model_summary$varcor
  # Specify the indices corresponding to the spline terms
  indices <- grep("^b[0-9]+$", rownames(model_summary$coefficients))
  index_max <- length(indices)
  
  coeffs_b <- model_summary$coefficients[,1][indices]
  sigma.matrix_b <- stats::vcov(model)[indices,indices]
  
  # Get number of unique (non-zero) exposure times
  exp_timepoints <- unique(data$Exposure[data$Exposure != 0])
  num_exp_timepoints <- length(exp_timepoints)
  max_exp_timepoint <- max(exp_timepoints)
  
  B <- matrix(NA, nrow=(J-1), ncol = nnode)
  for (i in 1:(J-1)) {
    for (j in 1:(nnode)) {
      B[i,j] <- ns_basis[i+1,j]
    }
  }
  
  coeffs_b_new <- as.numeric(B %*% coeffs_b)
  sigma.matrix <- B %*% sigma.matrix_b %*% t(B)
  model_rse <- vcovCR(model, type = rse_type)
  sigmarse.matrix <- model_rse[indices,indices]
  sigmarse.matrix <- B %*% sigmarse.matrix %*% t(B)
  
  
  A <- matrix(rep(1/num_exp_timepoints, num_exp_timepoints), nrow=1)
  model_est <- (A %*% coeffs_b_new)[1]
  model_se <- (sqrt(A %*% sigma.matrix %*% t(A)))[1,1]
  model_rse <- (sqrt(A %*% sigmarse.matrix %*% t(A)))[1,1]
  
  
  ### obtain CI depending on SS correction
  if (ss_correct == T){
    dof <- model_summary$ngrps[["Cluster"]] - 2
    modelci <- cal_confint_t(pe = model_est, df = dof, se = model_se)
    modelci_rse <- cal_confint_t(pe = model_est, df = dof, se = model_rse)
  } else {
    modelci <- model_est + c(-1.96,1.96) * model_se
    modelci_rse <- model_est + c(-1.96,1.96) * model_rse
  }
  
  list <- list(
    model = model,
    lte = "NA",
    reest = format_reest(model_reest),
    est = model_est,
    se = model_se,
    rse = model_rse,
    ci = modelci,
    rseci = modelci_rse,
    lte = "NA",
    lte_se = "NA",
    lte_rse = "NA",
    lteci = "NA",
    lterseci = "NA",
    converged = performance::check_convergence(model)[1],
    lme4_converged = ifelse(model@optinfo$conv$opt==0, TRUE, FALSE),
    messages = model@optinfo$conv$lme4$messages
  )
  return(list)
}



