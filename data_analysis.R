fit <- function(
    data,                  # standardized dataset
    family,                #
    rve_type,              # RVE type: MD was recommended 
    ss_correct,            # Whether small sample correction will be applied: T/F
    offset = FALSE,
    design                 # cs = cross-section; co = cohort 
    #re
){
  
  results <- list()
  
  if (offset == F) {f_offset <- ""} else {
    data$off <- as.numeric(as.character(data[[offset]]))
    f_offset <- " + offset(log(off))"
  }
  
  fre1 <- " + (1|Cluster) + (1|Cluster:Period) + (1|id_individual)"
  fre2 <- " + (1|Cluster) + (1|id_individual)"
  fre3 <- " + (1|Cluster) + (1|Cluster:Period)"
  fre4 <- " + (1|Cluster)"
  
  ################################################.
  #####       Immediate Treatment (IT)       #####
  ################################################
  
  formula1 <- paste0("Outcome ~ factor(Period) + Treatment", fre1, f_offset)
  formula2 <- paste0("Outcome ~ factor(Period) + Treatment", fre2, f_offset)
  formula3 <- paste0("Outcome ~ factor(Period) + Treatment", fre3, f_offset)
  formula4 <- paste0("Outcome ~ factor(Period) + Treatment", fre4, f_offset)
  
  if (design == "cs") {
    if(family == "gaussian"){
      model1 <- try(lme4::lmer(formula3, data = data), silent = TRUE)
      model2 <- try(lme4::lmer(formula4, data = data), silent = TRUE)
    } else {
      model1 <- try(glmer(formula3, family = family, data = data), silent = TRUE)
      model2 <- try(glmer(formula4, family = family, data = data), silent = TRUE)
    }
    
    
    if(is.null(model1) == T){
      itm1_result <- NULL
    } else {
      itm1_result <- get_coef(model1, rve_type, ss_correct)
    }
    
    if(is.null(model1) == T){
      itm2_result <- NULL
    } else {
      itm2_result <- get_coef(model2, rve_type, ss_correct)
    }
  } 
  if (design == "co"){
    if(family == "gaussian"){
      model1 <- try(lme4::lmer(formula1, data = data), silent = TRUE)
      model2 <- try(lme4::lmer(formula2, data = data), silent = TRUE)
      model3 <- try(lme4::lmer(formula4, data = data), silent = TRUE)
    } else {
      model1 <- try(glmer(formula1, family = family, data = data), silent = TRUE)
      model2 <- try(glmer(formula2, family = family, data = data), silent = TRUE)
      model3 <- try(glmer(formula4, family = family, data = data), silent = TRUE)
    }
    if(is.null(model1) == T){
      itm1_result <- NULL
    } else {
      itm1_result <- get_coef(model1, rve_type, ss_correct)
    }
    
    if(is.null(model2) == T){
      itm2_result <- NULL
    } else {
      itm2_result <- get_coef(model2, rve_type, ss_correct)
    }
    
    if(is.null(model3) == T){
      itm3_result <- NULL
    } else {
      itm3_result <- get_coef(model3, rve_type, ss_correct)
    }
  }
  
  ################################################.
  #####    Exposure Time Indicator (ETI)     #####
  ################################################.  
  formula1 <- paste0("Outcome ~ Period + factor(Exposure)", fre1, f_offset)
  formula2 <- paste0("Outcome ~ Period + factor(Exposure)", fre2, f_offset)
  formula3 <- paste0("Outcome ~ Period + factor(Exposure)", fre3, f_offset)
  formula4 <- paste0("Outcome ~ Period + factor(Exposure)", fre4, f_offset)
  
  if (design == "cs"){
    if(family == "gaussian"){
      model1 <- try(lme4::lmer(formula3, data = data), silent = TRUE)
      model2 <- try(lme4::lmer(formula4, data = data), silent = TRUE)
    } else {
      model1 <- try(glmer(formula3, family = family, data = data), silent = TRUE)
      model2 <- try(glmer(formula4, family = family, data = data), silent = TRUE)
    }
    
    if(is.null(model1) == T){
      etim1_result <- NULL
    } else {
      etim1_result <- get_eticoef(model1, rve_type, ss_correct)
    }
    
    if(is.null(model1) == T){
      etim2_result <- NULL
    } else {
      etim2_result <- get_eticoef(model2, rve_type, ss_correct)
    }
  }
  if (design == "co"){
    if(family == "gaussian"){
      model1 <- try(lme4::lmer(formula1, data = data), silent = TRUE)
      model2 <- try(lme4::lmer(formula2, data = data), silent = TRUE)
      model3 <- try(lme4::lmer(formula4, data = data), silent = TRUE)
    } else {
      model1 <- try(glmer(formula1, family = family, data = data), silent = TRUE)
      model2 <- try(glmer(formula2, family = family, data = data), silent = TRUE)
      model3 <- try(glmer(formula4, family = family, data = data), silent = TRUE)
    }
    if(is.null(model1) == T){
      etim1_result <- NULL
    } else {
      etim1_result <- get_eticoef(model1, rve_type, ss_correct)
    }
    
    if(is.null(model2) == T){
      etim2_result <- NULL
    } else {
      etim2_result <- get_eticoef(model2, rve_type, ss_correct)
    }
    
    if(is.null(model3) == T){
      etim3_result <- NULL
    } else {
      etim3_result <- get_eticoef(model3, rve_type, ss_correct)
    }
  }
  
  
  #####################################################.
  #####   Treatment Effect Heterogeneity (TEH)   #####
  #####################################################.
  formula1 <- paste0("Outcome ~ Period + Treatment + (0+Treatment|Exposure) ", fre1, f_offset)
  formula2 <- paste0("Outcome ~ Period + Treatment + (0+Treatment|Exposure) ", fre2, f_offset)
  formula3 <- paste0("Outcome ~ Period + Treatment + (0+Treatment|Exposure) ", fre3, f_offset)
  formula4 <- paste0("Outcome ~ Period + Treatment + (0+Treatment|Exposure) ", fre4, f_offset)
  
  
  
  if (design == "cs"){
    if(family == "gaussian"){
      model1 <- try(lme4::lmer(formula3, data = data), silent = TRUE)
      model2 <- try(lme4::lmer(formula4, data = data), silent = TRUE)
    } else {
      model1 <- try(glmer(formula3, family = family, data = data), silent = TRUE)
      model2 <- try(glmer(formula4, family = family, data = data), silent = TRUE)
    }
    
    if(is.null(model1) == T){
      tehm1_result <- NULL
    } else {
      tehm1_result <- get_tehcoef(model1, ss_correct)
    }
    
    if(is.null(model1) == T){
      tehm2_result <- NULL
    } else {
      tehm2_result <- get_tehcoef(model2, ss_correct)
    }
  }
  if (design == "co"){
    if(family == "gaussian"){
      model1 <- try(lme4::lmer(formula1, data = data), silent = TRUE)
      model2 <- try(lme4::lmer(formula2, data = data), silent = TRUE)
      model3 <- try(lme4::lmer(formula4, data = data), silent = TRUE)
    } else {
      model1 <- try(glmer(formula1, family = family, data = data), silent = TRUE)
      model2 <- try(glmer(formula2, family = family, data = data), silent = TRUE)
      model3 <- try(glmer(formula4, family = family, data = data), silent = TRUE)
    }
    if(is.null(model1) == T){
      tehm1_result <- NULL
    } else {
      tehm1_result <- get_tehcoef(model1, ss_correct)
    }
    
    if(is.null(model2) == T){
      tehm2_result <- NULL
    } else {
      tehm2_result <- get_tehcoef(model2, ss_correct)
    }
    
    if(is.null(model3) == T){
      tehm3_result <- NULL
    } else {
      tehm3_result <- get_tehcoef(model3, ss_correct)
    }
  }
  
  ############################################.
  #####    Natural Cubic Spline (NCS)    #####
  ############################################.
  
  
  ## Cubic spline
  data$Period <- as.numeric(as.character(data$Period))
  J <- length(unique(data$Exposure))
  nnode <-  ceiling((J)/2)
  ns_basis <- ns(c(0:(J-1)), knots = (1:(nnode-1))*(J-1)/(nnode))
  
  for (i in 1:(nnode)) {
    new_vector <- ns_basis[data$Exposure+1,i]  # Example: random numbers
    data[[paste0("b", i)]] <- new_vector
  }
  
  b_vars <- paste0("b", 1:(nnode))
  formula1 <- as.formula(paste("Outcome ~ factor(Period) + ", paste(b_vars, collapse = " + "), fre1 ))
  formula2 <- as.formula(paste("Outcome ~ factor(Period) + ", paste(b_vars, collapse = " + "), fre2 ))
  formula3 <- as.formula(paste("Outcome ~ factor(Period) + ", paste(b_vars, collapse = " + "), fre3 ))
  formula4 <- as.formula(paste("Outcome ~ factor(Period) + ", paste(b_vars, collapse = " + "), fre4 ))
  
  if (design == "cs"){
    if(family == "gaussian"){
      model1 <- try(lme4::lmer(formula3, data = data), silent = TRUE)
      model2 <- try(lme4::lmer(formula4, data = data), silent = TRUE)
    } else {
      model1 <- try(glmer(formula3, family = family, data = data), silent = TRUE)
      model2 <- try(glmer(formula4, family = family, data = data), silent = TRUE)
    }
    
    if(is.null(model1) == T){
      ncsm1_result <- NULL
    } else {
      ncsm1_result <- get_ncscoef(model1, data, rve_type, ss_correct, ns_basis, J, nnode)
    }
    
    if(is.null(model1) == T){
      ncsm2_result <- NULL
    } else {
      ncsm2_result <- get_ncscoef(model2, data, rve_type, ss_correct, ns_basis, J, nnode)
    }
  }
  if (design == "co"){
    if(family == "gaussian"){
      model1 <- try(lme4::lmer(formula1, data = data), silent = TRUE)
      model2 <- try(lme4::lmer(formula2, data = data), silent = TRUE)
      model3 <- try(lme4::lmer(formula4, data = data), silent = TRUE)
    } else {
      model1 <- try(glmer(formula1, family = family, data = data), silent = TRUE)
      model2 <- try(glmer(formula2, family = family, data = data), silent = TRUE)
      model3 <- try(glmer(formula4, family = family, data = data), silent = TRUE)
    }
    if(is.null(model1) == T){
      ncsm1_result <- NULL
    } else {
      ncsm1_result <- get_ncscoef(model1, data, rve_type, ss_correct, ns_basis, J, nnode)
    }
    
    if(is.null(model2) == T){
      ncsm2_result <- NULL
    } else {
      ncsm2_result <- get_ncscoef(model2, data, rve_type, ss_correct, ns_basis, J, nnode)
    }
    
    if(is.null(model3) == T){
      ncsm3_result <- NULL
    } else {
      ncsm3_result <- get_ncscoef(model3, data, rve_type, ss_correct, ns_basis, J, nnode)
    }
  }
  if (design == "cs"){
    results <- list(itm1=itm1_result, itm2=itm2_result, etim1 = etim1_result, etim2 = etim2_result, 
                    tehm1 = tehm1_result, tehm2 = tehm2_result, ncsm1 = ncsm1_result, ncsm2 = ncsm2_result)
  }
  if (design == "co"){
    results <- list(itm1=itm1_result, itm2=itm2_result, itm3=itm3_result, 
                    etim1 = etim1_result, etim2 = etim2_result, etim3 = etim3_result, 
                    tehm1 = tehm1_result, tehm2 = tehm2_result, tehm3 = tehm3_result,
                    ncsm1 = ncsm1_result, ncsm2 = ncsm2_result, ncsm3 = ncsm3_result)
  }
  return(results)
}

  


#testx <- fit(data = data, family = "gaussian", rve_type = "CR0", ss_correct=T, design = "co")
