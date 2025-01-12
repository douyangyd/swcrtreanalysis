# Analyze the Poisson outcome with cohort design
fit.data.binary.CO <- function(
    data,                  # standardized dataset
    rve_type,              # RVE type: MD was recommended 
    ss_correct             # Whether small sample correction will be applied: T/F
){
  results <- list()
  
  ## Immediate time model (IT)
  ### Model fitting
  model1 <- try(glmer(Outcome ~  Period + Treatment + (1|Cluster) + (1|Cluster:Period) + (1|id_individual), 
                      family = "poisson", data = data), silent = TRUE)
  
  # Check if model1 converged
  if (check_convergence(model1)[1] == F) {
    # Fit the next model if model1 failed to converge
    model2 <- try(glmer(Outcome ~  Period + Treatment + (1|Cluster) + (1|Cluster:Period), 
                        family = "poisson", data = data), silent = TRUE)
    
    # Check if model2 also failed to converge
    if (check_convergence(model2)[1] == F) {
      model3 <- try(glmer(Outcome ~  Period + Treatment + (1|Cluster), 
                          family = "poisson", data = data), silent = TRUE)
      if(check_convergence(model3)[1] == F) {
        ITconverge <- 0
        message("All models failed to converge.")
      } else {
        IT <- model3
        ITconverge <- 3
        message("model3 converged successfully.")
      }
    } else {
      IT <- model2
      ITconverge <- 2
      message("model2 converged successfully.")
    }
  } else {
    IT <- model1
    ITconverge <- 1
    message("model1 converged successfully.")
  }
  
  
  if(ITconverge == 0){
    IT_est <- IT_vcov <- ITRVE <- IT_coef <- ITRVE_trt <- ITCI <- ITCI_RVE <- NULL 
  } else {
    ### Obtain estimated coef and SE for the model
    IT_summary <- summary(IT)
    IT_est <- summary(IT)$coefficients
    IT_vcov <- vcov(IT)
    IT_coef <- IT_summary$coefficients["Treatment",]
    
    ### Obtain RVE for estimated trt effects
    ITRVE <- vcovCR(IT, type= rve_type)
    ITRVE_trt <- sqrt(ITRVE["Treatment", "Treatment"])
    
    ### Construct CI depending on SS correction
    if (ss_correct == T){
      df <- IT_summary$ngrps - 2
      ITCI <- cal_confint_t(pe = IT_coef[1], df = df, se = IT_coef[2])
      ITCI_RVE <- cal_confint_t(pe = IT_coef[1], df = df, se = ITRVE_trt)
    } else {
      ITCI <- c(IT_coef[1] - 1.96*IT_coef[2], IT_coef[1] + 1.96*IT_coef[2])
      ITCI_RVE <- c(IT_coef[1] - 1.96*ITRVE_trt, IT_coef[1] + 1.96*ITRVE_trt)
    }
  }
  
  ## Exposure time indicator model (ETI)
  ### Model fitting
  model1 <- try(glmer(Outcome ~  Period + Exposure + (1|Cluster) + (1|Cluster:Period) + (1|id_individual), 
                      family = "poisson", data = data), silent = TRUE)
  
  # Check if model1 converged
  if (check_convergence(model1)[1] == F) {
    # Fit the next model if model1 failed to converge
    model2 <- try(glmer(Outcome ~  Period + Exposure + (1|Cluster) + (1|Cluster:Period), 
                        family = "poisson", data = data), silent = TRUE)
    
    # Check if model2 also failed to converge
    if (check_convergence(model2)[1] == F) {
      model3 <- try(glmer(Outcome ~  Period + Exposure + (1|Cluster), 
                          family = "poisson", data = data), silent = TRUE)
      if(check_convergence(model3)[1] == F) {
        ETIconverge <- 0
        message("All models failed to converge.")
      } else {
        ETI <- model3
        ETIconverge <- 3
        message("model3 converged successfully.")
      }
    } else {
      ETI <- model2
      ETIconverge <- 2
      message("model2 converged successfully.")
    }
  } else {
    ETI <- model1
    ETIconverge <- 1
    message("model1 converged successfully.")
  }
  
  if(ETIconverge == 0){
    ETI_est <- ETI_vcov <- ETIRVE <- ETI_coef <- ETI_se <- ETIRVE_trt <- ETICI <- ETICI_RVE <- lrtest <- NULL 
  } else {
    lrtest <- anova(IT,ETI)
    ETI_summary <- summary(ETI)
    ETI_est <- ETI_summary$coefficients
    ETI_vcov <- vcov(ETI)
    
    ### number of periods
    index <- nlevels(data$Period)
    ### index for exposure time treatment effect in the model
    indices <- (index+1):dim(model.matrix(ETI))[2]
    
    ### Obtain estimates and SE 
    sigma.matrix <- vcov(ETI)[indices,indices]
    A <- (1/(index-1)) * matrix(rep(1,length(indices)), nrow=1)
    ETI_coef <- mean(ETI_summary$coefficients[indices,1])
    ETI_se <- as.numeric(sqrt(A %*% sigma.matrix %*% t(A)))
    
    ### Obtain RVE
    ETIRVE <- vcovCR(ETI, type= rve_type)
    sigmaRVE.matrix <- ETIRVE[indices,indices]
    ETIRVE_trt <- as.numeric(sqrt(A %*% sigmaRVE.matrix %*% t(A)))
    
    ### obtain CI depending on SS correction
    if (ss_correct == T){
      df <- ETI_summary$ngrps - 2
      ETICI <- cal_confint_t(pe = ETI_coef, df = df, se = ETI_se)
      ETICI_RVE <- cal_confint_t(pe = ETI_coef, df = df, se = ETIRVE_trt)
    } else {
      ETICI <- c(ETI_coef - 1.96*ETI_se, ETI_coef + 1.96*ETI_se)
      ETICI_RVE <- c(ETI_coef - 1.96*ETIRVE_trt, ETI_coef + 1.96* ETIRVE_trt)
    }
  }
  
  
  ## Treatment effect hetero (TEH)
  # Fit the first model
  model1 <- try(glmer(Outcome ~  Period + Treatment + (0 + Treatment|Exposure) + (1|Cluster) + (1|Cluster:Period) + (1|id_individual), 
                      family = "poisson", data = data), silent = TRUE)
  
  # Check if model1 converged
  if (check_convergence(model1)[1] == F) {
    # Fit the next model if model1 failed to converge
    model2 <- try(glmer(Outcome ~  Period + Treatment + (0 + Treatment|Exposurer) + (1|Cluster) + (1|Cluster:Period), 
                        family = "poisson", data = data), silent = TRUE)
    
    # Check if model2 also failed to converge
    if (check_convergence(model2)[1] == F) {
      model3 <- try(glmer(Outcome ~  Period + Treatment + (0 + Treatment|Exposure) + (1|Cluster), 
                          family = "poisson", data = data), silent = TRUE)
      if(check_convergence(model3)[1] == F) {
        TEHconverge <- 0
        message("All models failed to converge.")
      } else {
        TEH <- model3
        TEHconverge <- 3
        message("model3 converged successfully.")
      }
    } else {
      TEH <- model2
      TEHconverge <- 2
      message("model2 converged successfully.")
    }
  } else {
    TEH <- model1
    TEHconverge <- 1
    message("model1 converged successfully.")
  }
  
  #TEH <- glmer(Outcome ~  Period + Treatment + (Exposure|Cluster) + (1|Cluster:Period), 
  #             family = "poisson", data = data)
  if(TEHconverge == 0){
    TEH_est <- TEH_vcov <- TEHRVE <- TEH_coef <- TEHCI <- NULL 
  } else {
    TEH_summary <- summary(TEH)
    TEH_est <- TEH_summary$coefficients
    TEH_vcov <- vcov(TEH)
    TEH_coef <- TEH_summary$coefficients["Treatment",]
    
    ### Obtain RVE for estimated trt effects
    #TEHRVE <- vcovCR(TEH, type= rve_type)
    #TEHRVE_trt <- sqrt(TEHRVE["Treatment", "Treatment"])
    if (ss_correct == T){
      df <- TEH_summary$ngrps - 2
      TEHCI <- cal_confint_t(pe = TEH_coef[1], df = df, se = TEH_coef[2])
      #  TEHCI_RVE <- cal_confint_t(pe = TEH_coef[1], df = df, se = TEHRVE_trt)
    } else {
      TEHCI <- c(TEH_coef[1] - 1.96*TEH_coef[2], TEH_coef[1] + 1.96*TEH_coef[2])
      #  TEHCI_RVE <- c(TEH_coef[1] - 1.96*RVE_trt, TEH_coef[1] + 1.96*TEHRVE_trt)
    }
  }
  
  
  ## Cubic spline
  data$Exposure <- as.numeric(as.character(data$Exposure))
  data$Period <- as.numeric(as.character(data$Period))
  J <- length(unique(data$Period))
  nnode <-  ceiling((J)/2)
  ns_basis <- ns(c(0:(J-1)), knots = (1:(nnode-1))*(J-1)/(nnode))
  
  for (i in 1:(nnode)) {
    new_vector <- ns_basis[data$Exposure+1,i]  # Example: random numbers
    data[[paste0("b", i)]] <- new_vector
  }
  
  #data$b1 <- ns_basis[data$Exposure+1,1]
  #data$b2 <- ns_basis[data$Exposure+1,2]
  #data$b3 <- ns_basis[data$Exposure+1,3]
  #data$b4 <- ns_basis[data$Exposure+1,4]
  #formula <- Outcome ~ Period + b1 + b2 + b3 + b4 + (1|Cluster)
  
  b_vars <- paste0("b", 1:(nnode))
  formula1 <- as.formula(paste("Outcome ~ factor(Period) + ", paste(b_vars, collapse = " + "), "+ (1|Cluster) + (1|Cluster:Period) + (1|id_individual) "))
  formula2 <- as.formula(paste("Outcome ~ factor(Period) + ", paste(b_vars, collapse = " + "), "+ (1|Cluster) + (1|Cluster:Period)"))
  formula3 <- as.formula(paste("Outcome ~ factor(Period) + ", paste(b_vars, collapse = " + "), "+ (1|Cluster)"))
  
  model1 <- try(glmer(formula1, 
                      family = "poisson", data = data), silent = TRUE)
  
  # Check if model1 converged
  if (check_convergence(model1)[1] == F) {
    # Fit the next model if model1 failed to converge
    model2 <- try(glmer(formula2, 
                        family = "poisson", data = data), silent = TRUE)
    
    # Check if model2 also failed to converge
    if (check_convergence(model2)[1] == F) {
      model3 <- try(glmer(formula3, 
                          family = "poisson", data = data), silent = TRUE)
      if(check_convergence(model3)[1] == F) {
        NCSconverge <- 0
        message("All models failed to converge.")
      } else {
        model <- model3
        NCSconverge <- 3
        message("model3 converged successfully.")
      }
    } else {
      model <- model2
      NCSconverge <- 2
      message("model2 converged successfully.")
    }
  } else {
    model <- model1
    NCSconverge <- 1
    message("model1 converged successfully.")
  }
  
  
  if(NCSconverge == 0){
    NCS_est <- NCS_vcov <- theta_l_hat <- sigma_l_hat <- NCSRVE <- NCS_ATAE <- NCS_se <- NCSRVE_trt <- NCSCI <- NCSCI_RVE <- NULL 
  } else {
    NCS_summary <- summary(model)
    coeff_names <- names(summary(model)$coefficients[,1])
    b_hat <- as.numeric(summary(model)$coefficients[,1])
    sigma_b_hat <- vcov(model)
    indices <- c(1:length(coeff_names))[str_sub(coeff_names,1,1)=="b"]
    coeff_names <- coeff_names[indices]
    b_hat <- b_hat[indices]
    sigma_b_hat <- sigma_b_hat[indices,indices]
    sigma_b_hat <- as.matrix(sigma_b_hat)
    
    #RVE
    NCSRVE <- vcovCR(model, type= rve_type)
    NCSRVE_b_hat <- NCSRVE[indices,indices]
    B <- matrix(NA, nrow=(J-1), ncol = nnode)
    for (i in 1:(J-1)) {
      for (j in 1:(nnode)) {
        B[i,j] <- ns_basis[i+1,j]
      }
    }
    theta_l_hat <- as.numeric(B %*% b_hat)
    sigma_l_hat <- B %*% sigma_b_hat %*% t(B)
    NCSRVE <- B %*% NCSRVE_b_hat %*% t(B)
    A <- (1/(J-1)) * matrix(rep(1,length(theta_l_hat)), nrow=1)
    NCS_ATAE = (A %*% theta_l_hat)[1,1]
    NCS_se = sqrt(A %*% sigma_l_hat %*% t(A))[1,1]
    NCSRVE_trt = sqrt(A %*% NCSRVE %*% t(A))[1,1]
    
    
    ### obtain CI depending on SS correction
    if (ss_correct == T){
      df <- NCS_summary$ngrps - 2
      NCSCI <- cal_confint_t(pe = NCS_ATAE, df = df, se = NCS_se)
      NCSCI_RVE <- cal_confint_t(pe = NCS_ATAE, df = df, se = NCSRVE_trt)
    } else {
      NCSCI <- c(NCS_ATAE - 1.96*NCS_se, NCS_ATAE + 1.96*NCS_se)
      NCSCI_RVE <- c(NCS_ATAE - 1.96*NCSRVE_trt, NCS_ATAE + 1.96* NCSRVE_trt)
    }
  }
  
  ## return results
  results[["IT"]]$est <- IT_est
  results[["IT"]]$vcov <- IT_vcov
  results[["IT"]]$rve <- ITRVE
  results[["IT"]]$summary <- data.frame(TATE = IT_coef[1], SD = IT_coef[2], 
                                        RVE = ITRVE_trt, CI_L = ITCI[1], CI_U = ITCI[2], 
                                        CI_RVE_L = ITCI_RVE[1], CI_RVE_U = ITCI_RVE[2])
  results[["ETI"]]$est <- ETI_est
  results[["ETI"]]$vcov <- ETI_vcov
  results[["ETI"]]$rve <- ETIRVE
  results[["ETI"]]$summary <- data.frame(TATE = ETI_coef, SD = ETI_se, 
                                         RVE = ETIRVE_trt, CI_L = ETICI[1], CI_U = ETICI[2], 
                                         CI_RVE_L = ETICI_RVE[1], CI_RVE_U = ETICI_RVE[2])
  results[["ETI"]]$lrtest <- lrtest
  
  results[["NCS"]]$est <- theta_l_hat
  results[["NCS"]]$vcov <- sigma_l_hat
  results[["NCS"]]$rve <- NCSRVE
  results[["NCS"]]$summary <- data.frame(TATE = NCS_ATAE, SD = SE_NCS_ATAE, 
                                         RVE = NCSRVE_trt, CI_L = NCSCI[1], CI_U = NCSCI[2], 
                                         CI_RVE_L = NCSCI_RVE[1], CI_RVE_U = NCSCI_RVE[2])
  
  results[["TEH"]]$est <- TEH_est
  #results[["TEH"]]$vcov <- TEH_vcov
  #results[["TEH"]]$rve <- TEHRVE
  results[["TEH"]]$summary <- data.frame(TATE = TEH_coef[1], SD = TEH_coef[2], 
                                         CI_L = TEHCI[1], CI_U = TEHCI[2])
  
  results[["convergence"]] <- data.frame(IT=ITconverge, ETI=ETIconverge, TEH=TEHconverge, NCS=NCSconverge)
  return(results)
}
