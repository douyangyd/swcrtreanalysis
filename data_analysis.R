fit <- function(
    data,                  # standardized dataset
    family,                #
    rse_type,              # rse type: MD was recommended 
    ss_correct,            # Whether small sample correction will be applied: T/F
    offset = FALSE,
    design                 # cs = cross-section; co = cohort 
    #re
){
  
  # `steppdwedge` setup
  data2 <- data
  if (is.factor(data2$Period)) { data2$Period <- as.numeric(data2$Period) }
  if (design=="cs") {
    swdat <- load_data(time="Period", cluster_id="Cluster", treatment="Treatment",
                       outcome="Outcome", data=data2)
  } else if (design=="co") {
    swdat <- load_data(time="Period", cluster_id="Cluster", individual_id="id_individual",
                       treatment="Treatment", outcome="Outcome", data=data2)
  }
  
  results <- list()
  
  if (offset == F) {f_offset <- ""} else {
    data$off <- as.numeric(as.character(data[[offset]]))
    f_offset <- " + offset(log(off))"
  }
  
  ################################################.
  #####       Immediate Treatment (IT)       #####
  ################################################
  
  if (design == "cs") {
    model1 <- try(analyze(dat=swdat, family=family, re=c("clust","time")), silent=T)
    model2 <- try(analyze(dat=swdat, family=family, re=c("clust")), silent=T)
    model3 <- NULL
  } 
  if (design == "co"){
    model1 <- try(analyze(dat=swdat, family=family, re=c("clust","time","ind")), silent=T)
    model2 <- try(analyze(dat=swdat, family=family, re=c("clust","ind")), silent=T)
    model3 <- try(analyze(dat=swdat, family=family, re=c("clust")), silent=T)
  }
    if(is.null(model1) == T){
      itm1_result <- NULL
    } else {
      itm1_result <- get_coef(model1, rse_type, ss_correct)
    }
    
    if(is.null(model2) == T){
      itm2_result <- NULL
    } else {
      itm2_result <- get_coef(model2, rse_type, ss_correct)
    }
    
    if(is.null(model3) == T){
      itm3_result <- NULL
    } else {
      itm3_result <- get_coef(model3, rse_type, ss_correct)
    }
  
  
  ################################################.
  #####    Exposure Time Indicator (ETI)     #####
  ################################################.  

  if (design == "cs"){
    model1 <- try(analyze(dat=swdat, estimand_type="TATE", exp_time="ETI", family=family, re=c("clust","time")), silent=T)
    model2 <- try(analyze(dat=swdat, estimand_type="TATE", exp_time="ETI", family=family, re=c("clust")), silent=T)
    model3 <- NULL
  }
  if (design == "co"){
    model1 <- try(analyze(dat=swdat, estimand_type="TATE", exp_time="ETI", family=family, re=c("clust","time","ind")), silent=T)
    model2 <- try(analyze(dat=swdat, estimand_type="TATE", exp_time="ETI", family=family, re=c("clust","ind")), silent=T)
    model3 <- try(analyze(dat=swdat, estimand_type="TATE", exp_time="ETI", family=family, re=c("clust")), silent=T)
  }
    if(is.null(model1) == T){
      etim1_result <- NULL
    } else {
      etim1_result <- get_eticoef(model1, rse_type, ss_correct)
    }
    
    if(is.null(model2) == T){
      etim2_result <- NULL
    } else {
      etim2_result <- get_eticoef(model2, rse_type, ss_correct)
    }
    
    if(is.null(model3) == T){
      etim3_result <- NULL
    } else {
      etim3_result <- get_eticoef(model3, rse_type, ss_correct)
    }
  
  
  
  #####################################################.
  #####   Treatment Effect Heterogeneity (TEH)   #####
  #####################################################.

  if (design == "cs"){
    model1 <- try(analyze(dat=swdat, estimand_type="TATE", exp_time="TEH", family=family, re=c("clust","time")), silent=T)
    model2 <- try(analyze(dat=swdat, estimand_type="TATE", exp_time="TEH", family=family, re=c("clust")), silent=T)
  }
  if (design == "co"){
    model1 <- try(analyze(dat=swdat, estimand_type="TATE", exp_time="TEH", family=family, re=c("clust","time","ind")), silent=T)
    model2 <- try(analyze(dat=swdat, estimand_type="TATE", exp_time="TEH", family=family, re=c("clust","ind")), silent=T)
    model3 <- try(analyze(dat=swdat, estimand_type="TATE", exp_time="TEH", family=family, re=c("clust")), silent=T)
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
  
  
  ############################################.
  #####    Natural Cubic Spline (NCS)    #####
  ############################################.
  
  
  ## Cubic spline
  data$Period <- as.numeric(as.character(data$Period))
  J <- length(unique(data$Exposure))
  nnode <-  ceiling((J)/2)
  n_knots <- nnode + 1
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
    model1 <- try(analyze(dat=swdat, estimand_type="TATE", exp_time="NCS", family=family,
                          re=c("clust","time"), n_knots_exp=n_knots), silent=T)
    model2 <- try(analyze(dat=swdat, estimand_type="TATE", exp_time="NCS", family=family,
                          re=c("clust")), n_knots_exp=n_knots, silent=T)
    model3 <- NULL
    # if(family == "gaussian"){
    #   model1 <- try(lme4::lmer(formula3, data = data), silent = TRUE)
    #   model2 <- try(lme4::lmer(formula4, data = data), silent = TRUE)
    # } else {
    #   model1 <- try(glmer(formula3, family = family, data = data), silent = TRUE)
    #   model2 <- try(glmer(formula4, family = family, data = data), silent = TRUE)
    # }
  }
  if (design == "co"){
    model1 <- try(analyze(dat=swdat, estimand_type="TATE", exp_time="NCS", family=family,
                          re=c("clust","time","ind"), n_knots_exp=n_knots), silent=T)
    model2 <- try(analyze(dat=swdat, estimand_type="TATE", exp_time="NCS", family=family,
                          re=c("clust","ind"), n_knots_exp=n_knots), silent=T)
    model3 <- try(analyze(dat=swdat, estimand_type="TATE", exp_time="NCS", family=family,
                          re=c("clust"), n_knots_exp=n_knots), silent=T)
  }
    # if(family == "gaussian"){
    #   model1 <- try(lme4::lmer(formula1, data = data), silent = TRUE)
    #   model2 <- try(lme4::lmer(formula2, data = data), silent = TRUE)
    #   model3 <- try(lme4::lmer(formula4, data = data), silent = TRUE)
    # } else {
    #   model1 <- try(glmer(formula1, family = family, data = data), silent = TRUE)
    #   model2 <- try(glmer(formula2, family = family, data = data), silent = TRUE)
    #   model3 <- try(glmer(formula4, family = family, data = data), silent = TRUE)
    # }
    if(is.null(model1) == T){
      ncsm1_result <- NULL
    } else {
      ncsm1_result <- get_ncscoef(model1, data, rse_type, ss_correct, ns_basis, J, nnode)
    }
    
    if(is.null(model2) == T){
      ncsm2_result <- NULL
    } else {
      ncsm2_result <- get_ncscoef(model2, data, rse_type, ss_correct, ns_basis, J, nnode)
    }
    
    if(is.null(model3) == T){
      ncsm3_result <- NULL
    } else {
      ncsm3_result <- get_ncscoef(model3, data, rse_type, ss_correct, ns_basis, J, nnode)
    }
  
  if (design == "cs"){
    results <- list(family = family, itm1=itm1_result, itm2=itm2_result, etim1 = etim1_result, etim2 = etim2_result, 
                    tehm1 = tehm1_result, tehm2 = tehm2_result, ncsm1 = ncsm1_result, ncsm2 = ncsm2_result)
  }
  if (design == "co"){
    results <- list(family = family, itm1=itm1_result, itm2=itm2_result, itm3=itm3_result, 
                    etim1 = etim1_result, etim2 = etim2_result, etim3 = etim3_result, 
                    tehm1 = tehm1_result, tehm2 = tehm2_result, tehm3 = tehm3_result,
                    ncsm1 = ncsm1_result, ncsm2 = ncsm2_result, ncsm3 = ncsm3_result)
  }
  return(results)
}

#testx <- fit(data = data, family = "gaussian", rse_type = "CR0", ss_correct=T, design = "co")
