fit <- function(
    data,                  # standardized dataset
    family,                #
    rse_type,              # rse type: MD was recommended 
    ss_correct,            # Whether small sample correction will be applied: T/F
    offset = NULL,
    aggregate = F,
    design                 # cs = cross-section; co = cohort 
    #re
){
  
  # `steppdwedge` setup
  data2 <- na.omit(data)
  if (is.factor(data2$Period)) { data2$Period <- as.numeric(data2$Period) }
  if (aggregate == T) {
    swdat <<- load_data(
      time="Period", cluster_id="Cluster", treatment="Treatment",
      outcome = c("out_numerator", "out_denominator"), exposure_time="Exposure", data=data2
    )
  } else {
    if (design=="cs") {
      swdat <<- load_data(
        time="Period", cluster_id="Cluster", treatment="Treatment",
        outcome="Outcome", exposure_time="Exposure", data=data2
      )
    } else if (design=="co") {
      swdat <<- load_data(
        time="Period", cluster_id="Cluster", individual_id="id_individual",
        treatment="Treatment", outcome="Outcome", exposure_time="Exposure",
        data=data2
      )
    }
  }

  results <- list()
  
  if (is.null(offset) == T) {swdat$offset <- NULL} else {
    swdat$offset <- log(as.numeric(as.character(data2$offset)))
    swdat <- swdat[swdat$offset>=0,]
  }
  
  
  
  ################################################.
  #####       Immediate Treatment (IT)       #####
  ################################################
  
  if (design == "cs") {
    model1 <- try(analyze(dat=swdat, family=family, offset = offset, re=c("clust","time")), silent=T)
    model2 <- try(analyze(dat=swdat, family=family, offset = offset, re=c("clust")), silent=T)
    model3 <- NULL
  } 
  if (design == "co"){
    model1 <- try(analyze(dat=swdat, family=family, offset = offset, re=c("clust","time","ind")), silent=T)
    model2 <- try(analyze(dat=swdat, family=family, offset = offset, re=c("clust","ind")), silent=T)
    model3 <- try(analyze(dat=swdat, family=family, offset = offset, re=c("clust")), silent=T)
  }
  
  if(any(class(model1) %in% "try-error") ){
    itm1_result <- NULL
  } else {
    itm1_result <- get_coef(model1, rse_type, ss_correct)
  }
  
  if(any(class(model2) %in% "try-error") ){
    itm2_result <- NULL
  } else {
    itm2_result <- get_coef(model2, rse_type, ss_correct)
  }
  
  if(is.null(model3) | any(class(model3) %in% "try-error") ){
    itm3_result <- NULL
  } else {
    itm3_result <- get_coef(model3, rse_type, ss_correct)
  }
  
  
  ################################################.
  #####    Exposure Time Indicator (ETI)     #####
  ################################################.  
  
  if (design == "cs"){
    model1 <- try(analyze(dat=swdat, estimand_type="TATE", exp_time="ETI", family=family, offset = offset, re=c("clust","time")), silent=T)
    model2 <- try(analyze(dat=swdat, estimand_type="TATE", exp_time="ETI", family=family, offset = offset, re=c("clust")), silent=T)
    model3 <- NULL
  }
  if (design == "co"){
    model1 <- try(analyze(dat=swdat, estimand_type="TATE", exp_time="ETI", family=family, offset = offset, re=c("clust","time","ind")), silent=T)
    model2 <- try(analyze(dat=swdat, estimand_type="TATE", exp_time="ETI", family=family, offset = offset, re=c("clust","ind")), silent=T)
    model3 <- try(analyze(dat=swdat, estimand_type="TATE", exp_time="ETI", family=family, offset = offset, re=c("clust")), silent=T)
  }
  
  if(any(class(model1) %in% "try-error") ){
    etim1_result <- NULL
  } else {
    etim1_result <- get_eticoef(model1, rse_type, ss_correct)
  }
  
  if(any(class(model2) %in% "try-error") ){
    etim2_result <- NULL
  } else {
    etim2_result <- get_eticoef(model2, rse_type, ss_correct)
  }
  
  if(is.null(model3) | any(class(model3) %in% "try-error") ){
    etim3_result <- NULL
  } else {
    etim3_result <- get_eticoef(model3, rse_type, ss_correct)
  }
  
  
  
  #####################################################.
  #####   Treatment Effect Heterogeneity (TEH)   #####
  #####################################################.
  
  if (design == "cs"){
    model1 <- try(analyze(dat=swdat, estimand_type="TATE", exp_time="TEH", family=family, offset = offset, re=c("clust","time")), silent=T)
    model2 <- try(analyze(dat=swdat, estimand_type="TATE", exp_time="TEH", family=family, offset = offset, re=c("clust")), silent=T)
  }
  if (design == "co"){
    model1 <- try(analyze(dat=swdat, estimand_type="TATE", exp_time="TEH", family=family, offset = offset, re=c("clust","time","ind")), silent=T)
    model2 <- try(analyze(dat=swdat, estimand_type="TATE", exp_time="TEH", family=family, offset = offset, re=c("clust","ind")), silent=T)
    model3 <- try(analyze(dat=swdat, estimand_type="TATE", exp_time="TEH", family=family, offset = offset, re=c("clust")), silent=T)
  }
  
  if(any(class(model1) %in% "try-error") ){
    tehm1_result <- NULL
  } else {
    tehm1_result <- get_tehcoef(model1, ss_correct)
  }
  
  if(any(class(model2) %in% "try-error") ){
    tehm2_result <- NULL
  } else {
    tehm2_result <- get_tehcoef(model2, ss_correct)
  }
  
  if(is.null(model3) | any(class(model3) %in% "try-error") ){
    tehm3_result <- NULL
  } else {
    tehm3_result <- get_tehcoef(model3, ss_correct)
  }
  
  
  ############################################.
  #####    Natural Cubic Spline (NCS)    #####
  ############################################.
  
  
  ## Cubic spline
  if (max(swdat$exposure_time)==2){
    ncsm1_result <- etim1_result
    ncsm2_result <- etim2_result
    ncsm3_result <- etim3_result
  } else {
    n_knots <- ceiling(length(unique(swdat$exposure_time))/2)
    if (design == "cs"){
      model1 <- try(analyze(dat=swdat, estimand_type="TATE", exp_time="NCS", family=family,
                            offset = offset, re=c("clust","time"), n_knots_exp=n_knots), silent=T)
      model2 <- try(analyze(dat=swdat, estimand_type="TATE", exp_time="NCS", family=family,
                            offset = offset, re=c("clust"), n_knots_exp=n_knots), silent=T)
      model3 <- NULL
    }
    if (design == "co"){
      model1 <- try(analyze(dat=swdat, estimand_type="TATE", exp_time="NCS", family=family,
                            offset = offset, re=c("clust","time","ind"), n_knots_exp=n_knots), silent=T)
      model2 <- try(analyze(dat=swdat, estimand_type="TATE", exp_time="NCS", family=family,
                            offset = offset, re=c("clust","ind"), n_knots_exp=n_knots), silent=T)
      model3 <- try(analyze(dat=swdat, estimand_type="TATE", exp_time="NCS", family=family,
                            offset = offset, re=c("clust"), n_knots_exp=n_knots), silent=T)
    }
    
    if(any(class(model1) %in% "try-error") ){
      ncsm1_result <- NULL
    } else {
      ncsm1_result <- get_ncscoef(model1, swdat, n_knots, rse_type, ss_correct)
    }
    
    if(any(class(model2) %in% "try-error") ){
      ncsm2_result <- NULL
    } else {
      ncsm2_result <- get_ncscoef(model2, swdat, n_knots, rse_type, ss_correct)
    }
    
    if(is.null(model3) | any(class(model3) %in% "try-error") ){
      ncsm3_result <- NULL
    } else {
      ncsm3_result <- get_ncscoef(model3, swdat, n_knots, rse_type, ss_correct)
    }
  }
  
  
  
  ##########################################.
  #####    COnstruct results object    #####
  ##########################################.
  
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

