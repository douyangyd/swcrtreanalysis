# Identify different types of outcome variables
outcome_columns <- grep("out", names(data), value = TRUE)
outcome_columns <- outcome_columns[grep("con|bin|poiss", outcome_columns)]
off_columns <- grep("off", names(data), value = TRUE)
ind_columns <- grep("id_individual", names(data), value = TRUE)
out_poiss_vars <- grep("out_poiss", names(data), value = TRUE)  # Poisson outcomes
out_other_vars <- setdiff(grep("out_", names(data), value = TRUE), out_poiss_vars)  # Other outcomes
out_other_vars <- out_other_vars[grep("con|bin", out_other_vars)]

process_data <- function(data){
  # Create the list of datasets
  dataset_list <- list()
  
  # Process Poisson outcome variables
  for (outcome in out_poiss_vars) {
    suffix <- sub("out_poiss_", "", outcome)
    offset_var <- paste0("off_poiss_", suffix)
    #dataset_list$outcome <- dataset_list[[outcome]]
    
    if (offset_var %in% names(data)) {
      datacopy <- data
      datacopy$outcome <- datacopy[[outcome]]
      datacopy$offset <- datacopy[[offset_var]]
      datacopy <- datacopy[, !names(data) %in% c(outcome_columns,off_columns)]
      dataset_list[[outcome]] <- datacopy
    } else {
      datacopy <- data
      datacopy$outcome <- datacopy[[outcome]]
      datacopy <- datacopy[, !names(data) %in% c(outcome_columns,off_columns)]
      dataset_list[[outcome]] <- datacopy
    }
  }
  
  # Process other outcome variables (binary, continuous)
  for (outcome in out_other_vars) {
    datacopy <- data
    datacopy$outcome <- datacopy[[outcome]]
    datacopy <- datacopy[, !names(data) %in% c(outcome_columns,off_columns)]
    dataset_list[[outcome]] <- datacopy
  }
  
  # Print structure of the list
  #str(dataset_list)
  
  data.format <- function(dataset_list){
    data <- list()
    for (i in 1:length(dataset_list)){
      data[[i]] <- dplyr::rename(dataset_list[[i]], Outcome = outcome, Cluster = id_cluster, Treatment = trt, Period = time, Exposure = time_on_trt)
      data[[i]][,c("Cluster", "Period")] <- lapply(data[[i]][,c("Cluster", "Period")], factor)
    }
    return(data)
  }
  
  ## Return the list of the dataset
  data <- data.format(dataset_list)
  
  return(data)
}

data <- process_data(data)


family <- sapply(outcome_columns, function(x) {
  if (grepl("poiss", x)) return("poisson")
  if (grepl("bin", x)) return("binomial")
  if (grepl("con", x)) return("gaussian")
})
