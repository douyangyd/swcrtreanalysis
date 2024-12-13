# Data Processing
process_data <- function(data) {

### Identify outcome columns (assuming outcome columns are named 'outcome1', 'outcome2', etc.)
outcome_columns <- grep("out", names(data), value = TRUE)

### number of binary outcome
binary_columns <-grep("bin", names(data), value=TRUE)
### number of continuous outcome
gaussian_columns <-grep("con", names(data), value=TRUE)
### number of count outcome
poisson_columns <-grep("poiss", names(data), value=TRUE)


### Identify offset columns (only for count outcomes, for other outcomes, this is always 0)
offset_columns <- grep("off", names(data), value = TRUE)

### Identify indiivudal ID columns  
ind_columns <- grep("id_individual", names(data), value = TRUE)

### Create a list to hold copies of the dataset
dataset_copies <- list()

### Loop through each outcome column to create a copy of the dataset
### The outcome of this function is a list of datasets (length of the list equals to the number of columns)
### Each dataset contains one outcome with all other variables that we need for model fitting
### I choose this way to avoid the dynamic changes of outcome variable names when fit the model 
for (outcome in outcome_columns) {
  dataset_copy <- data
  dataset_copy$outcome <- dataset_copy[[outcome]]
  dataset_copy <- dataset_copy[, !names(dataset_copy) %in% outcome_columns]  # Drop original outcome columns
  dataset_copies[[outcome]] <- dataset_copy
}

### Display the copies
dataset_copies

### Data formulation
data.format <- function(dataset_copies){
  data <- list()
  for (i in 1:length(dataset_copies)){
    data[[i]] <- dplyr::rename(dataset_copies[[i]], Outcome = outcome, Cluster = id_cluster, Treatment = trt, Period = time, Exposure = time_on_trt)
    data[[i]][,c("Cluster", "Treatment", "Period", "Exposure")] <- lapply(data[[i]][,c("Cluster", "Treatment", "Period", "Exposure")], factor)
  }
  return(data)
}

## Return the list of the dataset
data <- data.format(dataset_copies)

return(data)

}
