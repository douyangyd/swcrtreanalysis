# Set configuration
cfg <- list(
  quiet = T,
  pkgs = c("dplyr", "msm", "lmerTest", "splines", "performance", "stringr",
           "sandwich", "clubSandwich", "WoodburyMatrix", "steppedwedge"),
  path = "G:/Shared drives/Stepped Wedge Data Files/Trials/",
  wd = "swcrtreanalysis",
  ind = 23 # This corresponds to the datasets in dir(cfg$path)
  # ind = c(1:74)
)

# Load packages
# To install `steppedwedge`, run the following:
#     devtools::install_github(repo="Avi-Kenny/steppedwedge", dependencies=T)
for (pkg in c(cfg$pkgs)) {
  if (cfg$quiet) {
    suppressMessages({ do.call("library", list(pkg)) })
  } else {
    do.call("library", list(pkg))
  }
}

# Source necessary files
setwd(cfg$wd)
files <- c("Data_processing.R", "data_analysis.R", "helpers.R")
for (file in files) { do.call("source", list(file)) }
source("G:/Shared drives/Stepped Wedge Data Files/Software/Robust Variance/vcovCRglmerMod.R") # Temporary

# Loop through datasets and run analysis
for (i in cfg$ind) {
  
  data <- read.csv(paste0(
    cfg$path, dir(cfg$path)[i], "/Standardized data/trial_data.csv"
  ))
  head(data)
  data <- process_data(data)
  
  test1 <- fit(
    data = data[[1]],
    family = "gaussian", 
    design = "cs",
    rve_type = "classic",
    ss_correct = T
  )
  
}
