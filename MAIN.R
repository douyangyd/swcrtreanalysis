
# Set configuration
cfg <- list(
  quiet = T,
  pkgs = c("dplyr", "msm", "lmerTest", "splines", "performance", "stringr",
           "sandwich", "clubSandwich", "WoodburyMatrix", "steppedwedge",
           "htmltools", "webshot", "DT", "here", "ggplot2"),
  #path = "G:/Shared drives/Stepped Wedge Data Files/Trials/",
  path = "~/Library/CloudStorage/GoogleDrive-derek6561@gmail.com/Shared drives/Stepped Wedge Data Files/Trials/",
  ind = c(10) # This corresponds to the datasets in dir(cfg$path)
  # ind = c(1:74)
)
# c(9,15:17,21,26:28)
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
#here::i_am("swcrtreanalysis/MAIN.R")
#setwd(paste0(here::here(), "/swcrtreanalysis"))
#setwd("G:/Shared drives/Stepped Wedge Data Files/Software/Data Analysis")
setwd("~/Library/CloudStorage/GoogleDrive-derek6561@gmail.com/Shared drives/Stepped Wedge Data Files/Software/Data Analysis")
files <- c("data_analysis.R", "helpers.R", "process results.R")
for (file in files) { do.call("source", list(file)) }
#source("G:/Shared drives/Stepped Wedge Data Files/Software/Robust Variance/vcovCRglmerMod.R") # Temporary
source("~/Library/CloudStorage/GoogleDrive-derek6561@gmail.com/Shared drives/Stepped Wedge Data Files/Software/Robust Variance/vcovCRglmerMod.R")
# Loop through datasets and run analysis
for (i in cfg$ind) {
  
  data <- read.csv(paste0(
    cfg$path, dir(cfg$path)[i], "/Standardized data/trial_data.csv"
  ))
  dataname <- print(dir(cfg$path)[i])
  head(data)
  source("~/Library/CloudStorage/GoogleDrive-derek6561@gmail.com/Shared drives/Stepped Wedge Data Files/software/Data Analysis/Data_processing.R")
  # source("~/OneDrive/SickKids/Method Project/TimedependentRVE/swcrtreanalysis2/data_processing_test.R") 
  # number of outcomes in this dataset
  length(outcome_columns)
  # type of outcomes for each outcome
  family
  # cs or co design
  design <- ifelse(length(ind_columns) == 0, "cs", "co")
  # create a Result dir under the trial folder (if not exist)
  drive_folder <- paste0(
    cfg$path, dir(cfg$path)[i], "/Results"
  )  # Windows
  # drive_folder <- "~/Google Drive/Results/"  # Mac/Linux
  
  # Ensure the folder exists
  if (!dir.exists(drive_folder)) dir.create(drive_folder, recursive = TRUE)
  
  
  for (j in 1:length(data)) {
    result <- fit(
      data = data[[j]],
      family = family[j], 
      design = design,
      rse_type = "MD",
      offset = if(is.null(data[[j]]$offset)==T){NULL}else{data[[j]]$offset},
      aggregate = F,
      ss_correct = T
    )
    table <- summary_table(result)
    plot <- plot_design(swdat)$design_plot
    file_name1 <- paste0(dataname, "_model", "_", outcome_columns[j], ".rds")
    file_name2 <- paste0(dataname, "_table", "_", outcome_columns[j], ".pdf")
    file_name3 <- paste0(dataname, "_plot", "_", outcome_columns[j], ".pdf")
    result$abbr <- paste(c("itm = immediate treatment model", 
                               "etim = exposure-time indicator model", 
                               "tehm = treatment-exposure heterogeneity model", 
                               "ncsm = natural cubic spline model",
                               "reest = random effect estimates", 
                               "est = estimated average treatment effect",
                               "se = standard error", 
                               "rse = robust standard error", 
                               "ci = confidence interval",
                               "lte = long-term treatment effect"))
    #browser()
    result$design <- design
    setwd(drive_folder)
    saveRDS(result, file_name1)
    html_file <- tempfile(fileext = ".html")
    saveWidget(table, html_file)
    webshot(html_file, file_name2)
    ggsave(file_name3, plot, width=8, height=11, units="in", dpi=400)
  }
}

