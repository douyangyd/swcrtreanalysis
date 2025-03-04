summary_table <- function(results){
  # Define the rows (metrics)
  rows <- c(
    "est",
    "se",
    "ci",
    "rse",
    "rseci",
    "lte",
    "lte_se",
    "lteci",
    "lte_rse",
    "lterseci",
    "converged",
    "lme4_converged",
    "messages",
    "reest"
  )
  
  if (length(results) == 9){
    table_data <- data.frame(
      Metric = rows,
      itm1 = NA,
      itm2 = NA,
      etim1 = NA,
      etim2 = NA,
      tehm1 = NA,
      tehm2 = NA,
      ncsm1 = NA,
      ncsm2 = NA,
      stringsAsFactors = FALSE
    )
  } else {
    table_data <- data.frame(
      Metric = rows,
      itm1 = NA,
      itm2 = NA,
      itm3 = NA,
      etim1 = NA,
      etim2 = NA,
      etim3 = NA,
      tehm1 = NA,
      tehm2 = NA,
      tehm3 = NA,
      ncsm1 = NA,
      ncsm2 = NA,
      ncsm3 = NA,
      stringsAsFactors = FALSE
    )
  }
  
  
  format_value <- function(value) {
    if (is.numeric(value) && length(value) == 1) {
      # Single numeric value
      return(format(value, digits = 3))
    } else if (is.numeric(value) && length(value) == 2) {
      # Confidence interval or range
      return(paste0("[", format(value[1], digits = 3), ", ", format(value[2], digits = 3), "]"))
    } else {
      # Assume it's already a string
      return(as.character(value))
    }
  }
  
  # Populate the table by extracting and formatting values from the lists
  for (col in names(results)[-1]) {
    table_data[[col]] <- sapply(rows, function(row) format_value(results[[col]][[row]]))
  }
  table_data <- table_data[,-1]
  
  if (length(results) == 9){
  colnames(table_data) <- c("Immediate Time (IT) Model 1", 
                            "Immediate Time (IT) Model 2",
                            "Exposure Time Indicator (ETI) Model 1", 
                            "Exposure Time Indicator (ETI) Model 2",
                            "Treatment Effect Heterogeneity (TEH) Model 1", 
                            "Treatment Effect Heterogeneity (TEH) Model 2",
                            "Natural Cubic Spline (NCS) Model 1",
                            "Natural Cubic Spline (NCS) Model 2")
  } else {
    colnames(table_data) <- c("Immediate Time (IT) Model 1", 
                              "Immediate Time (IT) Model 2",
                              "Immediate Time (IT) Model 3",
                              "Exposure Time Indicator (ETI) Model 1", 
                              "Exposure Time Indicator (ETI) Model 2",
                              "Exposure Time Indicator (ETI) Model 3",
                              "Treatment Effect Heterogeneity (TEH) Model 1", 
                              "Treatment Effect Heterogeneity (TEH) Model 2",
                              "Treatment Effect Heterogeneity (TEH) Model 3",
                              "Natural Cubic Spline (NCS) Model 1",
                              "Natural Cubic Spline (NCS) Model 2",
                              "Natural Cubic Spline (NCS) Model 3")
  }
  
  rownames(table_data) <- c("Estimate",
                            "Standard Error (SE)",
                            "95% Confidence Interval",
                            "Robust Standard Error",
                            "95% CI with Robust SE",
                            "Long-term treatment effect (LTE)",
                            "LTE SE",
                            "LTE 95% Confidence Interval",
                            "LTE Robust SE",
                            "LTE 95% CI with Robust SE",
                            "Convergence",
                            "lme4 Convergence",
                            "Model Message",
                            "Random Effect Estimates")
  
  caption_text <- if (ncol(table_data) == 8) {
    "Model 1 Random effect: Cluster and Cluster by time. Model 2 Random effect: Cluster only"
  } else {
    "Model 1 Random effect: Cluster, Cluster by time and indiviudal. Model 2 Random effect: Cluster, Cluster by time. Model 3 Random effect: Cluster only"
  }
  # Render as an interactive DT table
  
  
  datatable(table_data, 
            rownames = TRUE, 
            options = list(
              pageLength = nrow(table_data),
              dom = 't', # Only show the table without search or pagination
              autoWidth = TRUE
            ),
            caption = htmltools::tags$caption(
              style = 'caption-side: bottom; text-align: left;',
              caption_text
            ), escape = FALSE)
}
