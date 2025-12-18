library(readxl)
library(ComplexHeatmap)
library(circlize)
install.packages("viridis")  # only run once
library(viridis)
library(ggplot2)
library(tidyverse)

  # Load the data
  data <- read_excel("DeltaFdivideF2_batch2 combined with all batch 1.xlsx", sheet = "Sheet1")
  
  # Convert columns 2 to 81 into matrix (cell traces)
  df_matrix <- as.matrix(data[, 1:80])
  
  # Transpose so that: rows = cells, columns = timepoints
  df_matrix <- t(df_matrix)
  
  # Set cell names and time as column names
  rownames(df_matrix) <- colnames(data)[1:80]
  colnames(df_matrix) <- round(data[[1]], 1)  # use rounded time values for columns
  
  # Extract only finite values to avoid NA/Inf in color scaling
  finite_vals <- df_matrix[is.finite(df_matrix)]
  

  # focus on central 99.98% of the values for better contrast
  q <- quantile(finite_vals, probs = c(0.001, 0.999))
  col_fun <- colorRamp2(
    seq(q[1], q[2], length.out = 100),
    viridis(100)
  )
  
  top_third_mean <- function(x) {
    x <- x[is.finite(x)]  # remove NA/Inf
    n <- length(x)
    top_n <- ceiling(n / 3)
    mean(sort(x, decreasing = TRUE)[1:top_n])
  }
  
  # Apply to each row (cell)
  top_means <- apply(df_matrix, 1, top_third_mean)
  
  # Sort based on those means
  row_order <- order(top_means, decreasing = TRUE)
  df_matrix_ordered <- df_matrix[row_order, ]
  
  # Plot heatmap
  Heatmap(df_matrix_ordered,
          name = "ΔF/F",
          col = col_fun,
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_row_names = TRUE,
          show_column_names = FALSE,
          row_title = "ΔF/F",
          column_title = "Time (s)",
          column_title_side = "bottom")
  
  
 
