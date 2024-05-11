#!/usr/bin/env Rscript

# Load necessary libraries
suppressPackageStartupMessages({
  library(copula)
  library(stats)
  library(data.table)
  library(optparse)
})

# Set up command-line options using optparse
option_list <- list(
  make_option(c("-i", "--input_file_1"), type = "character", default = "", help = "Path to the first TSV file containing gene expression data."),
  make_option(c("-j", "--input_file_2"), type = "character", default = "", help = "Path to the second TSV file containing gene expression data."),
  make_option(c("-n", "--iterations"), type = "integer", default = 10, help = "Number of iterations for performance measurement."),
  make_option(c("-o", "--output_path"), type = "character", default = "", help = "Output directory path where the performance results will be saved.")
)

# Parse command-line arguments
parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

# Validate inputs
if (args$input_file_1 == "" || args$input_file_2 == "" || args$output_path == "") {
  print_help(parser)
  stop("Missing arguments, please provide all required inputs.", call. = FALSE)
}

# Function to compute the distance matrix
distance_matrix <- function(data1, data2) {
  dist_mat <- matrix(0, nrow = 2000, ncol = 2000)
  data_col1 <- matrix(0, nrow = ncol(data1) - 1, ncol = 2)
  data_col2 <- matrix(0, nrow = ncol(data2) - 1, ncol = 2)

  for (i in 1:1999) {
    for (j in (i + 1):2000) {
      data_col1[, 1] = t(data1[i, -1])
      data_col1[, 2] = t(data1[j, -1])
      data_col2[, 1] = t(data2[i, -1])
      data_col2[, 2] = t(data2[j, -1])

      u1 <- pobs(data_col1)
      u2 <- pobs(data_col2)

      ec1 <- C.n(u1, data_col1, smoothing = 'none', ties.method = 'average')
      ec2 <- C.n(u2, data_col2, smoothing = 'none', ties.method = 'average')

      p = ks.test(ec1, ec2)
      dist_mat[i, j] = p$statistic
    }
  }
  return(dist_mat)
}

# Load data
data1 <- fread(args$input_file_1, sep = "\t")
data2 <- fread(args$input_file_2, sep = "\t")

# Measure execution time over multiple runs
execution_times <- numeric(args$iterations)

# Initialize progress bar
pb <- txtProgressBar(min = 0, max = args$iterations, style = 3)
on.exit(close(pb))

for (i in 1:args$iterations) {
  cat(sprintf("Executing iteration %d...\n", i))
  start_time <- Sys.time()
  distance <- distance_matrix(data1, data2)
  execution_times[i] <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  setTxtProgressBar(pb, i)  # Update progress bar
}

# Save execution times to a CSV file
execution_data <- data.frame(Execution = 1:args$iterations, Time = execution_times)
output_csv <- file.path(args$output_path, "r_performance.csv")
fwrite(execution_data, output_csv)

cat(sprintf("Saved execution times to %s\n", output_csv))
