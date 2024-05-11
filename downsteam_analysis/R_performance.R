# Install necessary packages
# install.packages("copula")
# install.packages("data.table")  # For fast data reading/writing

library(copula)
library(stats)
library(data.table)

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
data1 <- fread("/Users/pranto/Desktop/BIONet/BRCA_normal.tsv", sep = "\t")
data2 <- fread("/Users/pranto/Desktop/BIONet/BRCA_tumor.tsv", sep = "\t")

# Measure execution time over multiple runs
execution_times <- numeric(10)

for (i in 1:10) {
  print(i)
  start_time <- Sys.time()
  distance <- distance_matrix(data1, data2)
  execution_times[i] <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
}

# Save execution times to a CSV file
execution_data <- data.frame(Execution = 1:10, Time = execution_times)
fwrite(execution_data, "/Users/pranto/Desktop/BIONet/R_execution_times.csv")

print("Saved execution times to /Users/pranto/Desktop/BIONet/R_execution_times.csv")
