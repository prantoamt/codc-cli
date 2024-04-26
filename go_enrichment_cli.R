#!/usr/bin/env Rscript

# Load necessary libraries
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
})

# Function to parse command line arguments
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  parsed_args <- list()
  for (arg in args) {
    key_value <- strsplit(arg, "=")[[1]]
    if (length(key_value) == 2) {
      parsed_args[[sub("^--", "", key_value[1])]] <- key_value[2]
    } else if (length(key_value) == 1) {
      parsed_args[[sub("^--", "", key_value[1])]] <- TRUE
    }
  }
  return(parsed_args)
}

# Parse the command line arguments
args <- parse_args()

# Check for help request
if ('help' %in% names(args)) {
  cat("Usage: go_enrichment_cli.R --input_file=FILE --output_path=PATH [OPTIONS]\n")
  cat("Options:\n")
  cat("  --input_file=FILE           Path to the input TSV file with gene network data.\n")
  cat("  --output_path=PATH          Path to save the output plot image.\n")
  cat("  --threshold=VALUE           Weight threshold to filter gene pairs (default 0.6).\n")
  cat("  --key_type=TYPE             The key type of the gene identifiers (default SYMBOL). Other types can be 'ENSEMBL', 'ENTREZID', etc., \n")
  cat("  --ontology=TYPE             Ontology to use: BP, CC, or MF (default CC).\n")
  cat("  --p_adjust_method=METHOD    Method for adjusting p-values (default BH).\n")
  cat("  --qvalue_cutoff=VALUE       Q-value cutoff for significant enrichment (default 0.05).\n")
  cat("  --show_category=NUMBER      Number of categories to show in the bar plot (default 10).\n")
  cat("  --help                      Show this help message and exit.\n")
  quit(status = 0)
}

input_file_path <- if (!is.null(args[['input_file']])) args[['input_file']] else stop("Error: --input_file argument is required.\n")
output_path <- if (!is.null(args[['output_path']])) args[['output_path']] else stop("Error: --output_path argument is required.\n")
weight_threshold <- if (!is.null(args[['threshold']])) as.numeric(args[['threshold']]) else 0.6
key_type <- if (!is.null(args[['key_type']])) args[['key_type']] else 'SYMBOL'
ontology <- if (!is.null(args[['ontology']])) args[['ontology']] else 'CC'
adjust_method <- if (!is.null(args[['p_adjust_method']])) args[['p_adjust_method']] else 'BH'
qvalue_cutoff <- if (!is.null(args[['qvalue_cutoff']])) as.numeric(args[['qvalue_cutoff']]) else 0.05
show_category <- if (!is.null(args[['show_category']])) as.integer(args[['show_category']]) else 10

# Initialize progress bar
pb <- txtProgressBar(min = 0, max = 5, style = 3)
on.exit(close(pb))

# Step 1: Reading data
cat("\t Reading data from input file...\n")
Sys.sleep(1)  # simulate delay
setTxtProgressBar(pb, 1)
data <- read.csv(input_file_path, header = TRUE, sep = "\t")

# Step 2: Filtering data
cat("\t Filtering data based on weight threshold...\n")
Sys.sleep(1)  # simulate delay
setTxtProgressBar(pb, 2)
filtered_data <- data[data$Weight >= weight_threshold, ]
genes_to_test <- unique(c(filtered_data$Target, filtered_data$Regulator))

# Step 3: Performing GO enrichment analysis
cat("\t Analysing...\n")
Sys.sleep(1)  # simulate delay
setTxtProgressBar(pb, 3)
ego <- enrichGO(gene = genes_to_test, OrgDb = org.Hs.eg.db, keyType = key_type, 
                ont = ontology, pAdjustMethod = adjust_method, qvalueCutoff = qvalue_cutoff,
                readable = TRUE)

# Step 4: Preparing to save the plot
cat("\t Saving the plot to output file...\n")
Sys.sleep(1)  # simulate delay
setTxtProgressBar(pb, 4)
if (!dir.exists(dirname(output_path))) {
  dir.create(dirname(output_path), recursive = TRUE)
}

# Step 5: Plotting data
cat("\t Plotting data...\n")
Sys.sleep(1)  # simulate delay
setTxtProgressBar(pb, 5)
cat("\t Finished...\n")
setwd(dirname(output_path))  # Set the working directory to the output path's directory
output_file <- file.path(output_path, paste0("ontology", ".png"))

# Customize plot colors, axis labels, etc.
png(output_file, width=2000, height=1600, res=300)
par(mar=c(6.1, 5.1, 7.1, 3.1))  # Adjust margins to fit labels and title
barplot_title <- paste("GO Enrichment Analysis -", ontology, "Ontology")
barplot(ego, showCategory = show_category)+ggtitle(barplot_title)

# Close the device and suppress the output
invisible(dev.off())
cat("\n \n Plot saved to:", output_file, "\n")