# dependencies.R

# Helper function to install CRAN and Bioconductor packages
install_packages <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      message(paste("Installing package:", pkg))
      if (pkg %in% bioconductor_packages) {
        BiocManager::install(pkg)
      } else {
        install.packages(pkg, dependencies = TRUE)
      }
    } else {
      message(paste("Package already installed:", pkg))
    }
  }
}

# Load or install BiocManager for installing Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  library(BiocManager)
}

# Define CRAN and Bioconductor packages
bioconductor_packages <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot")
cran_packages <- c("copula", "stats", "data.table", "optparse")

# Install packages using the helper function
install_packages(c(bioconductor_packages, cran_packages))

# Verify installations
message("Verifying package installations...")
required_packages <- c(bioconductor_packages, cran_packages)

# Loop through and verify
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    stop(paste("Package not installed:", pkg))
  }
}

message("All packages installed successfully.")
