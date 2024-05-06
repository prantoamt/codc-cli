# Stage 1: R Environment
FROM r-base:4.1.0 as stage1
WORKDIR /usr/src/app

# Install necessary system libraries
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libudunits2-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    ca-certificates

# Install R packages via BiocManager
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org/'); \
    BiocManager::install(c('clusterProfiler', 'org.Hs.eg.db', 'enrichplot'));"

# List installed things in R directory
RUN ls /usr/local/lib/R/site-library

# Stage 2: Python Environment
FROM python:3.11-slim as final-stage
WORKDIR /usr/src/app

# Copy R installation from the first stage
COPY --from=stage1 /usr/local /usr/local

# Make sure the R executable path is recognized
ENV PATH="/usr/local/bin:$PATH"

# Check Rscript availability
RUN which Rscript
RUN Rscript --version

# # Validate R installations
# RUN Rscript -e "if (!requireNamespace('clusterProfiler', quietly = TRUE)) stop('clusterProfiler not found', call. = FALSE);" \
#     && Rscript -e "if (!requireNamespace('org.Hs.eg.db', quietly = TRUE)) stop('org.Hs.eg.db not found', call. = FALSE);" \
#     && Rscript -e "if (!requireNamespace('enrichplot', quietly = TRUE)) stop('enrichplot not found', call. = FALSE);"

# Install system dependencies for Python and validation checks
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Python package manager
RUN pip install --no-cache-dir pdm

# Copy Python project files
COPY . .

# Install Python dependencies
RUN pdm install --prod


# Define the Docker entrypoint
ENTRYPOINT ["pdm", "run"]

# Default command
CMD ["cli", "--help"]
