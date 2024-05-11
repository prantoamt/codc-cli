# Stage 1: R Environment
FROM bioconductor/bioconductor_docker:RELEASE_3_17 as stage1
WORKDIR /usr/src/app

# Copy the dependencies file
COPY dependencies.R ./dependencies.R

# Install R packages specified in the dependencies.R script
RUN Rscript dependencies.R

# List installed things in R directory
RUN ls /usr/local/lib/R/site-library

# Stage 2: Python Environment
FROM python:3.11-slim as final-stage
WORKDIR /usr/src/app

# Copy R installation from the first stage
COPY --from=stage1 /usr/local/lib/R /usr/local/lib/R

# Set the environment variable for R home directory and include R binary in PATH
ENV R_HOME=/usr/local/lib/R
ENV PATH="${R_HOME}/bin:${PATH}"

# Check Rscript availability
RUN which Rscript
RUN Rscript --version

# Install system dependencies for Python
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
ENTRYPOINT ["pdm", "run", "cli"]

# Default command
CMD ["pdm", "run", "cli", "--help"]
