# Use Bioconductor base image
FROM bioconductor/bioconductor_docker:RELEASE_3_17

# Set working directory
WORKDIR /usr/src/app

# Install Python and system dependencies
RUN apt-get update && apt-get install -y \
    python3-pip \
    python3-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libblas3 \
    liblapack3 \
    libicu-dev \
    ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Update pip and install Python package manager (e.g., pip or Poetry)
RUN pip3 install --no-cache-dir --upgrade pip \
    && pip3 install --no-cache-dir pdm

# Copy R dependency installation script and run it
COPY dependencies.R ./dependencies.R
RUN Rscript dependencies.R

# Copy Python project files
COPY . .

# Install Python dependencies
RUN pdm install --prod

# Define the Docker entrypoint
ENTRYPOINT ["pdm", "run", "cli"]

# Default command
CMD ["pdm", "run", "cli", "--help"]
