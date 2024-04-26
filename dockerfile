# Use an official Python runtime as a parent image
FROM python:3.11-slim

# Set the working directory in the container
WORKDIR /usr/src/app

# Install system libraries required for R and common packages
RUN apt-get update && apt-get install -y \
    r-base \
    r-base-dev \
    gfortran \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy the current directory contents into the container at /usr/src/app
COPY . .

# Install Python packages using PDM
RUN pip install pdm
RUN pdm install --prod

# Install R packages
RUN R -e "install.packages(c('clusterProfiler', 'org.Hs.eg.db', 'enrichplot'), repos='http://cran.rstudio.com/')"


# Define an entry point for the container; this is the command-line tool
ENTRYPOINT ["pdm", "run", "dc_copula_cli.py"]
# If additional arguments are required, they can be provided in the CMD instruction.
CMD ["--help"]