# Use an official Python runtime as a parent image
FROM python:3.11-slim

# Set the working directory in the container
WORKDIR /usr/src/app

# Copy the current directory contents into the container at /usr/src/app
COPY . .

# Install PDM
RUN pip install pdm

# Install any needed packages specified in pyproject.toml
RUN pdm install --prod

# Define an entry point for the container; this is the command-line tool
ENTRYPOINT ["pdm", "run", "dc_copula_cli.py"]
# If additional arguments are required, they can be provided in the CMD instruction.
CMD ["--help"]
