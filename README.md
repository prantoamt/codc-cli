# CODC CLI

## Table of Contents
- [Brief Description](#brief-description)
- [Reference to the Publication](#reference-to-the-publication)
- [Available Commands](#available-commands)
- [Installation Instructions](#installation-instructions)
  - [Using Docker](#using-docker)
  - [Using Locally](#using-locally)
- [Execution Instructions Using Example Data](#execution-instructions-using-example-data)
  - [Using Docker](#using-docker-1)
  - [Using Locally](#using-locally-1)
- [Explanation of the Relevant Parameters](#explanation-of-the-relevant-parameters)
- [Input File Format Specification](#input-file-format-specification)
- [Output File Format Specification](#output-file-format-specification)
- [Explanation and Interpretation of the Output](#explanation-and-interpretation-of-the-output)
- [Recommended Hyperparameters by the Authors](#recommended-hyperparameters-by-the-authors)

## Brief Description
CODC is a Command-Line Interface (CLI) designed for analyzing gene expression data to identify differential coexpression between two conditions using a copula-based approach. 
This tool is reimplemented in python and R based on the author's R implementation due to calculation complexoty issue. Besides the coexpression calculation scripts, it also provides downstream and performace measurement scripts.

## Reference to the Publication
This tool implements the method proposed by Ray, S., Lall, S., & Bandyopadhyay, S. in "CODC: a Copula-based model to identify differential coexpression." You can read the study here: [npj Systems Biology and Applications (2020)](https://doi.org/10.1038/s41540-020-0137-9).

## Available Commands
The CLI includes commands for:
- Copula based differential co-expression calculation (`codc`)
- GO enrichment analysis (`go-enrichment`)
- Performance measurement of Python (`python-performance`)
- Performance measurement of R (`r-performance`)

This readme, explains Copula based differential co-expression calculation (`codc`).

## Installation Instructions
### Using Docker
```bash
docker build -t codc-tool .
```

### Using Locally
Install PDM (Python package manager) if not already installed:
```bash
pip install pdm
```
Then, install the packages using PDM:
```bash
pdm install
```

## Execution of codc Using Example Data
Before executing the commands below, go to the root directory of the CODC project. Outwise, you may have directory path mismatch.
The commands will output the `network.tsv`in `/tests/data/` directory
### Using Docker
```bash
docker run --rm -v /tests/data:/data codc-tool codc --input_file_1 /data/BRCA_normal.tsv --input_file_2 /data/BRCA_tumor.tsv --output_path /data
```

### Using Locally
```bash
pdm run codc --input_file_1 /tests/data/BRCA_normal.tsv --input_file_2 /tests/data/BRCA_tumor.tsv --output_path /tests/data/
```

## Explanation of the Relevant Parameters
#### `--input_file_1`
- **Description**: Path to the TSV file containing gene expression data for the first condition, often a disease state such as "tumor."
- **Required**: Yes
- **Example**: `--inputfile_1 /path/to/tumor_data.tsv`

#### `--input_file_2`
- **Description**: Path to the TSV file containing gene expression data for the second condition, typically a control or normal state.
- **Required**: Yes
- **Example**: `--inputfile_2 /path/to/normal_data.tsv`

#### `--output_path`
- **Description**: The directory where the output TSV file will be saved. This file will contain the computed differential coexpression network.
- **Required**: Yes
- **Example**: `--output_path /path/to/output`
- **Output Details**: The output is a TSV file named `network.tsv`, which includes columns for target gene, regulator gene, condition, and the weight of the coexpression difference.

#### `--ties_method`
- **Description**: Method to handle ties in data ranking within the pseudo-observations calculation.
- **Required**: No (default is "average")
- **Options**:
  - `average`: Average ranks of ties.
  - `max`: Use the maximum rank for ties.
- **Example**: `--ties_method max`

#### `--smoothing`
- **Description**: Specifies the smoothing technique applied to the empirical copula calculation.
- **Required**: No (default is "none")
- **Options**:
  - `none`: No smoothing applied.
  - `beta`: Use a beta smoothing approach.
  - `checkerboard`: Apply checkerboard smoothing.
- **Example**: `--smoothing beta`

#### `--ks_stat_method`
- **Description**: Determines the method used for computing the Kolmogorov-Smirnov statistic, which quantifies the differential coexpression.
- **Required**: No (default is "asymp")
- **Options**:
  - `asymp`: Use asymptotic properties of the KS statistic.
  - `auto`: Automatically determine the best method based on data characteristics.
  - `exact`: Compute an exact KS statistic.
- **Example**: `--ks_stat_method exact`

## Input File Format Specification
Input files must be in a tab-separated format with gene names in rows and sample IDs in columns. Example:

| Gene   | TCGA-A7-A0CE    | TCGA-A7-A0CH    |
|--------|-----------------|-----------------|
| ACTA1	| 6.872032023	   | 4.947203749     |
| MYL2	| 0.415445555	   | 0.0             |


## Output File Format Specification
The output `network.tsv` is a tab-separated file that includes:
- **Target**: Target gene of the edge.
- **Regulator**: Source gene of the edge.
- **Condition**: Describes the differential co-expression across conditions.
- **Weight**: Numerical value indicating the strength of the relationship.

Example output:

|  Target   |	Regulator   |  Condition                   |  Weight |
|-----------|--------------|------------------------------|---------|
|  MYL2    	|ACTA1	      |Diff Co-Exp of both Condition |	0.1111  |


## Explanation and Interpretation of the Output
The `network.tsv` output file lists gene pairs that are differentially coexpressed between two conditions, providing insights into gene interactions under different conditions.

## Recommended Hyperparameters by the Authors
There were no specific hyperparameters recommended by the authors. The default parameters used are based on typical settings derived from the author's R implementation:
- `ks_stat_method = asymp`
- `ties_method = average`
- `smoothing = none`