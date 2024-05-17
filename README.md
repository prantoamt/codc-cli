# CODC CLI: Command-line interface to calculate copula-based differential gene co-expression

![Cover Image](data/images/cover.png)

<details>
<summary> Table of Contents </summary>

- [Brief Description](#brief-description)
- [Reference to the Publication](#reference-to-the-publication)
- [Methodology](#methodology)
- [Available Commands](#available-commands)
- [Installation Instructions](#installation-instructions)
  - [Using Docker](#using-docker)
  - [Using Locally](#or-using-locally)
- [Execution Of CODC Using BRCA Data](#execution-of-codc-using-brca-data)
  - [Using Docker](#using-docker-1)
  - [Using Locally](#or-using-locally-1)
- [Explanation of the Relevant Parameters](#explanation-of-the-relevant-parameters)
- [Input File Format Specification](#input-file-format-specification)
- [Output File Format Specification](#output-file-format-specification)
- [Explanation and Interpretation of the Output](#explanation-and-interpretation-of-the-output)
- [Recommended Hyperparameters by the Authors](#recommended-hyperparameters-by-the-authors)

</details>

## Brief Description
The CODC CLI tool is designed for analyzing gene expression data to calculate differential co-expression using a copula-based approach. It is implemented in Python based on the [R implementation](https://github.com/Snehalikalall/CODC/blob/master/distance_mat_calculation.R) of Ray et al.  for enhanced performance, with support for parallel processing. The tool allows users to compute differential co-expression networks and provides additional commands for downstream analysis and performance measurement. Installation can be done via Docker or locally using PDM, a Python package manager. The tool expects input files in TSV format and outputs the co-expression network as a TSV file as well. 

## Reference to the Publication
This tool implements the method proposed by Ray, S., Lall, S., & Bandyopadhyay, S. in ["CODC: a Copula-based model to identify differential co-expression."](https://doi.org/10.1038/s41540-020-0137-9).

## Methodology
The methodology to compute the copula based differential co-expression and mathematical explaination is detailed [here](downstream-analysis/methodology.md)

## Available Commands
The CLI includes commands for:
- Copula based differential co-expression calculation (`codc`)
- [GO enrichment analysis (`go-enrichment`)](downstream-analysis/go-enrichment.md)
- [Performance measurement of Python script (`python-performance`)](downstream-analysis/performance-measure.md)
- [Performance measurement of R script (`r-performance`)](downstream-analysis/performance-measure.md)

This readme, explains Copula based differential co-expression calculation (`codc`).

## Installation Instructions

### Clone the repository and go to the project root dir

Before installing and running the CLI tool, you have to clone the repo and navigate
to the project's root directory.

```bash
git clone git@github.com:bionetslab/grn-benchmark.git && cd grn-benchmark/src/codc-cli-tool
```

### Using Docker
```bash
docker build -t codc-tool .
```

### OR Using Locally
Install PDM (Python package manager) if not already installed:
```bash
pip install pdm
```
Then, install the packages using PDM:
```bash
pdm install
```

## Execution of codc Using BRCA Data

The commands below will output the `network.tsv`in `./data/` directory

### Using Docker
```bash
docker run --rm -v ./data:/data codc-tool codc --input_file_1 /data/BRCA_normal.tsv --input_file_2 /data/BRCA_tumor.tsv --output_path /data --batch_size 100
```

### OR Using Locally
```bash
pdm run cli codc --input_file_1 ./data/BRCA_normal.tsv --input_file_2 ./data/BRCA_tumor.tsv --output_path ./data --batch_size 100
```

## Explanation of the Relevant Parameters

#### `--input_file_1`
- **Description**: Path to the TSV file containing gene expression data for the first condition.
- **Required**: Yes
- **Example**: `--inputfile_1 /path/to/condition1.tsv`

#### `--input_file_2`
- **Description**: Path to the TSV file containing gene expression data for the second condition.
- **Required**: Yes
- **Example**: `--inputfile_2 /path/to/condition2.tsv`

#### `--output_path`
- **Description**: The directory where the output TSV file will be saved. This file will contain the computed differential co-expression network based on copula approach.
- **Required**: Yes
- **Example**: `--output_path /path/to/output`
- **Output Details**: The output is a TSV file named `network.tsv`, which includes columns for target gene, regulator gene, condition, and the weight as the co-expression difference.

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
- **Description**: Determines the method used for computing the Kolmogorov-Smirnov statistic, which quantifies the differential co-expression.
- **Required**: No (default is "asymp")
- **Options**:
  - `asymp`: Use asymptotic properties of the KS statistic.
  - `auto`: Automatically determine the best method based on data characteristics.
  - `exact`: Compute an exact KS statistic.
- **Example**: `--ks_stat_method exact`

#### `--batch_size`
- **Description**: Determines how many pair of genes will be executed in each batch in parallel execution.
- **Required**: No (default is 100)
- **Example**: `--batch_size 100`

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

|  Target   |	Regulator    |  Condition                         |  Weight |
|-----------|--------------|------------------------------------|---------|
|  MYL2    	|ACTA1	       | Diff Co-Exp between both Condition |	0.1111  |


## Explanation and Interpretation of the Output
The `network.tsv` output file lists gene pairs that are differentially coexpressed between two conditions, providing insights into gene interactions under different conditions.

## Recommended Hyperparameters by the Authors
There were no specific hyperparameters recommended by the authors. The default parameters used are based on typical settings derived from the author's R implementation:
- `ks_stat_method = asymp`
- `ties_method = average`
- `smoothing = none`

-------------------------------------------------------

Md Badiuzzaman Pranto, Friedrich-Alexander-Universität, Erlangen-Nürnberg.