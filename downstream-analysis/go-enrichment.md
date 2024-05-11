# GO Enrichment Analysis

This document provides instructions for running the Gene Ontology (GO) enrichment analysis script using the R programming language or Docker. The script performs enrichment analysis on gene sets to identify significant biological terms associated with gene lists derived from differential coexpression analysis (network.tsv).

<details>
<summary> Table of Contents </summary>

- [Installation](#installation)
  - [Using Docker](#using-docker)
  - [Using Locally](#using-locally)
- [Usage](#usage)
  - [Command Line Arguments](#command-line-arguments)
  - [Running with Docker](#running-with-docker)
- [Output Description](#output-description)
- [Contributing](#contributing)

</details>

## Installation

### Using Docker

To use the script with Docker, first build the Docker image if not done already:

```bash
docker build -t codc-tool .
```

### Using Locally

Ensure R is installed on your system, then install the required R packages:

```R
Rscript dependencies.R
```

## Usage

### Command Line Arguments

The script accepts the following arguments:

- `--input_file`: Path to the input TSV file containing gene network data (network.tsv).
- `--output_path`: Directory path to save the output plot image.
- `--threshold`: Weight threshold to filter gene pairs. Default is 0.6.
- `--key_type`: The key type of the gene identifiers. Default is 'SYMBOL'.
- `--ontology`: Ontology to use: BP, CC, or MF. Default is 'CC'.
- `--p_adjust_method`: Method for adjusting p-values. Default is 'BH'.
- `--qvalue_cutoff`: Q-value cutoff for significant enrichment. Default is 0.05.
- `--show_category`: Number of categories to show in the bar plot. Default is 10.

### Running with Docker

To run the GO enrichment analysis using Docker, use the following command:

```bash
docker run --rm -v /path/to/data:/data codc-tool go-enrichment --input_file=/data/network.tsv --output_path=/data --threshold=0.5 --key_type=SYMBOL --ontology=BP --p_adjust_method=BH --qvalue_cutoff=0.05 --show_category=20
```

Replace `/path/to/data` with the actual path to your data directory, and adjust the command-line arguments as needed.

## Output Description

The script will output a plot image in PNG format showing the results of the GO enrichment analysis. This image includes a bar plot visualizing the significant GO terms associated with the gene list analyzed.