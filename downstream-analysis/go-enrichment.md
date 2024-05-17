# GO Enrichment Analysis

This document provides instructions for running the Gene Ontology (GO) enrichment analysis script using Docker. The script performs enrichment analysis on gene sets to identify significant biological terms associated with gene lists derived from differential coexpression analysis (network.tsv).

<details>
<summary> Table of Contents </summary>

- [Installation](#installation)
- [Usage](#usage)
  - [Running with Docker](#running-with-docker)
  - [Command Line Arguments](#command-line-arguments)
- [Output Description](#output-description)

</details>

## Installation

### Clone the repository and go to the project root dir

Before installing and running the CLI tool, you have to clone the repo and navigate
to the project's root directory.

```bash
git clone git@github.com:bionetslab/grn-benchmark.git && cd grn-benchmark/src/codc-cli-tool
```

To use the script with Docker, first build the Docker image if not done already:

```bash
docker build -t codc-tool .
```

## Usage

### Running with Docker

To run the GO enrichment analysis using Docker, use the following command:

```bash
docker run --rm -v ./data:/data codc-tool go-enrichment --input_file=/data/network.tsv --output_path=/data --threshold=0.6 --key_type=SYMBOL --ontology=BP --p_adjust_method=BH --qvalue_cutoff=0.05 --show_category=10
```

Replace `./data` with the path to your data directory where `network.tsv` exists if it is different than `./data`, and adjust the command-line arguments as needed.

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

## Output Description

The script will output a plot image in PNG format showing the results of the GO enrichment analysis. This image includes a bar plot visualizing the significant GO terms associated with the gene list analyzed.