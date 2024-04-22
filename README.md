# Gene-coexpress

This tool computes a matrix of differential coexpression scores for gene pairs across two phenotypic conditions (e.g., tumor vs. normal) using the Kolmogorov-Smirnov distance between their empirical copulas. It is designed to highlight differences in gene expression distributions, offering insights that are crucial for understanding phenotypic variations.

## Reference

This tool is based on the methods proposed in [CODC: a Copula-based model to identify differential coexpression.](https://doi.org/10.1038/s41540-020-0137-9). Please cite this paper if you use this tool in your research.

## Installation

The tool is implemented in Python and requires some external libraries. We recommend using `pdm` for managing dependencies and Python environments. Here’s how you can set it up:

1. Install PDM if you haven't already:
   ```bash
   pip install pdm
   ```

2. Clone the repository:
   ```bash
   git clone git@github.com:prantoamt/gene-coexpress.git
   cd gene-coexpress
   ```

3. Install dependencies using PDM:
   ```bash
   pdm install
   ```

## Usage

To run the tool, you need two TSV files containing gene expression data for the two conditions you are comparing. Use the following command to execute the script:

```bash
pdm run dc_copula_cli.py --inputfile_1 /path/to/your/file --inputfile_2 /path/to/your/file --output_path output/path
```

### Parameters

- `--inputfile_1`: Path to the TSV file containing gene expression data for the first condition (e.g., tumor).
- `--inputfile_2`: Path to the TSV file containing gene expression data for the second condition (e.g., normal).
- `--output_path`: Path where the output TSV file will be saved. The file name will be `network.tsv`.
- `--ties_method`: Method to handle ties in data ranking; options are 'average' and 'max'.
- `--smoothing`: Smoothing method for the empirical copula; options are 'none', 'beta', and 'checkerboard'.
- `--ks_stat_method`: Method to compute the Kolmogorov-Smirnov statistic; options include 'asymp', 'auto', and 'exact'.

## Input File Format

Each input file should be a TSV format containing:
- First column is named 'Gene' and contains the gene names
- All following columns are named after a sample/cell.

## Output File Format

The output is a TSV file named `network.tsv`, containing the differential coexpression network matrix:
- First column target: Target of the edge
- Second column regulator: Source of the edge
- Third column condition: Condition that the edge belongs to
- Fourth column weight: Weight of the edge