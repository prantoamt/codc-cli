# Gene-Coexpress

**Gene-coexpress** is a computational tool designed to analyze gene expression data and identify differential coexpression between two phenotypic conditions, such as tumor versus normal tissue. This analysis is based on the Kolmogorov-Smirnov distance measured between the empirical copulas of gene pairs, highlighting significant changes in gene expression distributions.

## Reference

This tool implements methods described in the paper ["CODC: a Copula-based model to identify differential coexpression"](https://doi.org/10.1038/s41540-020-0137-9). If you find this tool useful in your research, please cite this publication.

## Quick Start

### Installation

Gene-coexpress is implemented in Python and utilizes `pdm` for dependency management. To set up the tool, follow these steps:

1. **Install PDM**:
   ```bash
   pip install pdm
   ```
2. **Clone the repository and navigate to the directory**:
   ```bash
   git clone git@github.com:prantoamt/gene-coexpress.git
   cd gene-coexpress
   ```
3. **Install dependencies**:
   ```bash
   pdm install
   ```

### Basic Usage

Run the tool using the command below, specifying the paths to your input TSV files and the desired output path for the results:
```bash
pdm run dc_copula_cli.py --inputfile_1 /path/to/tumor_data.tsv --inputfile_2 /path/to/normal_data.tsv --output_path /path/to/output/network.tsv
```

For more detailed usage examples, parameter explanations, and tutorials, please see the [Usage Examples](https://github.com/prantoamt/gene-coexpress/wiki/Usage-Examples) section of our project Wiki.

## Documentation

For a comprehensive guide on setting up, running, and interpreting the output of Gene-Coexpress, refer to the [project Wiki](https://github.com/prantoamt/gene-coexpress/wiki).

## Support

If you encounter any problems or have suggestions, please open an issue on the [GitHub issue tracker](https://github.com/prantoamt/gene-coexpress/issues).