import numpy as np
import pandas as pd
import click

from analyzer import GeneExpressionAnalyzer
from copula.empirical_copula import EmpiricalCopula


@click.command()
@click.option(
    "--input_file_1",
    type=str,
    required=True,
    help="Path to the TSV file containing gene expression data.",
)
@click.option(
    "--input_file_2",
    type=str,
    required=True,
    help="Path to the TSV file containing gene expression data.",
)
@click.option(
    "--output_path",
    type=str,
    required=True,
    help="Output path to store the resulting TSV file containing the differential coexpression network. The file name will be network.tsv",
)
@click.option(
    "--ties_method",
    type=click.Choice(["average", "max"]),
    default="average",
    help="Method for ranking ties within pseudo-observations.",
)
@click.option(
    "--smoothing",
    type=click.Choice(["none", "beta", "checkerboard"]),
    default="none",
    help="Type of smoothing to apply to the empirical copula.",
)
@click.option(
    "--ks_stat_method",
    type=click.Choice(["asymp", "auto", "exact"]),
    default="asymp",
    help="Mode parameter for the ks_2samp function, which determines how the Kolmogorov-Smirnov statistic is computed.",
)
def main(
    input_file_1, input_file_2, output_path, ties_method, smoothing, ks_stat_method
):
    """
    This script computes a network of differential coexpression scores for gene pairs across two conditions,
    using the Kolmogorov-Smirnov distance between their empirical copulas.
    It assesses the similarity in joint gene expression distributions between phenotypes such as 'tumor' and 'normal'
    (Could be other as well). For detail explaination: https://github.com/prantoamt/gene-coexpress/wiki
    """
    # Loading data from TSV files
    df1 = pd.read_csv(input_file_1, delimiter="\t")
    df2 = pd.read_csv(input_file_2, delimiter="\t")

    # Initializing the EmpiricalCopula and GeneExpressionAnalyzer instances
    empirical_copula = EmpiricalCopula()
    analyzer = GeneExpressionAnalyzer(empirical_copula=empirical_copula)

    # Computing the network using the specified methods
    network_df = analyzer.compute_dc_copula_network(
        df1,
        df2,
        ties_method=ties_method,
        smoothing=smoothing,
        ks_stat_method=ks_stat_method,
    )

    # Saving the network to the specified output path
    output_path = f"{output_path}/network.tsv"
    network_df.to_csv(output_path, sep="\t", index=False, header=True)
    print(f"Saved the computed network to {output_path}")


if __name__ == "__main__":
    main()
