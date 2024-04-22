import numpy as np
import pandas as pd
import click

from analyzer import GeneExpressionAnalyzer
from copula.empirical_copula import EmpiricalCopula


@click.command()
@click.option(
    "--inputfile_1",
    type=str,
    required=True,
    help="Path to the TSV file containing the tumor gene expression data.",
)
@click.option(
    "--inputfile_2",
    type=str,
    required=True,
    help="Path to the TSV file containing the normal gene expression data.",
)
@click.option(
    "--output_path",
    type=str,
    required=True,
    help="Output path to save the resulting TSV file containing the distance matrix.",
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
    help="Mode parameter for the ks_2samp function.",
)
def main(inputfile_1, inputfile_2, output_path, ties_method, smoothing, ks_stat_method):
    """
    This script computes a matrix of differential coexpression scores for gene pairs across two conditions,
    using the Kolmogorov-Smirnov distance between their empirical copulas. It assesses the similarity
    in gene expression distributions between phenotypes such as 'tumor' and 'normal'.
    """
    df1 = pd.read_csv(inputfile_1, delimiter="\t")
    df2 = pd.read_csv(inputfile_2, delimiter="\t")

    empirical_copula = EmpiricalCopula()
    analyzer = GeneExpressionAnalyzer(empirical_copula=empirical_copula)

    distance_matrix = analyzer.compute_dc_copula_matrix(
        df1,
        df2,
        ties_method=ties_method,
        smoothing=smoothing,
        ks_stat_method=ks_stat_method,
    )

    np.savetxt(output_path, distance_matrix, delimiter="\t", fmt="%f")
    print(f"Saved the computed distance matrix to {output_path}")


if __name__ == "__main__":
    main()
