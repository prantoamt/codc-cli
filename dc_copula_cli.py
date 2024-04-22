import numpy as np
import pandas as pd
import click
from tqdm import tqdm
from scipy.stats import ks_2samp

# Assuming GeneExpressionAnalyzer and EmpiricalCopula are defined in other modules or earlier in the script
from analyzer import GeneExpressionAnalyzer
from Copula.empirical_copula import EmpiricalCopula


@click.command()
@click.option(
    "--tumor_dataset",
    type=str,
    required=True,
    help="Path to the CSV file containing the tumor gene expression data.",
)
@click.option(
    "--normal_dataset",
    type=str,
    required=True,
    help="Path to the CSV file containing the normal gene expression data.",
)
@click.option(
    "--output_path",
    type=str,
    required=True,
    help="Output path to save the resulting CSV file containing the distance matrix.",
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
def main(
    tumor_dataset, normal_dataset, output_path, ties_method, smoothing, ks_stat_method
):
    """
    This script computes a matrix of differential coexpression scores for gene pairs across two conditions,
    using the Kolmogorov-Smirnov distance between their empirical copulas. It assesses the similarity
    in gene expression distributions between phenotypes such as 'tumor' and 'normal'.
    """
    df1 = pd.read_csv(tumor_dataset)
    df2 = pd.read_csv(normal_dataset)

    empirical_copula = EmpiricalCopula()
    analyzer = GeneExpressionAnalyzer(empirical_copula=empirical_copula)

    distance_matrix = analyzer.compute_dc_copula_matrix(
        df1,
        df2,
        ties_method=ties_method,
        smoothing=smoothing,
        ks_stat_method=ks_stat_method,
    )

    np.savetxt(output_path, distance_matrix, delimiter=",", fmt="%f")
    print(f"Saved the computed distance matrix to {output_path}")


if __name__ == "__main__":
    main()
