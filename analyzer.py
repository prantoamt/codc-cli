import numpy as np
import pandas as pd
from tqdm import tqdm

from scipy.stats import ks_2samp


class GeneExpressionAnalyzer:
    TIES_AVERAGE = "average"
    TIES_MIN = "min"
    TIES_MAX = "max"
    TIES_DENSE = "dense"
    TIES_ORIGINAL = "ordinal"

    def __init__(self, empirical_copula):
        """
        Initializes the GeneExpressionAnalyzer with an instance of EmpiricalCopula.

        Args:
            empirical_copula (EmpiricalCopula): An instance of EmpiricalCopula used to compute pseudo-observations
                                                 and empirical copulas.
        """
        self.empirical_copula = empirical_copula

    def compute_dc_copula_network(
        self,
        df1: pd.DataFrame,
        df2: pd.DataFrame,
        ties_method: str = "average",
        smoothing: str = "none",
        ks_stat_method: str = "asymp",
    ) -> np.ndarray:
        """
        Computes a network of differential coexpression scores for gene pairs across two conditions.
        The score for each gene pair is calculated using the Kolmogorov-Smirnov distance between their empirical
        copulas, representing the degree of differential coexpression. This network is used to assess the similarity
        in gene expression distributions between 'tumor' and 'normal' phenotypes for example.

        Args:
            df1 (pd.DataFrame): Gene expression data for the first phenotype (e.g., tumor),
                                where rows are genes and columns are samples. Expects the first column
                                to contain non-numeric identifiers which will be ignored.
            df2 (pd.DataFrame): Gene expression data for the second phenotype (e.g., normal),
                                structured like df1.
            ties_method (str): Specifies the method for ranking ties within the pseudo-observations.
                            Options are 'average', 'min', 'max', 'dense', 'ordinal'. Default is 'average'.
            smoothing (str): Specifies the type of smoothing to apply to the empirical copula.
                            Options are 'none', 'beta', 'checkerboard'. Default is 'none'.
            ks_stat_method (str): method parameter for the ks_2samp function, determining how the
                                Kolmogorov-Smirnov statistic is computed. Default is 'asymp'.

        Returns:
            np.ndarray: A network where:
                        First column target: Target of the edge
                        Second column regulator: Source of the edge
                        Third column condition: Condition that the edge belongs to
                        Fourth column weight: Weight of the edge
        """
        # Extract gene names from the first column
        gene_names = df1.iloc[:, 0].values
        assert np.array_equal(
            gene_names, df2.iloc[:, 0].values
        ), "Gene lists must match!"

        # Initialize an empty List for the network
        rows_list = []

        # Extracting numeric data from the dataframes, assuming the first column is the header
        data1 = df1.iloc[:, 1:].values
        data2 = df2.iloc[:, 1:].values
        n_genes = min(data1.shape[0], data2.shape[0])

        # Print dataset summary
        print(f"Starting DC Copula coexpression calculation:")
        print(f" - Number of gene pairs to be analyzed: {n_genes*n_genes}")
        print(f" - Ties method: {ties_method}")
        print(f" - Smoothing technique: {smoothing}")
        print(f" - KS statistic mode: {ks_stat_method}")

        # Create a tqdm progress bar
        pbar = tqdm(
            total=(n_genes * (n_genes - 1)) // 2,
            desc="Computing distances",
            unit="pair",
        )

        for i in range(n_genes - 1):
            for j in range(i + 1, n_genes):
                gene_pair_data1 = np.vstack((data1[i, :], data1[j, :])).T
                gene_pair_data2 = np.vstack((data2[i, :], data2[j, :])).T

                u1 = self.empirical_copula.pseudo_observations(
                    gene_pair_data1, ties_method
                )
                u2 = self.empirical_copula.pseudo_observations(
                    gene_pair_data2, ties_method
                )

                ec1 = self.empirical_copula.empirical_copula(
                    u1, gene_pair_data1, ties_method, smoothing
                )
                ec2 = self.empirical_copula.empirical_copula(
                    u2, gene_pair_data2, ties_method, smoothing
                )

                ks_stat, _ = ks_2samp(ec1, ec2, method=ks_stat_method)

                # Append each pair as a separate row in the List
                rows_list.append(
                    {
                        "Target": gene_names[j],
                        "Regulator": gene_names[i],
                        "Condition": "Diff Co-Exp of both Condition",
                        "Weight": ks_stat,
                    }
                )

                pbar.update(1)  # Update the progress bar after each iteration

        pbar.close()  # Close the progress bar when done
        # Creating a DataFrame from the list of rows
        network_df = pd.DataFrame(rows_list)
        return network_df
