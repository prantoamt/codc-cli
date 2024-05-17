import pytest
import pandas as pd
import numpy as np
from analyzer import GeneExpressionAnalyzer
from copula.empirical_copula import EmpiricalCopula


@pytest.fixture
def setup_data():
    df1 = pd.DataFrame(
        {
            "Gene": ["Gene1", "Gene2", "Gene3"],
            "Sample1": [10, 25, 30],
            "Sample2": [4, 50, 60],
        }
    )

    df2 = pd.DataFrame(
        {
            "Gene": ["Gene1", "Gene2", "Gene3"],
            "Sample1": [3, 2, 1],
            "Sample2": [4, 2, 3],
        }
    )

    copula = EmpiricalCopula()
    analyzer = GeneExpressionAnalyzer(empirical_copula=copula)
    return df1, df2, analyzer


def test_output_format(setup_data):
    df1, df2, analyzer = setup_data
    network_df = analyzer.compute_dc_copula_network_parallel(df1, df2)
    assert isinstance(network_df, pd.DataFrame)
    assert set(network_df.columns) == {"Target", "Regulator", "Condition", "Weight"}


def test_output_values(setup_data):
    df1, df2, analyzer = setup_data
    network_df = analyzer.compute_dc_copula_network_parallel(df1, df2)
    assert len(network_df) == 3  # n(n-1)/2 pairs for n = 3
    assert all(network_df["Weight"] >= 0)  # Ensure all weights are non-negative


def test_error_on_mismatched_genes(setup_data):
    df1, df2, analyzer = setup_data
    df_mismatched = df2.copy()
    df_mismatched["Gene"] = ["Gene1", "Gene2", "Gene4"]  # 'Gene4' instead of 'Gene3'
    with pytest.raises(AssertionError):
        analyzer.compute_dc_copula_network_parallel(df1, df_mismatched)


@pytest.mark.parametrize("method", ["asymp", "auto", "exact"])
def test_ks_stat_method(setup_data, method):
    df1, df2, analyzer = setup_data
    network_df = analyzer.compute_dc_copula_network_parallel(
        df1, df2, ks_stat_method=method
    )
    assert len(network_df) == 3  # Check for the correct number of pairs processed


def test_compute_dc_copula_network_parallel_result():
    # Load your test data
    df1 = pd.read_csv("./tests/data/BRCA_normal_subset.tsv", sep="\t")
    df2 = pd.read_csv("./tests/data/BRCA_tumor_subset.tsv", sep="\t")
    empirical_copula = EmpiricalCopula()
    analyzer = GeneExpressionAnalyzer(empirical_copula=empirical_copula)

    # Perform the analysis
    network_df = analyzer.compute_dc_copula_network_parallel(
        df1, df2, ks_stat_method="asymp"
    )

    # Expected results prepared
    expected_network_df = pd.DataFrame(
        {
            "Target": [
                "MYL2",
                "ACTN2",
                "CKM",
                "MYL1",
                "MYH2",
                "ATP2A1",
                "KBTBD10",
                "NRAP",
                "NEB",
                "ACTN2",
                "CKM",
                "MYL1",
                "MYH2",
                "ATP2A1",
                "KBTBD10",
                "NRAP",
                "NEB",
                "CKM",
                "MYL1",
                "MYH2",
                "ATP2A1",
                "KBTBD10",
                "NRAP",
                "NEB",
                "MYL1",
                "MYH2",
                "ATP2A1",
                "KBTBD10",
                "NRAP",
                "NEB",
                "MYH2",
                "ATP2A1",
                "KBTBD10",
                "NRAP",
                "NEB",
                "ATP2A1",
                "KBTBD10",
                "NRAP",
                "NEB",
                "KBTBD10",
                "NRAP",
                "NEB",
                "NRAP",
                "NEB",
                "NEB",
            ],
            "Regulator": [
                "ACTA1",
                "ACTA1",
                "ACTA1",
                "ACTA1",
                "ACTA1",
                "ACTA1",
                "ACTA1",
                "ACTA1",
                "ACTA1",
                "MYL2",
                "MYL2",
                "MYL2",
                "MYL2",
                "MYL2",
                "MYL2",
                "MYL2",
                "MYL2",
                "ACTN2",
                "ACTN2",
                "ACTN2",
                "ACTN2",
                "ACTN2",
                "ACTN2",
                "ACTN2",
                "CKM",
                "CKM",
                "CKM",
                "CKM",
                "CKM",
                "CKM",
                "MYL1",
                "MYL1",
                "MYL1",
                "MYL1",
                "MYL1",
                "MYH2",
                "MYH2",
                "MYH2",
                "MYH2",
                "ATP2A1",
                "ATP2A1",
                "ATP2A1",
                "KBTBD10",
                "KBTBD10",
                "NRAP",
            ],
            "Condition": ["Diff Co-Exp of both Condition"] * 45,
            "Weight": [
                0.111111,
                0.111111,
                0.222222,
                0.000000,
                0.222222,
                0.222222,
                0.222222,
                0.222222,
                0.111111,
                0.222222,
                0.333333,
                0.555556,
                0.888889,
                0.111111,
                0.222222,
                0.333333,
                0.444444,
                0.222222,
                0.000000,
                0.111111,
                0.222222,
                0.222222,
                0.222222,
                0.222222,
                0.222222,
                0.444444,
                0.111111,
                0.333333,
                0.333333,
                0.333333,
                0.555556,
                0.000000,
                0.000000,
                0.555556,
                0.000000,
                0.222222,
                0.333333,
                0.555556,
                0.222222,
                0.333333,
                0.222222,
                0.222222,
                0.333333,
                0.222222,
                0.333333,
            ],
        }
    )

    # Sort dataframes by multiple columns to ensure consistent ordering for comparison
    network_df = network_df.sort_values(
        by=["Target", "Regulator", "Weight"]
    ).reset_index(drop=True)
    expected_network_df = expected_network_df.sort_values(
        by=["Target", "Regulator", "Weight"]
    ).reset_index(drop=True)

    # Use pandas testing assert function to compare dataframes
    pd.testing.assert_frame_equal(network_df, expected_network_df)
