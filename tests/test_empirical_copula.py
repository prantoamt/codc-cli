import numpy as np
import pandas as pd
import pytest
from copula.empirical_copula import (
    EmpiricalCopula,
)  # Make sure to import your class correctly


def test_check_evaluation_points_valid():
    copula = EmpiricalCopula()
    valid_points = np.array([0.0, 0.5, 1.0])
    # Should pass without any exception
    copula.check_evaluation_points(valid_points)


def test_check_evaluation_points_invalid():
    copula = EmpiricalCopula()
    invalid_points = np.array([-0.1, 1.1])
    # Expecting ValueError for invalid points
    with pytest.raises(ValueError) as excinfo:
        copula.check_evaluation_points(invalid_points)
    assert "'evaluation_points' must be in [0,1]." in str(excinfo.value)


def test_pseudo_observations_avg_ties():
    # Example data, mimicking the expected input format and data
    data = pd.read_csv("./tests/data/BRCA_normal_subset.tsv", sep="\t")
    # Exclude the first column from the data
    data = data.iloc[:, 1:].values  # This selects all columns except the first one
    expected_output = np.array(
        [
            [
                0.90909091,
                0.81818182,
                0.72727273,
                0.72727273,
                0.81818182,
                0.81818182,
                0.81818182,
                0.72727273,
                0.18181818,
            ],
            [
                0.22727273,
                0.27272727,
                0.22727273,
                0.22727273,
                0.22727273,
                0.18181818,
                0.27272727,
                0.18181818,
                0.27272727,
            ],
            [
                0.22727273,
                0.54545455,
                0.22727273,
                0.45454545,
                0.45454545,
                0.36363636,
                0.45454545,
                0.54545455,
                0.90909091,
            ],
            [
                0.54545455,
                0.27272727,
                0.5,
                0.54545455,
                0.63636364,
                0.54545455,
                0.54545455,
                0.36363636,
                0.63636364,
            ],
            [
                0.09090909,
                0.27272727,
                0.22727273,
                0.22727273,
                0.22727273,
                0.18181818,
                0.13636364,
                0.18181818,
                0.09090909,
            ],
            [
                0.36363636,
                0.27272727,
                0.5,
                0.22727273,
                0.22727273,
                0.18181818,
                0.13636364,
                0.18181818,
                0.45454545,
            ],
            [
                0.81818182,
                0.63636364,
                0.63636364,
                0.63636364,
                0.90909091,
                0.63636364,
                0.63636364,
                0.63636364,
                0.54545455,
            ],
            [
                0.63636364,
                0.90909091,
                0.90909091,
                0.90909091,
                0.72727273,
                0.90909091,
                0.90909091,
                0.90909091,
                0.72727273,
            ],
            [
                0.45454545,
                0.27272727,
                0.22727273,
                0.22727273,
                0.22727273,
                0.45454545,
                0.36363636,
                0.45454545,
                0.36363636,
            ],
            [
                0.72727273,
                0.72727273,
                0.81818182,
                0.81818182,
                0.54545455,
                0.72727273,
                0.72727273,
                0.81818182,
                0.81818182,
            ],
        ]
    )

    # Instantiate the EmpiricalCopula object
    copula = EmpiricalCopula()

    # Calculate pseudo-observations
    result = copula.pseudo_observations(
        data=data, ties_method=EmpiricalCopula.TIES_AVERAGE
    )

    # Assert the output is as expected
    np.testing.assert_array_almost_equal(
        result,
        expected_output,
        decimal=6,
        err_msg="Pseudo-observations average ties do not match expected values.",
    )


# def test_pseudo_observations_max_ties():
#     # Example data, mimicking the expected input format and data
#     data = pd.read_csv("copula/tests/data/BRCA_normal_subset.csv")
#     # Exclude the first column from the data
#     data = data.iloc[:, 1:].values  # This selects all columns except the first one
#     expected_output = np.array(
#         [
#             [
#                 0.90909091,
#                 0.81818182,
#                 0.72727273,
#                 0.72727273,
#                 0.81818182,
#                 0.81818182,
#                 0.81818182,
#                 0.72727273,
#                 0.18181818,
#             ],
#             [
#                 0.22727273,
#                 0.27272727,
#                 0.22727273,
#                 0.22727273,
#                 0.22727273,
#                 0.18181818,
#                 0.27272727,
#                 0.18181818,
#                 0.27272727,
#             ],
#             [
#                 0.22727273,
#                 0.54545455,
#                 0.22727273,
#                 0.45454545,
#                 0.45454545,
#                 0.36363636,
#                 0.45454545,
#                 0.54545455,
#                 0.90909091,
#             ],
#             [
#                 0.54545455,
#                 0.27272727,
#                 0.5,
#                 0.54545455,
#                 0.63636364,
#                 0.54545455,
#                 0.54545455,
#                 0.36363636,
#                 0.63636364,
#             ],
#             [
#                 0.09090909,
#                 0.27272727,
#                 0.22727273,
#                 0.22727273,
#                 0.22727273,
#                 0.18181818,
#                 0.13636364,
#                 0.18181818,
#                 0.09090909,
#             ],
#             [
#                 0.36363636,
#                 0.27272727,
#                 0.5,
#                 0.22727273,
#                 0.22727273,
#                 0.18181818,
#                 0.13636364,
#                 0.18181818,
#                 0.45454545,
#             ],
#             [
#                 0.81818182,
#                 0.63636364,
#                 0.63636364,
#                 0.63636364,
#                 0.90909091,
#                 0.63636364,
#                 0.63636364,
#                 0.63636364,
#                 0.54545455,
#             ],
#             [
#                 0.63636364,
#                 0.90909091,
#                 0.90909091,
#                 0.90909091,
#                 0.72727273,
#                 0.90909091,
#                 0.90909091,
#                 0.90909091,
#                 0.72727273,
#             ],
#             [
#                 0.45454545,
#                 0.27272727,
#                 0.22727273,
#                 0.22727273,
#                 0.22727273,
#                 0.45454545,
#                 0.36363636,
#                 0.45454545,
#                 0.36363636,
#             ],
#             [
#                 0.72727273,
#                 0.72727273,
#                 0.81818182,
#                 0.81818182,
#                 0.54545455,
#                 0.72727273,
#                 0.72727273,
#                 0.81818182,
#                 0.81818182,
#             ],
#         ]
#     )

#     # Instantiate the EmpiricalCopula object
#     copula = EmpiricalCopula()

#     # Calculate pseudo-observations
#     result = copula.pseudo_observations(data=data, ties_method=EmpiricalCopula.TIES_MAX)

#     # Assert the output is as expected
#     np.testing.assert_array_almost_equal(
#         result,
#         expected_output,
#         decimal=6,
#         err_msg="Pseudo-observations max ties do not match expected values.",
#     )


def test_empirical_distribution_function():
    # Load the data from the specified CSV file
    data = pd.read_csv("./tests/data/BRCA_normal_subset.tsv", sep="\t")

    # Exclude the first column from the data
    data = data.iloc[:, 1:].values  # This selects all columns except the first one
    expected_output = np.array([0.2, 0.2, 0.3, 0.4, 0.1, 0.2, 0.5, 0.6, 0.3, 0.5])
    # Instantiate the EmpiricalCopula object
    copula = EmpiricalCopula()

    # Calculate pseudo-observations
    pobs = copula.pseudo_observations(
        data=data, ties_method=EmpiricalCopula.TIES_AVERAGE
    )
    result = copula.empirical_distribution_function(pobs, pobs)
    # Assert the output is as expected
    np.testing.assert_array_almost_equal(
        result,
        expected_output,
        decimal=6,
        err_msg="empirical-distribution-function do not match expected values.",
    )


def test_empirical_copula_avg_ties():
    # Load the data from the specified CSV file
    data = pd.read_csv("./tests/data/BRCA_normal_subset.tsv", sep="\t")

    # Exclude the first column from the data
    data = data.iloc[:, 1:].values  # This selects all columns except the first one
    expected_output = np.array([0.2, 0.2, 0.3, 0.4, 0.1, 0.2, 0.5, 0.6, 0.3, 0.5])
    # Instantiate the EmpiricalCopula object
    copula = EmpiricalCopula()

    # Calculate pseudo-observations
    pobs = copula.pseudo_observations(
        data=data, ties_method=EmpiricalCopula.TIES_AVERAGE
    )
    result = copula.empirical_copula(evaluation_points=pobs, data=data)
    # Assert the output is as expected
    np.testing.assert_array_almost_equal(
        result,
        expected_output,
        decimal=6,
        err_msg="empirical-copula avg ties do not match expected values.",
    )


def test_empirical_copula_max_ties():
    # Load the data from the specified CSV file
    data = pd.read_csv("./tests/data/BRCA_normal_subset.tsv", sep="\t")

    # Exclude the first column from the data
    data = data.iloc[:, 1:].values  # This selects all columns except the first one
    expected_output = np.array([0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.6, 0.0, 0.5])
    # Instantiate the EmpiricalCopula object
    copula = EmpiricalCopula()

    # Calculate pseudo-observations
    pobs = copula.pseudo_observations(
        data=data, ties_method=EmpiricalCopula.TIES_AVERAGE
    )
    result = copula.empirical_copula(
        evaluation_points=pobs, data=data, ties_method=EmpiricalCopula.TIES_MAX
    )
    # Assert the output is as expected
    np.testing.assert_array_almost_equal(
        result,
        expected_output,
        decimal=6,
        err_msg="empirical-copula max ties do not match expected values.",
    )
