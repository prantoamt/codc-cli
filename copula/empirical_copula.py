import numpy as np
from scipy.stats import rankdata, beta
from typing import Optional


class EmpiricalCopula:
    """
    A class used to calculate empirical copula and related statistics for gene expression data or similar datasets.
    """

    TIES_AVERAGE = "average"
    TIES_MIN = "min"
    TIES_MAX = "max"
    TIES_DENSE = "dense"
    TIES_ORIGINAL = "ordinal"

    def __init__(self):
        pass  # Initialization can include setting up logging or other utilities if needed

    def check_evaluation_points(self, evaluation_points: np.ndarray) -> None:
        """
        Validates whether the evaluation points are within the valid range [0, 1]. This function ensures that
        all provided points are suitable for processes that require normalized data inputs, such as calculating
        empirical copulas.

        Args:
            evaluation_points (np.ndarray): An array of evaluation points to be checked.

        Raises:
            ValueError: If any of the evaluation points are outside the interval [0, 1].

        Examples:
            >>> check_evaluation_points(np.array([0.1, 0.5, 0.9]))
            None  # No exception is raised

            >>> check_evaluation_points(np.array([-0.1, 1.1]))
            ValueError: 'u' must be in [0,1].
        """
        if np.any(evaluation_points < 0) or np.any(evaluation_points > 1):
            raise ValueError("'evaluation_points' must be in [0,1].")

    def pseudo_observations(
        self,
        data: np.ndarray,
        ties_method: Optional[str] = TIES_AVERAGE,
    ) -> np.ndarray:
        """
        Converts data into pseudo-observations based on the ranks of the data points within each column.
        Pseudo-observations are computed as the rank of each data point divided by (n+1) where n is the
        number of data points in each column.

        Args:
            data (np.ndarray): A 2D array where each row represents a data point and each column a variable.
            ties_method (str): Method to use for ranking data points that have the same value. Accepts
                            'average', 'min', 'max', 'dense', and 'ordinal' as specified in scipy.stats.rankdata.

        Returns:
            np.ndarray: A 2D array of the same shape as `data` containing the pseudo-observations of the input data.
        """

        num_samples, num_variables = data.shape
        ranks = np.zeros_like(data, dtype=float)

        # Apply rankdata to each column individually (still a loop but using vectorized operations where possible)
        for i in range(num_variables):
            ranks[:, i] = rankdata(data[:, i], method=ties_method)

        # Compute pseudo-observations in a vectorized manner
        pseudo_observations = ranks / (num_samples + 1)
        return pseudo_observations

    def empirical_distribution_function(
        self, evaluation_points: np.ndarray, data: np.ndarray
    ) -> np.ndarray:
        """
        Computes the empirical distribution function (EDF) of a given dataset at specified evaluation points.

        The empirical distribution function at a point u_k is the proportion of data points in the dataset
        that are less than or equal to u_k in all dimensions.

        Args:
            evaluation_points (np.ndarray): An array of points where the EDF is to be evaluated. Each row
                                            corresponds to a point in the multidimensional data space.
            data (np.ndarray): The dataset array where each row is a data point and each column corresponds to
                            a dimension of the data space.

        Returns:
            np.ndarray: An array containing the empirical CDF values at each of the evaluation points.
        """

        num_samples = data.shape[0]  # Number of data points in the dataset
        empirical_cdf = np.zeros(
            len(evaluation_points)
        )  # Initialize the array for empirical CDF values

        # Calculate the empirical CDF for each evaluation point
        for point_index, eval_point in enumerate(evaluation_points):
            # Count how many data points are less than or equal to the evaluation point across all dimensions
            empirical_cdf[point_index] = (
                np.sum(np.all(data <= eval_point, axis=1)) / num_samples
            )

        return empirical_cdf

    def empirical_copula(
        self,
        evaluation_points: np.ndarray,
        data: np.ndarray,
        ties_method: str = "average",
        smoothing: Optional[str] = "none",
    ) -> np.ndarray:
        """
        Computes the empirical copula of a given dataset at specified evaluation points using
        different smoothing methods if required. This function transforms the data into pseudo-observations
        first and then computes the empirical distribution function or smoothed version based on the selected
        method.

        Args:
            evaluation_points (np.ndarray): Points at which the copula is to be evaluated. Each row corresponds
                                            to a point in the multidimensional unit cube [0, 1]^d.
            data (np.ndarray): A 2D array of data points where each row is an observation and each column is a variable.
            ties_method (str): Method to use for ranking data points that have the same value in the computation of
                            pseudo-observations. Options include 'average', 'min', 'max', 'dense', and 'ordinal'.
            smoothing (Optional[str]): Specifies the type of smoothing to apply to the empirical copula. Options are
                                    'none', 'beta', and 'checkerboard'. Defaults to 'none' which computes the plain
                                    empirical distribution function.

        Returns:
            np.ndarray: An array containing the empirical copula values at each of the evaluation points using the
                        selected smoothing method.

        Raises:
            ValueError: If the smoothing method provided is not supported.
        """
        # Check the validity of evaluation points
        self.check_evaluation_points(evaluation_points)

        # Convert data to pseudo-observations
        data_pseudo_observations = self.pseudo_observations(
            data, ties_method=ties_method
        )

        # Compute the empirical distribution function based on the selected smoothing method
        if smoothing == "none":
            return self.empirical_distribution_function(
                evaluation_points, data_pseudo_observations
            )
        elif smoothing == "beta":
            return self.beta_smoothed_edf(evaluation_points, data_pseudo_observations)
        elif smoothing == "checkerboard":
            return self.checkerboard_smoothing(
                evaluation_points, data_pseudo_observations
            )
        else:
            raise ValueError(f"Unsupported smoothing method: {smoothing}")

    def beta_smoothed_edf(
        self, evaluation_points: np.ndarray, data_matrix: np.ndarray
    ) -> np.ndarray:
        """
        Computes a smoothed empirical distribution function (EDF) for the input array `evaluation_points` based
        on the data matrix `data_matrix` using Beta distribution smoothing. Each dimension of `data_matrix` contributes
        to the smoothing independently.

        Args:
        evaluation_points (np.ndarray): An array of points where the smoothed EDF is evaluated. Each row represents
                                        a point and each column a dimension.
        data_matrix (np.ndarray): A data matrix where each row is an observation and each column is a dimension.

        Returns:
        np.ndarray: An array of smoothed EDF values for each point in `evaluation_points`.
        """

        # Number of observations and dimensions in the data_matrix
        num_observations, num_dimensions = data_matrix.shape

        # Array to store the smoothed EDF values for each point in `evaluation_points`
        smoothed_edf_values = np.zeros(len(evaluation_points))

        # Iterate over each point in `evaluation_points`
        for point_idx in range(len(evaluation_points)):
            # Array to store Beta CDF values for the current point across all dimensions
            dimension_cdfs = np.zeros(num_dimensions)

            # Process each dimension independently
            for dimension_idx in range(num_dimensions):
                # Compute ranks for the current dimension with 'average' tie-breaking
                ranks = rankdata(data_matrix[:, dimension_idx], method="average")

                # Alpha and Beta parameters for the Beta distributions
                alpha_params = ranks  # Alpha parameters based on ranks
                beta_params = (
                    num_observations - ranks + 1
                )  # Beta parameters based on the complement of ranks

                # Compute the Beta CDF at the current point's dimension value and take the mean across all samples
                dimension_cdfs[dimension_idx] = beta.cdf(
                    evaluation_points[point_idx, dimension_idx],
                    alpha_params,
                    beta_params,
                ).mean()

            # Average the Beta CDF values across all dimensions for the current point
            smoothed_edf_values[point_idx] = dimension_cdfs.mean()

        return smoothed_edf_values

    def checkerboard_smoothing(
        self,
        evaluation_points: np.ndarray,
        data_matrix: np.ndarray,
        num_blocks: int = 10,
    ) -> np.ndarray:
        """
        Computes a smoothed empirical distribution function (EDF) for the input array `evaluation_points` using a checkerboard
        partitioning approach on the data matrix `data_matrix`. Each dimension is divided into blocks (bins), and the count
        of points falling into each block is used to estimate the distribution.

        Args:
            evaluation_points (np.ndarray): Points at which the smoothed EDF is evaluated, shaped (m, d) where m is the
                                            number of points and d is the number of dimensions.
            data_matrix (np.ndarray): Data matrix where each row is an observation and each column is a dimension, shaped (n, d).
            num_blocks (int): Number of blocks to partition each dimension of the data space.

        Returns:
            np.ndarray: An array of smoothed EDF values for each point in `evaluation_points`, one per point.
        """

        num_samples, num_dimensions = data_matrix.shape
        block_size = 1.0 / num_blocks
        smoothed_edf_values = np.zeros_like(
            evaluation_points[:, 0]
        )  # Initialize smoothed EDF array for each evaluation point

        # Iterate over each dimension to perform block smoothing
        for dimension_idx in range(num_dimensions):
            block_counts = np.zeros(
                num_blocks
            )  # Array to count the number of data points in each block for this dimension

            # Count data points in each block
            for data_point in data_matrix[:, dimension_idx]:
                block_index = int(data_point // block_size)
                block_counts[block_index] += 1

            # Calculate the smoothed EDF values for each evaluation point in the current dimension
            for point_idx, eval_point in enumerate(evaluation_points[:, dimension_idx]):
                eval_block_index = int(eval_point // block_size)
                # Sum frequencies up to the block containing the evaluation point
                smoothed_edf_values[point_idx] += (
                    np.sum(block_counts[:eval_block_index]) / num_samples
                )

        # Normalize the summed probabilities by the number of dimensions to get the average
        return smoothed_edf_values / num_dimensions
