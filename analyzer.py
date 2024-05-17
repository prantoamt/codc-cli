from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import shared_memory
import os

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

    def compute_pairs(
        self,
        indices,
        shm_name_data1,
        shm_name_data2,
        gene_names,
        ties_method,
        smoothing,
        ks_stat_method,
        data1_shape,
        data2_shape,
        dtype,
    ):
        # Access shared memory
        existing_shm_data1 = shared_memory.SharedMemory(name=shm_name_data1)
        existing_shm_data2 = shared_memory.SharedMemory(name=shm_name_data2)

        # Create numpy arrays from the buffers of the shared memory
        np_data1 = np.ndarray(data1_shape, dtype=dtype, buffer=existing_shm_data1.buf)
        np_data2 = np.ndarray(data2_shape, dtype=dtype, buffer=existing_shm_data2.buf)

        results = []
        for i, j in indices:
            try:
                gene_pair_data1 = np.vstack((np_data1[:, i], np_data1[:, j])).T
                gene_pair_data2 = np.vstack((np_data2[:, i], np_data2[:, j])).T
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
                results.append(
                    {
                        "Target": gene_names[j],
                        "Regulator": gene_names[i],
                        "Condition": "Diff Co-Exp of both Condition",
                        "Weight": ks_stat,
                    }
                )
            except Exception as e:
                print(f"Error processing pair ({i}, {j}): {e}")

        # Clean up shared memory
        existing_shm_data1.close()
        existing_shm_data2.close()

        return results

    def compute_dc_copula_network_parallel(
        self,
        df1,
        df2,
        ties_method="average",
        smoothing="none",
        ks_stat_method="asymp",
        batch_size=100,
    ):
        gene_names = df1.iloc[:, 0].values
        assert np.array_equal(
            gene_names, df2.iloc[:, 0].values
        ), "Gene lists must match!"

        data1 = df1.iloc[:, 1:].values.T
        data2 = df2.iloc[:, 1:].values.T

        # Create shared memory
        shm_data1 = shared_memory.SharedMemory(create=True, size=data1.nbytes)
        shm_data2 = shared_memory.SharedMemory(create=True, size=data2.nbytes)
        # Create numpy arrays on the buffer of the shared memory
        np_data1 = np.ndarray(data1.shape, dtype=data1.dtype, buffer=shm_data1.buf)
        np_data2 = np.ndarray(data2.shape, dtype=data2.dtype, buffer=shm_data2.buf)
        np.copyto(np_data1, data1)
        np.copyto(np_data2, data2)

        n_genes = min(len(data1[0]), len(data2[0]))
        pairs = [(i, j) for i in range(n_genes - 1) for j in range(i + 1, n_genes)]
        batches = [pairs[i : i + batch_size] for i in range(0, len(pairs), batch_size)]

        results = pd.DataFrame()

        # Print dataset summary
        print(f"Starting DC Copula coexpression calculation:")
        print(f"-------------------------")
        print(f" - Number of gene pairs to be analyzed: {n_genes * (n_genes - 1) // 2}")
        print(f" - Batch size: {batch_size}")
        print(f" - Number of batches: {len(batches)}")
        print(f" - Ties method: {ties_method}")
        print(f" - Smoothing technique: {smoothing}")
        print(f" - KS statistic mode: {ks_stat_method}")
        print(f"-------------------------")

        with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
            futures = []

            print("\nQueueing tasks...")
            submission_progress = tqdm(
                total=len(batches), desc="Queueing gene pairs", unit="batch"
            )

            for batch in batches:
                futures.append(
                    executor.submit(
                        self.compute_pairs,
                        batch,
                        shm_data1.name,
                        shm_data2.name,
                        gene_names,
                        ties_method,
                        smoothing,
                        ks_stat_method,
                        data1.shape,
                        data2.shape,
                        data1.dtype,
                    )
                )
                submission_progress.update(1)
            submission_progress.close()

            print("\nProcessing gene pairs...")
            completion_progress = tqdm(
                as_completed(futures),
                total=len(futures),
                desc="Computing distances",
                unit="batch",
            )

            for future in completion_progress:
                batch_results = future.result()
                results = pd.concat(
                    [results, pd.DataFrame(batch_results)], ignore_index=True
                )
                completion_progress.update(1)
            completion_progress.close()

        # Clean up shared memory
        shm_data1.close()
        shm_data1.unlink()
        shm_data2.close()
        shm_data2.unlink()

        return results
