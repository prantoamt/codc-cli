import time
import pandas as pd
import csv
from analyzer import GeneExpressionAnalyzer
from copula.empirical_copula import EmpiricalCopula


# Configuration
input_file_1 = (
    "/Users/pranto/Desktop/BIONet/BRCA_normal.tsv"  # Replace with your input file path
)
input_file_2 = (
    "/Users/pranto/Desktop/BIONet/BRCA_normal.tsv"  # Replace with your input file path
)
ties_method = "average"  # Replace with your parameter
smoothing = "none"  # Replace with your parameter
ks_stat_method = "asymp"  # Replace with your parameter
output_dir = "/Users/pranto/Desktop/BIONet"  # Replace with your output directory
timings_csv = f"{output_dir}/python_execution_times.csv"

# Loading data from TSV files once
df1 = pd.read_csv(input_file_1, delimiter="\t")
df2 = pd.read_csv(input_file_2, delimiter="\t")

# Initializing the EmpiricalCopula and GeneExpressionAnalyzer instances
empirical_copula = EmpiricalCopula()
analyzer = GeneExpressionAnalyzer(empirical_copula=empirical_copula)

# Measuring execution times
execution_times = []
for _ in range(10):
    start_time = time.time()
    network_df = analyzer.compute_dc_copula_network(
        df1,
        df2,
        ties_method=ties_method,
        smoothing=smoothing,
        ks_stat_method=ks_stat_method,
    )
    elapsed_time = time.time() - start_time
    execution_times.append(elapsed_time)
    print(f"Execution {_ + 1}: {elapsed_time:.4f} seconds")

# Storing execution times in a CSV file
with open(timings_csv, "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["Execution", "Time"])
    for idx, exec_time in enumerate(execution_times, start=1):
        writer.writerow([idx, exec_time])

print(f"Saved execution times to {timings_csv}")
