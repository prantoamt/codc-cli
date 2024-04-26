import subprocess

import pandas as pd
import click

from analyzer import GeneExpressionAnalyzer
from copula.empirical_copula import EmpiricalCopula


class CustomFormatter(click.HelpFormatter):
    def write_text(self, text):
        if text:
            self.write_paragraph()
            with self.indentation():
                self.write_wrapped_text(text)


@click.group()
@click.pass_context
def cli(ctx):
    """
    \b
    Gene-Coexpress Command Line Interface
    -------------------------------------
    \b
    This CLI provides tools for analyzing gene expression data to identify differential coexpression
    between two phenotypes and for performing gene ontology (GO) enrichment analysis. The tools leverage
    statistical methods to assess changes in gene expression distributions and to link these changes with
    biological annotations.
    \b
    Documentation:
    For a detailed explanation of each command and options, use the '--help' option with the specific
    command. For more detailed theoretical background and methodological details, visit the project wiki at:
    https://github.com/prantoamt/gene-coexpress/wiki
    """
    ctx.help_option_names = ["--help"]
    ctx.max_content_width = 120  # Adjust the width to allow longer lines


@cli.command(
    "copula-diff-coexpress", short_help="Compute network of differential coexpression."
)
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
def copula_diff_coexpress(
    input_file_1, input_file_2, output_path, ties_method, smoothing, ks_stat_method
):
    """
    Compute a network of differential coexpression scores using the
    Kolmogorov-Smirnov distance between empirical copulas. This command
    helps in assessing the similarity in joint gene expression distributions
    between two conditions.
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


@cli.command("go-enrichment", short_help="Run the GO enrichment analysis.")
@click.option(
    "--input_file",
    type=str,
    required=True,
    help="Path to the input TSV file with gene network data.",
)
@click.option(
    "--output_path",
    type=str,
    required=True,
    help="Path where the output plot image will be saved.",
)
@click.option(
    "--threshold",
    type=float,
    default=0.6,
    help="Weight threshold to filter gene pairs.",
)
@click.option(
    "--key_type",
    type=str,
    default="SYMBOL",
    help="The key type of the gene identifiers (default SYMBOL). Other types can be 'ENSEMBL', 'ENTREZID', etc.",
)
@click.option(
    "--ontology", type=str, default="CC", help="Ontology to use: BP, CC, or MF."
)
@click.option(
    "--p_adjust_method", type=str, default="BH", help="Method for adjusting p-values."
)
@click.option(
    "--qvalue_cutoff",
    type=float,
    default=0.05,
    help="Q-value cutoff for significant enrichment.",
)
@click.option(
    "--show_category",
    type=int,
    default=10,
    help="Number of categories to show in the bar plot.",
)
def go_enrichment(
    input_file,
    output_path,
    threshold,
    key_type,
    ontology,
    p_adjust_method,
    qvalue_cutoff,
    show_category,
):
    """
    Perform gene ontology (GO) enrichment analysis on a set of genes to identify
    significant biological terms associated with gene lists derived from differential
    coexpression analysis. This analysis helps in understanding the biological functions
    that are potentially altered under different conditions.
    """

    # Define the path to the R script here to avoid specifying it every time
    R_SCRIPT_PATH = "./go_enrichment_cli.R"
    command = (
        f"Rscript {R_SCRIPT_PATH} "
        f"--input_file={input_file} "
        f"--output_path={output_path} "
        f"--threshold={threshold} "
        f"--key_type={key_type} "
        f"--ontology={ontology} "
        f"--p_adjust_method={p_adjust_method} "
        f"--qvalue_cutoff={qvalue_cutoff} "
        f"--show_category={show_category}"
    )

    try:
        # Start the subprocess and open a pipe to its stdout.
        with subprocess.Popen(
            command,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        ) as proc:
            # Read the output line by line as it becomes available.
            while True:
                output = proc.stdout.readline()
                if output == "" and proc.poll() is not None:
                    break
                if output:
                    print(output.strip(), flush=True)
            # Check for errors
            stderr = proc.stderr.read()
            if stderr:
                print(f"Error output: {stderr}", flush=True)
    except subprocess.CalledProcessError as e:
        click.echo(f"Failed to execute R script. Command tried: {command}")
        click.echo(f"Error output: {e.stderr}")


cli.add_command(copula_diff_coexpress)
cli.add_command(go_enrichment)

if __name__ == "__main__":
    cli()