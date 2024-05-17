# Methodology

The methodology adopted for identifying differential gene co-expression using a copula-based approach is systematically structured into computational steps that transform raw gene expression data into a quantifiable measure of differential co-expression. Here, I detail each computational step and the mathematical underpinnings that guide the analyses.

The source code of R's copula package, which do this calculation on the R implementation can be found [here](https://rdrr.io/cran/copula/src/R/empCopula.R)

## Conceptual Framework

### Pseudo-Observations

Converting raw data samples into pseudo-observations is a crucial step in copula-based analysis. This process involves ranking the data within each variable and normalizing these ranks to the unit interval $[0, 1)$, which prepares them for further analysis using empirical copulas. The implementations use the equation 1 for this purpose. There is no direct method in R which implements this, but they did it indirectly. The python implementation `pseudo_observations()` function can be found under `copula/empirical_copula.py` file.

$$
p_{ij} = \frac{\text{rank}(x_{ij})}{n + 1} [equation 1]
$$

Here, $p_{ij}$ is the pseudo-observation corresponding to $x_{ij}$. The function $\text{rank}(x_{ij})$ assigns a rank to each data point based on its value relative to other samples in the same variable, with adjustments based on the chosen ties method (e.g., 'average', 'min', 'max', 'dense', 'ordinal').

### Empirical Copula 

A copula is a statistical tool used to describe the dependency between random variables. The empirical copula specifically captures the joint distribution of variables by transforming their margins to uniform distributions on the interval $[0, 1]$.

$$
C(u, v) = \text{Prob}(U \leq u, V \leq v) [equation 2]
$$

where $U$ and $V$ are the uniform transformations (pseudo-observations) of the original variables $X$ and $Y$, representing gene expressions.

### Comparing Empirical Copulas with the KS Test

The Kolmogorov-Smirnov test is employed to compare the empirical copulas of a gene pair across two conditions using the equation, quantifying the maximum distance between their cumulative distribution functions.

$$
D = \sup_x |F_{1,n}(x) - F_{2,n}(x)| [equation 3]
$$

where $F_{1,n}$ and $F_{2,n}$ are the empirical distribution functions of the gene pair under the first and second conditions, respectively.

### Copula-based Differential Co-Expression

The equation 4 is proposed by the authors to calculate the final distance, which utilizes the previous equations in order to calculate the differential co-expression of a particular gene pair.

$$
\begin{split}
    DC_{Copula}(g_i, g_j) &= D(C(g_i^{condition 1}, g_j^{condition 1}), \\
                           &\quad C(g_i^{condition 2}, g_j^{condition 2})) [equation 4]
\end{split}
$$

where $g_i$ and $g_j$ are genes being compared, and $DC_{Copula}(g_i, g_j)$ is the differential coexpression score for the gene pair $(g_i, g_j)$.

## Computational Steps

#### Copula Based Differential Co-Expression

**Requirements:**
- Input file 1
- Input file 2
- Output path
- Ensure the correctness of input format and existence of the output directory

**Procedure:**
1. **Data Preparation:**
   - Verify that both matrices contain identical genes in the same order.

2. **Pairing and Transforming Gene Data:**
   - For each gene pair \(i\) and \(j\), construct two matrices for each condition by horizontally stacking their expression profiles.

3. **Generating Pseudo-Observations:**
   - Convert expression levels to a scale of \[0, 1\] by ranking the expression values and scaling these ranks using the specified equation 1.

4. **Empirical Copula Construction:**
   - Establish the joint distribution function using equation 2, for each gene pair in each condition based on the uniform pseudo-observations.

5. **Differential Coexpression Calculation:**
   - Apply the KS test to the empirical copulas from both conditions to determine the maximum divergence in their distribution behaviors. Use equation 3 on this purpose.

6. **Compiling and Reporting Results:**
   - Output a `network.tsv` file detailing the differential coexpression network, specifying target, regulator, condition, and weight for each gene pair.

**Return:**
- `network.tsv`
