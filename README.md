# diffprot
Annotate, compute, and visualize fold-changes and Z-scores for differential proteomics results.

### Functions
**`enrich()`**: Compute AUC XICs fold-changes (optional) and PSMs fold-changes/Z-scores for `_Proteins.txt` and `_PSMs.txt` assignment files output by Proteome Discoverer 2.3. Annotates results given an annotation table (optional). Annotation table must contain column with name ending in `gene_names_primary`.

**`psmplot()`**: Visualize control PSMs versus test PSMs on a log-log plot. Uses dataframe output by `enrich()` as input. Requires annotation column with name ending in `gene_names_primary`.

**`combine_reps()`**: Combine biological replicates, i.e., combines dataframes output by `enrich()`.

**`combine_exps()`**: Combines experiments, i.e., combines dataframes output by `combine_reps()`.

**`zplot()`**: Visualize PSM-based Z-scores for proteins biological replicates. Uses dataframe output by `combine_reps()` or `combine_exps()` as input. Requires annotation column with name ending in `gene_names_primary`. 

### Installation
Run the following code:
``` r
install.packages("devtools")  # if devtools not installed
devtools::install_github("rachaelcox/diffprot")
library(diffprot)
```
### Example workflow
