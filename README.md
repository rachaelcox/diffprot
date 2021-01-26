# diffprot
Annotate, compute, and visualize fold-changes and Z-scores for differential proteomics results.

### Functions
**`enrich()`**: Compute AUC XICs fold-changes (optional) and PSMs fold-changes/Z-scores for `_Proteins.txt` and `_PSMs.txt` assignment files output by Proteome Discoverer 2.3. Annotates results given an annotation table (optional). Annotation table must contain column named `gene_names_primary`.

**`psmplot()`**: Visualize control PSMs versus test PSMs on a log-log plot. Uses dataframe output by `enrich()` as input. Requires annotation column with name ending in `gene_names_primary`.

**`combine_reps()`**: Combine biological replicates, i.e., combines dataframes output by `enrich()`.

**`combine_exps()`**: Combines experiments, i.e., combines dataframes output by `combine_reps()`.

**`zplot()`**: Visualize PSM-based Z-scores for proteins biological replicates. Uses dataframe output by `combine_reps()` or `combine_exps()` as input. Requires annotation column with name ending in `gene_names_primary`. 

### Data preparation

The `enrich()` function requires a minimum of 3 files to work:

1. A `_Protein.txt` file from Proteome Discoverer.
2. A `_PSMs.txt` file from Proteome Discoverer.
3. A `_meta.txt` file that describes the relationship between control/test samples and their file IDs assigned by Proteome Discoverer.
4. Optional: an annotation file with columns `accession` and `gene_names_primary`.

##### Preparing Proteome Discoverer files

The first two files must be output from Proteome Discoverer in a specific way. **The control and test samples must be processed under a consensus workflow, such that each replicate experiment has its own .pdResult file.** For example, in an APMS experiment with two biological replicates where the test sample contains a GFP-Cadherin fusion and the control sample contains a GFP-only vector:

![Replicate 1](/data_prep_pics/consensus_assignment_b1.png)
![Replicate 2](/data_prep_pics/consensus_assignment_b2.png)

In the above screenshots, there are 4 technical replicates for both the test (GFPc-Cdh_1a-d) and control (GFP_1a-d) samples for biological replicate #1. All technical replicates, for both the control and test sample, are processed together into a single results file, `cadherin_b1_052920.pdResult`. Biological replicate #2 has 2 technical replicates for the test (GFP-c-Cdh_repeat_1a/b) and control (GFP_1a/b), which are processed together into another result file `cadherin_b1_052920.pdResult`.

Open each .pdResult file and export tab-delimited `_Proteins.txt` and `_PSMs.txt` files as shown in the screenshot below (remember to check "Generate R-Friendly Headers"):

![Export Proteome Discoverer files](/data_prep_pics/export.png)

Finally, move these files to your project directory for subsequent analysis in R. For two biological replicates, you should have 4 total files:

![Proteins.txt and PSMs.txt files for 2 biological replicates](/data_prep_pics/files.png)

##### Preparing the meta file

Proteome Discoverer will assign a unique file ID for each injection of each sample in your experiments, as shown in the second column of the screenshot below (F19-F30):

![Meta information for each sample](/data_prep_pics/meta_info.png)

The `enrich()` function needs to know which file IDs correspond to test and control samples, and which of those correspond to which biological replicate. This is achieved with a **tab-delimited** `meta.txt` file with three columns, `exp_name`, `type`, and `id`. For this example APMS experiment, the information above is converted to the meta file below:

### Installation
Run the following code:
``` r
install.packages("devtools")  # if devtools not installed
devtools::install_github("rachaelcox/diffprot")
library(diffprot)
```
### Example workflow
