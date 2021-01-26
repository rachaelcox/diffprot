# diffprot
Annotate, compute, and visualize fold-changes and Z-scores for differential proteomics results.

## Functions
**`enrich()`**: Compute AUC XICs fold-changes (optional) and PSMs fold-changes/Z-scores for `_Proteins.txt` and `_PSMs.txt` assignment files output by Proteome Discoverer 2.3. Annotates results given an annotation table (optional). Annotation table must contain column named `gene_names_primary`. Outputs a .csv file with results.

**`psmplot()`**: Visualize control PSMs versus test PSMs on a log-log plot. Uses dataframe output by `enrich()` as input. Requires annotation column with name ending in `gene_names_primary`. Outputs .png and .pdf files.

**`combine_reps()`**: Combine biological replicates, i.e., combines dataframes output by `enrich()`. Outputs a .csv file with results.

**`combine_exps()`**: Combines experiments, i.e., combines dataframes output by `combine_reps()`.

**`zplot()`**: Visualize PSM-based Z-scores for proteins biological replicates. Uses dataframe output by `combine_reps()` or `combine_exps()` as input. Requires annotation column with name ending in `gene_names_primary`. Outputs .png and .pdf files.

## Data preparation

The `enrich()` function requires a minimum of 3 files to work:

1. A `_Protein.txt` file from Proteome Discoverer.
2. A `_PSMs.txt` file from Proteome Discoverer.
3. A `_meta.txt` file that describes the relationship between control/test samples and their file IDs assigned by Proteome Discoverer.
4. Optional: an annotation file with columns `accession` and `gene_names_primary`.

#### Preparing Proteome Discoverer files

The first two files must be output from Proteome Discoverer in a specific way. **The control and test samples must be processed under a consensus workflow, such that each replicate experiment has its own .pdResult file.** For example, in an APMS experiment with two biological replicates where the test sample contains a GFP-Cadherin fusion and the control sample contains a GFP-only vector:

![Replicate 1](/data_prep/consensus_assignment_b1.PNG)
![Replicate 2](/data_prep/consensus_assignment_b2.PNG)

In the above screenshots, there are 4 technical replicates for both the test (GFPc-Cdh_1a-d) and control (GFP_1a-d) samples for biological replicate #1. All technical replicates, for both the control and test sample, are processed together into a single results file, `cadherin_b1_052920.pdResult`. Biological replicate #2 has 2 technical replicates for the test (GFP-c-Cdh_repeat_1a/b) and control (GFP_1a/b), which are processed together into another result file `cadherin_b1_052920.pdResult`.

Open each .pdResult file and export tab-delimited `_Proteins.txt` and `_PSMs.txt` files as shown in the screenshot below (remember to check "Generate R-Friendly Headers"):

![Export Proteome Discoverer files](/data_prep/export.PNG)

Finally, move these files to your project directory for subsequent analysis in R. For two biological replicates, you should have 4 total files:

![Proteins.txt and PSMs.txt files for 2 biological replicates](/data_prep/files.PNG)

#### Preparing the meta file

Proteome Discoverer will assign a unique file ID for each injection of each sample in your experiments, as shown in the second column of the screenshot below (F19-F30):

![Meta information for each sample](/data_prep/meta_info.PNG)

The `enrich()` function needs to know which file IDs correspond to test and control samples, and which of those correspond to which biological replicate. This is achieved with a **tab-delimited** `meta.txt` file with three columns, `exp_name`, `type`, and `id`. For this example APMS experiment, the information above is converted to the meta file below:

![Meta file format](/data_prep/meta_file.PNG)

I usually make the meta files in Vim, though any text editor will do. You can also make the file in Excel and export as tab-delimited. Whatever you do, it's a good idea to run `cat -T <filename>` to show tabs and make sure Excel hasn't added anything strange.

## Installation
Run the following code:
```r
install.packages("devtools")  # if devtools not installed
devtools::install_github("rachaelcox/diffprot")
library(diffprot)
```
## Example workflow
Compute fold-changes and z-scores for test versus control samples for each biological replicate using `enrich()`. This function outputs a .csv (openable in Excel) with statistics computed for each protein detected in either the test or control cases. In this APMS example, we are only interested in proteins positively enriched in our test case, so we need to set `one_sided = TRUE` to specify we want one-sided statistical tests.
```r
# compute differential protein abundance for biological replicate #1
cadherin_b1 <- enrich(exp_id = "cadherin_b1",
                      meta_file = "cadherin_meta.txt",
                      psm_file = "cadherin_b1_052920_Proteins.txt",
                      pd_file = "cadherin_b1_052920_PSMs.txt",
                      one_sided = TRUE)
```
We often use abstract accessions or grouped identifiers in search databases for assigning peptide mass spectrometry data. Sometimes these accessions are not useful, so this package will also annotate your data given a properly formatted annotation file (UniProt is a great source for this).
```r
# optionally annotate your data
cadherin_b1 <- enrich(exp_id = "cadherin_b1",
                      meta_file = "cadherin_meta.txt",
                      psm_file = "cadherin_b1_052920_Proteins.txt",
                      pd_file = "cadherin_b1_052920_PSMs.txt",
                      one_sided = TRUE,
                      annot_file = "xenla_annots.tab",
                      outfile_name = "cadherin_b1")
```
Ideally you have more than one biological replicate to power your results; if so, you can use `combine_reps()` to combine them, calculate a `joint_zscore` and resulting probability statistics. This function outputs another .csv file.
```r
# compute differential protein abundance for biological replicate #2
cadherin_b2 <- enrich(exp_id = "cadherin_b2",
                      meta_file = "cadherin_meta.txt",
                      psm_file = "cadherin_b2_052920_Proteins.txt",
                      pd_file = "cadherin_b2_052920_PSMs.txt",
                      one_sided = TRUE,
                      annot_file = "xenla_annots.tab",
                      outfile_name = "cadherin_b2")

# combine bio reps
cadherin_all <- combine_reps(rep1 = cadherin_b1,
                             rep2 = cadherin_b2,
                             one_sided = TRUE,
                             outfile_prefix = "cadherin")
```




