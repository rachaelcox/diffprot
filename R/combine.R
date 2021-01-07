#' Combine replicates for a proteomics experiment.
#'
#' Combine replicate proteomics experiments. Calculates mean PSMs across test and control runs. Calculates weighted z-scores using Stouffer's Z-score method. Calculates p-values and FDR on the weighted Z-score using the Benjamini-Hochberg procedure. See README for file formats and example usage on the \href{https://github.com/rachaelcox/enrichr}{enrichr GitHub Repository}.
#'
#' @param rep1 A data frame produced by `enrich_fxn()` corresponding to replicate 1.
#' @param rep2 A data frame produced by `enrich_fxn()` corresponding to replicate 2.
#' @param rep3 Optional. A data frame produced by `enrich_fxn()` corresponding to replicate 3.
#' @param one_sided Specify one-sided or two-sided statistical test; default = `FALSE` (two-sided).
#' @param outfile_prefix Optional. File prefix for output `.csv` files. Default = "reps_combined".
#' @return A data frame and `.csv` with PSM information across biological replicates combined into one table.
#' @export
combine_reps <- function(rep1, rep2, rep3, one_sided = FALSE, outfile_prefix){

  # combine replicates into one dataframe
  rep1 <- rep1 %>%
    dplyr::mutate(rep = "b1") %>%
    dplyr::select(-matches("abundance.*"), -matches("number_of.*"),
           -total_PSMs, -PSM_fc, -pval, -fdr_bh)

  rep2 <- rep2 %>%
    dplyr::mutate(rep = "b2") %>%
    dplyr::select(-matches("abundance.*"), -matches("number_of.*"),
           -total_PSMs, -PSM_fc, -pval, -fdr_bh)

  if(missing(outfile_prefix)){
    outfile_prefix = "reps_combined"
  }


  if(missing(rep3)){

    combined_df <- rep1 %>%
      bind_rows(rep2) %>%
      pivot_wider(names_from = rep,
                  values_from = c(ctrl_PSMs, exp_PSMs, PSM_zscore, PSM_log2fc))


  } else {

    rep3 <- rep3 %>%
      dplyr::mutate(rep = "b3") %>%
      dplyr::select(-matches("abundance.*"), -matches("number_of.*"),
             -total_PSMs, -PSM_fc, -pval, -fdr_bh)

    combined_df <- rep1 %>%
      bind_rows(list(rep2, rep3)) %>%
      pivot_wider(names_from = rep,
                  values_from = c(ctrl_PSMs, exp_PSMs, PSM_zscore, PSM_log2fc))

  }

  # define experiment, control, and zscore columns for each rep
  ecols <- grep('^exp_PSMs', names(combined_df), value = TRUE)
  ccols <- grep('^ctrl_PSMs', names(combined_df), value = TRUE)
  zcols <- grep('^PSM_zscore', names(combined_df), value = TRUE)

  print("Joining data frames...")
  print(glimpse(combined_df))

  # compute cross-replicate statistics
  print("Performing cross-replicate calculations...")
  combined_df <- combined_df %>%
    dplyr::mutate_if(is.numeric, replace_na, replace = 0) %>%
    dplyr::mutate(joint_zscore = rowSums(select(., zcols))/sqrt(length(zcols))) %>%
    dplyr::mutate(mean_ctrl_PSMs = rowMeans(select(., ccols))) %>%
    dplyr::mutate(mean_exp_PSMs = rowMeans(select(., ecols))) %>%
    dplyr::select(accession, matches("_b\\d"), matches("mean_"),
           joint_zscore, everything()) %>%
    dplyr::arrange(desc(joint_zscore))

  # calculate probabilities
  if(!one_sided){
    combined_df <- combined_df %>%
      dplyr::mutate(pval = pnorm(joint_zscore)) %>%
      dplyr::mutate(fdr_bh = p.adjust(pval, method = "BH", n = length(pval)))
  } else {
    combined_df <- combined_df %>%
      dplyr::mutate(pval = pnorm(joint_zscoree, lower.tail = FALSE)) %>%
      dplyr::mutate(fdr_bh = p.adjust(pval, method = "BH", n = length(pval)))
  }

  # label cutoffs
  combined_df <- combined_df %>%
    dplyr::mutate(conf_90 = dplyr::case_when(fdr_bh <= 0.10 ~ TRUE,
                                             TRUE ~ FALSE),
                  conf_95 = dplyr::case_when(fdr_bh <= 0.05 ~ TRUE,
                                             TRUE ~ FALSE),
                  conf_99 = dplyr::case_when(fdr_bh <= 0.01 ~ TRUE,
                                             TRUE ~ FALSE))

  readr::write_csv(combined_df, sprintf("%s_combined.csv", outfile_prefix))
  print(combined_df)
  return(combined_df)
}

#' Combine and compare proteomics experiments.
#'
#' Combine different proteomics experiments, i.e., an HA-tag pulldown and a TurboID-biotin pulldown. Compares joint z-scores and FDR across each experiment. See README for file formats and example usage on the \href{https://github.com/rachaelcox/enrichr}{enrichr GitHub Repository}.
#'
#' @param exp1 A data frame produced by `combine_reps()` corresponding to experiment 1.
#' @param exp2 A data frame produced by `combine_reps()` corresponding to experiment 2.
#' @param exp1_id Prefix that will annotate experiment 1 columns.
#' @param exp2_id Prefix that will annotate experiment 2 columns.
#' @param outfile_prefix Optional. File prefix for output `.csv` files. Default = "exps_combined".
#' @return A data frame and `.csv` with PSM information across biological replicates combined into one table.
#' @export
combine_exps <- function(exp1, exp2, exp1_id, exp2_id, outfile_prefix){

  cols <- grep('joint_zscore|fdr_bh', names(exp1), value = TRUE)

  if(length(cols) == 0){
    stop("It doesn't look like you're using the correct dataframes as input. \\
         Are you looking for `combine_reps()`? This function uses data frames output by `combine_reps()`. \\
         See `?combine_exps` and `?combine_reps` for more details.")
  }

  if(missing(outfile_prefix)){
    outfile_prefix = "exps_combined"
  }

  exp1_sig <- paste0("Significant in ", toupper(exp1_id))
  exp2_sig <- paste0("Significant in ", toupper(exp2_id))
  both_sig <- "Significant in both"
  no_sig <- "Below threshold in both"

  exp1_id <- paste0("_", exp1_id)
  exp2_id <- paste0("_", exp2_id)

  df1 <- exp1 %>%
    select(accession, joint_zscore, fdr_bh, gene_names_primary) %>%
    rename_with(.cols = all_of(cols), .fn = ~paste0(., exp1_id))

  df2 <- exp2 %>%
    select(accession, joint_zscore, fdr_bh, gene_names_primary) %>%
    rename_with(.cols = all_of(cols), .fn = ~paste0(., exp2_id))

  df1_fdr <- colnames(df1)[3]
  df1_z <- colnames(df1)[2]
  df2_fdr <- colnames(df2)[3]
  df2_z <- colnames(df2)[2]

  combined_df <- full_join(df1, df2, by = c("accession", "gene_names_primary")) %>%
    mutate(conf_90 = dplyr::case_when(get(df1_fdr) <= 0.10 & get(df2_fdr) <= 0.10 ~ both_sig,
                               get(df1_fdr) <= 0.10 & get(df2_fdr) > 0.10 ~ exp1_sig,
                               get(df1_fdr) > 0.10 & get(df2_fdr) <= 0.10 ~ exp2_sig,
                               get(df1_fdr) > 0.10 & get(df2_fdr) > 0.10 ~ no_sig),
           conf_95 = dplyr::case_when(get(df1_fdr) <= 0.05 & get(df2_fdr) <= 0.05 ~ both_sig,
                               get(df1_fdr) <= 0.05 & get(df2_fdr) > 0.05 ~ exp1_sig,
                               get(df1_fdr) > 0.05 & get(df2_fdr) <= 0.05 ~ exp2_sig,
                               get(df1_fdr) > 0.05 & get(df2_fdr) > 0.05 ~ no_sig),
           conf_99 = dplyr::case_when(get(df1_fdr) <= 0.01 & get(df2_fdr) <= 0.01 ~ both_sig,
                               get(df1_fdr) <= 0.01 & get(df2_fdr) > 0.01 ~ exp1_sig,
                               get(df1_fdr) > 0.01 & get(df2_fdr) <= 0.01 ~ exp2_sig,
                               get(df1_fdr) > 0.01 & get(df2_fdr) > 0.01 ~ no_sig)) %>%
    mutate(abs_zscore = abs(get(df1_z)) + abs(get(df2_z))) %>%
    arrange(desc(abs_zscore))

  print(combined_df)
  return(combined_df)
  write_csv(combined_df, sprintf("%s_combined.csv", outfile_prefix))
  }
