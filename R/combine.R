#' Combine replicate differential proteomics experiments.
#'
#' Combine replicate proteomics experiments. Calculates mean PSMs across test and control runs. Calculates weighted z-scores using Stouffer's Z-score method. See README for file formats and example usage on the \href{https://github.com/rachaelcox/enrichr}{enrichr GitHub Repository}.
#'
#' @param rep1 A data frame produced by `enrich_fxn()` corresponding to replicate 1.
#' @param rep2 A data frame produced by `enrich_fxn()` corresponding to replicate 2.
#' @param rep3 Optional. A data frame produced by `enrich_fxn()` corresponding to replicate 3.
#' @param outfile_prefix Optional. File prefix for output `.csv` files. Default = "output_combined."
#' @return A `.csv` with PSM information across biological replicates combined into one table.
#' @export
combine_reps <- function(rep1, rep2, rep3, outfile_prefix){

  rep1 <- rep1 %>%
    mutate(rep = "b1") %>%
    select(-matches("abundance.*"), -matches("number_of.*"),
           -total_PSMs, -PSM_fc)

  rep2 <- rep2 %>%
    mutate(rep = "b2") %>%
    select(-matches("abundance.*"), -matches("number_of.*"),
           -total_PSMs, -PSM_fc)


  if(missing(rep3)){

    combined_df <- rep1 %>%
      bind_rows(rep2) %>%
      pivot_wider(names_from = rep,
                  values_from = c(ctrl_PSMs, exp_PSMs, PSM_zscore, PSM_log2fc))


  } else {

    rep3 <- rep3 %>%
      mutate(rep = "b3") %>%
      select(-matches("abundance.*"), -matches("number_of.*"),
             -total_PSMs, -PSM_fc)

    combined_df <- rep1 %>%
      bind_rows(list(rep2, rep3)) %>%
      pivot_wider(names_from = rep,
                  values_from = c(ctrl_PSMs, exp_PSMs, PSM_zscore, PSM_log2fc))

  }

  ecols <- grep('^exp_PSMs', names(combined_df), value = TRUE)
  ccols <- grep('^ctrl_PSMs', names(combined_df), value = TRUE)
  zcols <- grep('^PSM_zscore', names(combined_df), value = TRUE)

  print("Joining data frames...")
  print(glimpse(combined_df))


  print("Performing cross-replicate calculations...")
  combined_df <- combined_df %>%
    mutate_if(is.numeric, replace_na, replace = 0) %>%
    mutate(joint_zscore = rowSums(select(., zcols))/sqrt(length(zcols))) %>%
    mutate(mean_ctrl_PSMs = rowMeans(select(., ccols))) %>%
    mutate(mean_exp_PSMs = rowMeans(select(., ecols))) %>%
    select(accession, matches("_b\\d"), matches("mean_"),
           joint_zscore, everything()) %>%
    #mutate_if(is.numeric, replace_na, replace = 0) %>%
    arrange(desc(joint_zscore)) %>%
    mutate(pval = pnorm(joint_zscore, lower.tail = FALSE)) %>%
    mutate(fdr_bh = p.adjust(pval, method = "BH", n = length(pval))) %>%
    mutate(conf_90 = case_when(fdr_bh <= 0.10 ~ TRUE,
                               TRUE ~ FALSE),
           conf_95 = case_when(fdr_bh <= 0.05 ~ TRUE,
                               TRUE ~ FALSE),
           conf_99 = case_when(fdr_bh <= 0.01 ~ TRUE,
                               TRUE ~ FALSE))

  write_csv(combined_df, sprintf("%s_combined.csv", outfile_prefix))
  print(combined_df)
  return(combined_df)
}
