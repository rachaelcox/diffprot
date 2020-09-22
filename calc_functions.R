format_psm_file <- function(exp_id, meta_file, psm_file){
  
  meta <- read_tsv(meta_file) %>% 
    filter(exp_name == exp_id)
  
  psm_df <- read_tsv(psm_file, col_types = list("Quan Info" = col_character())) %>%
  clean_names() %>%
    mutate(file_id = tolower(file_id)) %>% 
  separate_rows(protein_accessions, sep='; ') %>%
    filter(!grepl("CONTAMINANT", protein_accessions)) %>%
  group_by(file_id, protein_accessions) %>%
    tally() %>%
    ungroup() %>%
  left_join(meta, by=c("file_id" = "id")) %>%
  select(-exp_name, -file_id)
  print(psm_df)

  
  psm_df <- psm_df %>%
  group_by(protein_accessions, type) %>%
    summarise(PSMs = sum(n)) %>%
    ungroup() %>%
  pivot_wider(names_from = type, values_from = PSMs, values_fill = list(PSMs = 0)) %>%
  rename(ctrl_PSMs = ctrl, exp_PSMs = exp)
  
  print(psm_df)
  return(psm_df)
  
}

format_pd_file <- function(exp_id, meta_file, pd_file, psm_df){
  
    meta <- read_tsv(meta_file) %>% 
      filter(exp_name == exp_id)
    
    pd_df <- read_tsv(pd_file) %>%
      clean_names() %>%
      rename(number_PSMs = number_of_ps_ms) %>%
      filter(!grepl("CONTAMINANT", accession)) %>%
      select(accession, description, number_PSMs,
             number_of_unique_peptides, number_of_peptides,
             matches("abundance.*")) %>%
      rename(eggnogID=description) %>%
      left_join(psm_df, by=c("accession" = "protein_accessions")) %>%
      select(accession, eggnogID, ctrl_PSMs, exp_PSMs, everything())

    names(pd_df) <- names(pd_df) %>%
      gsub("abundance_", "", .) %>%
      gsub("_sample", "", .) %>%
      gsub("_control", "", .) %>%
      str_replace("^f", meta$type[match(.,meta$id)])
    
    return(pd_df)
}


impute_fxn <- function(pd_df){
      
    print(pd_df)
  
    pd_impute_df <- pd_df %>%
      mutate(comb_ctrl = rowSums(pd_df[, grep("^ctrl\\d.*", names(pd_df))])) %>% 
      mutate(comb_exp = rowSums(pd_df[, grep("^exp\\d.*", names(pd_df))])) %>%
      mutate(min_col = pmin(comb_ctrl, comb_exp, na.rm=TRUE))
    
    pd_impute_df[c("comb_ctrl","comb_exp")][is.na(pd_impute_df[c("comb_ctrl","comb_exp")])] = min(pd_impute_df$min_col, na.rm = TRUE)
    
    return(pd_impute_df)
}

enrichment_fxn <- function(pd_impute_df){
  
      pd_enriched <- pd_impute_df %>%
        # XIC fold change calcuations
        mutate(mean_abundance = (comb_exp+comb_ctrl)/2) %>%
        mutate(abundance_fc = case_when(comb_ctrl > comb_exp ~ -comb_ctrl/comb_exp,
                                       comb_ctrl < comb_exp ~ comb_exp/comb_ctrl)) %>%
        mutate(abundance_log2fc = log2((comb_exp/comb_ctrl))) %>%
        # PSM fold change calculations
        mutate(expPSMs_pseudo = exp_PSMs + 1,
               ctrlPSMs_pseudo = ctrl_PSMs + 1) %>%
        mutate(F0ctrl_pseudo = ctrlPSMs_pseudo / sum(ctrlPSMs_pseudo)) %>%
        mutate(F0exp_pseudo = expPSMs_pseudo / sum(expPSMs_pseudo)) %>%
        mutate(PSM_fc = case_when(F0ctrl_pseudo > F0exp_pseudo ~ -F0ctrl_pseudo/F0exp_pseudo,
                                  F0ctrl_pseudo < F0exp_pseudo ~ F0exp_pseudo/F0ctrl_pseudo)) %>%
        mutate(PSM_log2fc = log2((F0exp_pseudo/F0ctrl_pseudo))) %>%
        # PSM z-score calculations
        mutate(F0_ctrl = ctrl_PSMs / (sum(ctrl_PSMs))) %>%
        mutate(F0_exp = exp_PSMs / (sum(exp_PSMs))) %>%
        mutate(F1 = (ctrl_PSMs + exp_PSMs) / (sum(ctrl_PSMs) + sum(exp_PSMs))) %>%
        mutate(PSM_zscore = 
                 (F0_exp - F0_ctrl)/
                 sqrt((
                   (F1 * (1-F1) / sum(exp_PSMs))) + 
                     (F1 * (1-F1) / sum(exp_PSMs)))) %>%
        # format final dataframe
        arrange(desc(PSM_zscore)) %>%
        rename(total_PSMs = number_PSMs, ctrl_abundance = comb_ctrl,
               exp_abundance = comb_exp) %>% 
        select(accession, eggnogID, number_of_unique_peptides, number_of_peptides, 
               total_PSMs, ctrl_abundance, exp_abundance, mean_abundance, abundance_fc,
               abundance_log2fc, ctrl_PSMs, exp_PSMs, PSM_fc, PSM_log2fc, PSM_zscore)
      
      print(pd_enriched)
      
      return(pd_enriched)
}

abundance_fxn <- function(exp_id, meta_file, psm_file, pd_file, annot_file, outfile_name = "output.csv"){

    if(str_sub(outfile_name, -3, -1) != "csv"){
      outfile_name = paste0(outfile_name, ".csv")
    }
    
    # load in and format Proteome Discoverer files, meta file
    psm_df <- format_psm_file(exp_id, meta_file, psm_file)
    
    pd_df <- format_pd_file(exp_id, meta_file, pd_file, psm_df)

    # impute missing values
    pd_impute_df <- impute_fxn(pd_df)

    # sum technical replicates, calculate fold change and log ratios
    pd_enriched <- enrichment_fxn(pd_impute_df)

    # join on annotations if specified
    if(!missing(annot_file)){
       annot_table <- read_tsv(annot_file) %>%
          rename(eggnogID=ID)

       annot_table <- annot_table %>% 
         filter(!is.na(eggnogID)) %>%
         column_to_rownames(var = "eggnogID")
       
       # manual annotation cleaning -- remove for final package?
       annot_table["ENOG411D4D3","HUMAN_GeneNames_Primary"] <- "C16orf71"
       annot_table["ENOG411DGHA","HUMAN_GeneNames_Primary"] <- "DNAI2"
       annot_table["ENOG411D33H","HUMAN_GeneNames_Primary"] <- "H2B Histones"
       annot_table["ENOG411CPDK","HUMAN_GeneNames_Primary"] <- "Beta Tubulins"
       annot_table["ENOG411D33H","HUMAN_GeneNames_Primary"] <- "H2B Histones"
       
       annot_table <- annot_table %>%
         rownames_to_column(var = "eggnogID")
       
       pd_enriched <- pd_enriched %>%
        left_join(annot_table, by = "eggnogID") %>%
        mutate(XENLA_XenBase_GeneNames = gsub('NA, ', '', XENLA_XenBase_GeneNames))
       
    }

    write_csv(pd_enriched, outfile_name)
    
    return(pd_enriched)
    
}
