format_psm_file <- function(exp_id, meta_file, psm_file){

  meta <- readr::read_tsv(meta_file) %>%
    dplyr::filter(exp_name == exp_id)

  psm_df <- readr::read_tsv(psm_file, col_types = list("Quan Info" = col_character())) %>%
    janitor::clean_names() %>%
    dplyr::mutate(file_id = tolower(file_id)) %>%
    tidyr::separate_rows(protein_accessions, sep='; ') %>%
    dplyr::filter(!grepl("CONTAMINANT", protein_accessions)) %>%
    dplyr::group_by(file_id, protein_accessions) %>%
    dplyr::tally() %>%
    dplyr::ungroup() %>%
    dplyr::left_join(meta, by=c("file_id" = "id")) %>%
    dplyr::select(-exp_name, -file_id)

  print(psm_df)

  psm_df <- psm_df %>%
    dplyr::group_by(protein_accessions, type) %>%
    dplyr::summarise(PSMs = sum(n)) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_wider(names_from = type, values_from = PSMs,
                       values_fill = list(PSMs = 0)) %>%
    dplyr::rename(ctrl_PSMs = ctrl, exp_PSMs = exp)

  print(psm_df)
  return(psm_df)

}

format_pd_file <- function(exp_id, meta_file, pd_file, psm_df){

  meta <- readr::read_tsv(meta_file) %>%
    dplyr::filter(exp_name == exp_id)

  pd_df <- readr::read_tsv(pd_file) %>%
    janitor::clean_names() %>%
    dplyr::rename(number_PSMs = number_of_ps_ms) %>%
    dplyr::filter(!grepl("CONTAMINANT", accession)) %>%
    dplyr::select(accession, description, number_PSMs,
           number_of_unique_peptides, number_of_peptides,
           matches("abundance.*")) %>%
    dplyr::left_join(psm_df, by=c("accession" = "protein_accessions")) %>%
    dplyr::select(accession, description, ctrl_PSMs, exp_PSMs, everything())

  names(pd_df) <- names(pd_df) %>%
    gsub("abundance_", "", .) %>%
    gsub("_sample", "", .) %>%
    gsub("_control", "", .) %>%
    stringr::str_replace("^f", meta$type[match(.,meta$id)])

  return(pd_df)
}
