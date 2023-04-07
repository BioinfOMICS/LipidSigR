library(magrittr)
library(dplyr)
library(pathview)

# Load data and create example data sets
data(DE_exp_data)
data(DE_lipid_char_table)
data(DE_group_info)
data("DE_lipid_gene_path")
data("DE_pathway_gene_list")
group_info <- DE_group_info
DE.lipid.gene.path <- DE_lipid_gene_path
DE.pathway.gene.list <- DE_pathway_gene_list

# Start tests
test_that("pathview_function works", {
  exp_transform <- data_process(DE_exp_data, exclude_var_missing=TRUE,
      missing_pct_limit=50, replace_zero=TRUE, zero2what='min', xmin=0.5,
      replace_NA=TRUE, NA2what='min', ymin=0.5, pct_transform=TRUE,
      data_transform=TRUE, trans_type='log', centering=FALSE, scaling=FALSE)

  exp_transform_non_log <- data_process(DE_exp_data, exclude_var_missing=TRUE,
      missing_pct_limit=50, replace_zero=TRUE, zero2what='min', xmin=0.5,
      replace_NA=TRUE, NA2what='min', ymin=0.5, pct_transform=TRUE,
      data_transform=FALSE, trans_type='log', centering=FALSE, scaling=FALSE)
  lipid_char_filter <- DE_lipid_char_table %>%
    filter(feature %in% exp_transform$feature)
  DE_species_table_sig <- DE_species_2(exp_transform_non_log,
      data_transform = TRUE, group_info = group_info, paired = FALSE,
      test = 't.test', adjust_p_method = 'BH', sig_stat = 'p.adj',
      sig_pvalue = 0.05, sig_FC = 2, plotly=TRUE)$DE_species_table_sig
  char_var <- colnames(lipid_char_filter)[-1]
  enrich_char_table <- Enrichment (DE_species_table_sig,
      lipid_char_table = lipid_char_filter, char_var = char_var[1],
      sig_pvalue = 0.05, plotly=TRUE)$enrich_char_table
  sig_enrich_class <- enrich_char_table %>% filter(significant == 'YES') %>%
    distinct(characteristic) %>% .$characteristic
  dir.create(file.path(getwd(),"pathview_result"), recursive=TRUE)
  outPath <- file.path(getwd(), "pathview_result")
  pathview_function(sig_enrich_class, path = getwd(), DE.lipid.gene.path,
      DE.pathway.gene.list)
})
