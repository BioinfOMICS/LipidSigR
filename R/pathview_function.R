#' @title DE_pathview_function
#' @description The function is to produce the enrichment plot that assists users to determine whether significant lipid species are enriched in the categories of the selected characteristics. To obtain the enrichment plot, a vector of lipid classes and the pathway information from the database, such as KEGG, need to be provided. A differentially expressed analysis pathview plot and a data frame with key_name, DB, ID, path_id, path_name will be returned as a result.
#' @param lipid_class a vector comprising lipid class, such as PC, PE, TAG
#' @param lipid_gene_path A data frame of pathway information, which can be retrieved and reorganized from KEGG or other databases. It includes the related pathways of lipids, ID, pathway ID, pathway name, gene ID, and gene name.
#' @param path a character string, output path.
#' @param pathway_gene_list A list that comprises the genes included in each pathway for compiling the significant network of the selected lipid class. \emph{NOTE: the information can be retrieved from the KEGG database or others.}
#' @return Return a network plot of pathway analysis and a data frame revealing the details of the plot
#' \enumerate{
#' \item pathview plot save as .png
#' \item A data frame with key_name, DB, ID, path_id, path_name
#' }
#' @export
#' @examples
#' \donttest{
#' data("DE_exp_data")
#' data("DE_lipid_char_table")
#' data("DE_group_info")
#' data("DE.lipid.gene.path")
#' data("DE.pathway.gene.list")
#' exp_data <- DE_exp_data
#' lipid_char_table <- DE_lipid_char_table
#' group_info <- DE_group_info
#' DE.lipid.gene.path <- DE.lipid.gene.path
#' DE.pathway.gene.list <- DE.pathway.gene.list
#' exp_transform <- data_process(exp_data, exclude_var_missing=T, missing_pct_limit=50,
#'                               replace_zero=T, zero2what='min', xmin=0.5, replace_NA=T,
#'                               NA2what='min', ymin=0.5, pct_transform=T, data_transform=T,
#'                               trans_type='log', centering=F,  scaling=F)
#' exp_transform_non_log <- data_process(exp_data, exclude_var_missing=T, missing_pct_limit=50,
#'                                       replace_zero=T, zero2what='min', xmin=0.5,
#'                                       replace_NA=T, NA2what='min', ymin=0.5,
#'                                       pct_transform=T, data_transform=F, trans_type='log',
#'                                       centering=F,  scaling=F)
#' lipid_char_filter <- lipid_char_table %>% filter(feature %in% exp_transform$feature)
#' DE_species_table_sig <- DE_species_2(exp_transform_non_log, data_transform = T,
#'                                      group_info = group_info, paired = F,
#'                                      test = 't.test', adjust_p_method = 'BH',
#'                                      sig_stat = 'p.adj', sig_pvalue = 0.05,
#'                                      sig_FC = 2)$DE_species_table_sig
#' char_var <- colnames(lipid_char_filter)[-1]
#' enrich_char_table <- Enrichment (DE_species_table_sig, lipid_char_table = lipid_char_filter,
#'                                  char_var = char_var[1], sig_pvalue = 0.05)$enrich_char_table
#' sig_enrich_class <- enrich_char_table %>% filter(significant == 'YES') %>%
#'                     distinct(characteristic) %>% .$characteristic
#' dir.create(file.path(getwd(),"pathview_result"), recursive=T)
#' outPath <- file.path(getwd(), "pathview_result")
#' pathview_function(sig_enrich_class, path = getwd(), DE.lipid.gene.path, DE.pathway.gene.list)
#' }
pathview_function <- function(lipid_class, path, lipid_gene_path, pathway_gene_list){

  kegg_cpd <- lipid_gene_path[c(1,2,4,7,8)] %>% dplyr::filter(key_name%in%lipid_class) %>%
    dplyr::filter(DB=='KEGG') %>% unique()
  kegg_cpd <- kegg_cpd[stats::complete.cases(kegg_cpd),]
  kegg_path <- sort(table(stats::na.omit(kegg_cpd$path_id)),decreasing = T)
  pathview_path <- kegg_path %>% names() %>% stringr::str_sub(start = 5)


  cpd <- lipid_gene_path %>% dplyr::filter(key_name%in%lipid_class) %>%
    dplyr::filter(DB=='KEGG') %>% .$ID

  lipid_class_list <- unique(cpd)

  if (length(lipid_class_list)==0) {
    stop('Lipid class name not found')
  }
  pathview_cpd <- rep(1,length(lipid_class_list))

  names(pathview_cpd) <- lipid_class_list

  setwd(path)

  for(a in 1:length(kegg_path)){

    pathway_gene <- which(stringr::str_sub(names(pathway_gene_list),start = 5)==pathview_path[a])

    pathway_gene <- pathway_gene_list[[pathway_gene]] %>% stringr::str_sub(start = 10)

    pathview_gene <- rep(0, length(pathway_gene))
    names(pathview_gene) <- pathway_gene

    tryCatch({
      pathview::pathview(gene.data = pathview_gene, cpd.data=pathview_cpd, pathway.id  = pathview_path[a], species = "hsa",
               out.suffix = "pathview",discrete=list(cpd=T),
               limit = list(cpd = c(0,1)), bins = list(cpd = 1),
               mid = list(cpd = "red"),low = list(cpd = "white", gene='#bfffbf'),
               plot.col.key=F, kegg.native = T,same.layer=T) #kegg.dir='C:/Users/user/Desktop/KEGG/pathview/
    }, error = function(e){message(e)})

  }


  return(kegg_cpd)

}
