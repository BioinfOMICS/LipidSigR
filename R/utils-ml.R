###### ML_process
.ml_process <- function(processed_se, char, transform, type=c("normal", "transpose")){
    species_abundance <- .extract_df(processed_se, type="abundance")
    .check_imputation(species_abundance)
    if(isTRUE(.check_nor_negative(species_abundance[-1]))) {
        stop("Detect negative values in the input abundance data. We recommended choosing a different normalization method during data processing.")
    }
    group_info <- .extract_df(processed_se, type="group")
    if(any(char!="none")){
        char_data_all <- convert_sp2char(processed_se, transform='none')
        char_abundance_all <- .extract_df(char_data_all, type="abundance")
        char_abundance <-  char_abundance_all %>%
            dplyr::filter(
                sapply(paste0(char, "|"), function(prefix) startsWith(feature, prefix)) %>% rowSums() > 0)
        abundance_bind <- rbind(char_abundance, species_abundance)
    } else {
        abundance_bind <- species_abundance
    }
    ## transform data (ML_data[[1]])
    abundance_trans <- .transform(abundance_bind, transform=transform)
    if (type=="normal") {
        return(ML_data=abundance_trans)
    } else if (type=="transpose") {
        ## (ML_data[[2]]) : transpose, replace ":" as "__", " " as "_", "|" as "___"
        abundance_t <- abundance_trans %>%
            dplyr::mutate(feature = .lipid_name_replace(abundance_trans$feature, type="replace")) %>%
            # dplyr::mutate(feature = gsub(':','__',gsub(' ','_',feature))) %>%
            # dplyr::mutate(feature = gsub('\\|','___', feature)) %>%
            tibble::column_to_rownames(var="feature") %>%
            t() %>% as.data.frame() %>% tibble::rownames_to_column("sample_name")
        ## bind group by group_info
        abundance_with_group <- abundance_t %>%
            dplyr::left_join(group_info, by='sample_name') %>%
            #dplyr::select(-sample_name) %>%
            tibble::column_to_rownames(var="sample_name") %>%
            dplyr::select(group, dplyr::everything())
        return(ML_data=abundance_with_group)
    } else {
        stop("type must be 'normal' or 'transpose'.")
    }
}

### plot_ROC, plot_PR functions
# Function to calculate AUC statistics
# set feature=NULL for plotting tab
.auc_stats <- function(data, feature) {
    if (is.null(feature)){
        cv_auc <- data %>% unique() %>% dplyr::pull(auc)
    } else {
        cv_auc <- data %>%
            dplyr::filter(feature_num == feature) %>% dplyr::distinct(auc) %>%
            dplyr::pull(auc)
    }
    t_test_res <- stats::t.test(cv_auc, mu=0.5)
    list(
        lower95=round(t_test_res$conf.int[1], 3),
        upper95=round(t_test_res$conf.int[2], 3),
        AUC_pvalue=formatC(t_test_res$p.value, digits = 2, format = "e"),
        mean_AUC=round(t_test_res$estimate, 3)
    )
}
# Function to create AUC label
.auc_label <- function(stats) {
    stringr::str_c(
        'AUC=', as.character(stats$mean_AUC), ' (',as.character(stats$lower95), '-',
        as.character(stats$upper95),')')
    #sprintf("AUC=%.3f (%.3f-%.3f)", stats$mean_AUC, stats$lower95, stats$upper95)
}

# Function to pad feature number
.label_pad <- function(feature_num) {
    padding <- dplyr::case_when(
        nchar(feature_num) == 1 ~ 5,
        nchar(feature_num) == 2 ~ 4,
        TRUE ~ 3
    )
    stringr::str_pad(feature_num, width = padding)
}

.plot_roc_pr_tab <- function(data1, data2, type=c("ROC", "PR"), feature_num){
    total_feature <- sort(unique(as.character(data1$feature_num)))
    if (type=="ROC") {
        auc_col <- c("cv_fold", "feature_num", "sensitivity", "specificity")
    } else if (type=="PR") {
        auc_col <- c("cv_fold", "feature_num", "precision", "recall")
    }
    auc_data <- data1[c("ranking_method", "ml_method", "cv_fold",
                        "feature_num", paste0(type, "_AUC"))] %>%
        dplyr::rename(auc=paste0(type, "_AUC"))
    auc_res <- purrr::map(total_feature, ~.auc_stats(auc_data, .x))
    auc_chr <- purrr::map_chr(auc_res, .auc_label)
    auc_label <- purrr::map2_chr(total_feature, auc_chr, ~paste(.label_pad(.x), .y))
    auc_tab <- data2[auc_col] %>%
        dplyr::mutate(
            feature_num = as.character(feature_num),
            feature_num_label = auc_label[match(feature_num, total_feature)],
            feature_num_label1 = auc_chr[match(feature_num, total_feature)])
    plot_tab <- rbind(data1[auc_col], data2[c(auc_col)]) %>%
        dplyr::filter(feature_num==feature_num) %>%
        dplyr::mutate(cv=ifelse(cv_fold == 'mean', '1', '0'))
    return(list(auc_tab=auc_tab, plot_tab=plot_tab))
}
###
## lipid name replace ":" as "__", " " as "_", "|" as "___" or replace back
.lipid_name_replace <- function(list, type=c("replace", "revert")){
    if (type=="replace") {
        list <- stringr::str_replace_all(list, c("\\|" = "___", ":" = "__", " " = "_"))
    } else if (type=="revert") {
        list <- stringr::str_replace_all(list, c("___" = "\\|", "__" = ":", "_" = " "))
    } else {
        stop("type must be replace or revert.")
    }
    return(list)
}
