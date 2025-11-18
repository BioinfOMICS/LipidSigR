.check_valid_inputData <- function(
        abundance, goslin_lipid_char, group_info, se_type, paired_sample) {
    # input data not null & 0*0 data frame
    .check_df_null_0(abundance, data_name="Abundance data" )
    .check_df_null_0(goslin_lipid_char,
                     data_name="Goslin lipid characteristics table" )
    # if (isFALSE(.check_dataType(abundance, "data.frame", type="partial")) ){
    #     stop("Abundance data must be a data frame.")
    # }
    # if (isFALSE(.check_dataType(goslin_lipid_char, "data.frame", type="partial")) ){
    #     stop("Goslin lipid characteristics table must be a data frame.")
    # }
    # if group_info exist, check n_group & paired_sample parameter not NULL

    if (se_type %in% c("de_two", "de_multiple", "ml", "corr") ) {
        if (!is.null(group_info) ){
            if (all(dim(group_info) == c(0, 0)) ) {
                stop("Group information table must be a data frame and cannot be data frame with 0 columns and 0 rows.")
            }
            if (se_type=="de_two" && !is.logical(paired_sample)) {
                stop("paired_sample must be logical value for two group data of differential analysis.")
            } else if (se_type!="de_two" && !is.null(paired_sample)){
                stop("paired_sample should be NULL.")
            }
        }
    } else if (se_type=="profiling") {
        if (!is.null(group_info)) {
            warning("Detect group information data, which is optional for profiling analysis.")
        }
    } else {
        stop("se_type must be 'profiling', 'de_two', 'de_multiple', 'ml', or 'corr'." )
    }
}

.check_df_null_0 <- function(data, data_name) {
    if (is.null(data)) {
        stop(data_name, " cannot be null.")
    }
    if (all(dim(data) == c(0, 0)) ) {
        stop(data_name, " must be a data frame and cannot be data frame with 0 columns and 0 rows.")
    }
}

.check_group_info <- function(abundance, group_info, se_type, paired_sample){
    if (se_type=="de_two") {
        group_col <- c("sample_name", "label_name", "group", "pair")
        if (ncol(group_info)!=4) {
            stop("Group information table must contain columns of sample_name, label_name, group, and pair.")
        } else if (isFALSE(.check_colName(group_info, group_col, type="identical"))) {
            stop("Group information table column names must be arranged in order of sample_name, label_name, group, and pair.")
        }
    } else if (se_type=="de_multiple") {
        ## group_info column names correct
        group_col <- c("sample_name", "label_name", "group")
        if (ncol(group_info)!=3) {
            stop("Group information table must contain columns of sample_name, label_name, and group.")
        } else if (isFALSE(.check_colName(group_info, group_col, type="identical"))) {
            stop("Group information table column names must be arranged in order of sample_name, label_name, and group.")
        }
    } else if (se_type=="ml") {
        group_col <- c("sample_name", "group")
        if (ncol(group_info)!=2) {
            stop("Group information table must contain columns of sample_name and group.")
        } else if (isFALSE(.check_colName(group_info, group_col, type="identical"))) {
            stop("Group information table column names must be arranged in order of sample_name and group.")
        }
    } else if (se_type=="corr") {
        if (ncol(group_info) < 3) {
            stop("Group information table must include columns for sample_name and at least two clinical terms.")
        } else if (isFALSE(.check_colName(group_info[1], "sample_name", type="identical"))) {
            stop("The first column of group information table column must be 'sample_name'.")
        } else if (isFALSE( .check_corr_input(group_info) ) ) {
            stop("Except for the first column, sample_name, the group information table must include at least two columns with numeric values and no missing values.")
        }
    }
    ## abundance colnames match group_info sample_name
    if (isFALSE(.check_colName(
        abundance[-1], group_info$sample_name, type="identical"))) {
        stop("Sample names in group information table 'sample_name' column must be as same as the sample names in lipid abundance data.")
    }
    ## group_info sample names no duplicate
    if (isFALSE(.check_no_duplicates(group_info$sample_name))) {
        stop("Sample names in group information table must be unique.")
    }
    ## special check
    if (se_type %in% c("de_two", "de_multiple")) {
        ## no NA value in "sample_name", "label_name", "group"
        if (isFALSE(
            .check_no_NA(group_info, c("sample_name", "label_name", "group")))) {
            stop("NA values are not allowed in the 'sample_name', 'label_name', and 'group' columns.")
        }
        ## group number / samples
        if (se_type=="de_two"){
            if (length(unique(group_info$group))!=2) {
                stop("The column 'group' in group information table must contain 2 groups.")
            }
            ## paired
            .check_pair_sample(group_info, paired_sample)
        } else if (se_type=="de_multiple") {
            if (length(unique(group_info$group))<=2) {
                stop("The column 'group' in group information table must contain more than 2 groups.")
            } else {
                ## check for more than 2 sample in each group
                if (isFALSE(.check_multiGroup_sample(group_info)) ){
                    stop("In Group information table of multiple groups, each group must contains more than 2 samples.")
                }
            }
        }
    } else if (se_type=="ml") {
        if (isFALSE(
            .check_no_NA(group_info, c("sample_name", "group")))) {
            stop("NA values are not allowed in the 'sample_name' and 'group' columns.")
        }
        if (isFALSE(.check_numeric_01(group_info$group) )){
            stop("The 'group' column in the group information table for machine learning analysis must be numeric, with either 0 or 1 values.")
        }
        ## Each group must have more than 30 samples.
        if (isFALSE(all(
            .check_rowNum(group_info[group_info$group==0, ], 30),
            .check_rowNum(group_info[group_info$group==1, ], 30)  ))) {

        }
    }
}

.check_pair_sample <- function(group_info, paired_sample){
    ## has pair column
    if (isFALSE(paired_sample)) {
        if (isTRUE(all(is.na(group_info$pair))) ) {
            pass <- TRUE
        } else {
            pass <- FALSE
        }
    } else {
        ## have pair samples
        pair_num <- seq_len(nrow(group_info)/2)
        if (isTRUE(setequal(as.integer(unique(group_info$pair)), pair_num))) {
            pass <- TRUE
        } else {
            pass <- FALSE
        }
    }
    if (isFALSE(pass)) {
        stop(
        "Incorrect values in group information table 'pair' column. Each pair
        should be sequentially numbered from 1 to N without any missing, blank,
        or skipped numbers; otherwise, the value must be marked as NA.")
    }
}

.check_goslin_annotation <- function(goslin_annotation){
    ## The necessary column names
    goslin.colname <- c(
        "Normalized.Name", "Original.Name", "Grammar", "Adduct",
        "Lipid.Maps.Category", "Lipid.Maps.Main.Class", "Species.Name",
        "Molecular.Species.Name", "Sn.Position.Name", "Level", "Total.C",
        "Total.OH", "Total.DB", "Mass", "Sum.Formula", "FA1.Position",
        "FA1.C", "FA1.OH", "FA1.DB", "FA1.Bond.Type", "FA1.DB.Positions",
        "FA2.Position", "FA2.C", "FA2.OH", "FA2.DB", "FA2.Bond.Type",
        "FA2.DB.Positions", "LCB.Position", "LCB.C", "LCB.OH", "LCB.DB",
        "LCB.Bond.Type", "LCB.DB.Positions", "FA3.Position", "FA3.C", "FA3.OH",
        "FA3.DB", "FA3.Bond.Type", "FA3.DB.Positions", "FA4.Position", "FA4.C",
        "FA4.OH", "FA4.DB", "FA4.Bond.Type", "FA4.DB.Positions")
    ## Check for any unrecognized lipid names.
    if (any(goslin_annotation$Grammar %in% 'NOT_PARSEABLE')) {
        sub_goslin_annotation <- goslin_annotation %>%
            dplyr::filter(Grammar != 'NOT_PARSEABLE')
        num_rm_lipid <- nrow(goslin_annotation)-nrow(sub_goslin_annotation)
        warning('Identified and removed ', num_rm_lipid,
                ' unrecognized lipid names.')
        return(sub_goslin_annotation)
    }
    ## Check column name
    if (isFALSE(.check_colName(goslin_annotation, goslin.colname, type="partial"))) {
        stop('The `goslin_annotation` table needs to include the following ',
            'columns: ', paste0(goslin.colname, collapse=', '), '.')
    }
    ## Check lipid name
    if (isFALSE(.check_no_duplicates(goslin_annotation$Original.Name))) {
        stop('The `Original.Name` column in the `goslin_annotation` table ',
             'should not contain duplicated lipid names.')
    }
    return(goslin_annotation)
}

.check_abundance <- function(abundance, se_type){
    # the first column name is “feature” (if all characters change column name)
    if (isFALSE(.check_colName(abundance[1], "feature", type="identical")) ) {
        colnames(abundance)[1] <- "feature"
    }
    if(isFALSE(.check_all_character_value(abundance$feature)) ) {
        stop("The first column of abundance data must contain a list of lipid names (features).")
    }
    ## no duplicates & NA in features
    if (isFALSE(.check_no_NA(abundance, filter_col="feature")) ) {
        stop("Lipid names (features) in abundance data must not contain NAs.")
    }
    if (isFALSE(.check_no_duplicates(abundance$feature)) ) {
        stop("Lipid names (features) in abundance data must be unique.")
    }
    #abundance values are all numeric (help convert if all values are numbers)
    if (isFALSE(.check_all_numeric_character(abundance[-1])) ){
        warning("Converting character values of lipid abundance into numeric values. Note: character string will be convert as NA")
    }
    abundance[-1] <- sapply(abundance[-1], as.numeric)
    ## feature to row names
    abundance_mat <- abundance %>% tibble::column_to_rownames(var="feature")
    ## remove features with constant values (including all NAs and all zeros) + warnings
    abundance_rmR <- .rm_constant_feature(abundance_mat)
    if (nrow(abundance_rmR) < nrow(abundance_mat)) {
        rm_rowNum <- nrow(abundance_mat) - nrow(abundance_rmR)
        # at least two features
        if (isTRUE(.check_rowNum(abundance_mat, 2))) {
            warning("Remove ", rm_rowNum, "/", nrow(abundance), " features with constant values (such as all NAs or all zeros).")
        } else {
            ## after remove the features is not enough for analysis
            stop("Remove ", rm_rowNum, "/", nrow(abundance), " features with constant values (such as all NAs or all zeros). No features for further analysis.")
        }
    }
    ## remove samples with constant values (including all NAs and all zeros)  + warnings
    abundance_rmC <- .rm_constant_sample(abundance_rmR)
    if (ncol(abundance_rmC) < ncol(abundance_mat)) {
        rm_colName <- paste(colnames(abundance_mat)[!(colnames(abundance_mat) %in% colnames(abundance_rmC))], collapse = ", ")
        rm_colNum <- ncol(abundance_mat) - ncol(abundance_rmC)

        sample_limit <- switch(
            se_type,
            "profiling"=2, "de_two"=2, "de_multiple"=2, "ml"=60, "corr"=10)
        ## at least over sample limit
        if (isTRUE(.check_colNum(abundance_rmC, sample_limit))) {
            warning("Remove ", rm_colNum, "/", nrow(abundance), " samples with constant values (including all NAs and all zeros). Remove sample: ", rm_colName, ".")
        } else {
            ## after remove the samples is not enough for analysis
            stop("Remove ", rm_colNum, "/", nrow(abundance), " samples with constant values (including all NAs and all zeros). Remove sample: ", rm_colName, ". Not enough samples for further analysis.")
        }
    }
    abundance_processed <- abundance_rmC %>% tibble::rownames_to_column("feature")
    return(abundance_processed)
}

.check_inputSE <- function(se, metadata_list=NULL) {
    ## is SE
    if (isFALSE(inherits(se, "SummarizedExperiment") )) {
        stop("The input must be a SummarizedExperiment object construct by as_summarized_experiment function or output from upstream analysis function.")
    }
    if (all(dim(se) == c(0, 0) )) {
        stop("The input must be a SummarizedExperiment object construct by as_summarized_experiment function or output from upstream analysis function.")
    }
    if (any(
        length(SummarizedExperiment::assays(se))==0,
        length(SummarizedExperiment::rowData(se))==0,
        length(SummarizedExperiment::colData(se))==0) ) {
        stop("The input must be a SummarizedExperiment object construct by as_summarized_experiment function or output from upstream analysis function.")
    }
    ## correct slot
    slot_name <- c("assays", "colData", "elementMetadata")
    if (isFALSE(all(slot_name %in% slotNames(se)) )) {
        stop("Incorrect SummarizedExperiment structure, please construct by as_summarized_experiment function or output from upstream analysis function.")
    }
    ## at least two features (only check for the original SE)
    if (is.null(metadata_list)) {
        if (isFALSE(.check_rowNum(SummarizedExperiment::assay(se), 2)))(
            stop(
                "Input SummarizedExperiment contains fewer than two features. ",
                "At least two rows (features) are required for this analysis. ",
                "Please check your object or use `as_summarized_experiment()` to reconstruct it."
            )
        )
    }
    if (!is.null(metadata_list)) {
        for (i in seq(metadata_list)) {
            ## checking condition: not null
            if (is.null(S4Vectors::metadata(se)[[metadata_list[i]]])) {
                stop("Correct SummarizedExperiment metadata not found, please use output from upstream analysis function.")
            }
        }
    }
}

## type= c("identical", "all", "partial")
.check_colName <- function(df, expect_colName, type) {
    if (type=="identical") {
        identical(colnames(df), expect_colName)
    } else if (type=="all") {
        all(colnames(df) %in% expect_colName) &&
            length(colnames(df)) == length(expect_colName)
    } else if (type=="partial") {
        all(expect_colName %in% colnames(df))
    } else {
        stop("type should be identical, all, or partial")
    }
}

.check_colNum <- function(df, col_num){
    # >= col_num
    ifelse(ncol(df) >= col_num, TRUE, FALSE)
}

.check_rowNum <- function(df, row_num){
    # >= row_num
    ifelse(nrow(df) >= row_num, TRUE, FALSE)
}


.check_no_NA <- function(df, filter_col){
    data <- df[, filter_col]
    ifelse(anyNA(data), FALSE, TRUE)
}

.check_no_duplicates <- function(x) {
    ifelse(any(duplicated(x)), FALSE, TRUE)
}

## the class is character but value is numeric
.check_all_numeric_character <- function(df){
    all(sapply(df, function(col) {
        is_numeric <- all(grepl("^[+-]?\\d*\\.?\\d+$", col[!is.na(col)]))
        char_val <- any(grepl("[:;]", col[!is.na(col)]))
        return(is_numeric && !char_val)
    }) )
}

## only contain character value (number in character type not exist)
.check_all_character_value <- function(df){
    if (isTRUE(.check_all_character(df)) ){
        all(sapply(df, function(col) {
            non_numeric <- all(!grepl("^[+-]?\\d*\\.?\\d*$", col[!is.na(col)]))
            return(non_numeric)
        }) )
    } else {
        return(FALSE)
    }
}

## for ML group
.check_numeric_01 <- function(col){
    col <- as.numeric(col)
    ifelse(all(col %in% c(0, 1)), TRUE, FALSE)
}

.check_corr_input <- function(group_info) {
    # Get all columns except the first one
    numeric_cols <- group_info[, -1, drop=FALSE]
    # Check if there are at least 2 numeric columns
    num_numeric_cols <- sum(sapply(numeric_cols, is.numeric))
    if (num_numeric_cols < 2) {
        return (FALSE)
    }
    # Check for any NA values in numeric columns
    na_cols <- any(
        sapply(numeric_cols[, sapply(numeric_cols, is.numeric), drop=FALSE],
               function(x) any(is.na(x))))
    check_res <- ifelse(isTRUE(na_cols), FALSE, TRUE)
    return (check_res)
}

.check_all_character <- function(df){
    ifelse(all(sapply(df, is.character)), TRUE, FALSE)
}

.rm_constant_sample <- function(data) {
    constant_cols <- sapply(data, function(col) {
        all_NA0 <- all(col == 0 | is.na(col))
        constant_val <- length(unique(col[!is.na(col)])) <= 1 && length((col[!is.na(col)])) > 1
        any(constant_val, all_NA0)
    })
    data_filt <- data[, !constant_cols, drop = FALSE]
    return(data_filt=data_filt)
}

.rm_constant_feature <- function(data) {
    constant_rows <- apply(data, 1, function(row) {
        all_NA0 <- all(row == 0 | is.na(row))
        constant_val <- length(unique(row[!is.na(row)])) <= 1 && length((row[!is.na(row)])) > 1
        any(constant_val, all_NA0)
    })
    data_filt <- data[!constant_rows, , drop = FALSE]
    return(data_filt=data_filt)
}

## for analysis function to check correct group number
# how to use: identical("multiple", .check_nGroup(group_info))
.check_nGroup <- function(group_info) {
    detect_groupNum <- length(unique(group_info$group))
    if (detect_groupNum==2) {
        group <- "two"
    } else if (detect_groupNum > 2) {
        # check for more than 2 sample in each group
        if (isTRUE(.check_multiGroup_sample(group_info)) ){
            group <- "multiple"
        } else {
            stop("In Group information table of multiple groups, each group must contains more than 2 samples.")
        } # end sample check
    } else {
        stop("Group information table must have 2 groups or multiple groups.")
    }
    return(group)
}

.check_multiGroup_sample <- function(group_info){
    group_name <- unique(group_info$group)
    check_sample <- c()
    for (i in seq(group_name) ) {
        if (nrow(group_info[group_info$group==group_name[i] , ]) > 1) {
            check_sample <- append(check_sample, TRUE)
        } else {
            check_sample <- append(check_sample, FALSE)
        }
    }
    ifelse(all(check_sample), TRUE, FALSE)
}

## type=c("identical", "partial")
.check_dataType <- function(data, expect_type, type) {
    if (type=="identical") {
        identical(expect_type, class(data))
    } else if (type=="partial") {
        expect_type %in% class(data)
    } else {
        stop("type should be all or partial")
    }
}

.check_numeric_range <- function(var, min_val, max_val) {
    if (!is.null(min_val) && !is.null(max_val)) {
        ifelse( (min_val <= var) && (var<= max_val), TRUE, FALSE)
    } else if (is.null(min_val) && !is.null(max_val)) {
        ifelse(var <= max_val, TRUE, FALSE)
    } else if (!is.null(min_val) && is.null(max_val)) {
        ifelse(min_val <= var, TRUE, FALSE)
    } else {
        stop("At least min_val or max_val should be input")
    }
}

.check_imputation <- function(abundance){
    if (any(sapply(abundance[-1], function(x) any(x == 0 | is.na(x)))) ) {
        stop("Detect 0 values or NAs in abundance data. Please perform data imputation by the data_process function first.")
    }
}

## type=c('common', 'deChar', 'chain_db', 'ml')
.check_char <- function(se, char, type=c('common', 'deChar', 'chain_db', 'ml')){
    #data("lipidChar")
    abundance <- .extract_df(se, type='abundance') %>%
        tidyr::gather(sample_name, values, -feature)
    lipid_char <- .extract_df(se, type='lipid') %>%
        dplyr::select(feature, dplyr::all_of(lipidChar$characteristic)) %>%
        tidyr::gather(char, char_feature, -feature) %>%
        dplyr::filter(!is.na(char_feature))
    group_info <- .extract_df(se, type='group')

    data_info <- lipid_char %>%
        dplyr::left_join(abundance, by='feature', relationship='many-to-many') %>%
        dplyr::filter(!is.na(char_feature)) %>%
        dplyr::left_join(group_info, by='sample_name')

    ## for general situation
    common_char <- data_info %>%
        dplyr::group_by(char, char_feature, sample_name) %>%
        dplyr::summarise(sum(values, na.rm=TRUE), .groups='drop') %>%
        dplyr::distinct(char) %>%
        dplyr::left_join(lipidChar, by=c('char'='characteristic')) %>%
        dplyr::arrange(match(
            aspect, c('Lipid classification', 'Fatty acid properties',
                      'Physical or chemical properties', 'Cellular component',
                      'Function')), char)
    common_char_list <- common_char$char
    if (type=="common") {
        ifelse(char %in% common_char_list, TRUE, FALSE)
    } else if (type=="deChar") {
        sp2ratio <- convert_sp2ratio(se, transform='none')
        data_ratio <- S4Vectors::metadata(sp2ratio)[["char_list"]]
        de_char_list <- c(common_char_list, data_ratio)
        ifelse(char %in% de_char_list, TRUE, FALSE)
    } else if (type=="chain_db") {
        heatmap_char <- common_char %>%
            dplyr::filter(aspect != 'Fatty acid properties')
        heatmap_char_list <- heatmap_char$char
        ifelse(char %in% heatmap_char_list, TRUE, FALSE)
    } else if (type=="ml"){
        ml_char_list <- c("none", common_char_list[common_char_list %in% c("class", "Total.OH", "Total.DB", "FA.OH", "FA.DB", "FA.C")] )
        ifelse(char %in% ml_char_list, TRUE, FALSE)
    } else {
        stop("type must be common, deChar, chain_db, ml.")
    }
}

## check no contain negative value
.check_nor_negative <- function(data) {
    ifelse(any(sapply(data, function(x) any(x < 0))), TRUE, FALSE)
}

## check for correlation df column for split into condition/ adjusted table
.check_df_corrCol <- function(group_info, condition_col, adjusted_col, cor_type) {
    if (is.null(condition_col)) {
        stop(paste0("The condition_col parameter cannot be null. Select at least two columns from ",
                    paste(colnames(group_info)[-1], collapse = ", "), "." ))
    } else if (length(condition_col) <2) {
        stop(paste0("The condition_col must include at least two columns from ",
                    paste(colnames(group_info)[-1], collapse = ", "), "." ))
    } else if (isFALSE( all(condition_col %in% colnames(group_info)[-1])) ) {
        stop(paste0("The condition_col must be column names selected from ",
                    paste(colnames(group_info)[-1], collapse = ", "), "." ))
    } else {
        condition_table <- group_info %>%
            dplyr::select(sample_name, all_of(condition_col)) %>% as.data.frame()
        if (isFALSE(.check_corr_input(condition_table)) ) {
            stop("Except for the first column, sample_name, the group information table must include at least two columns with numeric values and no missing values.")
        }
    }

    if (cor_type=="lr") {
        if (!is.null(adjusted_col)) {
            if (isFALSE( all(adjusted_col %in% colnames(group_info)[-1])) ) {
                stop(paste0("The adjusted_col must be column names selected from ",
                            paste(colnames(group_info)[-1], collapse = ", "), "." ))
            } else if (any(condition_col %in% adjusted_col)) {
                stop("The condition_col and adjusted_col must not have overlapping columns.")
            }
        } else {
            message("No adjusted variables is selected.")
        }
    }
}

## Only for checking SE in deSp & deChar downstream function
.check_de_outputSE <- function(de_se, de_type=c("deSp", "deChar", "all")){
    if (isFALSE(inherits(de_se, "SummarizedExperiment") )) {
        stop("The input must be a SummarizedExperiment object construct by upstream analysis function, e.g. deSp_twoGroup, deChar_twoGroup...")
    }
    group_info <- .extract_df(de_se, type = "group")
    group_col <- c("sample_name", "label_name", "group")
    if (isFALSE(.check_colName(group_info, group_col, type="partial")) ) {
        stop("Incorrect SummarizedExperiment structure.
             Please use output from upstream analysis function (e.g. deSp_twoGroup, deChar_twoGroup...).")
    } else {
        if (length(unique(group_info$group))<2) {
            stop("Incorrect SummarizedExperiment structure.
                 Please use output from upstream analysis function (e.g. deSp_twoGroup, deChar_twoGroup...).")
        } else {
            analyze_group <- ifelse(length(unique(group_info$group))==2, "two", "multiple")
        }
    }
    type <- ifelse("all_deSp_result" %in% names(S4Vectors::metadata(de_se)), "deSp", "deChar")
    if (!(de_type %in% c("deSp", "deChar", "all")) ) {
        stop("de_type must be deSp, deChar, or all.")
    } else if (de_type!="all" && type!=de_type) {
        stop("Incorrect SummarizedExperiment structure.
             Please use output from upstream ", de_type, "_twoGroup or ", de_type, "_multiGroup", " function.")
    }
    result_list <- c(
        paste0("all_", type, "_result"), paste0("sig_", type, "_result"), "processed_abundance")
    stat_var <- c("significant", "p_cutoff")
    if (analyze_group=="two") {
        if (type=="deSp") {
            .check_inputSE(de_se, metadata_list=c(result_list, stat_var, "FC_cutoff", "transform"))
        } else if (type=="deChar") {
            .check_inputSE(de_se, metadata_list=c(result_list, "char", stat_var, "FC_cutoff"))
        } else {
            stop("Incorrect SummarizedExperiment structure.
                 Please use output from upstream analysis function (e.g. deSp_twoGroup, deChar_twoGroup...).")
        }
    } else {
        if (type=="deSp") {
            .check_inputSE(de_se, metadata_list=c(result_list, stat_var, "transform"))
        } else if (type=="deChar") {
            .check_inputSE(de_se, metadata_list=c(result_list, "char", "post_hoc_sig", "post_hoc_p_cutoff"))
        } else {
            stop("Incorrect SummarizedExperiment structure.
                 Please use output from upstream analysis function (e.g. deSp_multiGroup, deChar_multiGroup...).")
        }
    }
}

## Only for checking SE in ML downstream function
.check_ml_outputSE <- function(ml_se){
    var_list <- c("transform", "ranking_method", "ml_method")
    result_list <- c(
        "model_result", "confusion_matrix", "roc_result", "pr_result",
        "mean_roc_result", "feature_importance_result")
    list_res <- c("selected_features", "best_model", "best_model_feature")
    .check_inputSE(ml_se, metadata_list=c("char", var_list, "nfold", "feature_option", result_list, list_res) )
    all_valid <- all(
        vapply(var_list, function(res) is.character(S4Vectors::metadata(ml_se)[[res]]), logical(1)),
        vapply(result_list, function(res) is.data.frame(S4Vectors::metadata(ml_se)[[res]]), logical(1)),
        vapply(list_res, function(res) is.list(S4Vectors::metadata(ml_se)[[res]]), logical(1))
    )
    if (!all_valid) {
        stop("Incorrect SummarizedExperiment structure. Please use output from ml_model function.")
    }
    if (length(S4Vectors::metadata(ml_se)[["feature_option"]])!=7) {
        stop("Incorrect SummarizedExperiment structure. Please use output from ml_model function.")
    }
}
