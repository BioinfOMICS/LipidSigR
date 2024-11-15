#' @title convert_sp2ratio
#' @description This function is designed to convert species abundance into
#' specific ratio-based abundance.
#' @param processed_se A SummarizedExperiment object constructed by
#' \code{\link{as_summarized_experiment}} and processed by \code{\link{data_process}}.
#' @param transform Character. Method for ratio-based abundance transformation.
#' Allowed methods include "none" and "log2". Select 'none' to skip data transformation. Default is \code{'log2'}.
#' @return Return a SummarizedExperiment object.
#' @export
#' @examples
#' data("de_data_twoGroup")
#' processed_se_twoGroup <- data_process(
#'     se=de_data_twoGroup, exclude_missing=TRUE, exclude_missing_pct=70,
#'     replace_na_method='min', replace_na_method_ref=0.5,
#'     normalization='Percentage')
#' ratioAbund_twoGroup <- convert_sp2ratio(
#'     processed_se_twoGroup, transform='log2')
#'
#' data("se_multiGroup")
#' processed_se_multiGroup <- data_process(
#'     se=se_multiGroup, exclude_missing=TRUE, exclude_missing_pct=70,
#'     replace_na_method='min', replace_na_method_ref=0.5,
#'     normalization='Percentage')
#' ratioAbund_multiGroup <- convert_sp2ratio(
#'     processed_se_multiGroup, transform='log2')

convert_sp2ratio <- function(processed_se, transform=c('none', 'log2')){

    ## Check parameter
    .check_inputSE(processed_se, metadata_list=NULL)
    if(is.null(transform) | isFALSE(transform %in% c('none', 'log2'))){
        stop("The 'transform' parameter must be either 'none' or 'log2'.")
    }
    ether.ester=odd.even=LPL=LPE=LPS=LPG=LPI=LPA=LPC=PC.PE=TG.DG=NULL
    ## Extract SE
    speciesAbund <- .extract_df(processed_se, type='abundance')
    .check_imputation(speciesAbund)
    if(isTRUE(.check_nor_negative(speciesAbund[-1]))) {
        stop("Detect negative values in the input abundance data. We recommended choosing a different normalization method during data processing.")
    }
    speciesChar <- .extract_df(processed_se, type='lipid')
    ## Species to class
    classAbund <- speciesAbund %>%
        tidyr::gather(sample_name, values, -feature) %>%
        dplyr::left_join(speciesChar[, c('feature', 'class')], by='feature') %>%
        dplyr::group_by(class, sample_name) %>%
        dplyr::summarise(values=sum(values, na.rm=TRUE), .groups='drop') %>%
        tidyr::spread(sample_name, values) %>%
        dplyr::select('feature'='class', dplyr::everything())
    ## Ether/ester, LysoPL/PL, PC/PE, TG/DG
    if(nrow(classAbund) > 0){
        # Ether/ester
        ether.ester <- .etherEsterRatio(classAbund)
        # LysoPE/PE
        LPE <- .twoClassRatio(classAbund, 'LPE', 'PE', 'LPL')
        # LysoPS/PS
        LPS <- .twoClassRatio(classAbund, 'LPS', 'PS', 'LPL')
        # LysoPG/PG
        LPG <- .twoClassRatio(classAbund, 'LPG', 'PG', 'LPL')
        # LysoPI/PI
        LPI <- .twoClassRatio(classAbund, 'LPI', 'PI', 'LPL')
        # LysoPA/PA
        LPA <- .twoClassRatio(classAbund, 'LPA', 'PA', 'LPL')
        # LysoPC/PC
        LPC <- .twoClassRatio(classAbund, 'LPC', 'PC', 'LPL')
        # PC/PE
        PC.PE <- .twoClassRatio(classAbund, 'PC', 'PE', 'others')
        # TG/DG
        TG.DG <- .twoClassRatio(classAbund, 'TG', 'DG', 'others')
    }
    ## Odd/Even
    odd.even <- .oddEvenRatio(speciesAbund, speciesChar)
    ## LysoPL/PL
    LPL <- .lysoGpRatio(speciesAbund, speciesChar)

    ## Merge all ratios
    if(all(is.null(ether.ester), is.null(odd.even), is.null(LPL), is.null(LPE),
           is.null(LPS), is.null(LPG), is.null(LPI), is.null(LPA), is.null(LPC),
           is.null(PC.PE), is.null(TG.DG))) {
        warning(paste0('There are no ratio characteristics that can be ',
                       'converted in your dataset.'))
        return(NULL)
    }
    ratio.exp <- data.table::rbindlist(list(
        ether.ester, odd.even, LPL, LPE, LPS, LPG, LPI, LPA, LPC, PC.PE, TG.DG),
        use.names=TRUE, fill=TRUE) %>%
        dplyr::select(characteristic, feature,
                      dplyr::all_of(rownames(SummarizedExperiment::colData(processed_se)))) %>%
        as.data.frame()
    ratio.mat <- ratio.exp %>%
        dplyr::mutate(index=paste0(characteristic, '|', feature)) %>%
        dplyr::select(index, everything(), -characteristic, -feature) %>%
        tibble::column_to_rownames(var='index') %>% as.matrix()
    if(transform == 'log2'){
        ratio.exp[,-c(1, 2)] <- log2(ratio.exp[,-c(1, 2)])
        ratio.mat <- log2(ratio.mat)
    }
    ## Characteristic
    char.list <- unique(ratio.exp$characteristic)
    ## SE
    ratioExpr <- SummarizedExperiment::SummarizedExperiment(
        assays=list(transform_table=ratio.mat),
        rowData=S4Vectors::DataFrame(characteristic=ratio.exp$characteristic,
                          feature=ratio.exp$feature),
        colData=SummarizedExperiment::colData(processed_se),
        metadata=list(char_list=char.list,
                      ratio_abundance=ratio.exp))
    message(paste0('There are ', length(char.list), ' ratio characteristics ',
                   'that can be converted in your dataset.'))
    return(ratio_se=ratioExpr)
}

.etherEsterRatio <- function(expr.mat){
    # Detect ether lipid class
    ether.class <- expr.mat$feature[which(stringr::str_detect(expr.mat$feature, ' O-| P-'))]
    # Detect ester lipid class
    ester.class <- expr.mat$feature[which(expr.mat$feature %in% stringr::str_remove(ether.class, ' O-| P-'))]
    # Intersection
    pair.class <- intersect(ester.class, stringr::str_remove(ether.class, ' O-| P-'))
    if(length(pair.class)==0) return(NULL)
    # Reshape expression matrix
    reshape.expr <- expr.mat %>%
        tidyr::gather(sample_name, exp, -1) %>%
        dplyr::filter(feature %in% c(pair.class, paste0(pair.class, ' O-'),
                                     paste0(pair.class, ' P-'))) %>%
        dplyr::mutate(class=stringr::str_remove(feature, ' O-| P-'),
                      group=ifelse(stringr::str_detect(feature, ' O-| P-'), 'ether', 'ester')) %>%
        dplyr::select(-feature) %>%
        dplyr::group_by(sample_name, class, group) %>%
        dplyr::summarise(exp=sum(exp, na.rm=TRUE), .groups='drop') %>%
        tidyr::spread(group, exp) %>%
        dplyr::mutate(ratio=ether/ester,
                      feature=class,
                      characteristic='Chains Ether/Ester linked ratio')
    # Create ratio matrix
    ratio.mat <- reshape.expr %>%
        dplyr::select(characteristic, feature, sample_name, ratio) %>%
        tidyr::spread(sample_name, ratio)
    return(ratio.mat)
}

.oddEvenRatio <- function(speciesAbund, speciesChar){
    if(nrow(speciesAbund) == 0) return(NULL)
    # Identify lipid species contain either odd-chain or even-chain
    chain.type <- speciesChar %>%
        dplyr::filter(feature %in% speciesAbund$feature, !is.na(FA.C)) %>%
        tidyr::separate_rows(FA.C, sep='\\|') %>%
        dplyr::mutate(type=ifelse(as.numeric(FA.C) %% 2 == 1, 'odd', 'even'))
    if(length(unique(chain.type$type)) != 2) return(NULL)
    # Choose lipid classes include both odd-chain and even-chain
    candidate <- chain.type %>%
        dplyr::group_by(class, type) %>%
        dplyr::summarise(count=dplyr::n(), .groups='drop') %>%
        tidyr::spread(type, count) %>%
        dplyr::filter(!is.na(odd), !is.na(even))
    if(nrow(candidate) == 0) return(NULL)
    # Species to characteristics
    char.mat <- speciesAbund %>%
        dplyr::left_join(chain.type[, c('feature', 'class', 'FA.C', 'type')],
                         by='feature') %>%
        dplyr::filter(class %in% candidate$class) %>%
        dplyr::select(feature, class, FA.C, type, everything()) %>%
        tidyr::gather(sample_name, exp, -1:-4) %>%
        dplyr::group_by(class, sample_name, type) %>%
        dplyr::summarise(exp=sum(exp, na.rm=TRUE), .groups='drop') %>%
        dplyr::mutate(feature=paste0(class, '_', type)) %>%
        dplyr::select(feature, sample_name, exp) %>%
        tidyr::spread(sample_name, exp)
    # Reshape expression matrix
    reshape.expr <- char.mat %>%
        tidyr::gather(sample_name, exp, -1) %>%
        tidyr::separate(feature, into=c('class', 'type'), sep='_') %>%
        tidyr::spread(type, exp) %>%
        dplyr::mutate(ratio=odd/even,
                      feature=class,
                      characteristic='Chains odd/even ratio')
    # Create ratio matrix
    ratio.mat <- reshape.expr %>%
        dplyr::select(characteristic, feature, sample_name, ratio) %>%
        tidyr::spread(sample_name, ratio)
    return(ratio.mat)
}

.twoClassRatio <- function(expr.mat, riskClass, refClass, aspect=c('LPL', 'others')){
    sub.mat <- expr.mat %>%
        dplyr::filter(feature %in% c(riskClass, refClass))
    if(nrow(sub.mat) != 2) return(NULL)
    # Reshape expression matrix
    reshape.expr <- sub.mat %>%
        tidyr::gather(sample_name, exp, -1) %>%
        tidyr::spread(feature, exp) %>%
        dplyr::select(sample_name, 'risk'=all_of(riskClass), 'ref'=all_of(refClass)) %>%
        dplyr::mutate(ratio=risk/ref,
                      feature=paste0(riskClass, '/', refClass),
                      characteristic=ifelse(aspect == 'LPL',
                                            'Ratio of Lysophospholipids to Phospholipids',
                                            'Ratio of specific lipid class A to lipid class B'))
    # Create ratio matrix
    ratio.mat <- reshape.expr %>%
        dplyr::select(characteristic, feature, sample_name, ratio) %>%
        tidyr::spread(sample_name, ratio)
    return(ratio.mat)
}

.lysoGpRatio <- function(speciesAbund, speciesChar){
    if(nrow(speciesAbund) == 0) return(NULL)
    # Identify GP lipid species
    class.type <- speciesChar %>%
        dplyr::filter(
            feature %in% speciesAbund$feature, Lipid.Maps.Category == 'GP') %>%
        dplyr::mutate(
            type=ifelse(
                stringr::str_starts(class, 'LP')
                | class %in% c('LCDPDAG', 'LCL', 'LDMPE', 'LMMPE', 'LBPA'),
                'lyso', 'pl'))
    # Species to characteristics
    char.mat <- speciesAbund %>%
        dplyr::left_join(class.type[, c('feature', 'type')], by='feature') %>%
        dplyr::filter(!is.na(type)) %>%
        dplyr::select(feature, type, everything()) %>%
        tidyr::gather(sample_name, exp, -1:-2) %>%
        dplyr::group_by(sample_name, type) %>%
        dplyr::summarise(exp=sum(exp, na.rm=TRUE), .groups='drop') %>%
        dplyr::select('feature'='type', sample_name, exp) %>%
        tidyr::spread(sample_name, exp)
    if(nrow(char.mat) != 2) return(NULL)
    # Reshape expression matrix
    reshape.expr <- char.mat %>%
        tidyr::gather(sample_name, exp, -1) %>%
        tidyr::spread(feature, exp) %>%
        dplyr::mutate(ratio=lyso/pl,
                      feature='LPL/PL',
                      characteristic='Ratio of Lysophospholipids to Phospholipids')
    # Create ratio matrix
    ratio.mat <- reshape.expr %>%
        dplyr::select(characteristic, feature, sample_name, ratio) %>%
        tidyr::spread(sample_name, ratio)
    return(ratio.mat)
}
