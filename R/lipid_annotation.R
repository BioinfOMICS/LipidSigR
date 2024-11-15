#' @title lipid_annotation
#' @description This function combines adjusting the annotation table
#' returned by Goslin with adding mappings between lipids and the LION ontology
#' and mappings with other resource IDs. The returned lipid characteristics table
#' will be used for further analysis in LipidSigR.
#' @param goslin_annotation A data frame of lipid characteristics from the
#' \code{\link[rgoslin]{parseLipidNames}} function of \pkg{\link{rgoslin}} package.
#' @return Return a data frame of lipid characteristics table.
#' @export
#' @examples
#' library(dplyr)
#' data("abundance_twoGroup")
#' parse_lipid <- rgoslin::parseLipidNames(abundance_twoGroup$feature)
#' recognized_lipid <- parse_lipid$Original.Name[which(parse_lipid$Grammar != 'NOT_PARSEABLE')]
#' goslin_annotation <- parse_lipid %>%
#'     dplyr::filter(Original.Name %in% recognized_lipid)
#' lipid_annotation_table <- lipid_annotation(goslin_annotation)


lipid_annotation <- function(goslin_annotation){

    if(all(goslin_annotation$Grammar == 'NOT_PARSEABLE')){
        stop(paste(
            "None of the lipid names are parseable by Goslin. Please check if",
            "your lipid dialects conform to the formats supported by Goslin."))
    }
    ## SubFunction1: modifyGoslinResults
    goslin.char <- .modifyGoslinResults(goslin_annotation)
    ## SubFunction2: Reoeder Molecular species name
    goslin.char <- .reorderMolecularName(goslin.char)
    ## SubFunction3: lionIdConversion
    lion.res <- .lionIdConversion(goslin.char)
    lion.char <- lion.res$lion.char.table
    ## SubFunction4: faChainLengthCategory3
    chain.db <- .faUnsaturationCategory(goslin.char)
    ## SubFunction5: faChainLengthCategory3
    chain.length <- .faChainLengthCategory(goslin.char)
    ## SubFunction6: gatomIdConversion
    gatom.res <- .gatomIdConversion(goslin.char)
    gatom.id <- gatom.res$gatom.chebi
    ## SubFunction7: lipidmapsNetworkIdConversion
    lipidmaps.id <- .lipidmapsNetworkIdConversion(goslin.char)
    ## SubFunction8: ID mapping
    resource.id <- .idMapping(goslin.char)
    ## SubFunction9: Merge all results
    lipid.char.all <- .mergeResult(goslin.char, lion.char, chain.db, chain.length,
                                   gatom.id, lipidmaps.id, resource.id)
    lipid.char <- lipid.char.all %>%
        dplyr::select(1:13, 40:43, 14:39, 47:62, 45, 63:64, 46, 65:73) %>%
        as.data.frame()

    return(char_table=lipid.char)

}

##################
##### Goslin #####
##################
.modifyGoslinResults <- function(goslin_annotation){

    # Parseable Lipids
    goslin.res <- goslin_annotation %>% dplyr::mutate(Index=1:nrow(.))
    goslin.res[goslin.res == ''] <- NA
    goslin.res[goslin.res == 'NA'] <- NA
    # Remove adduct
    goslin.res %<>%
        dplyr::mutate(
            Species.Name=ifelse(!is.na(Adduct) & !is.na(Species.Name),
                                stringr::str_remove(Species.Name, '\\[.*'), Species.Name),
            Molecular.Species.Name=ifelse(!is.na(Adduct) & !is.na(Molecular.Species.Name),
                                          stringr::str_remove(Molecular.Species.Name, '\\[.*'),
                                          Molecular.Species.Name),
            Sn.Position.Name=ifelse(!is.na(Adduct) & !is.na(Sn.Position.Name),
                                    stringr::str_remove(Sn.Position.Name, '\\[.*'),
                                    Sn.Position.Name))
    # Subset ST lipids
    ST <- .goslinST(goslin.res)
    # Subset SE lipids
    SE <- .goslinSE(goslin.res)
    # Subset plasmalogen (XX P-) lipids
    plasmalogen <- .goslinPlasmalogen(goslin.res)
    # Other lipids
    others <- .goslinOthers(goslin.res)
    # Merging lipid species excluding ST
    merge.res <- data.table::rbindlist(list(SE, plasmalogen, others),
                                       use.names=TRUE, fill=TRUE)
    # Species level lipids
    if('SPECIES' %in% unique(merge.res$Level)){
        no.FA.lipid <- merge.res %>% dplyr::filter(Level == 'SPECIES')
    }else{
        no.FA.lipid <- NULL
    }
    # Create FA.chain.long
    FA.chain.long <- .goslinFaChainLong(merge.res)
    # Summarise FA, FA.C, FA.DB, FA.OH into one column
    FA.chain.summ <- .goslinFaSumm(FA.chain.long)
    # Spread FA.chain.long
    FA.chain.sp <- .goslinFaChainSpread(FA.chain.long)
    # FA.chain table
    FA.chain <- FA.chain.sp %>% dplyr::inner_join(FA.chain.summ, by='index')
    FA.chain <- data.table::rbindlist(list(FA.chain, no.FA.lipid),
                                      use.names=TRUE, fill=TRUE)
    FA.chain %<>% dplyr::select(
        'feature'='Original.Name', class, Lipid.Maps.Category, Species.Name,
        Molecular.Species.Name, 'Structural.Species.Name'='Sn.Position.Name',
        Level, Mass, Sum.Formula, 'Total.FA'='Total.Chain.Sum', Total.C,
        Total.DB, Total.OH, 'LCB'='LCB_SN', LCB2, 'LCB.C'='LCB_C',
        'LCB.DB'='LCB_DB', 'LCB.OH'='LCB_OH', LCB.Bond.Type, 'FA1'='FA1_SN',
        'FA1.C'='FA1_C', 'FA1.DB'='FA1_DB', 'FA1.OH'='FA1_OH', FA1.Bond.Type,
        'FA2'='FA2_SN', 'FA2.C'='FA2_C', 'FA2.DB'='FA2_DB', 'FA2.OH'='FA2_OH',
        FA2.Bond.Type, 'FA3'='FA3_SN', 'FA3.C'='FA3_C', 'FA3.DB'='FA3_DB',
        'FA3.OH'='FA3_OH', FA3.Bond.Type, 'FA4'='FA4_SN', 'FA4.C'='FA4_C',
        'FA4.DB'='FA4_DB', 'FA4.OH'='FA4_OH', FA4.Bond.Type, FA, FA.C, FA.DB,
        FA.OH, Index) %>%
        dplyr::arrange(Index)
    # Revised Lysophospholipids
    lyso <- .goslinLysophospholipid(FA.chain)
    # Revised goslin result
    goslin.m <- data.table::rbindlist(list(lyso, ST),
                                      use.names=TRUE, fill=TRUE) %>%
        dplyr::arrange(Index)
    goslin.m[goslin.m == ''] <- NA
    return(goslin.char=goslin.m)
}

.goslinST <- function(goslin.res){
    ST <- goslin.res %>%
        dplyr::filter(stringr::str_starts(Lipid.Maps.Main.Class, 'ST')) %>%
        dplyr::mutate(class=ifelse(Lipid.Maps.Main.Class == 'ST 27:1;O',
                                   'Cholesterol', Lipid.Maps.Main.Class)) %>%
        dplyr::select(
            'feature'='Original.Name', class, Lipid.Maps.Category, Species.Name,
            Molecular.Species.Name, 'Structural.Species.Name'='Sn.Position.Name',
            Level, Mass, Sum.Formula, Index)
    return(ST)
}

.goslinSE <- function(goslin.res){
    SE <- goslin.res %>%
        dplyr::filter(stringr::str_starts(Lipid.Maps.Main.Class, 'SE')) %>%
        dplyr::mutate(
            class=ifelse(Lipid.Maps.Main.Class == 'SE 27:1', 'CE',
                         Lipid.Maps.Main.Class),
            Total.Chain.Sum=stringr::str_remove(
                Species.Name, paste0(Lipid.Maps.Main.Class, '/')),
            Species.Name=ifelse(
                Lipid.Maps.Main.Class == 'SE 27:1',
                stringr::str_replace(Species.Name, 'SE 27:1\\/', 'CE '),
                Species.Name),
            Molecular.Species.Name=ifelse(
                Lipid.Maps.Main.Class == 'SE 27:1',
                stringr::str_replace(Molecular.Species.Name, 'SE 27:1\\/', 'CE '),
                Molecular.Species.Name),
            Sn.Position.Name=ifelse(
                Lipid.Maps.Main.Class == 'SE 27:1',
                stringr::str_replace(Sn.Position.Name, 'SE 27:1\\/', 'CE '),
                Sn.Position.Name))
    return(SE)
}

.goslinPlasmalogen <- function(goslin.res){
    plasmalogen <- goslin.res %>%
        dplyr::filter(stringr::str_detect(Normalized.Name, ' P-')) %>%
        dplyr::mutate(
            class=paste0(Lipid.Maps.Main.Class, ' P-'),
            Total.DB=Total.DB - 1,
            Species.Name=ifelse(
                Total.OH == 0, paste0(class, Total.C, ':', Total.DB),
                ifelse(Total.OH == 1, paste0(class, Total.C, ':', Total.DB, ';O'),
                       paste0(class, Total.C, ':', Total.DB, ';O', Total.OH))),
            Total.Chain.Sum=ifelse(
                Total.OH == 0, paste0(Total.C, ':', Total.DB),
                ifelse(Total.OH == 1, paste0(Total.C, ':', Total.DB, ';O'),
                       paste0(Total.C, ':', Total.DB, ';O', Total.OH))))
    return(plasmalogen)
}

.goslinOthers <- function(goslin.res){
    others <- goslin.res %>%
        dplyr::filter(stringr::str_starts(Lipid.Maps.Category, 'ST', negate=TRUE)) %>%
        dplyr::filter(stringr::str_detect(Normalized.Name, ' P-', negate=TRUE)) %>%
        tidyr::separate(
            Species.Name, into=c('class', 'Total.Chain.Sum'), sep=' ', remove=F) %>%
        dplyr::mutate(
            class=ifelse(
                stringr::str_detect(Total.Chain.Sum, '-'),
                paste0(class, ' ', stringr::str_sub(Total.Chain.Sum, 1, 2)),
                class),
            Total.Chain.Sum=ifelse(
                stringr::str_detect(Total.Chain.Sum, '-'),
                stringr::str_remove(Total.Chain.Sum,
                                    stringr::str_sub(Total.Chain.Sum, 1, 2)),
                Total.Chain.Sum))
    return(others)
}

.goslinFaChainLong <- function(merge.res){
    FA.chain.long <- merge.res %>%
        dplyr::mutate(SN=ifelse(
            Level == 'MOLECULAR_SPECIES',
            stringr::str_remove(Molecular.Species.Name, class),
            ifelse(!Level %in% c('SPECIES', 'MOLECULAR_SPECIES'),
                   stringr::str_remove(Sn.Position.Name, class), NA))) %>%
        dplyr::mutate(SN=stringr::str_remove_all(SN, ' ')) %>%
        dplyr::mutate(SN=ifelse(stringr::str_starts(SN, '\\/'),
                                stringr::str_remove(SN, '\\/'), SN)) %>%
        dplyr::mutate(SN=stringr::str_replace_all(SN, '\\/', '_'),
                      index=1:nrow(.)) %>%
        tidyr::separate_rows(SN, sep='_') %>%
        dplyr::filter(!is.na(SN) | SN != '0:0') %>%
        dplyr::group_by(index) %>%
        dplyr::mutate(label=paste0('FA', dplyr::row_number())) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(label=ifelse(
            Lipid.Maps.Category == 'SP',
            stringr::str_replace(label, 'FA1', 'LCB'), label)) %>%
        dplyr::mutate(label=ifelse(
            Lipid.Maps.Category == 'SP',
            stringr::str_replace(label, 'FA2', 'FA1'), label)) %>%
        dplyr::mutate(label=ifelse(
            Lipid.Maps.Category == 'GP' & stringr::str_starts(class, 'LP')
            & FA1.C == 0 & FA1.DB == 0, 'FA2', label)) %>%
        tidyr::separate(SN, into=c('CDB', 'OH'), sep=';', remove=F,
                        fill='right', extra='merge') %>%
        tidyr::separate(CDB, into=c('C', 'DB'), sep=':') %>%
        tidyr::separate(OH, into=c('O1', 'O2'), sep=';', fill='right', extra='drop') %>%
        dplyr::mutate(OH=ifelse(stringr::str_detect(O1, 'O'), O1, O2)) %>%
        dplyr::mutate(OH=ifelse(stringr::str_detect(OH, 'O'), OH, NA)) %>%
        dplyr::mutate(
            C=as.numeric(C),
            DB=as.numeric(DB),
            OH=ifelse(is.na(OH) & !is.na(C), 0,
                      ifelse(OH == 'O', 1,
                             ifelse(nchar(OH) > 1,
                                    as.numeric(stringr::str_remove(OH, 'O')), NA)))) %>%
        dplyr::select(-O1, -O2)
    return(FA.chain.long)
}

.goslinFaSumm <- function(FA.chain.long){
    FA.chain.summ <- FA.chain.long %>%
        dplyr::group_by(index) %>%
        dplyr::reframe(FA=paste0(SN, collapse='|'),
                       FA.C=paste0(C, collapse='|'),
                       FA.DB=paste0(DB, collapse='|'),
                       FA.OH=paste0(OH, collapse='|'))
    return(FA.chain.summ)
}

.goslinFaChainSpread <- function(FA.chain.long){
    FA.chain.sp <- FA.chain.long %>%
        tidyr::pivot_wider(names_from=label, values_from=c(SN, C, DB, OH),
                           names_glue="{label}_{.value}")
    if(!'LCB_SN' %in% colnames(FA.chain.sp)){
        FA.chain.sp %<>% dplyr::mutate(LCB_SN=NA, LCB_C=NA, LCB_DB=NA, LCB_OH=NA)
    }
    FA.chain.sp %<>%
        dplyr::mutate(LCB2=ifelse(stringr::str_detect(LCB_SN, ';O2'),
                                  paste0('d', stringr::str_remove(LCB_SN, ';O2')),
                                  LCB_SN)) %>%
        dplyr::mutate(LCB2=ifelse(stringr::str_detect(LCB2, ';O3'),
                                  paste0('t', stringr::str_remove(LCB2, ';O3')),
                                  LCB2)) %>%
        tidyr::separate(col=Total.Chain.Sum, into=c('CDB', 'OH.t'), sep=';',
                        remove=FALSE, extra='merge', fill='right') %>%
        tidyr::separate(OH.t, into=c('O1.t', 'O2.t'), sep=';', fill='right', extra='drop') %>%
        dplyr::mutate(OH.t=ifelse(stringr::str_detect(O1.t, 'O'), O1.t, O2.t)) %>%
        dplyr::mutate(OH.t=ifelse(stringr::str_detect(OH.t, 'O'), OH.t, NA)) %>%
        dplyr::mutate(Total.OH=ifelse(
            is.na(OH.t), 0, ifelse(OH.t == 'O', 1,
                                   as.numeric(stringr::str_remove(OH.t, 'O')))))
    if(!'FA1_SN' %in% colnames(FA.chain.sp)){
        FA.chain.sp %<>% dplyr::mutate(FA1_SN=NA, FA1_C=NA, FA1_DB=NA, FA1_OH=NA)
    }
    if(!'FA2_SN' %in% colnames(FA.chain.sp)){
        FA.chain.sp %<>% dplyr::mutate(FA2_SN=NA, FA2_C=NA, FA2_DB=NA, FA2_OH=NA)
    }
    if(!'FA3_SN' %in% colnames(FA.chain.sp)){
        FA.chain.sp %<>% dplyr::mutate(FA3_SN=NA, FA3_C=NA, FA3_DB=NA, FA3_OH=NA)
    }
    if(!'FA4_SN' %in% colnames(FA.chain.sp)){
        FA.chain.sp %<>% dplyr::mutate(FA4_SN=NA, FA4_C=NA, FA4_DB=NA, FA4_OH=NA)
    }
    return(FA.chain.sp)
}

.goslinLysophospholipid <- function(FA.chain){
    lyso <- FA.chain %>%
        dplyr::mutate(
            FA1=ifelse(is.na(FA1) & Lipid.Maps.Category == 'GP'
                       & stringr::str_starts(class, 'LP'), '0:0', FA1),
            FA1.C=ifelse(is.na(FA1.C) & Lipid.Maps.Category == 'GP'
                         & stringr::str_starts(class, 'LP'), 0, FA1.C),
            FA1.DB=ifelse(is.na(FA1.DB) & Lipid.Maps.Category == 'GP'
                          & stringr::str_starts(class, 'LP'), 0, FA1.DB),
            FA1.OH=ifelse(is.na(FA1.OH) & Lipid.Maps.Category == 'GP'
                          & stringr::str_starts(class, 'LP'), 0, FA1.OH),
            FA2=ifelse(is.na(FA2) & Lipid.Maps.Category == 'GP'
                       & stringr::str_starts(class, 'LP'), '0:0', FA2),
            FA2.C=ifelse(is.na(FA2.C) & Lipid.Maps.Category == 'GP'
                         & stringr::str_starts(class, 'LP'), 0, FA2.C),
            FA2.DB=ifelse(is.na(FA2.DB) & Lipid.Maps.Category == 'GP'
                          & stringr::str_starts(class, 'LP'), 0, FA2.DB),
            FA2.OH=ifelse(is.na(FA2.OH) & Lipid.Maps.Category == 'GP'
                          & stringr::str_starts(class, 'LP'), 0, FA2.OH))
    return(lyso)
}

##########################################
##### Reorder Molecular Species Name #####
##########################################
.reorderMolecularName <- function(goslin.char){
    separate.chain <- goslin.char %>%
        dplyr::filter(!is.na(Molecular.Species.Name)
                      & stringr::str_detect(Molecular.Species.Name, '_')) %>%
        dplyr::select(Index, Molecular.Species.Name) %>%
        dplyr::mutate(
            Molecular.Species.Name=ifelse(
                stringr::str_starts(Molecular.Species.Name, 'PE-N\\(FA '),
                stringr::str_replace(Molecular.Species.Name, 'PE-N\\(FA ', 'PE-N\\(FA_'),
                Molecular.Species.Name)) %>%
        dplyr::mutate(Molecular.Species.Name=ifelse(
            stringr::str_detect(Molecular.Species.Name, ' O-| P-'),
            stringr::str_replace(Molecular.Species.Name, ' ', '_'),
            Molecular.Species.Name)) %>%
        dplyr::mutate(Molecular.Species.Name=ifelse(
            stringr::str_detect(Molecular.Species.Name, '_O-'),
            stringr::str_replace(Molecular.Species.Name, '_O-', '_O- '),
            Molecular.Species.Name)) %>%
        dplyr::mutate(Molecular.Species.Name=ifelse(
            stringr::str_detect(Molecular.Species.Name, '_P-'),
            stringr::str_replace(Molecular.Species.Name, '_P-', '_P- '),
            Molecular.Species.Name)) %>%
        tidyr::separate(col=Molecular.Species.Name, into=c('class', 'chain'),
                        sep=' ', remove=FALSE, extra='merge') %>%
        dplyr::mutate(class=stringr::str_replace(class, '_', ' ')) %>%
        tidyr::separate_rows(chain, sep='_')
    reorder.chain <- separate.chain %>%
        tidyr::separate(col=chain, into=c('C', 'tmp'), sep=':') %>%
        tidyr::separate(col=tmp, into=c('DB', 'Fun'), sep=';', extra='merge', fill='right') %>%
        dplyr::mutate(
            C=as.numeric(C),
            DB=as.numeric(DB)) %>%
        dplyr::group_by(Index, Molecular.Species.Name, class) %>%
        dplyr::arrange(C, DB, Fun) %>%
        dplyr::mutate(
            new.chain=ifelse(is.na(Fun), paste0(C, ':', DB),
                             paste0(C, ':', DB, ';', Fun))) %>%
        dplyr::reframe(chain=paste0(new.chain, collapse='_')) %>%
        dplyr::mutate(new.molecular.name=paste0(class, ' ', chain)) %>%
        dplyr::mutate(new.molecular.name=ifelse(
            stringr::str_detect(new.molecular.name, ' O- '),
            stringr::str_replace(new.molecular.name, ' O- ', ' O-'),
            ifelse(stringr::str_detect(new.molecular.name, ' P- '),
                   stringr::str_replace(new.molecular.name, ' P- ', ' P-'),
                   new.molecular.name)))
    reorder.res <- goslin.char %>%
        dplyr::left_join(reorder.chain[, c('Index', 'new.molecular.name')], by='Index') %>%
        dplyr::mutate(Molecular.Species.Name=ifelse(
            !is.na(new.molecular.name), new.molecular.name, Molecular.Species.Name)) %>%
        dplyr::select(-new.molecular.name)
    return(reorder.res)
}

################
##### LION #####
################
.lionIdConversion <- function(goslin.char){
    class.sys <- lipidmaps_category
    lion.char <- lion_char
    lion.db <- lion_id_anno
    lipid2lion <- lion.db %>% dplyr::distinct(name, LION_ID)
    ## Modify Class Name
    goslin.char2 <- goslin.char %>%
        dplyr::mutate(LION.Class=ifelse(
            class == 'FA', 'FFA',
            ifelse(class == 'Cholesterol', 'cholesterol',
                   ifelse(Lipid.Maps.Category == 'GP'
                          & stringr::str_starts(class, 'LP'),
                          stringr::str_sub(class, 2, -1), class))))
    #### Sphingolipids ####
    SP <- goslin.char2 %>%
        dplyr::filter(Lipid.Maps.Category == 'SP') %>%
        dplyr::mutate(LION.abbr=ifelse(
            Level == 'SPECIES',
            paste0(LION.Class, '(', Total.C, ':', Total.DB, ')'),
            paste0(LION.Class, '(', LCB2, '/', FA1, ')'))) %>%
        dplyr::mutate(LION.abbr=stringr::str_remove(LION.abbr, '\\/NA'))
    #### Lysophospholipids ####
    LPL <- goslin.char2 %>%
        dplyr::filter(
            Lipid.Maps.Category == 'GP', stringr::str_starts(class, 'LP')) %>%
        dplyr::mutate(
            LION.abbr=ifelse(Level == 'SPECIES',
                             paste0(LION.Class, '(', Total.C, ':', Total.DB, ')'),
                             paste0(LION.Class, '(', FA1, '/', FA2, ')'))) %>%
        dplyr::mutate(LION.abbr=stringr::str_replace(LION.abbr, ' O-\\(', '\\(O-')) %>%
        dplyr::mutate(LION.abbr=stringr::str_replace(LION.abbr, ' P-\\(', '\\(P-'))
    #### Others lipids ####
    ## Sort FA
    sort.FA <- goslin.char2 %>%
        dplyr::filter(!feature %in% c(SP$feature, LPL$feature)) %>%
        dplyr::filter(Level != 'SPECIES') %>%
        dplyr::select(feature, Index, FA1, FA2, FA3, FA4) %>%
        tidyr::gather(FA, chain, -1:-2) %>%
        dplyr::filter(!is.na(chain)) %>%
        dplyr::mutate(chain2=stringr::str_remove_all(chain, ';.*')) %>%
        tidyr::separate(chain2, into=c('C', 'DB'), sep=':', remove=FALSE) %>%
        dplyr::mutate(C=as.numeric(C),
                      DB=as.numeric(DB)) %>%
        dplyr::arrange(Index, C, DB) %>%
        dplyr::group_by(Index) %>%
        dplyr::mutate(FA2=paste0('LION.FA', dplyr::row_number())) %>%
        dplyr::ungroup() %>%
        dplyr::select(-FA, -chain, -C, -DB) %>%
        tidyr::spread(FA2, chain2)
    if(!'LION.FA1' %in% colnames(sort.FA)) sort.FA %<>% dplyr::mutate(LION.FA1=NA)
    if(!'LION.FA2' %in% colnames(sort.FA)) sort.FA %<>% dplyr::mutate(LION.FA2=NA)
    if(!'LION.FA3' %in% colnames(sort.FA)) sort.FA %<>% dplyr::mutate(LION.FA3=NA)
    if(!'LION.FA4' %in% colnames(sort.FA)) sort.FA %<>% dplyr::mutate(LION.FA4=NA)

    others.abbr <- goslin.char2 %>%
        dplyr::filter(!feature %in% c(SP$feature, LPL$feature)) %>%
        dplyr::left_join(sort.FA, by=c('feature', 'Index')) %>%
        dplyr::mutate(
            LION.abbr=ifelse(
                Level == 'SPECIES',
                paste0(LION.Class, '(', Total.C, ':', Total.DB, ')'),
                paste0(LION.Class, '(', LION.FA1, '/', LION.FA2, '/',
                       LION.FA3, '/', LION.FA4, ')'))) %>%
        dplyr::mutate(LION.abbr=stringr::str_remove_all(LION.abbr, '\\/NA')) %>%
        dplyr::mutate(LION.abbr=stringr::str_remove_all(LION.abbr, '\\(NA\\)')) %>%
        dplyr::mutate(LION.abbr=stringr::str_replace(LION.abbr, ' O-\\(', '\\(O-')) %>%
        dplyr::mutate(LION.abbr=stringr::str_replace(LION.abbr, ' P-\\(', '\\(P-'))
    #### Map to LION ID and characteristics ####
    lion.id <- data.table::rbindlist(list(SP, LPL, others.abbr),
                                     use.names=TRUE, fill=TRUE) %>%
        dplyr::left_join(lipid2lion, by=c('LION.abbr'='name')) %>%
        dplyr::left_join(lion.char, by=c('LION.abbr'='Name', 'LION_ID')) %>%
        dplyr::left_join(class.sys, by=c('Lipid.Maps.Category'='goslin.category'))
    #### Selected column ####
    lion.final <- lion.id %>%
        dplyr::select(feature, LION.Class, LION.abbr, 'LION.ID'='LION_ID',
                      'Category'='lipidmaps.category', Main.Class, Sub.Class,
                      Cellular.Component, Function, Bond.type, Headgroup.Charge,
                      Lateral.Diffusion, Bilayer.Thickness, Intrinsic.Curvature,
                      Transition.Temperature, Index) %>%
        dplyr::arrange(Index)
    return(list(lion.char.table=lion.final,
                lion.all.results=lion.id))

}

####################################
##### FA.Unsaturation.Category #####
####################################
.faUnsaturationCategory <- function(goslin.char){
    cate <- goslin.char %>%
        dplyr::select(feature, Index, FA.DB) %>%
        tidyr::separate_rows(FA.DB, sep="\\|") %>%
        dplyr::mutate(FA.DB=as.numeric(FA.DB)) %>%
        dplyr::mutate(Category1=ifelse(
            !is.na(FA.DB) & FA.DB < 2, 'FA with less than 2 double bonds',
            ifelse(!is.na(FA.DB) & FA.DB >= 2, 'Polyunsaturated fatty acid', NA))) %>%
        dplyr::mutate(
            Category2=ifelse(
                !is.na(FA.DB) & FA.DB == 0, 'Saturated fatty acid',
                ifelse(!is.na(FA.DB) & FA.DB == 1, 'Monounsaturated fatty acid',
                       ifelse(!is.na(FA.DB) & FA.DB == 2, 'FA with 2 double bonds',
                              ifelse(!is.na(FA.DB) & FA.DB %in% c(3:5),
                                     'FA with more than 2 double bonds',
                                     ifelse(!is.na(FA.DB) & FA.DB > 5,
                                            'FA with more than 5 double bonds', NA)))))) %>%
        dplyr::group_by(feature, Index) %>%
        dplyr::reframe(FA.Unsaturation.Category1=paste0(Category1, collapse='|'),
                       FA.Unsaturation.Category2=paste0(Category2, collapse='|'))
    cate[cate == 'NA'] <- NA
    return(cate)
}
####################################
##### FA.Chain.Length.Category #####
####################################
.faChainLengthCategory <- function(goslin.char){
    cate <- goslin.char %>%
        dplyr::select(feature, Index, FA.C) %>%
        tidyr::separate_rows(FA.C, sep="\\|") %>%
        dplyr::mutate(FA.C=as.numeric(FA.C)) %>%
        dplyr::mutate(Category1=ifelse(
            !is.na(FA.C) & FA.C <= 18, 'FA with 18 carbons or less',
            ifelse(!is.na(FA.C) & FA.C > 18, 'FA with more than 18 carbons', NA))) %>%
        dplyr::mutate(
            Category2=ifelse(
                !is.na(FA.C) & FA.C < 13, 'FA with less than 13 carbons',
                ifelse(!is.na(FA.C) & FA.C %in% c(13:15), 'FA with 13-15 carbons',
                       ifelse(!is.na(FA.C) & FA.C %in% c(16:18),
                              'FA with 16-18 carbons',
                              ifelse(!is.na(FA.C) & FA.C %in% c(19:21),
                                     'FA with 19-21 carbons',
                                     ifelse(!is.na(FA.C) & FA.C %in% c(22:24),
                                            'FA with 22-24 carbons',
                                            ifelse(!is.na(FA.C) & FA.C > 24,
                                                   'FA with more than 24 carbons', NA))))))) %>%
        dplyr::mutate(
            Category3=ifelse(
                FA.C >= 2 & FA.C <= 6, 'Short-chain FA',
                ifelse(FA.C >= 7 & FA.C <= 12, 'Medium-chain FA',
                       ifelse(FA.C >= 13 & FA.C <= 21, 'Long-chain FA',
                              ifelse(FA.C >= 22, 'Very long-chain FA', NA))))) %>%
        dplyr::group_by(feature, Index) %>%
        dplyr::reframe(FA.Chain.Length.Category1=paste0(Category1, collapse='|'),
                       FA.Chain.Length.Category2=paste0(Category2, collapse='|'),
                       FA.Chain.Length.Category3=paste0(Category3, collapse='|'))
    cate[cate == 'NA'] <- NA
    return(cate)
}

#################
##### GATOM #####
#################
.gatomIdConversion <- function(goslin.char){
    #### Sphingolipids ####
    SP <- goslin.char %>%
        dplyr::filter(Lipid.Maps.Category == 'SP') %>%
        dplyr::mutate(GATOM.abbr=ifelse(Level == 'SPECIES', Species.Name,
                                        paste0(class, '(', LCB2, '/', FA1, ')'))) %>%
        dplyr::mutate(GATOM.abbr=stringr::str_remove_all(GATOM.abbr, '\\/NA')) %>%
        dplyr::mutate(GATOM.abbr=ifelse(
            class %in% c('ACer', 'SPB', 'SPBP', 'HexCer', 'Hex2Cer', 'Hex3Cer',
                         'SHexCer', 'SHex2Cer'), Species.Name, GATOM.abbr))
    #### Lysophospholipids ####
    LPL <- goslin.char %>%
        dplyr::filter(Lipid.Maps.Category == 'GP', stringr::str_starts(class, 'LP')) %>%
        dplyr::mutate(GATOM.abbr=paste0(class, '(', Total.FA, ')')) %>%
        dplyr::mutate(GATOM.abbr=stringr::str_replace(GATOM.abbr, ' O-\\(', '\\(O-')) %>%
        dplyr::mutate(GATOM.abbr=stringr::str_replace(GATOM.abbr, ' P-\\(', '\\(P-'))
    #### FA ####
    FA <- goslin.char %>%
        dplyr::filter(class == 'FA') %>%
        dplyr::mutate(GATOM.abbr=Species.Name)
    #### MG O-, DG O-, TG O- ####
    GLO <- goslin.char %>%
        dplyr::filter(class %in% c('MG O-', 'DG O-', 'TG O-')) %>%
        dplyr::mutate(GATOM.abbr=paste0(class, '(', Total.FA, ')')) %>%
        dplyr::mutate(GATOM.abbr=stringr::str_replace(GATOM.abbr, ' O-\\(', '\\(O-')) %>%
        dplyr::mutate(GATOM.abbr=stringr::str_replace(GATOM.abbr, ' P-\\(', '\\(P-'))
    #### GL OH ####
    GLOH <- goslin.char %>%
        dplyr::filter(Lipid.Maps.Category == 'GL', Total.OH > 0) %>%
        dplyr::mutate(GATOM.abbr=Species.Name)
    #### Other lipids ####
    rm.lipid <- c(SP$feature, LPL$feature, FA$feature, GLO$feature, GLOH$feature)
    others <- goslin.char %>%
        dplyr::filter(!feature %in% rm.lipid) %>%
        dplyr::mutate(GATOM.abbr=ifelse(
            Level == 'SPECIES', paste0(class, '(', Total.FA, ')'),
            ifelse(Level == 'MOLECULAR_SPECIES',
                   paste0(class, '(', FA1, '_', FA2, '_', FA3, '_', FA4, ')'),
                   paste0(class, '(', FA1, '/', FA2, '/', FA3, '/', FA4, ')')))) %>%
        dplyr::mutate(GATOM.abbr=stringr::str_remove_all(GATOM.abbr, '\\/NA')) %>%
        dplyr::mutate(GATOM.abbr=stringr::str_remove_all(GATOM.abbr, '\\_NA')) %>%
        dplyr::mutate(GATOM.abbr=stringr::str_remove_all(GATOM.abbr, '\\(NA\\)')) %>%
        dplyr::mutate(GATOM.abbr=stringr::str_replace(GATOM.abbr, ' O-\\(', '\\(O-')) %>%
        dplyr::mutate(GATOM.abbr=stringr::str_replace(GATOM.abbr, ' P-\\(', '\\(P-'))
    #### Merge all results ####
    chebi <- data.table::rbindlist(list(SP, LPL, FA, GLO, GLOH, others),
                       use.names=TRUE, fill=TRUE)
    #### Selected column ####
    gatom.final <- chebi %>%
        dplyr::select(feature, GATOM.abbr, Index) %>%
        dplyr::arrange(Index)
    return(list(gatom.chebi=gatom.final,
                gatom.all.results=chebi))
}

##########################################
##### LIPID MAPS reaction network ID #####
##########################################
.lipidmapsNetworkIdConversion <- function(goslin.char){
    lipidmaps.class <- networkNode
    lipidmaps.class %<>% dplyr::group_by(goslin) %>% dplyr::distinct(id, goslin)
    # Map to LIPID MAPS class
    class.id <- goslin.char %>%
        dplyr::left_join(lipidmaps.class, by=c('class'='goslin')) %>%
        dplyr::arrange(Index) %>%
        # Sphingosine & Sphinganine
        dplyr::mutate(id=ifelse(
            Species.Name == 'SPB 18:1;O2', 'Sphingosine',
            ifelse(Species.Name == 'SPB 18:0;O2', 'Sphinganine', id))) %>%
        dplyr::select(feature, 'LIPIDMAPS.reaction.abbr'='id', Index)
    return(class.id)
}

###############################
##### Resource ID mapping #####
###############################
.idMapping <- function(goslin.char){
    mapping.table <- resource_id_anno
    # Species level
    if('SPECIES' %in% unique(goslin.char$Level)){
        species.lipid <- goslin.char %>% dplyr::filter(Level == 'SPECIES')
        species.id <- mapping.table %>%
            dplyr::filter(Species.Name %in% species.lipid$Species.Name) %>%
            dplyr::select(-Molecular.Species.Name, -Sn.Position.Name,
                          -core, -main_class, -sub_class) %>%
            tidyr::gather(database, databaseid, -1) %>%
            dplyr::filter(!is.na(databaseid)) %>%
            dplyr::distinct() %>%
            dplyr::arrange(Species.Name, database, databaseid) %>%
            dplyr::group_by(Species.Name, database) %>%
            dplyr::reframe(databaseid=paste0(databaseid, collapse='|')) %>%
            tidyr::spread(database, databaseid)
        species.tab <- species.lipid %>%
            dplyr::left_join(species.id, by='Species.Name')
    }else{
        species.tab <- NULL
    }
    # Molecular level
    if('MOLECULAR_SPECIES' %in% unique(goslin.char$Level)){
        molecular.lipid <- goslin.char %>% dplyr::filter(Level == 'MOLECULAR_SPECIES')
        molecular.id <- mapping.table %>%
            dplyr::filter(Molecular.Species.Name %in% molecular.lipid$Molecular.Species.Name) %>%
            dplyr::select(-Species.Name, -Sn.Position.Name,
                          -core, -main_class, -sub_class) %>%
            tidyr::gather(database, databaseid, -1) %>%
            dplyr::filter(!is.na(databaseid)) %>%
            dplyr::distinct() %>%
            dplyr::arrange(Molecular.Species.Name, database, databaseid) %>%
            dplyr::group_by(Molecular.Species.Name, database) %>%
            dplyr::reframe(databaseid=paste0(databaseid, collapse='|')) %>%
            tidyr::spread(database, databaseid)
        molecular.tab <- molecular.lipid %>%
            dplyr::left_join(molecular.id, by='Molecular.Species.Name')
    }else{
        molecular.tab <- NULL
    }
    # Sn position level
    if(any(!unique(goslin.char$Level) %in% c('SPECIES', 'MOLECULAR_SPECIES'))){
        position.lipid <- goslin.char %>% dplyr::filter(!Level %in% c('SPECIES', 'MOLECULAR_SPECIES'))
        position.id <- mapping.table %>%
            dplyr::filter(Sn.Position.Name %in% position.lipid$Structural.Species.Name) %>%
            dplyr::select(-Species.Name, -Molecular.Species.Name,
                          -core, -main_class, -sub_class) %>%
            tidyr::gather(database, databaseid, -1) %>%
            dplyr::filter(!is.na(databaseid)) %>%
            dplyr::distinct() %>%
            dplyr::arrange(Sn.Position.Name, database, databaseid) %>%
            dplyr::group_by(Sn.Position.Name, database) %>%
            dplyr::reframe(databaseid=paste0(databaseid, collapse='|')) %>%
            tidyr::spread(database, databaseid)
        position.tab <- position.lipid %>%
            dplyr::left_join(position.id,
                             by=c('Structural.Species.Name'='Sn.Position.Name'))
    }else{
        position.tab <- NULL
    }
    # Merge result
    id.mapping.res <- data.table::rbindlist(list(species.tab, molecular.tab, position.tab),
                                            use.names=TRUE, fill=TRUE) %>%
        dplyr::arrange(Index)
    if(!'LIPIDMAPS' %in% colnames(id.mapping.res)) id.mapping.res$LIPIDMAPS <- NA
    if(!'SwissLipids' %in% colnames(id.mapping.res)) id.mapping.res$SwissLipids <- NA
    if(!'HMDB' %in% colnames(id.mapping.res)) id.mapping.res$HMDB <- NA
    if(!'ChEBI' %in% colnames(id.mapping.res)) id.mapping.res$ChEBI <- NA
    if(!'KEGG' %in% colnames(id.mapping.res)) id.mapping.res$KEGG <- NA
    if(!'MetaNetX' %in% colnames(id.mapping.res)) id.mapping.res$MetaNetX <- NA
    if(!'LipidBank' %in% colnames(id.mapping.res)) id.mapping.res$LipidBank <- NA
    if(!'PubChem' %in% colnames(id.mapping.res)) id.mapping.res$PubChem <- NA
    if(!'PlantFA' %in% colnames(id.mapping.res)) id.mapping.res$PlantFA <- NA
    id.mapping.res %<>% dplyr::select(
        feature, Index, 'LIPID.MAPS.ID'='LIPIDMAPS',
        'SwissLipids.ID'='SwissLipids', 'HMDB.ID'='HMDB', 'ChEBI.ID'='ChEBI',
        'KEGG.ID'='KEGG', 'LipidBank.ID'='LipidBank', 'PubChem.CID'='PubChem',
        'MetaNetX.ID'='MetaNetX', 'PlantFA.ID'='PlantFA')
    return(id.mapping.res)
}

############################
##### Merge all result #####
############################
.mergeResult <- function(goslin.char, lion.char, chain.db, chain.length,
                         gatom.id, lipidmaps.id, resource.id){
    lipid.char <- goslin.char %>%
        dplyr::inner_join(lion.char, by=c('feature', 'Index')) %>%
        dplyr::inner_join(chain.db, by=c('feature', 'Index')) %>%
        dplyr::inner_join(chain.length, by=c('feature', 'Index')) %>%
        dplyr::inner_join(gatom.id, by=c('feature', 'Index')) %>%
        dplyr::inner_join(lipidmaps.id, by=c('feature', 'Index')) %>%
        dplyr::inner_join(resource.id, by=c('feature', 'Index')) %>%
        dplyr::arrange(Index) %>%
        dplyr::select(-Index)
    lipid.char[lipid.char == ''] <- NA
    return(lipid.char)
}

