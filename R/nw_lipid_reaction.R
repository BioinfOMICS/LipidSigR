#' @title nw_lipid_reaction
#' @description This function returns node and edge tables that can be used to
#' generate a network using the network visualization tool. It is designed to
#' graphically represent significant lipid classes and species within lipid
#' biosynthesis pathways.
#' @param deSp_se A SummarizedExperiment object with results computed by \code{\link{deSp_twoGroup}}.
#' @param organism Character. The species to which the genes will be matched.
#' Allowed species are "human" and "mouse". Default is \code{'human'}.
#' @param show_sp Character. Determine how lipid species around the lipid class will be displayed.
#' Must be one of "all", "sigClass", and "none". Select "all" to show all species,
#' "sigClass" to show species in significant lipid classes, and "none" to not show any species.
#' Default is \code{'sigClass'}.
#' @param show_all_reactions Logical. If show_all_reactions=TURE, all the reactions will be showed.
#' Default is \code{FALSE}.
#' @param sp_significant Character. The p-value to be used for the statistically
#' significant of lipid species. Must be one of "pval" or "padj". Default is \code{'pval'}.
#' @param sp_p_cutoff Numeric. Significant level of lipid species. Default is \code{0.05}.
#' @param sp_FC_cutoff Numeric. Significance of the fold-change of lipid species. Default is \code{1}.
#' @param class_significant Character. The p-value to be used for the statistically
#' significant of lipid class. Must be one of "pval" or "padj". Default is \code{'pval'}.
#' @param class_p_cutoff Numeric. Significant level of lipid class. Default is \code{0.05}.
#' @param class_FC_cutoff Numeric. Significance of the fold-change of lipid class. Default is \code{1}.
#' @return Return a list of 4 tables.
#' \enumerate{
#' \item table_edge: a table of network edges.
#' \item table_node; a table of network nodes.
#' \item table_reaction: a table of reactions.
#' \item table_stat: a table of statistical results.
#' }
#' @export
#' @examples
#' data("de_data_twoGroup")
#' processed_se <- data_process(
#'     se=de_data_twoGroup, exclude_missing=TRUE, exclude_missing_pct=70,
#'     replace_na_method='min', replace_na_method_ref=0.5,
#'     normalization='Percentage', transform='log10')
#' deSp_se <- deSp_twoGroup(
#'     processed_se, ref_group='ctrl', test='t-test',
#'     significant='pval', p_cutoff=0.05, FC_cutoff=1, transform='log10')
#' network_table <- nw_lipid_reaction(
#'     deSp_se, organism='mouse', show_sp='sigClass', show_all_reactions=FALSE,
#'     sp_significant='pval', sp_p_cutoff=0.05, sp_FC_cutoff=1,
#'     class_significant='pval', class_p_cutoff=0.05, class_FC_cutoff=1)

nw_lipid_reaction <- function(
        deSp_se, organism=c('human', 'mouse'), show_sp=c('all', 'sigClass', 'none'),
        show_all_reactions=FALSE, sp_significant=c('pval', 'padj'), sp_p_cutoff=0.05,
        sp_FC_cutoff=1, class_significant=c('pval', 'padj'), class_p_cutoff=0.05,
        class_FC_cutoff=1){

    ## Check parameter
    .check_de_outputSE(deSp_se, de_type="deSp")
    if(is.null(organism) | isFALSE(organism %in% c('human', 'mouse'))){
        stop("The 'organism' parameter must be either 'human' or 'mouse'.")
    }
    if(is.null(show_sp) | isFALSE(show_sp %in% c('all', 'sigClass', 'none'))){
        stop("The 'show_sp' parameter must be one of 'all', 'sigClass' or 'none'.")
    }
    if(is.null(sp_significant) | isFALSE(sp_significant %in% c('pval', 'padj'))){
        stop("The 'sp_significant' parameter must be either 'pval' or 'padj'.")
    }
    if(is.null(class_significant) | isFALSE(class_significant %in% c('pval', 'padj'))){
        stop("The 'class_significant' parameter must be either 'pval' or 'padj'.")
    }
    if(!is.numeric(sp_p_cutoff) | isFALSE(.check_numeric_range(sp_p_cutoff, 0, 1))){
        stop("The 'sp_p_cutoff' parameter must be a numeric value between 0 and 1.")
    }
    if(!is.numeric(class_p_cutoff) | isFALSE(.check_numeric_range(class_p_cutoff, 0, 1))){
        stop("The 'class_p_cutoff' parameter must be a numeric value between 0 and 1.")
    }
    if(!is.numeric(sp_FC_cutoff) | isFALSE(sp_FC_cutoff >= 1)){
        stop("The minimum value for the 'sp_FC_cutoff' parameter is 1.")
    }
    if(!is.numeric(class_FC_cutoff) | isFALSE(class_FC_cutoff >= 1)){
        stop("The minimum value for the 'class_FC_cutoff' parameter is 1.")
    }
    ## Extract data from SE
    deSp <- S4Vectors::metadata(deSp_se)$all_deSp_result
    abundance <- S4Vectors::metadata(deSp_se)$processed_abundance
    lipid_char <- .extract_df(deSp_se, type = "lipid")
    group_info <- .extract_df(deSp_se, type = "group")
    ## Check two group or multiple group
    group_type <- .check_nGroup(group_info)
    if(group_type != 'two') stop('This function is used for two groups. If your data consists of multiple groups, please subset two of those groups for analysis.')
    ## Check nSample > 1
    nSample <- group_info %>% dplyr::group_by(group) %>%
        dplyr::summarise(n=dplyr::n(), .groups='drop')
    if(any(nSample$n < 2)){
        stop("This function requires each group to have at least two samples.")
    }
    transform <- S4Vectors::metadata(deSp_se)$transform
    ## Detect test from all_deSp_result
    test <- unique(deSp$method)
    ## Detect ref_group & paired_sample from group_info
    ref_group <- unique(group_info$original_group_name[which(group_info$group == 'ctrl')])
    paired_sample <- ifelse(all(is.na(group_info$pair)), FALSE, TRUE)
    ## Join all_deSp_result and lipid_char
    deLipid <- deSp %>% dplyr::left_join(lipid_char, by='feature')
    ## Get class abundance and node ID
    class_info <- .nw_classInfo(abundance, lipid_char)
    ## Lipid class abundance
    class_abund <- class_info$class.abund
    ## Network node ID
    node_id <- class_info$network.id
    #### DE class ####
    ## mean_ctrl, mean_exp, method, FC, log2FC
    class_mean_sd <- .mean_sd_twoGroup(abundance=class_abund, group_info)
    ## Data transformation
    class_abund <- .transform(abundance=class_abund, transform)
    if (length(group_info$group=="ctrl")==1 & length(group_info$group=="exp")==1) {
        res_table <- class_mean_sd %>% dplyr::select(-c("sd_ctrl", "sd_exp"))
        diff_table <- .sig_feature(res_table, significant='FC', p_cutoff=NULL, FC_cutoff)
    } else {
        ## statistic table
        class_stat_table <- .stat_twoGroup(abundance=class_abund, group_info, paired_sample, test)
        ## result table
        res_table <- class_mean_sd %>% dplyr::left_join(class_stat_table, by='feature') %>%
            as.data.frame()
        diff_table <- .sig_feature(res_table, significant=class_significant,
                                   p_cutoff=class_p_cutoff, FC_cutoff=class_FC_cutoff)
    }
    deClass <- diff_table$all_table
    deClass %<>% dplyr::left_join(node_id, by='feature')

    ## Significant lipids
    sig_species <- .sig_feature(stat_table=deLipid, significant=sp_significant,
                              p_cutoff=sp_p_cutoff, FC_cutoff=sp_FC_cutoff)[[2]]
    sig_species %<>% dplyr::filter(feature != class)
    if(nrow(sig_species) == 0) sig_species <- NULL
    ## Significant lipid classes
    sig_class <- .sig_feature(stat_table=deClass, significant=class_significant,
                              p_cutoff=class_p_cutoff, FC_cutoff=class_FC_cutoff)[[2]]
    if(nrow(sig_class) == 0) sig_class <- NULL
    ## Significant gene
    deGene <- NULL
    # gene_significant <- 'padj'
    # gene_p_cutoff <- 0.05
    # gene_FC_cutoff <- 2
    if(!is.null(deGene)){
        colnames(deGene)[which(colnames(deGene) == 'log2FoldChange')] <- 'log2FC'
        colnames(deGene)[which(colnames(deGene) == 'pvalue')] <- 'pval'
        sig_gene <- .sig_feature(stat_table=deGene, significant=gene_significant,
                                 p_cutoff=gene_p_cutoff, FC_cutoff=gene_FC_cutoff)[[2]]
        sig_gene %<>% filter(!is.na(gene_symbol))
        if(nrow(sig_gene) == 0) sig_gene <- NULL
    }else{
        sig_gene <- NULL
    }


    if(!is.null(sig_class)){
        ## Data specific edges
        edges <- .dataSpecificEdges(deClass, sig_species, sig_class, networkEdge,
                                    organism, show_all_reactions, show_sp, sig_gene)
        if(nrow(edges) == 0) return(warning('No reactions meet the criteria under these conditions; please relax your significant cutoffs.'))
        reaction_table <- edges %>%
            dplyr::filter(arrows == 'to') %>%
            dplyr::mutate(Path=paste0(from, ' -> ', to)) %>%
            tidyr::separate(
                col=title, into=c('Reaction', 'Enzyme', 'Gene',
                                  'Up-regulated gene', 'Down-regelated gene'),
                sep='</br>', fill='right') %>%
            dplyr::mutate(Reaction=stringr::str_remove(Reaction, 'Reaction: '),
                          Enzyme=stringr::str_remove(Enzyme, 'Enzyme: '),
                          Gene=stringr::str_remove(Gene, 'Gene: ')) %>%
            dplyr::mutate(Reaction=ifelse(Reaction == 'NA', NA, Reaction),
                          Enzyme=ifelse(Enzyme == 'NA', NA, Enzyme),
                          Gene=ifelse(Gene == 'NA', NA, Gene)) %>%
            dplyr::select(Path, Reaction, Enzyme, Gene)

        ## Data specific nodes
        nodes <- .dataSpecificNodes(
            edges, deLipid, deClass, sig_species, sig_class, networkNode,
            sp_significant, class_significant)
        lipid_stat <- nodes %>%
            dplyr::mutate(Level=ifelse(shape == 'square', 'Class', 'Species')) %>%
            tidyr::separate(col=title, into=c('Lipid', 'Log2FC', 'Significance'),
                            sep='</br>', fill='right') %>%
            tidyr::gather(column, value, 8:10) %>%
            dplyr::mutate(value=stringr::str_remove(value, '.*: ')) %>%
            tidyr::spread(column, value) %>%
            dplyr::mutate(Log2FC=suppressWarnings(as.numeric(Log2FC)),
                          Significance=suppressWarnings(as.numeric(Significance))) %>%
            dplyr::select(Level, Lipid, Log2FC, Significance) %>%
            dplyr::arrange(Level, Lipid)

        return(list(table_edge=edges,
                    table_node=nodes,
                    table_reaction=reaction_table,
                    table_stat=lipid_stat))
    }else{
        return(warning('No lipid classes meet these significant cutoffs; please adjust the thresholds you have set.'))
    }
}

.dataSpecificEdges <- function(deClass, sig_species, sig_class, networkEdge,
                               organism, show_all_reactions, show_sp, sig_gene){
    # Human/Mouse gene
    if(organism == 'human') networkEdge$gene <- networkEdge$reaction_gene_h
    if(organism == 'mouse') networkEdge$gene <- networkEdge$reaction_gene_m
    # Data specific lipid class subnetwork
    mapped.edge <- networkEdge %>%
        dplyr::left_join(deClass[,c('feature', 'networkId')], by=c('from'='networkId')) %>%
        dplyr::left_join(deClass[,c('feature', 'networkId')], by=c('to'='networkId')) %>%
        dplyr::mutate(from=ifelse(!is.na(feature.x), feature.x, from),
                      to=ifelse(!is.na(feature.y), feature.y, to)) %>%
        dplyr::filter(from %in% sig_class$feature | to %in% sig_class$feature)

    if(show_all_reactions == FALSE){
        remove.node <- mapped.edge %>%
            dplyr::mutate(from.sig=ifelse(from %in% sig_class$feature, 'yes', 'no'),
                          to.sig=ifelse(to %in% sig_class$feature, 'yes', 'no')) %>%
            dplyr::filter(from.sig == 'no' | to.sig == 'no') %>%
            dplyr::mutate(yes=ifelse(from.sig == 'yes', from, to),
                          no=ifelse(from.sig == 'no', from, to)) %>%
            dplyr::distinct(yes, no) %>%
            dplyr::group_by(no) %>%
            dplyr::summarise(n=dplyr::n()) %>%
            dplyr::filter(n == 1)
        mapped.edge %<>% dplyr::filter(!from %in% remove.node$no,
                                       !to %in% remove.node$no)
    }

    if(!is.null(sig_species)){
        species.edge <- .speciesEdges(mapped.edge, sig_species, sig_class, show_sp)
    }else{
        species.edge <- NULL
    }
    class.edge <- .classEdges(mapped.edge, sig_gene)
    edges <- data.table::rbindlist(list(class.edge, species.edge),
                                   use.names=TRUE, fill=TRUE)
    return(edges)
}

.speciesEdges <- function(mapped.edge, sig_species, sig_class, show_sp){
    # Data specific lipid species subnetwork
    mapped.species <- sig_species %>%
        dplyr::filter(class %in% c(mapped.edge$from, mapped.edge$to))
    if(show_sp == 'sigClass'){
        mapped.species %<>% dplyr::filter(class %in% sig_class$feature)
    }else if(show_sp == 'none'){
        mapped.species <- NULL
    }
    # Species edge table
    if(!is.null(mapped.species)){
        species.edge <- mapped.species %>%
            dplyr::mutate(label='', from=class, to=feature, title=feature,
                          color='#DDDDDD', dashes=FALSE, width=3, arrows='',
                          smooth=FALSE, shadow=FALSE) %>%
            dplyr::select(from, to, label, title, color, dashes, width, arrows,
                          smooth, shadow) %>%
            dplyr::distinct()
    }else{
        species.edge <- NULL
    }
    return(species.edge)
}

.classEdges <- function(mapped.edge, sig_gene){
    class.edge <- mapped.edge %>%
        dplyr::mutate(
            label=ifelse(is.na(EC.numbers), gene, EC.numbers),
            title=paste0('Reaction: ', Reactions, '</br>Enzyme: ', EC.numbers,
                         '</br>Gene: ', gene),
            dashes=FALSE, font.size=20, arrows='to', smooth=TRUE, shadow=FALSE)

    if(!is.null(sig_gene)){
        mapped.gene <- mapped.edge %>%
            tidyr::separate_rows(gene, sep=', ') %>%
            dplyr::filter(gene %in% sig_gene$gene_symbol)
        if(nrow(mapped.gene) > 0){
            mapped.gene %<>%
                dplyr::left_join(sig_gene, by=c('gene'='gene_symbol')) %>%
                dplyr::mutate(direction=ifelse(
                    log2FC > 0, 'Up', ifelse(log2FC < 0, 'Down', 'NS'))) %>%
                dplyr::group_by(id, Reactions, EC.numbers, direction) %>%
                dplyr::summarise(gene=paste0(gene, collapse=', '), .groups='drop') %>%
                tidyr::spread(direction, gene)
            if(!'Up' %in% colnames(mapped.gene)) mapped.gene$Up <- NA
            if(!'Down' %in% colnames(mapped.gene)) mapped.gene$Down <- NA

            class.edge %<>%
                dplyr::left_join(
                    mapped.gene, by=c("id", "Reactions", "EC.numbers")) %>%
                dplyr::mutate(
                    title=paste0(title, '</br>Up-regulated gene: ',
                                 ifelse(!is.na(Up), Up, NA),
                                 '</br>Down-regulated gene: ',
                                 ifelse(!is.na(Down), Down, NA)),
                    color=ifelse(
                        is.na(Up) & is.na(Down), '#888888',
                        ifelse(!is.na(Up) & is.na(Down), '#FF4500',
                               ifelse(is.na(Up) & !is.na(Down), '#4169E1', '#228B22'))),
                    width=ifelse(is.na(Up) & is.na(Down), 4,
                                 ifelse(!is.na(Up) & is.na(Down), 8,
                                        ifelse(is.na(Up) & !is.na(Down), 8, 8))))
        }else{
            class.edge %<>% dplyr::mutate(color='#FF8800', width=8)
        }
    }else{
        class.edge %<>% dplyr::mutate(color='#FF8800', width=8)
    }
    class.edge %<>% dplyr::select(from, to, label, title, color, dashes,
                                  font.size, width, arrows, smooth, shadow)
    return(class.edge)
}

.dataSpecificNodes <- function(
        edges, deLipid, deClass, sig_species, sig_class, networkNode,
        sp_significant, class_significant){

    lipid.class <- edges %>% dplyr::filter(arrows != '') %>%
        tidyr::gather(from_to, class, 1:2) %>% dplyr::distinct(class) %>% .$class
    lipid.species <- edges %>% dplyr::filter(arrows == '') %>%
        dplyr::distinct(to) %>% .$to

    if(length(lipid.species) > 0){
        species.node <- .speciesNodes(lipid.species, deLipid, sig_species, sp_significant)
    }else{
        species.node <- NULL
    }
    class.node <- .classNodes(lipid.class, deClass, sig_class, networkNode, class_significant, deLipid)

    nodes <- data.table::rbindlist(list(class.node, species.node),
                                   use.names=TRUE, fill=TRUE)
    return(nodes)

}

.speciesNodes <- function(lipid.species, deLipid, sig_species, sp_significant){

    species.node <- sig_species %>%
        dplyr::filter(feature %in% lipid.species) %>%
        dplyr::mutate(
            id=feature,
            label=ifelse(Level == 'SPECIES', Species.Name,
                         ifelse(Level == 'MOLECULAR_SPECIES',
                                Molecular.Species.Name, Structural.Species.Name)),
            group='species', shape='dot', size=10, font.size=16, shadow=FALSE,
            title=paste0('Lipid species: ', id, '</br>Log2FC: ',
                         format(log2FC, digits=3))) %>%
        dplyr::mutate(label=stringr::str_remove_all(stringr::str_remove(label, class), ' '))

    if('pval' %in% colnames(deLipid) | 'padj' %in% colnames(deLipid)){
        species.node %<>%
            dplyr::group_by(feature) %>%
            dplyr::mutate(
                p=ifelse(sp_significant == 'pval', pval, padj),
                size=scales::rescale(-log10(p), to=c(8, 20)),
                title=paste0(title, '</br>', sp_significant, ': ',
                             format(p, digits=3, scientific=TRUE))) %>%
            dplyr::ungroup()
    }

    species.node.color <- .speciesNodeColor(species.node)
    species.node %<>% dplyr::inner_join(species.node.color, by=c('id', 'log2FC')) %>%
        dplyr::select(id, label, group, shape, size, font.size, shadow, title, color)

}

.speciesNodeColor <- function(species.node){

    s.nodes.color <- species.node %>%
        dplyr::distinct(id, log2FC) %>%
        dplyr::mutate(color=NA) %>%
        dplyr::arrange(log2FC)

    # negative value
    s.nodes.color.neg <- s.nodes.color %>% dplyr::filter(log2FC < 0)
    if(nrow(s.nodes.color.neg) > 0){
        s.n_color.neg <- length(unique(s.nodes.color.neg$log2FC))
        s.unique_color.neg <- unique(s.nodes.color.neg$log2FC)
        for (i in 1:nrow(s.nodes.color.neg)){
            for (a in 1:s.n_color.neg) {
                if (s.nodes.color.neg$log2FC[i] == s.unique_color.neg[a]){
                    s.nodes.color.neg$color[i] <- colorRampPalette(c("#4169E1" , "#CCEEFF"))(s.n_color.neg)[a]
                }
            }
        }
    }else{
        s.nodes.color.neg <- NULL
    }
    # positive value
    s.nodes.color.pos <- s.nodes.color %>% dplyr::filter(log2FC > 0)
    if(nrow(s.nodes.color.pos) > 0){
        s.n_color.pos <- length(unique(s.nodes.color.pos$log2FC))
        s.unique_color.pos <- unique(s.nodes.color.pos$log2FC)
        for (i in 1:nrow(s.nodes.color.pos)) {
            for (a in 1:s.n_color.pos) {
                if (s.nodes.color.pos$log2FC[i] == s.unique_color.pos[a]){
                    s.nodes.color.pos$color[i] <- colorRampPalette(c("#FFE4E1", "#FF4500"))(s.n_color.pos)[a]
                }
            }
        }

    }else{
        s.nodes.color.pos <- NULL
    }

    s.nodes.color <- data.table::rbindlist(list(s.nodes.color.neg, s.nodes.color.pos),
                                           use.names=T, fill=T)
    return(s.nodes.color)

}

.classNodes <- function(lipid.class, deClass, sig_class, networkNode, class_significant,deLipid){

    deClass$p <- NA
    if(class_significant == 'pval') deClass$p <- deClass$pval
    if(class_significant == 'padj') deClass$p <- deClass$padj

    class.node <- networkNode %>%
        dplyr::left_join(deClass, by=c('id'='networkId')) %>%
        dplyr::mutate(id=ifelse(!is.na(feature), feature, id)) %>%
        dplyr::filter(id %in% lipid.class) %>%
        dplyr::distinct(id, log2FC, p) %>%
        dplyr::mutate(label=id, group='class', shape='square', size=20,
                      font.size=45, shadow=FALSE,
                      title=paste0('Lipid class: ', id, '</br>Log2FC: ',
                                   format(log2FC, digits=3)))

    if('pval' %in% colnames(deLipid) | 'padj' %in% colnames(deLipid)){
        class.node %<>%
            dplyr::mutate(size=ifelse(!id %in% sig_class$feature, NA, size)) %>%
            dplyr::mutate(
                size=ifelse(is.na(size), 20, scales::rescale(-log10(p), to=c(30, 40))),
                title=paste0(title, '</br>', class_significant, ': ',
                             format(p, digits=3, scientific=T)))
    }

    class.node.color <- .classNodeColor(class.node, sig_class)
    class.node %<>% dplyr::left_join(class.node.color, by=c('id', 'log2FC')) %>%
        dplyr::mutate(color=ifelse(is.na(color), '#CCCCCC', color)) %>%
        dplyr::select(id, label, group, shape, size, font.size, shadow, title, color)

}

.classNodeColor <- function(class.node, sig_class){

    c.nodes.color <- class.node %>%
        dplyr::filter(id %in% sig_class$feature) %>%
        dplyr::distinct(id, log2FC) %>%
        dplyr::mutate(color=NA) %>%
        dplyr::arrange(log2FC)

    ## negative value ##
    c.nodes.color.neg <- c.nodes.color %>% dplyr::filter(log2FC < 0)
    if(nrow(c.nodes.color.neg) > 0){
        c.n_color.neg <- length(unique(c.nodes.color.neg$log2FC))
        c.unique_color.neg <- unique(c.nodes.color.neg$log2FC)
        for (i in 1:nrow(c.nodes.color.neg)){
            for (a in 1:c.n_color.neg){
                if (c.nodes.color.neg$log2FC[i] == c.unique_color.neg[a]){
                    c.nodes.color.neg$color[i] <- colorRampPalette(c("#6A5ACD" , "#CCCCFF"))(c.n_color.neg)[a]
                }
            }
        }
    }else{
        c.nodes.color.neg <- NULL
    }

    ## positive value ##
    c.nodes.color.pos <- c.nodes.color %>% dplyr::filter(log2FC > 0)
    if(nrow(c.nodes.color.pos) > 0){
        c.n_color.pos <- length(unique(c.nodes.color.pos$log2FC))
        c.unique_color.pos <- unique(c.nodes.color.pos$log2FC)
        for (i in 1:nrow(c.nodes.color.pos)){
            for (a in 1:c.n_color.pos){
                if (c.nodes.color.pos$log2FC[i] == c.unique_color.pos[a]){
                    c.nodes.color.pos$color[i] <- colorRampPalette(c("#FFFACD", "#FFD700"))(c.n_color.pos)[a]
                }
            }
        }
    }else{
        c.nodes.color.pos <- NULL
    }

    c.nodes.color <- data.table::rbindlist(list(c.nodes.color.neg, c.nodes.color.pos),
                                           use.names=T, fill=T)

    return(c.nodes.color)

}


