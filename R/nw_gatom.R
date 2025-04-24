#' @title nw_gatom
#' @description This function returns node and edge tables that can be used to
#' generate a network using the network visualization tool. It isolates
#' significant subnetworks within a constructed metabolite-level network.
#' @param deSp_se A SummarizedExperiment object with results computed by \code{\link{deSp_twoGroup}}.
#' @param organism Character. The species to which the genes will be matched.
#' Allowed species are "human" and "mouse". Default is \code{'human'}.
#' @param n_lipid Numeric. The number of lipids to be scored positively determines
#' the size of the resulting module. The higher the number, the larger the module.
#' Default is \code{50}.
#' @param sp_significant Character. The p-value to be used for the statistically
#' significant of lipid species. Must be one of "pval" or "padj". Default is \code{'pval'}.
#' @param sp_p_cutoff Numeric. Significant level of lipid species. Default is \code{0.05}.
#' @param sp_FC_cutoff Numeric. Significance of the fold-change of lipid species. Default is \code{1}.
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
#' network_table <- nw_gatom(
#'     deSp_se, organism='mouse', n_lipid=50, sp_significant='pval',
#'     sp_p_cutoff=0.05, sp_FC_cutoff=1)

nw_gatom <- function(
        deSp_se, organism=c('human', 'mouse'), n_lipid=50,
        sp_significant=c('pval', 'padj'), sp_p_cutoff=0.05, sp_FC_cutoff=1){

    ## Check de SE
    .check_de_outputSE(deSp_se, de_type="deSp")

    ## Check parameter
    if(is.null(organism) | isFALSE(organism %in% c('human', 'mouse'))){
        stop("The 'organism' parameter must be either 'human' or 'mouse'.")
    }
    if(!is.numeric(n_lipid) | isFALSE(n_lipid >= 1)){
        stop("The minimum value for the 'n_lipid' parameter is 1.")
    }
    if(is.null(sp_significant) | isFALSE(sp_significant %in% c('pval', 'padj'))){
        stop("The 'sp_significant' parameter must be either 'pval' or 'padj'.")
    }
    if(!is.numeric(sp_p_cutoff) | isFALSE(.check_numeric_range(sp_p_cutoff, 0, 1))){
        stop("The 'sp_p_cutoff' parameter must be a numeric value between 0 and 1.")
    }
    if(!is.numeric(sp_FC_cutoff) | isFALSE(sp_FC_cutoff >= 1)){
        stop("The minimum value for the 'sp_FC_cutoff' parameter is 1.")
    }
    ## Extract data from SE
    deSp <- S4Vectors::metadata(deSp_se)$all_deSp_result
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

    deLipid <- deSp %>% dplyr::inner_join(lipid_char, by='feature')
    if(colnames(deLipid)[1] != 'GATOM.abbr'){
        deLipid %<>% dplyr::select(GATOM.abbr, log2FC, pval, padj, everything())
    }

    if(organism == 'human'){
        org.gatom.anno <- gatom_Hs_anno
    }else if(organism == 'mouse'){
        org.gatom.anno <- gatom_Mm_anno
    }

    ## Significant lipids
    if(sp_significant == 'pval') deLipid %<>% dplyr::mutate(p=pval)
    if(sp_significant == 'padj') deLipid %<>% dplyr::mutate(p=padj)
    deLipid %<>%
        dplyr::filter(abs(log2FC) > log2(sp_FC_cutoff), p < sp_p_cutoff) %>%
        dplyr::filter(!is.na(GATOM.abbr))
    if(nrow(deLipid) == 0) return(warning('No lipid species meet these significant cutoffs; please adjust the thresholds you have set.'))
    ## Significant genes
    deGene=n_gene=NULL
    if(!is.null(deGene)){
        if(gene_significant == 'pval') deGene %<>% dplyr::mutate(p=pvalue)
        if(gene_significant == 'padj') deGene %<>% dplyr::mutate(p=padj)
        deGene %<>%
            dplyr::filter(
                abs(log2FoldChange) > log2(gene_FC_cutoff), p < gene_p_cutoff) %>%
            dplyr::select('ID'='gene_symbol', 'pval'='pvalue',
                          'log2FC'='log2FoldChange', baseMean)
        if(nrow(deGene) == 0){
            deGene <- NULL
            return(warning('No genes meet these significant cutoffs; please adjust the thresholds you have set.'))
        }
    }


    ## Getting atom graph
    g <- gatom::makeMetabolicGraph(
        network=gatom_network , topology='metabolites', org.gatom.anno,
        gene.de=deGene, gene2reaction.extra=gatom_gene2reaction,
        met.db=gatom_lipid_db, met.de=deLipid)
    ## Scoring graph, obtaining an instance of SGMWCS problem
    gs <- gatom::scoreGraph(g, k.gene=n_gene, k.met=n_lipid)
    # Initialize an SMGWCS solver
    solver <- mwcsr::rnc_solver()
    # Finding a module
    res <- tryCatch(
        {
            mwcsr::solve_mwcsp(solver, gs)
        },
        error = function(e) {
            message(
                "[MWCSP solve failed]\n",
                "Error message: ", e$message, "\n"
            )
            NULL
        }
    )
    if (is.null(res)) {
        return(NULL)
    }

    m <- res$graph
    # node, edge table
    node.table <- data.table::data.table(igraph::as_data_frame(m, what="vertices"))
    if(!any(colnames(node.table) %in% 'feature')){
        node.table %<>%
            dplyr::left_join(deLipid[,c('feature', 'log2FC', 'pval', 'padj')],
                             by=c('log2FC', 'pval', 'padj')) %>%
            dplyr::distinct()
    }
    edge.table <- data.table::data.table(igraph::as_data_frame(m, what="edges"))
    if(all(is.na(node.table$log2FC)) | nrow(edge.table) == 0){
        return(warning('No moudles meet the criteria under these conditions; please relax your significant cutoffs.'))
    }

    nodes <- .gatomNetworkNode(node.table)
    edges <- .gatomNetworkEdge(edge.table, nodes, deGene)


    reaction_table <- edges %>% dplyr::select(Path, Reaction, Enzyme, 'Gene'='label')
    lipid_stat <- nodes %>%
        dplyr::select('Lipid'='lipid', 'Log2FC'='log2FC', 'p-value'='pval') %>%
        dplyr::arrange(dplyr::desc(Log2FC))

    return(list(table_edge=edges,
                table_node=nodes,
                table_reaction=reaction_table,
                table_stat=lipid_stat))
    # lm1 <- abbreviateLabels(m, orig.names=TRUE, abbrev.names=TRUE)
    # createShinyCyJSWidget(lm1)
}

.gatomNetworkNode <- function(node.table){

    node.table %<>%
        dplyr::mutate(id=name,
                      lipid=ifelse(is.na(feature), label, feature),
                      group=ifelse(!is.na(feature), "User's", "Others"),
                      label=stringr::str_wrap(lipid, width=35),
                      shape='dot',
                      size=ifelse(is.na(pval), 8,
                                  scales::rescale(-log10(pval), to=c(15, 30))),
                      font.size=20, shadow=FALSE,
                      title=paste0('Lipid: ', label, '</br>Group: ', group,
                                   ' lipids', '</br>CHEBI ID: ',
                                   stringr::str_remove(name, 'CHEBI:'),
                                   '</br>Log2FC: ', format(log2FC, digits=3),
                                   '</br>pval: ', format(pval, digits=3, scientific=T)))

    node.color <- .gatomNodeColor(node.table)

    nodes <- node.table %>%
        dplyr::left_join(node.color, by=c('id', 'log2FC')) %>%
        dplyr::mutate(color=ifelse(is.na(color), '#CCCCCC', color)) %>%
        dplyr::select(id, label, group, shape, size, font.size, shadow, title,
                      color, lipid, log2FC, pval)

    return(nodes)
}

.gatomNodeColor <- function(node.table){

    s.nodes.color <- node.table %>%
        dplyr::filter(!is.na(log2FC)) %>%
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

.gatomNetworkEdge <- function(edge.table, nodes, deGene){
    rhea_reaction$ID <- as.character(rhea_reaction$ID)
    edge.table %<>% dplyr::filter(from != to) %>%
        dplyr::left_join(rhea_reaction, by=c('reaction'='ID')) %>%
        dplyr::mutate(dashes=FALSE, font.size=20, width=2, smooth=FALSE,
                      shadow=FALSE, color='#DDDDDD') %>%
        dplyr::left_join(nodes[, c('id', 'lipid')], by=c('from'='id')) %>%
        dplyr::left_join(nodes[, c('id', 'lipid')], by=c('to'='id')) %>%
        dplyr::mutate(Path=paste0(lipid.x, ' -> ', lipid.y),
                      Enzyme=ifelse(enzyme == '-.-.-.-', NA, enzyme),
                      title=paste0(paste0('Reaction: ', DEFINITION,
                                          '</br>Rhea ID: ', reaction,
                                          '</br>Gene symbol: ', label,
                                          '</br>Enzyme: ', Enzyme))) %>%
        dplyr::select(from, to, 'from_label'='lipid.x', 'to_label'='lipid.y',
                      label, everything())

    if(!is.null(deGene)){
        edge.table %<>%
            dplyr::mutate(color=ifelse(is.na(log2FC), '#DDDDDD',
                                       ifelse(log2FC > 0, '#FF4500', '#4169E1')),
                          width=ifelse(is.na(log2FC), 2, 5),
                          title=paste0(title, '</br>Log2FC: ',
                                       format(log2FC, digits=3),
                                       '</br>pval: ', format(pval, digits =3,
                                                             scientific=T)))
    }else{
        edge.table %<>% dplyr::mutate(log2FC=NA)
    }

    edges <- edge.table %>%
        dplyr::select(from, to, from_label, to_label, label, dashes, font.size,
                      width, smooth, shadow, color, title, Path,
                      'Reaction'='DEFINITION', Enzyme, log2FC, pval)
    return(edges)
}


