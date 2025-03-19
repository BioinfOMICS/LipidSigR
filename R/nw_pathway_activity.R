#' @title nw_pathway_activity
#' @description This function returns node and edge tables that can be used to
#' generate a network using the network visualization tool. It computes flux
#' changes in the lipid reaction network, facilitating the identification of
#' active or suppressed pathways.
#' @param deSp_se A SummarizedExperiment object with results computed by \code{\link{deSp_twoGroup}}.
#' @param organism Character. The species to which the genes will be matched.
#' Allowed species are "human" and "mouse". Default is \code{'human'}.
#' @return Return a list of 4 tables.
#' \enumerate{
#' \item table_edge: a table of network edges.
#' \item table_node; a table of network nodes.
#' \item table_pathway_score: a table of pathway score
#' \item table_zScore: a table of z score.
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
#' network_table <- nw_pathway_activity(deSp_se, organism='mouse')

nw_pathway_activity <- function(deSp_se, organism=c('human', 'mouse')){

    ## Check parameter
    .check_de_outputSE(deSp_se, de_type="deSp")
    if(is.null(organism) | isFALSE(organism %in% c('human', 'mouse'))){
        stop("The 'organism' parameter must be either 'human' or 'mouse'.")
    }
    ## Extract data from SE
    abundance <- S4Vectors::metadata(deSp_se)$processed_abundance
    lipid_char <- .extract_df(deSp_se, type = "lipid")
    group_info <- .extract_df(deSp_se, type = "group")
    ## Check two group or multiple group
    group_type <- .check_nGroup(group_info)
    if(group_type != 'two') stop('This function is used for two groups. If your data consists of multiple groups, please subset two of those groups for analysis.')
    group_info %<>%
        dplyr::arrange(match(group, c('ctrl', 'exp')), pair)

    ## Check nSample > 1
    nSample <- group_info %>% dplyr::group_by(group) %>%
        dplyr::summarise(n=dplyr::n(), .groups='drop')
    if(any(nSample$n < 2)){
        stop("Comparing the means of two independent groups requires each group to have at least two samples.")
    }
    ## Detect ref_group & paired_sample from group_info
    ref_group <- unique(group_info$original_group_name[which(group_info$group == 'ctrl')])
    paired_sample <- ifelse(all(is.na(group_info$pair)), FALSE, TRUE)
    ## Get class abundance and node ID
    class_info <- .nw_classInfo(abundance, lipid_char)
    ## Lipid class abundance
    class_abund <- class_info$class.abund
    ## Network node ID
    node_id <- class_info$network.id
    ## Data specific node and edge
    data.node <- node_id %>% dplyr::inner_join(networkNode, by=c('networkId'='id'))
    if(nrow(data.node) == 0) return(warning('All lipid classes were not involved in the reference network; therefore, results could not be generated.'))
    data.edge <- networkEdge %>%
        dplyr::filter(from %in% data.node$networkId,
                      to %in% data.node$networkId) %>%
        dplyr::left_join(data.node[,c('feature', 'networkId')],
                         by=c('from'='networkId'), relationship="many-to-many") %>%
        dplyr::left_join(data.node[,c('feature', 'networkId')],
                         by=c('to'='networkId'), relationship="many-to-many") %>%
        dplyr::select('from'='feature.x', 'to'='feature.y', 'network.from'='from',
                      'network.to'='to', dplyr::everything()) %>%
        dplyr::distinct()
    if(organism == 'human') data.edge$gene <- data.edge$reaction_gene_h
    if(organism == 'mouse') data.edge$gene <- data.edge$reaction_gene_m

    zScore_table <- .pathwayActivityZScore(
        class_abund, group_info, edge_table=data.edge, paired_sample)
    if(all(is.na(zScore_table$p))) return(warning('Unable to perform t-test for all paths, thus results could not be generated.'))
    pathway_score <- .pathwayActivityPathwayScore(zScore_table, edge_table=data.edge)
    edges <- .pathwayActivityNetworkEdge(zScore_table)
    nodes <- .pathwayActivityNetworkNode(node_table=data.node)

    return(list(table_edge=edges,
                table_node=nodes,
                table_pathway_score=pathway_score$score,
                table_zScore=zScore_table))
}

.pathwayActivityZScore <- function(class_abund, group_info, edge_table, paired_sample){
    class_abund_gather <- class_abund %>% tidyr::gather(sample_name, abund, -feature)
    ## Fold change: to/from (product/substrate)
    FC <- edge_table %>%
        dplyr::select(from, to) %>%
        dplyr::left_join(class_abund_gather, by=c('from'='feature'),
                         relationship="many-to-many") %>%
        dplyr::left_join(class_abund_gather, by=c('to'='feature', 'sample_name')) %>%
        dplyr::mutate(FC=abund.y/abund.x) %>%
        dplyr::select(dplyr::everything(), 'from_abund'='abund.x', 'to_abund'='abund.y')
    ## T-test
    ttest <- FC %>%
        dplyr::left_join(group_info[,c('sample_name', 'group')], by='sample_name') %>%
        dplyr::group_by(from, to, group) %>%
        dplyr::summarise(list.FC=list(FC)) %>%
        tidyr::spread(group, list.FC)
    ttest %<>% dplyr::mutate(p=tryCatch({
        stats::t.test(unlist(exp), unlist(ctrl), paired=paired_sample,
                      alternative='greater')$p.value},
        error=function(e){ NA }))
    ## z-score
    z <- ttest %>% dplyr::mutate(z=stats::qnorm(1 - p)) %>%
        dplyr::select(-exp, -ctrl) %>%
        dplyr::left_join(edge_table[,c('from', 'to', 'Reactions', 'EC.numbers',
                                      'Crossreferences', 'gene')],
                         by=c('from', 'to')) %>% as.data.frame()
    return(z)
}

.pathwayActivityPathwayScore <- function(zScore_table, edge_table){
    # Convert the DataFrame to an igraph object
    ig.o <- igraph::graph_from_data_frame(edge_table[,1:2], directed=TRUE)
    # Initialize an empty DataFrame to store the paths
    path_df <- NULL
    # Loop through each pair of vertices to find all simple paths
    for(from_node in names(igraph::V(ig.o))){
        for(to_node in names(igraph::V(ig.o))){
            if(from_node != to_node){
                paths <- igraph::all_simple_paths(ig.o, from=from_node, to=to_node)
                if(length(paths) == 0) next
                for(i in 1:length(paths)){
                    path_name <- paste0(from_node, 'to', to_node, '.', i)
                    node_seq <- paste0(names(paths[[i]]), collapse=' -> ')
                    path_df <- rbind(path_df, data.frame(path_name=path_name,
                                                         nodes=node_seq,
                                                         stringsAsFactors=FALSE))
                }
            }
        }
    }
    ## Separate simple path
    split.path <- path_df %>%
        dplyr::mutate(path_length=stringr::str_count(nodes, ' -> ')) %>%
        tidyr::separate_rows(nodes, sep=' -> ') %>%
        dplyr::group_by(path_name) %>%
        dplyr::mutate(index1=dplyr::row_number(),
                      index2=ifelse(index1 == max(index1), NA, index1 + 1))
    ## Construct sub path
    sub.path <- split.path %>%
        dplyr::left_join(split.path[,c('path_name', 'nodes', 'index1')],
                         by=c('path_name', 'index2'='index1')) %>%
        dplyr::filter(!is.na(index2)) %>%
        dplyr::left_join(zScore_table, by=c('nodes.x'='from', 'nodes.y'='to'),
                         relationship="many-to-many") %>%
        dplyr::distinct(path_name, nodes.x, nodes.y, .keep_all=TRUE)
    ## Calculate pathway score
    path.score <- sub.path %>%
        dplyr::group_by(path_name) %>%
        dplyr::summarise(pathway_score=sum(z, na.rm=TRUE)/sqrt(mean(path_length)),
                         genes=paste0(gene, collapse=', ')) %>%
        tidyr::separate_rows(genes, sep=', ') %>%
        dplyr::group_by(path_name, pathway_score) %>%
        dplyr::distinct(genes, .keep_all=TRUE) %>%
        dplyr::mutate(N=dplyr::n(),
                      genes=ifelse(N == 1 & genes == 'NA', NA, genes)) %>%
        dplyr::filter(is.na(genes) | genes != 'NA') %>%
        dplyr::arrange(path_name, genes) %>%
        dplyr::summarise(genes=paste0(genes, collapse=', ')) %>%
        dplyr::inner_join(path_df, by='path_name') %>%
        dplyr::select(path_name, 'reaction_chain'='nodes', pathway_score, genes) %>%
        dplyr::arrange(dplyr::desc(pathway_score))
    return(list(split.path=sub.path,
                score=path.score))
}

.pathwayActivityNetworkEdge <- function(zScore_table){
    edges <- zScore_table %>%
        dplyr::mutate(
            label=format(z, digits=3), arrows='to', font.size=20, shadow=FALSE,
            width=scales::rescale(abs(z), to=c(1, 7)), dashes=FALSE, smooth=TRUE,
            title=paste0('Reaction: ', Reactions, '</br>Pathway score: ', label,
                         '</br>Enzyme: ', EC.numbers, '</br>Gene: ',
                         stringr::str_wrap(gene, 30)),
            color=ifelse(z > 1.645, '#FF4500',
                         ifelse(z <=1.645 & z > 0, '#FFD2D2',
                                ifelse(z < 0 & z >= -1.645, '#DDDDFF',
                                       ifelse(z < -1.645, '#4169E1', '#666666'))))) %>%
        dplyr::select(from, to, label, width, arrows, dashes, title, smooth,
                      shadow, color, font.size)
    return(edges)
}

.pathwayActivityNetworkNode <- function(node_table){
    nodes <- node_table %>%
        dplyr::mutate(id=feature, label=feature, size=30, shape='dot',
                      title=feature, color='#DDDDDD', border='#888888',
                      font.size=30, shadow=FALSE) %>%
        dplyr::select(id, label, size, shape, title, color, border,
                      font.size, shadow) %>%
        dplyr::distinct()
    return(nodes)
}



