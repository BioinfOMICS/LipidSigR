% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enrichment_lsea.R
\name{enrichment_lsea}
\alias{enrichment_lsea}
\title{enrichment_lsea}
\usage{
enrichment_lsea(
  deSp_se,
  char = NULL,
  rank_by = c("log2FC", "pval", "padj", "statistic"),
  significant = c("pval", "padj"),
  p_cutoff = 0.05,
  n_lipid = 2
)
}
\arguments{
\item{deSp_se}{The resulting SummarizedExperiment object from the differential
expression analysis function, such as \code{\link{deSp_twoGroup}} and
\code{\link{deSp_multiGroup}}.}

\item{char}{Character/NULL. A lipid characteristic selected from the common list
returned by \code{\link{list_lipid_char}}. Enter NULL to plot all characteristics.}

\item{rank_by}{Character. The method to rank lipids. Allowed rankings include
"log2FC," "pval," "padj," and "statistic". "log2FC" is only permitted for two-group data.
Default is \code{'statistic'}.}

\item{significant}{Character. The p-value to be used for the statistically
significant. Must be one of "pval" or "padj". Default is \code{'pval'}.}

\item{p_cutoff}{Numeric. The threshold to distinguish enriched lipid-sets
from not-enriched ones. Default is \code{0.05}.}

\item{n_lipid}{Numeric. The minimum number of lipids in a lipid set to be
included in enrichment. Default is \code{2}.}
}
\value{
Return a list of enrichment result, 1 interactive plot, 1 static plot, 1 table, and 2 lists.
\enumerate{
\item enrich_result: a table of enrichment result.
\item static_barPlot: a static bar plot with the top 10 significant up-regulated
and down-regulated terms for datasets involving two groups and the top 20 for multiple groups.
\item interactive_barPlot: an interactive bar plot with the top 10 significant up-regulated
and down-regulated terms for datasets involving two groups and the top 20 for multiple groups.
\item table_barPlot: table for plotting bar plots.
\item lipid_set: a list of lipid set.
\item ranked_list: lipid ranking list.
}
}
\description{
This function conducts Lipid Set Enrichment Analysis (LSEA), a
computational method for determining whether a predefined set of lipids shows
statistically significant, concordant differences among two or multiple
biological states (e.g., phenotypes).
}
\examples{
data("de_data_twoGroup")
processed_se_twoGroup <- data_process(
    se=de_data_twoGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage', transform='log10')
deSp_se_twoGroup <- deSp_twoGroup(
    processed_se_twoGroup, ref_group='ctrl', test='t-test',
    significant='pval', p_cutoff=0.05, FC_cutoff=1, transform='log10')
lsea_all_twoGroup <- enrichment_lsea(
    deSp_se_twoGroup, char=NULL, rank_by='statistic', significant='pval',
    p_cutoff=0.05)
char_list <- list_lipid_char(processed_se_twoGroup)$common_list
print(char_list)
lsea_one_twoGroup <- enrichment_lsea(
    deSp_se_twoGroup, char='class', rank_by='statistic', significant='pval',
    p_cutoff=0.05)

data("se_multiGroup")
processed_se_multiGroup <- data_process(
    se=se_multiGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage', transform='log10')
deSp_se_multiGroup <- deSp_multiGroup(
    processed_se_multiGroup, ref_group='ctrl', test='One-way ANOVA',
    significant='pval', p_cutoff=0.05, transform='log10')
lsea_all_multiGroup <- enrichment_lsea(
    deSp_se_multiGroup, char=NULL, rank_by='statistic', significant='pval',
    p_cutoff=0.05)
char_list <- list_lipid_char(processed_se_multiGroup)$common_list
print(char_list)
lsea_one_multiGroup <- enrichment_lsea(
    deSp_se_multiGroup, char='class', rank_by='statistic',
    significant='pval', p_cutoff=0.05)
}
