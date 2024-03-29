% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DE_species.R
\name{DE_species}
\alias{DE_species}
\title{DE_species}
\usage{
DE_species(
  exp_data_SE,
  data_transform = TRUE,
  paired = FALSE,
  test = "t.test",
  adjust_p_method = "BH",
  sig_stat = "p.adj",
  sig_pvalue = 0.05,
  sig_FC = 2
)
}
\arguments{
\item{exp_data_SE}{A SummarizedExperiment object contains information about 
various features, such as molecules, lipid class, etc., and their 
expression levels for each sample. The assay is a matrix representing the 
expression of lipid features across all samples. Missing values (NAs)  are 
allowed. The row data  corresponds to specific lipid features, such as class 
and total length. The first column's name must be "feature" (lipid species), 
and NAs are allowed for this data. The column data comprises sample names, 
sample labels, group names, and pair numbers that represent 'the pair' for 
conducting t-tests or Wilcoxon tests. NAs are allowed.}

\item{data_transform}{Logical. If data_transform = TRUE,
transform exp_data by log10.}

\item{paired}{Logical. If paired = TRUE,
data are paired samples. (default: FALSE)}

\item{test}{A character string indicating which method to be used for
comparing means. Allowed method include \bold{"t.test"} and
\bold{"wilcox.test"}. (default: "t.test")}

\item{adjust_p_method}{Correction method, a character string. One of
\bold{"holm"}, \bold{"hochberg"}, \bold{"hommel"}, \bold{"bonferroni"},
\bold{"BH"}, \bold{"BY"}, \bold{"fdr"}, \bold{"none"}, can be
abbreviated. (default: "BH")}

\item{sig_stat}{A character string indicating which p-value is to be used
for the statistically significant. One of \bold{"p.adj"} or
\bold{"p"}. (default: "p.adj")}

\item{sig_pvalue}{Numeric. Significant level. (default: 0.05)}

\item{sig_FC}{Numeric. Significance of the fold-change. (default: 2)}
}
\value{
Return a list with 2 SummarizedExperiment objects.
\enumerate{
\item exp_data_stat: an SE object comprises the results of
differential expression analysis, including fold change, p-value,
adjusted p-value.
\item exp_data_stat_sig: an SE object comprises the significant
results of differential expression analysis, including fold change,
p-value, adjusted p-value.
}
}
\description{
Compute differentially expressed analysis of
two groups (independent) based on the group Information to find
significant lipid species. To find differentially expressed lipid species,
a statistical method (t-test or Wilcoxon test) and the cut-offs for
significant lipid species need to be chosen, and the p-value will then be
adjusted by Benjamini-Hochberg or other methods.
}
\examples{
data("DE_data")
exp_data_SE <- DE_data
exp_transform_non_log <- data_process(exp_data_SE, exclude_var_missing=TRUE,
    missing_pct_limit=50, replace_zero=TRUE, zero2what='min', xmin=0.5, 
    replace_NA=TRUE, NA2what='min', ymin=0.5, pct_transform=TRUE,
    data_transform=FALSE, trans_type='log', centering=FALSE, scaling=FALSE)
DE_species(exp_data_SE=exp_transform_non_log, data_transform=FALSE, 
    paired=FALSE, test='wilcoxon test', adjust_p_method='BH',
    sig_stat='p.adj', sig_pvalue=0.05, sig_FC=2)
}
