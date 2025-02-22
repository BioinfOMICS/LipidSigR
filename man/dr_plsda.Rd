% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dimension_reduction.R
\name{dr_plsda}
\alias{dr_plsda}
\title{dr_plsda}
\usage{
dr_plsda(
  de_se,
  ncomp = 2,
  scaling = TRUE,
  clustering = c("kmeans", "kmedoids", "hclustering", "dbscan", "group_info"),
  cluster_num = 2,
  kmedoids_metric = NULL,
  distfun = NULL,
  hclustfun = NULL,
  eps = NULL,
  minPts = NULL
)
}
\arguments{
\item{de_se}{The resulting SummarizedExperiment object from the differential
expression analysis function, such as \code{\link{deSp_twoGroup}}, \code{\link{deSp_multiGroup}},
\code{\link{deChar_twoGroup}}, and \code{\link{deChar_multiGroup}}.}

\item{ncomp}{Numeric. The number of components to include in
the model. Default is \code{2}.}

\item{scaling}{Logical. If scaling = TRUE, each block is
standardized to zero means and unit variances. Default is \code{TRUE}.}

\item{clustering}{Character. The method to be used for clustering. Allowed method
include "kmeans", "kmedoids", "hclustering", "dbscan", "group_info". Default is \code{"kmeans"}.}

\item{cluster_num}{Numeric. A positive integer specifying the number of clusters.
The number must be between 1 and 10. Default is \code{2}.}

\item{kmedoids_metric}{Character. The metric to be used for calculating
dissimilarities between observations when choosing \code{"kmedoids"} as clustering method.
Must be one of "euclidean" and "manhattan". If "kmedoids" is not selected as
the clustering method, set the value to NULL.}

\item{distfun}{Character. The distance measure to be used when choosing \code{"hclustering"}
as clustering method. Allow method include "pearson", "kendall", "spearman",
"euclidean", "manhattan", "maximum", "canberra", "binary", and "minkowski".
If "hclustering" is not selected as the clustering method, set the value to NULL.}

\item{hclustfun}{Character. The agglomeration method to be used when choosing
\code{"hclustering"} as clustering method. This should be (an unambiguous abbreviation of)
one of "ward.D", "ward.D2", "single", "complete", "average" (=UPGMA),
"mcquitty" (= WPGMA), "median" (= WPGMC), or "centroid" (= UPGMC).
If "hclustering" is not selected as the clustering method, set the value to NULL.}

\item{eps}{Numeric. The size of the epsilon neighborhood when choosing \code{"dbscan"} as clustering method.
If "dbscan" is not selected as the clustering method, set the value to NULL.}

\item{minPts}{number of minimum points in the eps region (for core points) when choosing dbscan as clustering method.}
}
\value{
Return a list with 1 data frame, 1 interactive plot, and 1 static plot.
\enumerate{
\item plsda_result: A data frame of PLS-DA data.
\item table_plsda_loading: table for plotting PLS-DA loading plot.
\item interacitve_plsda & static_plsda: PLS-DA plot.
\item interactive_loadingPlot & static_loadingPlot: PLS-DA loading plot, display the variables that contribute to the definition of each component.
}
}
\description{
Partial least squares Discriminant Analysis (PLS-DA) is a dimensionality reduction
method that transforms data from a high-dimensional space into a low-dimensional
space while retaining the original data's essential properties.
}
\examples{
data("de_data_twoGroup")
processed_se <- data_process(
    de_data_twoGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage')
deSp_se <- deSp_twoGroup(processed_se, ref_group='ctrl', test='t-test',
    significant='pval', p_cutoff=0.05, FC_cutoff=1, transform='log10')
result_plsda <- dr_plsda(deSp_se, ncomp=2, scaling=TRUE, clustering='group_info',
    cluster_num=2, kmedoids_metric=NULL, distfun=NULL, hclustfun=NULL, eps=NULL, minPts=NULL)
}
