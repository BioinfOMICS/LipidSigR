% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ml_shap.R
\name{plot_shap_sample}
\alias{plot_shap_sample}
\title{plot_shap_sample}
\usage{
plot_shap_sample(shap_se, sample_id = 10)
}
\arguments{
\item{shap_se}{A SummarizedExperiment object with results computed by \code{\link{ml_shap}}.}

\item{sample_id}{Numeric. The number of samples to display for each feature. Default is \code{10}.}
}
\value{
Return 1 interactive plot, 1 static plot, and 1 table.
\enumerate{
\item interactive_sample_feature_importance & static_sample_feature_importance:
SHAP feature importance plot.
\item table_sample_feature_importance: table for plotting SHAP feature importance.
}
}
\description{
This function plots the SHAP feature importance results of user-selected sample.
}
\examples{
data("ml_se_sub")
shap_se <- ml_shap(ml_se_sub, feature_num=10, nsim=5)
sample_id_list <- unique(S4Vectors::metadata(shap_se)$shap_result$ID)
res <- plot_shap_sample(shap_se, sample_id=sample_id_list[10])
}
