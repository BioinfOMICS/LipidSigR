\name{group_info_multiOmics}
\alias{group_info_multiOmics}
\docType{data}
\title{
group_info_multiOmics
}
\description{
A data frame of grouping information for the samples in
\code{\link{abundance_lipidOmics}} and \code{\link{abundance_rna}}, used for
the multi-omics integration workflow (\code{\link{as_MultiAssayExperiment}}.
Samples come from a hepatocellular carcinoma (HCC) cohort of 98 patients
recruited between 2011 and 2013, and the \code{group} column ('tsne0'/'tsne1')
divides the 59 patients with paired lipidomic and transcriptomic profiling into
two subgroups based on t-SNE clustering of their lipidomic profiles.
}
\usage{data("group_info_multiOmics")}
\format{
  A data frame with 59 observations on the following 4 variables.
  \describe{
    \item{\code{sample_name}}{a character vector of de-identified sample IDs,
    matching the sample columns in \code{abundance_lipidOmics} and
    \code{abundance_rna}.}
    \item{\code{label_name}}{a character vector of sample labels (identical
    to \code{sample_name} for this dataset).}
    \item{\code{group}}{a character vector, one of 'tsne0' or 'tsne1',
    identifying which of the two t-SNE-derived lipidomic subgroups each
    sample belongs to.}
    \item{\code{pair}}{a logical vector. Not applicable for this
    (unpaired) dataset; all values are \code{NA}.}
  }
}
\source{
Ma WL, et al. Ether-lipids accumulation promotes hepatocellular carcinoma
progression linked to PPAR\enc{α}{alpha} deficiency. J Biomed Sci.
2025;32(1):89. \doi{10.1186/s12929-025-01178-y}
}
\examples{
data(group_info_multiOmics)
}
\keyword{datasets}
