\name{corr_group_info}
\alias{corr_group_info}
\docType{data}
\title{
corr_group_info
}
\description{
A data frame of clinical data for correlation analysis workflow.
The data was obtained from Bowler RP et al.'s study, "Plasma Sphingolipids
Associated with Chronic Obstructive Pulmonary Disease Phenotypes." Targeted mass
spectrometry identified 69 distinct plasma sphingolipid species in 129 current
and former smokers. This cohort was analyzed to explore the associations between
plasma sphingolipids and COPD subphenotypes, including airflow obstruction,
emphysema, and frequent exacerbations.
}
\usage{data("corr_group_info")}
\format{
  A data frame with 129 observations on the following 9 variables.
  \describe{
    \item{\code{sample_name}}{a character vector}
    \item{\code{FEV1_FVC}}{a numeric vector}
    \item{\code{Emphysema}}{a numeric vector}
    \item{\code{Exacerbations}}{a numeric vector}
    \item{\code{Age}}{a numeric vector}
    \item{\code{Sex}}{a numeric vector}
    \item{\code{Smoking}}{a numeric vector}
    \item{\code{BMI}}{a numeric vector}
    \item{\code{FEV1}}{a numeric vector}
  }
}
\source{
Bowler RP, Jacobson S, Cruickshank C, Hughes GJ, Siska C, Ory DS, Petrache I, Schaffer JE, Reisdorph N, Kechris K. Plasma sphingolipids associated with chronic obstructive pulmonary disease phenotypes. Am J Respir Crit Care Med. 2015;191(3):275-84.
}
\examples{
data(corr_group_info)
}
\keyword{datasets}
