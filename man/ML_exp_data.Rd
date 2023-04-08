\name{ML_exp_data}
\alias{ML_exp_data}
\docType{data}
\title{
ML_exp_data
}
\description{
A data frame of lipid expression data for machine learning analysis workflow.
The data was modified from the paper "The landscape of cancer cell line
metabolism" and divided cancer cell lines into sensitive or resistant to SCD
gene knockout evenly based on gene dependency scores (CERES). 89 lipid species
were used to predict 228 cancer cell lines with a label of 0 (sensitive) or 1
(resistant) by a binary classifier.
}
\usage{data("ML_exp_data")}
\value{
  A data frame with 89 observations on the following 229 variables.
  \describe{
    \item{\code{feature}}{a character vector}
    \item{\code{ACH_000973}}{a numeric vector}
    \item{\code{ACH_000070}}{a numeric vector}
    \item{\code{ACH_000411}}{a numeric vector}
    \item{\code{ACH_001306}}{a numeric vector}
    \item{\code{ACH_000137}}{a numeric vector}
    \item{\code{ACH_000697}}{a numeric vector}
    \item{\code{ACH_000557}}{a numeric vector}
    \item{\code{ACH_000838}}{a numeric vector}
    \item{\code{ACH_000802}}{a numeric vector}
    \item{\code{ACH_000228}}{a numeric vector}
    \item{\code{ACH_000402}}{a numeric vector}
    \item{\code{ACH_000832}}{a numeric vector}
    \item{\code{ACH_000518}}{a numeric vector}
    \item{\code{ACH_000856}}{a numeric vector}
    \item{\code{ACH_000516}}{a numeric vector}
    \item{\code{ACH_000264}}{a numeric vector}
    \item{\code{ACH_000713}}{a numeric vector}
    \item{\code{ACH_000329}}{a numeric vector}
    \item{\code{ACH_000120}}{a numeric vector}
    \item{\code{ACH_000458}}{a numeric vector}
    \item{\code{ACH_000920}}{a numeric vector}
    \item{\code{ACH_000253}}{a numeric vector}
    \item{\code{ACH_000350}}{a numeric vector}
    \item{\code{ACH_000256}}{a numeric vector}
    \item{\code{ACH_000278}}{a numeric vector}
    \item{\code{ACH_000608}}{a numeric vector}
    \item{\code{ACH_000055}}{a numeric vector}
    \item{\code{ACH_000095}}{a numeric vector}
    \item{\code{ACH_000243}}{a numeric vector}
    \item{\code{ACH_000211}}{a numeric vector}
    \item{\code{ACH_000233}}{a numeric vector}
    \item{\code{ACH_000258}}{a numeric vector}
    \item{\code{ACH_000435}}{a numeric vector}
    \item{\code{ACH_000330}}{a numeric vector}
    \item{\code{ACH_000308}}{a numeric vector}
    \item{\code{ACH_000936}}{a numeric vector}
    \item{\code{ACH_000706}}{a numeric vector}
    \item{\code{ACH_000295}}{a numeric vector}
    \item{\code{ACH_000585}}{a numeric vector}
    \item{\code{ACH_000906}}{a numeric vector}
    \item{\code{ACH_000487}}{a numeric vector}
    \item{\code{ACH_000633}}{a numeric vector}
    \item{\code{ACH_000082}}{a numeric vector}
    \item{\code{ACH_000098}}{a numeric vector}
    \item{\code{ACH_000738}}{a numeric vector}
    \item{\code{ACH_000047}}{a numeric vector}
    \item{\code{ACH_000746}}{a numeric vector}
    \item{\code{ACH_000485}}{a numeric vector}
    \item{\code{ACH_000374}}{a numeric vector}
    \item{\code{ACH_000859}}{a numeric vector}
    \item{\code{ACH_000150}}{a numeric vector}
    \item{\code{ACH_000840}}{a numeric vector}
    \item{\code{ACH_000467}}{a numeric vector}
    \item{\code{ACH_000012}}{a numeric vector}
    \item{\code{ACH_000454}}{a numeric vector}
    \item{\code{ACH_000190}}{a numeric vector}
    \item{\code{ACH_000954}}{a numeric vector}
    \item{\code{ACH_000941}}{a numeric vector}
    \item{\code{ACH_000996}}{a numeric vector}
    \item{\code{ACH_000831}}{a numeric vector}
    \item{\code{ACH_000984}}{a numeric vector}
    \item{\code{ACH_000542}}{a numeric vector}
    \item{\code{ACH_000847}}{a numeric vector}
    \item{\code{ACH_000393}}{a numeric vector}
    \item{\code{ACH_000014}}{a numeric vector}
    \item{\code{ACH_000799}}{a numeric vector}
    \item{\code{ACH_000274}}{a numeric vector}
    \item{\code{ACH_000054}}{a numeric vector}
    \item{\code{ACH_000986}}{a numeric vector}
    \item{\code{ACH_000552}}{a numeric vector}
    \item{\code{ACH_000926}}{a numeric vector}
    \item{\code{ACH_000118}}{a numeric vector}
    \item{\code{ACH_000882}}{a numeric vector}
    \item{\code{ACH_000915}}{a numeric vector}
    \item{\code{ACH_000569}}{a numeric vector}
    \item{\code{ACH_000476}}{a numeric vector}
    \item{\code{ACH_000734}}{a numeric vector}
    \item{\code{ACH_000324}}{a numeric vector}
    \item{\code{ACH_000237}}{a numeric vector}
    \item{\code{ACH_000132}}{a numeric vector}
    \item{\code{ACH_000584}}{a numeric vector}
    \item{\code{ACH_000928}}{a numeric vector}
    \item{\code{ACH_000993}}{a numeric vector}
    \item{\code{ACH_000711}}{a numeric vector}
    \item{\code{ACH_000653}}{a numeric vector}
    \item{\code{ACH_000404}}{a numeric vector}
    \item{\code{ACH_000507}}{a numeric vector}
    \item{\code{ACH_000167}}{a numeric vector}
    \item{\code{ACH_000259}}{a numeric vector}
    \item{\code{ACH_000773}}{a numeric vector}
    \item{\code{ACH_000293}}{a numeric vector}
    \item{\code{ACH_000969}}{a numeric vector}
    \item{\code{ACH_000862}}{a numeric vector}
    \item{\code{ACH_000250}}{a numeric vector}
    \item{\code{ACH_000714}}{a numeric vector}
    \item{\code{ACH_000576}}{a numeric vector}
    \item{\code{ACH_000541}}{a numeric vector}
    \item{\code{ACH_000622}}{a numeric vector}
    \item{\code{ACH_000858}}{a numeric vector}
    \item{\code{ACH_000265}}{a numeric vector}
    \item{\code{ACH_000631}}{a numeric vector}
    \item{\code{ACH_000524}}{a numeric vector}
    \item{\code{ACH_000855}}{a numeric vector}
    \item{\code{ACH_000809}}{a numeric vector}
    \item{\code{ACH_000865}}{a numeric vector}
    \item{\code{ACH_000824}}{a numeric vector}
    \item{\code{ACH_000627}}{a numeric vector}
    \item{\code{ACH_000471}}{a numeric vector}
    \item{\code{ACH_000769}}{a numeric vector}
    \item{\code{ACH_000819}}{a numeric vector}
    \item{\code{ACH_000252}}{a numeric vector}
    \item{\code{ACH_000390}}{a numeric vector}
    \item{\code{ACH_000787}}{a numeric vector}
    \item{\code{ACH_000477}}{a numeric vector}
    \item{\code{ACH_000796}}{a numeric vector}
    \item{\code{ACH_000621}}{a numeric vector}
    \item{\code{ACH_000884}}{a numeric vector}
    \item{\code{ACH_000910}}{a numeric vector}
    \item{\code{ACH_000935}}{a numeric vector}
    \item{\code{ACH_000881}}{a numeric vector}
    \item{\code{ACH_000601}}{a numeric vector}
    \item{\code{ACH_000362}}{a numeric vector}
    \item{\code{ACH_000851}}{a numeric vector}
    \item{\code{ACH_000335}}{a numeric vector}
    \item{\code{ACH_000045}}{a numeric vector}
    \item{\code{ACH_000938}}{a numeric vector}
    \item{\code{ACH_000294}}{a numeric vector}
    \item{\code{ACH_000510}}{a numeric vector}
    \item{\code{ACH_000921}}{a numeric vector}
    \item{\code{ACH_000589}}{a numeric vector}
    \item{\code{ACH_000869}}{a numeric vector}
    \item{\code{ACH_000766}}{a numeric vector}
    \item{\code{ACH_000448}}{a numeric vector}
    \item{\code{ACH_000021}}{a numeric vector}
    \item{\code{ACH_000747}}{a numeric vector}
    \item{\code{ACH_000781}}{a numeric vector}
    \item{\code{ACH_000521}}{a numeric vector}
    \item{\code{ACH_000841}}{a numeric vector}
    \item{\code{ACH_000311}}{a numeric vector}
    \item{\code{ACH_000785}}{a numeric vector}
    \item{\code{ACH_000481}}{a numeric vector}
    \item{\code{ACH_000924}}{a numeric vector}
    \item{\code{ACH_000092}}{a numeric vector}
    \item{\code{ACH_000800}}{a numeric vector}
    \item{\code{ACH_000463}}{a numeric vector}
    \item{\code{ACH_000395}}{a numeric vector}
    \item{\code{ACH_000403}}{a numeric vector}
    \item{\code{ACH_000427}}{a numeric vector}
    \item{\code{ACH_000911}}{a numeric vector}
    \item{\code{ACH_000116}}{a numeric vector}
    \item{\code{ACH_000336}}{a numeric vector}
    \item{\code{ACH_000013}}{a numeric vector}
    \item{\code{ACH_000159}}{a numeric vector}
    \item{\code{ACH_001151}}{a numeric vector}
    \item{\code{ACH_000947}}{a numeric vector}
    \item{\code{ACH_000663}}{a numeric vector}
    \item{\code{ACH_000022}}{a numeric vector}
    \item{\code{ACH_000606}}{a numeric vector}
    \item{\code{ACH_000732}}{a numeric vector}
    \item{\code{ACH_000735}}{a numeric vector}
    \item{\code{ACH_000937}}{a numeric vector}
    \item{\code{ACH_000307}}{a numeric vector}
    \item{\code{ACH_000320}}{a numeric vector}
    \item{\code{ACH_000169}}{a numeric vector}
    \item{\code{ACH_000041}}{a numeric vector}
    \item{\code{ACH_000960}}{a numeric vector}
    \item{\code{ACH_000791}}{a numeric vector}
    \item{\code{ACH_000261}}{a numeric vector}
    \item{\code{ACH_000505}}{a numeric vector}
    \item{\code{ACH_000943}}{a numeric vector}
    \item{\code{ACH_000473}}{a numeric vector}
    \item{\code{ACH_000242}}{a numeric vector}
    \item{\code{ACH_000410}}{a numeric vector}
    \item{\code{ACH_000839}}{a numeric vector}
    \item{\code{ACH_000782}}{a numeric vector}
    \item{\code{ACH_000376}}{a numeric vector}
    \item{\code{ACH_000099}}{a numeric vector}
    \item{\code{ACH_000748}}{a numeric vector}
    \item{\code{ACH_000017}}{a numeric vector}
    \item{\code{ACH_000361}}{a numeric vector}
    \item{\code{ACH_000822}}{a numeric vector}
    \item{\code{ACH_000665}}{a numeric vector}
    \item{\code{ACH_000363}}{a numeric vector}
    \item{\code{ACH_000260}}{a numeric vector}
    \item{\code{ACH_000366}}{a numeric vector}
    \item{\code{ACH_000811}}{a numeric vector}
    \item{\code{ACH_000974}}{a numeric vector}
    \item{\code{ACH_000932}}{a numeric vector}
    \item{\code{ACH_000549}}{a numeric vector}
    \item{\code{ACH_000302}}{a numeric vector}
    \item{\code{ACH_000483}}{a numeric vector}
    \item{\code{ACH_000623}}{a numeric vector}
    \item{\code{ACH_000907}}{a numeric vector}
    \item{\code{ACH_000517}}{a numeric vector}
    \item{\code{ACH_000500}}{a numeric vector}
    \item{\code{ACH_000683}}{a numeric vector}
    \item{\code{ACH_000537}}{a numeric vector}
    \item{\code{ACH_000460}}{a numeric vector}
    \item{\code{ACH_000338}}{a numeric vector}
    \item{\code{ACH_000365}}{a numeric vector}
    \item{\code{ACH_000953}}{a numeric vector}
    \item{\code{ACH_000677}}{a numeric vector}
    \item{\code{ACH_000820}}{a numeric vector}
    \item{\code{ACH_000651}}{a numeric vector}
    \item{\code{ACH_000018}}{a numeric vector}
    \item{\code{ACH_000571}}{a numeric vector}
    \item{\code{ACH_000720}}{a numeric vector}
    \item{\code{ACH_000318}}{a numeric vector}
    \item{\code{ACH_000488}}{a numeric vector}
    \item{\code{ACH_000917}}{a numeric vector}
    \item{\code{ACH_000408}}{a numeric vector}
    \item{\code{ACH_000452}}{a numeric vector}
    \item{\code{ACH_000397}}{a numeric vector}
    \item{\code{ACH_000387}}{a numeric vector}
    \item{\code{ACH_000146}}{a numeric vector}
    \item{\code{ACH_000885}}{a numeric vector}
    \item{\code{ACH_000459}}{a numeric vector}
    \item{\code{ACH_000317}}{a numeric vector}
    \item{\code{ACH_000495}}{a numeric vector}
    \item{\code{ACH_000364}}{a numeric vector}
    \item{\code{ACH_000579}}{a numeric vector}
    \item{\code{ACH_000522}}{a numeric vector}
    \item{\code{ACH_000545}}{a numeric vector}
    \item{\code{ACH_000484}}{a numeric vector}
    \item{\code{ACH_000827}}{a numeric vector}
    \item{\code{ACH_000332}}{a numeric vector}
    \item{\code{ACH_000469}}{a numeric vector}
    \item{\code{ACH_000570}}{a numeric vector}
  }
}
\source{
Li H, Ning S, Ghandi M, Kryukov GV, Gopal S, Deik A, et al. The landscape of
cancer cell line metabolism. Nature medicine. 2019;25(5):850.
}
\examples{
data(ML_exp_data)
}
\keyword{datasets}
