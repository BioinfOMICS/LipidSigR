\name{ml_abundance}
\alias{ml_abundance}
\docType{data}
\title{
ml_abundance
}
\description{
A data frame of lipid abundance data for machine learning analysis workflow.
The data was modified from the paper "The landscape of cancer cell line metabolism"
and divided cancer cell lines into sensitive or resistant to SCD gene knockout
evenly based on gene dependency scores (CERES).
}
\usage{data("ml_abundance")}
\format{
  A data frame with 88 observations on the following 229 variables.
  \describe{
    \item{\code{feature}}{a character vector}
    \item{\code{ACH000973}}{a numeric vector}
    \item{\code{ACH000070}}{a numeric vector}
    \item{\code{ACH000411}}{a numeric vector}
    \item{\code{ACH001306}}{a numeric vector}
    \item{\code{ACH000137}}{a numeric vector}
    \item{\code{ACH000697}}{a numeric vector}
    \item{\code{ACH000557}}{a numeric vector}
    \item{\code{ACH000838}}{a numeric vector}
    \item{\code{ACH000802}}{a numeric vector}
    \item{\code{ACH000228}}{a numeric vector}
    \item{\code{ACH000402}}{a numeric vector}
    \item{\code{ACH000832}}{a numeric vector}
    \item{\code{ACH000518}}{a numeric vector}
    \item{\code{ACH000856}}{a numeric vector}
    \item{\code{ACH000516}}{a numeric vector}
    \item{\code{ACH000264}}{a numeric vector}
    \item{\code{ACH000713}}{a numeric vector}
    \item{\code{ACH000329}}{a numeric vector}
    \item{\code{ACH000120}}{a numeric vector}
    \item{\code{ACH000458}}{a numeric vector}
    \item{\code{ACH000920}}{a numeric vector}
    \item{\code{ACH000253}}{a numeric vector}
    \item{\code{ACH000350}}{a numeric vector}
    \item{\code{ACH000256}}{a numeric vector}
    \item{\code{ACH000278}}{a numeric vector}
    \item{\code{ACH000608}}{a numeric vector}
    \item{\code{ACH000055}}{a numeric vector}
    \item{\code{ACH000095}}{a numeric vector}
    \item{\code{ACH000243}}{a numeric vector}
    \item{\code{ACH000211}}{a numeric vector}
    \item{\code{ACH000233}}{a numeric vector}
    \item{\code{ACH000258}}{a numeric vector}
    \item{\code{ACH000435}}{a numeric vector}
    \item{\code{ACH000330}}{a numeric vector}
    \item{\code{ACH000308}}{a numeric vector}
    \item{\code{ACH000936}}{a numeric vector}
    \item{\code{ACH000706}}{a numeric vector}
    \item{\code{ACH000295}}{a numeric vector}
    \item{\code{ACH000585}}{a numeric vector}
    \item{\code{ACH000906}}{a numeric vector}
    \item{\code{ACH000487}}{a numeric vector}
    \item{\code{ACH000633}}{a numeric vector}
    \item{\code{ACH000082}}{a numeric vector}
    \item{\code{ACH000098}}{a numeric vector}
    \item{\code{ACH000738}}{a numeric vector}
    \item{\code{ACH000047}}{a numeric vector}
    \item{\code{ACH000746}}{a numeric vector}
    \item{\code{ACH000485}}{a numeric vector}
    \item{\code{ACH000374}}{a numeric vector}
    \item{\code{ACH000859}}{a numeric vector}
    \item{\code{ACH000150}}{a numeric vector}
    \item{\code{ACH000840}}{a numeric vector}
    \item{\code{ACH000467}}{a numeric vector}
    \item{\code{ACH000012}}{a numeric vector}
    \item{\code{ACH000454}}{a numeric vector}
    \item{\code{ACH000190}}{a numeric vector}
    \item{\code{ACH000954}}{a numeric vector}
    \item{\code{ACH000941}}{a numeric vector}
    \item{\code{ACH000996}}{a numeric vector}
    \item{\code{ACH000831}}{a numeric vector}
    \item{\code{ACH000984}}{a numeric vector}
    \item{\code{ACH000542}}{a numeric vector}
    \item{\code{ACH000847}}{a numeric vector}
    \item{\code{ACH000393}}{a numeric vector}
    \item{\code{ACH000014}}{a numeric vector}
    \item{\code{ACH000799}}{a numeric vector}
    \item{\code{ACH000274}}{a numeric vector}
    \item{\code{ACH000054}}{a numeric vector}
    \item{\code{ACH000986}}{a numeric vector}
    \item{\code{ACH000552}}{a numeric vector}
    \item{\code{ACH000926}}{a numeric vector}
    \item{\code{ACH000118}}{a numeric vector}
    \item{\code{ACH000882}}{a numeric vector}
    \item{\code{ACH000915}}{a numeric vector}
    \item{\code{ACH000569}}{a numeric vector}
    \item{\code{ACH000476}}{a numeric vector}
    \item{\code{ACH000734}}{a numeric vector}
    \item{\code{ACH000324}}{a numeric vector}
    \item{\code{ACH000237}}{a numeric vector}
    \item{\code{ACH000132}}{a numeric vector}
    \item{\code{ACH000584}}{a numeric vector}
    \item{\code{ACH000928}}{a numeric vector}
    \item{\code{ACH000993}}{a numeric vector}
    \item{\code{ACH000711}}{a numeric vector}
    \item{\code{ACH000653}}{a numeric vector}
    \item{\code{ACH000404}}{a numeric vector}
    \item{\code{ACH000507}}{a numeric vector}
    \item{\code{ACH000167}}{a numeric vector}
    \item{\code{ACH000259}}{a numeric vector}
    \item{\code{ACH000773}}{a numeric vector}
    \item{\code{ACH000293}}{a numeric vector}
    \item{\code{ACH000969}}{a numeric vector}
    \item{\code{ACH000862}}{a numeric vector}
    \item{\code{ACH000250}}{a numeric vector}
    \item{\code{ACH000714}}{a numeric vector}
    \item{\code{ACH000576}}{a numeric vector}
    \item{\code{ACH000541}}{a numeric vector}
    \item{\code{ACH000622}}{a numeric vector}
    \item{\code{ACH000858}}{a numeric vector}
    \item{\code{ACH000265}}{a numeric vector}
    \item{\code{ACH000631}}{a numeric vector}
    \item{\code{ACH000524}}{a numeric vector}
    \item{\code{ACH000855}}{a numeric vector}
    \item{\code{ACH000809}}{a numeric vector}
    \item{\code{ACH000865}}{a numeric vector}
    \item{\code{ACH000824}}{a numeric vector}
    \item{\code{ACH000627}}{a numeric vector}
    \item{\code{ACH000471}}{a numeric vector}
    \item{\code{ACH000769}}{a numeric vector}
    \item{\code{ACH000819}}{a numeric vector}
    \item{\code{ACH000252}}{a numeric vector}
    \item{\code{ACH000390}}{a numeric vector}
    \item{\code{ACH000787}}{a numeric vector}
    \item{\code{ACH000477}}{a numeric vector}
    \item{\code{ACH000796}}{a numeric vector}
    \item{\code{ACH000621}}{a numeric vector}
    \item{\code{ACH000884}}{a numeric vector}
    \item{\code{ACH000910}}{a numeric vector}
    \item{\code{ACH000935}}{a numeric vector}
    \item{\code{ACH000881}}{a numeric vector}
    \item{\code{ACH000601}}{a numeric vector}
    \item{\code{ACH000362}}{a numeric vector}
    \item{\code{ACH000851}}{a numeric vector}
    \item{\code{ACH000335}}{a numeric vector}
    \item{\code{ACH000045}}{a numeric vector}
    \item{\code{ACH000938}}{a numeric vector}
    \item{\code{ACH000294}}{a numeric vector}
    \item{\code{ACH000510}}{a numeric vector}
    \item{\code{ACH000921}}{a numeric vector}
    \item{\code{ACH000589}}{a numeric vector}
    \item{\code{ACH000869}}{a numeric vector}
    \item{\code{ACH000766}}{a numeric vector}
    \item{\code{ACH000448}}{a numeric vector}
    \item{\code{ACH000021}}{a numeric vector}
    \item{\code{ACH000747}}{a numeric vector}
    \item{\code{ACH000781}}{a numeric vector}
    \item{\code{ACH000521}}{a numeric vector}
    \item{\code{ACH000841}}{a numeric vector}
    \item{\code{ACH000311}}{a numeric vector}
    \item{\code{ACH000785}}{a numeric vector}
    \item{\code{ACH000481}}{a numeric vector}
    \item{\code{ACH000924}}{a numeric vector}
    \item{\code{ACH000092}}{a numeric vector}
    \item{\code{ACH000800}}{a numeric vector}
    \item{\code{ACH000463}}{a numeric vector}
    \item{\code{ACH000395}}{a numeric vector}
    \item{\code{ACH000403}}{a numeric vector}
    \item{\code{ACH000427}}{a numeric vector}
    \item{\code{ACH000911}}{a numeric vector}
    \item{\code{ACH000116}}{a numeric vector}
    \item{\code{ACH000336}}{a numeric vector}
    \item{\code{ACH000013}}{a numeric vector}
    \item{\code{ACH000159}}{a numeric vector}
    \item{\code{ACH001151}}{a numeric vector}
    \item{\code{ACH000947}}{a numeric vector}
    \item{\code{ACH000663}}{a numeric vector}
    \item{\code{ACH000022}}{a numeric vector}
    \item{\code{ACH000606}}{a numeric vector}
    \item{\code{ACH000732}}{a numeric vector}
    \item{\code{ACH000735}}{a numeric vector}
    \item{\code{ACH000937}}{a numeric vector}
    \item{\code{ACH000307}}{a numeric vector}
    \item{\code{ACH000320}}{a numeric vector}
    \item{\code{ACH000169}}{a numeric vector}
    \item{\code{ACH000041}}{a numeric vector}
    \item{\code{ACH000960}}{a numeric vector}
    \item{\code{ACH000791}}{a numeric vector}
    \item{\code{ACH000261}}{a numeric vector}
    \item{\code{ACH000505}}{a numeric vector}
    \item{\code{ACH000943}}{a numeric vector}
    \item{\code{ACH000473}}{a numeric vector}
    \item{\code{ACH000242}}{a numeric vector}
    \item{\code{ACH000410}}{a numeric vector}
    \item{\code{ACH000839}}{a numeric vector}
    \item{\code{ACH000782}}{a numeric vector}
    \item{\code{ACH000376}}{a numeric vector}
    \item{\code{ACH000099}}{a numeric vector}
    \item{\code{ACH000748}}{a numeric vector}
    \item{\code{ACH000017}}{a numeric vector}
    \item{\code{ACH000361}}{a numeric vector}
    \item{\code{ACH000822}}{a numeric vector}
    \item{\code{ACH000665}}{a numeric vector}
    \item{\code{ACH000363}}{a numeric vector}
    \item{\code{ACH000260}}{a numeric vector}
    \item{\code{ACH000366}}{a numeric vector}
    \item{\code{ACH000811}}{a numeric vector}
    \item{\code{ACH000974}}{a numeric vector}
    \item{\code{ACH000932}}{a numeric vector}
    \item{\code{ACH000549}}{a numeric vector}
    \item{\code{ACH000302}}{a numeric vector}
    \item{\code{ACH000483}}{a numeric vector}
    \item{\code{ACH000623}}{a numeric vector}
    \item{\code{ACH000907}}{a numeric vector}
    \item{\code{ACH000517}}{a numeric vector}
    \item{\code{ACH000500}}{a numeric vector}
    \item{\code{ACH000683}}{a numeric vector}
    \item{\code{ACH000537}}{a numeric vector}
    \item{\code{ACH000460}}{a numeric vector}
    \item{\code{ACH000338}}{a numeric vector}
    \item{\code{ACH000365}}{a numeric vector}
    \item{\code{ACH000953}}{a numeric vector}
    \item{\code{ACH000677}}{a numeric vector}
    \item{\code{ACH000820}}{a numeric vector}
    \item{\code{ACH000651}}{a numeric vector}
    \item{\code{ACH000018}}{a numeric vector}
    \item{\code{ACH000571}}{a numeric vector}
    \item{\code{ACH000720}}{a numeric vector}
    \item{\code{ACH000318}}{a numeric vector}
    \item{\code{ACH000488}}{a numeric vector}
    \item{\code{ACH000917}}{a numeric vector}
    \item{\code{ACH000408}}{a numeric vector}
    \item{\code{ACH000452}}{a numeric vector}
    \item{\code{ACH000397}}{a numeric vector}
    \item{\code{ACH000387}}{a numeric vector}
    \item{\code{ACH000146}}{a numeric vector}
    \item{\code{ACH000885}}{a numeric vector}
    \item{\code{ACH000459}}{a numeric vector}
    \item{\code{ACH000317}}{a numeric vector}
    \item{\code{ACH000495}}{a numeric vector}
    \item{\code{ACH000364}}{a numeric vector}
    \item{\code{ACH000579}}{a numeric vector}
    \item{\code{ACH000522}}{a numeric vector}
    \item{\code{ACH000545}}{a numeric vector}
    \item{\code{ACH000484}}{a numeric vector}
    \item{\code{ACH000827}}{a numeric vector}
    \item{\code{ACH000332}}{a numeric vector}
    \item{\code{ACH000469}}{a numeric vector}
    \item{\code{ACH000570}}{a numeric vector}
  }
}
\source{
Li H, Ning S, Ghandi M, Kryukov GV, Gopal S, Deik A, et al. The landscape of
cancer cell line metabolism. Nature medicine. 2019;25(5):850.
}
\examples{
data(ml_abundance)
}
\keyword{datasets}
