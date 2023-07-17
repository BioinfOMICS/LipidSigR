# LipidSigR

## Installation

Before the installation you should first install the following packages from Bioconductor:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("LipidSigR", version = "devel")
``` 

Or install the latest development version (on GitHub) :

``` r
devtools::install_github("BioinfOMICS/LipidSigR")
``` 

## Description 
"LipidSigR" is an R package developed based on LipidSig web-based tool (http://www.chenglab.cmu.edu.tw/lipidsig/). This package integrates a comprehensive analysis for streamlined data mining of lipidomic datasets. We provide 4 main analysis workflows, which is "Profiling", "Differential expression", "Machine learning", and "Correlation". Each section provides unique aspects to analyze the lipidome profiling data based on different characteristics including lipid class, chain length, unsaturation, hydroxyl group, and fatty acid composition.
"Profiling" provides an overview of comprehensive analyses to efficiently examine data quality, the clustering of samples, the correlation between lipid species, and the composition of lipid characteristics. "Differential expression" integrates many useful lipid-focused analyses for identifying significant lipid species or lipid characteristics. "Machine learning" provides a broad variety of feature selection methods and classifiers to build binary classification models. Furthermore, the subsequent analyses can help users to evaluate the learning algorithm's performance and to explore important lipid-related variables. "Correlation" illustrates and compares the relationships between different clinical phenotypes and lipid features.

