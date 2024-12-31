# LipidSigR
<!-- badges: start -->
[![R >4.3](https://img.shields.io/badge/R-%3E4.3-success.svg)](https://www.r-project.org/) 
<a href='#devtools'>![installed with devtools](https://img.shields.io/badge/installed%20with-devtools-blueviolet.svg)</a>
<!-- badges: end -->

<font size=3> **"LipidSigR"** is an R package developed based on LipidSig web-based tool 
[https://lipidsig.bioinfomics.org/](https://lipidsig.bioinfomics.org/). 

This package integrates a comprehensive analysis for streamlined data mining of 
lipidomic datasets. We provide five main analysis workflows for analyzing 
two-group and multi-group data: **"Profiling,"** **"Differential Expression,"** 
**"Enrichment,"** **"Network,"** **"Machine learning,"** and **"Correlation."** 
Each section offers unique aspects to analyzing lipidome profiling data based on 
various characteristics, including lipid class, chain length, unsaturation, 
hydroxyl groups, and fatty acid composition. 
Please note that only two-group data can conduct the "Network" and "Machine learning" workflow.

## Installation
We assume that you have already installed the R program (see the R project at 
[http://www.r-project.org](http://www.r-project.org)  and are familiar with it. 
You need to have R 4.3.0 or a later version installed for running LipidSigR.

Our package is available at the github 
[https://github.com/BioinfOMICS/LipidSigR](https://github.com/BioinfOMICS/LipidSigR). 
Following are the instructions for installing our package.

```(r)
# Step 1: Install devtools
install.packages("devtools")
library(devtools)

# Step 2: Install BiocManager
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Step 3: Install LipidSigR
## Update repositories
options(repos = c(
    CRAN = "https://cloud.r-project.org/",
    BiocManager::repositories()))

## Install dependencies and package
devtools::install_github(
    "BioinfOMICS/LipidSigR", 
    dependencies = TRUE)
```

LipidSigR relies on several dependencies. If an error indicates a missing package, 
you can install the required packages using the commands below.
```(r)
# LipidSigR package depends on several packages, which can be installed using the below commands:
BiocManager::install(
    c('fgsea', 'gatom', 'mixOmics', 'S4Vectors', 'BiocGenerics', 
      'SummarizedExperiment', 'rgoslin'))

install.packages(
    c('devtools', 'magrittr', 'plotly', 'tidyverse', 'factoextra', 'ggthemes', 
      'ggforce', 'Hmisc', 'hwordcloud', 'iheatmapr', 'Rtsne', 'uwot', 
      'wordcloud', 'rsample', 'ranger', 'caret', 'yardstick', 'fastshap', 
      'SHAPforxgboost', 'visNetwork', 'tidygraph', 'ggraph'))

devtools::install_github("ctlab/mwcsr")
```

## Introduction
<font size=3> After installation, you're ready to start using LipidSigR. 
Based on functionality, LipidSigR functions can be categorized as tool functions 
and 5 analysis workflows. Below is a brief introduction to each section.
Please note that only two-group data can conduct the "Network" and "Machine learning" workflow. 

* **Tool function**: Tool functions are utility functions designed to enhance 
the convenience of conducting analyses. They include constructing input 
SummarizedExperiment objects, viewing output results, listing selectable
lipid characteristics, performing data processing, and more. 
Please read `vignette("1_tool_function")`.

* **Profiling**: The profiling workflow provides an overview of comprehensive 
analyses for you to efficiently examine data quality, the clustering of samples, 
the correlation between lipid characteristics, and the composition of lipid 
characteristics. Please read `vignette("2_profiling")`.

* **Differential expression**: The differential expression workflow integrates 
many useful lipid-focused analyses for identifying significant lipid species or 
lipid characteristics. Please read `vignette("3_de")`.

* **Enrichment**: The enrichment workflow provides two main approaches: 
'Over Representation Analysis (ORA)' and 'Lipid Set Enrichment Analysis (LSEA)' 
to illustrates significant lipid species enriched in the categories of lipid 
class and determine whether an a priori-defined set of lipids shows 
statistically significant, concordant differences between different biological states 
(e.g., phenotypes). Please read `vignette("4_enrichment")`.

* **Network**: The network workflow provides functions for generates input table 
for constructing pathway activity network, lipid reaction network and GATOM network.
*(NOTE: Only provides for two-group data.)* Please read `vignette("5_network")`.

* **Machine learning**: The machine learning workflow provides functions for 
building binary classification models and several following analyses to evaluate 
algorithm performance and identify critical lipid-related variables. *(NOTE: Only provides for two-group data.)*
Please read `vignette("6_ml")`.

* **Correlation**: The correlation workflow provides functions such as correlation 
coefficients and linear regression to analyze continuous clinical features 
correlating with lipid species or characteristics.
Please read `vignette("7_corr")`.

You can analyze data using the tool functions and the five workflow functions. 
If you need help getting started, try our case examples for two-group and 
multi-group data! Refer to `vignette("8_case_twoGroup")` for two-group data 
analysis and `vignette("9_case_multiGroup")` for multi-group data analysis. 
These case examples provide a complete tutorial, from package installation and 
input data preparation to data analysis and result visualization.

## Citation
<font size=3> You can cite the `LipidSigR` publication as follows:

> <font size=2> Chia-Hsin Liu, Pei-Chun Shen, Wen-Jen Lin, Hsiu-Cheng Liu, Meng-Hsin Tsai, 
Tzu-Ya Huang, I-Chieh Chen, Yo-Liang Lai, Yu-De Wang, Mien-Chie Hung, Wei-Chung Cheng, 
LipidSig 2.0: integrating lipid characteristic insights into advanced lipidomics data analysis, 
Nucleic Acids Research, Volume 52, Issue W1, 5 July 2024, Pages W390–W397, 
doi: [10.1093/nar/gkae335](https://doi.org/10.1093/nar/gkae335); PMID: 38709887.
