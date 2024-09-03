# LipidSigR
<!-- badges: start -->
[![R >4.2](https://img.shields.io/badge/R-%3E4.2-success.svg)](https://www.r-project.org/) 
<a href='#devtools'>![installed with devtools](https://img.shields.io/badge/installed%20with-devtools-blueviolet.svg)</a>
<!-- badges: end -->

<font size=3> **"LipidSigR"** is an R package developed based on LipidSig web-based tool 
[https://lipidsig.bioinfomics.org/](https://lipidsig.bioinfomics.org/). 

This package integrates a comprehensive analysis for streamlined data mining of 
lipidomic datasets. We provide four main analysis workflows for analyzing 
two-group and multi-group data: **"Profiling,"** **"Differential Expression,"** 
**"Enrichment,"** and **"Network."** Each section offers unique aspects to 
analyzing lipidome profiling data based on various characteristics, including 
lipid class, chain length, unsaturation, hydroxyl groups, and fatty acid composition. 
Please note that only two-group data can conduct the "Network" workflow.

> [!IMPORTANT]
> * #### For instructions and details on LipidSigR, please refer to https://lipidsig.bioinfomics.org/lipidsigr/.

## Installation
We assume that you have already installed the R program (see the R project at 
[http://www.r-project.org](http://www.r-project.org)  and are familiar with it. 
You need to have R 4.2.0 or a later version installed for running LipidSigR.

Our package is available at the github 
[https://github.com/BioinfOMICS/LipidSigR](https://github.com/BioinfOMICS/LipidSigR). 
There are 2 recommended ways to install our package.

### 1. Install the package directly from github by using the devtools package
```(r)
# Step 1: Install devtools
install.packages("devtools")
library(devtools)

# Step 2: Install LipidSigR
devtools::install_github("BioinfOMICS/LipidSigR")

# LipidSigR package depends on several packages, which can be installed using the below commands:
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c('fgsea', 'gatom', 'mixOmics', 'S4Vectors', 'SummarizedExperiment', 'rgoslin'))

install.packages(c('devtools', 'magrittr', 'plotly', 'tidyverse'))
devtools::install_github("ctlab/mwcsr")
```

### 2. Clone Github and install locally
```(r)
git clone https://github.com/BioinfOMICS/LipidSigR.git
R CMD build LipidSigR
R CMD INSTALL LipidSigR_0.7.0.tar.gz
```

## Citation
<font size=3> You can cite the `LipidSigR` publication as follows:

> <font size=2> Chia-Hsin Liu, Pei-Chun Shen, Wen-Jen Lin, Hsiu-Cheng Liu, Meng-Hsin Tsai, 
Tzu-Ya Huang, I-Chieh Chen, Yo-Liang Lai, Yu-De Wang, Mien-Chie Hung, Wei-Chung Cheng, 
LipidSig 2.0: integrating lipid characteristic insights into advanced lipidomics data analysis, 
Nucleic Acids Research, Volume 52, Issue W1, 5 July 2024, Pages W390–W397, 
doi: [10.1093/nar/gkae335](https://doi.org/10.1093/nar/gkae335); PMID: 38709887.
