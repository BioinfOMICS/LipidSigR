# LipidSigR
<!-- badges: start -->
[![R >4.2](https://img.shields.io/badge/R-%3E4.2-success.svg)](https://www.r-project.org/) 
<a href='#devtools'>![installed with devtools](https://img.shields.io/badge/installed%20with-devtools-blueviolet.svg)</a>
<!-- badges: end -->

<font size=3> **"LipidSigR"** is an R package developed based on LipidSig web-based tool 
[https://lipidsig.bioinfomics.org/](https://lipidsig.bioinfomics.org/). 

This package integrates a comprehensive analysis for streamlined data mining of 
lipidomic datasets. We provide 4 main analysis workflows, which is 
**"Profiling"**, **"Differential expression"**, **"Enrichment"**, and 
**"Network"**. Each section provides unique aspects to analyze the lipidome 
profiling data based on different characteristics including lipid class, 
chain length, unsaturation, hydroxyl group, and fatty acid composition.

> [!IMPORTANT]
> ### $\textcolor{purple}{A\ full\ tutorial\ is\ currently\ being\ prepared\. Please\ stay\ tuned\ for\ its\ release\!  }$

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
BiocManager::install(c('fgsea', 'gatom', 'mixOmics', 
                       'S4Vectors', 'SummarizedExperiment', 'rgoslin'))

install.packages(c('devtools', 'magrittr', 'plotly', 'tidyverse'))
devtools::install_github("ctlab/mwcsr")
```

### 2. Clone Github and install locally
```(r)
git clone https://github.com/BioinfOMICS/LipidSigR.git
R CMD build LipidSigR
R CMD INSTALL LipidSigR_0.7.0.tar.gz
```
## Introduction
<font size=3> After the installation, you are ready to start using LipidSigR. There are 4 
workflows provided for analysis and each of them is introduced briefly below. 

* **Profiling**: The profiling workflow provides an overview of comprehensive 
analyses for you to efficiently examine data quality, the clustering of samples, 
the correlation between lipid species, and the composition of lipid characteristics.

* **Differential expression**: The differential expression workflow integrates 
many useful lipid-focused analyses for identifying significant lipid species or 
lipid characteristics.

* **Enrichment**: The enrichment workflow provides two main approaches: 
'Over Representation Analysis (ORA)' and 'Lipid Set Enrichment Analysis (LSEA)' 
to illustrates significant lipid species enriched in the categories of lipid 
class and determine whether an a priori-defined set of lipids shows 
statistically significant, concordant differences between two biological states 
(e.g., phenotypes).

* **Network**: The network wrokflow provides functions for generates input table 
for constructing pathway activity network, lipid reaction network and GATOM network.

If you need help starting the analysis, you can begin by reading our case 
example for two-group and multi-group data.
$\textcolor{red}{(\ A\ full\ tutorial\ is\ currently\ being\ prepared\. Please\ stay\ tuned\ for\ its\ release\ !\) }$

## Citation
<font size=3> You can cite the `LipidSigR` publication as follows:

> <font size=2> Chia-Hsin Liu, Pei-Chun Shen, Wen-Jen Lin, Hsiu-Cheng Liu, Meng-Hsin Tsai, 
Tzu-Ya Huang, I-Chieh Chen, Yo-Liang Lai, Yu-De Wang, Mien-Chie Hung, Wei-Chung Cheng, 
LipidSig 2.0: integrating lipid characteristic insights into advanced lipidomics data analysis, 
Nucleic Acids Research, Volume 52, Issue W1, 5 July 2024, Pages W390–W397, 
doi: [10.1093/nar/gkae335](https://doi.org/10.1093/nar/gkae335); PMID: 38709887.
