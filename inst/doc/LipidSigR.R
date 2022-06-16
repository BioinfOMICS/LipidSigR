### R code from vignette source 'LipidSigR.Rnw'

###################################################
### code chunk number 1: style-Sweave
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: LipidSigR.Rnw:21-24
###################################################
knitr::opts_chunk$set(warning=FALSE)
knitr::opts_chunk$set(fig.width = 6, fig.height = 5)
knitr::opts_chunk$set(fig.align ="center")


###################################################
### code chunk number 3: install_Bioconductor (eval = FALSE)
###################################################
## if (!requireNamespace("BiocManager", quietly=TRUE))
##     install.packages("BiocManager")
## BiocManager::install()


###################################################
### code chunk number 4: install_package (eval = FALSE)
###################################################
## if (!requireNamespace("BiocManager", quietly=TRUE))
##     install.packages("BiocManager")
## BiocManager::install("LipidSigR")


###################################################
### code chunk number 5: load_package
###################################################
library(LipidSigR)


###################################################
### code chunk number 6: LipidSigR.Rnw:107-111
###################################################
library(plotly)
library(pathview)
library(SHAPforxgboost)
library(visNetwork)


###################################################
### code chunk number 7: load_profiling_data
###################################################
# clears all objects from workspace
rm(list = ls())

# lipid expression data
data("profiling_exp_data")
exp_data <- profiling_exp_data
head(exp_data[, 1:5], 5)

# lipid characteristics table (only use in Section 3.5)
data("profiling_lipid_char_table")
lipid_char_table <- profiling_lipid_char_table
head(lipid_char_table[, 1:4], 5)


###################################################
### code chunk number 8: data_process
###################################################
# lipid expression data
head(exp_data[, 1:5], 5)
# data processing of exp_data
exp_transform_table <- data_process(exp_data, exclude_var_missing=TRUE,
                                    missing_pct_limit=50,
                                    replace_zero=TRUE, zero2what='min',
                                    xmin=0.5, replace_NA=TRUE,
                                    NA2what='min', ymin=0.5,
                                    pct_transform=TRUE,
                                    data_transform=TRUE, trans_type='log',
                                    centering=FALSE, scaling=FALSE)
# exp_data after data processing
head(exp_transform_table[, 1:5], 5)


###################################################
### code chunk number 9: Profiling: cross-sample variability
###################################################
# conduct profiling
profiling_result <- exp_profiling(exp_data)


###################################################
### code chunk number 10: numExpLipid
###################################################
# view result: histograms (number of expressed lipids)
profiling_result$i.expr.lip


###################################################
### code chunk number 11: lipid_amount
###################################################
# view result: histogram (total amount of lipid)
profiling_result$i.p.amount


###################################################
### code chunk number 12: desity_exp_distri
###################################################
# view result: density plot of expression distribution
profiling_result$p.hist.value


###################################################
### code chunk number 13: Profiling: dimensionality reductionPCA
###################################################
# data processing of exp_data
exp_transform_table <- data_process(exp_data, exclude_var_missing=TRUE,
                                    missing_pct_limit=50,
                                    replace_zero=TRUE, zero2what='min',
                                    xmin=0.5, replace_NA=TRUE,
                                    NA2what='min', ymin=0.5,
                                    pct_transform=TRUE,
                                    data_transform=TRUE, trans_type='log',
                                    centering=FALSE, scaling=FALSE)
# conduct PCA
PCA_result <- PCA(exp_transform_table,
                  group_info = NULL, sig_feature = NULL,
                  scaling=TRUE, centering=TRUE, cluster_method='kmeans',
                  group_num=2, var1 = NULL, var2 = NULL,
                  insert_ref_group=NULL, ref_group=NULL,
                  n_PC=c(1,2), top_n_feature=10)

# view result: PCA prcomp
head(PCA_result[[1]], 1)


###################################################
### code chunk number 14: LipidSigR.Rnw:235-237
###################################################
  # view result: PCA plot
  PCA_result[[4]]


###################################################
### code chunk number 15: LipidSigR.Rnw:246-248
###################################################
# view result: scree plot of top 10 principle components
PCA_result[[5]]


###################################################
### code chunk number 16: LipidSigR.Rnw:259-263
###################################################
# view result: data frame of PCA rotated data
head(PCA_result[[2]][,1:5], 5)
# view result: data frame of PCA contribution table
head(PCA_result[[3]][,1:5], 5)


###################################################
### code chunk number 17: LipidSigR.Rnw:266-268
###################################################
# view result: correlation circle plot of PCA variables
PCA_result[[7]]


###################################################
### code chunk number 18: LipidSigR.Rnw:279-281
###################################################
# view result: bar plot of contribution of top 10 features
PCA_result[[6]]


###################################################
### code chunk number 19: Profiling: dimensionality reduction - t-SNE
###################################################
# data processing of exp_data
exp_transform_table <- data_process(exp_data, exclude_var_missing=TRUE,
                                    missing_pct_limit=50,
                                    replace_zero=TRUE, zero2what='min',
                                    xmin=0.5, replace_NA=TRUE,
                                    NA2what='min', ymin=0.5,
                                    pct_transform=TRUE,
                                    data_transform=TRUE, trans_type='log',
                                    centering=FALSE, scaling=FALSE)
# conduct t-SNE
tsne_result <- tsne(exp_transform_table, group_info = NULL,
                    sig_feature = NULL, pca=TRUE, perplexity=5,
                    max_iter=500, cluster_method='kmeans',
                    group_num=2, var1 = 'euclidean', var2 = NULL,
                    insert_ref_group = NULL, ref_group = NULL)

# view result: data frame of t-SNE data
head(tsne_result[[1]], 5)


###################################################
### code chunk number 20: LipidSigR.Rnw:313-315
###################################################
# view result: t-SNE plot
tsne_result[[2]]


###################################################
### code chunk number 21: Profiling: dimensionality reduction - UMAP
###################################################
# data processing of exp_data
exp_transform_table <- data_process(exp_data, trans_type='log',
                                    exclude_var_missing=TRUE,
                                    missing_pct_limit=50,
                                    replace_zero=TRUE, replace_NA=TRUE,
                                    zero2what='min', xmin=0.5,
                                    NA2what='min', ymin=0.5,
                                    pct_transform=TRUE, centering=FALSE,
                                    data_transform=TRUE, scaling=FALSE )
# conduct UMAP
UMAP_result <- UMAP(exp_transform_table, group_info=NULL,
                    sig_feature=NULL, n_neighbors=15,
                    scale=TRUE, metric='euclidean', group_num=2,
                    cluster_method='kmeans', var1=NULL, var2=NULL,
                    insert_ref_group=NULL, ref_group=NULL)

# view result: data frame of UMAP data
head(UMAP_result[[1]], 5)


###################################################
### code chunk number 22: LipidSigR.Rnw:350-352
###################################################
# view result: UMAP plot
UMAP_result[[2]]


###################################################
### code chunk number 23: Profiling: correlation heatmap
###################################################
# data processing of exp_data
exp_transform <- data_process(exp_data, exclude_var_missing=TRUE,
                              missing_pct_limit=50, replace_zero=TRUE,
                              zero2what='min', xmin=0.5, replace_NA=TRUE,
                              NA2what='min', ymin=0.5, pct_transform=TRUE,
                              data_transform=TRUE, trans_type='log',
                              centering=FALSE, scaling=FALSE)
# correlation calculation
corr_result <- corr_heatmap(exp_transform, corr_method="pearson",
                            distfun="maximum", hclustfun="average")

# view result: matrix of correlation coefficients
head(corr_result$sample_corr_coef[, 1:5], 5)

# view result: matrix of correlation p-value
head(corr_result$sample_corr_p[, 1:4], 5)

# view result: matrix of reorder sample correlation
head(corr_result$reorder_sample_corr_coef[, 1:3], 5)

# view result: matrix of correlation coefficients between lipids
head(corr_result$lipids_corr_coef[, 1:5], 5)

# view result: matrix of correlation p-value between lipids
head(corr_result$lipid_corr_p[, 1:4], 5)


###################################################
### code chunk number 24: LipidSigR.Rnw:390-392
###################################################
# view result: sample-sample heatmap
corr_result$sample_hm


###################################################
### code chunk number 25: LipidSigR.Rnw:402-404
###################################################
# view result: lipid-lipid correlations heatmap
corr_result$lipids_hm


###################################################
### code chunk number 26: Profiling: lipid characteristics
###################################################
# lipid characteristic
char_var <- colnames(lipid_char_table)[-1]
# calculate lipid expression of selected characteristic
compo_result <- exp_compo_by_lipidinfo(exp_data, lipid_char_table, char_var[1])


###################################################
### code chunk number 27: LipidSigR.Rnw:422-424
###################################################
# view result: bar plot
compo_result$p.barplot.p


###################################################
### code chunk number 28: LipidSigR.Rnw:432-434
###################################################
# view result: stacked horizontal bar chart
compo_result$p.compos


###################################################
### code chunk number 29: load_DE_data
###################################################
# clears all objects from workspace
rm(list = ls())

# lipid expression data
data("DE_exp_data")
exp_data <- DE_exp_data
head(exp_data[, 1:5], 5)

# lipid characteristics table
data("DE_lipid_char_table")
lipid_char_table <- DE_lipid_char_table
head(lipid_char_table[, 1:4], 5)

# group information table
data("DE_group_info")
group_info <- DE_group_info
head(group_info, 5)


###################################################
### code chunk number 30: DE_data_process
###################################################
# lipid expression data
head(exp_data[, 1:5], 5)
# data processing of exp_data
exp_transform_table <- data_process(exp_data, exclude_var_missing=TRUE,
                                    missing_pct_limit=50,
                                    replace_zero=TRUE, zero2what='min',
                                    xmin=0.5, replace_NA=TRUE,
                                    NA2what='min', ymin=0.5,
                                    pct_transform=TRUE,
                                    data_transform=TRUE, trans_type='log',
                                    centering=FALSE, scaling=FALSE)
head(exp_transform_table[, 1:5], 5)


###################################################
### code chunk number 31: DE_lipid species: differentially expressed analysis
###################################################
# data processing of exp_data (without log10 transformation)
exp_transform_non_log <- data_process(exp_data, exclude_var_missing=TRUE,
                                      missing_pct_limit=50, replace_zero=TRUE,
                                      zero2what='min', xmin=0.5,
                                      replace_NA=TRUE, NA2what='min',
                                      ymin=0.5, pct_transform=TRUE,
                                      data_transform=FALSE, trans_type='log',
                                      centering=FALSE, scaling=FALSE)
# conduct differentially expressed analysis of lipid species
DE_species_result <- DE_species_2(exp_transform_non_log,
                                  data_transform=TRUE,
                                  group_info = group_info,
                                  paired=FALSE, test='t.test',
                                  adjust_p_method='BH', sig_stat='p.adj',
                                  sig_pvalue=0.05, sig_FC=2)

# view result: data frame of lipid species
head(DE_species_result$DE_species_table_all[, 1:5], 5)

# view result: data frame of significant lipid species
head(DE_species_result$DE_species_table_sig [, 1:5], 5)


###################################################
### code chunk number 32: LipidSigR.Rnw:540-542
###################################################
# view result: lollipop chart
DE_species_result$DE_species_dotchart_sig


###################################################
### code chunk number 33: LipidSigR.Rnw:550-552
###################################################
# view result: MA plot
DE_species_result$DE_species_maplot


###################################################
### code chunk number 34: LipidSigR.Rnw:561-563
###################################################
# view result: MA plot
DE_species_result$DE_species_volcano


###################################################
### code chunk number 35: DE_lipid species: dimensionality reduction - PCA
###################################################
# data processing of exp_data
exp_transform_table <- data_process(exp_data, exclude_var_missing=TRUE,
                                    missing_pct_limit=50,
                                    replace_zero=TRUE, zero2what='min',
                                    xmin=0.5, replace_NA=TRUE,
                                    NA2what='min', ymin=0.5,
                                    pct_transform=TRUE,
                                    data_transform=TRUE,
                                    trans_type='log',
                                    centering=FALSE, scaling=FALSE)
# data processing of exp_data (without log10 transformation)
exp_transform_non_log <- data_process(exp_data,
                                      exclude_var_missing=TRUE,
                                      missing_pct_limit=50,
                                      replace_zero=TRUE,
                                      zero2what='min', xmin=0.5,
                                      replace_NA=TRUE, NA2what='min',
                                      ymin=0.5, pct_transform=TRUE,
                                      data_transform=FALSE,
                                      trans_type='log',
                                      centering=FALSE, scaling=FALSE)
# filter significant lipid
DE_species_table_sig <- DE_species_2(exp_transform_non_log,
                                     data_transform=TRUE,
                                     group_info = group_info,
                                     paired=FALSE, test='t.test',
                                     adjust_p_method='BH',
                                     sig_stat='p.adj', sig_pvalue=0.05,
                                     sig_FC=2)$DE_species_table_sig
# conduct PCA
DEspec_PCA <- PCA(exp_transform_table,
                  group_info = group_info,
                  sig_feature = DE_species_table_sig$feature,
                  scaling=TRUE, centering=TRUE, group_num=2,
                  cluster_method='kmeans', var1=NULL, var2=NULL,
                  insert_ref_group=NULL, ref_group=NULL,
                  n_PC=c(1,2), top_n_feature=10)

# view result: PCA prcomp
head(DEspec_PCA[[1]], 1)

# view result: data frame of PCA rotated data
head(DEspec_PCA[[2]][,1:4], 5)

# view result: data frame of PCA contribution table
head(DEspec_PCA[[3]][,1:3], 5)


###################################################
### code chunk number 36: LipidSigR.Rnw:633-635
###################################################
# view result: PCA plot
DEspec_PCA[[4]]


###################################################
### code chunk number 37: LipidSigR.Rnw:643-645
###################################################
# view result: scree plot
DEspec_PCA[[5]]


###################################################
### code chunk number 38: LipidSigR.Rnw:653-655
###################################################
# view result: bar plot of contribution of top 10 features
DEspec_PCA[[6]]


###################################################
### code chunk number 39: LipidSigR.Rnw:663-665
###################################################
# view result: correlation circle plot of variables
DEspec_PCA[[7]]


###################################################
### code chunk number 40: DE_lipid species: dimensionality reduction - t-SNE
###################################################
# data processing of exp_data
exp_transform_table <- data_process(exp_data, exclude_var_missing=TRUE,
                                    missing_pct_limit=50,
                                    replace_zero=TRUE, xmin=0.5,
                                    zero2what='min', replace_NA=TRUE,
                                    ymin=0.5, NA2what='min',
                                    pct_transform=TRUE,
                                    data_transform=TRUE,
                                    trans_type='log',
                                    centering=FALSE, scaling=FALSE)
# data processing of exp_data (without log10 transformation)
exp_transform_non_log <- data_process(exp_data,
                                      exclude_var_missing=TRUE,
                                      missing_pct_limit=50,
                                      replace_zero=TRUE,
                                      zero2what='min', xmin=0.5,
                                      replace_NA=TRUE, NA2what='min',
                                      ymin=0.5, pct_transform=TRUE,
                                      data_transform=FALSE,
                                      trans_type='log',
                                      centering=FALSE, scaling=FALSE)
# filter significant lipid
DE_species_table_sig <- DE_species_2(exp_transform_non_log,
                                     data_transform=TRUE,
                                     group_info = group_info,
                                     paired=FALSE, test='t.test',
                                     adjust_p_method='BH',
                                     sig_stat='p.adj', sig_pvalue=0.05,
                                     sig_FC=2)$DE_species_table_sig
# conduct t-SNE
DEspec_tsne <- tsne(exp_transform_table,
                    group_info = group_info,
                    sig_feature=DE_species_table_sig$feature,
                    pca=TRUE, perplexity=5, max_iter=500,
                    cluster_method='kmeans', group_num=2,
                    var1 = 'euclidean', var2 = NULL,
                    insert_ref_group=NULL, ref_group=NULL)

# view result: data frame of t-SNE data
head(DEspec_tsne[[1]], 5)


###################################################
### code chunk number 41: LipidSigR.Rnw:720-722
###################################################
# view result: t-SNE plot
DEspec_tsne[[2]]


###################################################
### code chunk number 42: DE_lipid species: dimensionality reduction - UMAP
###################################################
# data processing of exp_data
exp_transform_table <- data_process(exp_data, exclude_var_missing=TRUE,
                                    missing_pct_limit=50,
                                    replace_zero=TRUE, xmin=0.5,
                                    zero2what='min', replace_NA=TRUE,
                                    ymin=0.5, NA2what='min',
                                    pct_transform=TRUE,
                                    data_transform=TRUE,
                                    trans_type='log',
                                    centering=FALSE, scaling=FALSE)
# data processing of exp_data (without log10 transformation)
exp_transform_non_log <- data_process(exp_data,
                                      exclude_var_missing=TRUE,
                                      missing_pct_limit=50,
                                      replace_zero=TRUE,
                                      zero2what='min', xmin=0.5,
                                      replace_NA=TRUE, NA2what='min',
                                      ymin=0.5, pct_transform=TRUE,
                                      data_transform=FALSE,
                                      trans_type='log',
                                      centering=FALSE, scaling=FALSE)
# filter significant lipid
DE_species_table_sig <- DE_species_2(exp_transform_non_log,
                                     data_transform=TRUE,
                                     group_info = group_info,
                                     paired=FALSE, test='t.test',
                                     adjust_p_method='BH',
                                     sig_stat='p.adj', sig_pvalue=0.05,
                                     sig_FC=2)$DE_species_table_sig
# conduct UMAP
DEspec_UMAP <- UMAP(exp_transform_table, group_num=2,
                    group_info=group_info, n_neighbors=15,
                    sig_feature=DE_species_table_sig$feature,
                    metric='euclidean', scale=TRUE, var1=NULL,
                    cluster_method='kmeans', var2=NULL,
                    insert_ref_group=NULL, ref_group=NULL)

# view result: data frame of UMAP data
head(DEspec_UMAP[[1]], 5)


###################################################
### code chunk number 43: LipidSigR.Rnw:777-779
###################################################
# view result: UMAP plot
DEspec_UMAP[[2]]


###################################################
### code chunk number 44: DE_lipid species: dimensionality reduction - PLS-DA
###################################################
# data processing of exp_data
exp_transform_table <- data_process(exp_data, replace_zero=TRUE,
                                    exclude_var_missing=TRUE,
                                    missing_pct_limit=50,
                                    zero2what='min',
                                    xmin=0.5, replace_NA=TRUE,
                                    NA2what='min', ymin=0.5,
                                    pct_transform=TRUE,
                                    data_transform=TRUE,
                                    trans_type='log', centering=FALSE,
                                    scaling=FALSE)
# data processing of exp_data (without log10 transformation)
exp_transform_non_log <- data_process(exp_data,
                                      exclude_var_missing=TRUE,
                                      missing_pct_limit=50,
                                      replace_zero=TRUE,
                                      zero2what='min', xmin=0.5,
                                      replace_NA=TRUE, NA2what='min',
                                      ymin=0.5, pct_transform=TRUE,
                                      data_transform=FALSE,
                                      trans_type='log',
                                      centering=FALSE, scaling=FALSE)
# filter significant lipid
DE_species_table_sig <- DE_species_2(exp_transform_non_log,
                                     data_transform=TRUE,
                                     group_info = group_info,
                                     paired=FALSE, test='t.test',
                                     adjust_p_method='BH',
                                     sig_stat='p.adj', sig_pvalue=0.05,
                                     sig_FC=2)$DE_species_table_sig
# conduct PLSDA
DEspec_PLSDA <- PLSDA(exp_transform_table, scaling=TRUE,
                      group_info=group_info, ncomp=2,
                      sig_feature=DE_species_table_sig$feature,
                      cluster_method='group_info',
                      group_num=NULL, var1=NULL, var2=NULL,
                      insert_ref_group=NULL, ref_group=NULL)

# view result: data frame of sample variate
head(DEspec_PLSDA[[1]], 5)

# view result: data frame of sample loading
head(DEspec_PLSDA[[2]], 5)


###################################################
### code chunk number 45: LipidSigR.Rnw:835-837
###################################################
# view result: PLS-DA sample plot
DEspec_PLSDA[[3]]


###################################################
### code chunk number 46: LipidSigR.Rnw:847-849
###################################################
# view result: PLS-DA variable plot
DEspec_PLSDA[[4]]


###################################################
### code chunk number 47: DE_lipid species: hierarchical clustering
###################################################
# data processing of exp_data
exp_transform <- data_process(exp_data, exclude_var_missing=TRUE,
                              missing_pct_limit=50, replace_zero=TRUE,
                              zero2what='min', xmin=0.5, replace_NA=TRUE,
                              NA2what='min', ymin=0.5, pct_transform=TRUE,
                              data_transform=TRUE, trans_type='log',
                              centering=FALSE, scaling=FALSE)
# data processing of exp_data (without log10 transformation)
exp_transform_non_log <- data_process(exp_data,
                                      exclude_var_missing=TRUE,
                                      missing_pct_limit=50,
                                      replace_zero=TRUE, zero2what='min',
                                      xmin=0.5, replace_NA=TRUE,
                                      NA2what='min', ymin=0.5,
                                      pct_transform=TRUE,
                                      data_transform=FALSE,
                                      trans_type='log',
                                      centering=FALSE, scaling=FALSE)
# filter lipid_char_table according to exp_transform
lipid_char_filter <- lipid_char_table %>%
  filter(feature %in% exp_transform$feature)
# filter significant lipid
DE_species_table_sig <- DE_species_2(exp_transform_non_log,
                                     data_transform=TRUE,
                                     group_info = group_info,
                                     paired=FALSE, test='t.test',
                                     adjust_p_method='BH',
                                     sig_stat='p.adj', sig_pvalue=0.05,
                                     sig_FC=2)$DE_species_table_sig
# get lipid characteristics
char_var <- colnames(lipid_char_filter)[-1]
# conduct hierarchical clustering
DEspec_Hcluster<- Hclustering(exp_transform, DE_species_table_sig,
                              group_info, lipid_char_filter,
                              char_var = char_var[1],
                              distfun='pearson', hclustfun='complete')

# view result: matrix of all lipid species heatmap
head(DEspec_Hcluster$all.lipid.data[, 1:4], 5)

# view result: matrix of significant lipid species heatmap
head(DEspec_Hcluster$sig.lipid.data[, 1:4], 5)


###################################################
### code chunk number 48: LipidSigR.Rnw:909-911
###################################################
# view result: heatmap of all lipid species
DEspec_Hcluster$all.lipid


###################################################
### code chunk number 49: LipidSigR.Rnw:919-920
###################################################
# view result: heatmap of significant lipid species


###################################################
### code chunk number 50: DE_lipid species: characteristics analysis
###################################################
# data processing of exp_data
exp_transform <- data_process(exp_data, exclude_var_missing=TRUE,
                              missing_pct_limit=50,
                              replace_zero=TRUE, zero2what='min',
                              xmin=0.5, replace_NA=TRUE,
                              NA2what='min', pct_transform=TRUE,
                              data_transform=TRUE, trans_type='log',
                              ymin=0.5, centering=FALSE, scaling=FALSE)
# data processing of exp_data (without log10 transformation)
exp_transform_non_log <- data_process(exp_data, exclude_var_missing=TRUE,
                                      missing_pct_limit=50,
                                      replace_zero=TRUE,
                                      zero2what='min', xmin=0.5,
                                      replace_NA=TRUE, NA2what='min',
                                      ymin=0.5, pct_transform=TRUE,
                                      data_transform=FALSE,
                                      trans_type='log',
                                      centering=FALSE, scaling=FALSE)
# filter lipid_char_table according to exp_transform
lipid_char_filter <- lipid_char_table %>%
  filter(feature %in% exp_transform_non_log$feature)
# filter significant lipid
DE_species_table_sig <- DE_species_2(exp_transform_non_log,
                                     data_transform=TRUE,
                                     group_info = group_info,
                                     paired=FALSE, test='t.test',
                                     adjust_p_method='BH',
                                     sig_stat='p.adj', sig_pvalue=0.05,
                                     sig_FC=2)$DE_species_table_sig
# get lipid characteristics
char_var <- colnames(lipid_char_filter)[-1]
# conduct Sig_lipid_feature
sig_feature_result <- Sig_lipid_feature(DE_species_table_sig,
                                        lipid_char_filter,
                                        char_var[1], sig_FC=2)


###################################################
### code chunk number 51: LipidSigR.Rnw:971-973
###################################################
# view result: bar chart
sig_feature_result$barPlot


###################################################
### code chunk number 52: LipidSigR.Rnw:983-985
###################################################
# view result: lollipop plot
sig_feature_result$lolipop


###################################################
### code chunk number 53: LipidSigR.Rnw:994-996
###################################################
# view result: word cloud
sig_feature_result$word


###################################################
### code chunk number 54: DE_lipid species: enrichment
###################################################
# data processing of exp_data
exp_transform <- data_process(exp_data, exclude_var_missing=TRUE,
                              missing_pct_limit=50, replace_zero=TRUE,
                              zero2what='min', xmin=0.5, replace_NA=TRUE,
                              NA2what='min', ymin=0.5, pct_transform=TRUE,
                              data_transform=TRUE, trans_type='log',
                              centering=FALSE, scaling=FALSE)
# data processing of exp_data (without log10 transformation)
exp_transform_non_log <- data_process(exp_data,
                                      exclude_var_missing=TRUE,
                                      missing_pct_limit=50,
                                      replace_zero=TRUE, zero2what='min',
                                      xmin=0.5, replace_NA=TRUE,
                                      NA2what='min', ymin=0.5,
                                      pct_transform=TRUE, data_transform=FALSE,
                                      trans_type='log',
                                      centering=FALSE, scaling=FALSE)
# filter lipid_char_table according to exp_transform
lipid_char_filter <- lipid_char_table %>%
  filter(feature %in% exp_transform$feature)
# filter significant lipid
DE_species_table_sig <- DE_species_2(exp_transform_non_log,
                                     data_transform = TRUE,
                                     group_info = group_info,
                                     paired=FALSE, test='t.test',
                                     adjust_p_method='BH',
                                     sig_stat='p.adj', sig_pvalue=0.05,
                                     sig_FC=2)$DE_species_table_sig
# get lipid characteristics
char_var <- colnames(lipid_char_filter)[-1]
# conduct enrichment analysis
enrich_result <- Enrichment(DE_species_table_sig,
                            lipid_char_table = lipid_char_filter,
                            char_var=char_var[1], sig_pvalue=0.05)
# view result: summary data frame of enrichment
head(enrich_result$enrich_char_table[, 1:5], 5)


###################################################
### code chunk number 55: LipidSigR.Rnw:1051-1053
###################################################
# view result: enrichment plot
enrich_result$enrich_char_barplot


###################################################
### code chunk number 56: DE_lipid species: enrichment_pathview_data
###################################################
# data frame of pathway information
data("DE_lipid_gene_path")
DE.lipid.gene.path <- DE_lipid_gene_path
head(DE.lipid.gene.path[, 1:4], 5)

# list of genes included in each pathway
data("DE_pathway_gene_list")
DE.pathway.gene.list <- DE_pathway_gene_list
head(DE.pathway.gene.list, 2)


###################################################
### code chunk number 57: DE_lipid species: enrichment_pathview (eval = FALSE)
###################################################
## # data processing of exp_data
## exp_transform <- data_process(exp_data, exclude_var_missing=TRUE,
##                               missing_pct_limit=50, replace_zero=TRUE,
##                               zero2what='min', xmin=0.5, replace_NA=TRUE,
##                               ymin=0.5, NA2what='min', data_transform=TRUE,
##                               pct_transform=TRUE, trans_type='log',
##                               centering=FALSE, scaling=FALSE)
## # data processing of exp_data (without log10 transformation)
## exp_transform_non_log <- data_process(exp_data, replace_zero=TRUE,
##                                       exclude_var_missing=TRUE,
##                                       missing_pct_limit=50,
##                                       zero2what='min', xmin=0.5,
##                                       replace_NA=TRUE, ymin=0.5,
##                                       NA2what='min', centering=FALSE,
##                                       pct_transform=TRUE, scaling=FALSE,
##                                       data_transform=FALSE, trans_type='log')
## # filter lipid_char_table according to exp_transform
## lipid_char_filter <- lipid_char_table %>%
##   filter(feature %in% exp_transform$feature)
## # filter significant lipid
## DE_species_table_sig <- DE_species_2(exp_transform_non_log,
##                                      data_transform=TRUE,
##                                      group_info = group_info,
##                                      paired=FALSE, test = 't.test',
##                                      adjust_p_method = 'BH',
##                                      sig_stat = 'p.adj', sig_pvalue=0.05,
##                                      sig_FC=2)$DE_species_table_sig
## # get lipid characteristics
## char_var <- colnames(lipid_char_filter)[-1]
## # get enrich characteristic table
## enrich_char_table <- Enrichment(DE_species_table_sig,
##                                 lipid_char_table=lipid_char_filter,
##                                 char_var = char_var[1],
##                                 sig_pvalue=0.05)$enrich_char_table
## # filter "class" of significant characteristics
## sig_enrich_class <- enrich_char_table %>% filter(significant == 'YES') %>%
##   distinct(characteristic) %>% .$characteristic
## # set output path
## dir.create(file.path(getwd(),"pathview_result"), recursive=TRUE,
##            showWarnings=TRUE)
## outPath <- file.path(getwd(), "pathview_result")
## # conduct pathview analysis
## path_result <- pathview_function(sig_enrich_class, path=outPath,
##                                  DE.lipid.gene.path, DE.pathway.gene.list)


###################################################
### code chunk number 58: LipidSigR.Rnw:1130-1172
###################################################
# data processing of exp_data
exp_transform <- data_process(exp_data, exclude_var_missing=TRUE,
                              missing_pct_limit=50, replace_zero=TRUE,
                              zero2what='min', xmin=0.5, replace_NA=TRUE,
                              ymin=0.5, NA2what='min', data_transform=TRUE,
                              pct_transform=TRUE, trans_type='log',
                              centering=FALSE, scaling=FALSE)
# data processing of exp_data (without log10 transformation)
exp_transform_non_log <- data_process(exp_data, replace_zero=TRUE,
                                      exclude_var_missing=TRUE,
                                      missing_pct_limit=50,
                                      zero2what='min', xmin=0.5,
                                      replace_NA=TRUE, ymin=0.5,
                                      NA2what='min', centering=FALSE,
                                      pct_transform=TRUE, scaling=FALSE,
                                      data_transform=FALSE, trans_type='log')
# filter lipid_char_table according to exp_transform
lipid_char_filter <- lipid_char_table %>%
  filter(feature %in% exp_transform$feature)
# filter significant lipid
DE_species_table_sig <- DE_species_2(exp_transform_non_log,
                                     data_transform=TRUE,
                                     group_info = group_info,
                                     paired=FALSE, test = 't.test',
                                     adjust_p_method = 'BH',
                                     sig_stat = 'p.adj', sig_pvalue=0.05,
                                     sig_FC=2)$DE_species_table_sig
# get lipid characteristics
char_var <- colnames(lipid_char_filter)[-1]
# get enrich characteristic table
enrich_char_table <- Enrichment(DE_species_table_sig,
                                lipid_char_table=lipid_char_filter,
                                char_var = char_var[1],
                                sig_pvalue=0.05)$enrich_char_table
# filter "class" of significant characteristics
sig_enrich_class <- enrich_char_table %>% filter(significant == 'YES') %>%
  distinct(characteristic) %>% .$characteristic
# set output path
outPath <- file.path(getwd())
# conduct pathview analysis
path_result <- pathview_function(sig_enrich_class, path=outPath,
                                 DE.lipid.gene.path, DE.pathway.gene.list)


###################################################
### code chunk number 59: LipidSigR.Rnw:1175-1179
###################################################
# view result: data frame of pathway information
path_result
# view result: pathview plot
showimage::show_image(file.path(outPath,"hsa00565.pathview.png"))


###################################################
### code chunk number 60: DE_lipid characteristics: differentially expressed analysis
###################################################
# get lipid characteristics
char_var <- colnames(lipid_char_table)[-1]
# aggregated(sum) expression data by selected characteristics
exp_data_Spe2Char <- Species2Char(exp_data, lipid_char_table,
                                  char_var = char_var[4])
# data processing of exp_data (without log10 transformation)
exp_transform_non_log <- data_process(exp_data_Spe2Char,
                                      exclude_var_missing=TRUE,
                                      missing_pct_limit=50,
                                      replace_zero=TRUE,
                                      zero2what='min', xmin=0.5,
                                      replace_NA=TRUE, NA2what='min',
                                      ymin=0.5, pct_transform=TRUE,
                                      data_transform=FALSE,
                                      trans_type='log',
                                      centering=FALSE, scaling=FALSE)
# conduct deferentially expressed of lipid characters
DE_char_result <- DE_char_2(exp_transform_non_log, data_transform=TRUE,
                            group_info = group_info, paired=FALSE,
                            sig_pvalue=0.05, sig_FC=2,
                            insert_ref_group=NULL, ref_group=NULL)

# view result: data frame of expression
head(DE_char_result$DE_char_exp_data[, 1:5], 5)

# view result: data frame of statistical analysis
head(DE_char_result$DE_char_table_all[, 1:5], 5)

# view result: data frame with value of the continuous lipid characteristics
head(DE_char_result$DE_char_combined_table[, 1:5])

# view result: data frame with statistics of t.test
head(DE_char_result$DE_char_combine_result_table[, 1:5])


###################################################
### code chunk number 61: LipidSigR.Rnw:1232-1234
###################################################
# view result: bar plot of split_class
DE_char_result$DE_char_barplot


###################################################
### code chunk number 62: LipidSigR.Rnw:1244-1246
###################################################
# view result: sqrt-scaled bar plot of split_class
DE_char_result$DE_char_barplot_sqrt


###################################################
### code chunk number 63: LipidSigR.Rnw:1256-1258
###################################################
# view result: line plot of split_class
DE_char_result$DE_char_trendplot


###################################################
### code chunk number 64: LipidSigR.Rnw:1266-1268
###################################################
# view result: sqrt-scaled line plot of split_class
DE_char_result$DE_char_trendplot_sqrt


###################################################
### code chunk number 65: LipidSigR.Rnw:1277-1279
###################################################
# view result: box plot of split_class
DE_char_result$DE_char_boxplot


###################################################
### code chunk number 66: DE_lipid characteristics: differentially expressed analysis_Subgroup
###################################################
# get lipid characteristics
char_var <- colnames(lipid_char_table)[-1]
# subgroup deferentially expressed of lipid characters
DE.sub.char.2 <- DE_sub_char_2(exp_data, data_transform=TRUE,
                               lipid_char_table=lipid_char_table,
                               split_var = char_var[1],
                               char_var = char_var[4],
                               group_info = group_info,
                               paired=FALSE, sig_pvalue=0.05,
                               sig_FC=2, exclude_var_missing=TRUE,
                               missing_pct_limit=50,
                               replace_zero=TRUE, zero2what='min',
                               xmin=0.5, replace_NA=TRUE,
                               NA2what='min', ymin=0.5,
                               pct_transform=TRUE, trans_type='log',
                               centering=FALSE, scaling=FALSE)
# get class of characteristics
char.class <- unique(DE.sub.char.2[[2]][1])
# visualize subgroup deferentially expressed of lipid characters
sub_char_result <- DE_sub_char_plot_2(DE.sub.char.2[[2]],
                                      DE.sub.char.2[[3]],
                                      group_info=group_info,
                                      char_var=char_var[4],
                                      split_var=char_var[1],
                                      split_class=char.class[5,],
                                      insert_ref_group=NULL, ref_group=NULL)


###################################################
### code chunk number 67: LipidSigR.Rnw:1322-1324
###################################################
# view result: bar plot of split_class
sub_char_result[[1]]


###################################################
### code chunk number 68: LipidSigR.Rnw:1333-1335
###################################################
# view result: sqrt-scaled bar plot of split_class
sub_char_result[[4]]


###################################################
### code chunk number 69: LipidSigR.Rnw:1343-1345
###################################################
# view result: line plot of split_class
sub_char_result[[2]]


###################################################
### code chunk number 70: LipidSigR.Rnw:1355-1357
###################################################
# view result: sqrt-scaled line plot of split_class
sub_char_result[[5]]


###################################################
### code chunk number 71: LipidSigR.Rnw:1365-1367
###################################################
# view result: box plot of split_class
sub_char_result[[3]]


###################################################
### code chunk number 72: DE_lipid characteristics: dimensionality reduction - PCA
###################################################
# get lipid characteristics
char_var <- colnames(lipid_char_table)[-1]
# sum expression data by selected characteristics
exp_data_Spe2Char <- Species2Char(exp_data, lipid_char_table,
                                  char_var = char_var[4])
# data processing of exp_data_Spe2Char
exp_transform_class <- data_process(exp_data_Spe2Char,
                                    exclude_var_missing=TRUE,
                                    missing_pct_limit=50,
                                    replace_zero=TRUE,
                                    zero2what='NA', xmin=0.5,
                                    replace_NA=TRUE,
                                    NA2what='min', ymin=0.5,
                                    pct_transform=TRUE,
                                    data_transform=TRUE,
                                    trans_type='log',
                                    centering=FALSE, scaling=FALSE)
# data processing of exp_data (without log10 transformation)
exp_transform_non_log <- data_process(exp_data_Spe2Char,
                                      exclude_var_missing=TRUE,
                                      missing_pct_limit=50,
                                      replace_zero=TRUE,
                                      zero2what='min', xmin=0.5,
                                      replace_NA=TRUE, NA2what='min',
                                      ymin=0.5, pct_transform=TRUE,
                                      data_transform=FALSE,
                                      trans_type='log',
                                      centering=FALSE, scaling=FALSE)
# filter significant lipid characteristics
DE_char_table_sig <- DE_char_2(exp_transform_non_log,
                               data_transform=TRUE,
                               group_info = group_info,
                               paired=FALSE, sig_pvalue=0.05,
                               sig_FC=2)$DE_char_table_all %>%
                               filter(sig=="yes")
# conduct PCA
DEchar_PCA <- PCA(exp_transform_class, group_info = group_info,
                  sig_feature = DE_char_table_sig[,1],
                  scaling=TRUE, centering=TRUE,
                  cluster_method='kmeans',
                  group_num=2, var1 = NULL, var2 = NULL,
                  insert_ref_group=NULL, ref_group=NULL,
                  n_PC=c(1,2), top_n_feature=10)

# view result: PCA prcomp
head(DEchar_PCA[[1]], 2)

# view result: data frame of PCA rotated data
head(DEchar_PCA[[2]][,1:5], 5)

# view result: data frame of PCA contribution table
head(DEchar_PCA[[3]])


###################################################
### code chunk number 73: LipidSigR.Rnw:1438-1440
###################################################
# view result: PCA plot
DEchar_PCA[[4]]


###################################################
### code chunk number 74: LipidSigR.Rnw:1448-1450
###################################################
# view result: scree plot
DEchar_PCA[[5]]


###################################################
### code chunk number 75: LipidSigR.Rnw:1458-1460
###################################################
# view result: correlation circle plot of variables
DEchar_PCA[[7]]


###################################################
### code chunk number 76: LipidSigR.Rnw:1470-1472
###################################################
# view result: bar plot of features contribution
DEchar_PCA[[6]]


###################################################
### code chunk number 77: DE_lipid characteristics: dimensionality reduction - t-SNE
###################################################
# get lipid characteristics
char_var <- colnames(lipid_char_table)[-1]
# sum expression data by selected characteristics
exp_data_Spe2Char <- Species2Char(exp_data, lipid_char_table,
                                  char_var = char_var[4])
# data processing of exp_data_Spe2Char
exp_transform_class <- data_process(exp_data_Spe2Char,
                                    exclude_var_missing=TRUE,
                                    missing_pct_limit=50,
                                    replace_zero=TRUE,
                                    zero2what='NA', xmin=0.5,
                                    replace_NA=TRUE,
                                    NA2what='min', ymin=0.5,
                                    pct_transform=TRUE,
                                    data_transform=TRUE,
                                    trans_type='log',
                                    centering=FALSE, scaling=FALSE)
# data processing of exp_data (without log10 transformation)
exp_transform_non_log <- data_process(exp_data_Spe2Char,
                                      exclude_var_missing=TRUE,
                                      missing_pct_limit=50,
                                      replace_zero=TRUE,
                                      zero2what='min', xmin=0.5,
                                      replace_NA=TRUE, NA2what='min',
                                      ymin=0.5, pct_transform=TRUE,
                                      data_transform=FALSE,
                                      trans_type='log',
                                      centering=FALSE, scaling=FALSE)
# filter significant lipid characteristics
DE_char_table_sig <- DE_char_2(exp_transform_non_log,
                               data_transform=TRUE,
                               group_info = group_info,
                               paired=FALSE, sig_pvalue=0.05,
                               sig_FC=2)$DE_char_table_all %>%
                               filter(sig=="yes")
# conduct t-SNE
DEchar_tsne <- tsne(exp_transform_class, group_info = group_info,
                    sig_feature=DE_char_table_sig[,1],
                    pca=TRUE, perplexity=5, max_iter=500,
                    cluster_method='kmeans', group_num=2,
                    var1 = 'euclidean', var2=NULL,
                    insert_ref_group=NULL, ref_group=NULL)

# view result: data frame of t-SNE data
head(DEchar_tsne[[1]], 5)


###################################################
### code chunk number 78: LipidSigR.Rnw:1531-1533
###################################################
# view result: t-SNE plot
DEchar_tsne[[2]]


###################################################
### code chunk number 79: DE_lipid characteristics: dimensionality reduction - UMAP
###################################################
# get lipid characteristics
char_var <- colnames(lipid_char_table)[-1]
# sum expression data by selected characteristics
exp_data_Spe2Char <- Species2Char(exp_data, lipid_char_table,
                                  char_var = char_var[4])
# data processing of exp_data_Spe2Char
exp_transform_class <- data_process(exp_data_Spe2Char,
                                    exclude_var_missing=TRUE,
                                    missing_pct_limit=50,
                                    replace_zero=TRUE,
                                    zero2what='NA', xmin=0.5,
                                    replace_NA=TRUE,
                                    NA2what='min', ymin=0.5,
                                    pct_transform=TRUE,
                                    data_transform=TRUE,
                                    trans_type='log',
                                    centering=FALSE, scaling=FALSE)
# data processing of exp_data (without log10 transformation)
exp_transform_non_log <- data_process(exp_data_Spe2Char,
                                      exclude_var_missing=TRUE,
                                      missing_pct_limit=50,
                                      replace_zero=TRUE,
                                      zero2what='min', xmin=0.5,
                                      replace_NA=TRUE, NA2what='min',
                                      ymin=0.5, pct_transform=TRUE,
                                      data_transform=FALSE,
                                      trans_type='log',
                                      centering=FALSE, scaling=FALSE)
# filter significant lipid characteristics
DE_char_table_sig <- DE_char_2(exp_transform_non_log,
                               data_transform=TRUE,
                               group_info = group_info,
                               paired=FALSE, sig_pvalue=0.05,
                               sig_FC=2)$DE_char_table_all %>%
                               filter(sig=="yes")
# conduct UMAP
DEchar_UMAP <- UMAP(exp_transform_class, group_info = group_info,
                    sig_feature = DE_char_table_sig[,1],
                    n_neighbors=15, scale=TRUE, group_num=2,
                    metric = 'euclidean',
                    cluster_method = 'kmeans',
                    var1=NULL, var2=NULL,
                    insert_ref_group=NULL, ref_group=NULL)

# view result: data frame of UMAP data
head(DEchar_UMAP[[1]], 5)


###################################################
### code chunk number 80: LipidSigR.Rnw:1593-1595
###################################################
# view result: UMAP plot
DEchar_UMAP[[2]]


###################################################
### code chunk number 81: DE_lipid characteristics: dimensionality reduction - PLS-DA
###################################################
# get lipid characteristics
char_var <- colnames(lipid_char_table)[-1]
# sum expression data by selected characteristics
exp_data_Spe2Char <- Species2Char(exp_data, lipid_char_table,
                                  char_var = char_var[4])
# data processing of exp_data_Spe2Char
exp_transform_class <- data_process(exp_data_Spe2Char,
                                    exclude_var_missing=TRUE,
                                    missing_pct_limit=50,
                                    replace_zero=TRUE,
                                    zero2what='NA', xmin=0.5,
                                    replace_NA=TRUE,
                                    NA2what='min', ymin=0.5,
                                    pct_transform=TRUE,
                                    data_transform=TRUE,
                                    trans_type='log',
                                    centering=FALSE, scaling=FALSE)
# data processing of exp_data (without log10 transformation)
exp_transform_non_log <- data_process(exp_data_Spe2Char,
                                      exclude_var_missing=TRUE,
                                      missing_pct_limit=50,
                                      replace_zero=TRUE,
                                      zero2what='min', xmin=0.5,
                                      replace_NA=TRUE, NA2what='min',
                                      ymin=0.5, pct_transform=TRUE,
                                      data_transform=FALSE,
                                      trans_type='log',
                                      centering=FALSE, scaling=FALSE)
# filter significant lipid characteristics
DE_char_table_sig <- DE_char_2(exp_transform_non_log,
                               data_transform=TRUE,
                               group_info = group_info,
                               paired=FALSE, sig_pvalue=0.05,
                               sig_FC=2)$DE_char_table_all %>%
                               filter(sig=="yes")
# conduct PLSDA
DEchar_PLSDA <- PLSDA(exp_transform_class, group_info = group_info,
                      sig_feature=DE_char_table_sig[,1],
                      ncomp=2, scaling=TRUE, cluster_method='group_info',
                      group_num = NULL, var1=NULL, var2=NULL,
                      insert_ref_group=NULL, ref_group=NULL)

# view result: data frame of sample variate
head(DEchar_PLSDA[[1]], 5)

# view result: data frame of sample loading
head(DEchar_PLSDA[[2]])


###################################################
### code chunk number 82: LipidSigR.Rnw:1654-1656
###################################################
# view result: PLS-DA sample plot
DEchar_PLSDA[[3]]


###################################################
### code chunk number 83: LipidSigR.Rnw:1664-1666
###################################################
# view result: PLS-DA loading plot
DEchar_PLSDA[[4]]


###################################################
### code chunk number 84: DE_lipid characteristics: hierarchical clustering
###################################################
# get lipid characteristics
char_var <- colnames(lipid_char_table)[-1]
# sum expression data by selected characteristics
exp_data_Spe2Char <- Species2Char(exp_data, lipid_char_table,
                                  char_var = char_var[4])
# data processing of exp_data_Spe2Char
exp_transform_class <- data_process(exp_data_Spe2Char,
                                    exclude_var_missing=TRUE,
                                    missing_pct_limit=50, replace_zero=TRUE,
                                    replace_NA=TRUE,zero2what='NA',
                                    NA2what='min', pct_transform=TRUE,
                                    xmin=0.5, data_transform=TRUE,
                                    ymin=0.5, trans_type='log',
                                    centering=FALSE, scaling=FALSE)
# data processing of exp_data (without log10 transformation)
exp_transform_non_log <- data_process(exp_data_Spe2Char,
                                      exclude_var_missing=TRUE,
                                      missing_pct_limit=50,
                                      replace_zero=TRUE, zero2what='min',
                                      xmin=0.5, replace_NA=TRUE, NA2what='min',
                                      ymin=0.5, pct_transform=TRUE,
                                      data_transform=FALSE, trans_type='log',
                                      centering=FALSE, scaling=FALSE)
# filter significant lipid characteristics
DE_char_table_sig <- DE_char_2(exp_transform_non_log,
                               data_transform=TRUE, group_info = group_info,
                               paired=FALSE, sig_pvalue=0.05, sig_FC=2
                              )$DE_char_table_all %>% filter(sig=="yes")
# conduct hierarchical clustering
DEchar_Hcluster <- Hclustering(exp_transform_class, DE_char_table_sig,
                               group_info, lipid_char_table=NULL,
                               char_var=NULL, distfun='pearson',
                               hclustfun='complete')

# view result: matrix of all lipid species heatmap
head(DEchar_Hcluster$all.lipid.data[, 1:5], 5)

# view result: matrix of significant lipid species heatmap
head(DEchar_Hcluster$sig.lipid.data[, 1:5])


###################################################
### code chunk number 85: LipidSigR.Rnw:1721-1723
###################################################
# view result: heatmap of all lipid species
DEchar_Hcluster$all.lipid


###################################################
### code chunk number 86: LipidSigR.Rnw:1731-1733
###################################################
# view result: heatmap of significant lipid species
DEchar_Hcluster$sig.lipid


###################################################
### code chunk number 87: load_ML_data
###################################################
# clears all objects from workspace
rm(list = ls())

# lipid expression data
data("ML_exp_data")
exp_data <- ML_exp_data
head(exp_data[, 1:5], 5)

# lipid characteristics table
data("ML_lipid_char_table")
lipid_char_table <- ML_lipid_char_table
head(lipid_char_table[, 1:4], 5)

# condition table
data("ML_condition_table")
condition_table <- ML_condition_table
head(condition_table, 5)


###################################################
### code chunk number 88: ML_data_process
###################################################
# lipid expression data
head(exp_data[, 1:5], 5)
# data processing of exp_data
exp_transform_table <- data_process(exp_data, exclude_var_missing=TRUE,
                                    missing_pct_limit=50, replace_zero=TRUE,
                                    zero2what='min', xmin=0.5, replace_NA=TRUE,
                                    NA2what='min', ymin=0.5,
                                    pct_transform=TRUE,
                                    data_transform=TRUE, trans_type='log',
                                    centering=FALSE, scaling=FALSE)
# exp_data after data processing
head(exp_transform_table[, 1:5], 5)


###################################################
### code chunk number 89: ML: ROC curve
###################################################
# get lipid characteristics
char_var <- colnames(lipid_char_table)[-1]
# data processing of machine learning
ML_data <- ML_data_process(exp_data, group_info = condition_table,
                           lipid_char_table, char_var[1],
                           exclude_var_missing=TRUE, missing_pct_limit=50,
                           replace_zero=TRUE, zero2what='min', xmin=0.5,
                           replace_NA=TRUE, NA2what='min', ymin=0.5,
                           pct_transform=TRUE, data_transform=TRUE,
                           trans_type='log', centering=FALSE, scaling=FALSE)
# conduct machine learning
ML_output <- ML_final(ML_data[[2]],ranking_method='Random_forest',
                      ML_method='Random_forest', split_prop=0.3, nfold=10)
# ROC curves
ROC_result <- ROC_plot_all(ML_output[[3]], ML_output[[5]], feature_n=10)

# view result: ROC data frame of 10 features
head(ROC_result[[3]][, 1:5], 5)

# view result: data frame of ROC values
head(ROC_result[[1]][, 1:5], 5)


###################################################
### code chunk number 90: LipidSigR.Rnw:1833-1835
###################################################
# view result: ROC curve plot
ROC_result[[2]]


###################################################
### code chunk number 91: LipidSigR.Rnw:1845-1847
###################################################
# view result: average ROC curve plot of 10 features
ROC_result[[4]]


###################################################
### code chunk number 92: ML: PR curve
###################################################
# get lipid characteristics
char_var <- colnames(lipid_char_table)[-1]
## data processing of machine learning
ML_data <- ML_data_process(exp_data, group_info = condition_table,
                           lipid_char_table, char_var[1],
                           exclude_var_missing=TRUE, missing_pct_limit=50,
                           replace_zero=TRUE, zero2what='min', xmin=0.5,
                           replace_NA=TRUE, NA2what='min', ymin=0.5,
                           pct_transform=TRUE, data_transform=TRUE,
                           trans_type='log', centering=FALSE, scaling=FALSE)
# conduct machine learning
ML_output <- ML_final(ML_data[[2]],ranking_method='Random_forest',
                      ML_method='Random_forest', split_prop=0.3, nfold=10)
# PR curves
PR_result <- PR_plot_all(ML_output[[4]], ML_output[[5]], feature_n=10)

# view result: data frame of precision and recall values
head(PR_result[[1]][, 1:5], 5)

# view result: data frame of PR values
head(PR_result[[3]][, 1:5], 5)


###################################################
### code chunk number 93: LipidSigR.Rnw:1879-1881
###################################################
# view result: PR curve plot
PR_result[[2]]


###################################################
### code chunk number 94: LipidSigR.Rnw:1891-1893
###################################################
# view result: average PR curve plot of 10 features
PR_result[[4]]


###################################################
### code chunk number 95: ML: model performance
###################################################
# get lipid characteristics
char_var <- colnames(lipid_char_table)[-1]
# data processing of machine learning
ML_data <- ML_data_process(exp_data, group_info = condition_table,
                           lipid_char_table, char_var[1],
                           exclude_var_missing=TRUE, missing_pct_limit=50,
                           replace_zero=TRUE, zero2what='min', xmin=0.5,
                           replace_NA=TRUE, NA2what='min', ymin=0.5,
                           pct_transform=TRUE, data_transform=TRUE,
                           trans_type='log', centering=FALSE, scaling=FALSE)
# conduct machine learning
ML_output <- ML_final(ML_data[[2]], ranking_method='Random_forest',
                      ML_method='Random_forest', split_prop=0.3, nfold=10)
# conduct model evaluation
evaluate_result <- evalution_plot(ML_output[[2]], method='Accuracy')

# view result: data frame of model evaluation information
head(evaluate_result[[1]][, 1:5], 5)


###################################################
### code chunk number 96: LipidSigR.Rnw:1939-1941
###################################################
# view result: model performance plot
evaluate_result[[2]]


###################################################
### code chunk number 97: ML: predicted probability
###################################################
# get lipid characteristics
char_var <- colnames(lipid_char_table)[-1]
# data processing of machine learning
ML_data <- ML_data_process(exp_data, group_info = condition_table,
                           lipid_char_table, char_var[1],
                           exclude_var_missing=TRUE, missing_pct_limit=50,
                           replace_zero=TRUE, zero2what='min', xmin=0.5,
                           replace_NA=TRUE, NA2what='min', ymin=0.5,
                           pct_transform=TRUE, data_transform=TRUE,
                           trans_type='log', centering=FALSE, scaling=FALSE)
# conduct machine learning
ML_output <- ML_final(ML_data[[2]], ranking_method='Random_forest',
                      ML_method='Random_forest', split_prop=0.3, nfold=10)
# compute and visualize the average predicted probabilities
prob_result <- probability_plot(ML_output[[1]], feature_n=10)

# view result: data frame of confusion matrix
head(prob_result[[1]][, 1:5], 5)

# view result: data frame of predicted probability and labels
head(prob_result[[4]][, 1:5], 5)


###################################################
### code chunk number 98: LipidSigR.Rnw:1976-1978
###################################################
# view result: the distribution of predicted probabilities
prob_result[[2]]


###################################################
### code chunk number 99: LipidSigR.Rnw:1988-1990
###################################################
# view result: confusion matrix of sample number and proportion
prob_result[[3]]


###################################################
### code chunk number 100: ML: feature importance_algorithm-based
###################################################
# get lipid characteristics
char_var <- colnames(lipid_char_table)[-1]
# data processing of machine learning
ML_data <- ML_data_process(exp_data, group_info = condition_table,
                           lipid_char_table, char_var[1],
                           exclude_var_missing=TRUE, missing_pct_limit=50,
                           replace_zero=TRUE, zero2what='min', xmin=0.5,
                           replace_NA=TRUE, NA2what='min', ymin=0.5,
                           pct_transform=TRUE, data_transform=TRUE,
                           trans_type='log', centering=FALSE, scaling=FALSE)
# conduct machine learning
ML_output <- ML_final(ML_data[[2]], ranking_method='Random_forest',
                      ML_method='Random_forest', split_prop=0.3, nfold=10)
# compute and rank the contribution of each feature
feature_result <- feature_plot(ML_output[[6]], ML_output[[7]],
                               feature_n=10, nfold=10)

# view result: data frame of the selected frequency
head(feature_result[[1]][, 1:5], 5)

# view result: data frame of feature importance
head(feature_result[[3]][, 1:5], 5)


###################################################
### code chunk number 101: LipidSigR.Rnw:2030-2032
###################################################
# view result: selected frequency plot
feature_result[[2]]


###################################################
### code chunk number 102: LipidSigR.Rnw:2040-2042
###################################################
# view result: feature importance plot
feature_result[[4]]


###################################################
### code chunk number 103: ML: SHAP analysis
###################################################
# get lipid characteristics
char_var <- colnames(lipid_char_table)[-1]
# data processing of machine learning
ML_data <- ML_data_process(exp_data, group_info = condition_table,
                           lipid_char_table, char_var[1],
                           exclude_var_missing=TRUE, missing_pct_limit=50,
                           replace_zero=TRUE, zero2what='min', xmin=0.5,
                           replace_NA=TRUE, NA2what='min', ymin=0.5,
                           pct_transform=TRUE, data_transform=TRUE,
                           trans_type='log', centering=FALSE, scaling=FALSE)
# conduct machine learning
ML_output <- ML_final(ML_data[[2]], ranking_method='Random_forest',
                      ML_method='Random_forest', split_prop=0.3, nfold=10)
# conduct SHAP
SHAP_output <- SHAP(ML_data[[2]], best_model=ML_output[[8]],
                    best_model_feature=ML_output[[9]],
                    ML_method='Random_forest', feature_n=10, nsim=5)


###################################################
### code chunk number 104: LipidSigR.Rnw:2077-2079
###################################################
# view result: SHAP feature importance plot
SHAP_output[[3]]


###################################################
### code chunk number 105: LipidSigR.Rnw:2087-2089
###################################################
# view result: SHAP summary plot
SHAP_output[[4]]


###################################################
### code chunk number 106: ML: SHAP analysis_sample
###################################################
# get lipid characteristics
char_var <- colnames(lipid_char_table)[-1]
# data processing of machine learning
ML_data <- ML_data_process(exp_data, group_info = condition_table,
                           lipid_char_table, char_var[1],
                           exclude_var_missing=TRUE, pct_transform=TRUE,
                           missing_pct_limit=50, replace_zero=TRUE,
                           zero2what='min', NA2what='min',
                           xmin=0.5, replace_NA=TRUE, ymin=0.5,
                           data_transform=TRUE, centering=FALSE,
                           trans_type='log', scaling=FALSE)
# conduct machine learning
ML_output <- ML_final(ML_data[[2]], ranking_method='Random_forest',
                      ML_method='Random_forest', split_prop=0.3, nfold=10)
# conduct SHAP
SHAP_output <- SHAP(ML_data[[2]], best_model=ML_output[[8]],
                    best_model_feature=ML_output[[9]], nsim=5,
                    ML_method='Random_forest', feature_n=10)
# visualize SHAP feature importance of 10 samples
SHAP_sample_result <- SHAP_sample(SHAP_output[[2]], n_sample=10)


###################################################
### code chunk number 107: LipidSigR.Rnw:2124-2126
###################################################
# view result: SHAP feature importance plot of 10 samples
SHAP_sample_result


###################################################
### code chunk number 108: ML: SHAP analysis_forceplot
###################################################
# get lipid characteristics
char_var <- colnames(lipid_char_table)[-1]
# data processing of machine learning
ML_data <- ML_data_process(exp_data, group_info = condition_table,
                           lipid_char_table, char_var[1],
                           exclude_var_missing=TRUE, missing_pct_limit=50,
                           replace_zero=TRUE, zero2what='min', xmin=0.5,
                           replace_NA=TRUE, NA2what='min', ymin=0.5,
                           pct_transform=TRUE, data_transform=TRUE,
                           trans_type='log', centering=FALSE, scaling=FALSE)
# conduct machine learning
ML_output <- ML_final(ML_data[[2]], ranking_method='Random_forest',
                      ML_method='Random_forest', split_prop=0.3, nfold=10)
# conduct SHAP
SHAP_output <- SHAP(ML_data[[2]], best_model=ML_output[[8]],
                    best_model_feature=ML_output[[9]],
                    ML_method='Random_forest', feature_n=10, nsim=5)
# visualize each predictor’s attributions
SHAP_force_result <- SHAP_forceplot(SHAP_output[[1]], topN_feature=10,
                                    cluster_method="ward.D", group_num=10)

# view result: data frame of force plot information
head(SHAP_force_result[[1]][, 1:5], 5)


###################################################
### code chunk number 109: LipidSigR.Rnw:2164-2166
###################################################
# view result: SHAP force plot
SHAP_force_result[[2]]


###################################################
### code chunk number 110: ML: SHAP analysis_dependence_plot
###################################################
# get lipid characteristics
char_var <- colnames(lipid_char_table)[-1]
# data processing of machine learning
ML_data <- ML_data_process(exp_data, group_info = condition_table,
                           lipid_char_table, char_var[1],
                           exclude_var_missing=TRUE, missing_pct_limit=50,
                           replace_zero=TRUE, zero2what='min', xmin=0.5,
                           replace_NA=TRUE, NA2what='min', ymin=0.5,
                           pct_transform=TRUE, data_transform=TRUE,
                           trans_type='log', centering=FALSE, scaling=FALSE)
# conduct machine learning
ML_output <- ML_final(ML_data[[2]], ranking_method='Random_forest',
                      ML_method='Random_forest', split_prop=0.3, nfold=10)
# conduct SHAP
SHAP_output <- SHAP(ML_data[[2]], best_model=ML_output[[8]],
                    best_model_feature=ML_output[[9]], nsim=5,
                    ML_method='Random_forest', feature_n=10)
# visualize SHAP values against feature values for each variable
SHAP_depend_result <- SHAP_dependence_plot(SHAP_output[[2]], x="C38.6.PC",
                                           y="C38.6.PC", color_var="C38.6.PC")



###################################################
### code chunk number 111: LipidSigR.Rnw:2202-2204
###################################################
# view result: SHAP dependence plot
SHAP_depend_result


###################################################
### code chunk number 112: ML: network (eval = FALSE)
###################################################
## # get lipid characteristics
## char_var <- colnames(lipid_char_table)[-1]
## # data processing of machine learning
## ML_data <- ML_data_process(exp_data, group_info = condition_table,
##                            lipid_char_table, char_var[1],
##                            exclude_var_missing=TRUE, missing_pct_limit=50,
##                            replace_zero=TRUE, zero2what='min', xmin=0.5,
##                            replace_NA=TRUE, NA2what='min', ymin=0.5,
##                            pct_transform=TRUE, data_transform=TRUE,
##                            trans_type='log', centering=FALSE, scaling=FALSE)
## # conduct machine learning
## ML_output <- ML_final(ML_data[[2]], ranking_method='Random_forest',
##                       ML_method='Random_forest', split_prop=0.3, nfold=10)
## # select feature importance from model of ML_output
## model_net <- model_for_net(ML_data[[2]], ML_method='Random_forest',
##                            varimp_method='Algorithm-based', ML_output[[8]],
##                            ML_output[[9]], feature_num=10, nsim=5)
## # compute correlation coefficients and visualize correlation network
## cor_net_result <- cor_network(ML_data[[1]], lipid_char_table,
##                               model_net[[2]], model_net[[3]],
##                               cor_method='pearson', edge_cutoff=0)


###################################################
### code chunk number 113: LipidSigR.Rnw:2243-2245 (eval = FALSE)
###################################################
## # view result: the network of feature importance
## cor_net_result


###################################################
### code chunk number 114: load_correlation_data
###################################################
# clears all objects from workspace
rm(list = ls())

# lipid expression data
data("corr_exp_data")
exp_data <- corr_exp_data
head(exp_data[, 1:5], 5)

# lipid characteristics table
data("corr_lipid_char_table")
lipid_char_table <- corr_lipid_char_table
head(lipid_char_table, 5)

# condition table (clinical factor)
data("corr_condition_table")
condition_table <- corr_condition_table
head(condition_table,  5)

# adjusted table
data("corr_adjusted_table")
adjusted_table <- corr_adjusted_table
head(adjusted_table, 5)


###################################################
### code chunk number 115: corrrelation_data_process
###################################################
# lipid expression data
head(exp_data[, 1:5], 5)
# data processing of exp_data
exp_transform_table <- data_process(exp_data, exclude_var_missing=TRUE,
                                    missing_pct_limit=50, replace_zero=TRUE,
                                    zero2what='min', xmin=0.5,
                                    replace_NA=TRUE, NA2what='min',
                                    ymin=0.5, pct_transform=TRUE,
                                    data_transform=TRUE, trans_type='log',
                                    centering=FALSE, scaling=FALSE)
# exp_data after data processing
head(exp_transform_table[, 1:5], 5)


###################################################
### code chunk number 116: Correlation_lipid species: correlation
###################################################
# data processing of exp_data
exp_transform <- data_process(exp_data, exclude_var_missing=TRUE,
                              missing_pct_limit=50,
                              replace_zero=TRUE, zero2what='min',
                              xmin=0.5, replace_NA=TRUE,
                              NA2what='min', ymin=0.5,
                              pct_transform=TRUE, data_transform=TRUE,
                              trans_type='log', centering=FALSE, scaling=FALSE)
# compute correlation coefficient and visualize by heatmap
COspec_clinCor <- Clin_Cor_heatmap(exp_transform, condition_table,
                                   test = 'pearson', adjust_p_method = 'BH',
                                   sig_stat = 'p.adj', sig_pvalue=1,
                                   sig_cor_coef=0, heatmap_col='statistic',
                                   distfun='spearman', hclustfun='average')

# view result: data frame of clinical features and lipid species
head(COspec_clinCor$Cor_table_all[, 1:5], 5)

# view result: data frame of significant clinical features and lipid species
head(COspec_clinCor$Cor_table_sig[, 1:5], 5)

# view result: clinical features and lipid species correlation reorder matrix
head(COspec_clinCor$Cor_reorder_mat[, 1:2])


###################################################
### code chunk number 117: LipidSigR.Rnw:2364-2366
###################################################
# view result: heatmap of clinical features and lipid species
COspec_clinCor$Cor_table_plot


###################################################
### code chunk number 118: Correlation_lipid species: linear regression
###################################################
# data processing of exp_data
exp_transform <- data_process(exp_data, exclude_var_missing=TRUE,
                              missing_pct_limit=50, replace_zero=TRUE,
                              zero2what='min', xmin=0.5, replace_NA=TRUE,
                              NA2what='min', ymin=0.5, pct_transform=TRUE,
                              data_transform=TRUE, trans_type='log',
                              centering=FALSE, scaling=FALSE)
# compute linear regression and visualize by heatmap
COspec_clin_LR <- Clin_LR_heatmap(exp_transform, condition_table,
                                  adjusted_table, adjust_p_method = 'BH',
                                  sig_stat = 'p.adj', sig_pvalue = 1,
                                  distfun='spearman', hclustfun='centroid',
                                  heatmap_col='beta_coef')

# view result: data frame of statistical results
head(COspec_clin_LR$LR_table_all[, 1:4], 5)

# view result: data frame of significant statistical results
head(COspec_clin_LR$LR_table_sig[, 1:4], 5)

# view result: matrix of heatmap
head(COspec_clin_LR$LR_reorder_mat[, 1:2])


###################################################
### code chunk number 119: LipidSigR.Rnw:2403-2405
###################################################
# view result: heatmap of linear regression
COspec_clin_LR$LR_table_plot


###################################################
### code chunk number 120: Correlation_lipid characteristics: correlation
###################################################
# get lipid characteristics
char_var <- colnames(lipid_char_table)[-1]
# aggregated(sum) expression data by selected characteristics
exp_data_Spe2Char <- Species2Char(exp_data, lipid_char_table,
                                  char_var = char_var[1])
# data processing of exp_data_Spe2Char
exp_transform_class <- data_process(exp_data_Spe2Char,
                                    exclude_var_missing=TRUE,
                                    missing_pct_limit=50, replace_zero=TRUE,
                                    zero2what='NA', xmin=0.5, replace_NA=TRUE,
                                    NA2what='min', ymin=0.5,
                                    pct_transform=TRUE,
                                    data_transform=FALSE, trans_type='log',
                                    centering=FALSE, scaling=FALSE)
# compute correlation coefficient and visualize by heatmap
COchar_clinCor <- Clin_Cor_heatmap(exp_transform_class, condition_table,
                                   test = 'pearson', adjust_p_method = 'BH',
                                   sig_stat = 'p.adj', sig_pvalue=1,
                                   sig_cor_coef=0, heatmap_col='statistic',
                                   distfun='spearman', hclustfun='average')

# view result: data frame of clinical features and lipid characteristics
head(COchar_clinCor$Cor_table_all[, 1:5], 5)

# view result: data frame of significant clinical features and lipid characteristics
head(COchar_clinCor$Cor_table_sig[, 1:5], 5)

# view result: clinical features and lipid characteristics correlation reorder matrix
head(COchar_clinCor$Cor_reorder_mat[, 1:2])


###################################################
### code chunk number 121: LipidSigR.Rnw:2472-2474
###################################################
# view result: heatmap of clinical features and lipid characteristics
COchar_clinCor$Cor_table_plot


###################################################
### code chunk number 122: Correlation_lipid characteristics: linear regression
###################################################
# get lipid characteristics
char_var <- colnames(lipid_char_table)[-1]
# aggregated(sum) expression data by selected characteristics
exp_data_Spe2Char <- Species2Char(exp_data, lipid_char_table,
                                  char_var = char_var[1])
# data processing of exp_data_Spe2Char
exp_transform_class <- data_process(exp_data_Spe2Char,
                                    exclude_var_missing=TRUE,
                                    missing_pct_limit=50, replace_zero=TRUE,
                                    zero2what='NA', xmin=0.5, replace_NA=TRUE,
                                    NA2what='min', ymin=0.5,
                                    pct_transform=TRUE,
                                    data_transform=FALSE, trans_type='log',
                                    centering=FALSE, scaling=FALSE)
# compute linear regression and visualize by heatmap
COchar_clin_LR <- Clin_LR_heatmap(exp_transform_class, condition_table,
                                  adjusted_table, adjust_p_method = 'BH',
                                  sig_stat = 'p.adj', sig_pvalue = 1,
                                  distfun='spearman', hclustfun='centroid',
                                  heatmap_col='beta_coef')

# view result: data frame of statistical results
head(COchar_clin_LR$LR_table_all[, 1:4], 5)

# view result: data frame of significant statistical results
head(COchar_clin_LR$LR_table_sig[, 1:4], 5)

# view result: matrix of heatmap
head(COchar_clin_LR$LR_reorder_mat[, 1:2])


###################################################
### code chunk number 123: LipidSigR.Rnw:2519-2521
###################################################
# view result: heatmap of linear regression
COchar_clin_LR$LR_table_plot


###################################################
### code chunk number 124: sessionInfo
###################################################
toLatex(sessionInfo())


