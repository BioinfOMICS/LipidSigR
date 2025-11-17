# LipidSigR 1.0.4
## Minor bug fixes and improvements
- Fix the deprecated syntax in the tidyverse and ggplot2 packages.
- Resolve unit test errors and warnings.
- Change the static heatmap dependency package to ComplexHeatmap.
- Fixed an issue in the `heatmap_clustering()` function where the `char` parameter selection did not provide a colour bar on the side of the heatmap. The hover information has also been fixed.
- The clustering method in `heatmap_correlation()` has been corrected to ensure consistent sorting for interactive and static heatmaps.


# LipidSigR 1.0.3

## Minor bug fixes and improvements
- `char_2wayAnova()` Remove check imputation.
- `dr_pca()`, `dr_tsne()`, `dr_umap()`, and `dr_plsda()` Modify the description of the clustering parameter.
- `.deChar_plot_tab_multiGroup()`, and `.table_deChar_twoGroup()` Improved recognition of special characters in regular expressions.

# LipidSigR 1.0.2

## Minor bug fixes and improvements
- `ml_model()` The problem of lipid names containing semicolon not being able to be modeled has now been fixed.
- `dr_pca()`, `dr_tsne()`, `dr_umap()`, and `dr_plsda()` Added check cluster_num.
- `nw_pathway_activity()` add feedback when there were no reactions in your dataset.
- `nw_gatom()` The error that was appearing at times has now been fixed.
- `enrichment_lsea()` add feedback when there was no enrichment in lipid characteristics.
- `enrichment_ora()` The error that was appearing at times has now been fixed.

# LipidSigR 1.0.1

## Minor bug fixes and improvements
- `data_process()` add parameter transform and add return metadata in SummarizedExperiment.
- `convert_sp2char()` The error that was appearing at times has now been fixed.
- `heatmap_chain_db()` The error that was appearing at times has now been fixed.

# LipidSigR 1.0.0
- Update Correlation workflow

# LipidSigR 0.9.0
- Update Machine learning workflow

# LipidSigR 0.7.0
- Released on Github
