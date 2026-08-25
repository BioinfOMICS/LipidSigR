# LipidSigR 1.0.5

## New Features

### Multi-Omics Integration Workflow (new)
- `as_MultiAssayExperiment()` (new, `utils-MultiAssay.R`): combines a
  lipidomics `SummarizedExperiment` with a second, non-lipidomics
  `SummarizedExperiment` (e.g. RNA-seq) into one `MultiAssayExperiment`
  for downstream multi-omics analysis.
- `as_summarized_experiment()`: added a new `se_type = "rna"` option for
  building the second (non-lipidomics) SE consumed by 
  `as_MultiAssayExperiment()`. `group_info` format is auto-detected as either 
  the `de_two`-style (4 columns, paired design inferred from the `pair` column) 
  or `de_multiple`-style (3 columns) layout.

### Correlation Workflow - Three New Heatmap Functions
- `corr_DE_heatmap()` (new): differential-expression-based correlation
  heatmap between clinical features and lipid species/characteristics.
  Supports parametric (t-test/ANOVA) and non-parametric
  (Wilcoxon/Kruskal-Wallis) tests.
- `corr_logR_heatmap()` (new): logistic-regression-based correlation
  heatmap, with support for adjusted covariates; heatmap can be colored by
  ln(OR) or z-statistic.
- `corr_roc_heatmap()` (new): ROC-analysis-based correlation heatmap,
  computing AUC with 95% CI for each lipid species/characteristic against a
  binary clinical condition.

### Dimensionality Reduction
- `dr_pcoa()` (new): Principal Coordinates Analysis, via
  `vegan::vegdist`/`ape`, with multiple distance methods, eigenvalue
  corrections (none/lingoes/cailliez), and the existing clustering options
  (kmeans, kmedoids, hclustering, dbscan, group_info).
- `dr_splsda()` (new): Sparse PLS-DA with LASSO-style feature selection;
  component number and feature selection are tuned automatically via
  `mixOmics::tune.splsda`.
- `dr_oplsda()` (new): Orthogonal PLS-DA (via `ropls::opls`), separating
  predictive from orthogonal variation; the orthogonal component count is
  auto-determined by cross-validation.

### Data Processing (`data_process()`)
- Added `exclude_missing_method` (`"global"` or `"group"`) to control
  whether the missing-value filter is applied across all samples or
  per-group.
- Added new normalization methods: `CubicSplines` (`affy`), `CyclicLoess`
  (`limma`), `EigenMS` (`ProteoMM`), `linear`, `Mean`, `nonLinear`, `RLR`
  (`MASS`), and `TMM` (`edgeR`).
- Added new imputation methods: `BPCA`, `GMS`, `impSeq`, `Impseqrob`,
  `LLS`, `MinProb`, and `PI`.
- Added a `"log2"` option for `transform`.
- **Removed**: the `"IRMI"` imputation method (`VIM::irmi`) has been
  dropped - pipelines relying on `IRMI` will need to switch methods.

### Lipid Characteristic Conversion
- `convert_sp2char()`: when splitting `FA`/`FA.C`/`FA.DB`/`FA.OH`
  characteristics that contain `"|"`, abundance values are now divided by
  the number of split components so totals stay quantitatively correct.

### Differential Expression
- `deChar_twoGroup()` / `deSp_twoGroup()`: added `"limma"` as a new `test`
  option (via `limma::lmFit`/`eBayes`/`topTable`), alongside `"t-test"` and
  `"Wilcoxon test"`.
  
## Parameter / Default Changes
- Across most exported functions, parameter validation has been refactored
  to use `match.arg()`, and as part of this several functions' default
  argument values were also updated: `corr_cor_heatmap()`,
  `corr_lr_heatmap()`, `heatmap_clustering()`, `heatmap_correlation()`,
  `enrichment_lsea()`, `nw_lipid_reaction()`, `ml_model()`, `plot_ml_pr()`,
  `plot_ml_probability()`, and `plot_ml_evaluation()`.
- `ml_shap()`: the allowed range for `nsim` was also adjusted.

## Bug Fixes

- `ml_model()`: Fixed the train/test split assignment, and corrected the
  argument order passed to `caret::confusionMatrix()`. `plot_ml_probability()`
  received the same `confusionMatrix()` fix.
- `enrichment_ora()`: Fixed an issue in how the background lipid count for
  the enrichment test was calculated, which could affect resulting p-values.
- `utils-ml.R` / `plot_ml_probability()` / `plot_ml_feature()`: Fixed a
  filtering issue so that the `feature_num` selection is now applied as
  intended.
- `boxPlot_feature_multiGroup()` / `boxPlot_feature_twoGroup()`: Fixed an
  issue where the "square" transform option wasn't being applied correctly,
  and corrected a column reference used in the paired-test calculation.
- `char_association()`: Fixed a few calculation and plotting issues
  affecting the fold-change direction indicator, error bars, and hover
  text.
- `subChar_twoGroup()` / `deChar_twoGroup()`: Fixed an issue where the
  significance column wasn't being read correctly when flagging
  significant results.
- `heatmap_correlation()` / `heatmap_clustering()`: Fixed issues in the
  pre-clustering data checks that could lead to missing values or
  duplicate values not being handled correctly.
- `utils-heatmap.R`: Added clearer error handling for edge cases (e.g. all
  values equal to 0, or zero-variance correlation rows) that previously
  produced uninformative failures.
- Corrected several inaccurate error messages that referenced the wrong
  function name or group type.
- `ml_corr_network()`: Fixed a couple of minor display issues affecting
  node labels/shapes and the legend for smaller networks.

---

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
