# corr_logR_heatmap can handle error inputs

    Code
      corr_logR_heatmap(processed_se, char = "char", condition_col = c("FEV1_FVC",
        "Emphysema"), adjusted_col = NULL, side_color_char = NULL, significant = "pval",
      p_cutoff = 0.05, adjust_p_method = "BH", distfun = "spearman", hclustfun = "centroid",
      heatmap_col = "z_statistic", transform = "log10", type = "Char")
    Message
      No adjusted variables is selected.
    Condition
      Error in `corr_logR_heatmap()`:
      ! Wrong char input, you can view the available char list by list_lipid_char function.

---

    Code
      corr_logR_heatmap(processed_se, char = NULL, condition_col = c("FEV1_FVC",
        "Emphysema"), adjusted_col = NULL, side_color_char = NULL, significant = "p",
      p_cutoff = 0.05, adjust_p_method = "BH", distfun = "spearman", hclustfun = "centroid",
      heatmap_col = "z_statistic", transform = "log10", type = "Sp")
    Message
      No adjusted variables is selected.
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "pval", "padj"

---

    Code
      corr_logR_heatmap(processed_se, char = NULL, condition_col = c("FEV1_FVC",
        "Emphysema"), adjusted_col = c("Age", "Sex", "Smoking", "BMI", "FEV1"),
      side_color_char = NULL, significant = "pval", p_cutoff = 1.2, adjust_p_method = "BH",
      distfun = "spearman", hclustfun = "centroid", heatmap_col = "z_statistic",
      transform = "log10", type = "Sp")
    Condition
      Error in `corr_logR_heatmap()`:
      ! p_cutoff must be a numeric value between 0 and 1.

---

    Code
      corr_logR_heatmap(processed_se, char = NULL, condition_col = c("FEV1_FVC",
        "Emphysema"), adjusted_col = NULL, side_color_char = NULL, significant = "pval",
      p_cutoff = 0.05, adjust_p_method = "other", distfun = "spearman", hclustfun = "centroid",
      heatmap_col = "z_statistic", transform = "log10", type = "Sp")
    Message
      No adjusted variables is selected.
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", "none"

---

    Code
      corr_logR_heatmap(processed_se, char = NULL, condition_col = c("FEV1_FVC",
        "Emphysema"), adjusted_col = NULL, side_color_char = NULL, significant = "pval",
      p_cutoff = 0.05, adjust_p_method = "BH", distfun = "Spearman", hclustfun = "centroid",
      heatmap_col = "z_statistic", transform = "log10", type = "Sp")
    Message
      No adjusted variables is selected.
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "spearman", "pearson", "kendall"

---

    Code
      corr_logR_heatmap(processed_se, char = NULL, condition_col = c("FEV1_FVC",
        "Emphysema"), adjusted_col = NULL, side_color_char = NULL, significant = "pval",
      p_cutoff = 0.05, adjust_p_method = "BH", distfun = "spearman", hclustfun = "Centroid",
      heatmap_col = "z_statistic", transform = "log10", type = "Sp")
    Message
      No adjusted variables is selected.
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "average", "ward.D", "ward.D2", "single", "complete", "mcquitty", "median", "centroid"

---

    Code
      corr_logR_heatmap(processed_se, char = NULL, condition_col = c("FEV1_FVC",
        "Emphysema"), adjusted_col = NULL, side_color_char = NULL, significant = "pval",
      p_cutoff = 0.05, adjust_p_method = "BH", distfun = "spearman", hclustfun = "centroid",
      heatmap_col = "statistic", transform = "log10", type = "Sp")
    Message
      No adjusted variables is selected.
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "ln_OR", "z_statistic"

---

    Code
      corr_logR_heatmap(processed_se, char = NULL, condition_col = c("FEV1_FVC",
        "Emphysema"), adjusted_col = NULL, side_color_char = NULL, significant = "pval",
      p_cutoff = 0.05, adjust_p_method = "BH", distfun = "spearman", hclustfun = "centroid",
      heatmap_col = "z_statistic", transform = "log2", type = "Sp")
    Message
      No adjusted variables is selected.
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "log10", "square", "cube", "none"

---

    Code
      corr_logR_heatmap(processed_se, char = NULL, condition_col = c("FEV1_FVC",
        "Emphysema"), adjusted_col = NULL, side_color_char = NULL, significant = "pval",
      p_cutoff = 0.05, adjust_p_method = "BH", distfun = "spearman", hclustfun = "centroid",
      heatmap_col = "z_statistic", transform = "log10", type = "other")
    Message
      No adjusted variables is selected.
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "Sp", "Char"

---

    Code
      corr_logR_heatmap(processed_se, char = "class", condition_col = c("FEV1_FVC",
        "Emphysema"), adjusted_col = NULL, side_color_char = NULL, significant = "pval",
      p_cutoff = 0.05, adjust_p_method = "BH", distfun = "spearman", hclustfun = "centroid",
      heatmap_col = "z_statistic", transform = "log10", type = "Sp")
    Message
      No adjusted variables is selected.
    Condition
      Error in `corr_logR_heatmap()`:
      ! If type is set to 'Sp' for running lipid species correlation analysis, char must be set to NULL.

---

    Code
      corr_logR_heatmap(processed_se, char = NULL, condition_col = c("FEV1_FVC",
        "Emphysema"), adjusted_col = NULL, side_color_char = "char", significant = "pval",
      p_cutoff = 0.05, adjust_p_method = "BH", distfun = "spearman", hclustfun = "centroid",
      heatmap_col = "z_statistic", transform = "log10", type = "Sp")
    Message
      No adjusted variables is selected.
    Condition
      Error in `corr_logR_heatmap()`:
      ! Wrong side_color_char input, you can view the available side_color_char list by list_lipid_char function.

---

    Code
      corr_logR_heatmap(processed_se, char = NULL, condition_col = c("FEV1_FVC",
        "Emphysema"), adjusted_col = NULL, side_color_char = "char", significant = "pval",
      p_cutoff = 0.05, adjust_p_method = "BH", distfun = "spearman", hclustfun = "centroid",
      heatmap_col = "z_statistic", transform = "log10", type = "Char")
    Message
      No adjusted variables is selected.
    Condition
      Error in `corr_logR_heatmap()`:
      ! If type is set to 'Char' for running lipid characteristics correlation analysis, side_color_char must be set to NULL.

# corr_logR_heatmap can handle condition_col and adjusted_col inputs.

    Code
      corr_logR_heatmap(processed_se, char = NULL, condition_col = c("ABC",
        "Emphysema"), adjusted_col = NULL, side_color_char = NULL, significant = "pval",
      p_cutoff = 0.05, adjust_p_method = "BH", distfun = "spearman", hclustfun = "centroid",
      heatmap_col = "z_statistic", transform = "log10", type = "Sp")
    Condition
      Error in `.check_df_corrCol()`:
      ! The condition_col must be column names selected from FEV1_FVC, Emphysema, Exacerbations, Age, Sex, Smoking, BMI, FEV1.

---

    Code
      corr_logR_heatmap(processed_se, char = NULL, condition_col = c("Emphysema"),
      adjusted_col = NULL, side_color_char = NULL, significant = "pval", p_cutoff = 0.05,
      adjust_p_method = "BH", distfun = "spearman", hclustfun = "centroid",
      heatmap_col = "z_statistic", transform = "log10", type = "Sp")
    Condition
      Error in `.check_df_corrCol()`:
      ! The condition_col must include at least two columns from FEV1_FVC, Emphysema, Exacerbations, Age, Sex, Smoking, BMI, FEV1.

---

    Code
      corr_logR_heatmap(processed_se, char = NULL, condition_col = NULL,
        adjusted_col = NULL, side_color_char = NULL, significant = "pval", p_cutoff = 0.05,
        adjust_p_method = "BH", distfun = "spearman", hclustfun = "centroid",
        heatmap_col = "z_statistic", transform = "log10", type = "Sp")
    Condition
      Error in `.check_df_corrCol()`:
      ! The condition_col parameter cannot be null. Select at least two columns from FEV1_FVC, Emphysema, Exacerbations, Age, Sex, Smoking, BMI, FEV1.

---

    Code
      corr_logR_heatmap(corr_data_err, char = NULL, condition_col = c("FEV1_FVC",
        "Emphysema"), adjusted_col = NULL, side_color_char = NULL, significant = "pval",
      p_cutoff = 0.05, adjust_p_method = "BH", distfun = "spearman", hclustfun = "centroid",
      heatmap_col = "z_statistic", transform = "log10", type = "Sp")
    Condition
      Error in `.check_df_corrCol()`:
      ! Except for the first column, sample_name, the group information table must include at least two columns with numeric values and no missing values.

---

    Code
      corr_logR_heatmap(processed_se, char = NULL, condition_col = c("FEV1_FVC",
        "Emphysema"), adjusted_col = c("FEV1_FVC", "Age", "Sex", "Smoking", "BMI"),
      side_color_char = NULL, significant = "pval", p_cutoff = 0.05, adjust_p_method = "BH",
      distfun = "spearman", hclustfun = "centroid", heatmap_col = "z_statistic",
      transform = "log10", type = "Sp")
    Condition
      Error in `.check_df_corrCol()`:
      ! The condition_col and adjusted_col must not have overlapping columns.

---

    Code
      corr_logR_heatmap(processed_se, char = NULL, condition_col = c("FEV1_FVC",
        "Emphysema"), adjusted_col = c("ABC", "Sex", "Smoking", "BMI"),
      side_color_char = NULL, significant = "pval", p_cutoff = 0.05, adjust_p_method = "BH",
      distfun = "spearman", hclustfun = "centroid", heatmap_col = "z_statistic",
      transform = "log10", type = "Sp")
    Condition
      Error in `.check_df_corrCol()`:
      ! The adjusted_col must be column names selected from FEV1_FVC, Emphysema, Exacerbations, Age, Sex, Smoking, BMI, FEV1.

# corr_logR_heatmap can handle all adjust_p_method.

    Code
      corr_logR_heatmap(processed_se, char = NULL, condition_col = c("Sex", "Smoking"),
      adjusted_col = NULL, side_color_char = NULL, significant = "padj", p_cutoff = 1,
      adjust_p_method = adjust_p, distfun = "spearman", hclustfun = "centroid",
      heatmap_col = "z_statistic", transform = "log10", type = "Sp")
    Message
      No adjusted variables is selected.
    Condition
      Error in `.corr_heatmap()`:
      ! Insufficient data for plotting heatmap.

---

    Code
      corr_logR_heatmap(processed_se, char = NULL, condition_col = c("Sex", "Smoking"),
      adjusted_col = NULL, side_color_char = NULL, significant = "padj", p_cutoff = 1,
      adjust_p_method = adjust_p, distfun = "spearman", hclustfun = "centroid",
      heatmap_col = "z_statistic", transform = "log10", type = "Sp")
    Message
      No adjusted variables is selected.
    Condition
      Error in `.corr_heatmap()`:
      ! Insufficient data for plotting heatmap.

---

    Code
      corr_logR_heatmap(processed_se, char = NULL, condition_col = c("Sex", "Smoking"),
      adjusted_col = NULL, side_color_char = NULL, significant = "padj", p_cutoff = 1,
      adjust_p_method = adjust_p, distfun = "spearman", hclustfun = "centroid",
      heatmap_col = "z_statistic", transform = "log10", type = "Sp")
    Message
      No adjusted variables is selected.
    Condition
      Error in `.corr_heatmap()`:
      ! Insufficient data for plotting heatmap.

