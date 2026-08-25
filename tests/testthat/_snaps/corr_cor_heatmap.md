# corr_cor_heatmap can handle all error inputs.

    Code
      corr_cor_heatmap(processed_se, char = "char", condition_col = c("FEV1_FVC",
        "Emphysema", "Exacerbations"), side_color_char = NULL, correlation = "pearson",
      significant = "pval", p_cutoff = 1, adjust_p_method = "BH", cor_coef_cutoff = 0,
      distfun = "spearman", hclustfun = "average", heatmap_col = "statistic",
      transform = "log10", type = "Char")
    Condition
      Error in `corr_cor_heatmap()`:
      ! Wrong char input, you can view the available char list by list_lipid_char function.

---

    Code
      corr_cor_heatmap(processed_se, char = NULL, condition_col = c("FEV1_FVC",
        "Emphysema", "Exacerbations"), side_color_char = NULL, correlation = "other",
      significant = "pval", p_cutoff = 1, adjust_p_method = "BH", cor_coef_cutoff = 0,
      distfun = "spearman", hclustfun = "average", heatmap_col = "statistic",
      transform = "log10", type = "Sp")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "pearson", "spearman", "kendall"

---

    Code
      corr_cor_heatmap(processed_se, char = NULL, condition_col = c("FEV1_FVC",
        "Emphysema", "Exacerbations"), side_color_char = NULL, correlation = "pearson",
      significant = "p", p_cutoff = 1, adjust_p_method = "BH", cor_coef_cutoff = 0,
      distfun = "spearman", hclustfun = "average", heatmap_col = "statistic",
      transform = "log10", type = "Sp")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "pval", "padj"

---

    Code
      corr_cor_heatmap(processed_se, char = NULL, condition_col = c("FEV1_FVC",
        "Emphysema", "Exacerbations"), side_color_char = NULL, correlation = "pearson",
      significant = "pval", p_cutoff = NULL, adjust_p_method = "BH", cor_coef_cutoff = 0,
      distfun = "spearman", hclustfun = "average", heatmap_col = "statistic",
      transform = "log10", type = "Sp")
    Condition
      Error in `corr_cor_heatmap()`:
      ! p_cutoff must be a numeric value between 0 and 1.

---

    Code
      corr_cor_heatmap(processed_se, char = NULL, condition_col = c("FEV1_FVC",
        "Emphysema", "Exacerbations"), side_color_char = NULL, correlation = "pearson",
      significant = "pval", p_cutoff = 1, adjust_p_method = "other", cor_coef_cutoff = 0,
      distfun = "spearman", hclustfun = "average", heatmap_col = "statistic",
      transform = "log10", type = "Sp")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", "none"

---

    Code
      corr_cor_heatmap(processed_se, char = NULL, condition_col = c("FEV1_FVC",
        "Emphysema", "Exacerbations"), side_color_char = NULL, correlation = "pearson",
      significant = "pval", p_cutoff = 1, adjust_p_method = "BH", cor_coef_cutoff = 2,
      distfun = "spearman", hclustfun = "average", heatmap_col = "statistic",
      transform = "log10", type = "Sp")
    Condition
      Error in `corr_cor_heatmap()`:
      ! cor_coef_cutoff must be a numeric value between 0 and 1.

---

    Code
      corr_cor_heatmap(processed_se, char = NULL, condition_col = c("FEV1_FVC",
        "Emphysema", "Exacerbations"), side_color_char = NULL, correlation = "pearson",
      significant = "pval", p_cutoff = 1, adjust_p_method = "BH", cor_coef_cutoff = 0,
      distfun = NULL, hclustfun = "average", heatmap_col = "statistic", transform = "log10",
      type = "Sp")
    Condition
      Error in `if (distfun %in% c("pearson", "kendall", "spearman")) ...`:
      ! argument is of length zero

---

    Code
      corr_cor_heatmap(processed_se, char = NULL, condition_col = c("FEV1_FVC",
        "Emphysema", "Exacerbations"), side_color_char = NULL, correlation = "pearson",
      significant = "pval", p_cutoff = 1, adjust_p_method = "BH", cor_coef_cutoff = 0,
      distfun = "spearman", hclustfun = "Average", heatmap_col = "statistic",
      transform = "log10", type = "Sp")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "average", "ward.D", "ward.D2", "single", "complete", "mcquitty", "median", "centroid"

---

    Code
      corr_cor_heatmap(processed_se, char = NULL, condition_col = c("FEV1_FVC",
        "Emphysema", "Exacerbations"), side_color_char = NULL, correlation = "pearson",
      significant = "pval", p_cutoff = 1, adjust_p_method = "BH", cor_coef_cutoff = 0,
      distfun = "spearman", hclustfun = "average", heatmap_col = "t_statistic",
      transform = "log10", type = "Sp")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "statistic", "cor_coef"

---

    Code
      corr_cor_heatmap(processed_se, char = NULL, condition_col = c("FEV1_FVC",
        "Emphysema", "Exacerbations"), side_color_char = NULL, correlation = "pearson",
      significant = "pval", p_cutoff = 1, adjust_p_method = "BH", cor_coef_cutoff = 0,
      distfun = "spearman", hclustfun = "average", heatmap_col = "statistic",
      transform = "log2", type = "Sp")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "log10", "square", "cube", "none"

---

    Code
      corr_cor_heatmap(processed_se, char = NULL, condition_col = c("FEV1_FVC",
        "Emphysema", "Exacerbations"), side_color_char = NULL, correlation = "pearson",
      significant = "pval", p_cutoff = 1, adjust_p_method = "BH", cor_coef_cutoff = 0,
      distfun = "spearman", hclustfun = "average", heatmap_col = "statistic",
      transform = "log10", type = "other")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "Sp", "Char"

---

    Code
      corr_cor_heatmap(processed_se, char = "class", condition_col = c("FEV1_FVC",
        "Emphysema", "Exacerbations"), side_color_char = NULL, correlation = "pearson",
      significant = "pval", p_cutoff = 1, adjust_p_method = "BH", cor_coef_cutoff = 0,
      distfun = "spearman", hclustfun = "average", heatmap_col = "statistic",
      transform = "log10", type = "Sp")
    Condition
      Error in `corr_cor_heatmap()`:
      ! If type is set to 'Sp' for running lipid species correlation analysis, char must be set to NULL.

---

    Code
      corr_cor_heatmap(processed_se, char = NULL, condition_col = c("FEV1_FVC",
        "Emphysema", "Exacerbations"), side_color_char = "class", correlation = "pearson",
      significant = "pval", p_cutoff = 1, adjust_p_method = "BH", cor_coef_cutoff = 0,
      distfun = "spearman", hclustfun = "average", heatmap_col = "statistic",
      transform = "log10", type = "Char")
    Condition
      Error in `corr_cor_heatmap()`:
      ! If type is set to 'Char' for running lipid characteristics correlation analysis, side_color_char must be set to NULL.

---

    Code
      corr_cor_heatmap(processed_se, char = NULL, condition_col = c("FEV1_FVC",
        "Emphysema", "Exacerbations"), side_color_char = "char", correlation = "pearson",
      significant = "pval", p_cutoff = 1, adjust_p_method = "BH", cor_coef_cutoff = 0,
      distfun = "spearman", hclustfun = "average", heatmap_col = "statistic",
      transform = "log10", type = "Sp")
    Condition
      Error in `corr_cor_heatmap()`:
      ! Wrong side_color_char input, you can view the available side_color_char list by list_lipid_char function.

# corr_cor_heatmap can handle condition_col and adjusted_col inputs.

    Code
      corr_cor_heatmap(processed_se, char = NULL, condition_col = c("ABC",
        "Emphysema", "Exacerbations"), side_color_char = "class", correlation = "pearson",
      significant = "pval", p_cutoff = 1, adjust_p_method = "BH", cor_coef_cutoff = 0,
      distfun = "spearman", hclustfun = "average", heatmap_col = "statistic",
      transform = "log10", type = "Sp")
    Condition
      Error in `.check_df_corrCol()`:
      ! The condition_col must be column names selected from FEV1_FVC, Emphysema, Exacerbations, Age, Sex, Smoking, BMI, FEV1.

---

    Code
      corr_cor_heatmap(processed_se, char = NULL, condition_col = c("Exacerbations"),
      side_color_char = "class", correlation = "pearson", significant = "pval",
      p_cutoff = 1, adjust_p_method = "BH", cor_coef_cutoff = 0, distfun = "spearman",
      hclustfun = "average", heatmap_col = "statistic", transform = "log10", type = "Sp")
    Condition
      Error in `.check_df_corrCol()`:
      ! The condition_col must include at least two columns from FEV1_FVC, Emphysema, Exacerbations, Age, Sex, Smoking, BMI, FEV1.

---

    Code
      corr_cor_heatmap(processed_se, char = NULL, condition_col = NULL,
        side_color_char = "class", correlation = "pearson", significant = "pval",
        p_cutoff = 1, adjust_p_method = "BH", cor_coef_cutoff = 0, distfun = "spearman",
        hclustfun = "average", heatmap_col = "statistic", transform = "log10", type = "Sp")
    Condition
      Error in `.check_df_corrCol()`:
      ! The condition_col parameter cannot be null. Select at least two columns from FEV1_FVC, Emphysema, Exacerbations, Age, Sex, Smoking, BMI, FEV1.

---

    Code
      corr_cor_heatmap(corr_data_err, char = NULL, condition_col = c("FEV1_FVC",
        "Emphysema", "Exacerbations"), side_color_char = "class", correlation = "pearson",
      significant = "pval", p_cutoff = 1, adjust_p_method = "BH", cor_coef_cutoff = 0,
      distfun = "spearman", hclustfun = "average", heatmap_col = "statistic",
      transform = "log10", type = "Sp")
    Condition
      Error in `.check_df_corrCol()`:
      ! Except for the first column, sample_name, the group information table must include at least two columns with numeric values and no missing values.

