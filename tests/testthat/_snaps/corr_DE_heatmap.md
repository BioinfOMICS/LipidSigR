# corr_DE_heatmap can handle error inputs

    Code
      corr_DE_heatmap(processed_se, char = "char", condition_col = c("Sex", "Smoking"),
      side_color_char = NULL, significant = "pval", p_cutoff = 0.05, adjust_p_method = "BH",
      distfun = "spearman", hclustfun = "centroid", transform = "log10", type = "Char")
    Message
      No adjusted variables is selected.
    Condition
      Error in `corr_DE_heatmap()`:
      ! Wrong char input, you can view the available char list by list_lipid_char function.

---

    Code
      corr_DE_heatmap(processed_se, char = NULL, condition_col = c("Sex", "Smoking"),
      side_color_char = NULL, significant = "p", p_cutoff = 0.05, adjust_p_method = "BH",
      distfun = "spearman", hclustfun = "centroid", transform = "log10", type = "Sp")
    Message
      No adjusted variables is selected.
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "pval", "padj"

---

    Code
      corr_DE_heatmap(processed_se, char = NULL, condition_col = c("Sex", "Smoking"),
      side_color_char = NULL, significant = "pval", p_cutoff = 1.2, adjust_p_method = "BH",
      distfun = "spearman", hclustfun = "centroid", transform = "log10", type = "Sp")
    Message
      No adjusted variables is selected.
    Condition
      Error in `corr_DE_heatmap()`:
      ! p_cutoff must be a numeric value between 0 and 1.

---

    Code
      corr_DE_heatmap(processed_se, char = NULL, condition_col = c("Sex", "Smoking"),
      side_color_char = NULL, significant = "pval", p_cutoff = 0.05, adjust_p_method = "other",
      distfun = "spearman", hclustfun = "centroid", transform = "log10", type = "Sp")
    Message
      No adjusted variables is selected.
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", "none"

---

    Code
      corr_DE_heatmap(processed_se, char = NULL, condition_col = c("Sex", "Smoking"),
      side_color_char = NULL, significant = "pval", p_cutoff = 0.05, adjust_p_method = "BH",
      distfun = "Spearman", hclustfun = "centroid", transform = "log10", type = "Sp")
    Message
      No adjusted variables is selected.
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "spearman", "pearson", "kendall"

---

    Code
      corr_DE_heatmap(processed_se, char = NULL, condition_col = c("Sex", "Smoking"),
      side_color_char = NULL, significant = "pval", p_cutoff = 0.05, adjust_p_method = "BH",
      distfun = "spearman", hclustfun = "Centroid", transform = "log10", type = "Sp")
    Message
      No adjusted variables is selected.
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "average", "ward.D", "ward.D2", "single", "complete", "mcquitty", "median", "centroid"

---

    Code
      corr_DE_heatmap(processed_se, char = NULL, condition_col = c("Sex", "Smoking"),
      side_color_char = NULL, significant = "pval", p_cutoff = 0.05, adjust_p_method = "BH",
      distfun = "spearman", hclustfun = "centroid", transform = "log2", type = "Sp")
    Message
      No adjusted variables is selected.
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "log10", "square", "cube", "none"

---

    Code
      corr_DE_heatmap(processed_se, char = NULL, condition_col = c("Sex", "Smoking"),
      side_color_char = NULL, significant = "pval", p_cutoff = 0.05, adjust_p_method = "BH",
      distfun = "spearman", hclustfun = "centroid", transform = "log10", type = "other")
    Message
      No adjusted variables is selected.
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "Sp", "Char"

---

    Code
      corr_DE_heatmap(processed_se, char = "class", condition_col = c("Sex",
        "Smoking"), side_color_char = NULL, significant = "pval", p_cutoff = 0.05,
      adjust_p_method = "BH", distfun = "spearman", hclustfun = "centroid",
      transform = "log10", type = "Sp")
    Message
      No adjusted variables is selected.
    Condition
      Error in `corr_DE_heatmap()`:
      ! If type is set to 'Sp' for running lipid species correlation analysis, char must be set to NULL.

---

    Code
      corr_DE_heatmap(processed_se, char = NULL, condition_col = c("Sex", "Smoking"),
      side_color_char = "char", significant = "pval", p_cutoff = 0.05,
      adjust_p_method = "BH", distfun = "spearman", hclustfun = "centroid",
      transform = "log10", type = "Sp")
    Message
      No adjusted variables is selected.
    Condition
      Error in `corr_DE_heatmap()`:
      ! Wrong side_color_char input, you can view the available side_color_char list by list_lipid_char function.

---

    Code
      corr_DE_heatmap(processed_se, char = NULL, condition_col = c("Sex", "Smoking"),
      side_color_char = "char", significant = "pval", p_cutoff = 0.05,
      adjust_p_method = "BH", distfun = "spearman", hclustfun = "centroid",
      transform = "log10", type = "Char")
    Message
      No adjusted variables is selected.
    Condition
      Error in `corr_DE_heatmap()`:
      ! If type is set to 'Char' for running lipid characteristics correlation analysis, side_color_char must be set to NULL.

# corr_DE_heatmap can handle condition_col and adjusted_col inputs.

    Code
      corr_DE_heatmap(processed_se, char = NULL, condition_col = c("ABC", "Emphysema"),
      side_color_char = NULL, significant = "pval", p_cutoff = 0.05, adjust_p_method = "BH",
      distfun = "spearman", hclustfun = "centroid", transform = "log10", type = "Sp")
    Condition
      Error in `.check_df_corrCol()`:
      ! The condition_col must be column names selected from FEV1_FVC, Emphysema, Exacerbations, Age, Sex, Smoking, BMI, FEV1.

---

    Code
      corr_DE_heatmap(processed_se, char = NULL, condition_col = c("Emphysema"),
      side_color_char = NULL, significant = "pval", p_cutoff = 0.05, adjust_p_method = "BH",
      distfun = "spearman", hclustfun = "centroid", transform = "log10", type = "Sp")
    Condition
      Error in `.check_df_corrCol()`:
      ! The condition_col must include at least two columns from FEV1_FVC, Emphysema, Exacerbations, Age, Sex, Smoking, BMI, FEV1.

---

    Code
      corr_DE_heatmap(processed_se, char = NULL, condition_col = NULL,
        side_color_char = NULL, significant = "pval", p_cutoff = 0.05,
        adjust_p_method = "BH", distfun = "spearman", hclustfun = "centroid",
        transform = "log10", type = "Sp")
    Condition
      Error in `.check_df_corrCol()`:
      ! The condition_col parameter cannot be null. Select at least two columns from FEV1_FVC, Emphysema, Exacerbations, Age, Sex, Smoking, BMI, FEV1.

