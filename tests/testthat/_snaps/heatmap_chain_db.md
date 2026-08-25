# Test error conditions

    Code
      heatmap_chain_db(processed_se_twoGroup, char = "invalid_char", ref_group = "ctrl")
    Condition
      Error in `heatmap_chain_db()`:
      ! The 'char' parameter is not in the list of lipidcharacteristics. Please execute list_lipid_char() to viewthe valid lipid characteristics you can input.

---

    Code
      heatmap_chain_db(processed_se_twoGroup, char = "class", char_feature = "invalid_feature",
        ref_group = "ctrl")
    Condition
      Error in `heatmap_chain_db()`:
      ! The 'char_feature' parameter must be one of the valid char features.

---

    Code
      heatmap_chain_db(processed_se_twoGroup, char = "class", ref_group = "invalid_group")
    Condition
      Error in `heatmap_chain_db()`:
      ! ref_group must be one of the group names in the 'group' column of the group information table.

---

    Code
      heatmap_chain_db(processed_se_twoGroup, char = "class", ref_group = "ctrl",
        test = "One-way ANOVA")
    Condition
      Error in `heatmap_chain_db()`:
      ! 'One-way ANOVA' and 'Kruskal-Wallis test' are for multiple group comparisons. Please choose either 't-test' or 'Wilcoxon test' for two group analysis.

---

    Code
      heatmap_chain_db(processed_se_multiGroup, char = "class", ref_group = "ctrl",
        test = "t-test")
    Condition
      Error in `heatmap_chain_db()`:
      ! 't-test' and 'Wilcoxon test' are for two group comparisons. Please choose either 'One-way ANOVA' or 'Kruskal-Wallis test' for multiple group analysis.

---

    Code
      heatmap_chain_db(processed_se_twoGroup, char = "class", ref_group = "ctrl",
        significant = "invalid")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "pval", "padj"

---

    Code
      heatmap_chain_db(processed_se_twoGroup, char = "class", ref_group = "ctrl",
        p_cutoff = 2)
    Condition
      Error in `heatmap_chain_db()`:
      ! The 'p_cutoff' parameter must be a numeric value between 0 and 1.

---

    Code
      heatmap_chain_db(processed_se_twoGroup, char = "class", ref_group = "ctrl",
        FC_cutoff = 0)
    Condition
      Error in `heatmap_chain_db()`:
      ! The minimum value for the 'FC_cutoff' parameter is 1.

---

    Code
      heatmap_chain_db(processed_se_twoGroup, char = "class", ref_group = "ctrl",
        transform = "invalid")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "log10", "square", "cube", "none"

# Test edge cases

    Code
      heatmap_chain_db(empty_se, char = "class", ref_group = "ctrl")
    Condition
      Error in `.check_inputSE()`:
      ! The input must be a SummarizedExperiment object construct by as_summarized_experiment function or output from upstream analysis function.

---

    Code
      heatmap_chain_db(single_feature_se, char = "class", ref_group = "ctrl")
    Condition
      Error in `.check_inputSE()`:
      ! Input SummarizedExperiment contains fewer than two features. At least two rows (features) are required for this analysis. Please check your object or use `as_summarized_experiment()` to reconstruct it.

---

    Code
      heatmap_chain_db(processed_se_twoGroup, char = "class", test = "One-way ANOVA",
        ref_group = "ctrl", p_cutoff = 1e-10)
    Condition
      Error in `heatmap_chain_db()`:
      ! 'One-way ANOVA' and 'Kruskal-Wallis test' are for multiple group comparisons. Please choose either 't-test' or 'Wilcoxon test' for two group analysis.

