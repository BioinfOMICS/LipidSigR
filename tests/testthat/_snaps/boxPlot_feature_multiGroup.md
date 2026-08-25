# Test that message return for invalid input

    Code
      boxPlot_feature_multiGroup(processed_se = processed_se_multiGroup, feature = "ABC",
        ref_group = "ctrl", test = "One-way ANOVA", post_hoc_sig = "padj", transform = "log10")
    Condition
      Error in `boxPlot_feature_multiGroup()`:
      ! 'feature' must be one of the lipids in your lipid abundance data; please revise your input.

---

    Code
      boxPlot_feature_multiGroup(processed_se = processed_se_multiGroup, feature = "PE O- 16:1;0_16:1;0",
        ref_group = "AA1", test = "One-way ANOVA", post_hoc_sig = "padj", transform = "log10")
    Condition
      Error in `boxPlot_feature_multiGroup()`:
      ! ref_group must be one of the group names in the 'group' column of the group information table.

---

    Code
      boxPlot_feature_multiGroup(processed_se = processed_se_multiGroup, feature = "PE O- 16:1;0_16:1;0",
        ref_group = "ctrl", test = "One-way-ANOVA", post_hoc_sig = "padj", transform = "log10")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "One-way ANOVA", "Kruskal-Wallis test"

---

    Code
      boxPlot_feature_multiGroup(processed_se = processed_se_multiGroup, feature = "PE O- 16:1;0_16:1;0",
        ref_group = "ctrl", test = "Kruskal-Wallis test", post_hoc_sig = "p-val",
        transform = "log10")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "pval", "padj"

---

    Code
      boxPlot_feature_multiGroup(processed_se = processed_se_multiGroup, feature = "PE O- 16:1;0_16:1;0",
        ref_group = "ctrl", test = "One-way ANOVA", post_hoc_sig = "padj", transform = "log2")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "log10", "square", "cube", "none"

---

    Code
      boxPlot_feature_multiGroup(processed_se = processed_se_two, feature = "PE O- 16:1;0-18:2;0",
        ref_group = "ctrl", test = "One-way ANOVA", post_hoc_sig = "padj", transform = "log10")
    Condition
      Error in `boxPlot_feature_multiGroup()`:
      ! This function is used for multiple groups. If your data consists of two groups, please use boxPlot_feature_twoGroup() for analysis.

