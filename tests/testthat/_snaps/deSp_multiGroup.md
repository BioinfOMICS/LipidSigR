# deSp_multiGroup handles invalid inputs correctly

    Code
      deSp_multiGroup(processed_se_multiGroup, ref_group = "abc", test = "One-way ANOVA",
        significant = "pval", p_cutoff = 0.05, transform = "log10")
    Condition
      Error in `deSp_multiGroup()`:
      ! ref_group must be one of the group names in the 'group' column of the group information table.

---

    Code
      deSp_multiGroup(processed_se_multiGroup, ref_group = "ctrl", test = "one-way ANOVA",
        significant = "pval", p_cutoff = 0.05, transform = "log10")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "One-way ANOVA", "Kruskal-Wallis test"

---

    Code
      deSp_multiGroup(processed_se_multiGroup, ref_group = "ctrl", test = "One-way ANOVA",
        significant = "p-value", p_cutoff = 0.05, transform = "log10")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "pval", "padj"

---

    Code
      deSp_multiGroup(processed_se_multiGroup, ref_group = "ctrl", test = "One-way ANOVA",
        significant = "pval", p_cutoff = 2, transform = "log10")
    Condition
      Error in `deSp_multiGroup()`:
      ! The 'p_cutoff' parameter must be a numeric value between 0 and 1.

---

    Code
      deSp_multiGroup(processed_se_multiGroup, ref_group = "ctrl", test = "One-way ANOVA",
        significant = "pval", p_cutoff = 0.05, transform = "log2")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "log10", "square", "cube", "none"

---

    Code
      deSp_multiGroup(processed_se_two, ref_group = "ctrl", test = "One-way ANOVA",
        significant = "pval", p_cutoff = 0.05, transform = "log10")
    Condition
      Error in `deSp_multiGroup()`:
      ! This function is used for multiple groups. If your data consists of two groups, please use deSp_twoGroup function for analysis.

---

    Code
      deSp_multiGroup(se_subset, ref_group = "ctrl", test = "One-way ANOVA",
        significant = "pval", p_cutoff = 0.05, transform = "log10")
    Condition
      Error in `.check_nGroup()`:
      ! In Group information table of multiple groups, each group must contains more than 2 samples.

# deSp_multiGroup handles missing values correctly

    Code
      deSp_multiGroup(processed_se = mod_se, ref_group = "ctrl", test = "One-way ANOVA",
        significant = "pval", p_cutoff = 0.05, transform = "log10")
    Condition
      Error in `.check_imputation()`:
      ! Detect 0 values or NAs in abundance data. Please perform data imputation by the data_process function first.

