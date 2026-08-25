# deSp_twoGroup handles errors correctly

    Code
      deSp_twoGroup(processed_se, ref_group = "nonexistent")
    Condition
      Error in `deSp_twoGroup()`:
      ! ref_group must be one of the group names in the 'group' column of the group information table.

---

    Code
      deSp_twoGroup(processed_se, ref_group = "ctrl", test = "invalid_test")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "t-test", "Wilcoxon test", "limma"

---

    Code
      deSp_twoGroup(processed_se, ref_group = "ctrl", significant = "invalid")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "pval", "padj"

---

    Code
      deSp_twoGroup(processed_se, ref_group = "ctrl", p_cutoff = 2)
    Condition
      Error in `deSp_twoGroup()`:
      ! p_cutoff must be a numeric value between 0 and 1.

---

    Code
      deSp_twoGroup(processed_se, ref_group = "ctrl", FC_cutoff = 0)
    Condition
      Error in `deSp_twoGroup()`:
      ! FC_cutoff must be a numeric value >= 1.

---

    Code
      deSp_twoGroup(processed_se, ref_group = "ctrl", transform = "invalid")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "log10", "square", "cube", "none"

---

    Code
      deSp_twoGroup(processed_se_multi, ref_group = "ctrl", test = "t-test",
        significant = "padj", p_cutoff = 0.05, FC_cutoff = 1, transform = "log10")
    Condition
      Error in `deSp_twoGroup()`:
      ! This function is used for two groups. If your data consists of multi groups, please use deSp_multiGroup function for analysis.

# plot_deSp_twoGroup handles errors correctly

    Code
      plot_deSp_twoGroup("not a SummarizedExperiment object")
    Condition
      Error in `.check_de_outputSE()`:
      ! The input must be a SummarizedExperiment object construct by upstream analysis function, e.g. deSp_twoGroup, deChar_twoGroup...

