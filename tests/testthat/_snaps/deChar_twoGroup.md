# deChar_twoGroup handles invalid inputs

    Code
      deChar_twoGroup(processed_se, char = "Invalid", ref_group = "ctrl", test = "t-test",
        significant = "padj", p_cutoff = 0.05, FC_cutoff = 1, transform = "log10")
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `deChar_twoGroup()`:
      ! The 'char' parameter is not in the list of lipid characteristics. Please execute list_lipid_char() to view the valid lipid characteristics you can input.

---

    Code
      deChar_twoGroup(processed_se, char = "Total.C", ref_group = "invalid", test = "t-test",
        significant = "padj", p_cutoff = 0.05, FC_cutoff = 1, transform = "log10")
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `deChar_twoGroup()`:
      ! ref_group must be one of the group names in the 'group' column of the group information table.

---

    Code
      deChar_twoGroup(processed_se, char = "Total.C", ref_group = "ctrl", test = "invalid",
        significant = "padj", p_cutoff = 0.05, FC_cutoff = 1, transform = "log10")
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "t-test", "Wilcoxon test", "limma"

---

    Code
      deChar_twoGroup(processed_se, char = "Total.C", ref_group = "ctrl", test = "t-test",
        significant = "invalid", p_cutoff = 0.05, FC_cutoff = 1, transform = "log10")
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "pval", "padj"

---

    Code
      deChar_twoGroup(processed_se, char = "Total.C", ref_group = "ctrl", test = "t-test",
        significant = "padj", p_cutoff = 2, FC_cutoff = 1, transform = "log10")
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `deChar_twoGroup()`:
      ! p_cutoff must be a numeric value between 0 and 1.

---

    Code
      deChar_twoGroup(processed_se, char = "Total.C", ref_group = "ctrl", test = "t-test",
        significant = "padj", p_cutoff = 0.05, FC_cutoff = 0)
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `deChar_twoGroup()`:
      ! FC_cutoff must be a numeric value >= 1.

---

    Code
      deChar_twoGroup(processed_se, char = "Total.C", ref_group = "ctrl", test = "t-test",
        significant = "padj", p_cutoff = 1, FC_cutoff = 1, transform = "log2")
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "log10", "square", "cube", "none"

---

    Code
      deChar_twoGroup(processed_se, char = "Chains odd/even ratio", ref_group = "ctrl",
        test = "t-test", significant = "padj", p_cutoff = 1, FC_cutoff = 1,
        transform = "log10")
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "log2", "none"

---

    Code
      deChar_twoGroup(processed_se_multi, char = "Total.C", ref_group = "ctrl", test = "t-test",
        significant = "padj", p_cutoff = 1, FC_cutoff = 1, transform = "log10")
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `deChar_twoGroup()`:
      ! This function is used for two groups. If your data consists of multi groups, please use deChar_multiGroup function for analysis.

# plot_deChar_twoGroup handles invalid inputs

    Code
      plot_deChar_twoGroup(deChar_se)
    Condition
      Error in `plot_deChar_twoGroup()`:
      ! This function is used for two groups. If your data consists of multi groups, please use plot_deChar_multiGroup function for analysis.

