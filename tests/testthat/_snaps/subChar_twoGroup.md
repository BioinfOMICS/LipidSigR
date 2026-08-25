# subChar_twoGroup handles errors correctly

    Code
      subChar_twoGroup(processed_se, char = "Invalid", subChar = "class", ref_group = "ctrl",
        test = "t-test", significant = "pval", p_cutoff = 0.05, FC_cutoff = 1,
        transform = "log10")
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `subChar_twoGroup()`:
      ! The 'char' parameter must be a numeric lipid characteristic, such as 'Total.C'. Please execute list_lipid_char function to view the valid lipid characteristics you can input.

---

    Code
      subChar_twoGroup(processed_se, char = "Total.C", subChar = "Invalid",
        ref_group = "ctrl", test = "t-test", significant = "pval", p_cutoff = 0.05,
        FC_cutoff = 1, transform = "log10")
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `subChar_twoGroup()`:
      ! The 'subChar' parameter is not in the list of lipid characteristics. Please execute list_lipid_char function to view the valid lipid characteristics you can input.

---

    Code
      subChar_twoGroup(processed_se, char = "Total.C", subChar = "Total.C",
        ref_group = "ctrl", test = "t-test", significant = "pval", p_cutoff = 0.05,
        FC_cutoff = 1, transform = "log10")
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `subChar_twoGroup()`:
      ! The 'char' and 'subChar' parameters cannot be the same.

---

    Code
      subChar_twoGroup(processed_se, char = "Total.C", subChar = "class", ref_group = "control",
        test = "t-test", significant = "pval", p_cutoff = 0.05, FC_cutoff = 1,
        transform = "log10")
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `subChar_twoGroup()`:
      ! ref_group must be one of the group names in the 'group' column of the group information table.

---

    Code
      subChar_twoGroup(processed_se, char = "Total.C", subChar = "class", ref_group = "ctrl",
        test = "ttest", significant = "pval", p_cutoff = 0.05, FC_cutoff = 1,
        transform = "log10")
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "t-test", "Wilcoxon test"

---

    Code
      subChar_twoGroup(processed_se, char = "Total.C", subChar = "class", ref_group = "ctrl",
        test = "t-test", significant = "pvalue", p_cutoff = 0.05, FC_cutoff = 1,
        transform = "log10")
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "pval", "padj"

---

    Code
      subChar_twoGroup(processed_se, char = "Total.C", subChar = "class", ref_group = "ctrl",
        test = "t-test", significant = "pval", p_cutoff = 2, FC_cutoff = 1,
        transform = "log10")
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `subChar_twoGroup()`:
      ! The 'p_cutoff' parameter must be a numeric value between 0 and 1.

---

    Code
      subChar_twoGroup(processed_se, char = "Total.C", subChar = "class", ref_group = "ctrl",
        test = "t-test", significant = "pval", p_cutoff = 0.05, FC_cutoff = 1,
        transform = "log2")
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "log10", "cube", "square", "none"

---

    Code
      subChar_twoGroup(processed_se_multi, char = "Total.C", subChar = "class",
        ref_group = "ctrl", test = "t-test", significant = "pval", p_cutoff = 0.05,
        FC_cutoff = 1, transform = "log10")
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `subChar_twoGroup()`:
      ! This function is used for two groups. If your data consists of multi groups, please use subChar_multiGroup function for analysis.

# plot_subChar_twoGroup handles errors correctly

    Code
      plot_subChar_twoGroup(subChar_se, subChar_feature = "Invalid")
    Condition
      Error in `plot_subChar_twoGroup()`:
      ! The 'subChar_feature' parameter must be one of the valid subChar features.

