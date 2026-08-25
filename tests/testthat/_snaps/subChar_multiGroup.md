# subChar_multiGroup can handle invalid inputs.

    Code
      subChar_multiGroup(processed_se, char = "Invalid", subChar = "class",
        ref_group = "ctrl", post_hoc = "One-way ANOVA", post_hoc_sig = "pval",
        post_hoc_p_cutoff = 0.05, transform = "log10")
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `subChar_multiGroup()`:
      ! The 'char' parameter must be a numeric lipid characteristic, such as 'Total.C'. Please execute list_lipid_char function to view the valid lipid characteristics you can input.

---

    Code
      subChar_multiGroup(processed_se, char = "Total.C", subChar = "Total.C",
        ref_group = "ctrl", post_hoc = "One-way ANOVA", post_hoc_sig = "pval",
        post_hoc_p_cutoff = 0.05, transform = "log10")
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `subChar_multiGroup()`:
      ! The 'char' and 'subChar' parameters cannot be the same.

---

    Code
      subChar_multiGroup(processed_se, char = "Total.C", subChar = "Invalid",
        ref_group = "ctrl", post_hoc = "One-way ANOVA", post_hoc_sig = "pval",
        post_hoc_p_cutoff = 0.05, transform = "log10")
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `subChar_multiGroup()`:
      ! The 'subChar' parameter is not in the list of lipid characteristics. Please execute list_lipid_char function to view the valid lipid characteristics you can input.

---

    Code
      subChar_multiGroup(processed_se, char = "Total.C", subChar = "class",
        ref_group = "Invalid", post_hoc = "One-way ANOVA", post_hoc_sig = "pval",
        post_hoc_p_cutoff = 0.05, transform = "log10")
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `subChar_multiGroup()`:
      ! ref_group must be one of the group names in the 'group' column of the group information table.

---

    Code
      subChar_multiGroup(processed_se, char = "Total.C", subChar = "class",
        ref_group = "ctrl", post_hoc = "Invalid", post_hoc_sig = "pval",
        post_hoc_p_cutoff = 0.05, transform = "log10")
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "One-way ANOVA", "Kruskal-Wallis test"

---

    Code
      subChar_multiGroup(processed_se, char = "Total.C", subChar = "class",
        ref_group = "ctrl", post_hoc = "One-way ANOVA", post_hoc_sig = "Pvalue",
        post_hoc_p_cutoff = 0.05, transform = "log10")
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "pval", "padj"

---

    Code
      subChar_multiGroup(processed_se, char = "Total.C", subChar = "class",
        ref_group = "ctrl", post_hoc = "One-way ANOVA", post_hoc_sig = "pval",
        post_hoc_p_cutoff = 2, transform = "log10")
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `subChar_multiGroup()`:
      ! The 'post_hoc_p_cutoff' parameter must be a numeric value between 0 and 1.

---

    Code
      subChar_multiGroup(processed_se, char = "Total.C", subChar = "class",
        ref_group = "ctrl", post_hoc = "One-way ANOVA", post_hoc_sig = "pval",
        post_hoc_p_cutoff = 0.05, transform = "log2")
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "log10", "cube", "square", "none"

---

    Code
      subChar_multiGroup(processed_se_two, char = "Total.C", subChar = "class",
        ref_group = "ctrl", post_hoc = "One-way ANOVA", post_hoc_sig = "pval",
        post_hoc_p_cutoff = 0.05, transform = "log10")
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `subChar_multiGroup()`:
      ! This function is used for multiple groups. If your data consists of two groups, please use subChar_twoGroup function for analysis.

