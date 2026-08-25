# Test for incorrect input and error message

    Code
      deChar_se <- deChar_multiGroup(processed_se, char = "lipid", ref_group = "ctrl",
        post_hoc = "One-way ANOVA", post_hoc_sig = "pval", post_hoc_p_cutoff = 0.05,
        transform = "log10")
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `deChar_multiGroup()`:
      ! The 'char' parameter is not in the list of lipid characteristics. Please execute list_lipid_char function to view the valid lipid characteristics you can input.

---

    Code
      deChar_se <- deChar_multiGroup(processed_se, char = "class", ref_group = "ctrl",
        post_hoc = "one-way ANOVA", post_hoc_sig = "pval", post_hoc_p_cutoff = 0.05,
        transform = "log10")
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "One-way ANOVA", "Kruskal-Wallis test"

---

    Code
      deChar_se <- deChar_multiGroup(processed_se, char = "class", ref_group = "ctrl",
        post_hoc = "One-way ANOVA", post_hoc_sig = "p-value", post_hoc_p_cutoff = 0.05,
        transform = "log10")
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "pval", "padj"

---

    Code
      deChar_se <- deChar_multiGroup(processed_se, char = "class", ref_group = "ctrl",
        post_hoc = "One-way ANOVA", post_hoc_sig = "pval", post_hoc_p_cutoff = 2,
        transform = "log10")
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `deChar_multiGroup()`:
      ! The 'post_hoc_p_cutoff' parameter must be a numeric value between 0 and 1.

---

    Code
      deChar_se <- deChar_multiGroup(processed_se, char = "class", ref_group = "ctrl",
        post_hoc = "One-way ANOVA", post_hoc_sig = "pval", post_hoc_p_cutoff = 0.05,
        transform = "log2")
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `deChar_multiGroup()`:
      ! The 'char' is not a ratio-based characteristic. The 'transform' parameter must be one of 'none', 'log10', 'square', or 'cube'.

---

    Code
      deChar_se <- deChar_multiGroup(processed_se, char = "Chains Ether/Ester linked ratio",
        ref_group = "ctrl", post_hoc = "One-way ANOVA", post_hoc_sig = "pval",
        post_hoc_p_cutoff = 0.05, transform = "log10")
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `deChar_multiGroup()`:
      ! The 'char' is a ratio-based characteristic. The 'transform' parameter must be either 'none' or 'log2'.

---

    Code
      deChar_se <- deChar_multiGroup(processed_se_two, char = "class", ref_group = "ctrl",
        post_hoc = "One-way ANOVA", post_hoc_sig = "pval", post_hoc_p_cutoff = 0.05,
        transform = "log10")
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `deChar_multiGroup()`:
      ! This function is used for multiple groups. If your data consists of two groups, please use deChar_twoGroup function for analysis.

---

    Code
      deChar_se <- deChar_multiGroup(processed_se, char = "class", ref_group = "123",
        post_hoc = "One-way ANOVA", post_hoc_sig = "pval", post_hoc_p_cutoff = 0.05,
        transform = "log10")
    Message
      There are 4 ratio characteristics that can be converted in your dataset.
    Condition
      Error in `deChar_multiGroup()`:
      ! ref_group must be one of the group names in the 'group' column of the group information table.

