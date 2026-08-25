# Test under normal conditions: multiple groups

    Code
      enrichment_lsea(deSp_se = deSp_se_multiGroup, char = NULL, rank_by = "log2FC",
        significant = "pval", p_cutoff, n_lipid = 2)
    Condition
      Error in `enrichment_lsea()`:
      ! Log2 fold change (log2FC) ranking is not suitable for multiple group data. Please choose p-value (pval), adjusted p-value (padj), or statistic as the ranking criterion instead.

# Test for incorrect input & error message

    Code
      enrichment_lsea(deSp_se_twoGroup, char = "abc", rank_by = "statistic",
        significant = "pval", p_cutoff = 0.05)
    Condition
      Error in `enrichment_lsea()`:
      ! The 'char' parameter must be one of the listed lipid characteristics. Please use the list_lipid_char function to view which lipid characteristics can be analyzed.

---

    Code
      enrichment_lsea(deSp_se_twoGroup, char = NULL, rank_by = "Statistic",
        significant = "pval", p_cutoff = 0.05)
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "statistic", "log2FC", "pval", "padj"

---

    Code
      enrichment_lsea(deSp_se_twoGroup, char = NULL, rank_by = "statistic",
        significant = "Pval", p_cutoff = 0.05)
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "pval", "padj"

---

    Code
      enrichment_lsea(deSp_se_twoGroup, char = NULL, rank_by = "statistic",
        significant = "pval", p_cutoff = 1.5)
    Condition
      Error in `enrichment_lsea()`:
      ! The 'p_cutoff' parameter must be a numeric value between 0 and 1.

