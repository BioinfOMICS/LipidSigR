# Test for incorrect input and related error messages

    Code
      enrichment_ora(deSp_se_twoGroup, char = "lipid", significant = "pval",
        p_cutoff = 0.05)
    Condition
      Error in `enrichment_ora()`:
      ! The 'char' parameter must be one of the listed lipid characteristics. Please use the list_lipid_char function to view which lipid characteristics can be analyzed.

---

    Code
      enrichment_ora(deSp_se_twoGroup, char = "class", significant = "Pval",
        p_cutoff = 0.05)
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "pval", "padj"

---

    Code
      enrichment_ora(deSp_se_twoGroup, char = "class", significant = "pval",
        p_cutoff = 1.5)
    Condition
      Error in `enrichment_ora()`:
      ! The 'p_cutoff' parameter must be a numeric value between 0 and 1.

