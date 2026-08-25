# Test for incorrect input & error messages

    Code
      plot_enrichment_lsea(error_input, char = "class", char_feature = "PC")
    Condition
      Error in `plot_enrichment_lsea()`:
      ! Incorrect LSEA result object. Please use enrichment_lsea function to obtain the correct LSEA result object.

---

    Code
      plot_enrichment_lsea(deSp_se_twoGroup, char = "class", char_feature = "PC")
    Condition
      Error in `plot_enrichment_lsea()`:
      ! Incorrect LSEA result object. Please use enrichment_lsea function to obtain the correct LSEA result object.

---

    Code
      plot_enrichment_lsea(lsea_res, char = "lipid", char_feature = "PC")
    Condition
      Error in `plot_enrichment_lsea()`:
      ! The 'char' and 'char_feature' parameters must be included in the 'enrich_result' returned by enrichment_lsea function.

---

    Code
      plot_enrichment_lsea(lsea_res, char = "class", char_feature = "CP")
    Condition
      Error in `plot_enrichment_lsea()`:
      ! The 'char' and 'char_feature' parameters must be included in the 'enrich_result' returned by enrichment_lsea function.

