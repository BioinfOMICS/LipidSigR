# Test for incorrect input

    Code
      suppressWarnings(nw_gatom(deSp_se = deSp_se_twoGroup, organism = "invalid",
        n_lipid = 50, sp_significant = "pval", sp_p_cutoff = 0.05, sp_FC_cutoff = 1))
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "human", "mouse"

---

    Code
      suppressWarnings(nw_gatom(deSp_se = deSp_se_twoGroup, organism = "mouse",
        n_lipid = 0, sp_significant = "pval", sp_p_cutoff = 0.05, sp_FC_cutoff = 1))
    Condition
      Error in `nw_gatom()`:
      ! The minimum value for the 'n_lipid' parameter is 1.

---

    Code
      suppressWarnings(nw_gatom(deSp_se = deSp_se_twoGroup, organism = "mouse",
        n_lipid = 50, sp_significant = "Pval", sp_p_cutoff = 0.05, sp_FC_cutoff = 1))
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "pval", "padj"

---

    Code
      suppressWarnings(nw_gatom(deSp_se = deSp_se_twoGroup, organism = "mouse",
        n_lipid = 50, sp_significant = "pval", sp_p_cutoff = 1.5, sp_FC_cutoff = 1))
    Condition
      Error in `nw_gatom()`:
      ! The 'sp_p_cutoff' parameter must be a numeric value between 0 and 1.

---

    Code
      suppressWarnings(nw_gatom(deSp_se = deSp_se_twoGroup, organism = "mouse",
        n_lipid = 50, sp_significant = "pval", sp_p_cutoff = 0.05, sp_FC_cutoff = 0))
    Condition
      Error in `nw_gatom()`:
      ! The minimum value for the 'sp_FC_cutoff' parameter is 1.

---

    Code
      suppressWarnings(nw_gatom(deSp_se = deSp_se_multi, organism = "mouse", n_lipid = 50,
        sp_significant = "pval", sp_p_cutoff = 0.05, sp_FC_cutoff = 1))
    Condition
      Error in `nw_gatom()`:
      ! This function is used for two groups. If your data consists of multiple groups, please subset two of those groups for analysis.

---

    Code
      suppressWarnings(nw_gatom(deSp_se, organism = "mouse", n_lipid = 50,
        sp_significant = "pval", sp_p_cutoff = 0.05, sp_FC_cutoff = 1))
    Condition
      Error in `nw_gatom()`:
      ! This function requires each group to have at least two samples.

