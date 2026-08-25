# Test for incorrect input and error message

    Code
      nw_lipid_reaction(deSp_se = deSp_se_twoGroup, organism = "humna", show_sp = "sigClass",
        show_all_reactions = FALSE, sp_significant = "pval", sp_p_cutoff = 0.05,
        sp_FC_cutoff = 1, class_significant = "pval", class_p_cutoff = 0.05,
        class_FC_cutoff = 1)
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "human", "mouse"

---

    Code
      nw_lipid_reaction(deSp_se = deSp_se_twoGroup, organism = "human", show_sp = "invalid",
        show_all_reactions = FALSE, sp_significant = "pval", sp_p_cutoff = 0.05,
        sp_FC_cutoff = 1, class_significant = "pval", class_p_cutoff = 0.05,
        class_FC_cutoff = 1)
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "sigClass", "all", "none"

---

    Code
      nw_lipid_reaction(deSp_se = deSp_se_twoGroup, organism = "human", show_sp = "sigClass",
        show_all_reactions = FALSE, sp_significant = "Pval", sp_p_cutoff = 0.05,
        sp_FC_cutoff = 1, class_significant = "pval", class_p_cutoff = 0.05,
        class_FC_cutoff = 1)
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "pval", "padj"

---

    Code
      nw_lipid_reaction(deSp_se = deSp_se_twoGroup, organism = "human", show_sp = "sigClass",
        show_all_reactions = FALSE, sp_significant = "pval", sp_p_cutoff = 0.05,
        sp_FC_cutoff = 1, class_significant = "Pval", class_p_cutoff = 0.05,
        class_FC_cutoff = 1)
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "pval", "padj"

---

    Code
      nw_lipid_reaction(deSp_se = deSp_se_twoGroup, organism = "human", show_sp = "sigClass",
        show_all_reactions = FALSE, sp_significant = "pval", sp_p_cutoff = 1.1,
        sp_FC_cutoff = 2, class_significant = "pval", class_p_cutoff = 0.05,
        class_FC_cutoff = 1)
    Condition
      Error in `nw_lipid_reaction()`:
      ! The 'sp_p_cutoff' parameter must be a numeric value between 0 and 1.

---

    Code
      nw_lipid_reaction(deSp_se = deSp_se_twoGroup, organism = "human", show_sp = "sigClass",
        show_all_reactions = FALSE, sp_significant = "pval", sp_p_cutoff = 0.05,
        sp_FC_cutoff = 1, class_significant = "pval", class_p_cutoff = 1.5,
        class_FC_cutoff = 1)
    Condition
      Error in `nw_lipid_reaction()`:
      ! The 'class_p_cutoff' parameter must be a numeric value between 0 and 1.

---

    Code
      nw_lipid_reaction(deSp_se = deSp_se_twoGroup, organism = "human", show_sp = "sigClass",
        show_all_reactions = FALSE, sp_significant = "pval", sp_p_cutoff = 0.05,
        sp_FC_cutoff = 0, class_significant = "pval", class_p_cutoff = 0.05,
        class_FC_cutoff = 1)
    Condition
      Error in `nw_lipid_reaction()`:
      ! The minimum value for the 'sp_FC_cutoff' parameter is 1.

---

    Code
      nw_lipid_reaction(deSp_se = deSp_se_twoGroup, organism = "human", show_sp = "sigClass",
        show_all_reactions = FALSE, sp_significant = "pval", sp_p_cutoff = 0.05,
        sp_FC_cutoff = 1, class_significant = "pval", class_p_cutoff = 0.05,
        class_FC_cutoff = 0)
    Condition
      Error in `nw_lipid_reaction()`:
      ! The minimum value for the 'class_FC_cutoff' parameter is 1.

---

    Code
      nw_lipid_reaction(deSp_se_multi, organism = "human", show_sp = "sigClass",
        show_all_reactions = FALSE, sp_significant = "pval", sp_p_cutoff = 0.05,
        sp_FC_cutoff = 1, class_significant = "pval", class_p_cutoff = 0.05,
        class_FC_cutoff = 1)
    Condition
      Error in `nw_lipid_reaction()`:
      ! This function is used for two groups. If your data consists of multiple groups, please subset two of those groups for analysis.

---

    Code
      nw_lipid_reaction(deSp_se, organism = "human", show_sp = "sigClass",
        show_all_reactions = FALSE, sp_significant = "pval", sp_p_cutoff = 0.05,
        sp_FC_cutoff = 1, class_significant = "pval", class_p_cutoff = 0.05,
        class_FC_cutoff = 1)
    Condition
      Error in `nw_lipid_reaction()`:
      ! This function requires each group to have at least two samples.

