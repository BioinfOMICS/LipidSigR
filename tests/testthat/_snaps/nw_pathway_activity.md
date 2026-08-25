# Test for incorrect input data.

    Code
      nw_pathway_activity(deSp_se_multi, organism = "human")
    Condition
      Error in `nw_pathway_activity()`:
      ! This function is used for two groups. If your data consists of multiple groups, please subset two of those groups for analysis.

---

    Code
      nw_pathway_activity(deSp_se_twoGroup, organism = "humna")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "human", "mouse"

---

    Code
      suppressWarnings(nw_pathway_activity(deSp_se, organism = "human"))
    Condition
      Error in `nw_pathway_activity()`:
      ! Comparing the means of two independent groups requires each group to have at least two samples.

