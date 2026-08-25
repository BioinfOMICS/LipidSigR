# heatmap_clustering handles different char parameters

    Code
      heatmap_clustering(deSp_se, char = NULL, distfun = "pearson", hclustfun = "complete",
        type = "all")
    Condition
      Error in `heatmap_clustering()`:
      ! Wrong char input, you can view the available char list by list_lipid_char function.

# heatmap_clustering handles different distfun parameters

    Code
      heatmap_clustering(deSp_se, char = "class", distfun = "invalid", hclustfun = "complete",
        type = "all")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "pearson", "kendall", "spearman"

# heatmap_clustering handles different hclustfun parameters

    Code
      heatmap_clustering(deSp_se, char = "class", distfun = "pearson", hclustfun = "invalid",
        type = "all")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "complete", "ward.D", "ward.D2", "single", "average", "mcquitty", "median", "centroid"

# heatmap_clustering handles different type parameters

    Code
      heatmap_clustering(deSp_se, char = "class", distfun = "pearson", hclustfun = "complete",
        type = "invalid")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "sig", "all"

# heatmap_clustering handles edge cases

    Code
      suppressWarnings(heatmap_clustering(sub, char = "class", distfun = "pearson",
        hclustfun = "complete", type = "sig"))
    Condition
      Error in `heatmap_clustering()`:
      ! Insufficient number of lipids.

---

    Code
      suppressWarnings(heatmap_clustering(no_sig_se, char = "class", distfun = "pearson",
        hclustfun = "complete", type = "sig"))
    Condition
      Error in `heatmap_clustering()`:
      ! Insufficient number of significant lipids.

# heatmap_clustering can handle one row data condition

    Code
      suppressWarnings(heatmap_clustering(sub, char = "class", distfun = "pearson",
        hclustfun = "complete", type = "sig"))
    Condition
      Error in `heatmap_clustering()`:
      ! Insufficient number of lipids.

