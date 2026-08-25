# heatmap_correlation handles all char.

    Code
      heatmap_correlation(processed_se, char = char_list[i], transform = "log10",
      correlation = "pearson", distfun = "maximum", hclustfun = "average", type = "class")
    Condition
      Error in `heatmap_correlation()`:
      ! Not enough data for plotting the heatmap.

---

    Code
      heatmap_correlation(processed_se, char = char_list[i], transform = "log10",
      correlation = "pearson", distfun = "maximum", hclustfun = "average", type = "class")
    Condition
      Error in `heatmap_correlation()`:
      ! Not enough data for plotting the heatmap.

# heatmap_correlation handles invalid parameters correctly

    Code
      heatmap_correlation(processed_se, char = NULL, transform = "log10",
        correlation = "pearson", distfun = "maximum", hclustfun = "average", type = "class")
    Condition
      Error in `heatmap_correlation()`:
      ! For 'class' type, char cannot be NULL.

---

    Code
      heatmap_correlation(processed_se, char = 12345, transform = "log10",
        correlation = "pearson", distfun = "maximum", hclustfun = "average", type = "class")
    Condition
      Error in `heatmap_correlation()`:
      ! Wrong char input, you can view the available char list by list_lipid_char function.

---

    Code
      heatmap_correlation(processed_se, char = NULL, transform = "invalid",
        correlation = "pearson", distfun = "maximum", hclustfun = "average", type = "sample")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "log10", "cube", "square", "none"

---

    Code
      heatmap_correlation(processed_se, char = NULL, transform = "log10",
        correlation = "ABC", distfun = "maximum", hclustfun = "average", type = "sample")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "pearson", "spearman"

---

    Code
      heatmap_correlation(processed_se, char = NULL, transform = "log10",
        correlation = "pearson", distfun = "ward.D2", hclustfun = "average", type = "sample")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "maximum", "pearson", "spearman", "kendall", "euclidean", "manhattan", "canberra", "binary", "minkowski"

---

    Code
      heatmap_correlation(processed_se, char = NULL, transform = "log10",
        correlation = "pearson", distfun = "maximum", hclustfun = "class", type = "sample")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "average", "ward.D", "ward.D2", "single", "complete", "mcquitty", "median", "centroid"

---

    Code
      heatmap_correlation(processed_se, char = NULL, transform = "log10",
        correlation = "pearson", distfun = "maximum", hclustfun = "average", type = "lipid")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "sample", "class"

# heatmap_correlation handles unprocessed data correctly

    Code
      heatmap_correlation(profiling_data, char = NULL, transform = "log10",
        correlation = "pearson", distfun = "maximum", hclustfun = "average", type = "sample")
    Condition
      Error in `.check_imputation()`:
      ! Detect 0 values or NAs in abundance data. Please perform data imputation by the data_process function first.

