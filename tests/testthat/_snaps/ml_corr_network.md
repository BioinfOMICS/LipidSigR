# ml_corr_network can handle input errors.

    Code
      ml_corr_network(ml_se, feature_importance = "correlation", correlation = "pearson",
        edge_cutoff = 0, feature_num = 10, nsim = 5)
    Condition
      Error in `ml_corr_network()`:
      ! feature_importance must be "Algorithm-based" or "SHAP".

---

    Code
      ml_corr_network(ml_se, feature_importance = "Algorithm-based", correlation = NULL,
        edge_cutoff = 0, feature_num = 10, nsim = 5)
    Condition
      Error in `ml_corr_network()`:
      ! correlation must be one of "pearson", "kendall", or "spearman".

---

    Code
      ml_corr_network(ml_se, feature_importance = "Algorithm-based", correlation = "pearson",
        edge_cutoff = 0, feature_num = "10", nsim = 5)
    Condition
      Error in `ml_corr_network()`:
      ! feature_num must be a numeric value chosen from the following options: 2, 3, 5, 10, 20, 50, 88 .

---

    Code
      ml_corr_network(ml_se, feature_importance = "Algorithm-based", correlation = "pearson",
        edge_cutoff = 0, feature_num = 10, nsim = -1)
    Condition
      Error in `ml_corr_network()`:
      ! nsim must be a integer >= 0.

---

    Code
      ml_corr_network(ml_se, feature_importance = "Algorithm-based", correlation = "pearson",
        edge_cutoff = "0", feature_num = 10, nsim = 5)
    Condition
      Error in `ml_corr_network()`:
      ! edge_cutoff must be a numeric value between 0 and 1.

