# ml_shap & plot_ml_shap works correctly and can handle error inputs

    Code
      ml_shap(ml_se, feature_num = "200", nsim = 5)
    Condition
      Error in `ml_shap()`:
      ! feature_num must be a numeric value chosen from the following options: 2, 3, 5, 10, 20, 50, 88 .

---

    Code
      ml_shap(ml_se, feature_num = 10, nsim = -1)
    Condition
      Error in `ml_shap()`:
      ! nsim must be a positive integer.

# plot_shap_sample can return correct plots.

    Code
      plot_shap_sample(shap_se, sample_id = "10")
    Condition
      Error in `plot_shap_sample()`:
      ! sample_id must be one of the ID number in 'ID' column of shap_result.

# plot_shap_force can return correct plots.

    Code
      plot_shap_force(shap_se, top_feature = 10, cluster_method = "cluster",
        group_num = 10)
    Condition
      Error in `plot_shap_force()`:
      ! cluster_method must be one of the strings "ward.D", "ward.D2",
               "single", "complete", "average", "mcquitty", "median", "centroid".

---

    Code
      plot_shap_force(shap_se, top_feature = 10, cluster_method = "ward.D",
        group_num = "10")
    Condition
      Error in `plot_shap_force()`:
      ! group_num must be a numeric value between 0 and 228 1.

---

    Code
      plot_shap_force(shap_se, top_feature = "10", cluster_method = "ward.D",
        group_num = 10)
    Condition
      Error in `plot_shap_force()`:
      ! top_feature must be a positive value.

# plot_shap_dependence can return correct plots.

    Code
      plot_shap_dependence(shap_se, feature = "feature1", shap_feature = selected_feature[
        2], interaction_index = selected_feature[2])
    Condition
      Error in `plot_shap_dependence()`:
      ! feature must be one of the feature names in the 'variable' column of shap_result.

---

    Code
      plot_shap_dependence(shap_se, feature = selected_feature[1], shap_feature = NULL,
      interaction_index = selected_feature[2])
    Condition
      Error in `plot_shap_dependence()`:
      ! shap_feature must be one of the feature names in the 'variable' column of shap_result.

---

    Code
      plot_shap_dependence(shap_se, feature = selected_feature[1], shap_feature = selected_feature[
        2], interaction_index = "abd")
    Condition
      Error in `plot_shap_dependence()`:
      ! interaction_index must be one of the feature names in the 'variable' column of shap_result.

