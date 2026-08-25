# Input parameter validation works correctly

    Code
      ml_model(processed_se, char = "char", transform = "log10", ranking_method = "Random_forest",
        ml_method = "Random_forest", split_prop = 0.7, nfold = 3, alpha = NULL)
    Condition
      Error in `ml_model()`:
      ! Wrong char input, you can view the available char list by list_lipid_char function.

---

    Code
      ml_model(processed_se, char = NULL, transform = "log10", ranking_method = "Random_forest",
        ml_method = "Random_forest", split_prop = 0.7, nfold = 3, alpha = NULL)
    Condition
      Error in `ml_model()`:
      ! Wrong char input, you can view the available char list by list_lipid_char function.

---

    Code
      ml_model(processed_se, char = "none", transform = "log20", ranking_method = "Random_forest",
        ml_method = "Random_forest", split_prop = 0.7, nfold = 3, alpha = NULL)
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "log10", "cube", "square", "none"

---

    Code
      ml_model(processed_se, char = "none", transform = "log10", ranking_method = "Random_forest",
        ml_method = "Random_forest", split_prop = 0, nfold = 3, alpha = NULL)
    Condition
      Error in `ml_model()`:
      ! split_prop must be a numeric value between 0.5 and 1.

---

    Code
      ml_model(processed_se, char = "none", transform = "log10", ranking_method = "Random_forest",
        ml_method = "Random_forest", split_prop = 0.7, nfold = 0, alpha = NULL)
    Condition
      Error in `ml_model()`:
      ! nfold must be a numeric value >= 1.

---

    Code
      ml_model(processed_se, char = "none", transform = "log10", ranking_method = "Random_forest",
        ml_method = "ElasticNet", split_prop = 0.7, nfold = 3, alpha = -1)
    Condition
      Error in `ml_model()`:
      ! alpha must be a numeric value between 0 and 1.

---

    Code
      ml_model(processed_se, char = "none", transform = "log10", ranking_method = "Random_forest",
        ml_method = "ElasticNet", split_prop = 0.7, nfold = 3, alpha = 2)
    Condition
      Error in `ml_model()`:
      ! alpha must be a numeric value between 0 and 1.

# Ranking methods & machine learning methods work correctly for all conditions

    Code
      ml_model(processed_se, char = "none", transform = "log10", ranking_method = "invalid",
        ml_method = "Random_forest", split_prop = 0.7, nfold = 3, alpha = NULL)
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "Random_forest", "p_value", "pvalue_FC", "ROC", "SVM", "Lasso", "Ridge", "ElasticNet"

---

    Code
      ml_model(processed_se, char = "none", transform = "log10", ranking_method = "ElasticNet",
        ml_method = "Random_forest", split_prop = 0.7, nfold = 3, alpha = "1")
    Condition
      Error in `ml_model()`:
      ! alpha must be a numeric value between 0 and 1.

---

    Code
      ml_model(processed_se, char = "none", transform = "log10", ranking_method = "Random_forest",
        ml_method = "invalid", split_prop = 0.7, nfold = 3, alpha = NULL)
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "Random_forest", "SVM", "Lasso", "Ridge", "ElasticNet", "xgboost"

# Data transformation options work correctly

    Code
      ml_model(processed_se, char = "none", transform = "invalid", ranking_method = "Random_forest",
        ml_method = "Random_forest", split_prop = 0.7, nfold = 3, alpha = NULL)
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "log10", "cube", "square", "none"

