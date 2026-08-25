# plot_ml_evaluation handles wrong input

    Code
      plot_ml_evaluation(ml_se, eval_method = "invalid")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "Accuracy", "Sensitivity", "Specificity", "Pos Pred Value", "Neg Pred Value", "Precision", "Recall", "F1", "Prevalence", "Detection Rate", "Detection Prevalence", "Balanced Accuracy"

