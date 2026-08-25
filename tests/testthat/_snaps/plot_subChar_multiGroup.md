# plot_subChar_multiGroup handles errors correctly

    Code
      plot_subChar_multiGroup(subChar_se, subChar_feature = "lipid")
    Condition
      Error in `plot_subChar_multiGroup()`:
      ! The 'subChar_feature' parameter must be one of the valid subChar features.

---

    Code
      plot_subChar_multiGroup(invalid_se, subChar_feature = "CL")
    Condition
      Error in `.check_inputSE()`:
      ! Correct SummarizedExperiment metadata not found, please use output from upstream analysis function.

