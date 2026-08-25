# data_process function work correctly

    Code
      data_process(se, exclude_missing = FALSE, exclude_missing_pct = 70,
        replace_na_method = "none", replace_na_method_ref = 0.5, normalization = "none")
    Condition
      Error in `data_process()`:
      ! Detect 0 values or NAs in abundance data. Please select a data imputation method in replace_na_method other than 'none'.

# data_process handles errors correctly

    Code
      data_process(se, exclude_missing = "TRUE", exclude_missing_pct = 70,
        replace_na_method = "min", replace_na_method_ref = 0.5, normalization = "Percentage",
        transform = "log10")
    Condition
      Error in `data_process()`:
      ! exclude_missing must be a logical value.

---

    Code
      data_process(se, exclude_missing = TRUE, exclude_missing_pct = 101,
        replace_na_method = "min", replace_na_method_ref = 0.5, normalization = "Percentage",
        transform = "log10")
    Condition
      Error in `data_process()`:
      ! exclude_missing_pct must be a numeric value between 5 and 100.

---

    Code
      data_process(se, exclude_missing = TRUE, exclude_missing_pct = 70,
        replace_na_method = "invalid_method", replace_na_method_ref = 0.5,
        normalization = "Percentage", transform = "log10")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "min", "BPCA", "GMS", "KNN", "impSeq", "Impseqrob", "LLS", "mean", "median", "MinProb", "PPCA", "PI", "QRILC", "RandomForest", "SVD", "none"

---

    Code
      data_process(se, exclude_missing = TRUE, exclude_missing_pct = 70,
        replace_na_method = "min", replace_na_method_ref = "invalid_value",
        normalization = "Percentage", transform = "log10")
    Condition
      Error in `data_process()`:
      ! replace_na_method_ref must be a numeric value.

---

    Code
      data_process(se, exclude_missing = TRUE, exclude_missing_pct = 70,
        replace_na_method = "min", replace_na_method_ref = 0.5, normalization = "invalid_method")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "Percentage", "CubicSplines", "CyclicLoess", "EigenMS", "linear", "Mean", "Median", "nonLinear", "PQN", "Quantile", "RLR", "Sum", "TMM", "none"

---

    Code
      data_process(se, exclude_missing = TRUE, exclude_missing_pct = 70,
        replace_na_method = "min", replace_na_method_ref = 0.5, normalization = "Percentage",
        transform = "invalid_method")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "log10", "log2", "cube", "square", "none"

---

    Code
      data_process(se, exclude_missing = TRUE, exclude_missing_pct = 70,
        replace_na_method = "BPCA", replace_na_method_ref = 11, normalization = "Percentage",
        transform = "log10")
    Condition
      Error in `.missing_impute()`:
      ! replace_na_method_ref must be a numeric value between 1 and 10.

---

    Code
      data_process(se, exclude_missing = TRUE, exclude_missing_pct = 70,
        replace_na_method = "KNN", replace_na_method_ref = 11, normalization = "Percentage",
        transform = "log10")
    Condition
      Error in `.missing_impute()`:
      ! replace_na_method_ref must be a numeric value between 1 and 10.

---

    Code
      data_process(se, exclude_missing = TRUE, exclude_missing_pct = 70,
        replace_na_method = "LLS", replace_na_method_ref = 27, normalization = "Percentage",
        transform = "log10")
    Message
      The replace_na_method_ref must be less than the sample size, so it is replaced by the sample size minus one.
    Condition
      Warning in `.missing_impute()`:
      183 negative values were detected after imputation with the LLS method; these values were reverted to half-minute imputation.
    Output
      class: SummarizedExperiment 
      dim: 182 23 
      metadata(2): processed_abund transform
      assays(1): abundance
      rownames(182): Cer 38:1;2 Cer 40:1;2 ... TAG 58:7;0 TAG 58:8;0
      rowData names(72): feature class ... MetaNetX.ID PlantFA.ID
      colnames(23): control_01 control_02 ... hfref_patient_12
        hfref_patient_13
      colData names(4): sample_name label_name group pair

---

    Code
      data_process(se, exclude_missing = TRUE, exclude_missing_pct = 70,
        replace_na_method = "min", replace_na_method_ref = 1.1, normalization = "Percentage",
        transform = "log10")
    Condition
      Error in `.missing_impute()`:
      ! replace_na_method_ref must be a numeric value between 0.1 and 0.5.

---

    Code
      data_process(se, exclude_missing = TRUE, exclude_missing_pct = 70,
        replace_na_method = "PPCA", replace_na_method_ref = 11, normalization = "Percentage",
        transform = "log10")
    Condition
      Error in `.missing_impute()`:
      ! replace_na_method_ref must be a numeric value between 1 and 10.

---

    Code
      data_process(se, exclude_missing = TRUE, exclude_missing_pct = 70,
        replace_na_method = "QRILC", replace_na_method_ref = 1.1, normalization = "Percentage",
        transform = "log10")
    Condition
      Error in `.missing_impute()`:
      ! replace_na_method_ref must be a numeric value between 0.1 and 1.

---

    Code
      data_process(se, exclude_missing = TRUE, exclude_missing_pct = 70,
        replace_na_method = "SVD", replace_na_method_ref = 11, normalization = "Percentage",
        transform = "log10")
    Condition
      Error in `.missing_impute()`:
      ! replace_na_method_ref must be a numeric value between 1 and 10.

