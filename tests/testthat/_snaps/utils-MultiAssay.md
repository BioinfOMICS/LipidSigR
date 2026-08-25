# as_MultiAssayExperiment errors when there is no sample overlap

    Code
      as_MultiAssayExperiment(lipidSE, rnaSE_noOverlap)
    Condition
      Error in `as_MultiAssayExperiment()`:
      ! No common sample names found between lipidSE and rnaSE. Please check the colnames of both objects.

# as_MultiAssayExperiment errors when lipidSE has no group_info-derived colData (se_type = 'profiling')

    Code
      as_MultiAssayExperiment(lipidSE_profiling, rnaSE)
    Condition
      Error in `.check_mae_colData()`:
      ! lipidSE has no group_info-derived colData columns - it looks like it was built with se_type='profiling'. as_MultiAssayExperiment() requires lipidSE to be built with se_type='de_two' or se_type='de_multiple'.

# as_MultiAssayExperiment errors when rnaSE has no group_info-derived colData

    Code
      as_MultiAssayExperiment(lipidSE, rnaSE_bad)
    Condition
      Error in `.check_mae_colData()`:
      ! rnaSE has no group_info-derived colData columns. as_MultiAssayExperiment() requires rnaSE to be built with se_type='rna' and a de_two or de_multiple-style group_info.

# as_MultiAssayExperiment errors when lipidSE colData has no 'sample_name' column

    Code
      as_MultiAssayExperiment(lipidSE_noSampleName, rnaSE)
    Condition
      Error in `.check_mae_colData()`:
      ! colData of lipidSE must contain a 'sample_name' column.

# as_MultiAssayExperiment errors when rnaSE colData has no 'sample_name' column

    Code
      as_MultiAssayExperiment(lipidSE, rnaSE_noSampleName)
    Condition
      Error in `.check_mae_colData()`:
      ! colData of rnaSE must contain a 'sample_name' column.

# as_MultiAssayExperiment errors when lipidSE and rnaSE group_info structures differ

    Code
      as_MultiAssayExperiment(lipidSE, rnaSE_deMultiple)
    Condition
      Error in `.check_mae_colData()`:
      ! colData column structure differs between lipidSE and rnaSE (lipidSE: sample_name, label_name, group, pair; rnaSE: sample_name, label_name, group). This usually means the two objects were constructed with incompatible se_type values. Please ensure both SummarizedExperiment objects are constructed with a compatible se_type so their colData structures match.

# as_MultiAssayExperiment errors when common-sample group_info values differ between lipidSE and rnaSE

    Code
      as_MultiAssayExperiment(lipidSE, rnaSE_flipped)
    Condition
      Error in `.check_mae_colData()`:
      ! colData values differ between lipidSE and rnaSE for common samples. Please ensure both SummarizedExperiment objects were constructed from the same group_info table.

# as_MultiAssayExperiment can handle wrong input

    Code
      as_MultiAssayExperiment(data.frame(), rnaSE)
    Condition
      Error in `.check_inputSE()`:
      ! The input must be a SummarizedExperiment object construct by as_summarized_experiment function or output from upstream analysis function.

---

    Code
      as_MultiAssayExperiment(lipidSE, data.frame())
    Condition
      Error in `.check_inputSE()`:
      ! The input must be a SummarizedExperiment object construct by as_summarized_experiment function or output from upstream analysis function.

---

    Code
      as_MultiAssayExperiment(lipidSE, matrix(1:4, 2, 2))
    Condition
      Error in `.check_inputSE()`:
      ! The input must be a SummarizedExperiment object construct by as_summarized_experiment function or output from upstream analysis function.

---

    Code
      as_MultiAssayExperiment(NULL, rnaSE)
    Condition
      Error in `.check_inputSE()`:
      ! The input must be a SummarizedExperiment object construct by as_summarized_experiment function or output from upstream analysis function.

---

    Code
      as_MultiAssayExperiment(lipidSE, NULL)
    Condition
      Error in `.check_inputSE()`:
      ! The input must be a SummarizedExperiment object construct by as_summarized_experiment function or output from upstream analysis function.

