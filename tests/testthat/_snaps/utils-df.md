# as_summarized_experiment can handle NULL or 0*0 data frame error input.

    Code
      as_summarized_experiment(abundance = NULL, two_char_table, group_info = NULL,
        se_type = "profiling", paired_sample = FALSE)
    Condition
      Error in `.check_df_null_0()`:
      ! Abundance data cannot be null.

---

    Code
      as_summarized_experiment(abundance = data.frame(), two_char_table, group_info = NULL,
      se_type = "profiling", paired_sample = FALSE)
    Condition
      Error in `.check_df_null_0()`:
      ! Abundance data must be a data frame and cannot be data frame with 0 columns and 0 rows.

---

    Code
      as_summarized_experiment(two_abundance, goslin_annotation = data.frame(),
      group_info = NULL, se_type = "de_two", paired_sample = FALSE)
    Condition
      Error in `.check_df_null_0()`:
      ! Goslin lipid characteristics table must be a data frame and cannot be data frame with 0 columns and 0 rows.

---

    Code
      as_summarized_experiment(two_abundance, two_char_table, group_info = data.frame(),
      se_type = "de_two", paired_sample = FALSE)
    Condition
      Error in `.check_valid_inputData()`:
      ! Group information table must be a data frame and cannot be data frame with 0 columns and 0 rows.

---

    Code
      as_summarized_experiment(two_abundance, two_char_table, group_info = "ABC",
        se_type = "de_two", paired_sample = FALSE)
    Condition
      Error in `.check_valid_inputData()`:
      ! Group information table must be a data frame and cannot be data frame with 0 columns and 0 rows.

# as_summarized_experiment can handle group_info error input.

    Code
      as_summarized_experiment(two_abundance, two_char_table, group_info_twoGroup,
        se_type = "invalid", paired_sample = FALSE)
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "profiling", "de_two", "de_multiple", "ml", "corr", "rna"

---

    Code
      as_summarized_experiment(two_abundance, two_char_table, group_info_twoGroup,
        se_type = "de_two", paired_sample = "invalid")
    Condition
      Error in `.check_valid_inputData()`:
      ! paired_sample must be logical value for two group data of differential analysis.

---

    Code
      as_summarized_experiment(two_abundance, two_char_table, group_info_twoGroup,
        se_type = "de_multiple", paired_sample = NULL)
    Condition
      Error in `.check_group_info()`:
      ! Group information table must contain columns of sample_name, label_name, and group.

---

    Code
      as_summarized_experiment(multi_abundance, multi_char_table,
        group_info_multiGroup, se_type = "de_multiple", paired_sample = TRUE)
    Condition
      Error in `.check_valid_inputData()`:
      ! paired_sample should be NULL.

# as_summarized_experiment can handle group_info format error.

    Code
      as_summarized_experiment(two_abundance, two_char_table, group_info_multiGroup,
        se_type = "de_two", paired_sample = FALSE)
    Condition
      Error in `.check_group_info()`:
      ! Group information table must contain columns of sample_name, label_name, group, and pair.

---

    Code
      as_summarized_experiment(multi_abundance, multi_char_table, group_info_twoGroup,
        se_type = "de_multiple", paired_sample = NULL)
    Condition
      Error in `.check_group_info()`:
      ! Group information table must contain columns of sample_name, label_name, and group.

---

    Code
      as_summarized_experiment(two_abundance, two_char_table, wrong_name, se_type = "de_two",
        paired_sample = FALSE)
    Condition
      Error in `.check_group_info()`:
      ! Group information table column names must be arranged in order of sample_name, label_name, group, and pair.

---

    Code
      as_summarized_experiment(multi_abundance, multi_char_table, wrong_name,
        se_type = "de_multiple", paired_sample = NULL)
    Condition
      Error in `.check_group_info()`:
      ! Group information table column names must be arranged in order of sample_name, label_name, and group.

---

    Code
      as_summarized_experiment(two_abundance[, 1:10], two_char_table,
      group_info_twoGroup, se_type = "de_two", paired_sample = FALSE)
    Condition
      Error in `.check_group_info()`:
      ! Sample names in group information table 'sample_name' column must be as same as the sample names in lipid abundance data.

---

    Code
      as_summarized_experiment(two_abundance, two_char_table, wrong_group_info,
        se_type = "de_two", paired_sample = FALSE)
    Condition
      Error in `.check_group_info()`:
      ! NA values are not allowed in the 'sample_name', 'label_name', and 'group' columns.

---

    Code
      as_summarized_experiment(two_abundance, two_char_table, wrong_twoGroup,
        se_type = "de_two", paired_sample = FALSE)
    Condition
      Error in `.check_group_info()`:
      ! The column 'group' in group information table must contain 2 groups.

---

    Code
      as_summarized_experiment(multi_abundance, multi_char_table, wrong_multiGroup,
        se_type = "de_multiple", paired_sample = NULL)
    Condition
      Error in `.check_group_info()`:
      ! The column 'group' in group information table must contain more than 2 groups.

---

    Code
      as_summarized_experiment(multi_oneSample, multi_char_table, wrong_multiGroup,
        se_type = "de_multiple", paired_sample = NULL)
    Condition
      Error in `.check_group_info()`:
      ! In Group information table of multiple groups, each group must contains more than 2 samples.

---

    Code
      as_summarized_experiment(two_abundance, two_char_table, group_info_twoGroup,
        se_type = "de_two", paired_sample = TRUE)
    Condition
      Error in `.check_pair_sample()`:
      ! Incorrect values in group information table 'pair' column. Each pair
              should be sequentially numbered from 1 to N without any missing, blank,
              or skipped numbers; otherwise, the value must be marked as NA.

---

    Code
      as_summarized_experiment(two_abundance[, 1:23], two_char_table,
      group_info_twoGroup, se_type = "de_two", paired_sample = TRUE)
    Condition
      Error in `.check_pair_sample()`:
      ! Incorrect values in group information table 'pair' column. Each pair
              should be sequentially numbered from 1 to N without any missing, blank,
              or skipped numbers; otherwise, the value must be marked as NA.

---

    Code
      as_summarized_experiment(two_abundance[, 1:23], two_char_table,
      group_info_twoGroup, se_type = "de_two", paired_sample = TRUE)
    Condition
      Error in `.check_pair_sample()`:
      ! Incorrect values in group information table 'pair' column. Each pair
              should be sequentially numbered from 1 to N without any missing, blank,
              or skipped numbers; otherwise, the value must be marked as NA.

# as_summarized_experiment can handle wrong correlation group_info.

    Code
      as_summarized_experiment(corr_abundance, corr_char_table, group_info, se_type = "corr",
        paired_sample = NULL)
    Condition
      Error in `.check_group_info()`:
      ! Group information table must include columns for sample_name and at least two clinical terms.

---

    Code
      as_summarized_experiment(corr_abundance, corr_char_table, group_info, se_type = "corr",
        paired_sample = NULL)
    Condition
      Error in `.check_group_info()`:
      ! Except for the first column, sample_name, the group information table must include at least two columns with numeric values and no missing values.

---

    Code
      as_summarized_experiment(corr_abundance, corr_char_table, group_info, se_type = "corr",
        paired_sample = NULL)
    Condition
      Error in `.check_group_info()`:
      ! Except for the first column, sample_name, the group information table must include at least two columns with numeric values and no missing values.

# as_summarized_experiment can handle wrong input.

    Code
      as_summarized_experiment(two_char_table, two_abundance, group_info_twoGroup,
        se_type = "de_two", paired_sample = FALSE)
    Condition
      Error in `.check_goslin_annotation()`:
      ! The `goslin_annotation` table needs to include the following columns: Normalized.Name, Original.Name, Grammar, Adduct, Lipid.Maps.Category, Lipid.Maps.Main.Class, Species.Name, Molecular.Species.Name, Sn.Position.Name, Level, Total.C, Total.OH, Total.DB, Mass, Sum.Formula, FA1.Position, FA1.C, FA1.OH, FA1.DB, FA1.Bond.Type, FA1.DB.Positions, FA2.Position, FA2.C, FA2.OH, FA2.DB, FA2.Bond.Type, FA2.DB.Positions, LCB.Position, LCB.C, LCB.OH, LCB.DB, LCB.Bond.Type, LCB.DB.Positions, FA3.Position, FA3.C, FA3.OH, FA3.DB, FA3.Bond.Type, FA3.DB.Positions, FA4.Position, FA4.C, FA4.OH, FA4.DB, FA4.Bond.Type, FA4.DB.Positions.

---

    Code
      as_summarized_experiment(group_info_twoGroup, two_char_table, two_abundance,
        se_type = "de_two", paired_sample = FALSE)
    Condition
      Error in `.check_group_info()`:
      ! Group information table must contain columns of sample_name, label_name, group, and pair.

# list_lipid_char can handle wrong input

    Code
      list_lipid_char(data.frame())
    Condition
      Error in `.check_inputSE()`:
      ! The input must be a SummarizedExperiment object construct by as_summarized_experiment function or output from upstream analysis function.

---

    Code
      list_lipid_char(abundance_twoGroup)
    Condition
      Error in `.check_inputSE()`:
      ! The input must be a SummarizedExperiment object construct by as_summarized_experiment function or output from upstream analysis function.

---

    Code
      list_lipid_char(abundance_twoGroup$feature)
    Condition
      Error in `.check_inputSE()`:
      ! The input must be a SummarizedExperiment object construct by as_summarized_experiment function or output from upstream analysis function.

# extract_summarized_experiment can output correct result list

    Code
      extract_summarized_experiment(data.frame())
    Condition
      Error in `extract_summarized_experiment()`:
      ! The input must be a SummarizedExperiment object construct
                   by as_summarized_experiment function or ouput from upstream analysis function.

---

    Code
      extract_summarized_experiment(abundance_twoGroup)
    Condition
      Error in `extract_summarized_experiment()`:
      ! The input must be a SummarizedExperiment object construct
                   by as_summarized_experiment function or ouput from upstream analysis function.

---

    Code
      extract_summarized_experiment(abundance_twoGroup$feature)
    Condition
      Error in `extract_summarized_experiment()`:
      ! The input must be a SummarizedExperiment object construct
                   by as_summarized_experiment function or ouput from upstream analysis function.

