# Load data and create example data sets
data(profiling_lipid_char_table)

# Start tests
test_that("lipid_char_table_gather function works.", {
  char_var <- colnames(profiling_lipid_char_table)[-1]
  expect_error(lipid_char_table_gather(profiling_lipid_char_table,
                                       char_var = "char_var"))
  result_table <- lipid_char_table_gather(profiling_lipid_char_table,
                                          char_var = char_var[8])
  expect_s3_class(result_table, "data.frame")
})
