# Load data and create example data sets
data(profiling_exp_data)
data(profiling_lipid_char_table)
test_lipid_char <- data.frame(feature = profiling_lipid_char_table$feature,
                              class = profiling_lipid_char_table$totallength,
                              totallength = profiling_lipid_char_table$class,
                              totaldb = profiling_lipid_char_table$class,
                              totaloh = profiling_lipid_char_table$class)

# Start tests
test_that("exp_compo_by_lipidinfo works", {
  char_var <- colnames(profiling_lipid_char_table)[-1]
  expect_type(char_var, "character")
  # test At least 2 samples.
  expect_error(exp_compo_by_lipidinfo(profiling_exp_data[, 1:2],
     profiling_lipid_char_table, char_var[1], plotly=TRUE))
  # test lipid_char_table content of column 'class' must be characters
  expect_error(exp_compo_by_lipidinfo(profiling_exp_data,
      test_lipid_char, char_var[1], plotly=TRUE))
  # test lipid_char_table content of column 'totallength' must be numeric
  test_lipid_char$class <- profiling_lipid_char_table$class
  expect_error(exp_compo_by_lipidinfo(profiling_exp_data,
      test_lipid_char, char_var[1], plotly=TRUE))
  # test lipid_char_table content of column 'totaldb' must be numeric
  test_lipid_char$totallength <- profiling_lipid_char_table$totallength
  expect_error(exp_compo_by_lipidinfo(profiling_exp_data,
      test_lipid_char, char_var[1], plotly=TRUE))
  # test lipid_char_table content of column 'totaloh' must be numeric
  test_lipid_char$totaldb <- profiling_lipid_char_table$totaldb
  expect_error(exp_compo_by_lipidinfo(profiling_exp_data,
      test_lipid_char, char_var[1], plotly=TRUE))

  compo_result <- exp_compo_by_lipidinfo(profiling_exp_data,
                                         profiling_lipid_char_table,
                                         char_var[1],
                                         plotly=TRUE)
  expect_s3_class(compo_result$p.barplot.p, "plotly")
  expect_s3_class(compo_result$p.compos, "plotly")
  compo_ggplot <- exp_compo_by_lipidinfo(profiling_exp_data,
                                         profiling_lipid_char_table,
                                         char_var[1],
                                         plotly=FALSE)
  expect_s3_class(compo_ggplot$p.barplot.p, "ggplot")
  expect_s3_class(compo_ggplot$p.compos, "ggplot")
})
