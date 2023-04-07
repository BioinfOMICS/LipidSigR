# Load data and create example data sets
data(DE_exp_data)
data(DE_lipid_char_table)
test_lipid_char <- data.frame(feature = DE_lipid_char_table$feature,
                              class = DE_lipid_char_table$totallength,
                              totallength = DE_lipid_char_table$class,
                              totaldb = DE_lipid_char_table$class,
                              totaloh = DE_lipid_char_table$class)

# Start tests
test_that("Species2Char function works", {
  char_var <- colnames(DE_lipid_char_table)[-1]
  expect_type(char_var, "character")
  # test The lipids names of lipid_char_table table must same as exp_data
  expect_error(Species2Char(DE_exp_data, DE_lipid_char_table[, -1],
                            char_var = char_var[4]))
  # test The row number of lipid_char_table table must same as exp_data
  expect_error(Species2Char(DE_exp_data, DE_lipid_char_table[1:3, ],
                            char_var = char_var[4]))
  # test lipid_char_table content of column 'class' must be characters
  expect_error(Species2Char(DE_exp_data, test_lipid_char,
                            char_var = char_var[4]))
  test_lipid_char$class <- DE_lipid_char_table$class
  # test lipid_char_table content of column 'totallength' must be numeric
  expect_error(Species2Char(DE_exp_data, test_lipid_char,
                            char_var = char_var[4]))
  test_lipid_char$totallength <- DE_lipid_char_table$totallength
  # test lipid_char_table content of column 'totaldb' must be numeric
  expect_error(Species2Char(DE_exp_data, test_lipid_char,
                            char_var = char_var[4]))
  test_lipid_char$totaldb <- DE_lipid_char_table$totaldb
  # test lipid_char_table content of column 'totaloh' must be numeric
  expect_error(Species2Char(DE_exp_data, test_lipid_char,
                            char_var = char_var[4]))

  Species2Char(DE_exp_data, DE_lipid_char_table, char_var = char_var[8])

  exp_data_Spe2Char <- Species2Char(DE_exp_data, DE_lipid_char_table,
                                    char_var = char_var[4])
  expect_s3_class(exp_data_Spe2Char, "data.frame")
  expect_equal(ncol(exp_data_Spe2Char), ncol(DE_exp_data))
})
