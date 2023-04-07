# Load data and create example data sets
data(DE_exp_data)
data(DE_lipid_char_table)
data(DE_group_info)

# Start tests
test_that("DE_sub_char_2 functions works", {
  exp_transform_table <- data_process(DE_exp_data, exclude_var_missing=TRUE,
     missing_pct_limit=50, replace_zero=TRUE, zero2what='min', xmin=0.5,
     replace_NA=TRUE, NA2what='min', ymin=0.5, pct_transform=TRUE,
     data_transform=TRUE, trans_type='log', centering=FALSE, scaling=FALSE)
  lipid_char_filter <-
    DE_lipid_char_table[DE_lipid_char_table$feature %in%
                          exp_transform_table$feature, ]
  char_var <- colnames(lipid_char_filter)[-1]
  DE.sub.char.2 <- DE_sub_char_2(DE_exp_data, data_transform=TRUE,
      lipid_char_table=DE_lipid_char_table, split_var = char_var[1],
      char_var = char_var[4], group_info = DE_group_info, paired=FALSE,
      sig_pvalue=0.05, sig_FC=2, exclude_var_missing=TRUE, missing_pct_limit=50,
      replace_zero=TRUE, zero2what='min', xmin=0.5, replace_NA=TRUE,
      NA2what='min', ymin=0.5, pct_transform=TRUE, trans_type='log',
      centering=FALSE, scaling=FALSE)
  expect_s3_class(DE.sub.char.2[[1]], "data.frame")
  expect_s3_class(DE.sub.char.2[[2]], "data.frame")
  expect_s3_class(DE.sub.char.2[[3]], "data.frame")
  expect_s3_class(DE.sub.char.2[[4]], "data.frame")
  char.class <- unique(DE.sub.char.2[[2]][1])
  sub_char_result <- DE_sub_char_plot_2(DE.sub.char.2[[2]], DE.sub.char.2[[3]],
      group_info=DE_group_info, char_var=char_var[4], split_var=char_var[1],
      split_class=char.class[5,], insert_ref_group=NULL, ref_group=NULL,
      plotly=TRUE)
  expect_s3_class(sub_char_result[[1]], "plotly")
  expect_s3_class(sub_char_result[[2]], "plotly")
  expect_s3_class(sub_char_result[[3]], "plotly")
  expect_s3_class(sub_char_result[[4]], "plotly")
  expect_s3_class(sub_char_result[[5]], "plotly")
  sub_char_ggplot <- DE_sub_char_plot_2(DE.sub.char.2[[2]], DE.sub.char.2[[3]],
      group_info=DE_group_info, char_var=char_var[4], split_var=char_var[1],
      split_class=char.class[5,], insert_ref_group=NULL, ref_group=NULL,
      plotly=FALSE)
  expect_s3_class(sub_char_ggplot[[1]], "ggplot")
  expect_s3_class(sub_char_ggplot[[2]], "ggplot")
  expect_s3_class(sub_char_ggplot[[3]], "ggplot")
  expect_s3_class(sub_char_ggplot[[4]], "ggplot")
  expect_s3_class(sub_char_ggplot[[5]], "ggplot")
})
