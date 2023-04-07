# Load data and create example data sets
data(profiling_exp_data)

# Start tests
test_that("exp_profiling function works.", {
  expect_s3_class(profiling_exp_data, "data.frame")
  expect_error(exp_profiling(profiling_exp_data[, 1], plotly=TRUE))
  result_plotly <- exp_profiling(profiling_exp_data, plotly=TRUE)
  expect_s3_class(result_plotly$i.expr.lip, "plotly")
  expect_s3_class(result_plotly$i.p.amount, "plotly")
  expect_s3_class(result_plotly$p.hist.value, "plotly")
  result_ggplot <- exp_profiling(profiling_exp_data, plotly=FALSE)
  expect_s3_class(result_ggplot$i.expr.lip, "ggplot")
  expect_s3_class(result_ggplot$i.p.amount, "ggplot")
  expect_s3_class(result_ggplot$p.hist.value, "ggplot")
})
