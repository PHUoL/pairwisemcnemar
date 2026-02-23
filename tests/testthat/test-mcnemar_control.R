test_that("mcnemar_control returns comparisons", {
  dat <- mcnemar_example_long(seed = 1)
  fit <- mcnemar_control(dat, "id", "condition", "outcome", control = "Control")
  expect_true(nrow(fit$results) >= 1)
})
