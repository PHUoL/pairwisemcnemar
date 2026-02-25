test_that("mcnemar_control returns comparisons", {
  dat <- mcnemar_example_long(seed = 1)
  fit <- mcnemar_control(dat, "id", "condition", "outcome", control = "Control")
  expect_true(nrow(fit$results) >= 1)
})

test_that("midp produces comparisons and p-values", {
  dat <- mcnemar_example_long(seed = 1)
  fit <- mcnemar_control(dat, "id", "condition", "outcome",
                         control = "Control", method = "midp")
  expect_gt(nrow(fit$results), 0)
  expect_true(all(is.na(fit$results$Z)))
  expect_true(all(!is.na(fit$results$p_value)))
})