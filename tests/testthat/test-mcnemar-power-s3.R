test_that("mcnemar_power() returns expected class and core columns", {
  dat <- mcnemar_example_long(seed = 1)

  fit <- mcnemar_power(
    data = dat,
    id = "id",
    condition = "condition",
    outcome = "outcome",
    control = "Control",
    methods = c("asymptotic", "midp", "exact_cond"),
    p_adjust = "holm"
  )

  expect_s3_class(fit, "mcnemarPower")
  expect_true(is.list(fit))
  expect_true("results" %in% names(fit))
  expect_true("settings" %in% names(fit))

  expect_true(is.data.frame(fit$results))
  expect_true(all(c(
    "control", "treatment", "method_label",
    "a", "b", "c", "d", "n",
    "pb_hat", "pc_hat", "PowerS", "N80S", "N90S"
  ) %in% names(fit$results)))
})

test_that("print and summary methods dispatch for mcnemarPower", {
  dat <- mcnemar_example_long(seed = 1)

  fit <- mcnemar_power(
    data = dat,
    id = "id",
    condition = "condition",
    outcome = "outcome",
    control = "Control",
    methods = c("asymptotic"),
    p_adjust = "holm"
  )

  out <- capture.output(ret <- print(fit))
  expect_s3_class(ret, "mcnemarPower")
  expect_true(length(out) > 0)

  s <- summary(
    fit,
    data = dat,
    comparisons = c("TrtA", "TrtB"),
    scenarios = c("indep"),
    nsim = 50,
    seed = 123
  )

  expect_s3_class(s, "summary.mcnemarPower")
  expect_true("single_test" %in% names(s))
  expect_true("adjusted_power" %in% names(s))
  expect_true("settings" %in% names(s))
  out <- capture.output(ret <- print(s))
  expect_s3_class(ret, "summary.mcnemarPower")
  expect_true(length(out) > 0)
})

test_that("mcnemar_power_summary() mirrors summary()", {
  dat <- mcnemar_example_long(seed = 1)

  fit <- mcnemar_power(
    data = dat,
    id = "id",
    condition = "condition",
    outcome = "outcome",
    control = "Control",
    methods = c("asymptotic"),
    p_adjust = "holm"
  )

  s1 <- summary(
    fit,
    data = dat,
    comparisons = c("TrtA"),
    scenarios = c("indep"),
    nsim = 25,
    seed = 123
  )

  s2 <- mcnemar_power_summary(
    fit,
    data = dat,
    comparisons = c("TrtA"),
    scenarios = c("indep"),
    nsim = 25,
    seed = 123
  )

  expect_s3_class(s2, "summary.mcnemarPower")
  expect_equal(class(s1), class(s2))
  expect_equal(names(s1), names(s2))
})
