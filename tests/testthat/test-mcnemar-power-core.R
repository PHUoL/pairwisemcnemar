test_that("mcnemar_power() returns expected class and core structure", {
  # Generate a deterministic toy dataset from the package example generator.
  dat <- mcnemar_example_long(seed = 1)
  
  # Fit a simple mcnemarPower object using three methods so the object
  # structure is exercised across a small method set.
  fit <- mcnemar_power(
    data = dat,
    id = "id",
    condition = "condition",
    outcome = "outcome",
    control = "Control",
    methods = c("asymptotic", "midp", "exact_cond"),
    p_adjust = "holm"
  )
  
  # Check the top-level S3 class and object structure.
  expect_s3_class(fit, "mcnemarPower")
  expect_true(is.list(fit))
  expect_true("results" %in% names(fit))
  expect_true("settings" %in% names(fit))
  
  # Check that the main result component is tabular and settings is a list.
  expect_true(is.data.frame(fit$results))
  expect_true(is.list(fit$settings))
})

test_that("mcnemar_power() results tibble contains required columns", {
  # Use a simple asymptotic-only fit for a lightweight structural test.
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
  
  # Check that the current implementation exposes the documented columns.
  expect_true(all(c(
    "control", "treatment", "method_label",
    "a", "b", "c", "d", "n",
    "pb_hat", "pc_hat", "PowerS", "N80S", "N90S"
  ) %in% names(fit$results)))
})

test_that("mcnemar_power() numerical outputs satisfy basic invariants", {
  # Use a modest n_max to keep the test quick while still exercising N80S/N90S.
  dat <- mcnemar_example_long(seed = 1)
  
  fit <- mcnemar_power(
    data = dat,
    id = "id",
    condition = "condition",
    outcome = "outcome",
    control = "Control",
    methods = c("asymptotic"),
    p_adjust = "holm",
    n_max = 120
  )
  
  res <- fit$results
  
  # Power is a probability, so it must lie in [0, 1].
  expect_true(all(res$PowerS >= 0 & res$PowerS <= 1))
  
  # Estimated discordant probabilities must also lie in [0, 1].
  expect_true(all(res$pb_hat >= 0 & res$pb_hat <= 1))
  expect_true(all(res$pc_hat >= 0 & res$pc_hat <= 1))
  
  # The total discordant probability cannot exceed 1.
  expect_true(all(res$pb_hat + res$pc_hat <= 1 + 1e-12))
  
  # If both sample-size targets are finite, 90% power cannot be achieved
  # at a smaller n than 80% power.
  ok <- is.finite(res$N80S) & is.finite(res$N90S)
  if (any(ok)) {
    expect_true(all(res$N90S[ok] >= res$N80S[ok]))
  }
})

test_that("mcnemar_power() respects treatment subset selection", {
  # Restricting 'treatments' should restrict the rows in the output.
  dat <- mcnemar_example_long(seed = 1)
  
  fit <- mcnemar_power(
    data = dat,
    id = "id",
    condition = "condition",
    outcome = "outcome",
    control = "Control",
    treatments = c("TrtA", "TrtC"),
    methods = c("asymptotic"),
    p_adjust = "holm"
  )
  
  expect_equal(sort(unique(fit$results$treatment)), sort(c("TrtA", "TrtC")))
})

test_that("print.mcnemarPower() runs invisibly", {
  # Printing should dispatch and return invisibly without error.
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
})

test_that("mcnemar_power() validates key arguments", {
  # Use one deterministic dataset and intentionally trigger several
  # documented input-validation errors.
  dat <- mcnemar_example_long(seed = 1)
  
  # Bad column name for id.
  expect_error(
    mcnemar_power(
      data = dat,
      id = "bad_id",
      condition = "condition",
      outcome = "outcome",
      control = "Control"
    )
  )
  
  # Invalid method name.
  expect_error(
    mcnemar_power(
      data = dat,
      id = "id",
      condition = "condition",
      outcome = "outcome",
      control = "Control",
      methods = "not_a_method"
    )
  )
  
  # Invalid n_max.
  expect_error(
    mcnemar_power(
      data = dat,
      id = "id",
      condition = "condition",
      outcome = "outcome",
      control = "Control",
      n_max = 1
    )
  )
  
  # Invalid target power outside (0, 1).
  expect_error(
    mcnemar_power(
      data = dat,
      id = "id",
      condition = "condition",
      outcome = "outcome",
      control = "Control",
      targets = c(0.8, 1.2)
    )
  )
  
  # Invalid outcome_levels length.
  expect_error(
    mcnemar_power(
      data = dat,
      id = "id",
      condition = "condition",
      outcome = "outcome",
      control = "Control",
      outcome_levels = c("No", "Yes", "Maybe")
    )
  )
})