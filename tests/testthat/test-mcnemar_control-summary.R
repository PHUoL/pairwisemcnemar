# Helper: deterministic fake mcnemarControl object for summary tests
make_fake_mcnemar_fit <- function() {
  results <- tibble::tibble(
    control = rep("Control", 4),
    treatment = c("A", "B", "C", "D"),
    n = c(30, 30, 30, 30),
    discordant = c(12, 10, 14, 8),
    method_label = rep("asymptotic", 4),
    Z = c(2.30, 2.05, 1.90, 0.80),
    # raw p-value order: A < B < C < D
    p_value = c(0.020, 0.040, 0.060, 0.200),
    # adjusted p-value order: B < C < A < D
    p_adjusted = c(0.080, 0.040, 0.060, 0.200),
    p_adjust_label = rep("holm", 4)
  )
  
  settings <- list(
    id = "id",
    condition = "condition",
    outcome = "outcome",
    control = "Control",
    treatments = c("A", "B", "C", "D"),
    method = "asymptotic",
    p_adjust = "holm",
    ci = 0.95,
    outcome_levels = c("No", "Yes"),
    drop_na = TRUE,
    gamma = 1e-4,
    num_pi_values = 1000L
  )
  
  structure(
    list(results = results, settings = settings),
    class = "mcnemarControl"
  )
}

test_that("summary.mcnemarControl returns class 'summary.mcnemarControl'", {
  dat <- mcnemar_example_long(seed = 1)
  fit <- mcnemar_control(
    dat,
    id = "id",
    condition = "condition",
    outcome = "outcome",
    control = "Control"
  )
  
  s <- summary(fit)
  
  expect_s3_class(s, "summary.mcnemarControl")
})

test_that("summary.mcnemarControl contains overview counts", {
  dat <- mcnemar_example_long(seed = 1)
  fit <- mcnemar_control(
    dat,
    id = "id",
    condition = "condition",
    outcome = "outcome",
    control = "Control"
  )
  
  s <- summary(fit)
  
  expect_true(is.list(s$overview))
  
  expect_true(all(c(
    "control", "method", "p_adjust", "ci", "alpha",
    "comparisons", "valid_p", "zero_discordant",
    "significant_raw", "significant_adjusted"
  ) %in% names(s$overview)))
  
  expect_equal(s$overview$control, fit$settings$control)
  expect_equal(s$overview$method, fit$settings$method)
  expect_equal(s$overview$p_adjust, fit$settings$p_adjust)
  expect_equal(s$overview$comparisons, nrow(fit$results))
})

test_that("summary.mcnemarControl contains significant-treatment vectors", {
  dat <- mcnemar_example_long(seed = 1)
  fit <- mcnemar_control(
    dat,
    id = "id",
    condition = "condition",
    outcome = "outcome",
    control = "Control"
  )
  
  s <- summary(fit)
  
  expect_true(is.character(s$significant_raw))
  expect_true(is.character(s$significant_adjusted))
  
  if (length(s$significant_raw) > 0) {
    expect_true(all(s$significant_raw %in% fit$results$treatment))
  }
  if (length(s$significant_adjusted) > 0) {
    expect_true(all(s$significant_adjusted %in% fit$results$treatment))
  }
})

test_that("summary.mcnemarControl contains discordant diagnostics", {
  dat <- mcnemar_example_long(seed = 1)
  fit <- mcnemar_control(
    dat,
    id = "id",
    condition = "condition",
    outcome = "outcome",
    control = "Control"
  )
  
  s <- summary(fit)
  
  expect_true(is.list(s$discordant))
  expect_true(all(c(
    "min", "q1", "median", "mean", "q3", "max",
    "low_n", "threshold"
  ) %in% names(s$discordant)))
})

test_that("summary.mcnemarControl contains ordered comparisons table", {
  dat <- mcnemar_example_long(seed = 1)
  fit <- mcnemar_control(
    dat,
    id = "id",
    condition = "condition",
    outcome = "outcome",
    control = "Control"
  )
  
  s <- summary(fit)
  
  expect_true("ordered" %in% names(s))
  expect_s3_class(s$ordered, "data.frame")
  
  expect_true(all(c(
    "control", "treatment", "n", "discordant",
    "method_label", "Z", "p_value", "p_adjusted", "p_adjust_label"
  ) %in% names(s$ordered)))
  
  expect_equal(nrow(s$ordered), nrow(fit$results))
})

test_that("summary.mcnemarControl default ordered table is ascending by p_adjusted", {
  fit <- make_fake_mcnemar_fit()
  
  s <- summary(fit)
  
  # Lowest adjusted p-value first: B (0.04), C (0.06), A (0.08), D (0.20)
  expect_equal(as.character(s$ordered$treatment), c("B", "C", "A", "D"))
  
  # Verify nondecreasing p_adjusted
  expect_true(all(diff(s$ordered$p_adjusted) >= 0))
})

test_that("summary.mcnemarControl alpha changes number of significant treatments", {
  fit <- make_fake_mcnemar_fit()
  
  s_small <- summary(fit, alpha = 0.05)
  s_large <- summary(fit, alpha = 0.10)
  
  # At alpha=0.05, adjusted significant = B only
  expect_equal(sort(s_small$significant_adjusted), "B")
  expect_equal(s_small$overview$significant_adjusted, 1)
  
  # At alpha=0.10, adjusted significant = A, B, C
  expect_equal(sort(s_large$significant_adjusted), c("A", "B", "C"))
  expect_equal(s_large$overview$significant_adjusted, 3)
  
  expect_true(s_large$overview$significant_adjusted >= s_small$overview$significant_adjusted)
  expect_true(s_large$overview$significant_raw >= s_small$overview$significant_raw)
})

test_that("summary.mcnemarControl sort_by = 'p_value' changes ordering behavior", {
  fit <- make_fake_mcnemar_fit()
  
  s_adj <- summary(fit, sort_by = "p_adjusted")
  s_raw <- summary(fit, sort_by = "p_value")
  
  # By adjusted p-value ascending: B, C, A, D
  expect_equal(as.character(s_adj$ordered$treatment), c("B", "C", "A", "D"))
  
  # By raw p-value ascending: A, B, C, D
  expect_equal(as.character(s_raw$ordered$treatment), c("A", "B", "C", "D"))
  
  expect_false(identical(as.character(s_adj$ordered$treatment),
                         as.character(s_raw$ordered$treatment)))
})

test_that("summary.mcnemarControl contains notes and next_step", {
  dat <- mcnemar_example_long(seed = 1)
  fit <- mcnemar_control(
    dat,
    id = "id",
    condition = "condition",
    outcome = "outcome",
    control = "Control"
  )
  
  s <- summary(fit)
  
  expect_true(is.character(s$notes))
  expect_true(is.character(s$next_step))
  expect_equal(length(s$next_step), 1L)
})

test_that("summary.mcnemarControl works for exact methods with Z = NA", {
  testthat::skip_if_not_installed("contingencytables")
  
  dat <- mcnemar_example_long(seed = 1)
  fit <- mcnemar_control(
    dat,
    id = "id",
    condition = "condition",
    outcome = "outcome",
    control = "Control",
    method = "exact_cond"
  )
  
  s <- summary(fit)
  
  expect_s3_class(s, "summary.mcnemarControl")
  expect_true(any(grepl("Z is reported as NA", s$notes)))
})

test_that("print.summary.mcnemarControl prints ordered comparisons without error", {
  dat <- mcnemar_example_long(seed = 1)
  fit <- mcnemar_control(
    dat,
    id = "id",
    condition = "condition",
    outcome = "outcome",
    control = "Control"
  )
  
  s <- summary(fit)
  
  out <- capture.output(print(s))
  
  expect_true(any(grepl("^mcnemarControl summary$", out)))
  expect_true(any(grepl("^Overview$", out)))
  expect_true(any(grepl("^Significant treatments$", out)))
  expect_true(any(grepl("^Discordant diagnostics$", out)))
  expect_true(any(grepl("^Ordered comparisons$", out)))
  expect_true(any(grepl("^Next step$", out)))
})