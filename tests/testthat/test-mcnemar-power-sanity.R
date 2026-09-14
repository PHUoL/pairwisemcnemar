test_that("PowerS exact_cond matches exact2x2 on mcnemar_example_long(n = 40, seed = 2)", {
  skip_if_not_installed("exact2x2")
  
  dat <- mcnemar_example_long(n = 40, seed = 2)
  
  fit <- mcnemar_power(
    data = dat,
    id = "id",
    condition = "condition",
    outcome = "outcome",
    control = "Control",
    methods = c("exact_cond"),
    p_adjust = "holm"
  )
  
  res <- fit$results
  expect_true(nrow(res) >= 1)
  
  for (i in seq_len(nrow(res))) {
    expect_equal(res$n[i], 40)
    
    ref_power <- exact2x2::powerPaired2x2(
      pb = res$pb_hat[i],
      pc = res$pc_hat[i],
      npairs = res$n[i],
      alternative = "two.sided"
    )$power
    
    # Use a slightly looser tolerance to allow tiny floating-point differences
    # between two valid implementations of the same exact conditional power.
    expect_equal(res$PowerS[i], ref_power, tolerance = 1e-7)
    
    expect_true(is.finite(res$PowerS[i]))
    expect_true(res$PowerS[i] >= 0 && res$PowerS[i] <= 1)
  }
})

test_that("Asymptotic TrtB sanity check gives strong power on seed = 2 toy data", {
  # In the toy generator, TrtB has the largest built-in positive shift.
  # This makes it the strongest treatment for a basic sanity check.
  dat <- mcnemar_example_long(n = 40, seed = 2)
  
  fit <- mcnemar_power(
    data = dat,
    id = "id",
    condition = "condition",
    outcome = "outcome",
    control = "Control",
    treatments = "TrtB",
    methods = c("asymptotic"),
    p_adjust = "holm",
    n_max = 100
  )
  
  res <- fit$results
  expect_equal(nrow(res), 1)
  expect_equal(res$treatment, "TrtB")
  expect_equal(res$n, 40)
  
  # TrtB should have high single-test power and finite 80%/90% sample sizes.
  expect_true(res$PowerS > 0.9)
  expect_true(is.finite(res$N80S))
  expect_true(is.finite(res$N90S))
  expect_true(res$N90S >= res$N80S)
})

test_that("N90S may be NA when target is not reached by chosen n_max", {
  # A weaker treatment and a smaller n_max may fail to reach 90% power,
  # in which case N90S is expected to be NA.
  dat <- mcnemar_example_long(n = 40, seed = 2)
  
  fit <- mcnemar_power(
    data = dat,
    id = "id",
    condition = "condition",
    outcome = "outcome",
    control = "Control",
    treatments = "TrtC",
    methods = c("asymptotic"),
    p_adjust = "holm",
    n_max = 100
  )
  
  res <- fit$results
  expect_equal(nrow(res), 1)
  expect_equal(res$treatment, "TrtC")
  
  # Either N90S is NA (target not reached) or, if finite, it must still
  # be at least as large as N80S.
  expect_true(is.na(res$N90S) || res$N90S >= res$N80S)
})

test_that("Sample-size targets are monotone whenever both are finite", {
  # Run a multi-method fit and verify that 90% sample-size targets are
  # never smaller than 80% targets.
  dat <- mcnemar_example_long(n = 40, seed = 2)
  
  fit <- mcnemar_power(
    data = dat,
    id = "id",
    condition = "condition",
    outcome = "outcome",
    control = "Control",
    methods = c("asymptotic", "midp", "exact_cond"),
    p_adjust = "holm",
    n_max = 150
  )
  
  res <- fit$results
  ok <- is.finite(res$N80S) & is.finite(res$N90S)
  expect_true(all(res$N90S[ok] >= res$N80S[ok]))
})

test_that("TrtB is at least as strong as one smaller-shift treatment under asymptotic PowerS", {
  # This is a broad sanity test, not a strict ranking theorem.
  # TrtA and TrtC have smaller generator shifts than TrtB.
  dat <- mcnemar_example_long(n = 40, seed = 2)
  
  fit <- mcnemar_power(
    data = dat,
    id = "id",
    condition = "condition",
    outcome = "outcome",
    control = "Control",
    methods = c("asymptotic"),
    p_adjust = "holm",
    n_max = 150
  )
  
  asym <- fit$results[fit$results$method_label == "asymptotic", , drop = FALSE]
  expect_true(all(c("TrtA", "TrtB", "TrtC") %in% asym$treatment))
  
  pw <- stats::setNames(asym$PowerS, asym$treatment)
  
  # TrtB should be at least as large as one of the weaker-shift treatments.
  expect_true(pw[["TrtB"]] >= min(pw[c("TrtA", "TrtC")]))
})