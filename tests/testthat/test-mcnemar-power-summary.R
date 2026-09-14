test_that("summary.mcnemarPower() returns expected class and top-level structure", {
  # Fit a simple asymptotic power object, then summarize under dep_rand.
  dat <- mcnemar_example_long(seed = 1)
  
  fit <- mcnemar_power(
    data = dat,
    id = "id",
    condition = "condition",
    outcome = "outcome",
    control = "Control",
    methods = c("asymptotic"),
    p_adjust = "holm",
    n_max = 100
  )
  
  s <- summary(
    fit,
    data = dat,
    comparisons = c("TrtA", "TrtB"),
    scenarios = c("dep_rand"),
    nsim = 50,
    seed = 123
  )
  
  # Check the summary object class and main components.
  expect_s3_class(s, "summary.mcnemarPower")
  expect_true("single_test" %in% names(s))
  expect_true("adjusted_power" %in% names(s))
  expect_true("settings" %in% names(s))
  expect_true("comparisons" %in% names(s))
  expect_true("individual_treatment" %in% names(s))
  expect_true("scenarios" %in% names(s))
  expect_true("nsim" %in% names(s))
})

test_that("individual_treatment defaults to comparisons", {
  # In the current implementation, if individual_treatment is NULL,
  # it defaults to the values supplied in comparisons.
  dat <- mcnemar_example_long(seed = 1)
  
  fit <- mcnemar_power(
    data = dat,
    id = "id",
    condition = "condition",
    outcome = "outcome",
    control = "Control",
    methods = c("asymptotic"),
    p_adjust = "holm",
    n_max = 100
  )
  
  s <- summary(
    fit,
    data = dat,
    comparisons = c("TrtA", "TrtB"),
    scenarios = c("dep_rand"),
    nsim = 20,
    seed = 1
  )
  
  expect_equal(sort(s$individual_treatment), sort(c("TrtA", "TrtB")))
})

test_that("summary scenario block contains expected dep_rand components", {
  # Check the structure of the dep_rand scenario block.
  dat <- mcnemar_example_long(seed = 1)
  
  fit <- mcnemar_power(
    data = dat,
    id = "id",
    condition = "condition",
    outcome = "outcome",
    control = "Control",
    methods = c("asymptotic"),
    p_adjust = "holm",
    n_max = 100
  )
  
  s <- summary(
    fit,
    data = dat,
    comparisons = c("TrtA", "TrtB"),
    scenarios = c("dep_rand"),
    nsim = 50,
    seed = 123
  )
  
  ss <- s$adjusted_power$asymptotic$dep_rand
  
  # individual_adjusted should be a tibble/data frame with one row per treatment.
  expect_true(is.data.frame(ss$individual_adjusted))
  expect_true(all(c("scenario", "treatment", "PowerAdj") %in% names(ss$individual_adjusted)))
  
  # subset_any, subset_all, and all_treatments should be scalar probabilities.
  expect_true(is.numeric(ss$subset_any) && length(ss$subset_any) == 1)
  expect_true(is.numeric(ss$subset_all) && length(ss$subset_all) == 1)
  expect_true(is.numeric(ss$all_treatments) && length(ss$all_treatments) == 1)
  
  # Logical probability ordering sanity:
  # any-subset should be at least as likely as all-subset,
  # and usually at least as likely as all-treatments.
  expect_true(ss$subset_any >= ss$subset_all)
  expect_true(ss$subset_any >= ss$all_treatments || isTRUE(all.equal(ss$subset_any, ss$all_treatments)))
})

test_that("summary sample_size contains four separate tibbles", {
  # The current implementation stores sample-size summaries in a list of
  # four separate tibbles: individual, any_subset, all_subset, all_treatments.
  dat <- mcnemar_example_long(seed = 1)
  
  fit <- mcnemar_power(
    data = dat,
    id = "id",
    condition = "condition",
    outcome = "outcome",
    control = "Control",
    methods = c("asymptotic"),
    p_adjust = "holm",
    n_max = 100
  )
  
  s <- summary(
    fit,
    data = dat,
    comparisons = c("TrtA", "TrtB"),
    scenarios = c("dep_rand"),
    nsim = 30,
    seed = 123
  )
  
  ss <- s$adjusted_power$asymptotic$dep_rand$sample_size
  
  expect_true(is.list(ss))
  expect_true(all(c("individual", "any_subset", "all_subset", "all_treatments") %in% names(ss)))
  
  expect_true(is.data.frame(ss$individual))
  expect_true(is.data.frame(ss$any_subset))
  expect_true(is.data.frame(ss$all_subset))
  expect_true(is.data.frame(ss$all_treatments))
  
  # Each sample-size tibble should use the current column contract.
  expect_true(all(c("scenario", "target_power", "treatment", "sample_size") %in% names(ss$individual)))
  expect_true(all(c("scenario", "target_power", "treatment", "sample_size") %in% names(ss$any_subset)))
  expect_true(all(c("scenario", "target_power", "treatment", "sample_size") %in% names(ss$all_subset)))
  expect_true(all(c("scenario", "target_power", "treatment", "sample_size") %in% names(ss$all_treatments)))
})

test_that("summary.mcnemarPower() respects explicit individual_treatment override", {
  # If individual_treatment is supplied explicitly, only that treatment
  # should appear in the individual sample-size tibble.
  dat <- mcnemar_example_long(seed = 1)
  
  fit <- mcnemar_power(
    data = dat,
    id = "id",
    condition = "condition",
    outcome = "outcome",
    control = "Control",
    methods = c("asymptotic"),
    p_adjust = "holm",
    n_max = 100
  )
  
  s <- summary(
    fit,
    data = dat,
    comparisons = c("TrtA", "TrtB"),
    individual_treatment = "TrtB",
    scenarios = c("dep_rand"),
    nsim = 30,
    seed = 123
  )
  
  expect_equal(s$individual_treatment, "TrtB")
  
  ind_tbl <- s$adjusted_power$asymptotic$dep_rand$sample_size$individual
  expect_true(all(ind_tbl$treatment == "TrtB"))
})

test_that("print.summary.mcnemarPower() runs invisibly", {
  # The summary print method should dispatch and return invisibly.
  dat <- mcnemar_example_long(seed = 1)
  
  fit <- mcnemar_power(
    data = dat,
    id = "id",
    condition = "condition",
    outcome = "outcome",
    control = "Control",
    methods = c("asymptotic"),
    p_adjust = "holm",
    n_max = 100
  )
  
  s <- summary(
    fit,
    data = dat,
    comparisons = c("TrtA", "TrtB"),
    scenarios = c("dep_rand"),
    nsim = 25,
    seed = 1
  )
  out <- capture.output(ret <- print(s))
  expect_s3_class(ret, "summary.mcnemarPower")
  expect_true(length(out) > 0)
})

test_that("mcnemar_power_summary() mirrors summary()", {
  # The convenience wrapper should return the same kind of object as summary().
  dat <- mcnemar_example_long(seed = 1)
  
  fit <- mcnemar_power(
    data = dat,
    id = "id",
    condition = "condition",
    outcome = "outcome",
    control = "Control",
    methods = c("asymptotic"),
    p_adjust = "holm",
    n_max = 100
  )
  
  s1 <- summary(
    fit,
    data = dat,
    comparisons = c("TrtA", "TrtB"),
    scenarios = c("dep_rand"),
    nsim = 25,
    seed = 123
  )
  
  s2 <- mcnemar_power_summary(
    fit,
    data = dat,
    comparisons = c("TrtA", "TrtB"),
    scenarios = c("dep_rand"),
    nsim = 25,
    seed = 123
  )
  
  expect_s3_class(s2, "summary.mcnemarPower")
  expect_equal(class(s1), class(s2))
  expect_equal(names(s1), names(s2))
})