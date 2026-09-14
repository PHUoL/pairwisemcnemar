test_that("mcnemar_table_summary returns expected object structure", {
  x <- data.frame(
    a = c(10L, 12L),
    b = c(4L, 5L),
    c = c(2L, 1L),
    d = c(20L, 18L),
    cohens_g = c(0.17, 0.33),
    ci_low = c(0.01, 0.10),
    ci_high = c(0.30, 0.45),
    stringsAsFactors = FALSE,
    row.names = c("Trt1", "Trt2")
  )
  attr(x, "mcnemar_table_meta") <- list(
    control = "Control",
    methodCI = "marginal",
    ci_level = 0.95
  )
  
  s <- mcnemar_table_summary(x)
  
  expect_s3_class(s, "summary.mcnemar_table")
  expect_true(all(c(
    "overview", "table", "paired_probabilities",
    "ratio_probabilities", "odds_ratios", "original"
  ) %in% names(s)))
  
  expect_null(s$paired_probabilities)
  expect_null(s$ratio_probabilities)
  expect_null(s$odds_ratios)
  
  expect_identical(s$overview$control, "Control")
  expect_false(isTRUE(s$overview$pairedP))
  expect_false(isTRUE(s$overview$ratioP))
  expect_false(isTRUE(s$overview$OR))
})

test_that("paired_probabilities and ratio_probabilities are NULL unless requested", {
  skip_if_not_installed("effectsize")
  
  dat <- mcnemar_example_long(n = 40, seed = 12)
  fit <- mcnemar_control(dat, "id", "condition", "outcome", control = "Control")
  tab <- mcnemar_table(fit, dat, methodCI = "marginal")
  
  s <- mcnemar_table_summary(tab)
  
  expect_null(s$paired_probabilities)
  expect_null(s$ratio_probabilities)
  expect_null(s$odds_ratios)
})

test_that("mcnemar_table_summary recovers treatment names from rownames", {
  x <- data.frame(
    a = 10L, b = 4L, c = 2L, d = 20L,
    stringsAsFactors = FALSE,
    row.names = "TrtA"
  )
  attr(x, "mcnemar_table_meta") <- list(
    control = "Control",
    methodCI = "none",
    ci_level = NA_real_
  )
  
  s <- mcnemar_table_summary(x)
  
  expect_identical(s$table$treatment, "TrtA")
  expect_identical(s$table$control, "Control")
  expect_identical(s$table$N, 36L)
})

test_that("mcnemar_table_summary uses control override when supplied", {
  x <- data.frame(
    a = 10L, b = 4L, c = 2L, d = 20L,
    stringsAsFactors = FALSE,
    row.names = "TrtA"
  )
  attr(x, "mcnemar_table_meta") <- list(
    control = "OriginalControl",
    methodCI = "none",
    ci_level = 0.95
  )
  
  s <- mcnemar_table_summary(x, control = "OverrideControl")
  
  expect_identical(s$overview$control, "OverrideControl")
  expect_identical(s$table$control, "OverrideControl")
})

test_that("mcnemar_table_summary carries adjusted ci columns from x when present", {
  x <- data.frame(
    a = 10L, b = 4L, c = 2L, d = 20L,
    cohens_g = 0.17,
    ci_low = 0.01,
    ci_high = 0.30,
    ci_low_adj = -0.05,
    ci_high_adj = 0.40,
    stringsAsFactors = FALSE,
    row.names = "TrtA"
  )
  attr(x, "mcnemar_table_meta") <- list(
    control = "Control",
    methodCI = "simultaneous",
    ci_level = 0.95
  )
  
  s <- mcnemar_table_summary(x)
  
  expect_true(all(c("ci_low_adj", "ci_high_adj") %in% names(s$table)))
  expect_true(isTRUE(s$overview$adjusted_ci))
})

test_that("mcnemar_table_summary include_flags collapses legacy sim_degenerate columns", {
  x <- data.frame(
    a = c(10L, 12L),
    b = c(4L, 5L),
    c = c(2L, 1L),
    d = c(20L, 18L),
    sim_degenerate_sim = c(TRUE, FALSE),
    sim_degenerate_logit = c(FALSE, TRUE),
    recommend_marginal = c(TRUE, FALSE),
    stringsAsFactors = FALSE,
    row.names = c("Trt1", "Trt2")
  )
  attr(x, "mcnemar_table_meta") <- list(
    control = "Control",
    methodCI = "simultaneous_logit",
    ci_level = 0.95
  )
  
  s <- mcnemar_table_summary(x, include_flags = TRUE)
  
  expect_true("sim_degenerate" %in% names(s$table))
  expect_identical(s$table$sim_degenerate, c(TRUE, TRUE))
  expect_identical(s$table$recommend_marginal, c(TRUE, FALSE))
})

test_that("mcnemar_table_summary rounds numeric columns when digits is supplied", {
  x <- data.frame(
    a = 10L, b = 4L, c = 2L, d = 20L,
    cohens_g = 0.1666667,
    ci_low = 0.0123456,
    ci_high = 0.3012345,
    stringsAsFactors = FALSE,
    row.names = "TrtA"
  )
  attr(x, "mcnemar_table_meta") <- list(
    control = "Control",
    methodCI = "marginal",
    ci_level = 0.95
  )
  
  s <- mcnemar_table_summary(x, digits = 2)
  
  expect_identical(s$table$cohens_g, 0.17)
  expect_identical(s$table$ci_low, 0.01)
  expect_identical(s$table$ci_high, 0.30)
})

test_that("summary.mcnemar_table wrapper dispatches to mcnemar_table_summary", {
  skip_if_not_installed("effectsize")
  
  dat <- mcnemar_example_long(n = 40, seed = 1)
  fit <- mcnemar_control(dat, "id", "condition", "outcome", control = "Control")
  tab <- mcnemar_table(fit, dat, methodCI = "marginal")
  
  s1 <- mcnemar_table_summary(tab)
  s2 <- summary(tab)
  
  expect_s3_class(s2, "summary.mcnemar_table")
  expect_identical(names(s1), names(s2))
  expect_equal(as.data.frame(s1$table), as.data.frame(s2$table))
  expect_equal(s1$overview, s2$overview)
})

test_that("mcnemar_table_summary validates logical arguments and x columns", {
  x <- data.frame(a = 1L, b = 1L, c = 1L, stringsAsFactors = FALSE)
  
  expect_error(
    mcnemar_table_summary(x),
    "must contain columns a, b, c, and d"
  )
  
  x2 <- data.frame(a = 1L, b = 1L, c = 1L, d = 1L, stringsAsFactors = FALSE)
  
  expect_error(mcnemar_table_summary(x2, include_flags = NA), "include_flags")
  expect_error(mcnemar_table_summary(x2, pairedP = NA), "pairedP")
  expect_error(mcnemar_table_summary(x2, ratioP = NA), "ratioP")
  expect_error(mcnemar_table_summary(x2, OR = NA), "OR")
})

test_that("mcnemar_table_summary validates R and digits", {
  x <- data.frame(a = 1L, b = 1L, c = 1L, d = 1L, stringsAsFactors = FALSE)
  
  expect_error(
    mcnemar_table_summary(x, R = 99),
    "R.*>= 100"
  )
  
  expect_error(
    mcnemar_table_summary(x, digits = -1),
    "digits.*non-negative"
  )
})

test_that("mcnemar_table_summary overview reports adjusted_ci correctly", {
  x <- data.frame(
    a = 10L, b = 4L, c = 2L, d = 20L,
    ci_low = 0.01,
    ci_high = 0.30,
    stringsAsFactors = FALSE,
    row.names = "TrtA"
  )
  attr(x, "mcnemar_table_meta") <- list(
    control = "Control",
    methodCI = "marginal",
    ci_level = 0.95
  )
  
  s <- mcnemar_table_summary(x)
  expect_false(isTRUE(s$overview$adjusted_ci))
  
  x$ci_low_adj <- -0.05
  x$ci_high_adj <- 0.40
  s2 <- mcnemar_table_summary(x)
  expect_true(isTRUE(s2$overview$adjusted_ci))
})