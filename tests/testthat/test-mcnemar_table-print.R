test_that("print.mcnemar_table prints header for non-bare output", {
  dat <- mcnemar_example_long(n = 40, seed = 1)
  fit <- mcnemar_control(dat, "id", "condition", "outcome", control = "Control")
  
  tab <- mcnemar_table(fit, dat)  # default now includes cohens_g -> non-bare output
  
  expect_true(inherits(tab, "mcnemar_table"))
  
  out <- capture.output(rv <- print(tab))
  
  expect_true(any(grepl("^mcnemar_table$", out)))
  expect_true(any(grepl("^ Control:", out)))
  expect_true(any(grepl("^ Treatments \\(K\\):", out)))
  expect_identical(rv, tab)
})

test_that("bare mcnemar_table output does not use custom print header", {
  dat <- mcnemar_example_long(n = 40, seed = 1)
  fit <- mcnemar_control(dat, "id", "condition", "outcome", control = "Control")
  
  tab <- mcnemar_table(fit, dat, methodCI = "none", include_g = FALSE)
  
  expect_false(inherits(tab, "mcnemar_table"))
  
  out <- capture.output(print(tab))
  expect_false(any(grepl("^mcnemar_table$", out)))
})

test_that("print.summary.mcnemar_table prints overview and main table", {
  x <- data.frame(
    a = 10L, b = 4L, c = 2L, d = 20L,
    cohens_g = 0.17,
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
  
  out <- capture.output(rv <- print(s))
  
  expect_true(any(grepl("^mcnemar_table summary$", out)))
  expect_true(any(grepl("Control: Control", out, fixed = TRUE)))
  expect_true(any(grepl("Comparisons: 1", out, fixed = TRUE)))
  expect_true(any(grepl("Paired probabilities: FALSE", out, fixed = TRUE)))
  expect_true(any(grepl("Ratio of proportions: FALSE", out, fixed = TRUE)))
  expect_true(any(grepl("Conditional odds ratios: FALSE", out, fixed = TRUE)))
  expect_identical(rv, s)
})

test_that("print.summary.mcnemar_table prints paired, ratio, and OR sections when present", {
  skip_if_not_installed("effectsize")
  skip_if_not_installed("contingencytables")
  
  dat <- mcnemar_example_long(n = 50, seed = 2)
  fit <- mcnemar_control(dat, "id", "condition", "outcome", control = "Control")
  tab <- mcnemar_table(fit, dat, methodCI = "marginal")
  
  s <- mcnemar_table_summary(
    tab,
    pairedP = TRUE,
    pairedPCI = "Wald",
    ratioP = TRUE,
    ratioCI = "Wald",
    OR = TRUE,
    ORCI = "Wilson"
  )
  
  out <- capture.output(print(s))
  
  expect_true(any(grepl("Paired probabilities: TRUE", out, fixed = TRUE)))
  expect_true(any(grepl("pairedPCI: Wald", out, fixed = TRUE)))
  expect_true(any(grepl("^Difference between marginal proportions$", out)))
  
  expect_true(any(grepl("Ratio of proportions: TRUE", out, fixed = TRUE)))
  expect_true(any(grepl("ratioCI: Wald", out, fixed = TRUE)))
  expect_true(any(grepl("^Ratio of proportions$", out)))
  
  expect_true(any(grepl("Conditional odds ratios: TRUE", out, fixed = TRUE)))
  expect_true(any(grepl("ORCI: Wilson", out, fixed = TRUE)))
  expect_true(any(grepl("^Conditional odds ratios$", out)))
})

test_that("print.summary.mcnemar_table reports OR p_adjust when OR branch is used", {
  skip_if_not_installed("effectsize")
  skip_if_not_installed("contingencytables")
  
  dat <- mcnemar_example_long(n = 50, seed = 3)
  fit <- mcnemar_control(dat, "id", "condition", "outcome", control = "Control")
  tab <- mcnemar_table(fit, dat, methodCI = "marginal")
  
  s <- mcnemar_table_summary(
    tab,
    OR = TRUE,
    ORCI = "Wald",
    p_adjust = "bonferroni"
  )
  
  out <- capture.output(print(s))
  
  expect_true(any(grepl("Conditional odds ratios: TRUE", out, fixed = TRUE)))
  expect_true(any(grepl("ORCI: Wald", out, fixed = TRUE)))
  expect_true(any(grepl("p_adjust: bonferroni", out, fixed = TRUE)))
})