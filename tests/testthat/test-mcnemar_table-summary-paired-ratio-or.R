test_that("pairedP single method returns Diff and ci_low/ci_high", {
  skip_if_not_installed("effectsize")
  skip_if_not_installed("contingencytables")
  
  dat <- mcnemar_example_long(n = 50, seed = 1)
  fit <- mcnemar_control(dat, "id", "condition", "outcome", control = "Control")
  tab <- mcnemar_table(fit, dat, methodCI = "marginal")
  
  s <- mcnemar_table_summary(tab, pairedP = TRUE, pairedPCI = "Wald")
  
  expect_false(is.null(s$paired_probabilities))
  expect_true(all(c("control", "treatment", "Diff", "ci_low", "ci_high") %in% names(s$paired_probabilities)))
  expect_equal(nrow(s$paired_probabilities), nrow(tab))
})

test_that("pairedP multiple methods keeps raw Diff and adds method-specific columns", {
  skip_if_not_installed("effectsize")
  skip_if_not_installed("contingencytables")
  
  dat <- mcnemar_example_long(n = 50, seed = 2)
  fit <- mcnemar_control(dat, "id", "condition", "outcome", control = "Control")
  tab <- mcnemar_table(fit, dat, methodCI = "marginal")
  
  s <- mcnemar_table_summary(
    tab,
    pairedP = TRUE,
    pairedPCI = c("Wald", "Newcombe")
  )
  
  expect_true(all(c(
    "Diff",
    "ci_low_Wald", "ci_high_Wald",
    "ci_low_Newcombe", "ci_high_Newcombe"
  ) %in% names(s$paired_probabilities)))
  
  raw_diff <- with(s$table, (b - c) / N)
  expect_equal(s$paired_probabilities$Diff, raw_diff)
})

test_that("pairedP bonferroni adds adjusted columns for single method", {
  skip_if_not_installed("effectsize")
  skip_if_not_installed("contingencytables")
  
  dat <- mcnemar_example_long(n = 50, seed = 3)
  fit <- mcnemar_control(dat, "id", "condition", "outcome", control = "Control")
  tab <- mcnemar_table(fit, dat, methodCI = "marginal")
  
  s <- mcnemar_table_summary(
    tab,
    pairedP = TRUE,
    pairedPCI = "Wald",
    p_adjust = "bonferroni"
  )
  
  expect_true(all(c("ci_low_adj", "ci_high_adj") %in% names(s$paired_probabilities)))
})

test_that("pairedP multiple methods with bonferroni adds method-specific adjusted columns", {
  skip_if_not_installed("effectsize")
  skip_if_not_installed("contingencytables")
  
  dat <- mcnemar_example_long(n = 50, seed = 11)
  fit <- mcnemar_control(dat, "id", "condition", "outcome", control = "Control")
  tab <- mcnemar_table(fit, dat, methodCI = "marginal")
  
  s <- mcnemar_table_summary(
    tab,
    pairedP = TRUE,
    pairedPCI = c("Wald", "Newcombe"),
    p_adjust = "bonferroni"
  )
  
  expect_true(all(c(
    "ci_low_adj_Wald", "ci_high_adj_Wald",
    "ci_low_adj_Newcombe", "ci_high_adj_Newcombe"
  ) %in% names(s$paired_probabilities)))
})

test_that("ratioP returns expected columns and control/treatment proportions", {
  skip_if_not_installed("effectsize")
  skip_if_not_installed("contingencytables")
  
  dat <- mcnemar_example_long(n = 50, seed = 4)
  fit <- mcnemar_control(dat, "id", "condition", "outcome", control = "Control")
  tab <- mcnemar_table(fit, dat, methodCI = "marginal")
  
  s <- mcnemar_table_summary(tab, ratioP = TRUE, ratioCI = "Wald")
  
  expect_false(is.null(s$ratio_probabilities))
  expect_true(all(c("control", "treatment", "Ratio", "ci_low", "ci_high") %in% names(s$ratio_probabilities)))
  
  expected_control <- with(s$table, (c + d) / N)
  expected_treatment <- with(s$table, (b + d) / N)
  
  expect_equal(s$ratio_probabilities$control, expected_control)
  expect_equal(s$ratio_probabilities$treatment, expected_treatment)
})

test_that("ratioP bonferroni adds adjusted columns", {
  skip_if_not_installed("effectsize")
  skip_if_not_installed("contingencytables")
  
  dat <- mcnemar_example_long(n = 50, seed = 5)
  fit <- mcnemar_control(dat, "id", "condition", "outcome", control = "Control")
  tab <- mcnemar_table(fit, dat, methodCI = "marginal")
  
  s <- mcnemar_table_summary(
    tab,
    ratioP = TRUE,
    ratioCI = "BonettPrice",
    p_adjust = "bonferroni"
  )
  
  expect_true(all(c("ci_low_adj", "ci_high_adj") %in% names(s$ratio_probabilities)))
})

test_that("ratioP reports treatment-over-control direction and correct marginal probability columns", {
  skip_if_not_installed("contingencytables")
  
  x <- data.frame(
    a = 10L, b = 8L, c = 2L, d = 20L,
    stringsAsFactors = FALSE,
    row.names = "TrtA"
  )
  attr(x, "mcnemar_table_meta") <- list(
    control = "Control",
    methodCI = "none",
    ci_level = 0.95
  )
  
  s <- mcnemar_table_summary(x, ratioP = TRUE, ratioCI = "Wald")
  
  expect_gt(s$ratio_probabilities$Ratio, 1)
  expect_equal(s$ratio_probabilities$control, (2 + 20) / (10 + 8 + 2 + 20))
  expect_equal(s$ratio_probabilities$treatment, (8 + 20) / (10 + 8 + 2 + 20))
})

test_that("OR=TRUE returns odds_ratios tibble with expected columns", {
  skip_if_not_installed("effectsize")
  skip_if_not_installed("contingencytables")
  
  dat <- mcnemar_example_long(n = 50, seed = 6)
  fit <- mcnemar_control(dat, "id", "condition", "outcome", control = "Control")
  tab <- mcnemar_table(fit, dat, methodCI = "marginal")
  
  s <- mcnemar_table_summary(tab, OR = TRUE, ORCI = "Wilson")
  
  expect_false(is.null(s$odds_ratios))
  expect_true(all(c("control", "treatment", "OR", "ci_low", "ci_high") %in% names(s$odds_ratios)))
  expect_true(isTRUE(s$overview$OR))
  expect_identical(s$overview$ORCI, "Wilson")
})

test_that("all ORCI methods are accepted and return OR output", {
  skip_if_not_installed("effectsize")
  skip_if_not_installed("contingencytables")
  
  dat <- mcnemar_example_long(n = 50, seed = 7)
  fit <- mcnemar_control(dat, "id", "condition", "outcome", control = "Control")
  tab <- mcnemar_table(fit, dat, methodCI = "marginal")
  
  methods <- c(
    "Wilson", "ClopperPearson_midP", "ClopperPearson",
    "Blaker", "Wald", "Wald_Laplace"
  )
  
  for (m in methods) {
    s <- mcnemar_table_summary(tab, OR = TRUE, ORCI = m)
    expect_equal(nrow(s$odds_ratios), nrow(tab))
    expect_true("OR" %in% names(s$odds_ratios))
  }
})

test_that("OR bonferroni adds adjusted columns", {
  skip_if_not_installed("effectsize")
  skip_if_not_installed("contingencytables")
  
  dat <- mcnemar_example_long(n = 50, seed = 8)
  fit <- mcnemar_control(dat, "id", "condition", "outcome", control = "Control")
  tab <- mcnemar_table(fit, dat, methodCI = "marginal")
  
  s <- mcnemar_table_summary(
    tab,
    OR = TRUE,
    ORCI = "Wald",
    p_adjust = "bonferroni"
  )
  
  expect_true(all(c("ci_low_adj", "ci_high_adj") %in% names(s$odds_ratios)))
})

test_that("bootstrap_maxt adds adjusted columns and optional degeneracy flags", {
  skip_if_not_installed("boot")
  skip_if_not_installed("effectsize")
  skip_if_not_installed("contingencytables")
  
  dat <- mcnemar_example_long(n = 50, seed = 9)
  fit <- mcnemar_control(dat, "id", "condition", "outcome", control = "Control")
  tab <- mcnemar_table(fit, dat, methodCI = "marginal")
  
  s_paired <- mcnemar_table_summary(
    tab,
    pairedP = TRUE,
    pairedPCI = "Wald",
    p_adjust = "bootstrap_maxt",
    fit = fit,
    data = dat,
    R = 200,
    seed = 1,
    include_flags = TRUE
  )
  expect_true(all(c("ci_low_adj", "ci_high_adj") %in% names(s_paired$paired_probabilities)))
  expect_true("pairedP_sim_degenerate" %in% names(s_paired$paired_probabilities))
  
  s_ratio <- mcnemar_table_summary(
    tab,
    ratioP = TRUE,
    ratioCI = "Wald",
    p_adjust = "bootstrap_maxt",
    fit = fit,
    data = dat,
    R = 200,
    seed = 1,
    include_flags = TRUE
  )
  expect_true(all(c("ci_low_adj", "ci_high_adj") %in% names(s_ratio$ratio_probabilities)))
  expect_true("ratioP_sim_degenerate" %in% names(s_ratio$ratio_probabilities))
  
  s_or <- mcnemar_table_summary(
    tab,
    OR = TRUE,
    ORCI = "Wald_Laplace",
    p_adjust = "bootstrap_maxt",
    fit = fit,
    data = dat,
    R = 200,
    seed = 1,
    include_flags = TRUE
  )
  expect_true(all(c("ci_low_adj", "ci_high_adj") %in% names(s_or$odds_ratios)))
  expect_true("OR_sim_degenerate" %in% names(s_or$odds_ratios))
})

test_that("bootstrap_maxt errors when fit or data are missing for paired/ratio/OR branches", {
  skip_if_not_installed("effectsize")
  skip_if_not_installed("contingencytables")
  skip_if_not_installed("boot")
  
  dat <- mcnemar_example_long(n = 50, seed = 10)
  fit <- mcnemar_control(dat, "id", "condition", "outcome", control = "Control")
  tab <- mcnemar_table(fit, dat, methodCI = "marginal")
  
  expect_error(
    mcnemar_table_summary(
      tab,
      pairedP = TRUE,
      p_adjust = "bootstrap_maxt"
    ),
    "fit|data|bootstrap_maxt"
  )
  
  expect_error(
    mcnemar_table_summary(
      tab,
      ratioP = TRUE,
      p_adjust = "bootstrap_maxt"
    ),
    "fit|data|bootstrap_maxt"
  )
  
  expect_error(
    mcnemar_table_summary(
      tab,
      OR = TRUE,
      ORCI = "Wald",
      p_adjust = "bootstrap_maxt"
    ),
    "fit|data|bootstrap_maxt"
  )
})

test_that("paired difference direction follows treatment-minus-control sign", {
  skip_if_not_installed("contingencytables")
  
  x <- data.frame(
    a = 10L, b = 8L, c = 2L, d = 20L,
    stringsAsFactors = FALSE,
    row.names = "TrtA"
  )
  attr(x, "mcnemar_table_meta") <- list(
    control = "Control",
    methodCI = "none",
    ci_level = 0.95
  )
  
  s_pos <- mcnemar_table_summary(x, pairedP = TRUE, pairedPCI = "Wald")
  expect_gt(s_pos$paired_probabilities$Diff, 0)
  
  x2 <- x
  x2$b <- 2L
  x2$c <- 8L
  
  s_neg <- mcnemar_table_summary(x2, pairedP = TRUE, pairedPCI = "Wald")
  expect_lt(s_neg$paired_probabilities$Diff, 0)
})

test_that("OR direction follows b/c ordering after reciprocal correction", {
  skip_if_not_installed("contingencytables")
  
  x <- data.frame(
    a = 10L, b = 8L, c = 2L, d = 20L,
    stringsAsFactors = FALSE,
    row.names = "TrtA"
  )
  attr(x, "mcnemar_table_meta") <- list(
    control = "Control",
    methodCI = "none",
    ci_level = 0.95
  )
  
  s_pos <- mcnemar_table_summary(x, OR = TRUE, ORCI = "Wald")
  expect_gt(s_pos$odds_ratios$OR, 1)
  
  x2 <- x
  x2$b <- 2L
  x2$c <- 8L
  
  s_neg <- mcnemar_table_summary(x2, OR = TRUE, ORCI = "Wald")
  expect_lt(s_neg$odds_ratios$OR, 1)
})