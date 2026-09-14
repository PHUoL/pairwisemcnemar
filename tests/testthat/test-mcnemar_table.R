test_that("mcnemar_table default output includes cohens_g and treatment names in rownames", {
  dat <- mcnemar_example_long(n = 40, seed = 1)
  fit <- mcnemar_control(
    dat, "id", "condition", "outcome",
    control = "Control", method = "asymptotic"
  )
  
  tab <- mcnemar_table(fit, dat)
  
  expect_s3_class(tab, "data.frame")
  expect_true(inherits(tab, "mcnemar_table"))
  expect_true(all(c("a", "b", "c", "d", "cohens_g") %in% colnames(tab)))
  expect_false(any(grepl("^ci_", colnames(tab))))
  
  expect_equal(nrow(tab), length(fit$settings$treatments))
  expect_setequal(rownames(tab), fit$settings$treatments)
  expect_true(all(tab$a >= 0 & tab$b >= 0 & tab$c >= 0 & tab$d >= 0))
})

test_that("mcnemar_table include_g=FALSE with methodCI='none' returns bare a,b,c,d only", {
  dat <- mcnemar_example_long(n = 40, seed = 1)
  fit <- mcnemar_control(
    dat, "id", "condition", "outcome",
    control = "Control", method = "asymptotic"
  )
  
  tab <- mcnemar_table(fit, dat, methodCI = "none", include_g = FALSE)
  
  expect_s3_class(tab, "data.frame")
  expect_false(inherits(tab, "mcnemar_table"))
  expect_identical(colnames(tab), c("a", "b", "c", "d"))
  expect_setequal(rownames(tab), fit$settings$treatments)
})

test_that("mcnemar_table include_diff adds diff = b - c", {
  dat <- mcnemar_example_long(n = 40, seed = 1)
  fit <- mcnemar_control(
    dat, "id", "condition", "outcome",
    control = "Control", method = "asymptotic"
  )
  
  tab <- mcnemar_table(fit, dat, include_diff = TRUE)
  
  expect_true(all(c("a", "b", "c", "d", "diff") %in% colnames(tab)))
  expect_equal(tab$diff, tab$b - tab$c)
})

test_that("mcnemar_table include_treatment adds treatment column and default rownames", {
  dat <- mcnemar_example_long(n = 40, seed = 1)
  fit <- mcnemar_control(
    dat, "id", "condition", "outcome",
    control = "Control", method = "asymptotic"
  )
  
  tab <- mcnemar_table(
    fit, dat,
    include_treatment = TRUE,
    include_diff = TRUE
  )
  
  expect_true(all(c("treatment", "a", "b", "c", "d", "diff") %in% colnames(tab)))
  expect_equal(nrow(tab), length(fit$settings$treatments))
  expect_setequal(tab$treatment, fit$settings$treatments)
  expect_equal(rownames(tab), as.character(seq_len(nrow(tab))))
  expect_equal(tab$diff, tab$b - tab$c)
})

test_that("mcnemar_table include_discordant adds discordant = b + c", {
  dat <- mcnemar_example_long(n = 40, seed = 1)
  fit <- mcnemar_control(dat, "id", "condition", "outcome", control = "Control")
  
  tab <- mcnemar_table(fit, dat, include_discordant = TRUE)
  
  expect_true("discordant" %in% colnames(tab))
  expect_equal(tab$discordant, tab$b + tab$c)
})

test_that("mcnemar_table include_g=TRUE adds cohens_g when methodCI='none' and no CI columns", {
  dat <- mcnemar_example_long(n = 40, seed = 1)
  fit <- mcnemar_control(
    dat, "id", "condition", "outcome",
    control = "Control", method = "asymptotic"
  )
  
  tab <- mcnemar_table(fit, dat, methodCI = "none", include_g = TRUE)
  
  expect_true("cohens_g" %in% colnames(tab))
  expect_false(any(grepl("^ci_", colnames(tab))))
})

test_that("mcnemar_table errors if methodCI != 'none' but include_g is FALSE", {
  dat <- mcnemar_example_long(n = 40, seed = 1)
  fit <- mcnemar_control(
    dat, "id", "condition", "outcome",
    control = "Control", method = "asymptotic"
  )
  
  expect_error(
    mcnemar_table(fit, dat, methodCI = "marginal", include_g = FALSE),
    "include_g must be TRUE|CIs imply Cohen's g"
  )
})

test_that("mcnemar_table methodCI='marginal' returns cohens_g and ci_low/ci_high", {
  skip_if_not_installed("effectsize")
  
  dat <- mcnemar_example_long(n = 40, seed = 1)
  fit <- mcnemar_control(
    dat, "id", "condition", "outcome",
    control = "Control", method = "asymptotic"
  )
  
  tab <- mcnemar_table(fit, dat, methodCI = "marginal")
  
  expect_true(all(c("cohens_g", "ci_low", "ci_high") %in% colnames(tab)))
  expect_false(any(c("ci_low_adj", "ci_high_adj") %in% colnames(tab)))
})

test_that("mcnemar_table methodCI='simultaneous' returns marginal and adjusted CI columns", {
  skip_if_not_installed("boot")
  skip_if_not_installed("effectsize")
  
  dat <- mcnemar_example_long(n = 40, seed = 1)
  fit <- mcnemar_control(
    dat, "id", "condition", "outcome",
    control = "Control", method = "asymptotic"
  )
  
  tab <- mcnemar_table(
    fit, dat,
    methodCI = "simultaneous",
    R = 200,
    seed = 1,
    verbose = FALSE
  )
  
  expect_true(all(c(
    "cohens_g", "ci_low", "ci_high",
    "ci_low_adj", "ci_high_adj",
    "sim_degenerate", "recommend_marginal"
  ) %in% colnames(tab)))
})

test_that("mcnemar_table methodCI='simultaneous_logit' returns bounded adjusted CI columns when finite", {
  skip_if_not_installed("boot")
  skip_if_not_installed("effectsize")
  
  dat <- mcnemar_example_long(n = 40, seed = 1)
  fit <- mcnemar_control(dat, "id", "condition", "outcome", control = "Control")
  
  tab <- mcnemar_table(
    fit, dat,
    methodCI = "simultaneous_logit",
    R = 200,
    seed = 1
  )
  
  expect_true(all(c(
    "cohens_g", "ci_low", "ci_high",
    "ci_low_adj", "ci_high_adj",
    "sim_degenerate", "recommend_marginal"
  ) %in% colnames(tab)))
  
  ok <- is.finite(tab$ci_low_adj) & is.finite(tab$ci_high_adj)
  if (any(ok)) {
    expect_true(all(tab$ci_low_adj[ok] >= -0.5 & tab$ci_high_adj[ok] <= 0.5))
  } else {
    succeed("No finite adjusted logit-scale CIs to check; all were degenerate/NA.")
  }
})

test_that("mcnemar_table rejects unsupported methodCI values from older API", {
  dat <- mcnemar_example_long(n = 40, seed = 1)
  fit <- mcnemar_control(dat, "id", "condition", "outcome", control = "Control")
  
  expect_error(mcnemar_table(fit, dat, methodCI = "both"))
  expect_error(mcnemar_table(fit, dat, methodCI = "simultaneous_both"))
  expect_error(mcnemar_table(fit, dat, methodCI = "all"))
})

test_that("mcnemar_table validates fit class and data type", {
  dat <- mcnemar_example_long(n = 40, seed = 1)
  fit <- mcnemar_control(dat, "id", "condition", "outcome", control = "Control")
  
  expect_error(
    mcnemar_table(list(), dat),
    "fit.*mcnemar_control|class 'mcnemarControl'"
  )
  
  expect_error(
    mcnemar_table(fit, as.list(dat)),
    "data.*must be a data.frame"
  )
})

test_that("mcnemar_table errors when required columns are missing", {
  dat <- mcnemar_example_long(n = 40, seed = 1)
  fit <- mcnemar_control(dat, "id", "condition", "outcome", control = "Control")
  
  dat_bad <- dat
  dat_bad$outcome <- NULL
  
  expect_error(
    mcnemar_table(fit, dat_bad),
    "missing required columns"
  )
})

test_that("mcnemar_table errors when duplicated id-condition rows are present", {
  dat <- mcnemar_example_long(n = 40, seed = 1)
  fit <- mcnemar_control(dat, "id", "condition", "outcome", control = "Control")
  
  dat_dup <- rbind(dat, dat[1, , drop = FALSE])
  
  expect_error(
    mcnemar_table(fit, dat_dup),
    "Multiple rows per id-condition detected"
  )
})

test_that("mcnemar_table validates ci_level and digits", {
  dat <- mcnemar_example_long(n = 40, seed = 1)
  fit <- mcnemar_control(dat, "id", "condition", "outcome", control = "Control")
  
  expect_error(
    mcnemar_table(fit, dat, methodCI = "marginal", ci_level = 1),
    "ci_level.*\\(0,1\\)"
  )
  
  expect_error(
    mcnemar_table(fit, dat, digits = -1),
    "digits.*non-negative"
  )
})

test_that("mcnemar_table non-bare output stores metadata and bare output does not", {
  dat <- mcnemar_example_long(n = 40, seed = 1)
  fit <- mcnemar_control(dat, "id", "condition", "outcome", control = "Control")
  
  tab_nonbare <- mcnemar_table(fit, dat, include_discordant = TRUE)
  meta_nonbare <- attr(tab_nonbare, "mcnemar_table_meta", exact = TRUE)
  deg_nonbare <- attr(tab_nonbare, "degenerate_treatments", exact = TRUE)
  
  expect_true(inherits(tab_nonbare, "mcnemar_table"))
  expect_true(is.list(meta_nonbare))
  expect_identical(meta_nonbare$control, "Control")
  expect_true(all(c("treatments", "methodCI", "ci_level", "R", "seed", "digits") %in% names(meta_nonbare)))
  expect_true(is.character(deg_nonbare) || length(deg_nonbare) == 0)
  
  tab_bare <- mcnemar_table(fit, dat, methodCI = "none", include_g = FALSE)
  expect_false(inherits(tab_bare, "mcnemar_table"))
  expect_null(attr(tab_bare, "mcnemar_table_meta", exact = TRUE))
})

test_that("mcnemar_table no-discordance case yields NA cohens_g", {
  dat <- data.frame(
    id = rep(1:10, each = 2),
    condition = rep(c("Control", "TreatmentA"), times = 10),
    outcome = rep("No", 20),
    stringsAsFactors = FALSE
  )
  
  fit <- mcnemar_control(
    dat,
    "id",
    "condition",
    "outcome",
    control = "Control",
    outcome_levels = c("No", "Yes")
  )
  
  tab <- mcnemar_table(fit, dat, methodCI = "none", include_g = TRUE)
  
  expect_true(all(tab$b == 0))
  expect_true(all(tab$c == 0))
  expect_true(all(tab$b + tab$c == 0))
  expect_true(all(is.na(tab$cohens_g)))
})