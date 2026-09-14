test_that("default p_adjust collapses to 'holm' when not supplied", {
  dat <- mcnemar_example_long(seed = 1)
  fit <- mcnemar_control(
    dat, id = "id", condition = "condition", outcome = "outcome",
    control = "Control"  # no p_adjust, no method change
  )
  expect_identical(fit$settings$p_adjust, "holm")
})

test_that("p_adjust='discrete_holm' is accepted only for method='exact_cond'", {
  testthat::skip_if_not_installed("contingencytables")
  
  dat <- mcnemar_example_long(seed = 123)
  
  # OK case
  fit_ok <- mcnemar_control(
    dat, id = "id", condition = "condition", outcome = "outcome",
    control = "Control", method = "exact_cond", p_adjust = "discrete_holm",
    quiet = TRUE
  )
  expect_s3_class(fit_ok, "mcnemarControl")
  expect_identical(fit_ok$settings$method, "exact_cond")
  expect_identical(fit_ok$settings$p_adjust, "discrete_holm")
  
  # Block all other methods
  for (m in c("asymptotic", "cc", "midp", "exact_uncond")) {
    expect_error(
      mcnemar_control(
        dat, id = "id", condition = "condition", outcome = "outcome",
        control = "Control", method = m, p_adjust = "discrete_holm", quiet = TRUE
      ),
      regexp = "p_adjust = 'discrete_holm' is available only when method = 'exact_cond'",
      fixed = TRUE
    )
  }
  
  # Default p_adjust collapses to "holm"
  fit_def <- mcnemar_control(
    dat, id = "id", condition = "condition", outcome = "outcome",
    control = "Control"  # no p_adjust supplied
  )
  expect_identical(fit_def$settings$p_adjust, "holm")
})

test_that("discrete_holm adjusted p-values are monotone and valid (no domination asserted)", {
  testthat::skip_if_not_installed("contingencytables")
  
  # Small multi-treatment dataset to produce a family of p-values
  dat <- mcnemar_example_long(n = 60, seed = 42)
  
  # exact_cond + discrete Holm (method under test)
  fit_disc <- mcnemar_control(
    dat,
    id = "id",
    condition = "condition",
    outcome = "outcome",
    control = "Control",
    method = "exact_cond",
    p_adjust = "discrete_holm",
    quiet = TRUE
  )
  
  # exact_cond + standard Holm (baseline; used only for a smoke comparison)
  fit_holm <- mcnemar_control(
    dat,
    id = "id",
    condition = "condition",
    outcome = "outcome",
    control = "Control",
    method = "exact_cond",
    p_adjust = "holm",
    quiet = TRUE
  )
  
  # Align rows by treatment for element-wise comparisons and stable ordering
  df_disc <- fit_disc$results[order(fit_disc$results$treatment), , drop = FALSE]
  df_holm <- fit_holm$results[order(fit_holm$results$treatment), , drop = FALSE]
  expect_identical(as.character(df_disc$treatment), as.character(df_holm$treatment))
  
  padj_disc <- df_disc$p_adjusted
  padj_holm <- df_holm$p_adjusted
  praw_disc <- df_disc$p_value
  
  # 1) Valid bounds: adjusted p in [0,1]
  idx <- which(is.finite(padj_disc))
  expect_true(all(padj_disc[idx] >= -1e-12 & padj_disc[idx] <= 1 + 1e-12))
  
  # 2) Monotonicity in the step-down order (raw p ascending)
  ord <- order(praw_disc, na.last = TRUE)
  step_vals <- padj_disc[ord]
  step_vals <- step_vals[is.finite(step_vals)]
  # non-decreasing (allow tiny numerical noise)
  expect_true(all(diff(step_vals) >= -1e-12))
  
  # 3) Optional smoke checks versus standard Holm (no domination asserted):
  #    - both vectors same length and aligned
  expect_length(padj_disc, length(padj_holm))
  expect_identical(names(padj_disc), NULL) # ensure no name-induced surprises
  
  #    - both in [0,1]
  idxH <- which(is.finite(padj_holm))
  expect_true(all(padj_holm[idxH] >= -1e-12 & padj_holm[idxH] <= 1 + 1e-12))
  
  # Note: We intentionally DO NOT assert discrete Holm <= Holm element-wise,
  # since discrete Holm is not guaranteed to uniformly dominate standard Holm.
})

test_that("mcnemar_control returns comparisons", {
  dat <- mcnemar_example_long(seed = 1)
  fit <- mcnemar_control(dat, "id", "condition", "outcome", control = "Control")
  
  expect_s3_class(fit, "mcnemarControl")
  expect_true(nrow(fit$results) >= 1)
  
  # Expected columns (after removing Cohen's g + CI from mcnemar_control)
  expect_true(all(
    c("control", "treatment", "n", "discordant",
      "method_label", "Z", "p_value", "p_adjusted", "p_adjust_label") %in% names(fit$results)
  ))
  
  # Ensure effect size / CI columns are NOT produced here anymore
  expect_false(any(c("cohens_g", "ci_low", "ci_high", "ci_label") %in% names(fit$results)))
  
  # Settings still contain default CI level (for mcnemar_table)
  expect_true(is.list(fit$settings))
  expect_true("ci" %in% names(fit$settings))
})

test_that("discrete_holm matches hand-computation on m in {3,4,5}", {
  testthat::skip_if_not_installed("contingencytables")
  
  # ------------------------------------------------------------------
  # Construct a tiny long-format dataset where:
  #   - Control outcome is "No" for all subjects
  #   - Treatment outcomes give discordant totals m = b in {3,4,5}
  #   - Therefore exact conditional two-sided p-values are:
  #       m = 3 -> 0.25
  #       m = 4 -> 0.125
  #       m = 5 -> 0.0625
  # ------------------------------------------------------------------
  N <- 12
  ids <- sprintf("S%02d", 1:N)
  
  ctrl <- tibble::tibble(
    id = ids,
    condition = "Control",
    outcome = "No"
  )
  
  T1 <- tibble::tibble(
    id = ids,
    condition = "T1",
    outcome = ifelse(seq_along(ids) <= 3, "Yes", "No")  # m = 3
  )
  
  T2 <- tibble::tibble(
    id = ids,
    condition = "T2",
    outcome = ifelse(seq_along(ids) <= 4, "Yes", "No")  # m = 4
  )
  
  T3 <- tibble::tibble(
    id = ids,
    condition = "T3",
    outcome = ifelse(seq_along(ids) <= 5, "Yes", "No")  # m = 5
  )
  
  dat <- dplyr::bind_rows(ctrl, T1, T2, T3)
  dat$outcome <- factor(dat$outcome, levels = c("No", "Yes"))
  
  # ------------------------------------------------------------------
  # Fit: discrete Holm and standard Holm
  # ------------------------------------------------------------------
  fit_disc <- mcnemar_control(
    dat,
    id = "id",
    condition = "condition",
    outcome = "outcome",
    control = "Control",
    method = "exact_cond",
    p_adjust = "discrete_holm",
    quiet = TRUE
  )
  
  fit_holm <- mcnemar_control(
    dat,
    id = "id",
    condition = "condition",
    outcome = "outcome",
    control = "Control",
    method = "exact_cond",
    p_adjust = "holm",
    quiet = TRUE
  )
  
  # Align by treatment
  dd <- fit_disc$results[order(fit_disc$results$treatment), , drop = FALSE]
  dh <- fit_holm$results[order(fit_holm$results$treatment), , drop = FALSE]
  
  expect_identical(as.character(dd$treatment), c("T1", "T2", "T3"))
  expect_identical(as.character(dh$treatment), c("T1", "T2", "T3"))
  
  # Sanity-check raw p-values and discordant totals
  expect_equal(dd$p_value, c(0.25, 0.125, 0.0625), tolerance = 1e-12)
  expect_equal(dd$discordant, c(3, 4, 5))
  
  # ------------------------------------------------------------------
  # Hand-compute the discrete Holm adjusted p-values
  # ------------------------------------------------------------------
  
  # Exact two-sided binomial McNemar p-value
  p_two_sided <- function(x, m) {
    p_lo <- stats::pbinom(x, m, 0.5, lower.tail = TRUE)
    p_hi <- stats::pbinom(x - 1L, m, 0.5, lower.tail = FALSE)
    min(2 * min(p_lo, p_hi), 1)
  }
  
  # q_i(t) = P_H0(p_i <= t | m)
  q_discrete <- function(t, m) {
    if (!is.finite(t) || t <= 0) return(0)
    if (t >= 1) return(1)
    
    xs <- 0:m
    pvals <- vapply(xs, p_two_sided, numeric(1L), m = m)
    probs <- stats::dbinom(xs, m, 0.5)
    
    sum(probs[pvals <= t])
  }
  
  p_raw <- dd$p_value
  m_vec <- dd$discordant
  
  # Order by raw p-value ascending (step-down order)
  ord <- order(p_raw, na.last = TRUE)
  p_ord <- p_raw[ord]
  m_ord <- m_vec[ord]
  
  s <- numeric(length(p_ord))
  for (j in seq_along(p_ord)) {
    s[j] <- sum(
      vapply(
        m_ord[j:length(m_ord)],
        function(m) q_discrete(p_ord[j], m),
        numeric(1L)
      )
    )
    s[j] <- min(s[j], 1)
  }
  
  adj_ord <- cummax(s)
  
  # Map back to treatment order
  inv <- integer(length(ord))
  inv[ord] <- seq_along(ord)
  adj_hand <- adj_ord[inv]
  
  # ------------------------------------------------------------------
  # Assertions
  # ------------------------------------------------------------------
  
  # In step-down order (smallest raw p first), expected values:
  #   0.0625, 0.125, 0.25
  expect_equal(adj_ord, c(0.0625, 0.125, 0.25), tolerance = 1e-12)
  
  # Check equality to implementation
  expect_equal(dd$p_adjusted, adj_hand, tolerance = 1e-12)
  
  # Monotone in step-down order
  expect_true(all(diff(adj_ord) >= -1e-12))
  
  # Valid bounds
  expect_true(all(dd$p_adjusted >= -1e-12 & dd$p_adjusted <= 1 + 1e-12))
  
  # Optional: in this specific micro-case, discrete Holm should be <= Holm
  expect_true(all(dd$p_adjusted <= dh$p_adjusted + 1e-12))
})

test_that("midp produces comparisons and p-values", {
  dat <- mcnemar_example_long(seed = 1)
  fit <- mcnemar_control(dat, "id", "condition", "outcome",
                         control = "Control", method = "midp")
  
  expect_gt(nrow(fit$results), 0)
  
  # For midp, Z is undefined in our output, and p_value is the mid-p value
  expect_true(all(is.na(fit$results$Z)))
  expect_true(all(!is.na(fit$results$p_value)))
  
  # Expected columns remain the same
  expect_true(all(
    c("control", "treatment", "n", "discordant",
      "method_label", "Z", "p_value", "p_adjusted", "p_adjust_label") %in% names(fit$results)
  ))
  
  # Ensure effect size / CI columns are NOT produced here anymore
  expect_false(any(c("cohens_g", "ci_low", "ci_high", "ci_label") %in% names(fit$results)))
})

test_that("mcnemar_control computes adjusted p-values", {
  dat <- mcnemar_example_long(seed = 1)
  fit <- mcnemar_control(dat, "id", "condition", "outcome", control = "Control")
  
  expect_true(all(!is.na(fit$results$p_adjusted)))
  expect_true(all(fit$results$p_adjusted >= 0 & fit$results$p_adjusted <= 1))
})

test_that("exact_cond produces comparisons and p-values; Z is NA", {
  testthat::skip_if_not_installed("contingencytables")
  
  dat <- mcnemar_example_long(seed = 1)
  fit <- mcnemar_control(dat, "id", "condition", "outcome",
                         control = "Control", method = "exact_cond")
  
  expect_gt(nrow(fit$results), 0)
  expect_true(all(is.na(fit$results$Z)))
  expect_true(all(!is.na(fit$results$p_value)))
  expect_true(all(fit$results$p_value >= 0 & fit$results$p_value <= 1))
})

test_that("exact_uncond produces comparisons and p-values; Z is NA", {
  testthat::skip_if_not_installed("contingencytables")
  
  dat <- mcnemar_example_long(seed = 1)
  
  # Keep this small for unit tests (unconditional exact can be slower)
  fit <- mcnemar_control(dat, "id", "condition", "outcome",
                         control = "Control", method = "exact_uncond",
                         gamma = 1e-4, num_pi_values = 50L)
  
  expect_gt(nrow(fit$results), 0)
  expect_true(all(is.na(fit$results$Z)))
  expect_true(all(!is.na(fit$results$p_value)))
  expect_true(all(fit$results$p_value >= 0 & fit$results$p_value <= 1))
})

test_that("p_adjust='discrete_holm' is accepted only for method='exact_cond'", {
  testthat::skip_if_not_installed("contingencytables")
  
  dat <- mcnemar_example_long(seed = 123)
  
  # OK: exact_cond + discrete_holm should run
  expect_silent({
    fit_ok <- mcnemar_control(
      dat, id = "id", condition = "condition", outcome = "outcome",
      control = "Control", method = "exact_cond", p_adjust = "discrete_holm",
      quiet = TRUE
    )
  })
  expect_s3_class(fit_ok, "mcnemarControl")
  expect_identical(fit_ok$settings$method, "exact_cond")
  expect_identical(fit_ok$settings$p_adjust, "discrete_holm")
  
  # NOT OK: asymptotic + discrete_holm should error at argument validation
  expect_error(
    mcnemar_control(
      dat, id = "id", condition = "condition", outcome = "outcome",
      control = "Control", method = "asymptotic", p_adjust = "discrete_holm",
      quiet = TRUE
    ),
    regexp = "p_adjust = 'discrete_holm' is available only when method = 'exact_cond'",
    fixed  = TRUE
  )
  
  # Sanity: other p_adjust values should be accepted for non-exact methods
  expect_silent(
    mcnemar_control(
      dat, id = "id", condition = "condition", outcome = "outcome",
      control = "Control", method = "asymptotic", p_adjust = "holm",
      quiet = TRUE
    )
  )
})

test_that("p_adjust='discrete_holm' is accepted only for method='exact_cond'", {
  testthat::skip_if_not_installed("contingencytables")
  
  dat <- mcnemar_example_long(seed = 123)
  
  # OK: exact_cond + discrete_holm should run
  expect_silent({
    fit_ok <- mcnemar_control(
      dat, id = "id", condition = "condition", outcome = "outcome",
      control = "Control", method = "exact_cond", p_adjust = "discrete_holm",
      quiet = TRUE
    )
  })
  expect_s3_class(fit_ok, "mcnemarControl")
  expect_identical(fit_ok$settings$method, "exact_cond")
  expect_identical(fit_ok$settings$p_adjust, "discrete_holm")
  
  # NOT OK: asymptotic + discrete_holm should error at argument validation
  expect_error(
    mcnemar_control(
      dat, id = "id", condition = "condition", outcome = "outcome",
      control = "Control", method = "asymptotic", p_adjust = "discrete_holm",
      quiet = TRUE
    ),
    regexp = "p_adjust = 'discrete_holm' is available only when method = 'exact_cond'",
    fixed  = TRUE
  )
  
  # Sanity: other p_adjust values should be accepted for non-exact methods
  expect_silent(
    mcnemar_control(
      dat, id = "id", condition = "condition", outcome = "outcome",
      control = "Control", method = "asymptotic", p_adjust = "holm",
      quiet = TRUE
    )
  )
})

test_that("p_adjust='BY' is accepted and returns valid adjusted p-values", {
  dat <- mcnemar_example_long(seed = 1)
  
  # Fast methods that do not require contingencytables exact routines
  for (m in c("asymptotic", "cc", "midp")) {
    fit <- mcnemar_control(
      dat,
      id = "id",
      condition = "condition",
      outcome = "outcome",
      control = "Control",
      method = m,
      p_adjust = "BY",
      quiet = TRUE
    )
    
    expect_s3_class(fit, "mcnemarControl")
    expect_identical(fit$settings$p_adjust, "BY")
    expect_identical(fit$settings$method, m)
    expect_true(all(!is.na(fit$results$p_adjusted)))
    expect_true(all(fit$results$p_adjusted >= 0 & fit$results$p_adjusted <= 1))
  }
  
  testthat::skip_if_not_installed("contingencytables")
  
  # exact_cond
  fit_ec <- mcnemar_control(
    dat,
    id = "id",
    condition = "condition",
    outcome = "outcome",
    control = "Control",
    method = "exact_cond",
    p_adjust = "BY",
    quiet = TRUE
  )
  expect_s3_class(fit_ec, "mcnemarControl")
  expect_identical(fit_ec$settings$p_adjust, "BY")
  expect_identical(fit_ec$settings$method, "exact_cond")
  expect_true(all(!is.na(fit_ec$results$p_adjusted)))
  expect_true(all(fit_ec$results$p_adjusted >= 0 & fit_ec$results$p_adjusted <= 1))
  
  # exact_uncond (keep small for unit-test speed)
  fit_eu <- mcnemar_control(
    dat,
    id = "id",
    condition = "condition",
    outcome = "outcome",
    control = "Control",
    method = "exact_uncond",
    p_adjust = "BY",
    gamma = 1e-4,
    num_pi_values = 50L,
    quiet = TRUE
  )
  expect_s3_class(fit_eu, "mcnemarControl")
  expect_identical(fit_eu$settings$p_adjust, "BY")
  expect_identical(fit_eu$settings$method, "exact_uncond")
  expect_true(all(!is.na(fit_eu$results$p_adjusted)))
  expect_true(all(fit_eu$results$p_adjusted >= 0 & fit_eu$results$p_adjusted <= 1))
})