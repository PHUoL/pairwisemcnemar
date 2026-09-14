#' Compare a control with one or more treatments using multiple McNemar tests
#'
#' @description
#' Runs paired 2x2 McNemar-type tests comparing a \emph{control} condition with
#' one or more \emph{treatment} conditions with data in \strong{long-format}.
#'
#' @details
#' \code{mcnemar_control()} performs a family of \strong{pairwise matched
#' comparisons} between a control condition and one or more treatment
#' conditions. For each treatment, the function constructs a paired 2x2 table
#' using the binary outcome under control and the binary outcome under that
#' treatment for the same subjects, and then applies the requested
#' McNemar-type test. In each comparison, inference is driven by the
#' discordant pairs.
#'
#' In \code{mcnemar_control()} the outcome levels supplied in \code{outcome_levels}
#' are \code{"No"} treated as failure and \code{"Yes"} as success. The table is
#' created with:
#'
#' * rows = control outcomes
#' * columns = treatment outcomes
#'
#' producing a table in the order:
#'
#' \preformatted{
#'               Treatment: No   Treatment: Yes
#' Control: No          a                b
#' Control: Yes         c                d
#' }
#'
#' In this layout, the \emph{discordant} pairs are the upper-right cell \code{b}
#' (control = No, treatment = Yes) and lower-left cell \code{c}
#' (control = Yes, treatment = No). These are the only cells used by all
#' McNemar-type tests.
#'
#' \strong{Interpreting the McNemar test statistic (\code{Z}) across methods:}
#'
#' Five different methods for the McNemar test can be selected, calling
#' functions from the \pkg{contingencytables} package.
#'
#' \itemize{
#' \item \code{method = "asymptotic"}: \code{Z} is signed because it is
#' computed from the difference of discordant counts (proportional to
#' \code{b - c}).
#'
#' \item \code{method = "cc"}: the continuity-corrected statistic uses
#' \code{abs(b - c)} and therefore returns a non-negative \code{Z}, which
#' should be interpreted as a magnitude-only statistic.
#'
#' \item \code{method = "midp"}: the mid-p test returns a p-value;
#' \code{Z} is not defined and is returned as \code{NA}.
#'
#' \item \code{method = "exact_cond"} and \code{method = "exact_uncond"}:
#' exact methods return p-values only; \code{Z} is not defined and is
#' returned as \code{NA}.
#' }
#'
#' In all cases, the default hypothesis test is two-sided.
#'
#' When multiple treatments are compared against the same control, the
#' resulting p-values form a multiple-testing family. The adjusted p-values
#' returned by \code{mcnemar_control()} control multiplicity across the set
#' of control-versus-treatment hypotheses. \code{p_adjust = "holm"} is the
#' default for planned control-versus-treatment comparisons by providing
#' family-wise error control.
#'
#' @details
#' \strong{Discrete Holm for exact McNemar (optional)}
#'
#' When \code{method = "exact_cond"}, you may set \code{p_adjust = "discrete_holm"}
#' to apply a McNemar-specific discrete Holm step-down adjustment. This procedure
#' exploits the exact conditional null distribution of each test's two-sided
#' p-value, which—conditional on the number of discordant pairs \eqn{m=b+c}—is
#' induced by \eqn{B(m, 0.5)}. At step \eqn{j} of the step-down sequence, the
#' method compares the ordered raw p-value \eqn{p_{(j)}} to the smallest
#' threshold \eqn{t} such that \eqn{\sum_{k=j}^{m} q_{(k)}(t) \le \alpha}, where
#' \eqn{q_{(k)}(t) = \Pr_{H_0}\{p_{(k)} \le t \mid m_{(k)}\}} is computed from the
#' exact Binomial distribution. The adjusted p-values are the minimal \eqn{\alpha}
#' at which each comparison would be rejected under this step-down rule, capped at
#' 1 and enforced to be non-decreasing in rank.
#'
#' \emph{When to use.} The discrete Holm option is most beneficial when many
#' comparisons have small numbers of discordant pairs (sparse or highly discrete
#' settings), where it can be materially less conservative than standard Holm yet
#' still achieves strong FWER control. Because it relies on the exact conditional
#' Binomial model, it is offered only for \code{method = "exact_cond"}.
#'
#' @details
#' \strong{Bootstrap multiplicity adjustment for asymptotic McNemar (optional)}
#'
#' When \code{method = "asymptotic"}, you may set
#' \code{p_adjust = "bootstrap"} to apply a correlation-aware bootstrap
#' multiplicity adjustment for the family of control-versus-treatment McNemar
#' comparisons. This option is intended to align with the bootstrap alternative
#' described by Westfall, Troendle, and Pennello (2010) for multiple McNemar
#' tests in multivariate binary data settings, where several paired comparisons
#' are based on the same subjects.
#'
#' The procedure uses the asymptotic McNemar \code{Z} statistics, resamples
#' subjects with replacement from the complete multivariate binary profiles, and
#' calibrates adjusted p-values using the bootstrap distribution of the maximum
#' absolute centered test statistic across the treatment family. The resulting
#' adjusted p-values are approximate, familywise-error-rate-oriented, and
#' account for the correlation structure induced by repeated observations on the
#' same subjects.
#'
#' \emph{When to use.} The bootstrap option is most useful when
#' control-versus-treatment McNemar comparisons are correlated because all
#' treatments are measured on the same subjects, and a dependence-aware
#' multiplicity adjustment is desired.
#'
#' \emph{Requirements.} Because this bootstrap adjustment is designed for the
#' multivariate repeated-subject setting, it uses only subjects with complete
#' observed outcomes across the full control-plus-treatments family included in
#' the adjustment.
#'
#' @details
#' \strong{Resampling adjustment for mid-\eqn{p} McNemar (optional)}
#'
#' When \code{method = "midp"}, you may set \code{p_adjust = "resampling"} to
#' apply a subject-level resampling-based single-step minP adjustment for the
#' family of control-versus-treatment mid-\eqn{p} McNemar tests. This option is
#' intended as a dependence-aware multiplicity companion to \code{method = "midp"}
#' in repeated-subject settings where all treatment comparisons are based on the
#' same subjects.
#'
#' The procedure uses the raw mid-\eqn{p} values from each treatment-versus-control
#' comparison, resamples complete subject profiles with replacement from the set of
#' subjects with non-missing outcomes across the full control-plus-treatments family,
#' recomputes the family of mid-\eqn{p} values in each resample, and calibrates
#' adjusted p-values from the resampling distribution of the minimum raw
#' mid-\eqn{p} value across the treatment family.
#'
#' \emph{When to use.} The resampling option is most useful when
#' control-versus-treatment comparisons are correlated because all treatments are
#' measured on the same subjects, and a dependence-aware multiplicity adjustment is
#' preferred over standard methods such as Holm.
#'
#' \emph{Requirements.} Because the resampling adjustment is family-based, it uses
#' only subjects with complete observed outcomes across the full control-plus-
#' treatments family included in the adjustment.
#'
#' @details
#' \strong{Appropriate use}
#'
#' \code{mcnemar_control()} is \strong{not} a full repeated-measures
#' model. McNemar's test is a paired binary test of marginal homogeneity for
#' two related responses; it does not explicitly model period effects,
#' treatment order, sequence, carryover, or treatment-by-period interactions.
#' As a result, when the same subjects are observed under control and several
#' treatments, the results should be interpreted as \strong{pairwise marginal
#' differences under the observed study design}, not necessarily as treatment
#' effects isolated from time/order effects.
#'
#' The function is most naturally interpreted when:
#' \itemize{
#' \item each subject contributes one binary response per condition,
#' \item control and treatment responses are measured on the same subjects,
#' \item the order of conditions is randomized or counterbalanced across
#' subjects,
#' \item and carryover is negligible or minimized by design.
#' }
#'
#' In repeated-measures studies where all subjects receive conditions in the
#' same fixed order (for example, control first and treatments afterwards),
#' treatment effects may be confounded with period/order effects. In that
#' setting, the pairwise McNemar comparisons returned by
#' \code{mcnemar_control()} should be interpreted as \strong{matched marginal
#' differences under the observed presentation order}, not necessarily as
#' treatment effects fully separated from time, sequence, or carryover.
#'
#' If the order of control and treatments is randomized or counterbalanced
#' across subjects, the same pairwise control-versus-treatment comparisons
#' remain valid but become more interpretable, because treatment is no longer
#' perfectly tied to a fixed sequence position. Randomization or
#' counterbalancing therefore improves the design-based interpretation of the
#' McNemar workflow, although \code{mcnemar_control()} still does not
#' explicitly estimate period, sequence, or carryover effects.
#'
#' When carryover is scientifically plausible, washout periods can help reduce
#' residual effects from earlier conditions before later conditions are
#' measured. A washout period does not guarantee the absence of carryover, but
#' it can reduce contamination of later responses by earlier treatments and
#' thereby improve the interpretability of repeated-condition comparisons.
#'
#' Randomization, counterbalancing, and washout improve the design-based
#' validity of the multiple pairwise McNemar tests, but they do not convert
#' \code{mcnemar_control()} into a full repeated-binary model. If order,
#' period, sequence, or carryover effects need to be estimated directly, users
#' should consider a repeated-binary framework such as Cochran's Q for a
#' global repeated-condition test, or model-based approaches such as
#' generalized estimating equations (GEE), generalized linear mixed models
#' (GLMM), or conditional logistic regression.
#'
#' Note: Cohen's g and confidence intervals are \strong{not} computed by
#' \code{mcnemar_control()}. Use \code{\link{mcnemar_table}()} to obtain the
#' underlying matched 2x2 cell counts and, if requested, Cohen's g and
#' confidence intervals for each control-versus-treatment comparison.
#'
#' @references
#' Edwards AL (1948). Note on the ``correction for continuity'' in testing the
#' significance of the difference between correlated proportions.
#' \emph{Psychometrika}, 13(3), 185--187.
#'
#' Fagerland MW, Lydersen S, Laake P (2013). The McNemar test for binary
#' matched-pairs data: mid-p and asymptotic are better than exact conditional.
#' \emph{BMC Medical Research Methodology}, 13, 91.
#'
#' Fagerland MW, Lydersen S, Laake P (2017). Statistical Analysis of
#' Contingency Tables. Chapman & Hall/CRC.
#'
#' Fay MP (2020). “Exact McNemar's Test and Matching Confidence Intervals,”
#' \emph{exact2x2} vignette, CRAN.
#'
#' McNemar Q (1947). Note on the sampling error of the difference between
#' correlated proportions or percentages. \emph{Psychometrika}, 12(2), 153--157.
#'
#' Westfall PH, Troendle JF, Pennello G (2010). Multiple McNemar Tests.
#' \emph{Biometrics}, 66(4):1185--1191.
#'
#' @param data A data frame in long format.
#' @param id Name of the subject identifier column (character).
#' @param condition Name of the condition column (character) containing control
#' and treatments.
#' @param outcome Name of the binary outcome column (character).
#' @param control The control level (value in \code{condition}).
#' @param treatments Optional vector of treatment levels. Default: all levels
#' except \code{control}.
#' @param method McNemar method: one of \code{"asymptotic"} (default),
#' \code{"cc"}, \code{"midp"}, \code{"exact_cond"}, or
#' \code{"exact_uncond"}.
#' @param p_adjust Multiple-testing adjustment: \code{"holm"} (default),
#' \code{"BH"}, \code{"BY"}, \code{"hochberg"}, \code{"discrete_holm"},
#' \code{"bootstrap"}, or \code{"resampling"}. \code{"discrete_holm"} is
#' available only when \code{method = "exact_cond"}. \code{"bootstrap"} is
#' available only when \code{method = "asymptotic"}. \code{"resampling"} is
#' available only when \code{method = "midp"}.
#' See \link[stats]{p.adjust}.
#' @param ci Confidence level stored in \code{settings} for downstream confidence
#' interval calculations in \code{\link{mcnemar_table}()}. Default is \code{0.95}.
#' @param outcome_levels Optional length-2 character vector giving the order of
#' the two outcome levels.
#' @param drop_na Logical; if \code{TRUE}, drop pairs missing control or
#' treatment outcome.
#' @param quiet Logical; if \code{FALSE}, message or warn for common edge cases.
#' @param gamma Numeric; Berger and Boos adjustment parameter for
#' \code{method = "exact_uncond"} that is used by
#' \code{\link[contingencytables]{McNemar_exact_unconditional_test_paired_2x2}()}.
#' Default is \code{1e-4}.
#' @param num_pi_values Integer; number of values in the nuisance-parameter grid
#' for \code{method = "exact_uncond"} that is used by
#' \code{\link[contingencytables]{McNemar_exact_unconditional_test_paired_2x2}()}.
#' Default is \code{1000L}.
#' @param R Integer; number of bootstrap / resampling replicates used when
#' \code{p_adjust = "bootstrap"} or \code{p_adjust = "resampling"}.
#' Default is \code{5000L}.
#' @param seed Optional integer seed for reproducible bootstrap / resampling
#' multiplicity adjustment when \code{p_adjust = "bootstrap"} or
#' \code{p_adjust = "resampling"}.
#'
#' @return
#' An object of class \code{"mcnemarControl"} with elements:
#' \itemize{
#' \item \code{results}: a data frame containing one row per
#' treatment-versus-control comparison, including the test statistic
#' (if defined), raw p-value, and multiplicity-adjusted p-value.
#' \item \code{settings}: a list containing the analysis settings used to
#' construct the comparisons.
#' }
#'
#' @examples
#' dat <- mcnemar_example_long(seed = 1)
#' fit <- mcnemar_control(
#'   dat,
#'   id = "id",
#'   condition = "condition",
#'   outcome = "outcome",
#'   control = "Control"
#' )
#'
#' fit
#' summary(fit)
#'
#' # Asymptotic McNemar with bootstrap multiplicity adjustment
#' fit_boot <- mcnemar_control( 
#' dat, 
#' id = "id", 
#' condition = "condition", 
#' outcome = "outcome", 
#' control = "Control", 
#' method = "asymptotic", 
#' p_adjust = "bootstrap", 
#' R = 1000L, 
#' seed = 1 
#' ) 
#' 
#' # Mid-p McNemar with resampling-based minP multiplicity adjustment 
#' fit_midp_resamp <- mcnemar_control( 
#' dat, 
#' id = "id", 
#' condition = "condition", 
#' outcome = "outcome", 
#' control = "Control", 
#' method = "midp", 
#' p_adjust = "resampling", 
#' R = 1000L, 
#' seed = 1 
#' ) 
#'
#' @seealso
#' \itemize{
#' \item \code{\link[rstatix]{pairwise_mcnemar_test}()} in the
#' \pkg{rstatix} package for generic \strong{all-pairs} post hoc McNemar
#' comparisons among multiple related conditions. Compared with that workflow,
#' \code{mcnemar_control()} is intended for planned control-versus-treatment
#' comparisons rather than every possible pair of repeated conditions.
#'
#' \item \code{\link[rstatix]{cochran_qtest}()} in the
#' \pkg{rstatix} package for a Cochran's Q test when the primary aim
#' is a single omnibus repeated-binary test across three or more matched
#' conditions before pairwise follow-up comparisons.
#'
#' \item Repeated-binary modeling approaches such as generalized estimating
#' equations (GEE), generalized linear mixed models (GLMM), or conditional
#' logistic regression when treatment order, period effects, sequence,
#' carryover, or subject-level correlation need to be modeled explicitly.
#' }
#'
#' @export
#' @importFrom stats p.adjust pbinom dbinom
#' @importFrom rlang .data
mcnemar_control <- function(data,
                            id,
                            condition,
                            outcome,
                            control,
                            treatments = NULL,
                            method = c("asymptotic", "cc", "midp", "exact_cond", "exact_uncond"),
                            p_adjust = c("holm", "BH", "BY", "hochberg", "discrete_holm", "bootstrap", "resampling"),
                            ci = 0.95,
                            outcome_levels = NULL,
                            drop_na = TRUE,
                            quiet = FALSE,
                            gamma = 1e-4,
                            num_pi_values = 1000L,
                            R = 5000L, 
                            seed = NULL) {
  
  
  method <- match.arg(method)
  p_adjust <- match.arg(p_adjust) # collapses default vector to "holm"
  
  if (p_adjust == "discrete_holm" && method != "exact_cond") {
    stop("p_adjust = 'discrete_holm' is available only when method = 'exact_cond'.")
  }
  
  if (p_adjust == "bootstrap" && method != "asymptotic") {
    stop("p_adjust = 'bootstrap' is available only when method = 'asymptotic'.")
  }
  
  if (p_adjust == "resampling" && method != "midp") {
    stop("p_adjust = 'resampling' is available only when method = 'midp'.")
  }
  
  if (!is.numeric(R) || length(R) != 1L || is.na(R) || R < 100L) {
    stop("`R` must be a single integer >= 100.") 
  } 
  R <- as.integer(R) 
  if (!is.null(seed)) { 
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) { 
      stop("`seed` must be NULL or a single numeric value.") 
    } 
    seed <- as.integer(seed) 
  } 
  
  if (p_adjust == "hochberg" && !quiet) {
    warning("Hochberg controls FWER under independence/positive dependence; interpret with care when tests are correlated.",
            call. = FALSE)
  }
  if (p_adjust == "BH" && !quiet) {
    message("BH controls the false discovery rate (FDR), whereas Holm/Hochberg control the family-wise error rate (FWER).")
  }
  if (p_adjust == "BY" && !quiet) {
    message("BY controls the false discovery rate (FDR) under arbitrary dependence, but is usually more conservative than BH.")
  }
  
  stopifnot(is.data.frame(data))
  
  for (nm in c(id, condition, outcome)) {
    if (!is.character(nm) || length(nm) != 1L) {
      stop("`id`, `condition`, and `outcome` must be single character column names.")
    }
    if (!nm %in% names(data)) {
      stop(sprintf("Column '%s' not found in `data`.", nm))
    }
  }
  
  # Coerce and validate levels
  cond_vals <- unique(as.character(data[[condition]]))
  control <- as.character(control)
  if (!control %in% cond_vals) {
    stop("`control` level not found in `condition` column.")
  }
  
  if (is.null(treatments)) treatments <- setdiff(cond_vals, control)
  treatments <- as.character(treatments)
  treatments <- treatments[!is.na(treatments)]
  if (length(treatments) < 1L) stop("No treatment levels found (or provided).")
  
  if (is.null(outcome_levels)) {
    outcome_levels <- unique(data[[outcome]])
    outcome_levels <- outcome_levels[!is.na(outcome_levels)]
    if (length(outcome_levels) != 2L) {
      stop("`outcome` must have exactly 2 non-missing levels (or supply `outcome_levels`).")
    }
  } else {
    if (length(outcome_levels) != 2L) stop("`outcome_levels` must be length 2.")
  }
  
  # validate exact_uncond args early
  if (!is.numeric(gamma) || length(gamma) != 1L || is.na(gamma) || gamma < 0) {
    stop("`gamma` must be a single non-negative number.")
  }
  if (length(num_pi_values) != 1L || is.na(num_pi_values) || num_pi_values < 2) {
    stop("`num_pi_values` must be a single integer >= 2.")
  }
  num_pi_values <- as.integer(num_pi_values)
  
  # reshape to wide
  df <- data[as.character(data[[condition]]) %in% c(control, treatments), , drop = FALSE]
  key <- paste(df[[id]], df[[condition]], sep = "__")
  if (any(duplicated(key))) {
    stop("Multiple rows per id-condition detected. Aggregate first or ensure uniqueness.")
  }
  
  wide <- tidyr::pivot_wider(
    df,
    id_cols = dplyr::all_of(id),
    names_from = dplyr::all_of(condition),
    values_from = dplyr::all_of(outcome)
  )
  
  available_conditions <- setdiff(names(wide), id)
  treatments <- intersect(treatments, available_conditions)
  if (length(treatments) < 1L) {
    stop(
      "No valid comparisons were produced. ",
      "After reshaping, available condition columns are: ",
      paste(available_conditions, collapse = ", "),
      ". Check that `control`/`treatments` match the values in the `condition` column."
    )
  }
  
  build_table <- function(x, y, levels) {
    x <- factor(x, levels = levels)
    y <- factor(y, levels = levels)
    tab <- table(x, y)
    if (!all(dim(tab) == c(2L, 2L))) {
      tab2 <- matrix(0L, 2L, 2L, dimnames = list(levels, levels))
      tab2[rownames(tab), colnames(tab)] <- tab
      tab <- tab2
    }
    tab
  }
  
  # Robust p-value extractor (contingencytables objects are list-like)
  extract_pvalue <- function(res) {
    common <- c("P", "p.value", "p_value", "pvalue", "p", "pval", "p_val", "two_sided_P")
    for (nm in common) {
      if (!is.null(res[[nm]]) && is.numeric(res[[nm]])) {
        val <- as.numeric(res[[nm]][1])
        if (is.finite(val) && val >= 0 && val <= 1) return(val)
      }
    }
    
    flat <- unlist(res, recursive = TRUE, use.names = TRUE)
    if (length(flat) == 0 || is.null(names(flat))) return(NA_real_)
    
    nm <- names(flat)
    is_prob <- is.finite(flat) & is.numeric(flat) & flat >= 0 & flat <= 1
    nm_ok <- grepl("p", nm, ignore.case = TRUE) &
      grepl("value|val|prob", nm, ignore.case = TRUE)
    
    cand <- flat[is_prob & nm_ok]
    if (length(cand) > 0) return(as.numeric(cand[[1]]))
    
    cand2 <- flat[is_prob]
    if (length(cand2) == 1) return(as.numeric(cand2[[1]]))
    
    NA_real_
  }
  
  run_mcnemar <- function(tab, method) {
    b <- tab[1, 2]
    c <- tab[2, 1]
    disc <- b + c
    
    if (disc == 0L) {
      if (!quiet) message("Zero discordant cells (b + c = 0): returning p = 1; Z undefined.")
      return(list(Z = NA_real_, p = 1, discordant = disc))
    }
    
    if (method == "asymptotic") {
      res <- contingencytables::McNemar_asymptotic_test_paired_2x2(tab)
      return(list(Z = as.numeric(res$Z), p = as.numeric(res$P), discordant = disc))
    }
    
    if (method == "cc") {
      res <- contingencytables::McNemar_asymptotic_test_CC_paired_2x2(tab)
      return(list(Z = as.numeric(res$Z), p = as.numeric(res$P), discordant = disc))
    }
    
    if (method == "midp") {
      res <- contingencytables::McNemar_midP_test_paired_2x2(tab)
      return(list(Z = NA_real_, p = as.numeric(res$midP), discordant = disc))
    }
    
    if (method == "exact_cond") {
      res <- contingencytables::McNemar_exact_cond_test_paired_2x2(tab)
      pval <- extract_pvalue(res) # exact conditional (two-sided)
      return(list(Z = NA_real_, p = pval, discordant = disc))
    }
    
    if (method == "exact_uncond") {
      res <- contingencytables::McNemar_exact_unconditional_test_paired_2x2(
        tab, gamma = gamma, num_pi_values = num_pi_values
      )
      pval <- extract_pvalue(res)
      return(list(Z = NA_real_, p = pval, discordant = disc))
    }
    
    stop("Unknown method.")
  }
  
  rows <- list()
  
  for (trt in treatments) {
    x <- wide[[control]]
    y <- wide[[trt]]
    
    if (drop_na) {
      keep <- !is.na(x) & !is.na(y)
      x <- x[keep]
      y <- y[keep]
    }
    
    tab <- build_table(x, y, outcome_levels)
    test <- run_mcnemar(tab, method)
    
    rows[[length(rows) + 1L]] <- tibble::tibble(
      control = control,
      treatment = trt,
      n = sum(tab),
      discordant = test$discordant,
      method_label = method,
      Z = as.numeric(test$Z),
      p_value = as.numeric(test$p),
      p_adjust_label = p_adjust
    )
  }
  
  results <- dplyr::bind_rows(rows)
  
  # ---- p-value adjustment ----
  if (p_adjust == "discrete_holm") {
    # custom discrete Holm for exact_cond
    pa <- .discrete_holm_mcnemar_exact_cond(
      p = results$p_value,
      m_disc = results$discordant
    )
    results$p_adjusted <- pa
    
  } else if (p_adjust == "bootstrap") {
    # WTP-style bootstrap adjustment for asymptotic McNemar tests 
    boot_adj <- .wtp_bootstrap_mcnemar_asymptotic( 
      wide = wide, 
      control = control, 
      treatments = treatments, 
      outcome_levels = outcome_levels, 
      R = R, 
      seed = seed, 
      quiet = quiet 
    ) 
    results$p_adjusted <- as.numeric(boot_adj$p_adjusted[results$treatment]) 
  } else if (p_adjust == "resampling") { 
    # resampling-based minP adjustment for raw mid-p McNemar p-values 
    resamp_adj <- .midp_resampling_mcnemar( 
      wide = wide, 
      control = control, 
      treatments = treatments, 
      outcome_levels = outcome_levels, 
      R = R, 
      seed = seed, 
      quiet = quiet 
    ) 
    results$p_adjusted <- as.numeric(resamp_adj$p_adjusted[results$treatment]) 
  } else { 
    results$p_adjusted <- stats::p.adjust(results$p_value, method = p_adjust)
  }
  
  cols <- c(
    "control", "treatment", "n", "discordant",
    "method_label", "Z", "p_value", "p_adjusted", "p_adjust_label")
  results <- dplyr::select(results, dplyr::all_of(cols))
  
  out <- list(
    results = results,
    settings = list(
      id = id,
      condition = condition,
      outcome = outcome,
      control = control,
      treatments = treatments,
      method = method,
      p_adjust = p_adjust,
      ci = ci,
      outcome_levels = outcome_levels,
      drop_na = drop_na,
      gamma = gamma,
      num_pi_values = num_pi_values, 
      R = R, 
      seed = seed 
    )
  )
  
  class(out) <- "mcnemarControl"
  out
}

# ----- Helpers for discrete Holm (exact_cond only) -----

# Standard two-sided exact McNemar/binomial p-value at count x out of m
.p_two_sided_binom <- function(x, m) {
  p_lo <- stats::pbinom(x, m, 0.5, lower.tail = TRUE)
  p_hi <- stats::pbinom(x - 1L, m, 0.5, lower.tail = FALSE)
  p <- 2 * pmin(p_lo, p_hi)
  pmin(p, 1)
}

# q_i(t) = P_H0( p_i <= t | discordant m )
.q_discrete_mcnemar <- function(t, m) {
  if (!is.finite(t) || t <= 0) return(0)
  if (t >= 1) return(1)
  xs <- 0:m
  pvals <- .p_two_sided_binom(xs, m)
  w <- stats::dbinom(xs, m, 0.5)
  sum(w[pvals <= t])
}

# Westfall–Troendle discrete Holm adjusted p-values for exact McNemar
# p : numeric vector of raw exact conditional two-sided p-values
# m_disc : integer vector of discordant totals (b+c) per comparison
.discrete_holm_mcnemar_exact_cond <- function(p, m_disc) {
  stopifnot(length(p) == length(m_disc))
  n <- length(p)
  out <- rep(NA_real_, n)
  ok <- which(is.finite(p) & is.finite(m_disc))
  if (length(ok) == 0) return(out)
  
  ord <- order(p[ok], na.last = TRUE)
  p_ord <- p[ok][ord]
  m_ord <- as.integer(m_disc[ok][ord])
  
  # step-down: s_j = sum_{k=j..n} q_k( p_(j) ); adj_j = max_{i<=j} s_i
  s <- rep(NA_real_, length(p_ord))
  for (j in seq_along(p_ord)) {
    qsum <- 0
    for (k in j:length(p_ord)) {
      qsum <- qsum + .q_discrete_mcnemar(p_ord[j], m_ord[k])
    }
    s[j] <- min(qsum, 1)
  }
  adj <- cummax(s)
  # map back to original order
  out_idx <- integer(length(ok)); out_idx[ord] <- seq_along(ord)
  out_vec <- rep(NA_real_, length(ok)); out_vec[out_idx] <- adj
  out[ok] <- pmin(out_vec, 1)
  out
}

# WTP-style bootstrap multiplicity adjustment for asymptotic McNemar tests
# using complete multivariate binary subject profiles.
.wtp_bootstrap_mcnemar_asymptotic <- function(wide,
                                              control,
                                              treatments,
                                              outcome_levels,
                                              R = 5000L,
                                              seed = NULL,
                                              quiet = FALSE) {
  needed_cols <- c(control, treatments)
  cc <- stats::complete.cases(wide[, needed_cols, drop = FALSE])
  
  if (!all(cc) && !quiet) {
    message("Bootstrap adjustment uses ", sum(cc), " complete subject profiles out of ",
            nrow(wide), " after requiring non-missing control and all treatment outcomes.")
  }
  
  wide2 <- wide[cc, , drop = FALSE]
  n_subj <- nrow(wide2)
  if (n_subj < 2L) {
    stop("Need at least 2 complete subjects for p_adjust = 'bootstrap'.")
  }
  
  .coerce01 <- function(v, levels2) {
    v <- factor(v, levels = levels2)
    if (anyNA(v)) return(NULL)
    as.integer(v == levels2[2L])
  }
  
  for (nm in needed_cols) {
    y01 <- .coerce01(wide2[[nm]], outcome_levels)
    if (is.null(y01)) {
      stop("Could not coerce outcomes to 0/1 for bootstrap adjustment. Check outcome levels and missingness.")
    }
    wide2[[nm]] <- y01
  }
  
  # Asymptotic signed McNemar Z statistic:
  # delta_hat / sqrt(Var_hat(delta_hat)),
  # where delta_hat = mean(D), D = treatment - control in {-1,0,1}.
  stat_fun <- function(dd) {
    ctrl <- dd[[control]]
    vapply(treatments, function(tr) {
      d <- dd[[tr]] - ctrl
      n <- length(d)
      if (n < 2L) return(NA_real_)
      
      delta_hat <- mean(d)
      varD_hat <- mean(d^2) - delta_hat^2
      se_hat <- sqrt(varD_hat / n)
      
      if (!is.finite(se_hat) || se_hat <= 0) {
        return(NA_real_)
      }
      delta_hat / se_hat
    }, numeric(1))
  }
  
  z_obs <- stat_fun(wide2)
  
  if (!is.null(seed)) set.seed(seed)
  
  boot_z <- matrix(NA_real_, nrow = R, ncol = length(treatments),
                   dimnames = list(NULL, treatments))
  
  for (r in seq_len(R)) {
    idx <- sample.int(n_subj, size = n_subj, replace = TRUE)
    boot_z[r, ] <- stat_fun(wide2[idx, , drop = FALSE])
  }
  
  # Center the bootstrap distribution around the observed statistics and
  # use the max absolute centered statistic as the familywise reference.
  z_center <- sweep(boot_z, 2, z_obs, "-")
  
  ok_rows <- apply(z_center, 1, function(z) all(is.finite(z)))
  if (!any(ok_rows)) {
    stop("No valid bootstrap replicates available for p_adjust = 'bootstrap'.")
  }
  
  max_abs <- apply(abs(z_center[ok_rows, , drop = FALSE]), 1, max, na.rm = TRUE)
  
  p_adj <- vapply(seq_along(z_obs), function(j) {
    zj <- abs(z_obs[j])
    if (!is.finite(zj)) return(NA_real_)
    mean(max_abs >= zj)
  }, numeric(1))
  
  names(p_adj) <- treatments
  
  list(
    p_adjusted = p_adj,
    observed_z = stats::setNames(z_obs, treatments),
    n_subjects = n_subj,
    n_boot = length(max_abs)
  )
}

# Resampling-based single-step minP adjustment for raw mid-p McNemar p-values
# using complete multivariate binary subject profiles.
.midp_resampling_mcnemar <- function(wide,
                                     control,
                                     treatments,
                                     outcome_levels,
                                     R = 5000L,
                                     seed = NULL,
                                     quiet = FALSE) {
  needed_cols <- c(control, treatments)
  cc <- stats::complete.cases(wide[, needed_cols, drop = FALSE])
  
  if (!all(cc) && !quiet) {
    message("Resampling adjustment uses ", sum(cc), " complete subject profiles out of ",
            nrow(wide), " after requiring non-missing control and all treatment outcomes.")
  }
  
  wide2 <- wide[cc, , drop = FALSE]
  n_subj <- nrow(wide2)
  if (n_subj < 2L) {
    stop("Need at least 2 complete subjects for p_adjust = 'resampling'.")
  }
  
  .coerce01 <- function(v, levels2) {
    v <- factor(v, levels = levels2)
    if (anyNA(v)) return(NULL)
    as.integer(v == levels2[2L])
  }
  
  for (nm in needed_cols) {
    y01 <- .coerce01(wide2[[nm]], outcome_levels)
    if (is.null(y01)) {
      stop("Could not coerce outcomes to 0/1 for resampling adjustment. Check outcome levels and missingness.")
    }
    wide2[[nm]] <- y01
  }
  
  # observed raw mid-p values from discordant counts
  midp_from_counts <- function(b, c) {
    m <- b + c
    if (!is.finite(m) || m <= 0) return(1)
    
    p_lo <- stats::pbinom(b, m, 0.5, lower.tail = TRUE)
    p_hi <- stats::pbinom(b - 1L, m, 0.5, lower.tail = FALSE)
    p_exact <- 2 * pmin(p_lo, p_hi)
    p_exact <- min(p_exact, 1)
    
    p_obs <- stats::dbinom(b, m, 0.5)
    min(1, p_exact - p_obs)
  }
  
  p_obs <- vapply(treatments, function(tr) {
    ctrl <- wide2[[control]]
    trt <- wide2[[tr]]
    b <- sum(ctrl == 0 & trt == 1)
    c <- sum(ctrl == 1 & trt == 0)
    midp_from_counts(b, c)
  }, numeric(1))
  
  # Subject-level sign-flip resampling of pairwise discordance indicators.
  # For each subject i, one Rademacher sign is drawn and applied across the
  # full treatment family to preserve subject-level dependence structure.
  if (!is.null(seed)) set.seed(seed)
  
  minp_star <- rep(NA_real_, R)
  
  d_mat <- sweep(as.matrix(wide2[, treatments, drop = FALSE]), 1, wide2[[control]], "-")
  
  for (r in seq_len(R)) {
    sgn <- sample(c(-1L, 1L), size = n_subj, replace = TRUE)
    d_star <- d_mat
    
    for (j in seq_len(ncol(d_star))) {
      nz <- d_star[, j] != 0
      d_star[nz, j] <- sgn[nz] * d_star[nz, j]
    }
    
    p_star <- vapply(seq_len(ncol(d_star)), function(j) {
      dj <- d_star[, j]
      b_star <- sum(dj == 1L)
      c_star <- sum(dj == -1L)
      midp_from_counts(b_star, c_star)
    }, numeric(1))
    
    minp_star[r] <- min(p_star, na.rm = TRUE)
  }
  
  p_adj <- vapply(seq_along(p_obs), function(j) {
    pj <- p_obs[j]
    if (!is.finite(pj)) return(NA_real_)
    mean(minp_star <= pj, na.rm = TRUE)
  }, numeric(1))
  
  names(p_adj) <- treatments
  
  list(
    p_adjusted = p_adj,
    p_observed = stats::setNames(p_obs, treatments),
    n_subjects = n_subj,
    n_resampling = R
  )
}

#' Summary for `mcnemar_control()` results
#'
#' @description
#' Produces an aggregated summary of an object returned by
#' \code{\link{mcnemar_control}()}. The summary is intended to complement the
#' comparison-level output printed by \code{print.mcnemarControl()} by providing
#' a higher-level overview of the analysis.
#'
#' @details
#' \code{mcnemar_control_summary()} summarizes the family of pairwise
#' control-versus-treatment McNemar comparisons stored in a
#' \code{"mcnemarControl"} object. In contrast to \code{print(fit)}, which
#' displays the comparison-level results table, \code{summary(fit)} returns a
#' structured object of class \code{"summary.mcnemarControl"} containing:
#'
#' \itemize{
#' \item an aggregated \strong{overview} of the analysis,
#' \item the treatments significant at the raw and adjusted levels,
#' \item \strong{discordant-pair diagnostics},
#' \item \strong{ordered comparisons} ranked by adjusted p-value by default,
#' \item method-specific \strong{notes},
#' \item and a suggested \strong{next step} pointing to
#' \code{\link{mcnemar_table}} with one of the currently supported CI modes
#' (\code{"marginal"}, \code{"simultaneous"}, or
#' \code{"simultaneous_logit"}) for inspection of cell counts, Cohen's
#' \eqn{g}, and confidence intervals.
#' }
#'
#' By default, the ordered comparisons are sorted with the \strong{lowest
#' adjusted p-value first}, because this places the strongest pairwise evidence
#' against the control at the top of the summary output. The ordering can be
#' changed to use raw p-values instead via \code{sort_by = "p_value"}.
#'
#' The discordant diagnostics are descriptive only. The
#' \code{low_discordant_threshold} argument is used to count how many
#' control-versus-treatment comparisons have relatively few discordant pairs
#' (\code{b + c}). This threshold is a heuristic summary aid and should not be
#' interpreted as a formal rule for selecting one McNemar test variant over
#' another.
#'
#' Like \code{\link{mcnemar_control}()} itself, this summary is based on
#' \strong{pairwise matched marginal comparisons}. It does not explicitly model
#' repeated-measures features such as period effects, treatment order, sequence,
#' or carryover. In repeated-measures studies with fixed condition order, the
#' summary should therefore be interpreted as describing pairwise matched
#' differences under the observed design rather than treatment effects fully
#' separated from order or carryover effects.
#'
#' The recommended downstream workflow no longer uses combined CI modes;
#' instead, users should select one of the supported \code{\link{mcnemar_table}()}
#' CI methods according to whether pointwise or multiplicity-adjusted
#' intervals are desired.
#'
#' @param x An object returned by \code{\link{mcnemar_control}()}.
#' @param alpha Significance threshold used for summary counts and treatment
#' listings. Default is \code{0.05}.
#' @param sort_by Which column to use when ordering comparisons:
#' \code{"p_adjusted"} (default) or \code{"p_value"}.
#' @param digits Digits for rounding numeric values shown in summary output.
#' Default is \code{4}.
#' @param low_discordant_threshold Integer threshold used for descriptive
#' summary diagnostics. Comparisons with discordant count \code{b + c} below
#' this value are counted as having low discordance. Default is \code{10}.
#' @param object An object returned by \code{\link{mcnemar_control}()}.
#' @param ... Additional arguments passed through to the summary constructor or
#' ignored by the print method.
#'
#' @return
#' \code{mcnemar_control_summary()} and \code{summary.mcnemarControl()} return
#' an object of class \code{"summary.mcnemarControl"} containing:
#' \itemize{
#' \item \code{overview}: a list with analysis-level counts, including
#' control label, method, p-value adjustment method, confidence level stored
#' in the fit, number of comparisons, number of valid p-values, number of
#' zero-discordant comparisons, and counts of significant raw and adjusted
#' comparisons.
#' \item \code{significant_raw}: character vector of treatment labels with
#' raw p-values below \code{alpha}.
#' \item \code{significant_adjusted}: character vector of treatment labels
#' with adjusted p-values below \code{alpha}.
#' \item \code{discordant}: a list summarizing the discordant counts across
#' comparisons (minimum, quartiles, median, mean, maximum, and the number of
#' comparisons below \code{low_discordant_threshold}).
#' \item \code{ordered}: a data frame of all control-versus-treatment
#' comparisons ordered by \code{sort_by}.
#' \item \code{notes}: character vector of method- or data-specific summary
#' notes.
#' \item \code{next_step}: a character string suggesting a reporting workflow
#' using \code{\link{mcnemar_table}()} with one of the supported CI modes
#' (\code{"marginal"}, \code{"simultaneous"}, or
#' \code{"simultaneous_logit"}).
#' \item \code{results}: the original comparison-level results table from the
#' \code{"mcnemarControl"} object.
#' }
#'
#' \code{print.summary.mcnemarControl()} prints the summary object in a
#' human-readable report format and returns the object invisibly.
#'
#' @examples
#' dat <- mcnemar_example_long(seed = 1)
#'
#' fit <- mcnemar_control(
#'   dat,
#'   id = "id",
#'   condition = "condition",
#'   outcome = "outcome",
#'   control = "Control"
#' )
#'
#' # Construct summary object
#' s <- summary(fit)
#'
#' # Print aggregated summary
#' s
#'
#' # Order comparisons by raw p-value instead of adjusted p-value
#' summary(fit, sort_by = "p_value")
#'
#' # Inspect cell counts, Cohen's g, and marginal confidence intervals
#' tab_marg <- mcnemar_table(
#'   fit,
#'   dat,
#'   methodCI = "marginal"
#' )
#'
#' # Inspect cell counts, Cohen's g, and simultaneous adjusted confidence
#' # intervals on the raw g-scale
#' tab_sim <- mcnemar_table(
#'   fit,
#'   dat,
#'   methodCI = "simultaneous",
#'   R = 1000,
#'   seed = 1
#' )
#'
#' # Inspect cell counts, Cohen's g, and simultaneous adjusted confidence
#' # intervals based on the logit scale
#' tab_logit <- mcnemar_table(
#'   fit,
#'   dat,
#'   methodCI = "simultaneous_logit",
#'   R = 1000,
#'   seed = 1
#' )
#'
#' @seealso
#' \itemize{
#' \item \code{\link{mcnemar_control}()} for the comparison-level matched test
#' workflow against a control condition.
#' \item \code{\link{mcnemar_table}()} for matched 2x2 cell counts, Cohen's g,
#' marginal confidence intervals, and simultaneous adjusted confidence
#' intervals corresponding to the pairwise control-versus-treatment
#' comparisons.
#' }
#'
#' @export
mcnemar_control_summary <- function(x,
                                    alpha = 0.05,
                                    sort_by = c("p_adjusted", "p_value"),
                                    digits = 4,
                                    low_discordant_threshold = 10,
                                    ...) {
  stopifnot(inherits(x, "mcnemarControl"))
  
  sort_by <- match.arg(sort_by)
  res <- x$results
  s <- x$settings
  
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number in (0,1).")
  }
  
  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) || digits < 0) {
    stop("`digits` must be a single non-negative number.")
  }
  digits <- as.integer(digits)
  
  if (!is.numeric(low_discordant_threshold) || length(low_discordant_threshold) != 1L ||
      is.na(low_discordant_threshold) || low_discordant_threshold < 0) {
    stop("`low_discordant_threshold` must be a single non-negative integer.")
  }
  low_discordant_threshold <- as.integer(low_discordant_threshold)
  
  n_comp <- nrow(res)
  
  valid_p <- sum(is.finite(res$p_value))
  zero_discordant <- sum(res$discordant == 0, na.rm = TRUE)
  
  sig_raw_idx <- which(is.finite(res$p_value) & res$p_value < alpha)
  sig_adj_idx <- which(is.finite(res$p_adjusted) & res$p_adjusted < alpha)
  
  sig_raw <- as.character(res$treatment[sig_raw_idx])
  sig_adj <- as.character(res$treatment[sig_adj_idx])
  
  discordant_vals <- res$discordant
  discordant_summary <- list(
    min = if (all(is.na(discordant_vals))) NA_real_ else min(discordant_vals, na.rm = TRUE),
    q1 = if (all(is.na(discordant_vals))) NA_real_ else as.numeric(stats::quantile(discordant_vals, 0.25, na.rm = TRUE, names = FALSE)),
    median = if (all(is.na(discordant_vals))) NA_real_ else stats::median(discordant_vals, na.rm = TRUE),
    mean = if (all(is.na(discordant_vals))) NA_real_ else mean(discordant_vals, na.rm = TRUE),
    q3 = if (all(is.na(discordant_vals))) NA_real_ else as.numeric(stats::quantile(discordant_vals, 0.75, na.rm = TRUE, names = FALSE)),
    max = if (all(is.na(discordant_vals))) NA_real_ else max(discordant_vals, na.rm = TRUE),
    low_n = sum(discordant_vals < low_discordant_threshold, na.rm = TRUE),
    threshold = low_discordant_threshold
  )
  
  # Ordered comparisons: lowest selected p-value first (ascending)
  ord <- order(res[[sort_by]], res$p_value, na.last = TRUE, decreasing = FALSE)
  ordered_tbl <- res[ord, , drop = FALSE]
  
  # Round for display only
  ordered_tbl_print <- ordered_tbl
  num_cols <- intersect(c("Z", "p_value", "p_adjusted", "n", "discordant"), names(ordered_tbl_print))
  for (nm in num_cols) {
    if (is.numeric(ordered_tbl_print[[nm]])) {
      ordered_tbl_print[[nm]] <- round(ordered_tbl_print[[nm]], digits)
    }
  }
  
  notes <- character(0)
  
  if (s$method %in% c("midp", "exact_cond", "exact_uncond")) {
    notes <- c(notes, sprintf("Z is reported as NA for method = '%s'.", s$method))
  }
  if (s$method == "exact_uncond") {
    notes <- c(notes, sprintf(
      "Unconditional exact method used gamma = %s and num_pi_values = %s.",
      format(s$gamma, trim = TRUE), s$num_pi_values
    ))
  }
  if (s$p_adjust == "bootstrap") {
    notes <- c(notes, sprintf( 
      "Bootstrap multiplicity adjustment used for asymptotic McNemar tests (R = %s%s).", 
      s$R, 
      if (!is.null(s$seed) && !is.na(s$seed)) { 
        paste0(", seed = ", s$seed) 
      } else { 
        "" 
      } 
    )) 
  } 
  if (s$p_adjust == "resampling") {
    notes <- c(notes, sprintf( 
      "Resampling-based minP multiplicity adjustment used for mid-p McNemar tests (R = %s%s).", 
      s$R, 
      if (!is.null(s$seed) && !is.na(s$seed)) { 
        paste0(", seed = ", s$seed) 
      } else { 
        "" 
      }
    ))
  }
  if (s$p_adjust == "BH") {
    notes <- c(notes, "BH controls the false discovery rate (FDR), not the family-wise error rate (FWER).")
  }
  if (s$p_adjust == "BY") {
    notes <- c(notes, "BY controls the false discovery rate (FDR) under arbitrary dependence, but is usually more conservative than BH.")
  }
  if (s$p_adjust == "hochberg") {
    notes <- c(notes, "Hochberg controls FWER under independence/positive dependence assumptions.")
  }
  if (zero_discordant > 0) {
    notes <- c(notes, sprintf("%d comparison(s) had zero discordant pairs.", zero_discordant))
  }
  if (valid_p < n_comp) {
    notes <- c(notes, sprintf("%d comparison(s) had missing p-values.", n_comp - valid_p))
  }
  
  next_step <- paste0(
    "Use mcnemar_table(fit, data, methodCI = \"marginal\"), ",
    "mcnemar_table(fit, data, methodCI = \"simultaneous\"), or ",
    "mcnemar_table(fit, data, methodCI = \"simultaneous_logit\") ",
    "to inspect cell counts, Cohen's g, and confidence intervals."
  )
  
  overview <- list(
    control = s$control,
    method = s$method,
    p_adjust = s$p_adjust,
    ci = s$ci,
    alpha = alpha,
    comparisons = n_comp,
    valid_p = valid_p,
    zero_discordant = zero_discordant,
    significant_raw = length(sig_raw),
    significant_adjusted = length(sig_adj)
  )
  
  out <- list(
    overview = overview,
    significant_raw = sig_raw,
    significant_adjusted = sig_adj,
    discordant = discordant_summary,
    ordered = ordered_tbl_print,
    notes = unique(notes),
    next_step = next_step,
    results = res
  )
  
  class(out) <- "summary.mcnemarControl"
  out
}

#' @rdname mcnemar_control_summary
#' @export
summary.mcnemarControl <- function(object,
                                   alpha = 0.05,
                                   sort_by = c("p_adjusted", "p_value"),
                                   digits = 4,
                                   low_discordant_threshold = 10,
                                   ...) {
  mcnemar_control_summary(
    x = object,
    alpha = alpha,
    sort_by = sort_by,
    digits = digits,
    low_discordant_threshold = low_discordant_threshold,
    ...
  )
}


#' @rdname mcnemar_control_summary
#'
#' @details
#' The \code{print.summary.mcnemarControl()} method prints the aggregated summary
#' returned by \code{summary(fit)} or \code{mcnemar_control_summary(fit)}. It
#' includes the overview section, lists of significant treatments, discordant
#' diagnostics, ordered comparisons, method-specific notes, and the recommended
#' next step.
#'
#' For the ordered-comparisons tibble, significant values in the columns
#' \code{p_value} and \code{p_adjusted} can optionally be highlighted in red.
#'
#' @param digits Number of digits used for display formatting of p-value columns
#'   in the ordered-comparisons tibble. Default is \code{4}.
#' @param color Logical; if \code{TRUE} (default in interactive sessions),
#'   significant values in the \code{p_value} and \code{p_adjusted} columns of
#'   the ordered-comparisons tibble are displayed in red when ANSI color output
#'   appears to be supported. If \code{FALSE}, the tibble is printed without
#'   color highlighting.
#'
#' @export
print.summary.mcnemarControl <- function(x, digits = 4, alpha = NULL, color = interactive(), ...) {
  ov <- x$overview
  
  if (is.null(alpha)) {
    if (!is.null(ov$alpha) && is.finite(ov$alpha)) {
      alpha <- ov$alpha
    } else if (!is.null(ov$ci) && is.finite(ov$ci)) {
      alpha <- 1 - ov$ci
    } else {
      alpha <- 0.05
    }
  }
  
  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) || digits < 0) {
    stop("`digits` must be a single non-negative number.")
  }
  
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number in (0,1).")
  }
  
  if (!is.logical(color) || length(color) != 1L || is.na(color)) {
    stop("`color` must be TRUE or FALSE.")
  }
  
  cat("mcnemarControl summary\n")
  cat(sprintf(" Control: %s\n", ov$control))
  cat(sprintf(" Method: %s\n", ov$method))
  cat(sprintf(" p-adjust: %s\n", ov$p_adjust))
  cat(sprintf(" Alpha: %.3f\n", ov$alpha))
  cat(sprintf(" Comparisons: %d\n\n", ov$comparisons))
  
  cat("Overview\n")
  cat(sprintf(" Valid p-values: %d / %d\n", ov$valid_p, ov$comparisons))
  cat(sprintf(" Zero-discordant comparisons: %d\n", ov$zero_discordant))
  cat(sprintf(" Significant (raw): %d\n", ov$significant_raw))
  cat(sprintf(" Significant (adjusted): %d\n\n", ov$significant_adjusted))
  
  cat("Significant treatments\n")
  if (length(x$significant_adjusted) > 0) {
    cat(" After adjustment: ",
        paste(x$significant_adjusted, collapse = ", "),
        "\n", sep = "")
  } else {
    cat(" After adjustment: None\n")
  }
  
  if (length(x$significant_raw) > 0) {
    cat(" Raw p-value < alpha: ",
        paste(x$significant_raw, collapse = ", "),
        "\n\n", sep = "")
  } else {
    cat(" Raw p-value < alpha: None\n\n")
  }
  
  cat("Discordant diagnostics\n")
  cat(sprintf(" min = %s\n", format(x$discordant$min, trim = TRUE)))
  cat(sprintf(" Q1 = %s\n", format(x$discordant$q1, trim = TRUE)))
  cat(sprintf(" median = %s\n", format(x$discordant$median, trim = TRUE)))
  cat(sprintf(" mean = %s\n", format(x$discordant$mean, trim = TRUE)))
  cat(sprintf(" Q3 = %s\n", format(x$discordant$q3, trim = TRUE)))
  cat(sprintf(" max = %s\n", format(x$discordant$max, trim = TRUE)))
  cat(sprintf(" Comparisons with discordant < %d: %d\n\n",
              x$discordant$threshold, x$discordant$low_n))
  
  cat("Ordered comparisons\n")
  if (nrow(x$ordered) > 0) {
    ord_tbl <- x$ordered
    
    # Rebuild from original numeric results when available
    if (!is.null(x$results) && is.data.frame(x$results)) {
      key_ord <- paste(ord_tbl$control, ord_tbl$treatment, sep = "\r")
      key_res <- paste(x$results$control, x$results$treatment, sep = "\r")
      m <- match(key_ord, key_res)
      
      if (all(!is.na(m))) {
        ord_tbl <- x$results[m, , drop = FALSE]
      }
    }
    
    ord_tbl <- .decorate_mcnemar_tbl(
      df = ord_tbl,
      alpha = alpha,
      digits = as.integer(digits),
      color = color
    )
    
    # Keep tibble-style printing
    print(ord_tbl, n = nrow(ord_tbl), ...)
  } else {
    cat(" None\n")
  }
  
  cat("\n")
  
  if (length(x$notes) > 0) {
    cat("Notes\n")
    for (msg in x$notes) {
      cat(" - ", msg, "\n", sep = "")
    }
    cat("\n")
  }
  
  cat("Next step\n")
  cat(" ", x$next_step, "\n", sep = "")
  invisible(x)
}

# ----- Internal helpers for tibble-aware colored p-value display -----

# Conservative ANSI-color support check
.supports_color <- function() {
  if (nzchar(Sys.getenv("NO_COLOR", unset = ""))) {
    return(FALSE)
  }
  
  if (identical(Sys.getenv("TERM", unset = ""), "dumb")) {
    return(FALSE)
  }
  
  interactive()
}

# Apply red styling only when enabled and supported
.colorize_text <- function(text, color = TRUE) {
  if (!isTRUE(color) || !.supports_color()) {
    return(text)
  }
  
  vapply(text, cli::col_red, character(1), USE.NAMES = FALSE)
}

# Constructor for a custom p-value vector class used only for tibble printing
.new_mcnemar_pval <- function(x, alpha = 0.05, color = interactive(), digits = 4L) {
  if (!is.numeric(x)) {
    stop("`x` must be numeric.")
  }
  
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number in (0,1).")
  }
  
  if (!is.logical(color) || length(color) != 1L || is.na(color)) {
    stop("`color` must be TRUE or FALSE.")
  }
  
  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) || digits < 0) {
    stop("`digits` must be a single non-negative number.")
  }
  
  vctrs::new_vctr(
    x,
    alpha = alpha,
    color = color,
    digits = as.integer(digits),
    class = "mcnemar_pval"
  )
}

# Format method used by pillar/tibble when rendering the column
#' @export
format.mcnemar_pval <- function(x, ..., digits = NULL) {
  vals <- vctrs::vec_data(x)
  
  alpha <- attr(x, "alpha", exact = TRUE)
  color <- attr(x, "color", exact = TRUE)
  digs  <- attr(x, "digits", exact = TRUE)
  
  if (!is.null(digits)) {
    digs <- digits
  }
  
  if (is.null(alpha)) alpha <- 0.05
  if (is.null(color)) color <- interactive()
  if (is.null(digs))  digs  <- 4L
  
  digs <- as.integer(digs)
  
  out <- rep(NA_character_, length(vals))
  
  ok <- is.finite(vals)
  out[ok] <- format(round(vals[ok], digs), nsmall = digs, trim = TRUE)
  
  idx_red <- ok & vals < alpha
  if (any(idx_red)) {
    out[idx_red] <- .colorize_text(out[idx_red], color = color)
  }
  
  out[is.nan(vals)]      <- "NaN"
  out[is.infinite(vals)] <- as.character(vals[is.infinite(vals)])
  
  out
}

# Tell tibble/pillar how to render this class in a tibble column
#' @exportS3Method pillar::pillar_shaft
pillar_shaft.mcnemar_pval <- function(x, ...) {
  pillar::new_pillar_shaft_simple(
    format(x),
    align = "right"
  )
}

# Keep the tibble header looking like a numeric/double column
#' @exportS3Method vctrs::vec_ptype_abbr
vec_ptype_abbr.mcnemar_pval <- function(x, ...) {
  "dbl"
}

#' @exportS3Method vctrs::vec_ptype_full
vec_ptype_full.mcnemar_pval <- function(x, ...) {
  "double"
}


# ----- Internal helper for tibble-aware p-value display -----
# Decorates p-value columns for tibble/pillar printing while leaving the
# underlying results structure otherwise unchanged.
.decorate_mcnemar_tbl <- function(df, alpha = 0.05, digits = 4L, color = interactive()) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data.frame.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number in (0,1).")
  }
  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) || digits < 0) {
    stop("`digits` must be a single non-negative number.")
  }
  if (!is.logical(color) || length(color) != 1L || is.na(color)) {
    stop("`color` must be TRUE or FALSE.")
  }
  
  out <- tibble::as_tibble(df)
  
  # Wrap p-value columns in the custom vector class so tibble/pillar will
  # format them using format.mcnemar_pval() and pillar_shaft.mcnemar_pval().
  p_cols <- intersect(c("p_value", "p_adjusted"), names(out))
  for (nm in p_cols) {
    if (is.numeric(out[[nm]])) {
      out[[nm]] <- .new_mcnemar_pval(
        out[[nm]],
        alpha = alpha,
        color = color,
        digits = as.integer(digits)
      )
    }
  }
  
  out
}

#' @rdname mcnemar_control
#'
#' @description
#' Prints the comparison-level results stored in an object returned by
#' \code{mcnemar_control()}. This is the compact row-by-row console display of
#' control-versus-treatment McNemar comparisons.
#'
#' @param x An object returned by \code{mcnemar_control()}.
#' @param digits Number of digits used for display formatting of p-value columns.
#' Default is \code{4}.
#' @param alpha Optional significance threshold used only for display
#' highlighting. If \code{NULL} (default), \code{alpha = 1 - x$settings$ci}.
#' @param color Logical; if \code{TRUE} (default in interactive sessions),
#' significant values in the \code{p_value} and \code{p_adjusted} columns of
#' the comparison-level results tibble are displayed in red when ANSI color
#' output appears to be supported. If \code{FALSE}, the tibble is printed
#' without color highlighting.
#' @param ... Additional arguments passed to tibble printing.
#'
#' @return
#' The object \code{x}, returned invisibly.
#'
#' @seealso
#' \itemize{
#' \item \code{\link{summary.mcnemarControl}()} for an aggregated summary view.
#' \item \code{\link{mcnemar_control_summary}()} for constructing the summary
#' object explicitly.
#' \item \code{\link{mcnemar_table}()} for matched 2x2 tables, Cohen's g, and
#' optional confidence intervals corresponding to the printed comparisons.
#' }
#'
#' @export
print.mcnemarControl <- function(x, digits = 4, alpha = NULL, color = interactive(), ...) {
  if (!inherits(x, "mcnemarControl")) {
    stop("`x` must be an object of class 'mcnemarControl'.")
  }
  if (is.null(x$settings) || !is.list(x$settings)) {
    stop("`x$settings` is missing or invalid.")
  }
  if (is.null(x$results) || !is.data.frame(x$results)) {
    stop("`x$results` is missing or invalid.")
  }
  
  s <- x$settings
  
  cat("mcnemarControl results\n")
  cat(sprintf(" Control: %s\n", s$control))
  cat(sprintf(" Method: %s\n", s$method))
  cat(sprintf(" p-adjust: %s\n", s$p_adjust))
  cat(sprintf(" CI level: %.3f\n", s$ci))
  cat(sprintf(" Comparisons: %d\n\n", nrow(x$results)))
  
  if (is.null(alpha)) {
    alpha <- 1 - s$ci
  }
  
  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) || digits < 0) {
    stop("`digits` must be a single non-negative number.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number in (0,1).")
  }
  if (!is.logical(color) || length(color) != 1L || is.na(color)) {
    stop("`color` must be TRUE or FALSE.")
  }
  
  res <- .decorate_mcnemar_tbl(
    df = x$results,
    alpha = alpha,
    digits = as.integer(digits),
    color = color
  )
  
  # Keep tibble-style printing
  print(res, n = nrow(res), ...)
  
  invisible(x)
}