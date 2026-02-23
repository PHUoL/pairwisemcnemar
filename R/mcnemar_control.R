#' Compare a control with one or more treatments using McNemar tests
#'
#' @description
#' Runs paired 2x2 McNemar tests comparing a *control* condition with one or more
#' *treatment* conditions from **long-format** data.
#'
#' @param data A data frame in long format.
#' @param id Name of the subject identifier column (character).
#' @param condition Name of the condition column (character) containing control and treatments.
#' @param outcome Name of the binary outcome column (character).
#' @param control The control level (value in `condition`).
#' @param treatments Optional vector of treatment levels. Default: all levels except `control`.
#' @param method McNemar method: one of `"asymptotic"` (default), `"cc"`, `"midp"`.
#' @param p_adjust Multiple-testing adjustment: `"holm"` (default), `"BH"`, or `"hochberg"`.
#' @param ci Confidence level for Cohen's *g* CI. Default 0.95.
#' @param outcome_levels Optional length-2 character vector giving outcome levels order.
#' @param drop_na Logical; if `TRUE`, drop pairs missing control or treatment outcome.
#' @param quiet Logical; if `FALSE`, message/warn for common edge cases.
#'
#' @return An object of class `"mcnemarControl"` with elements `results` and `settings`.
#'
#' @examples
#' dat <- mcnemar_example_long(seed = 1)
#' fit <- mcnemar_control(dat, id = "id", condition = "condition", outcome = "outcome",
#'                        control = "Control")
#' summary(fit)
#'
#' @export
#' @importFrom stats p.adjust
#' @importFrom rlang .data
mcnemar_control <- function(data,
                            id,
                            condition,
                            outcome,
                            control,
                            treatments = NULL,
                            method = c("asymptotic", "cc", "midp"),
                            p_adjust = c("holm", "BH", "hochberg"),
                            ci = 0.95,
                            outcome_levels = NULL,
                            drop_na = TRUE,
                            quiet = FALSE) {

  method <- match.arg(method)
  p_adjust <- match.arg(p_adjust)

  if (p_adjust == "hochberg" && !quiet) {
    warning("Hochberg controls FWER under independence/positive dependence; interpret with care when tests are correlated.", call. = FALSE)
  }
  if (p_adjust == "BH" && !quiet) {
    message("BH controls the false discovery rate (FDR), whereas Holm/Hochberg control the family-wise error rate (FWER).")
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

  # Coerce condition values and requested levels to character for robust matching
  cond_vals <- unique(as.character(data[[condition]]))
  control <- as.character(control)

  if (!control %in% cond_vals) {
    stop("`control` level not found in `condition` column.")
  }

  if (is.null(treatments)) {
    treatments <- setdiff(cond_vals, control)
  }
  treatments <- as.character(treatments)
  treatments <- treatments[!is.na(treatments)]

  if (length(treatments) < 1L) {
    stop("No treatment levels found (or provided).")
  }

  if (is.null(outcome_levels)) {
    outcome_levels <- unique(data[[outcome]])
    outcome_levels <- outcome_levels[!is.na(outcome_levels)]
    if (length(outcome_levels) != 2L) {
      stop("`outcome` must have exactly 2 non-missing levels (or supply `outcome_levels`).")
    }
  } else {
    if (length(outcome_levels) != 2L) stop("`outcome_levels` must be length 2.")
  }

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
  # Keep only treatments that actually exist as wide columns
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

  run_mcnemar <- function(tab, method) {
    b <- tab[1, 2]
    c <- tab[2, 1]
    disc <- b + c

    if (disc == 0L) {
      if (!quiet) message("Zero discordant cells (b + c = 0): returning p = 1; Z undefined.")
      if (method == "midp") return(list(Z = NA_real_, P = NA_real_, midP = 1, discordant = disc))
      return(list(Z = 0, P = 1, midP = NA_real_, discordant = disc))
    }

    if (method == "asymptotic") {
      res <- contingencytables::McNemar_asymptotic_test_paired_2x2(tab)
      return(list(Z = res$Z, P = res$P, midP = NA_real_, discordant = disc))
    }
    if (method == "cc") {
      res <- contingencytables::McNemar_asymptotic_test_CC_paired_2x2(tab)
      return(list(Z = res$Z, P = res$P, midP = NA_real_, discordant = disc))
    }
    if (method == "midp") {
      res <- contingencytables::McNemar_midP_test_paired_2x2(tab)
      return(list(Z = NA_real_, P = NA_real_, midP = res$P, discordant = disc))
    }
    stop("Unknown method.")
  }

  run_g <- function(tab, ci) {
    out <- tryCatch(effectsize::cohens_g(tab, ci = ci), error = function(e) NULL)
    if (is.null(out)) return(list(g = 0, g_low = NA_real_, g_high = NA_real_))
    list(g = out$Cohens_g[1], g_low = out$CI_low[1], g_high = out$CI_high[1])
  }

  rows <- list()
  for (trt in treatments) {
    x <- wide[[control]]
    y <- wide[[trt]]

    if (drop_na) {
      keep <- !is.na(x) & !is.na(y)
      x <- x[keep]; y <- y[keep]
    }

    tab <- build_table(x, y, outcome_levels)
    test <- run_mcnemar(tab, method)
    p_raw <- if (method == "midp") test$midP else test$P
    g <- run_g(tab, ci)

    rows[[length(rows) + 1L]] <- tibble::tibble(
      control = control,
      treatment = trt,
      n = sum(tab),
      discordant = test$discordant,
      method_label = method,
      Z = as.numeric(test$Z),
      p_value = as.numeric(p_raw),
      cohens_g = as.numeric(g$g),
      ci_low = as.numeric(g$g_low),
      ci_high = as.numeric(g$g_high),
      p_adjust_label = p_adjust,
      ci_label = paste0(ci * 100, "%")
    )
  }

  results <- dplyr::bind_rows(rows)
  results$p_adjusted <- stats::p.adjust(results$p_value, method = p_adjust)

    
    cols <- c(
      "control", "treatment", "n", "discordant",
      "method_label",
      "Z", "p_value", "p_adjusted",
      "cohens_g", "ci_low", "ci_high",
      "p_adjust_label", "ci_label"
    )
    
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
      drop_na = drop_na
    )
  )
  class(out) <- "mcnemarControl"
  out
}

#' Summary table for `mcnemar_control()` results
#'
#' @param x An object returned by \code{mcnemar_control()}.
#' @param digits Digits for rounding numeric columns.
#' @param ... Unused.
#'
#' @return A tibble summarising each treatment vs control comparison.
#' @export
mcnemar_control_summary <- function(x, digits = 4, ...) {
  stopifnot(inherits(x, "mcnemarControl"))
  res <- x$results
  num_cols <- c("Z", "p_value", "p_adjusted", "cohens_g", "ci_low", "ci_high")
  for (nm in intersect(num_cols, names(res))) res[[nm]] <- round(res[[nm]], digits)
  res
}

#' @export
summary.mcnemarControl <- function(object, ...) {
  mcnemar_control_summary(object, ...)
}

#' @export
print.mcnemarControl <- function(x, ...) {
  cat("mcnemarControl results\n")
  s <- x$settings
  cat(sprintf("  Control: %s\n", s$control))
  cat(sprintf("  Method: %s\n", s$method))
  cat(sprintf("  p-adjust: %s\n", s$p_adjust))
  cat(sprintf("  CI level: %.3f\n", s$ci))
  cat(sprintf("  Comparisons: %d\n\n", nrow(x$results)))
  print(mcnemar_control_summary(x, ...), row.names = FALSE)
  invisible(x)
}
