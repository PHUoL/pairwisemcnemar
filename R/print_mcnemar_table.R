#' @rdname mcnemar_table_summary
#' @export
print.summary.mcnemar_table <- function(x, ...) {
  ov <- x$overview
  cat("mcnemar_table summary\n")
  cat(sprintf(" Control: %s\n", ov$control))
  cat(sprintf(" Comparisons: %d\n", ov$comparisons))
  
  if (!is.null(ov$methodCI) && length(ov$methodCI) == 1L && !is.na(ov$methodCI)) {
    cat(sprintf(" methodCI: %s\n", ov$methodCI))
  }
  if (!is.null(ov$ci_level) && length(ov$ci_level) == 1L && !is.na(ov$ci_level)) {
    cat(sprintf(" CI level: %s\n", ov$ci_level))
  }
  if (!is.null(ov$adjusted_ci) && length(ov$adjusted_ci) == 1L) {
    cat(sprintf(" Adjusted CI columns present: %s\n",
                if (isTRUE(ov$adjusted_ci)) "TRUE" else "FALSE"))
  }
  
  cat(sprintf(" Include flags: %s\n",
              if (isTRUE(ov$include_flags)) "TRUE" else "FALSE"))
  
  cat(sprintf(" Paired probabilities: %s\n",
              if (isTRUE(ov$pairedP)) "TRUE" else "FALSE"))
  if (isTRUE(ov$pairedP)) {
    cat(sprintf(" pairedPCI: %s\n", paste(ov$pairedPCI, collapse = ", ")))
    cat(sprintf(" p_adjust: %s\n", ov$p_adjust))
  }
  
  cat(sprintf(" Ratio of proportions: %s\n",
              if (isTRUE(ov$ratioP)) "TRUE" else "FALSE"))
  if (isTRUE(ov$ratioP)) {
    cat(sprintf(" ratioCI: %s\n", ov$ratioCI))
    cat(sprintf(" p_adjust: %s\n", ov$p_adjust))
  }
  
  cat(sprintf(" Conditional odds ratios: %s\n",
              if (isTRUE(ov$OR)) "TRUE" else "FALSE"))
  if (isTRUE(ov$OR)) {
    cat(sprintf(" ORCI: %s\n", ov$ORCI))
    cat(sprintf(" p_adjust: %s\n", ov$p_adjust))
  }
  
  cat("\n")
  print(x$table, n = nrow(x$table), ...)
  
  cat("\n")
  if (!is.null(x$paired_probabilities)) {
    cat("Difference between marginal proportions\n")
    print(x$paired_probabilities, n = nrow(x$paired_probabilities), ...)
    cat("\n")
  }
  
  if (!is.null(x$ratio_probabilities)) {
    cat("Ratio of proportions\n")
    print(x$ratio_probabilities, n = nrow(x$ratio_probabilities), ...)
    cat("\n")
  }
  
  if (!is.null(x$odds_ratios)) {
    cat("Conditional odds ratios\n")
    print(x$odds_ratios, n = nrow(x$odds_ratios), ...)
    cat("\n")
  }
  
  invisible(x)
}

#' @export
print.mcnemar_table <- function(x, ...) {
  meta <- attr(x, "mcnemar_table_meta", exact = TRUE)
  deg  <- attr(x, "degenerate_treatments", exact = TRUE)
  
  # ---- Header block ----
  cat("mcnemar_table\n")
  
  if (is.list(meta)) {
    cat(" Control: ", meta$control, "\n", sep = "")
    cat(" Treatments (K): ", length(meta$treatments), "\n", sep = "")
    cat(" methodCI: ", meta$methodCI, "\n", sep = "")
    
    if (!is.na(meta$ci_level)) {
      cat(" CI level: ", meta$ci_level, "\n", sep = "")
    }
    if (!is.na(meta$R)) {
      cat(" R (bootstrap): ", meta$R, "\n", sep = "")
    }
    if (!is.na(meta$seed) && !is.null(meta$seed)) {
      cat(" seed: ", meta$seed, "\n", sep = "")
    }
    cat(" digits: ", meta$digits, "\n", sep = "")
    
    # Explain CI columns under the revised contract
    if (identical(meta$methodCI, "marginal")) {
      cat(" CI columns: ci_low / ci_high = marginal\n", sep = "")
    }
    if (identical(meta$methodCI, "simultaneous")) {
      cat(" CI columns: ci_low / ci_high = marginal; ci_low_adj / ci_high_adj = simultaneous adjusted (raw g-scale)\n", sep = "")
    }
    if (identical(meta$methodCI, "simultaneous_logit")) {
      cat(" CI columns: ci_low / ci_high = marginal; ci_low_adj / ci_high_adj = simultaneous adjusted (logit scale)\n", sep = "")
    }
  }
  
  # ---- Recommendation line (only when there is at least one degenerate treatment) ----
  if (is.character(deg) && length(deg) > 0) {
    cat(
      " Marginal CI recommended for treatment/s = ",
      paste(deg, collapse = ", "),
      "\n",
      sep = ""
    )
  }
  
  cat("\n")
  
  # Print table itself
  NextMethod()
  
  invisible(x)
}