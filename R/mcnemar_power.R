# ============================================================
# mcnemar_power.R
# Single-file, package-ready implementation
# ============================================================

#' Calculate single-test and multiplicity-adjusted power for control-vs-treatment McNemar comparisons
#'
#' @description
#' Builds paired 2x2 control-versus-treatment tables from long-format binary data,
#' then estimates \strong{single-test McNemar power} for each treatment and each
#' selected McNemar-type test method. Power is computed by \strong{complete enumeration}
#' of all possible paired tables with \eqn{N} pairs, weighting each table by its
#' multinomial probability under the estimated alternative probabilities and summing
#' the probability of all tables whose p-values are below the nominal significance
#' level.
#'
#' The fitted object stores a comparison-level tibble in \code{results} and the
#' analysis options in \code{settings}. A dedicated summary method extends the
#' analysis to \strong{multiplicity-adjusted power} under the p-value adjustment
#' rule chosen by \code{p_adjust}, including scenario-based sensitivity analyses
#' for independence, dependence with randomized/counterbalanced order, and
#' dependence with fixed order.
#'
#' @details
#' \strong{Output:}
#'
#' The printed output displays the key columns from \code{results}, including the
#' paired table counts, estimated discordant probabilities, the single-test power
#' \code{PowerS}, and single-test sample-size predictions \code{N80S} and
#' \code{N90S}. If \code{NA} appears in \code{N80S} or \code{N90S}, the
#' requested target was not reached by the chosen \code{n_max}.
#'
#' \strong{Relationship to \code{mcnemar_control()}:}
#'
#' The data reshaping, paired-table orientation, and interpretation of the discordant
#' cells follow the same conventions as \code{mcnemar_control()}. Specifically,
#' the paired 2x2 table is constructed with:
#'
#' \itemize{
#' \item \strong{rows} = control outcomes
#' \item \strong{columns} = treatment outcomes
#' }
#'
#' producing the table:
#'
#' \preformatted{
#'               Treatment: No   Treatment: Yes
#' Control: No          a               b
#' Control: Yes         c               d
#' }
#'
#' where \code{b} and \code{c} are the discordant cells.
#'
#' \strong{Single-test power (\code{PowerS}):}
#'
#' For each treatment-versus-control comparison and each requested test method,
#' \code{mcnemar_power()} estimates the discordant-pair probabilities
#' \eqn{\hat p_b = b/N} and \eqn{\hat p_c = c/N}. These estimated probabilities
#' define the alternative model for a paired 2x2 table with \eqn{N} pairs.
#' Single-test power is then computed by exact complete enumeration:
#'
#' \deqn{Power = \sum_{b,c} Pr(B=b, C=c) I\{p\text{-value}(b,c) \le \alpha\}}
#'
#' where \eqn{(B,C,N-B-C)} follow a trinomial distribution with probabilities
#' \eqn{(\hat p_b, \hat p_c, 1-\hat p_b-\hat p_c)}.
#'
#' \strong{Single-test sample size predictions:}
#'
#' The columns \code{N80S} and \code{N90S} report the smallest integer sample
#' sizes whose enumerated single-test power reaches the requested targets
#' (default 80% and 90%), keeping the estimated discordant probabilities fixed.
#' 
#' \strong{Special note for \code{method = "exact_uncond"}:}
#'
#' For \code{method = "exact_uncond"}, single-test power is computed differently
#' from the other McNemar methods. The asymptotic, continuity-corrected,
#' mid-P, and exact-conditional methods are evaluated from the discordant-pair
#' counts \code{b} and \code{c}. By contrast, the exact-unconditional method
#' is evaluated from the full paired 2x2 table \code{(a, b, c, d)} using
#' \code{contingencytables::McNemar_exact_unconditional_test_paired_2x2()}.
#'
#' This full-table exact-unconditional calculation is computationally more
#' expensive than the other methods, especially when sample-size targets
#' \code{N80S} and \code{N90S} are requested. The sample-size search for
#' \code{"exact_uncond"} repeatedly recomputes exact power over increasing
#' candidate sample sizes up to \code{n_max}, and each evaluation also depends
#' on the nuisance-parameter grid controlled by \code{num_pi_values}.
#'
#' As a result, \code{method = "exact_uncond"} may run noticeably slower than
#' \code{"asymptotic"}, \code{"cc"}, \code{"midp"}, or \code{"exact_cond"},
#' particularly for larger \code{n_max} or larger \code{num_pi_values}. 
#' Consider starting with a preliminary run for \code{method = "exact_uncond"} with lower 
#' \code{num_pi_values} than default and choose \code{n_max} that is close
#' what has been estimated from one of the other methods.
#'
#' \strong{Multiplicity-adjusted power in \code{summary()}:}
#'
#' \code{summary.mcnemarPower()} extends the single-test analysis to the multiple
#' comparison setting induced by comparing one control with multiple treatments.
#' It estimates individual adjusted power for each treatment, any-subset power
#' for a user-specified set of treatments, all-subset power for that set, and
#' all-treatment power for the full family of comparisons.
#'
#' @param data A data frame in long format.
#' @param id Name of the subject identifier column (character).
#' @param condition Name of the condition column (character) containing control
#'   and treatment levels.
#' @param outcome Name of the binary outcome column (character).
#' @param control The control level (value in \code{condition}).
#' @param treatments Optional vector of treatment levels. Default: all levels
#'   except \code{control}.
#' @param methods Character vector of McNemar methods to evaluate. Allowed values:
#'   \code{"asymptotic"}, \code{"cc"}, \code{"midp"}, \code{"exact_cond"},
#'   and \code{"exact_uncond"}.
#' The first four methods are evaluated from the discordant-pair counts
#' \code{b} and \code{c}. The \code{"exact_uncond"} method is evaluated from
#' the full paired 2x2 table \code{(a, b, c, d)} and may therefore be much
#' slower, especially when sample-size targets are also requested.   
#' @param p_adjust Multiple-testing adjustment rule to be stored with the fit and
#'   used by \code{summary()}. One of \code{"holm"} (default), \code{"BH"},
#'   \code{"BY"}, \code{"hochberg"}, or \code{"discrete_holm"}.
#' @param alpha Nominal significance level used for single-test power calculations.
#'   Default is \code{0.05}.
#' @param ci Confidence level stored in \code{settings} for consistency with the
#'   \code{mcnemar_control()} API. Default is \code{0.95}.
#' @param outcome_levels Optional length-2 character vector giving the order of
#'   the two binary outcome levels.
#' @param order Optional name of a column giving within-subject condition order.
#'   If \code{NULL}, row order within each subject is used when order-sensitive
#'   simulations are required in \code{summary()}.
#' @param drop_na Logical; if \code{TRUE}, drop subject pairs missing control or
#'   treatment outcomes. Default is \code{TRUE}.
#' @param quiet Logical; if \code{FALSE}, issue informative messages/warnings for
#'   common edge cases. Default is \code{FALSE}.
#' @param gamma Numeric; Berger and Boos adjustment parameter for
#'   \code{method = "exact_uncond"}. Default is \code{1e-4}.
#' @param num_pi_values Integer; number of grid values for the nuisance parameter
#' in the exact-unconditional method. Default is \code{1000L}. Larger values may
#' improve resolution of the nuisance-parameter search but also increase
#' computation time.
#' @param strict Logical; if \code{FALSE} (default), two-sided single-test power
#'   excludes rejections in the wrong direction when \eqn{\hat p_c > \hat p_b}.
#'   If \code{TRUE}, any two-sided rejection contributes to power.
#' @param targets Numeric vector of target single-test powers for sample-size
#'   prediction. The first two values are reported as \code{N80S} and \code{N90S}
#'   by default. Default is \code{c(0.80, 0.90)}.
#' @param n_max Maximum sample size to search when predicting single-test sample
#'   size. Default is \code{150L}. For \code{method = "exact_uncond"}, larger values 
#'   of \code{n_max} can increase computation time substantially because exact power 
#'   is recomputed repeatedly at each candidate sample size until the requested target 
#'   is reached or the search limit is exceeded.   
#' @param x An object of class \code{"summary.mcnemarPower"} or related
#' output returned by \code{\link{mcnemar_power}()} or its summary method,
#' used by the corresponding print method.
#'
#' @param digits Optional integer giving the number of digits to print.
#' For print methods, defaults to the value stored in the object when available.
#'
#' @param ... Additional arguments passed to or ignored by S3 methods.
#'
#' @return
#' An object of class \code{"mcnemarPower"}.
#'
#' The returned object contains:
#' \itemize{
#'   \item \code{results}: A tibble containing one row per treatment-versus-control
#'   comparison and per McNemar method.
#'   \item \code{settings}: A list of analysis settings used to construct the fit.
#' }
#'
#' The \code{results} tibble contains:
#' \itemize{
#'   \item \code{control}: The control condition label.
#'   \item \code{treatment}: The treatment condition label.
#'   \item \code{method_label}: The McNemar method used for the single-test power calculation.
#'   \item \code{a}, \code{b}, \code{c}, \code{d}: The paired 2x2 cell counts,
#'   where rows correspond to control outcomes and columns correspond to treatment outcomes.
#'   \item \code{n}: The total number of matched pairs.
#'   \item \code{pb_hat}: The estimated discordant-pair probability \eqn{\hat p_b = b/n}.
#'   \item \code{pc_hat}: The estimated discordant-pair probability \eqn{\hat p_c = c/n}.
#'   \item \code{PowerS}: The estimated single-test power for the specified McNemar
#'   method at the observed sample size \code{n}.
#'   \item \code{N80S}: The smallest integer sample size whose estimated single-test
#'   power reaches 80%, if found within \code{n_max}; otherwise \code{NA}.
#'   \item \code{N90S}: The smallest integer sample size whose estimated single-test
#'   power reaches 90%, if found within \code{n_max}; otherwise \code{NA}.
#' }
#'
#' \strong{Interpreting \code{b} and \code{c}:}
#'
#' The paired table is formed with \strong{rows = control outcomes} and
#' \strong{columns = treatment outcomes}. In this layout:
#' \itemize{
#' \item \code{b} = control = No, treatment = Yes
#' \item \code{c} = control = Yes, treatment = No
#' }
#' so \code{b > c} indicates more pairs improving under treatment than worsening
#' under treatment.
#'
#' \strong{Interpreting \code{N80S} and \code{N90S}:}
#'
#' If \code{NA} appears in \code{N80S} or \code{N90S}, the requested target
#' power was not reached by the chosen value of \code{n_max}. In that case,
#' users may increase \code{n_max} and rerun the function.
#'
#' @seealso \code{\link{summary.mcnemarPower}}, \code{\link{mcnemar_power_summary}},
#'   \code{\link{mcnemar_control}}
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
#' Fagerland MW, Lydersen S, Laake P (2014). Recommended tests and confidence
#' intervals for paired binomial proportions. \emph{Statistics in Medicine},
#' 33(16), 2850--2875.
#'
#' Fay MP (2020). Exact McNemar's Test and Matching Confidence Intervals,
#' \emph{exact2x2} vignette.
#'
#' Westfall PH, Troendle JF, Pennello G (2010). Multiple McNemar Tests.
#' \emph{Biometrics}, 66(4), 1185--1191.
#'
#' @examples
#' dat <- mcnemar_example_long(seed = 1)
#'
#' fit <- mcnemar_power(
#'   data = dat,
#'   id = "id",
#'   condition = "condition",
#'   outcome = "outcome",
#'   control = "Control",
#'   methods = c("asymptotic", "midp", "exact_cond")
#' )
#'
#' fit
#'
#' s <- summary(
#'   fit,
#'   data = dat,
#'   comparisons = c("TrtA", "TrtB"),
#'   scenarios = c("indep", "dep_rand", "dep_fixed"),
#'   nsim = 200,
#'   seed = 1
#' )
#'
#' s
#'
#' @export
mcnemar_power <- function(data,
                          id,
                          condition,
                          outcome,
                          control,
                          treatments = NULL,
                          methods = c("asymptotic", "cc", "midp", "exact_cond", "exact_uncond"),
                          p_adjust = c("holm", "BH", "BY", "hochberg", "discrete_holm"),
                          alpha = 0.05,
                          ci = 0.95,
                          outcome_levels = NULL,
                          order = NULL,
                          drop_na = TRUE,
                          quiet = FALSE,
                          gamma = 1e-4,
                          num_pi_values = 1000L,
                          strict = FALSE,
                          targets = c(0.80, 0.90),
                          n_max = 150L) {
  
  methods <- unique(as.character(methods))
  valid_methods <- c("asymptotic", "cc", "midp", "exact_cond", "exact_uncond")
  if (!all(methods %in% valid_methods)) {
    stop("All `methods` must be one or more of: ", paste(valid_methods, collapse = ", "))
  }
  
  p_adjust <- match.arg(p_adjust)
  
  stopifnot(is.data.frame(data))
  for (nm in c(id, condition, outcome)) {
    if (!is.character(nm) || length(nm) != 1L) {
      stop("`id`, `condition`, and `outcome` must be single character column names.")
    }
    if (!nm %in% names(data)) {
      stop(sprintf("Column '%s' not found in `data`.", nm))
    }
  }
  if (!is.null(order)) {
    if (!is.character(order) || length(order) != 1L || !order %in% names(data)) {
      stop("`order` must be NULL or a single column name present in `data`.")
    }
  }
  
  control <- as.character(control)
  cond_vals <- unique(as.character(data[[condition]]))
  if (!control %in% cond_vals) stop("`control` level not found in `condition` column.")
  
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
  
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number in (0,1).")
  }
  if (!is.numeric(ci) || length(ci) != 1L || is.na(ci) || ci <= 0 || ci >= 1) {
    stop("`ci` must be a single number in (0,1).")
  }
  if (!is.numeric(gamma) || length(gamma) != 1L || is.na(gamma) || gamma < 0) {
    stop("`gamma` must be a single non-negative number.")
  }
  if (length(num_pi_values) != 1L || is.na(num_pi_values) || num_pi_values < 2L) {
    stop("`num_pi_values` must be a single integer >= 2.")
  }
  num_pi_values <- as.integer(num_pi_values)
  
  if (!is.numeric(targets) || any(is.na(targets)) || any(targets <= 0 | targets >= 1)) {
    stop("`targets` must contain one or more numbers in (0,1).")
  }
  if (length(targets) < 2L && !quiet) {
    warning("`targets` has fewer than two values; N80S/N90S will use NA for missing targets.", call. = FALSE)
  }
  if (length(n_max) != 1L || is.na(n_max) || n_max < 2L) {
    stop("`n_max` must be a single integer >= 2.")
  }
  n_max <- as.integer(n_max)
  
  wide <- .build_wide_for_mcnemar_power(
    data = data,
    id = id,
    condition = condition,
    outcome = outcome,
    control = control,
    treatments = treatments
  )
  
  available_conditions <- setdiff(names(wide), id)
  treatments <- intersect(treatments, available_conditions)
  if (length(treatments) < 1L) {
    stop(
      "No valid treatment columns after reshaping. Available condition columns are: ",
      paste(available_conditions, collapse = ", "),
      "."
    )
  }
  
  order_info <- .extract_order_structure_mcnemar_power(
    data = data,
    id = id,
    condition = condition,
    control = control,
    treatments = treatments,
    order = order
  )
  
  cache_env <- new.env(parent = emptyenv())
  out_rows <- list()
  row_idx <- 0L
  
  for (trt in treatments) {
    x <- wide[[control]]
    y <- wide[[trt]]
    
    if (drop_na) {
      keep <- !is.na(x) & !is.na(y)
      x <- x[keep]
      y <- y[keep]
    }
    
    tab <- .build_table_2x2_mcnemar_power(x, y, outcome_levels = outcome_levels)
    
    a <- as.integer(tab[1, 1])
    b <- as.integer(tab[1, 2])
    c <- as.integer(tab[2, 1])
    d <- as.integer(tab[2, 2])
    n <- as.integer(sum(tab))
    
    pa_hat <- a / n
    pb_hat <- b / n
    pc_hat <- c / n
    pd_hat <- d / n
    
    for (mth in methods) {
      row_idx <- row_idx + 1L
      
      if (identical(mth, "exact_uncond")) {
        pw <- .power_single_enum_exact_uncond_full_mcnemar_power(
          pa = pa_hat,
          pb = pb_hat,
          pc = pc_hat,
          pd = pd_hat,
          n = n,
          alpha = alpha,
          gamma = gamma,
          num_pi_values = num_pi_values,
          strict = strict,
          cache_env = cache_env
        )
        
        n80 <- if (length(targets) >= 1L) {
          .solve_n_single_exact_uncond_mcnemar_power(
            pa = pa_hat,
            pb = pb_hat,
            pc = pc_hat,
            pd = pd_hat,
            target_power = targets[1],
            alpha = alpha,
            n_min = 2L,
            n_max = n_max,
            gamma = gamma,
            num_pi_values = num_pi_values,
            strict = strict,
            cache_env = cache_env,
            progress = TRUE,
            progress_every = 5L
          )
        } else {
          NA_integer_
        }
        
        n90 <- if (length(targets) >= 2L) {
          .solve_n_single_exact_uncond_mcnemar_power(
            pa = pa_hat,
            pb = pb_hat,
            pc = pc_hat,
            pd = pd_hat,
            target_power = targets[2],
            alpha = alpha,
            n_min = if (is.finite(n80)) n80 else 2L,
            n_max = n_max,
            gamma = gamma,
            num_pi_values = num_pi_values,
            strict = strict,
            cache_env = cache_env,
            progress = TRUE,
            progress_every = 2L
          )
        } else {
          NA_integer_
        }
        
      } else {
        pw <- .power_single_enum_mcnemar_power(
          pb = pb_hat,
          pc = pc_hat,
          n = n,
          method = mth,
          alpha = alpha,
          gamma = gamma,
          num_pi_values = num_pi_values,
          strict = strict,
          cache_env = cache_env
        )
        
        n80 <- if (length(targets) >= 1L) {
          .solve_n_single_mcnemar_power(
            pb = pb_hat,
            pc = pc_hat,
            target_power = targets[1],
            method = mth,
            alpha = alpha,
            n_min = 2L,
            n_max = n_max,
            gamma = gamma,
            num_pi_values = num_pi_values,
            strict = strict,
            cache_env = cache_env
          )
        } else {
          NA_integer_
        }
        
        n90 <- if (length(targets) >= 2L) {
          .solve_n_single_mcnemar_power(
            pb = pb_hat,
            pc = pc_hat,
            target_power = targets[2],
            method = mth,
            alpha = alpha,
            n_min = 2L,
            n_max = n_max,
            gamma = gamma,
            num_pi_values = num_pi_values,
            strict = strict,
            cache_env = cache_env
          )
        } else {
          NA_integer_
        }
      }
      
      out_rows[[row_idx]] <- tibble::tibble(
        control = control,
        treatment = trt,
        method_label = mth,
        a = a, b = b, c = c, d = d,
        n = n,
        pb_hat = pb_hat,
        pc_hat = pc_hat,
        PowerS = pw,
        N80S = n80,
        N90S = n90
      )
    }
  }
  
  results <- dplyr::bind_rows(out_rows)
  
  out <- list(
    results = results,
    settings = list(
      id = id,
      condition = condition,
      outcome = outcome,
      control = control,
      treatments = treatments,
      methods = methods,
      p_adjust = p_adjust,
      alpha = alpha,
      ci = ci,
      outcome_levels = outcome_levels,
      order = order,
      order_info = order_info,
      drop_na = drop_na,
      quiet = quiet,
      gamma = gamma,
      num_pi_values = num_pi_values,
      strict = strict,
      targets = targets,
      n_max = n_max
    )
  )
  
  class(out) <- "mcnemarPower"
  out
}

# ---- internal helpers ---------------------------------------------------------
.build_wide_for_mcnemar_power <- function(data, id, condition, outcome, control, treatments) {
  df <- data[as.character(data[[condition]]) %in% c(control, treatments), , drop = FALSE]
  key <- paste(df[[id]], df[[condition]], sep = "__")
  if (anyDuplicated(key)) {
    stop("Multiple rows per id-condition detected. Aggregate first or ensure uniqueness.")
  }
  
  tidyr::pivot_wider(
    df,
    id_cols = dplyr::all_of(id),
    names_from = dplyr::all_of(condition),
    values_from = dplyr::all_of(outcome)
  )
}

.build_table_2x2_mcnemar_power <- function(x, y, outcome_levels) {
  x <- factor(x, levels = outcome_levels)
  y <- factor(y, levels = outcome_levels)
  tab <- table(x, y)
  if (!all(dim(tab) == c(2L, 2L))) {
    tab2 <- matrix(0L, 2L, 2L, dimnames = list(outcome_levels, outcome_levels))
    tab2[rownames(tab), colnames(tab)] <- tab
    tab <- tab2
  }
  tab
}

.extract_order_structure_mcnemar_power <- function(data, id, condition, control, treatments, order = NULL) {
  keep <- as.character(data[[condition]]) %in% c(control, treatments)
  df <- data[keep, , drop = FALSE]
  
  if (!is.null(order)) {
    ord <- order(df[[id]], df[[order]], seq_len(nrow(df)))
  } else {
    ord <- order(df[[id]], seq_len(nrow(df)))
  }
  df <- df[ord, , drop = FALSE]
  
  split_seq <- split(as.character(df[[condition]]), df[[id]])
  split_seq <- lapply(split_seq, function(z) unname(z[!is.na(z)]))
  split_seq <- split_seq[lengths(split_seq) > 0]
  
  seq_strings <- vapply(split_seq, function(z) paste(z, collapse = " -> "), character(1))
  seq_tab <- sort(table(seq_strings), decreasing = TRUE)
  
  unique_sequences <- unique(unname(seq_strings))
  seq_lookup <- lapply(unique_sequences, function(ss) strsplit(ss, " -> ", fixed = TRUE)[[1]])
  names(seq_lookup) <- unique_sequences
  
  modal_sequence <- if (length(seq_tab) > 0L) names(seq_tab)[1] else NA_character_
  
  list(
    subject_sequences = split_seq,
    sequence_table = seq_tab,
    sequence_lookup = seq_lookup,
    modal_sequence = modal_sequence,
    unique_sequences = names(seq_tab)
  )
}

.pval_mcnemar_asymptotic_mcnemar_power <- function(b, c) {
  m <- b + c
  if (m == 0) return(1)
  z2 <- (b - c)^2 / m
  stats::pchisq(z2, df = 1, lower.tail = FALSE)
}

.pval_mcnemar_cc_mcnemar_power <- function(b, c) {
  m <- b + c
  if (m == 0) return(1)
  z2 <- (max(abs(b - c) - 1, 0))^2 / m
  stats::pchisq(z2, df = 1, lower.tail = FALSE)
}

.pval_mcnemar_exact_cond_mcnemar_power <- function(b, c) {
  m <- b + c
  if (m == 0) return(1)
  x <- min(b, c)
  p <- 2 * stats::pbinom(x, size = m, prob = 0.5)
  min(1, p)
}

.pval_mcnemar_midp_mcnemar_power <- function(b, c) {
  m <- b + c
  if (m == 0) return(1)
  x <- min(b, c)
  p_exact <- 2 * stats::pbinom(x, size = m, prob = 0.5)
  p_mid <- p_exact - stats::dbinom(x, size = m, prob = 0.5)
  min(1, max(0, p_mid))
}

.pval_mcnemar_exact_uncond_mcnemar_power <- function(a, b, c, d,
                                                     gamma = 1e-4,
                                                     num_pi_values = 1000L,
                                                     cache_env = NULL) {
  key <- paste0(
    "exact_uncond|", a, "|", b, "|", c, "|", d, "|",
    gamma, "|", as.integer(num_pi_values)
  )
  
  if (!is.null(cache_env) && exists(key, envir = cache_env, inherits = FALSE)) {
    return(get(key, envir = cache_env, inherits = FALSE))
  }
  
  vals <- c(a, b, c, d)
  if (!all(is.finite(vals)) || any(vals < 0)) {
    stop("`a`, `b`, `c`, and `d` must be non-negative finite counts.")
  }
  
  a <- as.integer(a)
  b <- as.integer(b)
  c <- as.integer(c)
  d <- as.integer(d)
  
  # Boundary case: no discordant pairs
  # McNemar-style null is trivially not rejected here, and the conditional helper
  # already returns 1 in the analogous case b + c == 0.
  if ((b + c) == 0L) {
    pval <- 1
    if (!is.null(cache_env)) assign(key, pval, envir = cache_env)
    return(pval)
  }
  
  tab <- matrix(c(a, b,
                  c, d), nrow = 2L, byrow = TRUE)
  
  # Some boundary full tables can make the contingencytables routine return
  # an internal NA test statistic. Catch those safely.
  res <- tryCatch(
    contingencytables::McNemar_exact_unconditional_test_paired_2x2(
      tab,
      gamma = gamma,
      num_pi_values = as.integer(num_pi_values)
    ),
    error = function(e) NULL
  )
  
  # If the exact unconditional routine fails on a boundary configuration,
  # treat the p-value as NA so the caller can skip that state cleanly.
  if (is.null(res)) {
    pval <- NA_real_
    if (!is.null(cache_env)) assign(key, pval, envir = cache_env)
    return(pval)
  }
  
  # Prefer documented top-level field
  pval <- NA_real_
  if (!is.null(res$P) && is.numeric(res$P)) {
    val <- as.numeric(res$P[1])
    if (is.finite(val) && val >= 0 && val <= 1) {
      pval <- val
    }
  }
  
  # Fallbacks for robustness
  if (!is.finite(pval)) {
    common <- c("P", "p.value", "p_value", "pvalue", "p", "pval", "p_val")
    for (nm in common) {
      if (!is.null(res[[nm]]) && is.numeric(res[[nm]])) {
        val <- as.numeric(res[[nm]][1])
        if (is.finite(val) && val >= 0 && val <= 1) {
          pval <- val
          break
        }
      }
    }
  }
  
  if (!is.finite(pval)) {
    flat <- unlist(res, recursive = TRUE, use.names = TRUE)
    if (length(flat) > 0L && !is.null(names(flat))) {
      flat_num <- suppressWarnings(as.numeric(flat))
      nm <- names(flat)
      
      idx <- which(
        is.finite(flat_num) &
          flat_num >= 0 & flat_num <= 1 &
          grepl("(^P$|p.value|p_value|pvalue|^p$|pval|p_val)", nm, ignore.case = TRUE)
      )
      
      if (length(idx) > 0L) {
        pval <- flat_num[idx[1]]
      } else {
        idx2 <- which(is.finite(flat_num) & flat_num >= 0 & flat_num <= 1)
        if (length(idx2) == 1L) {
          pval <- flat_num[idx2[1]]
        }
      }
    }
  }
  
  if (!is.null(cache_env)) assign(key, pval, envir = cache_env)
  pval
}

.pval_mcnemar_method_mcnemar_power <- function(a = NULL, b, c, d = NULL, method,
                                               gamma = 1e-4,
                                               num_pi_values = 1000L,
                                               cache_env = NULL) {
  switch(
    method,
    asymptotic = .pval_mcnemar_asymptotic_mcnemar_power(b, c),
    cc = .pval_mcnemar_cc_mcnemar_power(b, c),
    midp = .pval_mcnemar_midp_mcnemar_power(b, c),
    exact_cond = .pval_mcnemar_exact_cond_mcnemar_power(b, c),
    exact_uncond = .pval_mcnemar_exact_uncond_mcnemar_power(
      a = a, b = b, c = c, d = d,
      gamma = gamma,
      num_pi_values = num_pi_values,
      cache_env = cache_env
    ),
    stop("Unknown method: ", method)
  )
}

.power_single_enum_exact_uncond_full_mcnemar_power <- function(pa, pb, pc, pd, n,
                                                               alpha = 0.05,
                                                               gamma = 1e-4,
                                                               num_pi_values = 1000L,
                                                               strict = FALSE,
                                                               cache_env = NULL) {
  probs <- c(pa, pb, pc, pd)
  
  if (!is.numeric(probs) ||
      length(probs) != 4L ||
      any(is.na(probs)) ||
      any(probs < 0) ||
      abs(sum(probs) - 1) > 1e-12) {
    stop("`pa`, `pb`, `pc`, and `pd` must be non-negative probabilities summing to 1.")
  }
  
  n <- as.integer(n)
  if (n < 1L) stop("n must be >= 1.")
  
  lf  <- lfactorial(0:n)
  lpa <- if (pa == 0) -Inf else log(pa)
  lpb <- if (pb == 0) -Inf else log(pb)
  lpc <- if (pc == 0) -Inf else log(pc)
  lpd <- if (pd == 0) -Inf else log(pd)
  
  out <- 0
  
  for (a in 0:n) {
    max_bc <- n - a
    
    if (strict) {
      # full search
      for (b in 0:max_bc) {
        for (c in 0:(max_bc - b)) {
          d <- n - a - b - c
          
          if ((b + c) == 0L) next
          
          term_a <- if (a == 0L) 0 else a * lpa
          term_b <- if (b == 0L) 0 else b * lpb
          term_c <- if (c == 0L) 0 else c * lpc
          term_d <- if (d == 0L) 0 else d * lpd
          
          logpr <- lf[n + 1L] - lf[a + 1L] - lf[b + 1L] - lf[c + 1L] - lf[d + 1L] +
            term_a + term_b + term_c + term_d
          
          pr <- exp(logpr)
          if (!is.finite(pr) || pr == 0) next
          
          pval <- .pval_mcnemar_exact_uncond_mcnemar_power(
            a = a, b = b, c = c, d = d,
            gamma = gamma,
            num_pi_values = num_pi_values,
            cache_env = cache_env
          )
          
          if (!is.finite(pval)) next
          if (pval <= alpha) out <- out + pr
        }
      }
      
    } else {
      # exact-safe reduced search: only states with b >= c
      for (c in 0:max_bc) {
        b_min <- c
        b_max <- max_bc - c
        if (b_min > b_max) next
        
        for (b in b_min:b_max) {
          d <- n - a - b - c
          
          if ((b + c) == 0L) next
          
          term_a <- if (a == 0L) 0 else a * lpa
          term_b <- if (b == 0L) 0 else b * lpb
          term_c <- if (c == 0L) 0 else c * lpc
          term_d <- if (d == 0L) 0 else d * lpd
          
          logpr <- lf[n + 1L] - lf[a + 1L] - lf[b + 1L] - lf[c + 1L] - lf[d + 1L] +
            term_a + term_b + term_c + term_d
          
          pr <- exp(logpr)
          if (!is.finite(pr) || pr == 0) next
          
          pval <- .pval_mcnemar_exact_uncond_mcnemar_power(
            a = a, b = b, c = c, d = d,
            gamma = gamma,
            num_pi_values = num_pi_values,
            cache_env = cache_env
          )
          
          if (!is.finite(pval)) next
          if (pval <= alpha) out <- out + pr
        }
      }
    }
  }
  
  out
}

.power_single_enum_mcnemar_power <- function(pb, pc, n,
                                             method = c("asymptotic", "cc", "midp", "exact_cond"),
                                             alpha = 0.05,
                                             gamma = 1e-4,
                                             num_pi_values = 1000L,
                                             strict = FALSE,
                                             cache_env = NULL) {
  method <- match.arg(method)
  
  if (!is.numeric(pb) || !is.numeric(pc) ||
      length(pb) != 1L || length(pc) != 1L ||
      is.na(pb) || is.na(pc) || pb < 0 || pc < 0 || pb + pc > 1) {
    stop("pb and pc must be single probabilities with pb >= 0, pc >= 0, pb + pc <= 1.")
  }
  
  n <- as.integer(n)
  if (n < 1L) stop("n must be >= 1.")
  
  prest <- 1 - pb - pc
  
  # precompute exact-safe pieces
  lf <- lfactorial(0:n)
  lpb <- if (pb == 0) -Inf else log(pb)
  lpc <- if (pc == 0) -Inf else log(pc)
  lprst <- if (prest == 0) -Inf else log(prest)
  
  out <- 0
  
  for (b in 0:n) {
    for (c in 0:(n - b)) {
      rest <- n - b - c
      
      # exact-safe pruning
      if (!strict && b < c) next
      if ((b + c) == 0L) next
      
      term_b <- if (b == 0L) 0 else b * lpb
      term_c <- if (c == 0L) 0 else c * lpc
      term_r <- if (rest == 0L) 0 else rest * lprst
      
      logpr <- lf[n + 1L] - lf[b + 1L] - lf[c + 1L] - lf[rest + 1L] +
        term_b + term_c + term_r
      
      pr <- exp(logpr)
      if (!is.finite(pr) || pr == 0) next
      
      pval <- .pval_mcnemar_method_mcnemar_power(
        b = b,
        c = c,
        method = method,
        gamma = gamma,
        num_pi_values = num_pi_values,
        cache_env = cache_env
      )
      
      if (!is.finite(pval)) next
      
      if (pval <= alpha) {
        out <- out + pr
      }
    }
  }
  
  out
}

.solve_n_single_exact_uncond_mcnemar_power <- function(pa, pb, pc, pd,
                                                       target_power,
                                                       alpha = 0.05,
                                                       n_min = 2L,
                                                       n_max = 2000L,
                                                       gamma = 1e-4,
                                                       num_pi_values = 1000L,
                                                       strict = FALSE,
                                                       cache_env = NULL,
                                                       progress = FALSE,
                                                       progress_every = 1L,
                                                       show_notice = TRUE,
                                                       notice_fun = message) {
  
  target_power <- as.numeric(target_power)
  
  if (length(target_power) != 1L ||
      is.na(target_power) ||
      target_power <= 0 ||
      target_power >= 1) {
    stop("target_power must be a single number in (0,1).")
  }
  
  n_min <- as.integer(n_min)
  n_max <- as.integer(n_max)
  progress_every <- max(1L, as.integer(progress_every))
  
  pb_obj <- NULL
  t0 <- proc.time()[["elapsed"]]
  
  if (isTRUE(progress) && isTRUE(show_notice)) {
    notice_fun(
      paste0(
        "Computing sample-size targets for `exact_uncond` may take noticeably longer than other McNemar methods. ",
        "Unlike the asymptotic, continuity-corrected, mid-P, and exact-conditional methods, ",
        "A progress bar is shown for this method because those repeated evaluations can take time."
      )
    )
  }
  
  if (isTRUE(progress)) {
    pb_obj <- utils::txtProgressBar(min = n_min, max = n_max, style = 3)
    on.exit(close(pb_obj), add = TRUE)
  }
  
  for (n in seq.int(n_min, n_max)) {
    pw <- .power_single_enum_exact_uncond_full_mcnemar_power(
      pa = pa, pb = pb, pc = pc, pd = pd, n = n,
      alpha = alpha,
      gamma = gamma,
      num_pi_values = num_pi_values,
      strict = strict,
      cache_env = cache_env
    )
    
    if (isTRUE(progress)) {
      utils::setTxtProgressBar(pb_obj, n)
      if (((n - n_min + 1L) %% progress_every) == 0L) {
        elapsed <- proc.time()[["elapsed"]] - t0
        message(sprintf(
          " exact_uncond target=%.2f | n=%d/%d | current power=%.4f | elapsed=%.1fs",
          target_power, n, n_max, pw, elapsed
        ))
      }
    }
    
    if (is.finite(pw) && pw >= target_power) {
      return(as.integer(n))
    }
  }
  
  NA_integer_
}


.solve_n_single_mcnemar_power <- function(pb, pc,
                                          target_power,
                                          method,
                                          alpha = 0.05,
                                          n_min = 2L,
                                          n_max = 2000L,
                                          gamma = 1e-4,
                                          num_pi_values = 1000L,
                                          strict = FALSE,
                                          cache_env = NULL) {
  if (identical(method, "exact_uncond")) {
    stop("Use .solve_n_single_exact_uncond_mcnemar_power() for method = 'exact_uncond'.")
  }
  target_power <- as.numeric(target_power)
  if (length(target_power) != 1L || is.na(target_power) || target_power <= 0 || target_power >= 1) {
    stop("target_power must be a single number in (0,1).")
  }
  
  for (n in seq.int(as.integer(n_min), as.integer(n_max))) {
    pw <- .power_single_enum_mcnemar_power(
      pb = pb, pc = pc, n = n,
      method = method,
      alpha = alpha,
      gamma = gamma,
      num_pi_values = num_pi_values,
      strict = strict,
      cache_env = cache_env
    )
    if (is.finite(pw) && pw >= target_power) return(as.integer(n))
  }
  NA_integer_
}

.q_discrete_mcnemar_power <- function(t, m) {
  if (!is.finite(t) || t <= 0) return(0)
  if (t >= 1) return(1)
  xs <- 0:m
  pvals <- vapply(xs, function(x) .pval_mcnemar_exact_cond_mcnemar_power(x, m - x), numeric(1))
  w <- stats::dbinom(xs, size = m, prob = 0.5)
  sum(w[pvals <= t])
}

.discrete_holm_mcnemar_power <- function(p, m_disc) {
  stopifnot(length(p) == length(m_disc))
  n <- length(p)
  out <- rep(NA_real_, n)
  ok <- which(is.finite(p) & is.finite(m_disc))
  if (length(ok) == 0L) return(out)
  
  ord <- order(p[ok], na.last = TRUE)
  p_ord <- p[ok][ord]
  m_ord <- as.integer(m_disc[ok][ord])
  
  s <- rep(NA_real_, length(p_ord))
  for (j in seq_along(p_ord)) {
    qsum <- 0
    for (k in j:length(p_ord)) {
      qsum <- qsum + .q_discrete_mcnemar_power(p_ord[j], m_ord[k])
    }
    s[j] <- min(qsum, 1)
  }
  
  adj <- cummax(s)
  inv_ord <- integer(length(ord))
  inv_ord[ord] <- seq_along(ord)
  out_vec <- rep(NA_real_, length(ok))
  out_vec[inv_ord] <- adj
  out[ok] <- pmin(out_vec, 1)
  out
}

.adjust_p_family_mcnemar_power <- function(p, m_disc, p_adjust) {
  if (p_adjust == "discrete_holm") {
    .discrete_holm_mcnemar_power(p, m_disc)
  } else {
    stats::p.adjust(p, method = p_adjust)
  }
}

.sim_one_replicate_indep_mcnemar_power <- function(pair_info,
                                                   method,
                                                   p_adjust,
                                                   alpha,
                                                   gamma,
                                                   num_pi_values,
                                                   cache_env = NULL) {
  pvals <- numeric(nrow(pair_info))
  disc <- numeric(nrow(pair_info))
  
  for (j in seq_len(nrow(pair_info))) {
    n <- pair_info$n[j]
    pb <- pair_info$pb_hat[j]
    pc <- pair_info$pc_hat[j]
    draw <- as.vector(stats::rmultinom(1, size = n, prob = c(pb, pc, 1 - pb - pc)))
    b <- draw[1]
    c <- draw[2]
    pvals[j] <- .pval_mcnemar_method_mcnemar_power(
      b = b, c = c,
      method = method,
      gamma = gamma,
      num_pi_values = num_pi_values,
      cache_env = cache_env
    )
    disc[j] <- b + c
  }
  
  padj <- .adjust_p_family_mcnemar_power(pvals, disc, p_adjust)
  list(p = pvals, padj = padj, disc = disc)
}

.sim_one_replicate_dep_rand_mcnemar_power <- function(wide_data,
                                                      order_info,
                                                      control,
                                                      treatments,
                                                      method,
                                                      p_adjust,
                                                      alpha,
                                                      gamma,
                                                      num_pi_values,
                                                      rho = 0.30,
                                                      cache_env = NULL) {
  sigma_u <- sqrt(rho / max(1e-8, 1 - rho))
  all_cond <- c(control, treatments)
  
  marg <- stats::setNames(rep(NA_real_, length(all_cond)), all_cond)
  for (nm in all_cond) {
    z <- wide_data[[nm]]
    if (all(is.na(z))) {
      marg[[nm]] <- 0.5
    } else if (is.numeric(z) || is.logical(z)) {
      marg[[nm]] <- mean(as.numeric(z) > 0, na.rm = TRUE)
    } else {
      lv <- levels(as.factor(z))
      marg[[nm]] <- mean(as.character(z) == lv[min(2, length(lv))], na.rm = TRUE)
    }
  }
  
  nsub <- nrow(wide_data)
  u <- stats::rnorm(nsub, mean = 0, sd = sigma_u)
  
  if (length(order_info$sequence_table) > 0L) {
    seq_names <- names(order_info$sequence_table)
    seq_prob <- as.numeric(order_info$sequence_table) / sum(order_info$sequence_table)
    draw_seq <- sample(seq_names, size = nsub, replace = TRUE, prob = seq_prob)
    seq_list <- lapply(draw_seq, function(ss) order_info$sequence_lookup[[ss]])
  } else {
    seq_list <- replicate(nsub, sample(all_cond, length(all_cond), replace = FALSE), simplify = FALSE)
  }
  invisible(seq_list)
  
  sim <- vector("list", length(all_cond))
  names(sim) <- all_cond
  for (nm in all_cond) {
    eta <- stats::qlogis(pmin(pmax(marg[[nm]], 1e-6), 1 - 1e-6)) + u
    pr <- stats::plogis(eta)
    sim[[nm]] <- stats::rbinom(nsub, size = 1, prob = pr)
  }
  
  pvals <- numeric(length(treatments))
  disc <- numeric(length(treatments))
  for (j in seq_along(treatments)) {
    trt <- treatments[j]
    x <- sim[[control]]
    y <- sim[[trt]]
    b <- sum(x == 0 & y == 1)
    c <- sum(x == 1 & y == 0)
    pvals[j] <- .pval_mcnemar_method_mcnemar_power(
      b = b, c = c,
      method = method,
      gamma = gamma,
      num_pi_values = num_pi_values,
      cache_env = cache_env
    )
    disc[j] <- b + c
  }
  
  padj <- .adjust_p_family_mcnemar_power(pvals, disc, p_adjust)
  list(p = pvals, padj = padj, disc = disc)
}

.sim_one_replicate_dep_fixed_mcnemar_power <- function(wide_data,
                                                       control,
                                                       treatments,
                                                       method,
                                                       p_adjust,
                                                       alpha,
                                                       gamma,
                                                       num_pi_values,
                                                       rho = 0.30,
                                                       order_effect = 0.25,
                                                       fixed_order = NULL,
                                                       cache_env = NULL) {
  sigma_u <- sqrt(rho / max(1e-8, 1 - rho))
  if (is.null(fixed_order)) fixed_order <- c(control, treatments)
  fixed_order <- as.character(fixed_order)
  
  all_cond <- unique(c(control, treatments))
  marg <- stats::setNames(rep(NA_real_, length(all_cond)), all_cond)
  for (nm in all_cond) {
    z <- wide_data[[nm]]
    if (all(is.na(z))) {
      marg[[nm]] <- 0.5
    } else if (is.numeric(z) || is.logical(z)) {
      marg[[nm]] <- mean(as.numeric(z) > 0, na.rm = TRUE)
    } else {
      lv <- levels(as.factor(z))
      marg[[nm]] <- mean(as.character(z) == lv[min(2, length(lv))], na.rm = TRUE)
    }
  }
  
  nsub <- nrow(wide_data)
  u <- stats::rnorm(nsub, mean = 0, sd = sigma_u)
  sim <- vector("list", length(all_cond))
  names(sim) <- all_cond
  
  pos_map <- stats::setNames(seq_along(fixed_order), fixed_order)
  pos_map <- pos_map[names(pos_map) %in% all_cond]
  
  for (nm in all_cond) {
    jj <- if (nm %in% names(pos_map)) pos_map[[nm]] else match(nm, all_cond)
    base <- stats::qlogis(pmin(pmax(marg[[nm]], 1e-6), 1 - 1e-6))
    eta <- base + u + order_effect * (jj - 1)
    pr <- stats::plogis(eta)
    sim[[nm]] <- stats::rbinom(nsub, size = 1, prob = pr)
  }
  
  pvals <- numeric(length(treatments))
  disc <- numeric(length(treatments))
  for (j in seq_along(treatments)) {
    trt <- treatments[j]
    x <- sim[[control]]
    y <- sim[[trt]]
    b <- sum(x == 0 & y == 1)
    c <- sum(x == 1 & y == 0)
    pvals[j] <- .pval_mcnemar_method_mcnemar_power(
      b = b, c = c,
      method = method,
      gamma = gamma,
      num_pi_values = num_pi_values,
      cache_env = cache_env
    )
    disc[j] <- b + c
  }
  
  padj <- .adjust_p_family_mcnemar_power(pvals, disc, p_adjust)
  list(p = pvals, padj = padj, disc = disc)
}

#' Summary for `mcnemar_power()` results
#'
#' @description
#' Produces a structured summary of an object returned by `mcnemar_power()`,
#' including multiplicity-adjusted power estimates and scenario-based sample-size
#' predictions under independence, dependence with randomized/counterbalanced
#' order, and dependence with fixed order.
#'
#' @param object An object returned by `mcnemar_power()`.
#' @param data The original long-format data used for `mcnemar_power()`.
#' @param comparisons Optional character vector of treatment levels defining the
#'   subset for subset-power summaries. Default: all fitted treatments.
#' @param scenarios Character vector selecting one or more of `"indep"`,
#'   `"dep_rand"`, and `"dep_fixed"`.
#' @param alpha Optional significance threshold for adjusted-power summaries.
#'   Default is the `alpha` stored in the fit.
#' @param nsim Number of Monte Carlo simulation replicates per scenario and method.
#' @param rho Dependence parameter for the built-in latent random-intercept model.
#' @param order_effect Size of the additive fixed-order effect applied by the
#'   `"dep_fixed"` scenario.
#' @param fixed_order Optional character vector specifying the fixed treatment order.
#' @param targets Numeric vector of target powers for scenario-based sample-size prediction.
#' @param power_types Character vector selecting one or more sample-size targets
#'   among `"individual"`, `"any_subset"`, `"all_subset"`, and `"all_treatments"`.
#' @param individual_treatment Treatment label(s) used when `"individual"`
#'   sample size is requested.
#' @param sim_fun Optional user-supplied replicate simulator.
#' @param parallel Logical; if `TRUE`, enable parallel replicate generation where supported.
#' @param n_cores Number of cores used if `parallel = TRUE`.
#' @param seed Optional random seed for reproducibility.
#' @param n_max Maximum sample size searched for scenario-based sample-size estimation.
#' @param digits Number of digits used in the print method.
#' @param x An object of class \code{"summary.mcnemarPower"} to be printed.
#' @param ... Additional arguments passed to `sim_fun`, if supplied.
#'
#' @return
#' An object of class \code{"summary.mcnemarPower"}.
#'
#' The returned object contains:
#' \itemize{
#'   \item \code{single_test}: The original comparison-level results table from
#'   the \code{"mcnemarPower"} fit.
#'
#'   \item \code{adjusted_power}: A nested list of multiplicity-adjusted power
#'   summaries by McNemar method and scenario.
#'
#'   \item \code{settings}: The original settings stored in the fitted object.
#'
#'   \item \code{comparisons}: The subset of treatment labels used for
#'   subset-based family power summaries.
#'
#'   \item \code{individual_treatment}: The treatment label(s) used for
#'   \code{"individual"} sample-size summaries.
#'
#'   \item \code{scenarios}: The evaluated scenario(s), such as
#'   \code{"indep"}, \code{"dep_rand"}, and \code{"dep_fixed"}.
#'
#'   \item \code{nsim}: The number of Monte Carlo replicates used in the
#'   summary simulation.
#' }
#'
#' For each McNemar method and scenario, the \code{adjusted_power} component
#' contains:
#' \itemize{
#'   \item \code{individual_adjusted}: A tibble with one row per fitted treatment
#'   and the corresponding individual adjusted power after multiplicity adjustment.
#'
#'   \item \code{subset_any}: The estimated probability that at least one treatment
#'   in \code{comparisons} is significant after adjustment.
#'
#'   \item \code{subset_all}: The estimated probability that all treatments in
#'   \code{comparisons} are significant after adjustment.
#'
#'   \item \code{all_treatments}: The estimated probability that all fitted
#'   treatments are significant after adjustment.
#'
#'   \item \code{sample_size}: A list containing up to four tibbles for
#'   scenario-based sample-size summaries:
#'   \code{individual}, \code{any_subset}, \code{all_subset}, and
#'   \code{all_treatments}.
#' }
#'
#' If \code{comparisons = c("TrtA", "TrtB")}, then:
#' \itemize{
#'   \item \code{subset_any}: The probability that at least one of TrtA or TrtB
#'   is significant.
#'
#'   \item \code{subset_all}: The probability that both TrtA and TrtB are
#'   significant.
#'
#'   \item \code{all_treatments}: The probability that all fitted treatments
#'   (for example TrtA, TrtB, and TrtC) are significant.
#' }
#'
#' \strong{Interpreting \code{individual_adjusted}:}
#'
#' The \code{individual_adjusted} tibble reports adjusted power separately for
#' each fitted treatment, even when \code{comparisons} defines only a subset
#' for the subset-based family power summaries.
#' 
#' \strong{Scenario definitions:}
#'
#' The \code{scenarios} argument controls how multiplicity-adjusted power is
#' simulated in \code{summary.mcnemarPower()}.
#'
#' \itemize{
#' \item \code{"indep"} treats the control-versus-treatment McNemar comparisons
#' as independent. For each treatment, discordant counts are simulated from the
#' estimated paired-binomial model without adding subject-level dependence across
#' treatments.
#'
#' \item \code{"dep_rand"} represents a repeated-measures setting in which
#' outcomes from the same subject are dependent, but condition order is
#' randomized or counterbalanced. The simulation uses a shared subject-level
#' random effect to induce dependence and samples the condition sequence from the
#' observed subject-level order structure in the supplied long-format data.
#'
#' \item \code{"dep_fixed"} represents a repeated-measures setting in which the
#' control is presented first and the treatments follow a fixed sequence. The
#' simulation uses a shared subject-level random effect to induce dependence and
#' adds a position-specific \code{order_effect} across the fixed order.
#' }
#'
#' The \code{"dep_rand"} and \code{"dep_fixed"} scenarios are
#' \strong{sensitivity analyses}. This means they do not attempt to estimate a
#' repeated-binary model directly from the observed data. Instead, they examine
#' how the multiplicity-adjusted power results may change when additional
#' assumptions are introduced about dependence among repeated outcomes from the
#' same subject and, where relevant, about the order in which conditions are
#' presented.
#'
#' In \code{"dep_rand"}, the analysis assumes that outcomes measured on the same
#' subject are correlated, but that the condition order is randomized or
#' counterbalanced. In \code{"dep_fixed"}, the analysis assumes both
#' within-subject dependence and a fixed presentation order, so that later
#' conditions may differ systematically from earlier ones through the
#' \code{order_effect} parameter. These scenarios therefore provide a
#' simulation-based assessment of robustness: they show whether the power
#' conclusions remain similar, become weaker, or become stronger when plausible
#' dependence and order structures are imposed.
#'
#' Because these scenarios are not fitted repeated-binary models, they should not
#' be interpreted as direct estimates of period effects, sequence effects,
#' carry-over, or subject-level random effects from the observed study. Their role
#' is instead to help users judge how sensitive the power calculations are to
#' assumptions about the repeated-measures structure of the data.
#'
#' @seealso \code{mcnemar_power_summary}, \code{mcnemar_power}
#'
#' @export
summary.mcnemarPower <- function(object,
                                 data,
                                 comparisons = NULL,
                                 scenarios = c("indep", "dep_rand", "dep_fixed"),
                                 alpha = NULL,
                                 nsim = 2000,
                                 rho = 0.30,
                                 order_effect = 0.25,
                                 fixed_order = NULL,
                                 targets = c(0.80, 0.90),
                                 power_types = c("individual", "any_subset", "all_subset", "all_treatments"),
                                 individual_treatment = NULL,
                                 sim_fun = NULL,
                                 parallel = FALSE,
                                 n_cores = 1L,
                                 seed = NULL,
                                 n_max = 1000L,
                                 digits = 4,
                                 ...) {
  
  stopifnot(inherits(object, "mcnemarPower"))
  stopifnot(is.data.frame(data))
  if (!is.null(seed)) set.seed(seed)
  
  s <- object$settings
  res <- object$results
  
  if (is.null(alpha)) alpha <- s$alpha
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number in (0,1).")
  }
  if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) || nsim < 10) {
    stop("`nsim` must be a single integer >= 10.")
  }
  nsim <- as.integer(nsim)
  
  scenarios <- unique(as.character(scenarios))
  valid_scen <- c("indep", "dep_rand", "dep_fixed")
  if (!all(scenarios %in% valid_scen)) stop("All `scenarios` must be in: ", paste(valid_scen, collapse = ", "))
  
  power_types <- unique(as.character(power_types))
  valid_pt <- c("individual", "any_subset", "all_subset", "all_treatments")
  if (!all(power_types %in% valid_pt)) stop("All `power_types` must be in: ", paste(valid_pt, collapse = ", "))
  
  if (!is.numeric(targets) || any(is.na(targets)) || any(targets <= 0 | targets >= 1)) {
    stop("`targets` must contain one or more numbers in (0,1).")
  }
  if (length(n_max) != 1L || is.na(n_max) || n_max < 2L) stop("`n_max` must be a single integer >= 2.")
  n_max <- as.integer(n_max)
  
  if (!is.numeric(rho) || length(rho) != 1L || is.na(rho) || rho < 0 || rho >= 1) {
    stop("`rho` must be a single number in [0,1).")
  }
  if (!is.numeric(order_effect) || length(order_effect) != 1L || is.na(order_effect)) {
    stop("`order_effect` must be a single numeric value.")
  }
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) stop("`parallel` must be TRUE or FALSE.")
  if (length(n_cores) != 1L || is.na(n_cores) || n_cores < 1L) stop("`n_cores` must be a single integer >= 1.")
  n_cores <- as.integer(n_cores)
  if (!is.null(sim_fun) && !is.function(sim_fun)) stop("`sim_fun` must be NULL or a function.")
  
  wide <- .build_wide_for_mcnemar_power(
    data = data,
    id = s$id,
    condition = s$condition,
    outcome = s$outcome,
    control = s$control,
    treatments = s$treatments
  )
  
  order_info <- .extract_order_structure_mcnemar_power(
    data = data,
    id = s$id,
    condition = s$condition,
    control = s$control,
    treatments = s$treatments,
    order = s$order
  )
  
  if (is.null(comparisons)) comparisons <- s$treatments
  comparisons <- intersect(as.character(comparisons), s$treatments)
  if (length(comparisons) < 1L) stop("No valid `comparisons` remained after matching against fitted treatments.")
  
  if (is.null(individual_treatment)) {
    individual_treatment <- comparisons
  } else {
    individual_treatment <- intersect(as.character(individual_treatment), s$treatments)
  }
  
  if (length(individual_treatment) < 1L) {
    stop("`individual_treatment` must contain one or more fitted treatment levels.")
  }
  
  if (is.null(fixed_order)) {
    fixed_order <- if (isTRUE(nzchar(order_info$modal_sequence))) {
      order_info$sequence_lookup[[order_info$modal_sequence]]
    } else {
      c(s$control, s$treatments)
    }
  }
  fixed_order <- as.character(fixed_order)
  
  methods <- unique(res$method_label)
  summary_by_method <- list()
  cache_env <- new.env(parent = emptyenv())
  
  .run_reps <- function(FUN, nsim, parallel = FALSE, n_cores = 1L) {
    idx <- seq_len(nsim)
    if (isTRUE(parallel) && n_cores > 1L && .Platform$OS.type != "windows") {
      parallel::mclapply(idx, FUN, mc.cores = n_cores)
    } else {
      lapply(idx, FUN)
    }
  }
  
  for (mth in methods) {
    pair_info <- res[res$method_label == mth, , drop = FALSE]
    pair_info <- pair_info[match(s$treatments, pair_info$treatment), , drop = FALSE]
    
    one_rep_builtin <- function(scenario, wide_data_local, order_info_local) {
      if (scenario == "indep") {
        return(.sim_one_replicate_indep_mcnemar_power(
          pair_info = pair_info,
          method = mth,
          p_adjust = s$p_adjust,
          alpha = alpha,
          gamma = s$gamma,
          num_pi_values = s$num_pi_values,
          cache_env = cache_env
        ))
      }
      if (scenario == "dep_rand") {
        return(.sim_one_replicate_dep_rand_mcnemar_power(
          wide_data = wide_data_local,
          order_info = order_info_local,
          control = s$control,
          treatments = s$treatments,
          method = mth,
          p_adjust = s$p_adjust,
          alpha = alpha,
          gamma = s$gamma,
          num_pi_values = s$num_pi_values,
          rho = rho,
          cache_env = cache_env
        ))
      }
      if (scenario == "dep_fixed") {
        return(.sim_one_replicate_dep_fixed_mcnemar_power(
          wide_data = wide_data_local,
          control = s$control,
          treatments = s$treatments,
          method = mth,
          p_adjust = s$p_adjust,
          alpha = alpha,
          gamma = s$gamma,
          num_pi_values = s$num_pi_values,
          rho = rho,
          order_effect = order_effect,
          fixed_order = fixed_order,
          cache_env = cache_env
        ))
      }
      stop("Unknown scenario: ", scenario)
    }
    
    scen_out <- list()
    
    for (scen in scenarios) {
      rep_fun <- function(i) {
        if (is.null(sim_fun)) {
          one_rep_builtin(scen, wide, order_info)
        } else {
          sim_fun(
            scenario = scen,
            method = mth,
            wide_data = wide,
            settings = s,
            comparisons = comparisons,
            alpha = alpha,
            rho = rho,
            order_effect = order_effect,
            fixed_order = fixed_order,
            order_info = order_info,
            ...
          )
        }
      }
      
      reps <- .run_reps(rep_fun, nsim = nsim, parallel = parallel, n_cores = n_cores)
      sig_mat <- matrix(FALSE, nrow = nsim, ncol = length(s$treatments), dimnames = list(NULL, s$treatments))
      disc_mat <- matrix(NA_real_, nrow = nsim, ncol = length(s$treatments), dimnames = list(NULL, s$treatments))
      
      for (r in seq_len(nsim)) {
        rr <- reps[[r]]
        if (is.null(rr$padj)) stop("Each simulated replicate must return a `padj` component.")
        padj <- as.numeric(rr$padj)
        if (length(padj) != length(s$treatments)) {
          stop("Each simulated replicate must return `padj` with one value per fitted treatment.")
        }
        sig_mat[r, ] <- padj < alpha
        if (!is.null(rr$disc) && length(rr$disc) == length(s$treatments)) disc_mat[r, ] <- rr$disc
      }
      
      indiv_power <- colMeans(sig_mat)
      subset_idx <- which(colnames(sig_mat) %in% comparisons)
      indiv_idx <- which(colnames(sig_mat) %in% individual_treatment)
      
      any_subset <- mean(rowSums(sig_mat[, subset_idx, drop = FALSE]) >= 1)
      all_subset <- mean(rowSums(sig_mat[, subset_idx, drop = FALSE]) == length(subset_idx))
      all_treat <- mean(rowSums(sig_mat) == ncol(sig_mat))
      
      indiv_tbl <- tibble::tibble(
        scenario = scen,
        treatment = s$treatments,
        PowerAdj = indiv_power
      )
      
      estimate_n_target <- function(target,
                                    type = c("individual", "any_subset", "all_subset", "all_treatments"),
                                    treatment_name = NULL) {
        type <- match.arg(type)
        
        for (ntry in seq.int(2L, n_max)) {
          rep_fun_n <- function(i) {
            idx <- sample(seq_len(nrow(wide)), size = ntry, replace = TRUE)
            wide_n <- wide[idx, , drop = FALSE]
            order_n <- order_info
            
            if (length(order_info$subject_sequences) > 0L) {
              ids_all <- names(order_info$subject_sequences)
              draw_ids <- sample(ids_all, size = ntry, replace = TRUE)
              order_n$subject_sequences <- order_info$subject_sequences[draw_ids]
              if (length(order_n$subject_sequences) > 0L) {
                seq_strings_n <- vapply(order_n$subject_sequences, function(z) paste(z, collapse = " -> "), character(1))
                order_n$sequence_table <- sort(table(seq_strings_n), decreasing = TRUE)
                order_n$unique_sequences <- names(order_n$sequence_table)
              }
            }
            
            if (is.null(sim_fun)) {
              switch(
                scen,
                indep = .sim_one_replicate_indep_mcnemar_power(
                  pair_info = transform(pair_info, n = ntry),
                  method = mth,
                  p_adjust = s$p_adjust,
                  alpha = alpha,
                  gamma = s$gamma,
                  num_pi_values = s$num_pi_values,
                  cache_env = cache_env
                ),
                dep_rand = .sim_one_replicate_dep_rand_mcnemar_power(
                  wide_data = wide_n,
                  order_info = order_n,
                  control = s$control,
                  treatments = s$treatments,
                  method = mth,
                  p_adjust = s$p_adjust,
                  alpha = alpha,
                  gamma = s$gamma,
                  num_pi_values = s$num_pi_values,
                  rho = rho,
                  cache_env = cache_env
                ),
                dep_fixed = .sim_one_replicate_dep_fixed_mcnemar_power(
                  wide_data = wide_n,
                  control = s$control,
                  treatments = s$treatments,
                  method = mth,
                  p_adjust = s$p_adjust,
                  alpha = alpha,
                  gamma = s$gamma,
                  num_pi_values = s$num_pi_values,
                  rho = rho,
                  order_effect = order_effect,
                  fixed_order = fixed_order,
                  cache_env = cache_env
                ),
                stop("Unknown scenario: ", scen)
              )
            } else {
              sim_fun(
                scenario = scen,
                method = mth,
                wide_data = wide_n,
                settings = s,
                comparisons = comparisons,
                alpha = alpha,
                rho = rho,
                order_effect = order_effect,
                fixed_order = fixed_order,
                order_info = order_n,
                ...
              )
            }
          }
          
          reps_n <- .run_reps(rep_fun_n, nsim = nsim, parallel = parallel, n_cores = n_cores)
          sig_try <- matrix(FALSE, nrow = nsim, ncol = length(s$treatments),
                            dimnames = list(NULL, s$treatments))
          
          for (r in seq_len(nsim)) {
            rr <- reps_n[[r]]
            sig_try[r, ] <- as.numeric(rr$padj) < alpha
          }
          
          achieved <- switch(
            type,
            individual = {
              if (is.null(treatment_name)) stop("`treatment_name` must be supplied when type = 'individual'.")
              this_idx <- which(colnames(sig_try) == treatment_name)
              mean(sig_try[, this_idx])
            },
            any_subset = mean(rowSums(sig_try[, subset_idx, drop = FALSE]) >= 1),
            all_subset = mean(rowSums(sig_try[, subset_idx, drop = FALSE]) == length(subset_idx)),
            all_treatments = mean(rowSums(sig_try) == ncol(sig_try))
          )
          
          if (is.finite(achieved) && achieved >= target) return(as.integer(ntry))
        }
        
        NA_integer_
      }
      
      # ----------------------------------------------------------
      # Build separate sample-size tibbles for each power type
      # ----------------------------------------------------------
      
      # 1. individual
      individual_rows <- list()
      kk_ind <- 0L
      if ("individual" %in% power_types) {
        for (tt in targets) {
          for (trt_ind in individual_treatment) {
            kk_ind <- kk_ind + 1L
            individual_rows[[kk_ind]] <- tibble::tibble(
              scenario = scen,
              target_power = tt,
              treatment = trt_ind,
              sample_size = estimate_n_target(
                target = tt,
                type = "individual",
                treatment_name = trt_ind
              )
            )
          }
        }
      }
      individual_tbl <- if (length(individual_rows) > 0L) {
        dplyr::bind_rows(individual_rows)
      } else {
        tibble::tibble(
          scenario = character(0),
          target_power = numeric(0),
          treatment = character(0),
          sample_size = integer(0)
        )
      }
      
      # 2. any_subset
      any_subset_rows <- list()
      kk_any <- 0L
      if ("any_subset" %in% power_types) {
        for (tt in targets) {
          kk_any <- kk_any + 1L
          any_subset_rows[[kk_any]] <- tibble::tibble(
            scenario = scen,
            target_power = tt,
            treatment = NA_character_,
            sample_size = estimate_n_target(
              target = tt,
              type = "any_subset"
            )
          )
        }
      }
      any_subset_tbl <- if (length(any_subset_rows) > 0L) {
        dplyr::bind_rows(any_subset_rows)
      } else {
        tibble::tibble(
          scenario = character(0),
          target_power = numeric(0),
          treatment = character(0),
          sample_size = integer(0)
        )
      }
      
      # 3. all_subset
      all_subset_rows <- list()
      kk_allsub <- 0L
      if ("all_subset" %in% power_types) {
        for (tt in targets) {
          kk_allsub <- kk_allsub + 1L
          all_subset_rows[[kk_allsub]] <- tibble::tibble(
            scenario = scen,
            target_power = tt,
            treatment = NA_character_,
            sample_size = estimate_n_target(
              target = tt,
              type = "all_subset"
            )
          )
        }
      }
      all_subset_tbl <- if (length(all_subset_rows) > 0L) {
        dplyr::bind_rows(all_subset_rows)
      } else {
        tibble::tibble(
          scenario = character(0),
          target_power = numeric(0),
          treatment = character(0),
          sample_size = integer(0)
        )
      }
      
      # 4. all_treatments
      all_treatments_rows <- list()
      kk_alltrt <- 0L
      if ("all_treatments" %in% power_types) {
        for (tt in targets) {
          kk_alltrt <- kk_alltrt + 1L
          all_treatments_rows[[kk_alltrt]] <- tibble::tibble(
            scenario = scen,
            target_power = tt,
            treatment = NA_character_,
            sample_size = estimate_n_target(
              target = tt,
              type = "all_treatments"
            )
          )
        }
      }
      all_treatments_tbl <- if (length(all_treatments_rows) > 0L) {
        dplyr::bind_rows(all_treatments_rows)
      } else {
        tibble::tibble(
          scenario = character(0),
          target_power = numeric(0),
          treatment = character(0),
          sample_size = integer(0)
        )
      }
      
      # ----------------------------------------------------------
      # Store separate tibbles in the scenario output
      # ----------------------------------------------------------
      scen_out[[scen]] <- list(
        individual_adjusted = indiv_tbl,
        subset_any = any_subset,
        subset_all = all_subset,
        all_treatments = all_treat,
        sample_size = list(
          individual = individual_tbl,
          any_subset = any_subset_tbl,
          all_subset = all_subset_tbl,
          all_treatments = all_treatments_tbl
        ),
        mean_discordant = if (all(is.na(disc_mat))) {
          rep(NA_real_, ncol(disc_mat))
        } else {
          colMeans(disc_mat, na.rm = TRUE)
        }
      )
    }
    
    summary_by_method[[mth]] <- scen_out
  }
  
  out <- list(
    single_test = object$results,
    adjusted_power = summary_by_method,
    settings = object$settings,
    comparisons = comparisons,
    individual_treatment = individual_treatment,
    scenarios = scenarios,
    nsim = nsim,
    rho = rho,
    order_effect = order_effect,
    fixed_order = fixed_order,
    digits = as.integer(digits)
  )
  class(out) <- "summary.mcnemarPower"
  out
}

#' Convenience wrapper for `summary.mcnemarPower()`
#'
#' @description
#' Returns the same object as `summary(x, ...)` for an object produced by
#' `mcnemar_power()`. This wrapper mirrors the API style of
#' `mcnemar_control_summary()` and can be used when an explicit summary
#' function name is preferred over S3 dispatch.
#'
#' @param object An object returned by `mcnemar_power()`.
#' @param ... Additional arguments passed to `summary.mcnemarPower()`.
#'
#' @return
#' An object of class \code{"summary.mcnemarPower"}.
#'
#' This is a convenience wrapper around \code{summary()} for objects returned by
#' \code{mcnemar_power()}. It returns the same summary object as
#' \code{summary(object, ...)}.
#'
#' @details
#' Use \code{mcnemar_power_summary()} when an explicit summary function name is
#' preferred over S3 dispatch. The returned object contains the same elements as
#' \code{summary.mcnemarPower()}, including the single-test results table and the
#' scenario-based multiplicity-adjusted power summaries.
#'
#' @seealso \code{summary.mcnemarPower}, \code{mcnemar_control_summary}
#'
#' @export
mcnemar_power_summary <- function(object, ...) {
  summary(object, ...)
}

#' @rdname mcnemar_power
#' @export
print.mcnemarPower <- function(x, digits = 4, ...) {
  cat("mcnemarPower results\n")
  s <- x$settings
  cat(sprintf(" Control: %s\n", s$control))
  cat(sprintf(" Methods: %s\n", paste(s$methods, collapse = ", ")))
  cat(sprintf(" alpha: %.3f\n", s$alpha))
  cat(sprintf(" p-adjust (used in summary): %s\n", s$p_adjust))
  cat(sprintf(" n_max for single-test sample size: %d\n", s$n_max))
  cat(" Note: If NA is shown in N80S and/or N90S for a treatment, the target power\n")
  cat("       was not reached by the chosen n_max. Consider increasing n_max.\n")
  cat(sprintf(" Comparisons x methods: %d\n\n", nrow(x$results)))
  
  res <- x$results
  
  # round display-only numeric columns
  for (nm in intersect(c("pb_hat", "pc_hat", "PowerS"), names(res))) {
    if (is.numeric(res[[nm]])) {
      res[[nm]] <- round(res[[nm]], digits)
    }
  }
  
  keep_cols <- c(
    "control", "treatment", "method_label",
    "a", "b", "c", "d", "n",
    "pb_hat", "pc_hat", "PowerS",
    "N80S", "N90S"
  )
  keep_cols <- intersect(keep_cols, names(res))
  res <- res[, keep_cols, drop = FALSE]
  
  print(res, row.names = FALSE)
  invisible(x)
}

#' @rdname summary.mcnemarPower
#' @export
print.summary.mcnemarPower <- function(x, digits = .mcnemar_power_or(x$digits, 4L), ...) {
  cat("mcnemarPower summary\n")
  cat(sprintf(" p-adjust rule: %s\n", x$settings$p_adjust))
  cat(sprintf(" alpha: %.3f\n", x$settings$alpha))
  cat(sprintf(" nsim: %d\n", x$nsim))
  cat(sprintf(" subset for family power: %s\n", paste(x$comparisons, collapse = ", ")))
  cat(sprintf(
    " individual treatment target(s): %s\n\n",
    paste(x$individual_treatment, collapse = ", ")
  ))
  
  cat("Single-test results\n")
  tmp <- x$single_test
  num_cols <- intersect(c("pb_hat", "pc_hat", "PowerS"), names(tmp))
  for (nm in num_cols) tmp[[nm]] <- round(tmp[[nm]], digits)
  print(tmp, row.names = FALSE)
  cat("\n")
  
  for (mth in names(x$adjusted_power)) {
    cat("--------------------------------------------------\n")
    cat(sprintf("Method: %s\n", mth))
    for (scen in names(x$adjusted_power[[mth]])) {
      ss <- x$adjusted_power[[mth]][[scen]]
      cat(sprintf(" Scenario: %s\n", scen))
      tbl <- ss$individual_adjusted
      tbl$PowerAdj <- round(tbl$PowerAdj, digits)
      print(tbl, row.names = FALSE)
      cat(sprintf("  Any of subset significant: %.*f\n", digits, ss$subset_any))
      cat(sprintf("  All of subset significant: %.*f\n", digits, ss$subset_all))
      cat(sprintf("  All treatments significant: %.*f\n", digits, ss$all_treatments))
      cat("  Sample-size targets\n")
      
      if (!is.null(ss$sample_size$individual) && nrow(ss$sample_size$individual) > 0) {
        cat("   individual\n")
        print(ss$sample_size$individual, row.names = FALSE)
      }
      
      if (!is.null(ss$sample_size$any_subset) && nrow(ss$sample_size$any_subset) > 0) {
        cat("   any_subset\n")
        print(ss$sample_size$any_subset, row.names = FALSE)
      }
      
      if (!is.null(ss$sample_size$all_subset) && nrow(ss$sample_size$all_subset) > 0) {
        cat("   all_subset\n")
        print(ss$sample_size$all_subset, row.names = FALSE)
      }
      
      if (!is.null(ss$sample_size$all_treatments) && nrow(ss$sample_size$all_treatments) > 0) {
        cat("   all_treatments\n")
        print(ss$sample_size$all_treatments, row.names = FALSE)
      }
      
      cat("\n")
    }
  }
  invisible(x)
}

.mcnemar_power_or <- function(x, y) if (is.null(x)) y else x
