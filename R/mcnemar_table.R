#' Extract paired 2x2 contingency tables (a, b, c, d) with optional Cohen's g
#' and confidence intervals
#'
#' @description
#' After fitting \code{\link{mcnemar_control}()}, reconstruct the paired 2x2
#' table for each planned control-versus-treatment comparison and return the
#' matched cell counts. Optionally, return Cohen's \eqn{g} and one or more
#' confidence-interval representations for each treatment-versus-control
#' comparison.
#'
#' @details
#' \code{mcnemar_table()} is a reporting and effect-size companion to
#' \code{\link{mcnemar_control}()}. It takes the family of planned
#' control-versus-treatment McNemar comparisons encoded in a
#' \code{"mcnemarControl"} object, reconstructs the paired 2x2 table for each
#' treatment versus the common control, and returns the cell counts in the
#' control-row / treatment-column orientation:
#'
#' \preformatted{
#'                Treatment: No   Treatment: Yes
#' Control: No            a                b
#' Control: Yes           c                d
#' }
#'
#' In this layout, the discordant cells are \code{b} (control = No,
#' treatment = Yes) and \code{c} (control = Yes, treatment = No). Cohen's
#' \eqn{g} is computed from the discordant-pair success proportion as
#' \deqn{g = \frac{b}{b + c} - 0.5,}
#' and is therefore positive when treatment-favoring discordance dominates
#' (\code{b > c}), negative when control-favoring discordance dominates
#' (\code{c > b}), and zero when the discordant directions balance.
#'
#' By default (\code{methodCI = "none"}), no confidence intervals are returned
#' and only \code{cohens_g} is added (unless \code{include_g = FALSE} is set
#' explicitly).
#'
#' If confidence intervals are requested for Cohen's \eqn{g}, four CI modes are
#' available:
#'
#' \itemize{
#' \item \code{methodCI = "none"} returns no confidence intervals and includes
#'   \code{cohens_g} only.
#'
#' \item \code{methodCI = "marginal"} returns marginal (pointwise)
#'   confidence intervals in columns \code{ci_low} and \code{ci_high}. This
#'   matches the behavior of \code{\link[effectsize]{cohens_g}()} from the
#'   \pkg{effectsize} package, which forms a confidence interval for
#'   \eqn{P = g + 0.5} using \code{\link[stats]{prop.test}()} and shifts the
#'   limits back to the \eqn{g}-scale by subtracting 0.5.
#'   
#' \item \code{"simultaneous"}: return marginal confidence intervals in
#' \code{ci_low} and \code{ci_high}, and simultaneous adjusted confidence
#' intervals on the raw-\eqn{g} scale in \code{ci_low_adj} and
#' \code{ci_high_adj}. The adjusted intervals are obtained using a
#' subject-level bootstrap max-\eqn{|t|} procedure across the treatment family.
#' Because the adjustment is performed directly on the raw effect-size scale,
#' the resulting simultaneous limits are not constrained to the admissible
#' range of Cohen's \eqn{g} and may extend beyond \eqn{[-0.5, 0.5]}. No post
#' hoc truncation is applied. Use this option when you want the simultaneous
#' max-\eqn{|t|} intervals exactly as computed on the raw \eqn{g}-scale. 
#' 
#' \item \code{"simultaneous_logit"}: return marginal confidence intervals in
#' \code{ci_low} and \code{ci_high}, and simultaneous adjusted confidence
#' intervals in \code{ci_low_adj} and \code{ci_high_adj} obtained by applying
#' the same bootstrap max-\eqn{|t|} procedure on the logit scale of
#' \eqn{p = g + 0.5 = b/(b+c)}, followed by back-transformation to the
#' \eqn{g}-scale. Because the interval is constructed on a bounded probability
#' scale and then transformed back, the returned simultaneous limits remain
#' within the admissible range \eqn{[-0.5, 0.5]}. Use this option when
#' support-respecting simultaneous intervals are preferred, especially near the
#' boundary. This mirrors the general bounded-scale logic used by
#' \code{\link[effectsize]{cohens_g}()} for marginal intervals. 
#' 
#' }
#'
#' In contrast to the marginal intervals, both \code{"simultaneous"} and
#' \code{"simultaneous_logit"} are familywise procedures: they are intended for
#' the setting where one control is compared with several treatments and all
#' Cohen's \eqn{g} intervals are to be interpreted jointly rather than one at a
#' time. As a result, these adjusted intervals are typically wider than the
#' marginal intervals, reflecting multiplicity adjustment across the treatment
#' family.
#'
#' \strong{Degenerate simultaneous confidence intervals.}
#' In boundary/support cases (for example, one discordant direction is never
#' observed), nonparametric subject-level bootstraps can yield zero bootstrap
#' standard errors and point-mass bootstrap distributions, leading to degenerate
#' simultaneous intervals. No smoothing is applied. When this occurs,
#' \code{sim_degenerate} and \code{recommend_marginal} are returned, and
#' downstream printing methods may recommend \code{methodCI = "marginal"} for
#' affected treatments.
#'
#' \strong{Design-based interpretation.}
#' Like \code{\link{mcnemar_control}()}, this function summarizes pairwise
#' matched marginal differences under the observed repeated-condition design. It
#' does not explicitly model period effects, sequence, treatment order, or
#' carryover. When all subjects receive conditions in the same fixed order, the
#' resulting treatment-versus-control summaries should therefore be interpreted
#' as matched marginal differences under the observed presentation order rather
#' than treatment effects fully separated from time or order.
#'
#' @param fit An object of class \code{"mcnemarControl"} returned by
#'   \code{\link{mcnemar_control}()}.
#' @param data The original long-format \code{data.frame} used to fit
#'   \code{fit}. It must contain the columns specified in
#'   \code{fit$settings$id}, \code{fit$settings$condition}, and
#'   \code{fit$settings$outcome}, with unique \code{id}-\code{condition} rows.
#' @param include_treatment Logical; if \code{TRUE}, include a \code{treatment}
#'   column. If \code{FALSE}, treatment names are stored as row names. Default
#'   is \code{FALSE}.
#' @param include_diff Logical; if \code{TRUE}, include \code{diff = b - c}.
#'   Default is \code{FALSE}.
#' @param include_discordant Logical; if \code{TRUE}, include
#'   \code{discordant = b + c}. This equals the McNemar discordant-pairs count
#'   used in the test. Default is \code{FALSE}.
#' @param methodCI Character string specifying which confidence-interval method
#'   to use for Cohen's \eqn{g}. Supported values are:
#'   \itemize{
#'   \item \code{"none"}: no confidence intervals; return \code{cohens_g} only.
#'   \item \code{"marginal"}: return marginal confidence intervals in columns
#'     \code{ci_low} and \code{ci_high}.
#'   \item \code{"simultaneous"}: return marginal confidence intervals in
#'     \code{ci_low} and \code{ci_high}, and simultaneous adjusted confidence
#'     intervals on the raw-\eqn{g} scale in \code{ci_low_adj} and
#'     \code{ci_high_adj}.
#'   \item \code{"simultaneous_logit"}: return marginal confidence intervals in
#'     \code{ci_low} and \code{ci_high}, and simultaneous adjusted confidence
#'     intervals based on the logit scale in \code{ci_low_adj} and
#'     \code{ci_high_adj}.
#'   }
#' @param include_g Logical or \code{NULL}. If \code{NULL} (default), Cohen's
#'   \eqn{g} is included for all \code{methodCI} choices. If \code{FALSE} and
#'   \code{methodCI != "none"}, an error is thrown because confidence intervals
#'   imply inclusion of \code{cohens_g}.
#' @param ci_level Confidence level for confidence intervals. If \code{NULL}
#'   (default), uses \code{fit$settings$ci}. Only required / validated when CIs
#'   are requested (\code{methodCI != "none"}).
#' @param R Number of bootstrap resamples for simultaneous CI methods. Default
#'   is \code{5000}.
#' @param seed Optional integer seed for bootstrap reproducibility
#'   (simultaneous methods only).
#' @param na_policy Handling of bootstrap replicates with missing / undefined
#'   components: \code{"omit_replicates"} (default) drops replicates used for
#'   max-\eqn{t} when any required component is non-finite; \code{"pairwise"}
#'   uses pairwise standard deviations but still restricts max-\eqn{t}
#'   calibration to replicates with all required components finite.
#' @param two_sided Logical; two-sided simultaneous CIs (recommended). Default
#'   is \code{TRUE}.
#' @param verbose Logical; if \code{TRUE}, print messages from bootstrap
#'   routines and warnings about degeneracy. Default is \code{FALSE}.
#' @param digits Integer; number of decimal places used when returning
#'   \code{cohens_g} and CI columns. Default is \code{2}.
#'
#' @return
#' A \code{data.frame} containing \code{a}, \code{b}, \code{c}, and \code{d}
#' plus optional columns depending on the arguments. Depending on settings, the
#' result may include:
#'
#' \itemize{
#' \item \code{treatment} if \code{include_treatment = TRUE},
#' \item \code{diff = b - c} if \code{include_diff = TRUE},
#' \item \code{discordant = b + c} if \code{include_discordant = TRUE},
#' \item \code{cohens_g},
#' \item \code{ci_low} and \code{ci_high} for marginal confidence intervals,
#' \item \code{ci_low_adj} and \code{ci_high_adj} for simultaneous adjusted
#'   confidence intervals when \code{methodCI = "simultaneous"} or
#'   \code{methodCI = "simultaneous_logit"},
#' \item degeneracy and recommendation flags such as \code{sim_degenerate} and
#'   \code{recommend_marginal} for simultaneous CI methods.
#' }
#'
#' If \code{include_treatment = FALSE}, treatment names are stored as row names.
#' For non-bare outputs, the returned object may also carry metadata in a
#' \code{"mcnemar_table_meta"} attribute and may have class
#' \code{"mcnemar_table"} to support custom printing.
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
#' # Cell counts + Cohen's g only
#' tab_g <- mcnemar_table(
#'   fit,
#'   dat,
#'   methodCI = "none"
#' )
#'
#' # Marginal (pointwise) confidence intervals
#' tab_marg <- mcnemar_table(
#'   fit,
#'   dat,
#'   methodCI = "marginal"
#' )
#'
#' # Simultaneous adjusted confidence intervals on the raw g-scale
#' # together with marginal intervals
#' tab_sim <- mcnemar_table(
#'   fit,
#'   dat,
#'   methodCI = "simultaneous",
#'   R = 1000,
#'   seed = 1
#' )
#'
#' # Simultaneous adjusted confidence intervals on the logit scale
#' # together with marginal intervals
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
#' \item \code{\link{mcnemar_control}()} for the upstream control-versus-treatment
#'   McNemar testing workflow that defines the family of paired comparisons used
#'   here.
#' \item \code{\link{mcnemar_table_summary}()} and
#'   \code{summary.mcnemar_table()} for standardized summary tibbles and
#'   optional paired-probability summaries derived from
#'   \code{mcnemar_table()} output.
#' \item \code{\link[effectsize]{cohens_g}()} in the \pkg{effectsize} package for
#'   single-comparison Cohen's \eqn{g} intervals used in the marginal branch.
#' }
#'
#' @references
#' Cohen J (1988). \emph{Statistical Power Analysis for the Behavioral Sciences}
#'   (2nd ed.). Routledge.
#'
#' Fagerland MW, Lydersen S, Laake P (2017). \emph{Statistical Analysis of
#'   Contingency Tables}. Chapman & Hall/CRC.
#'
#' Kuchibhotla A, Kolassa J, Kuffner T (2022). Post-selection inference.
#'   \emph{Annual Review of Statistics and Its Application}, 9(1), 505--527.
#'
#' @export
mcnemar_table <- function(fit, data,
                          include_treatment = FALSE,
                          include_diff = FALSE,
                          include_discordant = FALSE,
                          methodCI = c("none", "marginal", "simultaneous", "simultaneous_logit"),
                          include_g = NULL,
                          ci_level = NULL,
                          R = 5000,
                          seed = NULL,
                          na_policy = c("omit_replicates", "pairwise"),
                          two_sided = TRUE,
                          verbose = FALSE,
                          digits = 2) {
  
  methodCI <- match.arg(methodCI)
  na_policy <- match.arg(na_policy)
  
  # digits validation
  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) || digits < 0) {
    stop("`digits` must be a single non-negative number.")
  }
  digits <- as.integer(digits)
  
  # include_g defaults to TRUE for all methodCI choices
  if (is.null(include_g)) {
    include_g <- TRUE
  } else {
    include_g <- isTRUE(include_g)
  }
  
  # if any CI method is requested, g must be included
  if (methodCI != "none" && !include_g) {
    stop("If methodCI != 'none', then include_g must be TRUE (CIs imply Cohen's g). ",
         "Set include_g = TRUE or use methodCI = 'none'.")
  }
  
  # --- Basic validation ---
  if (!inherits(fit, "mcnemarControl")) {
    stop("`fit` must be an object returned by mcnemar_control() (class 'mcnemarControl').")
  }
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame in long format (same data used to fit).")
  }
  
  s <- fit$settings
  req <- c(s$id, s$condition, s$outcome)
  missing_cols <- setdiff(req, names(data))
  if (length(missing_cols) > 0) {
    stop("`data` is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # CI level only relevant when methodCI != "none"
  if (is.null(ci_level)) ci_level <- s$ci
  if (methodCI != "none") {
    if (!is.numeric(ci_level) || length(ci_level) != 1L || ci_level <= 0 || ci_level >= 1) {
      stop("`ci_level` must be a single number in (0,1).")
    }
  }
  
  control <- as.character(s$control)
  treatments <- as.character(s$treatments)
  outcome_levels <- s$outcome_levels
  drop_na <- isTRUE(s$drop_na)
  
  # Helper: build 2x2 table with fixed level order, fill absent cells with zeros
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
  
  # Marginal CI via effectsize::cohens_g(tab)
  get_g_marginal_ci <- function(tab, ci_level) {
    if (!requireNamespace("effectsize", quietly = TRUE)) {
      stop("Package 'effectsize' is required for methodCI including 'marginal'. Please install it.")
    }
    out <- suppressWarnings(effectsize::cohens_g(tab, ci = ci_level))
    list(
      g = as.numeric(out$Cohens_g[1]),
      low = as.numeric(out$CI_low[1]),
      high = as.numeric(out$CI_high[1])
    )
  }
  
  # outcome -> 0/1 using fit$outcome_levels order (1 = second level)
  .coerce_outcome_01 <- function(x, levels2) {
    x <- factor(x, levels = levels2)
    if (anyNA(x)) return(NULL)
    as.integer(x == levels2[2])
  }
  
  # ------------------ Simultaneous max-|t| on g-scale ------------------
  sim_ci_g_scale <- function(dat_long) {
    if (!requireNamespace("boot", quietly = TRUE)) {
      stop("Package 'boot' is required for simultaneous CI methods. Please install it.")
    }
    
    d <- dat_long[, c(s$id, s$condition, s$outcome), drop = FALSE]
    d <- d[stats::complete.cases(d), , drop = FALSE]
    d[[s$condition]] <- as.factor(d[[s$condition]])
    d[[s$outcome]] <- factor(d[[s$outcome]], levels = outcome_levels)
    
    wide <- tidyr::pivot_wider(
      d,
      names_from = dplyr::all_of(s$condition),
      values_from = dplyr::all_of(s$outcome)
    )
    
    needed_cols <- c(s$id, control, treatments)
    missing_cols2 <- setdiff(needed_cols, names(wide))
    if (length(missing_cols2) > 0) {
      stop("Missing columns after pivot: ", paste(missing_cols2, collapse = ", "))
    }
    
    cc <- stats::complete.cases(wide[, needed_cols, drop = FALSE])
    if (!all(cc)) {
      if (verbose) {
        message("Dropping ", sum(!cc), " subjects with incomplete profiles (simultaneous g-scale).")
      }
      wide <- wide[cc, , drop = FALSE]
    }
    
    n_ids <- nrow(wide)
    if (n_ids < 2) stop("Need at least 2 subjects after filtering for simultaneous CI.")
    
    for (nm in c(control, treatments)) {
      y01 <- .coerce_outcome_01(wide[[nm]], outcome_levels)
      if (is.null(y01)) {
        stop("Missing outcomes after coercion; check outcome_levels and missingness.")
      }
      wide[[nm]] <- y01
    }
    
    g_hat <- vapply(treatments, function(tr) {
      ctrl <- wide[[control]]
      trt <- wide[[tr]]
      b <- sum(ctrl == 0 & trt == 1)
      c <- sum(ctrl == 1 & trt == 0)
      if ((b + c) == 0) return(NA_real_)
      b / (b + c) - 0.5
    }, numeric(1))
    
    boot_stat <- function(data0, indices) {
      dd <- data0[indices, , drop = FALSE]
      ctrl <- dd[[control]]
      vapply(treatments, function(tr) {
        trt <- dd[[tr]]
        b <- sum(ctrl == 0 & trt == 1)
        c <- sum(ctrl == 1 & trt == 0)
        if ((b + c) == 0) return(NA_real_)
        b / (b + c) - 0.5
      }, numeric(1))
    }
    
    if (!is.null(seed)) set.seed(seed)
    boot_out <- boot::boot(data = wide, statistic = boot_stat, R = R)
    Tmat <- boot_out$t
    
    if (na_policy == "omit_replicates") {
      ok <- apply(Tmat, 1, function(row) all(is.finite(row)))
      Tmat2 <- Tmat[ok, , drop = FALSE]
      if (verbose) {
        message("Valid bootstrap replicates (g-scale): ", nrow(Tmat2), " / ", nrow(Tmat))
      }
    } else {
      Tmat2 <- Tmat
    }
    
    se <- apply(Tmat2, 2, stats::sd, na.rm = TRUE)
    cols_use <- which(is.finite(se) & se > 0 & is.finite(g_hat))
    rows_for_maxt <- rep(FALSE, nrow(Tmat2))
    
    if (length(cols_use) > 0) {
      rows_for_maxt <- apply(
        Tmat2[, cols_use, drop = FALSE],
        1,
        function(row) all(is.finite(row))
      )
    }
    
    if (verbose) {
      message("Bootstrap SE summary (g-scale): min=", min(se, na.rm = TRUE),
              " median=", stats::median(se, na.rm = TRUE),
              " max=", max(se, na.rm = TRUE),
              " zeros=", sum(se == 0, na.rm = TRUE),
              " NA=", sum(!is.finite(se)))
      message("Usable rows for max-t (g-scale): ", sum(rows_for_maxt), " / ", nrow(Tmat2))
    }
    
    crit <- 0
    if (length(cols_use) > 0 && sum(rows_for_maxt) > 0) {
      Z <- sweep(Tmat2[rows_for_maxt, cols_use, drop = FALSE], 2, g_hat[cols_use], "-")
      for (j in seq_along(cols_use)) {
        sj <- se[cols_use[j]]
        if (!is.finite(sj) || sj <= 0) {
          Z[, j] <- 0
        } else {
          Z[, j] <- Z[, j] / sj
        }
      }
      maxabs <- apply(Z, 1, function(z) max(abs(z), na.rm = TRUE))
      maxabs <- maxabs[is.finite(maxabs)]
      if (length(maxabs) > 0L) {
        alpha <- 1 - ci_level
        crit <- stats::quantile(maxabs, probs = 1 - alpha, names = FALSE, type = 8, na.rm = TRUE)
      }
    }
    
    ci_low <- g_hat - crit * se
    ci_high <- g_hat + crit * se
    
    degenerate <- (!is.finite(se)) | (se == 0) | (!is.finite(g_hat))
    ci_low[!is.finite(g_hat)] <- NA_real_
    ci_high[!is.finite(g_hat)] <- NA_real_
    
    list(
      ci_low = stats::setNames(ci_low, treatments),
      ci_high = stats::setNames(ci_high, treatments),
      se = stats::setNames(se, treatments),
      degenerate = stats::setNames(degenerate, treatments),
      crit = crit,
      n_subjects = n_ids
    )
  }
  
  # ------------------ Simultaneous max-|t| on logit(p) scale ------------------
  sim_ci_logit_scale <- function(dat_long) {
    if (!requireNamespace("boot", quietly = TRUE)) {
      stop("Package 'boot' is required for simultaneous CI methods. Please install it.")
    }
    
    d <- dat_long[, c(s$id, s$condition, s$outcome), drop = FALSE]
    d <- d[stats::complete.cases(d), , drop = FALSE]
    d[[s$condition]] <- as.factor(d[[s$condition]])
    d[[s$outcome]] <- factor(d[[s$outcome]], levels = outcome_levels)
    
    wide <- tidyr::pivot_wider(
      d,
      names_from = dplyr::all_of(s$condition),
      values_from = dplyr::all_of(s$outcome)
    )
    
    needed_cols <- c(s$id, control, treatments)
    missing_cols2 <- setdiff(needed_cols, names(wide))
    if (length(missing_cols2) > 0) {
      stop("Missing columns after pivot: ", paste(missing_cols2, collapse = ", "))
    }
    
    cc <- stats::complete.cases(wide[, needed_cols, drop = FALSE])
    if (!all(cc)) {
      if (verbose) {
        message("Dropping ", sum(!cc), " subjects with incomplete profiles (simultaneous logit).")
      }
      wide <- wide[cc, , drop = FALSE]
    }
    
    n_ids <- nrow(wide)
    if (n_ids < 2) stop("Need at least 2 subjects after filtering for simultaneous CI.")
    
    for (nm in c(control, treatments)) {
      y01 <- .coerce_outcome_01(wide[[nm]], outcome_levels)
      if (is.null(y01)) {
        stop("Missing outcomes after coercion; check outcome_levels and missingness.")
      }
      wide[[nm]] <- y01
    }
    
    p_hat <- vapply(treatments, function(tr) {
      ctrl <- wide[[control]]
      trt <- wide[[tr]]
      b <- sum(ctrl == 0 & trt == 1)
      c <- sum(ctrl == 1 & trt == 0)
      if ((b + c) == 0) return(NA_real_)
      b / (b + c)
    }, numeric(1))
    
    eta_hat <- rep(NA_real_, length(treatments))
    names(eta_hat) <- treatments
    finite_mask <- is.finite(p_hat) & (p_hat > 0) & (p_hat < 1)
    eta_hat[finite_mask] <- stats::qlogis(p_hat[finite_mask])
    
    # observed boundary or no-discordant -> degenerate and excluded from max-t calibration
    degenerate_obs <- !finite_mask
    
    boot_stat <- function(data0, indices) {
      dd <- data0[indices, , drop = FALSE]
      ctrl <- dd[[control]]
      vapply(treatments, function(tr) {
        trt <- dd[[tr]]
        b <- sum(ctrl == 0 & trt == 1)
        c <- sum(ctrl == 1 & trt == 0)
        if ((b + c) == 0) return(NA_real_)
        p <- b / (b + c)
        if (p <= 0 || p >= 1) return(NA_real_)
        stats::qlogis(p)
      }, numeric(1))
    }
    
    if (!is.null(seed)) set.seed(seed)
    boot_out <- boot::boot(data = wide, statistic = boot_stat, R = R)
    Tmat <- boot_out$t
    
    se_eta <- apply(Tmat, 2, stats::sd, na.rm = TRUE)
    cols_use <- which(!degenerate_obs & is.finite(se_eta) & se_eta > 0)
    rows_for_maxt <- rep(FALSE, nrow(Tmat))
    
    if (length(cols_use) > 0) {
      rows_for_maxt <- apply(
        Tmat[, cols_use, drop = FALSE],
        1,
        function(row) all(is.finite(row))
      )
    }
    
    if (verbose) {
      message("Bootstrap SE (eta/logit) summary: min=", min(se_eta, na.rm = TRUE),
              " median=", stats::median(se_eta, na.rm = TRUE),
              " max=", max(se_eta, na.rm = TRUE),
              " zeros=", sum(se_eta == 0, na.rm = TRUE),
              " NA=", sum(!is.finite(se_eta)))
      message("Usable rows for max-t (logit): ", sum(rows_for_maxt), " / ", nrow(Tmat))
      if (any(degenerate_obs)) {
        message("Degenerate (boundary/no-discordant) treatments (logit): ",
                paste(treatments[degenerate_obs], collapse = ", "))
      }
    }
    
    crit <- 0
    if (length(cols_use) > 0 && sum(rows_for_maxt) > 0) {
      Z <- sweep(Tmat[rows_for_maxt, cols_use, drop = FALSE], 2, eta_hat[cols_use], "-")
      for (j in seq_along(cols_use)) {
        sj <- se_eta[cols_use[j]]
        if (!is.finite(sj) || sj <= 0) {
          Z[, j] <- 0
        } else {
          Z[, j] <- Z[, j] / sj
        }
      }
      maxabs <- apply(Z, 1, function(z) max(abs(z), na.rm = TRUE))
      maxabs <- maxabs[is.finite(maxabs)]
      if (length(maxabs) > 0L) {
        alpha <- 1 - ci_level
        crit <- stats::quantile(maxabs, probs = 1 - alpha, names = FALSE, type = 8, na.rm = TRUE)
      }
    }
    
    ci_eta_low <- rep(NA_real_, length(treatments))
    ci_eta_high <- rep(NA_real_, length(treatments))
    names(ci_eta_low) <- names(ci_eta_high) <- treatments
    ci_eta_low[!degenerate_obs] <- eta_hat[!degenerate_obs] - crit * se_eta[!degenerate_obs]
    ci_eta_high[!degenerate_obs] <- eta_hat[!degenerate_obs] + crit * se_eta[!degenerate_obs]
    
    # back-transform to g scale
    ci_g_low <- rep(NA_real_, length(treatments))
    ci_g_high <- rep(NA_real_, length(treatments))
    names(ci_g_low) <- names(ci_g_high) <- treatments
    
    ci_g_low[!degenerate_obs] <- stats::plogis(ci_eta_low[!degenerate_obs]) - 0.5
    ci_g_high[!degenerate_obs] <- stats::plogis(ci_eta_high[!degenerate_obs]) - 0.5
    
    # collapse to boundary if observed p is 0 or 1 (no smoothing)
    is_p0 <- is.finite(p_hat) & p_hat == 0
    is_p1 <- is.finite(p_hat) & p_hat == 1
    ci_g_low[is_p0] <- -0.5
    ci_g_high[is_p0] <- -0.5
    ci_g_low[is_p1] <- 0.5
    ci_g_high[is_p1] <- 0.5
    
    degenerate <- degenerate_obs | (!is.finite(se_eta)) | (se_eta == 0)
    
    list(
      ci_low = stats::setNames(ci_g_low, treatments),
      ci_high = stats::setNames(ci_g_high, treatments),
      se = stats::setNames(se_eta, treatments),
      degenerate = stats::setNames(degenerate, treatments),
      crit = crit,
      n_subjects = n_ids
    )
  }
  
  # --- Filter to relevant conditions and enforce uniqueness ---
  cond_chr <- as.character(data[[s$condition]])
  df <- data[cond_chr %in% c(control, treatments), , drop = FALSE]
  
  key <- paste(df[[s$id]], df[[s$condition]], sep = "__")
  if (any(duplicated(key))) {
    stop(
      "Multiple rows per id-condition detected in `data`. ",
      "This helper requires unique id-condition rows (same as mcnemar_control())."
    )
  }
  
  # Wide for cell reconstruction
  wide <- tidyr::pivot_wider(
    df,
    id_cols = dplyr::all_of(s$id),
    names_from = dplyr::all_of(s$condition),
    values_from = dplyr::all_of(s$outcome)
  )
  
  needs_marginal <- methodCI %in% c("marginal", "simultaneous", "simultaneous_logit")
  needs_sim_g <- identical(methodCI, "simultaneous")
  needs_sim_logit <- identical(methodCI, "simultaneous_logit")
  
  # --- Build base rows: a,b,c,d (+ optional) + cohens_g + marginal CIs when needed ---
  out_rows <- vector("list", length(treatments))
  
  for (i in seq_along(treatments)) {
    trt <- treatments[i]
    
    x <- wide[[control]]
    y <- wide[[trt]]
    
    if (drop_na) {
      keep <- !is.na(x) & !is.na(y)
      x <- x[keep]
      y <- y[keep]
    }
    
    tab <- build_table(x, y, outcome_levels)
    
    a <- as.integer(tab[1, 1])
    b <- as.integer(tab[1, 2])
    c <- as.integer(tab[2, 1])
    d_ <- as.integer(tab[2, 2])
    
    diff <- b - c
    discordant <- b + c
    
    if (isTRUE(include_treatment)) {
      res_i <- data.frame(
        treatment = trt,
        a = a, b = b, c = c, d = d_,
        stringsAsFactors = FALSE
      )
    } else {
      res_i <- data.frame(
        a = a, b = b, c = c, d = d_,
        stringsAsFactors = FALSE
      )
      rownames(res_i) <- trt
    }
    
    if (isTRUE(include_diff)) res_i$diff <- diff
    if (isTRUE(include_discordant)) res_i$discordant <- discordant
    
    if (include_g) {
      # default g from discordants
      res_i$cohens_g <- if (discordant == 0) NA_real_ else b / discordant - 0.5
      
      # marginal CI if requested (and overwrite g with effectsize g for consistency)
      if (needs_marginal) {
        gm <- get_g_marginal_ci(tab, ci_level)
        if (is.finite(gm$g)) res_i$cohens_g <- gm$g
        res_i$ci_low <- gm$low
        res_i$ci_high <- gm$high
      }
    }
    
    out_rows[[i]] <- res_i
  }
  
  res <- do.call(rbind, out_rows)
  if (isTRUE(include_treatment)) rownames(res) <- NULL
  
  # --- Add simultaneous adjusted CIs ---
  rec_marg <- rep(FALSE, nrow(res))
  
  # raw g-scale simultaneous adjusted CI
  if (needs_sim_g) {
    sim_g <- sim_ci_g_scale(data)
    
    if (isTRUE(include_treatment)) {
      res$ci_low_adj <- sim_g$ci_low[res$treatment]
      res$ci_high_adj <- sim_g$ci_high[res$treatment]
      res$sim_degenerate <- as.logical(sim_g$degenerate[res$treatment])
    } else {
      res$ci_low_adj <- sim_g$ci_low[rownames(res)]
      res$ci_high_adj <- sim_g$ci_high[rownames(res)]
      res$sim_degenerate <- as.logical(sim_g$degenerate[rownames(res)])
    }
    rec_marg <- rec_marg | res$sim_degenerate
  }
  
  # logit-scale simultaneous adjusted CI
  if (needs_sim_logit) {
    sim_l <- sim_ci_logit_scale(data)
    
    if (isTRUE(include_treatment)) {
      res$ci_low_adj <- sim_l$ci_low[res$treatment]
      res$ci_high_adj <- sim_l$ci_high[res$treatment]
      res$sim_degenerate <- as.logical(sim_l$degenerate[res$treatment])
    } else {
      res$ci_low_adj <- sim_l$ci_low[rownames(res)]
      res$ci_high_adj <- sim_l$ci_high[rownames(res)]
      res$sim_degenerate <- as.logical(sim_l$degenerate[rownames(res)])
    }
    rec_marg <- rec_marg | res$sim_degenerate
  }
  
  # --- Recommendation flags (kept for programmatic use; printing handled by print method) ---
  if (needs_sim_g || needs_sim_logit) {
    res$recommend_marginal <- rec_marg
  }
  
  # --- Store header metadata + degenerate treatments for print.mcnemar_table ---
  is_bare <- setequal(names(res), c("a", "b", "c", "d"))
  
  if (!is_bare) {
    attr(res, "mcnemar_table_meta") <- list(
      control = control,
      treatments = treatments,
      methodCI = methodCI,
      ci_level = if (methodCI == "none") NA_real_ else ci_level,
      R = if (methodCI %in% c("simultaneous", "simultaneous_logit")) R else NA_integer_,
      seed = if (methodCI %in% c("simultaneous", "simultaneous_logit")) seed else NA_integer_,
      digits = digits
    )
    
    if ("recommend_marginal" %in% names(res) && any(res$recommend_marginal %in% TRUE, na.rm = TRUE)) {
      idx <- which(res$recommend_marginal %in% TRUE)
      deg_trt <- if ("treatment" %in% names(res)) {
        as.character(res$treatment[idx])
      } else {
        rownames(res)[idx]
      }
      attr(res, "degenerate_treatments") <- unique(deg_trt)
    } else {
      attr(res, "degenerate_treatments") <- character(0)
    }
    
    class(res) <- c("mcnemar_table", class(res))
  }
  
  # --- Round effect size and CI columns for display ---
  cols_round <- c(
    "cohens_g",
    "ci_low", "ci_high",
    "ci_low_adj", "ci_high_adj"
  )
  cols_present <- intersect(cols_round, names(res))
  for (nm in cols_present) {
    if (is.numeric(res[[nm]])) {
      res[[nm]] <- round(res[[nm]], digits = digits)
    }
  }
  
  res
}

#' Summary tibble for \code{mcnemar_table()} results
#'
#' @description
#' Standardizes the output of \code{\link{mcnemar_table}()} into a tibble-like
#' summary with a consistent column layout for matched cell counts, Cohen's
#' \eqn{g}, marginal confidence intervals, and (when present) adjusted
#' simultaneous confidence intervals. Optionally, it can also compute:
#' \itemize{
#' \item a second tibble of differences between marginal proportions and their confidence intervals;
#' \item a third tibble of ratios of marginal proportions and their confidence intervals; and
#' \item a fourth tibble of conditional odds ratios and their confidence intervals.
#' }
#'
#' @details
#' The summary constructor is designed for reporting and inspection. It takes
#' the output returned by \code{\link{mcnemar_table}()} and builds a
#' standardized summary table. If the original table used row names for
#' treatments (\code{include_treatment = FALSE}), treatment names are recovered
#' from the row names.
#'
#' In the current \code{mcnemar_table()} workflow, the primary confidence
#' interval columns are standardized by construction:
#'
#' \itemize{
#' \item \code{ci_low} and \code{ci_high} represent the marginal (pointwise)
#' confidence intervals for Cohen's \eqn{g}.
#' \item \code{ci_low_adj} and \code{ci_high_adj}, when present, represent the
#' adjusted simultaneous confidence intervals returned by
#' \code{mcnemar_table()} for \code{methodCI = "simultaneous"} or
#' \code{methodCI = "simultaneous_logit"}.
#' }
#'
#' If \code{include_flags = TRUE}, the summary tibble also includes
#' \code{sim_degenerate} (when available) and
#' \code{recommend_marginal} (when available). If multiple simultaneous
#' degeneracy flags are present in an older object (for example
#' \code{sim_degenerate_sim} and \code{sim_degenerate_logit}), they are
#' collapsed into a single logical \code{sim_degenerate} column.
#'
#' If \code{pairedP = TRUE}, a second tibble is created containing differences
#' between marginal proportions and confidence intervals computed from the
#' paired 2x2 tables using one or more methods from the
#' \pkg{contingencytables} package. The estimate column in this tibble is named
#' \code{Diff}. Supported values of \code{pairedPCI} are:
#' \itemize{
#' \item \code{"Wald"}: Wald confidence interval for the difference between paired proportions;
#' \item \code{"cc"}: continuity-corrected Wald confidence interval for the difference between paired proportions;
#' \item \code{"BonettPrice"}: Bonett-Price Wald-type confidence interval for the difference between paired proportions;
#' \item \code{"AgrestiMin"}: Agresti-Min Wald-type confidence interval for the difference between paired proportions;
#' \item \code{"Newcombe"}: Newcombe square-and-add confidence interval for the difference between paired proportions;
#' \item \code{"Tango"}: Tango asymptotic score confidence interval for the difference between paired proportions.
#' }
#'
#'
#' When exactly one \code{pairedPCI} method is selected, the displayed
#' \code{Diff} estimate is taken from the estimate returned by that
#' \pkg{contingencytables} method, after alignment to the
#' treatment-minus-control direction used by this summary. When multiple
#' \code{pairedPCI} methods are selected, the main \code{Diff} column retains
#' the direct paired-difference estimate \eqn{(b-c)/N}, and method-specific
#' confidence-interval columns are added for each requested method.
#'
#'
#' If \code{ratioP = TRUE}, a third tibble is created containing the ratio of
#' marginal proportions and confidence intervals computed from the same paired
#' 2x2 tables. Supported values of \code{ratioCI} are:
#' \itemize{
#' \item \code{"BonettPrice"}: Bonett-Price hybrid Wilson score confidence interval for the ratio of paired proportions;
#' \item \code{"Asymptotic"}: Tang asymptotic score confidence interval for the ratio of paired proportions;
#' \item \code{"Mover"}: MOVER Wilson score confidence interval for the ratio of paired proportions;
#' \item \code{"BonettPriceCC"}: continuity-corrected Bonett-Price hybrid Wilson score confidence interval for the ratio of paired proportions;
#' \item \code{"Wald"}: Wald confidence interval for the ratio of paired proportions.
#' }
#'
#'
#' In the current table orientation used by \code{mcnemar_table()}, the
#' marginal success probabilities are
#' \eqn{P(\mathrm{Control} = \mathrm{Yes}) = (c+d)/N} and
#' \eqn{P(\mathrm{Treatment} = \mathrm{Yes}) = (b+d)/N}. The third tibble
#' therefore reports \code{control = (c+d)/N} and
#' \code{treatment = (b+d)/N}. The estimate column in this third tibble is
#' named \code{Ratio}. To make the selected paired-ratio method correspond to
#' the treatment-over-control success ratio, the paired 2x2 table is
#' reoriented before calling the \pkg{contingencytables} ratio procedure so
#' that the returned estimate is interpreted as
#' \eqn{P(\mathrm{Treatment} = \mathrm{Yes}) /
#' P(\mathrm{Control} = \mathrm{Yes})}. The displayed \code{Ratio} and its
#' confidence limits are therefore the method-specific estimate and interval
#' for the Treatment Yes / Control Yes ratio under that reoriented table.
#'
#'
#' If \code{OR = TRUE}, a fourth tibble is created containing conditional odds
#' ratios and confidence intervals computed from the same paired 2x2 tables
#' using a selected method from the \pkg{contingencytables} package. The
#' estimate column in this tibble is named \code{OR}. Supported values of
#' \code{ORCI} are:
#' \itemize{
#' \item \code{"Wilson"}: transformed Wilson score confidence interval for the conditional odds ratio;
#' \item \code{"ClopperPearson_midP"}: transformed Clopper-Pearson mid-P confidence interval for the conditional odds ratio;
#' \item \code{"ClopperPearson"}: transformed Clopper-Pearson exact confidence interval for the conditional odds ratio;
#' \item \code{"Blaker"}: transformed Blaker exact confidence interval for the conditional odds ratio;
#' \item \code{"Wald"}: Wald confidence interval for the conditional odds ratio;
#' \item \code{"Wald_Laplace"}: Wald confidence interval for the conditional odds ratio with Laplace adjustment.
#' }
#'
#' In the table orientation used by \code{mcnemar_table()},
#' \preformatted{
#'                 Treatment: No   Treatment: Yes
#' Control: No            a                b
#' Control: Yes           c                d
#' }
#' the discordant cells are \code{b} (control = No, treatment = Yes) and
#' \code{c} (control = Yes, treatment = No). The conditional odds ratio
#' reported by the fourth tibble is aligned to this summary's direction, so the
#' displayed estimate and confidence interval correspond to the discordant-cell
#' ordering \eqn{b/c}. If the orientation used internally by a selected
#' \pkg{contingencytables} function returns the reciprocal direction, the
#' estimate and confidence interval are inverted before being shown in the
#' tibble so that all reported conditional odds ratios follow the same
#' treatment-favoring versus control-favoring convention as the second and
#' third tibbles.
#'
#' If \code{p_adjust = "bonferroni"}, each selected paired CI method is
#' recomputed at \eqn{\alpha / K}, where \eqn{K} is the number of
#' treatment-versus-control comparisons, and adjusted CI columns are added.
#' This multiplicity adjustment may be used when \code{pairedP = TRUE},
#' \code{ratioP = TRUE}, and/or \code{OR = TRUE}.
#'
#' If \code{p_adjust = "bootstrap_maxt"}, simultaneous adjusted CIs are
#' computed using a subject-level bootstrap max-\eqn{t} procedure; this
#' requires the original \code{fit} and \code{data} arguments so the repeated
#' subject-level outcomes can be reconstructed. For paired differences, the
#' max-\eqn{t} calibration is applied on the paired-difference scale. For
#' ratios of marginal proportions, the max-\eqn{t} calibration is applied on
#' the log-ratio scale and then back-transformed to the ratio scale. For
#' conditional odds ratios, the max-\eqn{t} calibration is applied on the
#' log-odds-ratio scale and then back-transformed to the odds-ratio scale.
#'
#' The summary object is intended for reporting and downstream printing; it does
#' not re-fit McNemar tests or recompute Cohen's \eqn{g}. Instead, it organizes
#' the information already returned by \code{\link{mcnemar_table}()} into a more
#' stable rectangular format and, if requested, augments it with summaries of
#' differences and/or ratios between marginal proportions and/or conditional
#' odds ratios.
#'
#' @section S3 methods:
#' \itemize{
#' \item \code{summary.mcnemar_table()} dispatches to
#' \code{mcnemar_table_summary()} and returns the standardized summary object.
#' \item \code{print.summary.mcnemar_table()} formats and prints that summary
#' object, including the overview metadata, the main summary tibble, the
#' second tibble of differences between marginal proportions when present, the
#' third tibble of ratios of marginal proportions when present, and the fourth
#' tibble of conditional odds ratios when present.
#' }
#'
#' @param x An object returned by \code{\link{mcnemar_table}()}.
#' @param object An object returned by \code{\link{mcnemar_table}()}.
#' Used by the S3 method \code{summary.mcnemar_table()}.
#' @param include_flags Logical; if \code{TRUE}, keep the standardized
#' \code{sim_degenerate} column (when available) and
#' \code{recommend_marginal} column (when available). Default is
#' \code{FALSE}.
#' @param pairedP Logical; if \code{TRUE}, include a second tibble containing
#' differences between marginal proportions and confidence intervals. Default
#' is \code{FALSE}.
#' @param pairedPCI Character vector specifying one or more CI methods for the
#' differences between marginal proportions. Allowed values are
#' \code{"BonettPrice"} (default), \code{"AgrestiMin"},
#' \code{"Newcombe"}, \code{"Tango"}, \code{"Wald"}, and \code{"cc"}.
#' Multiple methods may be supplied, e.g.
#' \code{pairedPCI = c("BonettPrice", "Newcombe")}.
#' @param ratioP Logical; if \code{TRUE}, include a third tibble containing
#' ratios of marginal proportions and confidence intervals. Default is
#' \code{FALSE}.
#' @param ratioCI Character string specifying the confidence-interval method
#' for the ratio of marginal proportions. Allowed values are
#' \code{"BonettPrice"} (default), \code{"Asymptotic"},
#' \code{"Mover"}, \code{"BonettPriceCC"}, and \code{"Wald"}.
#' @param OR Logical; if \code{TRUE}, include a fourth tibble containing
#' conditional odds ratios and confidence intervals. Default is
#' \code{FALSE}.
#' @param ORCI Character string specifying the confidence-interval method
#' for the conditional odds ratio. Allowed values are
#' \code{"Wilson"} (default),
#' \code{"ClopperPearson_midP"},
#' \code{"ClopperPearson"},
#' \code{"Blaker"},
#' \code{"Wald"},
#' and \code{"Wald_Laplace"}.
#' @param p_adjust Multiplicity adjustment for the optional paired-confidence
#' interval tibbles: \code{"none"} (default), \code{"bonferroni"}, or
#' \code{"bootstrap_maxt"}. This adjustment may be used when
#' \code{pairedP = TRUE}, \code{ratioP = TRUE}, and/or \code{OR = TRUE}.
#' @param fit Optional object returned by
#' \code{\link[=mcnemar_control]{mcnemar_control}()}. Required only when
#' \code{p_adjust = "bootstrap_maxt"} and at least one of
#' \code{pairedP}, \code{ratioP}, or \code{OR} is \code{TRUE}.
#' @param data Optional original long-format data used to fit \code{fit}.
#' Required only when \code{p_adjust = "bootstrap_maxt"} and at least one of
#' \code{pairedP}, \code{ratioP}, or \code{OR} is \code{TRUE}.
#' @param R Number of bootstrap resamples used when
#' \code{p_adjust = "bootstrap_maxt"}. Default is \code{5000}.
#' @param seed Optional integer seed used when
#' \code{p_adjust = "bootstrap_maxt"}.
#' @param na_policy Handling of non-finite bootstrap replicates for
#' \code{p_adjust = "bootstrap_maxt"}: \code{"omit_replicates"} (default) or
#' \code{"pairwise"}.
#' @param control Optional character override for the control label. By default,
#' the function uses the control stored in the \code{"mcnemar_table_meta"}
#' attribute when available.
#' @param digits Optional non-negative integer. If supplied, round
#' \code{cohens_g}, \code{ci_low}, \code{ci_high}, and (when present)
#' \code{ci_low_adj} and \code{ci_high_adj} in the main summary tibble, and
#' numeric columns in the second, third, and fourth tibbles, to this many
#' digits.
#' @param ... Additional arguments ignored by the summary constructor or passed
#' through to printing methods.
#'
#' @return
#' \code{mcnemar_table_summary()} and \code{summary.mcnemar_table()} return an
#' object of class \code{"summary.mcnemar_table"} with elements:
#' \itemize{
#' \item \code{overview}: a list containing summary metadata, including control
#' label, CI method, confidence level, whether adjusted CI columns were
#' present, and settings for the optional second, third, and fourth tibbles.
#' \item \code{table}: a tibble-like object with standardized columns for
#' \code{control}, \code{treatment}, \code{a}, \code{b}, \code{c},
#' \code{d}, \code{N}, \code{cohens_g}, \code{ci_low}, and
#' \code{ci_high}; and, when available, \code{ci_low_adj} and
#' \code{ci_high_adj}.
#' \item \code{paired_probabilities}: either \code{NULL} or a tibble of
#' differences between marginal proportions and confidence intervals; when a
#' single \code{pairedPCI} method is selected, the displayed \code{Diff}
#' estimate is taken from that method's returned estimate after alignment to
#' the treatment-minus-control direction, whereas with multiple
#' \code{pairedPCI} methods the main \code{Diff} column retains the direct
#' paired-difference estimate and method-specific confidence-interval columns
#' are added.
#' \item \code{ratio_probabilities}: either \code{NULL} or a tibble containing
#' the control and treatment marginal success probabilities, together with the
#' treatment-over-control ratio and its confidence interval as returned by the
#' selected paired-ratio method after reorienting the paired 2x2 table so that
#' the reported \code{Ratio} is interpreted as
#' \eqn{P(\mathrm{Treatment} = \mathrm{Yes}) /
#' P(\mathrm{Control} = \mathrm{Yes})}.
#' \item \code{odds_ratios}: either \code{NULL} or a tibble containing the
#' control label, treatment label, the conditional odds ratio \code{OR}, and
#' its confidence interval. The displayed conditional odds ratio and its
#' confidence limits are aligned to the discordant-cell direction
#' \eqn{b/c}, where \code{b} is the upper-right cell
#' (control = No, treatment = Yes) and \code{c} is the lower-left cell
#' (control = Yes, treatment = No). When \code{p_adjust} is
#' \code{"bonferroni"} or \code{"bootstrap_maxt"}, adjusted CI columns
#' \code{ci_low_adj} and \code{ci_high_adj} are also included.
#' \item \code{original}: the original \code{\link{mcnemar_table}()} result.
#' }
#'
#' \code{print.summary.mcnemar_table()} prints the summary object in a compact,
#' human-readable format and returns it invisibly.
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
#' # Marginal CI summary
#' tab_marg <- mcnemar_table(
#'   fit,
#'   dat,
#'   methodCI = "marginal"
#' )
#' s_marg <- mcnemar_table_summary(tab_marg)
#' s_marg
#'
#' # Simultaneous adjusted CI summary (raw g-scale)
#' tab_sim <- mcnemar_table(
#'   fit,
#'   dat,
#'   methodCI = "simultaneous",
#'   R = 1000,
#'   seed = 1
#' )
#' s_sim <- mcnemar_table_summary(tab_sim)
#' s_sim
#'
#' # The summary table automatically includes ci_low_adj / ci_high_adj
#' # when they are present in the original mcnemar_table() output
#' s_sim$table
#'
#' # Include degeneracy / recommendation flags when present
#' s_flags <- mcnemar_table_summary(
#'   tab_sim,
#'   include_flags = TRUE
#' )
#'
#' # Add differences between marginal proportions and their confidence intervals
#' s_paired <- mcnemar_table_summary(
#'   tab_sim,
#'   pairedP = TRUE,
#'   pairedPCI = c("BonettPrice", "Newcombe")
#' )
#'
#' # Add ratios of marginal proportions and their confidence intervals
#' s_ratio <- mcnemar_table_summary(
#'   tab_sim,
#'   ratioP = TRUE,
#'   ratioCI = "BonettPrice"
#' )
#'
#' # Add differences between marginal proportions with simultaneous adjusted CIs
#' s_boot <- mcnemar_table_summary(
#'   tab_sim,
#'   pairedP = TRUE,
#'   p_adjust = "bootstrap_maxt",
#'   fit = fit,
#'   data = dat,
#'   R = 1000,
#'   seed = 1,
#'   na_policy = "omit_replicates"
#' )
#'
#' # Add ratios of marginal proportions with simultaneous adjusted CIs
#' s_ratio_boot <- mcnemar_table_summary(
#'   tab_sim,
#'   ratioP = TRUE,
#'   ratioCI = "BonettPrice",
#'   p_adjust = "bootstrap_maxt",
#'   fit = fit,
#'   data = dat,
#'   R = 1000,
#'   seed = 1,
#'   na_policy = "omit_replicates"
#' )
#'
#' # Add conditional odds ratios and their confidence intervals
#' s_or <- mcnemar_table_summary(
#'   tab_sim,
#'   OR = TRUE,
#'   ORCI = "Wilson"
#' )
#'
#' # Add conditional odds ratios with Bonferroni-adjusted confidence intervals
#' s_or_bonf <- mcnemar_table_summary(
#'   tab_sim,
#'   OR = TRUE,
#'   ORCI = "Blaker",
#'   p_adjust = "bonferroni"
#' )
#'
#' # Add conditional odds ratios with simultaneous adjusted CIs
#' s_or_boot <- mcnemar_table_summary(
#'   tab_sim,
#'   OR = TRUE,
#'   ORCI = "Wald_Laplace",
#'   p_adjust = "bootstrap_maxt",
#'   fit = fit,
#'   data = dat,
#'   R = 1000,
#'   seed = 1,
#'   na_policy = "omit_replicates"
#' )
#'
#' # S3 summary method dispatches to the same constructor
#' summary(tab_sim)
#'
#' @seealso
#' \itemize{
#' \item \code{\link{mcnemar_table}()} for reconstruction of paired 2x2 tables,
#' Cohen's \eqn{g}, marginal confidence intervals, and optional adjusted
#' simultaneous confidence intervals.
#' \item \code{\link[effectsize]{cohens_g}()} in the \pkg{effectsize} package for
#' single-comparison Cohen's \eqn{g} intervals used in the marginal branch.
#' }
#'
#' @references
#' Fagerland MW, Lydersen S, Laake P (2017). \emph{Statistical Analysis of
#' Contingency Tables}. Chapman & Hall/CRC.
#'
#' Kuchibhotla A, Kolassa J, Kuffner T (2022). Post-selection inference.
#' \emph{Annual Review of Statistics and Its Application}, 9(1), 505--527.
#' 
#' @export
mcnemar_table_summary <- function(x,
                                  include_flags = FALSE,
                                  pairedP = FALSE,
                                  pairedPCI = c("BonettPrice", "AgrestiMin", "Newcombe", "Tango", "Wald", "cc"),
                                  ratioP = FALSE,
                                  ratioCI = c("BonettPrice", "Asymptotic", "Mover", "BonettPriceCC", "Wald"),
                                  OR = FALSE,
                                  ORCI = c("Wilson", "ClopperPearson_midP", "ClopperPearson", "Blaker", "Wald", "Wald_Laplace"),
                                  p_adjust = c("none", "bonferroni", "bootstrap_maxt"),
                                  fit = NULL,
                                  data = NULL,
                                  R = 5000,
                                  seed = NULL,
                                  na_policy = c("omit_replicates", "pairwise"),
                                  control = NULL,
                                  digits = NULL,
                                  ...) {
  
  stopifnot(is.data.frame(x))
  
  p_adjust  <- match.arg(p_adjust)
  na_policy <- match.arg(na_policy)
  ratioCI   <- match.arg(ratioCI)
  ORCI      <- match.arg(ORCI)
  
  if (!is.logical(include_flags) || length(include_flags) != 1L || is.na(include_flags)) {
    stop("`include_flags` must be TRUE or FALSE.")
  }
  if (!is.logical(pairedP) || length(pairedP) != 1L || is.na(pairedP)) {
    stop("`pairedP` must be TRUE or FALSE.")
  }
  if (!is.logical(ratioP) || length(ratioP) != 1L || is.na(ratioP)) {
    stop("`ratioP` must be TRUE or FALSE.")
  }
  if (!is.logical(OR) || length(OR) != 1L || is.na(OR)) {
    stop("`OR` must be TRUE or FALSE.")
  }
  
  allowed_pairedPCI <- c("Wald", "cc", "BonettPrice", "AgrestiMin", "Newcombe", "Tango")
  pairedPCI <- match.arg(pairedPCI, choices = allowed_pairedPCI, several.ok = TRUE)
  
  if (!is.numeric(R) || length(R) != 1L || is.na(R) || R < 100) {
    stop("`R` must be a single number >= 100.")
  }
  R <- as.integer(R)
  
  meta <- attr(x, "mcnemar_table_meta", exact = TRUE)
  
  # control label
  if (is.null(control)) {
    control <- if (is.list(meta) && !is.null(meta$control)) {
      as.character(meta$control)
    } else {
      NA_character_
    }
  } else {
    control <- as.character(control)[1]
  }
  
  # treatment label from column if present, otherwise from row names
  treatment <- if ("treatment" %in% names(x)) {
    as.character(x$treatment)
  } else {
    rn <- rownames(x)
    if (is.null(rn)) rep(NA_character_, nrow(x)) else rn
  }
  
  # required cell-count columns
  req_counts <- c("a", "b", "c", "d")
  miss_counts <- setdiff(req_counts, names(x))
  if (length(miss_counts) > 0) {
    stop("`x` must contain columns a, b, c, and d. Missing: ",
         paste(miss_counts, collapse = ", "))
  }
  
  # Primary CI pair:
  # Prefer revised output convention (ci_low / ci_high).
  # Keep legacy fallbacks for older objects.
  has_pair <- function(df, low, high) {
    all(c(low, high) %in% names(df))
  }
  
  ci_pair <- if (has_pair(x, "ci_low", "ci_high")) {
    c("ci_low", "ci_high")
  } else if (has_pair(x, "ci_low_marg", "ci_high_marg")) {
    c("ci_low_marg", "ci_high_marg")
  } else if (has_pair(x, "ci_low_sim", "ci_high_sim")) {
    c("ci_low_sim", "ci_high_sim")
  } else if (has_pair(x, "ci_low_logit", "ci_high_logit")) {
    c("ci_low_logit", "ci_high_logit")
  } else {
    c(NA_character_, NA_character_)
  }
  
  ci_low  <- if (!is.na(ci_pair[1])) as.numeric(x[[ci_pair[1]]]) else rep(NA_real_, nrow(x))
  ci_high <- if (!is.na(ci_pair[2])) as.numeric(x[[ci_pair[2]]]) else rep(NA_real_, nrow(x))
  
  # Adjusted CI pair (revised convention)
  has_adjusted_ci <- has_pair(x, "ci_low_adj", "ci_high_adj")
  ci_low_adj  <- if (has_adjusted_ci) as.numeric(x$ci_low_adj)  else rep(NA_real_, nrow(x))
  ci_high_adj <- if (has_adjusted_ci) as.numeric(x$ci_high_adj) else rep(NA_real_, nrow(x))
  
  a  <- as.integer(x$a)
  b  <- as.integer(x$b)
  c_ <- as.integer(x$c)
  d  <- as.integer(x$d)
  N  <- a + b + c_ + d
  
  out_tbl <- tibble::tibble(
    control   = rep(control, nrow(x)),
    treatment = treatment,
    a = a,
    b = b,
    c = c_,
    d = d,
    N = N,
    cohens_g = if ("cohens_g" %in% names(x)) as.numeric(x$cohens_g) else NA_real_,
    ci_low = ci_low,
    ci_high = ci_high
  )
  
  # Automatically carry adjusted CI columns when present
  if (has_adjusted_ci) {
    out_tbl$ci_low_adj  <- ci_low_adj
    out_tbl$ci_high_adj <- ci_high_adj
  }
  
  # Flags
  if (isTRUE(include_flags)) {
    sim_cols <- intersect(
      c("sim_degenerate", "sim_degenerate_sim", "sim_degenerate_logit"),
      names(x)
    )
    if (length(sim_cols) > 0) {
      sim_deg <- Reduce(`|`, lapply(sim_cols, function(nm) {
        z <- as.logical(x[[nm]])
        z[is.na(z)] <- FALSE
        z
      }))
      out_tbl$sim_degenerate <- sim_deg
    }
    if ("recommend_marginal" %in% names(x)) {
      out_tbl$recommend_marginal <- as.logical(x$recommend_marginal)
    }
  }
  
  # ---------- paired-probability helpers ----------
  alpha <- 0.05
  if (is.list(meta) && !is.null(meta$ci_level) && is.finite(meta$ci_level)) {
    alpha <- 1 - as.numeric(meta$ci_level)
  }
  
  .extract_single_number <- function(obj, direct_names = character(0), regex_fallback = NULL) {
    for (nm in direct_names) {
      if (!is.null(obj[[nm]]) && is.numeric(obj[[nm]])) {
        val <- as.numeric(obj[[nm]][1])
        if (is.finite(val)) return(val)
      }
    }
    
    flat0 <- unlist(obj, recursive = TRUE, use.names = TRUE)
    if (length(flat0) == 0L || is.null(names(flat0))) {
      return(NA_real_)
    }
    
    nm   <- names(flat0)
    flat <- suppressWarnings(as.numeric(flat0))
    ok   <- is.finite(flat)
    
    if (!is.null(regex_fallback)) {
      idx <- which(ok & grepl(regex_fallback, nm, ignore.case = TRUE))
      if (length(idx) > 0L) return(flat[idx[1]])
    }
    
    NA_real_
  }
  
  .extract_ct_ci <- function(res) {
    est <- .extract_single_number(
      res,
      direct_names = c("estimate", "Estimate", "est", "delta", "Delta", "diff", "Difference", "or", "OR", "ratio", "Ratio"),
      regex_fallback = "estimate|^est$|delta|difference|diff|\\bor\\b|odds.?ratio|ratio"
    )
    low <- .extract_single_number(
      res,
      direct_names = c("CI_low", "Lower_limit", "lower_limit", "LCL", "lower", "Lower"),
      regex_fallback = "ci.*low|low.*ci|lower|lcl"
    )
    high <- .extract_single_number(
      res,
      direct_names = c("CI_high", "Upper_limit", "upper_limit", "UCL", "upper", "Upper"),
      regex_fallback = "ci.*high|high.*ci|upper|ucl"
    )
    c(estimate = est, ci_low = low, ci_high = high)
  }
  
  .pairedPCI_fun <- function(method) {
    switch(
      method,
      "Wald"        = contingencytables::Wald_CI_diff_paired_2x2,
      "cc"          = contingencytables::Wald_CI_diff_CC_paired_2x2,
      "BonettPrice" = contingencytables::Wald_CI_BonettPrice_paired_2x2,
      "AgrestiMin"  = contingencytables::Wald_CI_AgrestiMin_paired_2x2,
      "Newcombe"    = contingencytables::Newcombe_square_and_add_CI_paired_2x2,
      "Tango"       = contingencytables::Tango_asymptotic_score_CI_paired_2x2,
      stop("Unknown pairedPCI method: ", method)
    )
  }
  
  .ratioPCI_fun <- function(method) {
    switch(
      method,
      "BonettPrice"   = contingencytables::BonettPrice_hybrid_Wilson_score_CI_paired_2x2,
      "Asymptotic"    = contingencytables::Tang_asymptotic_score_CI_paired_2x2,
      "Mover"         = contingencytables::MOVER_Wilson_score_CI_paired_2x2,
      "BonettPriceCC" = contingencytables::BonettPrice_hybrid_Wilson_score_CI_CC_paired_2x2,
      "Wald"          = contingencytables::Wald_CI_ratio_paired_2x2,
      stop("Unknown ratioCI method: ", method)
    )
  }
  
  .orPCI_fun <- function(method) {
    switch(
      method,
      "Wilson"               = contingencytables::Transformed_Wilson_score_CI_paired_2x2,
      "ClopperPearson_midP"  = contingencytables::Transformed_Clopper_Pearson_midP_CI_paired_2x2,
      "ClopperPearson"       = contingencytables::Transformed_Clopper_Pearson_exact_CI_paired_2x2,
      "Blaker"               = contingencytables::Transformed_Blaker_exact_CI_paired_2x2,
      "Wald"                 = contingencytables::Wald_CI_OR_paired_2x2,
      "Wald_Laplace"         = contingencytables::Wald_CI_OR_Laplace_paired_2x2,
      stop("Unknown ORCI method: ", method)
    )
  }
  
  .sign_correct_ci <- function(raw_est, raw_low, raw_high, b, c, N) {
    delta_expected <- if (is.finite(N) && N > 0) (b - c) / N else NA_real_
    need_flip <- is.finite(raw_est) &&
      is.finite(delta_expected) &&
      delta_expected != 0 &&
      sign(raw_est) != sign(delta_expected)
    
    if (isTRUE(need_flip)) {
      c(
        estimate = -raw_est,
        ci_low   = -raw_high,
        ci_high  = -raw_low
      )
    } else {
      c(
        estimate = raw_est,
        ci_low   = raw_low,
        ci_high  = raw_high
      )
    }
  }
  
  # Reciprocal correction for paired ORs so reported OR aligns with b/c
  # under the control-row / treatment-column orientation:
  # b = upper-right (control No, treatment Yes)
  # c = lower-left (control Yes, treatment No)
  .reciprocal_correct_or_ci <- function(raw_est, raw_low, raw_high, b, c) {
    expected_or <- if (c > 0) {
      b / c
    } else if (b > 0 && c == 0) {
      Inf
    } else if (b == 0 && c > 0) {
      0
    } else {
      NA_real_
    }
    
    need_flip <- FALSE
    
    if (is.finite(expected_or) && expected_or > 0 &&
        is.finite(raw_est) && raw_est > 0) {
      d_direct <- abs(log(raw_est) - log(expected_or))
      d_flip   <- abs(log(1 / raw_est) - log(expected_or))
      need_flip <- is.finite(d_direct) && is.finite(d_flip) && (d_flip < d_direct)
    } else if (is.infinite(expected_or) && is.finite(raw_est) && raw_est == 0) {
      need_flip <- TRUE
    } else if (expected_or == 0 && is.infinite(raw_est)) {
      need_flip <- TRUE
    }
    
    if (isTRUE(need_flip)) {
      c(
        estimate = 1 / raw_est,
        ci_low   = 1 / raw_high,
        ci_high  = 1 / raw_low
      )
    } else {
      c(
        estimate = raw_est,
        ci_low   = raw_low,
        ci_high  = raw_high
      )
    }
  }
  
  .paired_maxt_ci <- function(fit, data, treatments, alpha, R, seed, na_policy) {
    if (!inherits(fit, "mcnemarControl")) {
      stop("`fit` must be an object returned by `mcnemar_control()` when p_adjust = 'bootstrap_maxt'.")
    }
    if (!is.data.frame(data)) {
      stop("`data` must be the original long-format data when p_adjust = 'bootstrap_maxt'.")
    }
    if (!requireNamespace("boot", quietly = TRUE)) {
      stop("Package 'boot' is required for p_adjust = 'bootstrap_maxt'. Please install it.")
    }
    
    s <- fit$settings
    req <- c(s$id, s$condition, s$outcome)
    miss <- setdiff(req, names(data))
    if (length(miss) > 0) {
      stop("`data` is missing required columns for bootstrap_maxt: ",
           paste(miss, collapse = ", "))
    }
    
    dat0 <- data[as.character(data[[s$condition]]) %in% c(s$control, treatments), req, drop = FALSE]
    dat0 <- dat0[stats::complete.cases(dat0), , drop = FALSE]
    key <- paste(dat0[[s$id]], dat0[[s$condition]], sep = "__")
    if (anyDuplicated(key)) {
      stop("`data` contains multiple rows per id-condition; bootstrap_maxt requires unique id-condition rows.")
    }
    
    wide <- tidyr::pivot_wider(
      dat0,
      id_cols = dplyr::all_of(s$id),
      names_from = dplyr::all_of(s$condition),
      values_from = dplyr::all_of(s$outcome)
    )
    
    needed_cols <- c(s$id, s$control, treatments)
    miss2 <- setdiff(needed_cols, names(wide))
    if (length(miss2) > 0) {
      stop("Missing condition columns after pivot for bootstrap_maxt: ",
           paste(miss2, collapse = ", "))
    }
    
    cc <- stats::complete.cases(wide[, needed_cols, drop = FALSE])
    wide <- wide[cc, , drop = FALSE]
    if (nrow(wide) < 2L) {
      stop("Need at least 2 complete subjects for p_adjust = 'bootstrap_maxt'.")
    }
    
    .coerce01 <- function(v, levels2) {
      v <- factor(v, levels = levels2)
      if (anyNA(v)) return(NULL)
      as.integer(v == levels2[2])
    }
    
    for (nm in c(s$control, treatments)) {
      y01 <- .coerce01(wide[[nm]], s$outcome_levels)
      if (is.null(y01)) {
        stop("Could not coerce outcomes to 0/1 for bootstrap_maxt; check outcome levels and missingness.")
      }
      wide[[nm]] <- y01
    }
    
    delta_hat <- vapply(treatments, function(tr) {
      ctrl <- wide[[s$control]]
      trt  <- wide[[tr]]
      nsub <- length(ctrl)
      if (nsub == 0L) return(NA_real_)
      b0 <- sum(ctrl == 0 & trt == 1)
      c0 <- sum(ctrl == 1 & trt == 0)
      (b0 - c0) / nsub
    }, numeric(1))
    
    boot_stat <- function(data0, indices) {
      dd   <- data0[indices, , drop = FALSE]
      ctrl <- dd[[s$control]]
      nsub <- length(ctrl)
      vapply(treatments, function(tr) {
        trt <- dd[[tr]]
        b0  <- sum(ctrl == 0 & trt == 1)
        c0  <- sum(ctrl == 1 & trt == 0)
        (b0 - c0) / nsub
      }, numeric(1))
    }
    
    if (!is.null(seed)) set.seed(seed)
    boot_out <- boot::boot(data = wide, statistic = boot_stat, R = R)
    Tmat <- boot_out$t
    
    if (na_policy == "omit_replicates") {
      ok_rows <- apply(Tmat, 1, function(z) all(is.finite(z)))
      Tmat2   <- Tmat[ok_rows, , drop = FALSE]
    } else {
      Tmat2 <- Tmat
    }
    
    se <- apply(Tmat2, 2, stats::sd, na.rm = TRUE)
    cols_use <- which(is.finite(se) & se > 0 & is.finite(delta_hat))
    rows_for_maxt <- rep(FALSE, nrow(Tmat2))
    if (length(cols_use) > 0L) {
      rows_for_maxt <- apply(Tmat2[, cols_use, drop = FALSE], 1, function(z) all(is.finite(z)))
    }
    
    crit <- 0
    if (length(cols_use) > 0L && sum(rows_for_maxt) > 0L) {
      Z <- sweep(Tmat2[rows_for_maxt, cols_use, drop = FALSE], 2, delta_hat[cols_use], "-")
      for (j in seq_along(cols_use)) {
        sj <- se[cols_use[j]]
        if (!is.finite(sj) || sj <= 0) {
          Z[, j] <- 0
        } else {
          Z[, j] <- Z[, j] / sj
        }
      }
      
      maxabs <- apply(Z, 1, function(z) max(abs(z), na.rm = TRUE))
      maxabs <- maxabs[is.finite(maxabs)]
      if (length(maxabs) > 0L) {
        crit <- stats::quantile(maxabs, probs = 1 - alpha, names = FALSE, type = 8, na.rm = TRUE)
      }
    }
    
    ci_low  <- delta_hat - crit * se
    ci_high <- delta_hat + crit * se
    degenerate <- (!is.finite(delta_hat)) | (!is.finite(se)) | (se == 0)
    
    list(
      Diff = stats::setNames(delta_hat, treatments),
      ci_low = stats::setNames(ci_low, treatments),
      ci_high = stats::setNames(ci_high, treatments),
      se = stats::setNames(se, treatments),
      degenerate = stats::setNames(degenerate, treatments),
      crit = crit
    )
  }
  
  .paired_maxt_ratio_ci <- function(fit, data, treatments, alpha, R, seed, na_policy) {
    if (!inherits(fit, "mcnemarControl")) {
      stop("`fit` must be an object returned by `mcnemar_control()` when p_adjust = 'bootstrap_maxt'.")
    }
    if (!is.data.frame(data)) {
      stop("`data` must be the original long-format data when p_adjust = 'bootstrap_maxt'.")
    }
    if (!requireNamespace("boot", quietly = TRUE)) {
      stop("Package 'boot' is required for p_adjust = 'bootstrap_maxt'. Please install it.")
    }
    
    s <- fit$settings
    req <- c(s$id, s$condition, s$outcome)
    miss <- setdiff(req, names(data))
    if (length(miss) > 0) {
      stop("`data` is missing required columns for bootstrap_maxt: ",
           paste(miss, collapse = ", "))
    }
    
    dat0 <- data[as.character(data[[s$condition]]) %in% c(s$control, treatments), req, drop = FALSE]
    dat0 <- dat0[stats::complete.cases(dat0), , drop = FALSE]
    key <- paste(dat0[[s$id]], dat0[[s$condition]], sep = "__")
    if (anyDuplicated(key)) {
      stop("`data` contains multiple rows per id-condition; bootstrap_maxt requires unique id-condition rows.")
    }
    
    wide <- tidyr::pivot_wider(
      dat0,
      id_cols = dplyr::all_of(s$id),
      names_from = dplyr::all_of(s$condition),
      values_from = dplyr::all_of(s$outcome)
    )
    
    needed_cols <- c(s$id, s$control, treatments)
    miss2 <- setdiff(needed_cols, names(wide))
    if (length(miss2) > 0) {
      stop("Missing condition columns after pivot for bootstrap_maxt: ",
           paste(miss2, collapse = ", "))
    }
    
    cc <- stats::complete.cases(wide[, needed_cols, drop = FALSE])
    wide <- wide[cc, , drop = FALSE]
    if (nrow(wide) < 2L) {
      stop("Need at least 2 complete subjects for p_adjust = 'bootstrap_maxt'.")
    }
    
    .coerce01 <- function(v, levels2) {
      v <- factor(v, levels = levels2)
      if (anyNA(v)) return(NULL)
      as.integer(v == levels2[2])
    }
    
    for (nm in c(s$control, treatments)) {
      y01 <- .coerce01(wide[[nm]], s$outcome_levels)
      if (is.null(y01)) {
        stop("Could not coerce outcomes to 0/1 for bootstrap_maxt; check outcome levels and missingness.")
      }
      wide[[nm]] <- y01
    }
    
    ratio_hat <- vapply(treatments, function(tr) {
      ctrl <- wide[[s$control]]
      trt  <- wide[[tr]]
      p_ctrl <- mean(ctrl == 1)
      p_trt  <- mean(trt == 1)
      if (!is.finite(p_ctrl) || !is.finite(p_trt) || p_ctrl <= 0 || p_trt <= 0) return(NA_real_)
      p_trt / p_ctrl
    }, numeric(1))
    
    eta_hat <- rep(NA_real_, length(treatments))
    names(eta_hat) <- treatments
    ok_obs <- is.finite(ratio_hat) & (ratio_hat > 0)
    eta_hat[ok_obs] <- log(ratio_hat[ok_obs])
    
    boot_stat <- function(data0, indices) {
      dd   <- data0[indices, , drop = FALSE]
      ctrl <- dd[[s$control]]
      vapply(treatments, function(tr) {
        trt <- dd[[tr]]
        p_ctrl <- mean(ctrl == 1)
        p_trt  <- mean(trt == 1)
        if (!is.finite(p_ctrl) || !is.finite(p_trt) || p_ctrl <= 0 || p_trt <= 0) return(NA_real_)
        log(p_trt / p_ctrl)
      }, numeric(1))
    }
    
    if (!is.null(seed)) set.seed(seed)
    boot_out <- boot::boot(data = wide, statistic = boot_stat, R = R)
    Tmat <- boot_out$t
    
    if (na_policy == "omit_replicates") {
      ok_rows <- apply(Tmat, 1, function(z) all(is.finite(z)))
      Tmat2 <- Tmat[ok_rows, , drop = FALSE]
    } else {
      Tmat2 <- Tmat
    }
    
    se <- apply(Tmat2, 2, stats::sd, na.rm = TRUE)
    cols_use <- which(is.finite(se) & se > 0 & is.finite(eta_hat))
    rows_for_maxt <- rep(FALSE, nrow(Tmat2))
    if (length(cols_use) > 0L) {
      rows_for_maxt <- apply(Tmat2[, cols_use, drop = FALSE], 1, function(z) all(is.finite(z)))
    }
    
    crit <- 0
    if (length(cols_use) > 0L && sum(rows_for_maxt) > 0L) {
      Z <- sweep(Tmat2[rows_for_maxt, cols_use, drop = FALSE], 2, eta_hat[cols_use], "-")
      for (j in seq_along(cols_use)) {
        sj <- se[cols_use[j]]
        if (!is.finite(sj) || sj <= 0) {
          Z[, j] <- 0
        } else {
          Z[, j] <- Z[, j] / sj
        }
      }
      
      maxabs <- apply(Z, 1, function(z) max(abs(z), na.rm = TRUE))
      maxabs <- maxabs[is.finite(maxabs)]
      if (length(maxabs) > 0L) {
        crit <- stats::quantile(maxabs, probs = 1 - alpha, names = FALSE, type = 8, na.rm = TRUE)
      }
    }
    
    ci_eta_low  <- eta_hat - crit * se
    ci_eta_high <- eta_hat + crit * se
    ci_low  <- exp(ci_eta_low)
    ci_high <- exp(ci_eta_high)
    
    degenerate <- (!is.finite(eta_hat)) | (!is.finite(se)) | (se == 0)
    ci_low[degenerate]  <- NA_real_
    ci_high[degenerate] <- NA_real_
    
    list(
      Ratio = stats::setNames(ratio_hat, treatments),
      ci_low = stats::setNames(ci_low, treatments),
      ci_high = stats::setNames(ci_high, treatments),
      se = stats::setNames(se, treatments),
      degenerate = stats::setNames(degenerate, treatments),
      crit = crit
    )
  }
  
  .paired_maxt_or_ci <- function(fit, data, treatments, alpha, R, seed, na_policy) {
    if (!inherits(fit, "mcnemarControl")) {
      stop("`fit` must be an object returned by `mcnemar_control()` when p_adjust = 'bootstrap_maxt'.")
    }
    if (!is.data.frame(data)) {
      stop("`data` must be the original long-format data when p_adjust = 'bootstrap_maxt'.")
    }
    if (!requireNamespace("boot", quietly = TRUE)) {
      stop("Package 'boot' is required for p_adjust = 'bootstrap_maxt'. Please install it.")
    }
    
    s <- fit$settings
    req <- c(s$id, s$condition, s$outcome)
    miss <- setdiff(req, names(data))
    if (length(miss) > 0) {
      stop("`data` is missing required columns for bootstrap_maxt: ",
           paste(miss, collapse = ", "))
    }
    
    dat0 <- data[as.character(data[[s$condition]]) %in% c(s$control, treatments), req, drop = FALSE]
    dat0 <- dat0[stats::complete.cases(dat0), , drop = FALSE]
    key <- paste(dat0[[s$id]], dat0[[s$condition]], sep = "__")
    if (anyDuplicated(key)) {
      stop("`data` contains multiple rows per id-condition; bootstrap_maxt requires unique id-condition rows.")
    }
    
    wide <- tidyr::pivot_wider(
      dat0,
      id_cols = dplyr::all_of(s$id),
      names_from = dplyr::all_of(s$condition),
      values_from = dplyr::all_of(s$outcome)
    )
    
    needed_cols <- c(s$id, s$control, treatments)
    miss2 <- setdiff(needed_cols, names(wide))
    if (length(miss2) > 0) {
      stop("Missing condition columns after pivot for bootstrap_maxt: ",
           paste(miss2, collapse = ", "))
    }
    
    cc <- stats::complete.cases(wide[, needed_cols, drop = FALSE])
    wide <- wide[cc, , drop = FALSE]
    if (nrow(wide) < 2L) {
      stop("Need at least 2 complete subjects for p_adjust = 'bootstrap_maxt'.")
    }
    
    .coerce01 <- function(v, levels2) {
      v <- factor(v, levels = levels2)
      if (anyNA(v)) return(NULL)
      as.integer(v == levels2[2])
    }
    
    for (nm in c(s$control, treatments)) {
      y01 <- .coerce01(wide[[nm]], s$outcome_levels)
      if (is.null(y01)) {
        stop("Could not coerce outcomes to 0/1 for bootstrap_maxt; check outcome levels and missingness.")
      }
      wide[[nm]] <- y01
    }
    
    or_hat <- vapply(treatments, function(tr) {
      ctrl <- wide[[s$control]]
      trt  <- wide[[tr]]
      b0 <- sum(ctrl == 0 & trt == 1)
      c0 <- sum(ctrl == 1 & trt == 0)
      if (b0 <= 0 || c0 <= 0) return(NA_real_)
      b0 / c0
    }, numeric(1))
    
    eta_hat <- rep(NA_real_, length(treatments))
    names(eta_hat) <- treatments
    ok_obs <- is.finite(or_hat) & (or_hat > 0)
    eta_hat[ok_obs] <- log(or_hat[ok_obs])
    
    boot_stat <- function(data0, indices) {
      dd   <- data0[indices, , drop = FALSE]
      ctrl <- dd[[s$control]]
      vapply(treatments, function(tr) {
        trt <- dd[[tr]]
        b0  <- sum(ctrl == 0 & trt == 1)
        c0  <- sum(ctrl == 1 & trt == 0)
        if (b0 <= 0 || c0 <= 0) return(NA_real_)
        log(b0 / c0)
      }, numeric(1))
    }
    
    if (!is.null(seed)) set.seed(seed)
    boot_out <- boot::boot(data = wide, statistic = boot_stat, R = R)
    Tmat <- boot_out$t
    
    if (na_policy == "omit_replicates") {
      ok_rows <- apply(Tmat, 1, function(z) all(is.finite(z)))
      Tmat2 <- Tmat[ok_rows, , drop = FALSE]
    } else {
      Tmat2 <- Tmat
    }
    
    se <- apply(Tmat2, 2, stats::sd, na.rm = TRUE)
    cols_use <- which(is.finite(se) & se > 0 & is.finite(eta_hat))
    rows_for_maxt <- rep(FALSE, nrow(Tmat2))
    if (length(cols_use) > 0L) {
      rows_for_maxt <- apply(Tmat2[, cols_use, drop = FALSE], 1, function(z) all(is.finite(z)))
    }
    
    crit <- 0
    if (length(cols_use) > 0L && sum(rows_for_maxt) > 0L) {
      Z <- sweep(Tmat2[rows_for_maxt, cols_use, drop = FALSE], 2, eta_hat[cols_use], "-")
      for (j in seq_along(cols_use)) {
        sj <- se[cols_use[j]]
        if (!is.finite(sj) || sj <= 0) {
          Z[, j] <- 0
        } else {
          Z[, j] <- Z[, j] / sj
        }
      }
      
      maxabs <- apply(Z, 1, function(z) max(abs(z), na.rm = TRUE))
      maxabs <- maxabs[is.finite(maxabs)]
      if (length(maxabs) > 0L) {
        crit <- stats::quantile(maxabs, probs = 1 - alpha, names = FALSE, type = 8, na.rm = TRUE)
      }
    }
    
    ci_eta_low  <- eta_hat - crit * se
    ci_eta_high <- eta_hat + crit * se
    ci_low  <- exp(ci_eta_low)
    ci_high <- exp(ci_eta_high)
    
    degenerate <- (!is.finite(eta_hat)) | (!is.finite(se)) | (se == 0)
    ci_low[degenerate]  <- NA_real_
    ci_high[degenerate] <- NA_real_
    
    list(
      OR = stats::setNames(or_hat, treatments),
      ci_low = stats::setNames(ci_low, treatments),
      ci_high = stats::setNames(ci_high, treatments),
      se = stats::setNames(se, treatments),
      degenerate = stats::setNames(degenerate, treatments),
      crit = crit
    )
  }
  
  # ---------- second tibble: differences between marginal proportions ----------
  paired_tbl <- NULL
  if (isTRUE(pairedP)) {
    K <- nrow(x)
    Diff_raw <- ifelse(N > 0, (b - c_) / N, NA_real_)
    
    paired_tbl <- tibble::tibble(
      control   = rep(control, nrow(x)),
      treatment = treatment,
      Diff = Diff_raw
    )
    
    # marginal pairedPCI CIs + method-matched estimate
    for (m in pairedPCI) {
      fn <- .pairedPCI_fun(m)
      
      vals <- lapply(seq_len(nrow(x)), function(i) {
        tab <- matrix(c(a[i], b[i], c_[i], d[i]), nrow = 2L, byrow = TRUE)
        res_ci <- fn(tab, alpha = alpha)
        ext <- .extract_ct_ci(res_ci)
        corr <- .sign_correct_ci(
          raw_est  = unname(ext["estimate"]),
          raw_low  = unname(ext["ci_low"]),
          raw_high = unname(ext["ci_high"]),
          b = b[i],
          c = c_[i],
          N = N[i]
        )
        stats::setNames(
          c(unname(corr["estimate"]), unname(corr["ci_low"]), unname(corr["ci_high"])),
          c("Diff", "ci_low", "ci_high")
        )
      })
      
      mat <- do.call(rbind, vals)
      mat <- as.matrix(mat)
      
      if (length(pairedPCI) == 1L) {
        paired_tbl$Diff    <- as.numeric(mat[, "Diff"])
        paired_tbl$ci_low  <- as.numeric(mat[, "ci_low"])
        paired_tbl$ci_high <- as.numeric(mat[, "ci_high"])
      } else {
        paired_tbl[[paste0("ci_low_", m)]]  <- as.numeric(mat[, "ci_low"])
        paired_tbl[[paste0("ci_high_", m)]] <- as.numeric(mat[, "ci_high"])
      }
    }
    
    if (p_adjust == "bonferroni") {
      alpha_adj <- alpha / K
      for (m in pairedPCI) {
        fn <- .pairedPCI_fun(m)
        
        vals_adj <- lapply(seq_len(nrow(x)), function(i) {
          tab <- matrix(c(a[i], b[i], c_[i], d[i]), nrow = 2L, byrow = TRUE)
          res_ci <- fn(tab, alpha = alpha_adj)
          ext <- .extract_ct_ci(res_ci)
          corr <- .sign_correct_ci(
            raw_est  = unname(ext["estimate"]),
            raw_low  = unname(ext["ci_low"]),
            raw_high = unname(ext["ci_high"]),
            b = b[i],
            c = c_[i],
            N = N[i]
          )
          stats::setNames(
            c(unname(corr["ci_low"]), unname(corr["ci_high"])),
            c("ci_low_adj", "ci_high_adj")
          )
        })
        
        mat_adj <- do.call(rbind, vals_adj)
        mat_adj <- as.matrix(mat_adj)
        
        if (length(pairedPCI) == 1L) {
          paired_tbl$ci_low_adj  <- as.numeric(mat_adj[, "ci_low_adj"])
          paired_tbl$ci_high_adj <- as.numeric(mat_adj[, "ci_high_adj"])
        } else {
          paired_tbl[[paste0("ci_low_adj_", m)]]  <- as.numeric(mat_adj[, "ci_low_adj"])
          paired_tbl[[paste0("ci_high_adj_", m)]] <- as.numeric(mat_adj[, "ci_high_adj"])
        }
      }
    }
    
    if (p_adjust == "bootstrap_maxt") {
      maxt <- .paired_maxt_ci(
        fit = fit,
        data = data,
        treatments = treatment,
        alpha = alpha,
        R = R,
        seed = seed,
        na_policy = na_policy
      )
      paired_tbl$ci_low_adj  <- as.numeric(maxt$ci_low[treatment])
      paired_tbl$ci_high_adj <- as.numeric(maxt$ci_high[treatment])
      if (isTRUE(include_flags)) {
        paired_tbl$pairedP_sim_degenerate <- as.logical(maxt$degenerate[treatment])
      }
    }
  }
  
  # ---------- third tibble: ratio of marginal proportions ----------
  ratio_tbl <- NULL
  if (isTRUE(ratioP)) {
    K <- nrow(x)
    
    # In this table orientation:
    # P(Control = Yes) = (c + d) / N
    # P(Treatment = Yes) = (b + d) / N
    ctrl_prop <- ifelse(N > 0, (c_ + d) / N, NA_real_)
    trt_prop  <- ifelse(N > 0, (b + d) / N, NA_real_)
    ratio_raw <- ifelse(ctrl_prop > 0, trt_prop / ctrl_prop, NA_real_)
    
    fn_ratio <- .ratioPCI_fun(ratioCI)
    
    vals_ratio <- lapply(seq_len(nrow(x)), function(i) {
      # Reoriented for ratio methods so that:
      # row 1 total = P(Treatment = Yes) = b + d
      # col 1 total = P(Control = Yes) = c + d
      tab <- matrix(c(d[i], b[i],
                      c_[i], a[i]), nrow = 2L, byrow = TRUE)
      
      res_ci <- fn_ratio(tab, alpha = alpha)
      ext <- .extract_ct_ci(res_ci)
      est_i  <- unname(ext["estimate"])
      low_i  <- unname(ext["ci_low"])
      high_i <- unname(ext["ci_high"])
      
      raw_i <- if (is.finite(c_[i] + d[i]) && (c_[i] + d[i]) > 0) {
        (b[i] + d[i]) / (c_[i] + d[i])
      } else {
        NA_real_
      }
      
      if (is.finite(est_i) && is.finite(raw_i) && raw_i > 0 && est_i > 0) {
        if (isTRUE(all.equal(est_i, 1 / raw_i, tolerance = 1e-8))) {
          est_i <- 1 / est_i
          low_old  <- low_i
          high_old <- high_i
          low_i  <- 1 / high_old
          high_i <- 1 / low_old
        }
      }
      
      stats::setNames(
        c(est_i, low_i, high_i),
        c("Ratio", "ci_low", "ci_high")
      )
    })
    
    mat_ratio <- do.call(rbind, vals_ratio)
    mat_ratio <- as.matrix(mat_ratio)
    
    ratio_tbl <- tibble::tibble(
      control   = ctrl_prop,
      treatment = trt_prop,
      Ratio = as.numeric(mat_ratio[, "Ratio"]),
      ci_low = as.numeric(mat_ratio[, "ci_low"]),
      ci_high = as.numeric(mat_ratio[, "ci_high"])
    )
    
    if (p_adjust == "bonferroni") {
      alpha_adj <- alpha / K
      
      vals_ratio_adj <- lapply(seq_len(nrow(x)), function(i) {
        tab <- matrix(c(d[i], b[i],
                        c_[i], a[i]), nrow = 2L, byrow = TRUE)
        
        res_ci <- fn_ratio(tab, alpha = alpha_adj)
        ext <- .extract_ct_ci(res_ci)
        est_i  <- unname(ext["estimate"])
        low_i  <- unname(ext["ci_low"])
        high_i <- unname(ext["ci_high"])
        
        raw_i <- if (is.finite(c_[i] + d[i]) && (c_[i] + d[i]) > 0) {
          (b[i] + d[i]) / (c_[i] + d[i])
        } else {
          NA_real_
        }
        
        if (is.finite(est_i) && is.finite(raw_i) && raw_i > 0 && est_i > 0) {
          if (isTRUE(all.equal(est_i, 1 / raw_i, tolerance = 1e-8))) {
            low_old  <- low_i
            high_old <- high_i
            low_i  <- 1 / high_old
            high_i <- 1 / low_old
          }
        }
        
        stats::setNames(
          c(low_i, high_i),
          c("ci_low_adj", "ci_high_adj")
        )
      })
      
      mat_ratio_adj <- do.call(rbind, vals_ratio_adj)
      mat_ratio_adj <- as.matrix(mat_ratio_adj)
      ratio_tbl$ci_low_adj  <- as.numeric(mat_ratio_adj[, "ci_low_adj"])
      ratio_tbl$ci_high_adj <- as.numeric(mat_ratio_adj[, "ci_high_adj"])
    }
    
    if (p_adjust == "bootstrap_maxt") {
      maxt_ratio <- .paired_maxt_ratio_ci(
        fit = fit,
        data = data,
        treatments = treatment,
        alpha = alpha,
        R = R,
        seed = seed,
        na_policy = na_policy
      )
      ratio_tbl$ci_low_adj  <- as.numeric(maxt_ratio$ci_low[treatment])
      ratio_tbl$ci_high_adj <- as.numeric(maxt_ratio$ci_high[treatment])
      if (isTRUE(include_flags)) {
        ratio_tbl$ratioP_sim_degenerate <- as.logical(maxt_ratio$degenerate[treatment])
      }
    }
  }
  
  # ---------- fourth tibble: conditional odds ratios ----------
  or_tbl <- NULL
  if (isTRUE(OR)) {
    K <- nrow(x)
    fn_or <- .orPCI_fun(ORCI)
    
    vals_or <- lapply(seq_len(nrow(x)), function(i) {
      # Keep native table orientation used throughout this summary:
      #      Treatment: No   Treatment: Yes
      # C No      a               b
      # C Yes     c               d
      # and then direction-correct if the package's OR orientation is reciprocal.
      tab <- matrix(c(a[i], b[i],
                      c_[i], d[i]), nrow = 2L, byrow = TRUE)
      
      res_ci <- fn_or(tab, alpha = alpha)
      ext <- .extract_ct_ci(res_ci)
      
      corr <- .reciprocal_correct_or_ci(
        raw_est  = unname(ext["estimate"]),
        raw_low  = unname(ext["ci_low"]),
        raw_high = unname(ext["ci_high"]),
        b = b[i],
        c = c_[i]
      )
      
      stats::setNames(
        c(unname(corr["estimate"]),
          unname(corr["ci_low"]),
          unname(corr["ci_high"])),
        c("OR", "ci_low", "ci_high")
      )
    })
    
    mat_or <- do.call(rbind, vals_or)
    mat_or <- as.matrix(mat_or)
    
    or_tbl <- tibble::tibble(
      control   = rep(control, nrow(x)),
      treatment = treatment,
      OR = as.numeric(mat_or[, "OR"]),
      ci_low = as.numeric(mat_or[, "ci_low"]),
      ci_high = as.numeric(mat_or[, "ci_high"])
    )
    
    if (p_adjust == "bonferroni") {
      alpha_adj <- alpha / K
      
      vals_or_adj <- lapply(seq_len(nrow(x)), function(i) {
        tab <- matrix(c(a[i], b[i],
                        c_[i], d[i]), nrow = 2L, byrow = TRUE)
        
        res_ci <- fn_or(tab, alpha = alpha_adj)
        ext <- .extract_ct_ci(res_ci)
        
        corr <- .reciprocal_correct_or_ci(
          raw_est  = unname(ext["estimate"]),
          raw_low  = unname(ext["ci_low"]),
          raw_high = unname(ext["ci_high"]),
          b = b[i],
          c = c_[i]
        )
        
        stats::setNames(
          c(unname(corr["ci_low"]), unname(corr["ci_high"])),
          c("ci_low_adj", "ci_high_adj")
        )
      })
      
      mat_or_adj <- do.call(rbind, vals_or_adj)
      mat_or_adj <- as.matrix(mat_or_adj)
      or_tbl$ci_low_adj  <- as.numeric(mat_or_adj[, "ci_low_adj"])
      or_tbl$ci_high_adj <- as.numeric(mat_or_adj[, "ci_high_adj"])
    }
    
    if (p_adjust == "bootstrap_maxt") {
      maxt_or <- .paired_maxt_or_ci(
        fit = fit,
        data = data,
        treatments = treatment,
        alpha = alpha,
        R = R,
        seed = seed,
        na_policy = na_policy
      )
      or_tbl$ci_low_adj  <- as.numeric(maxt_or$ci_low[treatment])
      or_tbl$ci_high_adj <- as.numeric(maxt_or$ci_high[treatment])
      
      if (isTRUE(include_flags)) {
        or_tbl$OR_sim_degenerate <- as.logical(maxt_or$degenerate[treatment])
      }
    }
  }
  
  # Rounding
  if (!is.null(digits)) {
    if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) || digits < 0) {
      stop("`digits` must be a single non-negative number.")
    }
    digits <- as.integer(digits)
    
    for (nm in intersect(c("cohens_g", "ci_low", "ci_high", "ci_low_adj", "ci_high_adj"), names(out_tbl))) {
      if (is.numeric(out_tbl[[nm]])) {
        out_tbl[[nm]] <- round(out_tbl[[nm]], digits)
      }
    }
    
    if (!is.null(paired_tbl)) {
      num_cols <- names(paired_tbl)[vapply(paired_tbl, is.numeric, logical(1))]
      for (nm in num_cols) {
        paired_tbl[[nm]] <- round(paired_tbl[[nm]], digits)
      }
    }
    
    if (!is.null(ratio_tbl)) {
      num_cols <- names(ratio_tbl)[vapply(ratio_tbl, is.numeric, logical(1))]
      for (nm in num_cols) {
        ratio_tbl[[nm]] <- round(ratio_tbl[[nm]], digits)
      }
    }
    
    if (!is.null(or_tbl)) {
      num_cols <- names(or_tbl)[vapply(or_tbl, is.numeric, logical(1))]
      for (nm in num_cols) {
        or_tbl[[nm]] <- round(or_tbl[[nm]], digits)
      }
    }
  }
  
  overview <- list(
    control = control,
    comparisons = nrow(out_tbl),
    methodCI = if (is.list(meta) && !is.null(meta$methodCI)) meta$methodCI else NA_character_,
    ci_level = if (is.list(meta) && !is.null(meta$ci_level)) meta$ci_level else NA_real_,
    adjusted_ci = has_adjusted_ci,
    include_flags = isTRUE(include_flags),
    pairedP = isTRUE(pairedP),
    pairedPCI = pairedPCI,
    ratioP = isTRUE(ratioP),
    ratioCI = ratioCI,
    OR = isTRUE(OR),
    ORCI = ORCI,
    p_adjust = p_adjust
  )
  
  out <- list(
    overview = overview,
    table = out_tbl,
    paired_probabilities = paired_tbl,
    ratio_probabilities = ratio_tbl,
    odds_ratios = or_tbl,
    original = x
  )
  class(out) <- "summary.mcnemar_table"
  out
}

#' @rdname mcnemar_table_summary
#' @export
summary.mcnemar_table <- function(object,
                                  include_flags = FALSE,
                                  pairedP = FALSE,
                                  pairedPCI = c("BonettPrice", "AgrestiMin", "Newcombe", "Tango", "Wald", "cc"),
                                  ratioP = FALSE,
                                  ratioCI = c("BonettPrice", "Asymptotic", "Mover", "BonettPriceCC", "Wald"),
                                  OR = FALSE,
                                  ORCI = c("Wilson", "ClopperPearson_midP", "ClopperPearson", "Blaker", "Wald", "Wald_Laplace"),
                                  p_adjust = c("none", "bonferroni", "bootstrap_maxt"),
                                  fit = NULL,
                                  data = NULL,
                                  R = 5000,
                                  seed = NULL,
                                  na_policy = c("omit_replicates", "pairwise"),
                                  control = NULL,
                                  digits = NULL,
                                  ...) {
  
  mcnemar_table_summary(
    x = object,
    include_flags = include_flags,
    pairedP = pairedP,
    pairedPCI = pairedPCI,
    ratioP = ratioP,
    ratioCI = ratioCI,
    OR = OR,
    ORCI = ORCI,
    p_adjust = p_adjust,
    fit = fit,
    data = data,
    R = R,
    seed = seed,
    na_policy = na_policy,
    control = control,
    digits = digits,
    ...
  )
}