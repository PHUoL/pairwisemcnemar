#' Example long-format paired binary data
#'
#' @description
#' Returns a toy dataset in long format suitable for \code{mcnemar_control()}.
#' Each `id` has a binary `outcome` for the `Control` and three treatments.
#'
#' @param n Number of subjects.
#' @param seed Random seed.
#'
#' @return A tibble with columns `id`, `condition`, `outcome`.
#' @export
mcnemar_example_long <- function(n = 40, seed = 1) {
  set.seed(seed)
  id <- seq_len(n)
  conditions <- c("Control", "TrtA", "TrtB", "TrtC")

  base_p <- stats::plogis(stats::rnorm(n, 0, 0.7))
  shifts <- c(Control = 0, TrtA = 0.30, TrtB = 0.10, TrtC = -0.05)

  out <- lapply(conditions, function(cc) {
    p <- pmin(pmax(base_p + shifts[[cc]], 0.02), 0.98)
    tibble::tibble(
      id = id,
      condition = cc,
      outcome = ifelse(stats::runif(n) < p, "Yes", "No")
    )
  })

  dat <- dplyr::bind_rows(out)
  dat <- dplyr::mutate(
    dat,
    condition = factor(.data$condition, levels = conditions),
    outcome   = factor(.data$outcome, levels = c("No", "Yes"))
  )
  dat
}
