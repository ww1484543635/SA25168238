#' Two-sample Cramér–von Mises test (Eq. 10.14) via permutation
#'
#' Implements the two-sample Cramér–von Mises statistic exactly as in Eq. (10.14),
#' using pooled ranks with average ties, and computes a right-tail permutation p-value
#' by re-labeling the pooled sample.
#'
#' @details
#' This is a distribution-free, rank-based two-sample test. It is useful when comparing
#' two samples without assuming normality. The permutation procedure approximates the
#' null distribution by repeatedly shuffling group labels.
#'
#' Practical guidelines:
#' \itemize{
#'   \item Use \code{B >= 999} for quick checks; \code{B >= 4999} for more stable p-values.
#'   \item If there are many ties, the rank-based statistic still works but the p-value can be noisier.
#' }
#'
#' @param x,y Numeric vectors (two samples).
#' @param B Integer, number of permutations.
#' @param seed Integer seed for reproducibility. Use \code{NULL} to not set seed.
#' @param return_perm Logical, whether to return the permutation statistics.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{stat}: observed CvM statistic (Eq. 10.14).
#'   \item \code{p_value}: right-tail permutation p-value.
#'   \item \code{perm_stat}: permutation statistics of length \code{B} (only if \code{return_perm = TRUE}).
#' }
#'
#' @examples
#' set.seed(2025)
#' x <- rnorm(20)
#' y <- rnorm(25, 0.5)
#' cvm_perm_test_1014(x, y, B = 999)
#'
#' @export
cvm_perm_test_1014 <- function(x, y, B = 2999L, seed = 2025L, return_perm = FALSE) {
  x <- as.numeric(x); y <- as.numeric(y)
  n <- length(x); m <- length(y)
  if (n < 1L || m < 1L) stop("x and y must be non-empty.")
  B <- as.integer(B)
  if (B < 199L) warning("B is small; p-value may be noisy.")

  if (!is.null(seed)) set.seed(as.integer(seed))

  # (10.14) statistic using pooled ranks with average ties
  cvm_w2_1014 <- function(x, y) {
    n <- length(x); m <- length(y)
    z <- c(x, y)
    r_all <- rank(z, ties.method = "average")
    r_x <- r_all[seq_len(n)]
    r_y <- r_all[n + seq_len(m)]

    U  <- n * sum((r_x - (1:n))^2) + m * sum((r_y - (1:m))^2)
    U/(n*m*(n+m)) - (4*n*m - 1)/(6*(n+m))
  }

  stat0 <- as.numeric(cvm_w2_1014(x, y))

  z <- c(x, y)
  N <- n + m
  perm <- numeric(B)

  for (b in seq_len(B)) {
    ix <- sample.int(N, N, replace = FALSE)
    perm[b] <- cvm_w2_1014(z[ix[1:n]], z[ix[(n+1):N]])
  }

  # right-tail permutation p-value (as in HW)
  p <- (1 + sum(perm >= stat0)) / (B + 1)

  out <- list(stat = stat0, p_value = p)
  if (isTRUE(return_perm)) out$perm_stat <- perm
  out
}
