#' Random-walk Metropolis sampler for theta (Rcpp version)
#'
#' Random-walk Metropolis sampler for \eqn{\theta \in (0,1)} in the genetic linkage model
#' (Exercise 11.8). Multinomial probabilities are
#' \eqn{(1/2 + \theta/4,\ (1-\theta)/4,\ (1-\theta)/4,\ \theta/4)}.
#' Proposal: \eqn{\theta' = \theta + U(-w,w)}. Proposals outside (0,1) are rejected.
#'
#' @details
#' This function calls a compiled C++ routine (\code{rw_mh_theta_cpp_impl}) for speed.
#' The returned draws are post burn-in and optionally thinned.
#'
#' Practical guidelines:
#' \itemize{
#'   \item Tune \code{w} so that \code{accept_rate} is roughly 0.2--0.6.
#'   \item Use a burn-in (e.g., 500--2000) before computing posterior summaries.
#'   \item Thinning is optional; only use it when you need to reduce storage.
#' }
#'
#' @param k Integer vector of length 4 (counts).
#' @param m Integer, total iterations (including burn-in).
#' @param w Numeric, proposal half-width.
#' @param x0 Numeric, initial value in (0,1).
#' @param burn Integer, burn-in iterations to drop from the start.
#' @param thin Integer, keep every \code{thin}-th draw after burn-in.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{chain}: numeric vector of posterior draws after burn-in/thinning.
#'   \item \code{accept_rate}: acceptance rate computed on the full chain (before thinning).
#' }
#'
#' @examples
#' k <- c(125, 18, 20, 34)
#' out <- rw_mh_theta_cpp(k, m = 5000, w = 0.10, x0 = 0.50, burn = 500, thin = 2)
#' mean(out$chain)
#' out$accept_rate
#'
#' @export
rw_mh_theta_cpp <- function(k, m = 10000L, w = 0.10, x0 = 0.50,
                           burn = 0L, thin = 1L) {
  if (!exists("rw_mh_theta_cpp_impl", mode = "function")) {
    stop("Compiled function rw_mh_theta_cpp_impl not found. Did you run usethis::use_rcpp() and Rcpp::compileAttributes()?")
  }
  k <- as.integer(k)
  m <- as.integer(m)
  burn <- as.integer(burn)
  thin <- as.integer(thin)

  if (length(k) != 4L) stop("k must be length 4.")
  if (m < 2L) stop("m must be >= 2.")
  if (burn < 0L || burn >= m) stop("burn must be in [0, m-1].")
  if (thin < 1L) stop("thin must be >= 1.")
  if (!is.finite(w) || w <= 0) stop("w must be positive.")
  if (!is.finite(x0) || x0 <= 0 || x0 >= 1) stop("x0 must be in (0,1).")

  res <- rw_mh_theta_cpp_impl(k, m, as.numeric(w), as.numeric(x0))
  x <- res$chain
  if (burn > 0L) x <- x[(burn+1L):length(x)]
  if (thin > 1L) x <- x[seq(1L, length(x), by = thin)]

  list(chain = x, accept_rate = as.numeric(res$accept_rate))
}
