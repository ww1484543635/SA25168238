#' HMC sampler for theta (Rcpp)
#'
#' This function generates samples of \eqn{\theta \in (0,1)} for the genetic linkage
#' example (Exercise 11.8). The posterior is proportional to a multinomial likelihood
#' with probabilities
#' \eqn{(1/2 + \theta/4,\ (1-\theta)/4,\ (1-\theta)/4,\ \theta/4)}.
#'
#' Compared with a simple random-walk MH, HMC uses gradient information to propose
#' larger moves with a high acceptance rate. Internally we use a logit transform
#' \eqn{z=\log(\theta/(1-\theta))} so the algorithm runs on the real line.
#'
#' @param k Integer vector of length 4, the observed counts \code{(k1,k2,k3,k4)}.
#' @param m Integer, total number of iterations (including burn-in).
#' @param w Numeric, HMC step size (epsilon). Smaller values are more stable but slower.
#' @param x0 Numeric, initial value of theta in (0,1).
#' @param L Integer, number of leapfrog steps per iteration.
#' @param burn Integer, number of burn-in iterations removed from the start.
#' @param thin Integer, keep every \code{thin}-th draw after burn-in.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{chain}: numeric vector of posterior draws after burn-in and thinning.
#'   \item \code{accept_rate}: acceptance rate of HMC proposals.
#' }
#'
#' @examples
#' k <- c(125, 18, 20, 34)
#' out <- rw_mh_theta_cpp(k, m = 2000, w = 0.03, L = 15, burn = 200, thin = 2)
#' mean(out$chain)
#' out$accept_rate
#'
#' @export

rw_mh_theta_cpp <- function(k, m = 10000L, w = 0.02, x0 = 0.50, L = 10L,
                            burn = 0L, thin = 1L) {
  if (!exists("rw_mh_theta_cpp_impl", mode = "function")) {
    stop("Compiled function rw_mh_theta_cpp_impl not found. Did you run usethis::use_rcpp() and Rcpp::compileAttributes()?")
  }
  k <- as.integer(k)
  m <- as.integer(m)
  burn <- as.integer(burn)
  thin <- as.integer(thin)
  L <- as.integer(L)
  if (L < 1L) stop("L must be >= 1.")

  if (length(k) != 4L) stop("k must be length 4.")
  if (m < 2L) stop("m must be >= 2.")
  if (burn < 0L || burn >= m) stop("burn must be in [0, m-1].")
  if (thin < 1L) stop("thin must be >= 1.")
  if (!is.finite(w) || w <= 0) stop("w must be positive.")
  if (!is.finite(x0) || x0 <= 0 || x0 >= 1) stop("x0 must be in (0,1).")

  res <- rw_mh_theta_cpp_impl(k, m, as.numeric(w), as.numeric(x0), as.integer(L))
  x <- res$chain
  if (burn > 0L) x <- x[(burn+1L):length(x)]
  if (thin > 1L) x <- x[seq(1L, length(x), by = thin)]

  list(chain = x, accept_rate = as.numeric(res$accept_rate))
}
