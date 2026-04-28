#' Numerical Constants for TAC/WSC/MPC Metric Equivalence
#'
#' Computes the explicit finite-sample constants K_TAC(gamma, q), K_MPC(gamma)
#' and the unified two-sided envelope (c_L, c_U) from Theorems 5-8 of the
#' companion paper.
#'
#' @param gamma Mixing decay rate; must satisfy \code{gamma > 1}.
#' @param q Moment order used in Ibragimov-Rio covariance bound; must satisfy
#'   \code{q > 2} and \code{gamma * (1 - 2/q) > 1}.
#' @param n Sample size used for the two-sided envelope.
#' @param C_alpha Strong-mixing constant (default 1).
#'
#' @return A list with components \code{K_TAC}, \code{K_MPC}, \code{c_L},
#'   \code{c_U}, \code{b_optimal}, \code{m_optimal}, and \code{regime}
#'   (a short label describing which method is favored).
#'
#' @examples
#' htdv_equivalence_constants(gamma = 2, q = 6, n = 60)
#'
#' @export
htdv_equivalence_constants <- function(gamma, q = 6, n = 100, C_alpha = 1) {
  .assert_positive_scalar(gamma, "gamma")
  .assert_positive_scalar(q, "q")
  .assert_positive_scalar(n, "n")
  .assert_positive_scalar(C_alpha, "C_alpha")
  if (gamma <= 1) {
    stop("'gamma' must be > 1 for summable mixing.", call. = FALSE)
  }
  if (q <= 2) {
    stop("'q' must be > 2.", call. = FALSE)
  }
  s <- gamma * (1 - 2 / q)
  if (s <= 1) {
    stop("Condition gamma*(1 - 2/q) > 1 violated; increase gamma or q.",
         call. = FALSE)
  }
  K_TAC <- 1 + 16 * C_alpha^(1 - 2 / q) * .zeta_riemann(s)
  K_MPC <- (gamma / (gamma - 1))^2 * (1 + 4 / (gamma - 1))
  c_U <- sqrt(max(K_TAC, K_MPC))
  penalty <- 8 * C_alpha^(1 - 2 / q) / (gamma - 1) *
    n^(-(gamma - 1) / (2 * gamma + 1))
  c_L <- sqrt(max(0, 1 - penalty))
  b_opt <- max(1, round(n^(1 / (2 * gamma + 1))))
  m_opt <- max(1, floor(n / b_opt))
  regime <- if (c_L < 0.5) {
    "pre-asymptotic: intervals must be widened"
  } else if (K_TAC < K_MPC) {
    "TAC favored"
  } else {
    "MPC favored"
  }
  list(K_TAC = K_TAC, K_MPC = K_MPC,
       c_L = c_L, c_U = c_U,
       b_optimal = b_opt, m_optimal = m_opt,
       gamma = gamma, q = q, n = n,
       regime = regime)
}

.zeta_riemann <- function(s) {
  if (s <= 1) return(Inf)
  sum(seq_len(10000L)^(-s)) +
    stats::integrate(function(x) x^(-s), lower = 10000, upper = Inf)$value
}
