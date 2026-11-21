# Self-check for step 1 --perContextTau (and --slopeCov): the C++ Sigma product and Sigma^-1 solves
# (Woodbury and PCG) match a dense reference with one slope tau per context (and intercept-slope covariances).
# Run after installing CASTIE (e.g. in the container): Rscript tutorial/check_per_context_tau.R
suppressMessages(library(CASTIE))
set.seed(1)
did <- rep(1:30, sample(3:8, 30, TRUE))
N <- length(did)
E <- scale(cbind(rnorm(N), rbinom(N, 1, 0.5)))
w <- runif(N, 0.2, 2)
b <- rnorm(N)
same <- outer(did, did, "==")

CASTIE:::set_g_omp_num_threads(1)
CASTIE:::set_useGRMtoFitNULL(FALSE)
CASTIE:::set_I_mat_inR(did)
CASTIE:::set_EMat(E)

check <- function(tau, perContext, cov = FALSE) {
  CASTIE:::set_perContextTau(perContext)
  CASTIE:::set_slopeCov(cov)
  CASTIE:::set_num_Kmat(length(tau))
  k <- ncol(E)
  Tm <- diag(c(tau[2], if (perContext) tau[2 + 1:k] else rep(tau[3], k)))
  if (cov) Tm[1, -1] <- Tm[-1, 1] <- tau[2 + k + 1:k]
  Z <- cbind(1, E)
  S <- diag(tau[1] / w) + (Z %*% Tm %*% t(Z)) * same
  Sb <- drop(S %*% b)
  stopifnot(max(abs(CASTIE:::getCrossprod_multiV_eMat(b, w, tau, FALSE) - Sb)) < 1e-4 * max(abs(Sb)))
  x <- solve(S, b)
  for (pcg in c(FALSE, TRUE)) {
    CASTIE:::set_usePCG(pcg)
    xc <- CASTIE:::getPCG1ofSigmaAndVector_multiV(w, tau, b, 500, 1e-7, FALSE)
    stopifnot(max(abs(xc - x)) < 1e-3 * max(abs(x)))
  }
  CASTIE:::set_usePCG(FALSE)
}
check(c(1, 0.3, 0.2), FALSE)      # shared slope tau (default)
check(c(1, 0.3, 0.5, 0.05), TRUE) # one slope tau per context
check(c(1, 0.3, 0, 0.4), TRUE)    # a context with no slope variance
check(c(1, 0.3, 0.5, 0.05, 0.2, -0.05), TRUE, TRUE)       # intercept-slope covariances of both signs
check(c(1, 0.3, 0.5, 0.01, 0.38, -0.003), TRUE, TRUE)     # near-singular T (correlation 0.98, small tau2)
check(c(1, 0.3, 0, 0.4, 0, 0.3), TRUE, TRUE)             # covariance with a zero slope variance elsewhere
cat("check_per_context_tau: OK\n")
