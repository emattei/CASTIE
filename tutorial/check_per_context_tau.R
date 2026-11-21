# Self-check for step 1 --perContextTau: the C++ Sigma product and Sigma^-1 solves
# (Woodbury and PCG) match a dense reference with one slope tau per context.
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

check <- function(tau, perContext) {
  CASTIE:::set_perContextTau(perContext)
  CASTIE:::set_num_Kmat(length(tau))
  slope <- if (perContext) tau[-(1:2)] else rep(tau[3], ncol(E))
  S <- diag(tau[1] / w) + tau[2] * same + (E %*% (slope * t(E))) * same
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
cat("check_per_context_tau: OK\n")
