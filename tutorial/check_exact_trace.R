# Self-check for step 1 --exactTrace: the exact donor-block traces tr(P K_k) of GetTrace_multiV_eMat match a
# dense reference, with a shared slope tau and with one slope tau per context (some taus 0).
# Run after installing CASTIE (e.g. in the container): Rscript tutorial/check_exact_trace.R
suppressMessages(library(CASTIE))
set.seed(1)
did <- rep(1:30, sample(3:8, 30, TRUE))
N <- length(did)
E <- scale(cbind(rnorm(N), rbinom(N, 1, 0.5)))
X <- cbind(1, rnorm(N), E)
w <- runif(N, 0.2, 2)
same <- outer(did, did, "==")

CASTIE:::set_g_omp_num_threads(1)
CASTIE:::set_useGRMtoFitNULL(FALSE)
CASTIE:::set_I_mat_inR(did)
CASTIE:::set_EMat(E)
CASTIE:::set_usePCG(FALSE)

check <- function(tau, perContext) {
  CASTIE:::set_perContextTau(perContext)
  CASTIE:::set_num_Kmat(length(tau))
  slope <- if (perContext) tau[-(1:2)] else rep(tau[3], ncol(E))
  S <- diag(tau[1] / w) + tau[2] * same + (E %*% (slope * t(E))) * same
  SX <- solve(S, X)
  covm <- solve(crossprod(X, SX))
  P <- solve(S) - SX %*% covm %*% t(SX)
  kern <- c(list(same), lapply(seq_len(ncol(E)), function(j) outer(E[, j], E[, j]) * same))
  if (!perContext) kern <- list(kern[[1]], Reduce(`+`, kern[-1]))
  ref <- sapply(kern, function(K) sum(P * K)) # tr(P K) for symmetric K
  CASTIE:::set_exactTrace(TRUE)
  ex <- CASTIE:::GetTrace_multiV_eMat(SX, X, w, tau, c(1L, rep(0L, length(tau) - 1)), covm, 30, 500, 1e-7, 0.0025, FALSE)
  CASTIE:::set_exactTrace(FALSE)
  stopifnot(max(abs(ex - ref) / abs(ref)) < 1e-4)
}
check(c(1, 0.3, 0.2), FALSE)      # shared slope tau (default)
check(c(1, 0.3, 0.5, 0.05), TRUE) # one slope tau per context
check(c(1, 0.3, 0, 0.4), TRUE)    # a context with no slope variance
check(c(1, 0, 0, 0), TRUE)        # all variances 0
cat("check_exact_trace: OK\n")
