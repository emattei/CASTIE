# Check a fitted step 1 model (count trait, eMat, repeated cells per donor, no GRM) against the REML equations,
# recomputed in plain R with exact donor-block traces: at the constrained optimum the AI-REML step over the
# components > 0 is ~0 (and boundary components have score <= 0).
#   Rscript tutorial/check_null_fit.R <outputPrefix.rda>...
# The full-system step is what an AI step that keeps boundary components would take (see --activeSetAI).
for (f in commandArgs(TRUE)) {
  e <- new.env(); load(f, envir = e); m <- e$modglmm
  mu <- as.vector(m$fitted.values); n <- length(mu); th <- m$theta
  X <- if (is.null(m$X)) m$obj.noK$XXVX_inv %*% m$obj.noK$XVX else m$X # older null models lack X
  y <- as.vector(m$y); off <- if (is.null(m$offset)) 0 else as.vector(m$offset)
  Yw <- log(mu) - off + (y - mu) / mu                   # working response, Poisson log link
  W <- mu; U <- cbind(1, m$eMat); K <- ncol(U)
  d <- match(m$sampleID, unique(m$sampleID)); D <- max(d); G <- diag(th[-1], K)
  A <- array(0, c(K, K, D))
  for (a in 1:K) for (b in a:K) { s <- rowsum(W * U[, a] * U[, b], d)[, 1]; A[a, b, ] <- s; A[b, a, ] <- s }
  H <- array(0, c(K, K, D)); tr1 <- numeric(K)
  for (j in 1:D) { H[, , j] <- solve(diag(K) + G %*% A[, , j]); tr1 <- tr1 + diag(A[, , j] %*% H[, , j]) }
  sig_inv <- function(B) {
    B <- as.matrix(B); WB <- W * B; nc <- ncol(B)
    C <- lapply(1:K, function(a) rowsum(U[, a] * WB, d)); R <- array(0, c(D, K, nc))
    for (j in 1:D) R[j, , ] <- H[, , j] %*% G %*% do.call(rbind, lapply(C, function(Ca) Ca[j, , drop = FALSE]))
    corr <- matrix(0, n, nc); for (a in 1:K) corr <- corr + U[, a] * matrix(R[d, a, ], n, nc)
    WB - W * corr
  }
  SX <- sig_inv(X); covm <- solve(crossprod(X, SX))
  PY <- drop(sig_inv(Yw) - SX %*% (covm %*% crossprod(SX, Yw)))
  tr <- sapply(1:K, function(k) { M <- rowsum(U[, k] * SX, d); tr1[k] - sum((M %*% covm) * M) })
  ypapy <- sapply(1:K, function(k) sum(rowsum(U[, k] * PY, d)^2))
  score <- ypapy - tr
  P <- function(b) drop(sig_inv(b) - SX %*% (covm %*% crossprod(SX, b)))
  KPY <- sapply(1:K, function(k) U[, k] * rowsum(U[, k] * PY, d)[d, 1])   # K_k P Y, n x K
  AI <- crossprod(KPY, sig_inv(KPY) - SX %*% (covm %*% crossprod(SX, KPY)))
  free <- th[-1] > 0                                                         # boundary components stay at 0
  full <- solve(AI, score); cat(sprintf("  full-system step (CASTIE): %s\n", paste(signif(full, 3), collapse = " ")))
  step <- numeric(K); if (any(free)) step[free] <- solve(AI[free, free, drop = FALSE], score[free])
  cat(sprintf("%s\n  theta %s\n  trace %s\n  score %s\n  AI-REML step %s  (relative to theta %s)\n",
              basename(dirname(dirname(f))), paste(signif(th[-1], 4), collapse = " "),
              paste(round(tr, 1), collapse = " "), paste(round(score, 2), collapse = " "),
              paste(signif(step, 3), collapse = " "),
              paste(ifelse(free, sprintf("%+.4f", step / th[-1]), "boundary"), collapse = " ")))
}
