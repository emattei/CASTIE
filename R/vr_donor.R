# Donor-level variance-ratio quantities for extractVarianceRatio_multiV when exactTrace = TRUE.
#
# Sigma is block-diagonal by donor: Sigma_d = diag(tau0 / W_d) + U_d Gm U_d' with U = [1, eMat] and
# Gm = diag(tau1, slope taus). A marker's genotype is constant within a donor (G0 = g[b]), so every
# per-marker quantity of the variance-ratio loop is a linear or quadratic form in the donor genotype g:
#   1_d' Sigma_d^-1 1_d etc. = Q_d = A_d - A_d (I + Gm A_d)^-1 Gm A_d,  A_d = U_d' diag(W_d / tau0) U_d
#   X' Sigma^-1 (g[b] * e)    = rowsum(Sigma^-1 X * e, b)' g
# and the X-adjustments of the no-random-effect model expand through donor sums of the noK matrices. The
# N-length work is done once per model here; each marker then costs O(D p) instead of O(N p).
# Sigma^-1 (G - X c) = Sigma^-1 G - Sigma^-1 X c, and P = Sigma^-1 - Sigma^-1 X A^-1 X' Sigma^-1 has P X = 0,
# so var1 (= G~' P G~) needs only the unadjusted donor genotype.
.vr_donor_setup <- function(obj, W, tauVecNew, b, Sigma_iX, XtSigma_iX_inv) {
  eMat <- obj$eMat
  nE <- ncol(eMat)
  k <- 1 + nE
  D <- max(b)
  U <- cbind(1, eMat)
  wt <- as.numeric(W) / tauVecNew[1] # model vectors can be stored as N x 1 matrices
  Gm <- c(tauVecNew[2], if (length(tauVecNew) == 2 + nE) tauVecNew[2 + seq_len(nE)] else rep(tauVecNew[3], nE))
  dsum <- function(x) rowsum(x, b, reorder = TRUE) # donor sums; b is 1..D
  A <- array(0, c(D, k, k))
  for (a in 1:k) for (cc in a:k) A[, a, cc] <- A[, cc, a] <- dsum(wt * U[, a] * U[, cc])[, 1]
  Q <- array(0, c(D, k, k))
  for (d in 1:D) {
    Ad <- A[d, , ]
    Q[d, , ] <- Ad - Ad %*% solve(diag(k) + Gm * Ad, Gm * Ad)
  }
  noK <- obj$obj.noK
  B <- noK$XXVX_inv # G~ = G0 - B (XV G0)
  tXV <- t(noK$XV)
  mu <- as.numeric(obj$fitted.values)
  m <- (if (obj$traitType == "binary") mu * (1 - mu) else mu) * as.numeric(obj$varWeights) # var2null weights
  r <- as.numeric(obj$residuals) * as.numeric(obj$varWeights)
  list(
    b = b, nE = nE, D = D, p = ncol(Sigma_iX), Ainv = XtSigma_iX_inv, eMat = eMat, B = B,
    q00 = Q[, 1, 1], qjj = Q[, -1, -1, drop = FALSE], qj0 = Q[, -1, 1, drop = FALSE],
    RX = dsum(Sigma_iX), RXE = lapply(seq_len(nE), function(j) dsum(Sigma_iX * eMat[, j])),
    RV = dsum(tXV), RVE = lapply(seq_len(nE), function(j) dsum(tXV * eMat[, j])),
    Rr = dsum(r)[, 1], Rer = dsum(eMat * r), Br = colSums(B * r), BrD = dsum(B * r),
    Mm = dsum(m)[, 1], MmE2 = dsum(eMat^2 * m), BmD = dsum(B * m),
    BmDE = lapply(seq_len(nE), function(j) dsum(B * (m * eMat[, j]))), BmB = crossprod(B, B * m)
  )
}

# The per-marker quantities of the variance-ratio loop for donor genotype g (cell genotype g[b]).
# var2sparseGE needs Sigma_V^-1 and stays with the caller.
.vr_donor_marker <- function(s, g) {
  g <- as.numeric(g) # Get_OneSNP_Geno* can return a D x 1 matrix
  g2 <- g * g
  u <- crossprod(s$RX, g)[, 1]
  var1 <- sum(g2 * s$q00) - sum(u * (s$Ainv %*% u))
  c0 <- crossprod(s$RV, g)[, 1] # XV G0
  S <- sum(g * s$Rr) - sum(c0 * s$Br)
  var2null <- sum(g2 * s$Mm) - 2 * sum(g * (s$BmD %*% c0)) + sum(c0 * (s$BmB %*% c0))
  S_donor <- g * s$Rr - (s$BrD %*% c0)[, 1]
  var1GE <- I_21 <- var2nullGE <- numeric(s$nE)
  R_donor <- matrix(0, s$D, s$nE)
  cE <- matrix(0, ncol(s$RV), s$nE)
  for (j in seq_len(s$nE)) {
    cj <- crossprod(s$RVE[[j]], g)[, 1] # XV (G0 * e_j)
    uj <- crossprod(s$RXE[[j]], g)[, 1]
    Ainv_uj <- s$Ainv %*% uj
    var1GE[j] <- sum(g2 * s$qjj[, j, j]) - sum(uj * Ainv_uj)
    I_21[j] <- sum(g2 * s$qj0[, j, 1]) - sum(u * Ainv_uj)
    var2nullGE[j] <- sum(g2 * s$MmE2[, j]) - 2 * sum(g * (s$BmDE[[j]] %*% cj)) + sum(cj * (s$BmB %*% cj))
    R_donor[, j] <- g * s$Rer[, j] - (s$BrD %*% cj)[, 1]
    cE[, j] <- cj
  }
  c_coeff <- I_21 / var1
  Q_donor <- R_donor - outer(S_donor, c_coeff)
  centered <- sweep(Q_donor, 2, colSums(Q_donor) / s$D, `-`)
  list(
    var1 = var1, S = S, var2null = var2null, var1GE = var1GE, var2nullGE = var2nullGE,
    varModelGEcond = var1GE - I_21^2 / var1, varSWGEcond = colSums(centered^2) * s$D / (s$D - s$p), cE = cE
  )
}
