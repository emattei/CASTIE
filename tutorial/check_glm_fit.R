# castie_glm_fit against stats::glm.fit: same iterations and QR, coefficients equal up to rounding.
# Rscript tutorial/check_glm_fit.R [path to R/glm_fit_chol.R]   (default: the installed CASTIE)
a <- commandArgs(TRUE)
if (length(a)) source(a[1]) else castie_glm_fit <- CASTIE::castie_glm_fit
set.seed(1)
n <- 70000
lib <- factor(sample(20, n, TRUE))
d <- data.frame(age = rnorm(n), tss = rnorm(n), lib = lib, off = log(rgamma(n, 5, 1 / 2000)))
eta <- -8 + 0.2 * d$age + 0.3 * d$tss + rnorm(20, 0, 0.3)[lib] + d$off
d$y <- rpois(n, exp(eta))
d$dup <- d$age * 2 # rank-deficient design: must take the QR path throughout
check <- function(f, fam = poisson()) {
  t0 <- system.time(a <- glm(f, data = d, offset = off, family = fam))[3]
  t1 <- system.time(b <- glm(f, data = d, offset = off, family = fam, method = castie_glm_fit))[3]
  rel <- function(u, v) max(abs(u - v) / pmax(abs(u), 1e-300), na.rm = TRUE)
  stopifnot(a$iter == b$iter, a$converged == b$converged, a$rank == b$rank, identical(is.na(coef(a)), is.na(coef(b))),
            rel(coef(a), coef(b)) < 1e-9, rel(fitted(a), fitted(b)) < 1e-9, rel(a$deviance, b$deviance) < 1e-12,
            max(abs(a$qr$qr - b$qr$qr)) < 1e-8, identical(a$qr$pivot, b$qr$pivot), rel(a$effects, b$effects) < 1e-6)
  cat(sprintf("%-28s iter %d  coef rel %.1e  dev rel %.1e  glm %.2fs  castie %.2fs\n", deparse(f), a$iter,
              rel(coef(a), coef(b)), rel(a$deviance, b$deviance), t0, t1))
}
check(y ~ age + tss + lib)
check(y ~ age + tss + lib + dup)
check(I(y > 0) ~ age + tss + lib, binomial())
cat("castie_glm_fit OK\n")
