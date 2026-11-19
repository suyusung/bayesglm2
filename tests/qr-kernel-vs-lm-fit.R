## Bit-for-bit pin: our pivoted-QR kernel must agree with lm.fit's
## coefficients and rank exactly (same dqrls, same tol = 1e-7).
stopifnot(require(bayesglm2))

mk <- function(seed, n, p, rank = NULL) {
    set.seed(seed)
    X <- matrix(rnorm(n * p), n, p)
    if (!is.null(rank)) {
        Q <- qr.Q(qr(matrix(rnorm(n * rank), n, rank)))
        X <- Q %*% matrix(rnorm(rank * p), rank, p)   # exactly rank `rank`
    }
    colnames(X) <- paste0("x", seq_len(ncol(X)))
    y <- drop(X %*% (seq_len(ncol(X)) * 0.5)) + rnorm(n)
    list(X = X, y = as.vector(y))
}

cases <- list(mk(1, 50, 3), mk(2, 60, 5), mk(3, 80, 4, rank = 2),
              mk(4, 40, 6, rank = 3), mk(5, 30, 1))

for (cs in cases) {
    mine <- bayesglm2:::bglm_qr_solve(cs$X, cs$y)
    ref  <- stats::lm.fit(cs$X, cs$y)
    stopifnot(identical(mine$rank, ref$rank))
    stopifnot(isTRUE(all.equal(unname(mine$coefficients),
                               unname(ref$coefficients), tol = 0)))
    stopifnot(identical(dimnames(mine$qr$qr), dimnames(ref$qr$qr)))
}

# integer input must be coerced, not rejected
Xi <- matrix(1:8, 4, 2); yi <- c(1, 3, 5, 7)
mine <- bayesglm2:::bglm_qr_solve(Xi, yi)
stopifnot(isTRUE(all.equal(unname(mine$coefficients),
                           unname(stats::lm.fit(Xi, yi)$coefficients),
                           tol = 0)))

cat("qr-kernel-vs-lm-fit: PASS\n")
