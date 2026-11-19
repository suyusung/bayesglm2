/* =====================================================================
 * bglm_qrls -- the linear-algebra core of the bayesglm2 engine.
 *
 * Pivoted-QR least squares: min || X b - y ||_2 with column pivoting
 * and tolerance-based rank detection.  This is the function the old
 * engine reached through lm.fit(); the engine calls it directly.  It
 * drives dqrls -- the same netlib LINPACK routine that R's own
 * Cdqrls (the .Call behind lm.fit, src/library/stats/src/lm.c)
 * drives, with the same argument pattern, pivot initialisation and
 * work-array size -- so rank decisions and coefficients agree with
 * lm.fit by construction (pinned bit-for-bit by
 * tests/qr-kernel-vs-lm-fit.R).
 *
 * Returns list(coef, effects, residuals, rank, qr, qraux, pivot):
 *   coef       length p, original column order; aliased (unpivoted)
 *              columns are exactly zero
 *   effects    Q'y, length n (lm.fit's "effects" at full rank)
 *   residuals  y - X b, length n
 *   rank       numerical rank (leading |r_kk| > tol*|r_11|)
 *   qr         n x p: R in the upper triangle, Householder below
 *   qraux      the p tau factors
 *   pivot      1-based column pivot
 *
 * W (nullable): per-row weights.  When non-NULL the solve runs on
 * diag(w) %*% X, diag(w) %*% y with the scaling fused into the
 * work-array copy -- the same double multiplications, in one pass,
 * no intermediate n x p allocation.
 * ===================================================================== */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>    /* F77_NAME(dqrls): the lm.fit backend */
#include <string.h>

SEXP bglm_qrls(SEXP X, SEXP y, SEXP Tol, SEXP W)
{
    if (!isMatrix(X) || !isReal(X) || !isReal(y) || !isReal(Tol))
        error("'X', 'y' and 'Tol' must be double");
    int n = nrows(X), p = ncols(X), ny = 1;
    if (length(y) != n) error("'y' must have length nrow(X)");
    if (n < p) error("the design must have at least as many rows as columns");
    if (n == 0 || p == 0) error("empty design");
    double tol = REAL(Tol)[0];
    if (!R_FINITE(tol) || tol <= 0.0) error("'Tol' must be a positive finite number");
    int weighted = !isNull(W);
    if (weighted && (!isReal(W) || length(W) != n))
        error("'w', if given, must be a double vector of length nrow(X)");
    double *wp = weighted ? REAL(W) : NULL;

    for (R_xlen_t i = 0; i < XLENGTH(X); i++)
        if (!R_FINITE(REAL(X)[i])) error("NA/NaN/Inf in 'X'");
    for (R_xlen_t i = 0; i < XLENGTH(y); i++)
        if (!R_FINITE(REAL(y)[i])) error("NA/NaN/Inf in 'y'");
    if (weighted)
        for (R_xlen_t i = 0; i < XLENGTH(W); i++)
            if (!R_FINITE(wp[i])) error("NA/NaN/Inf in 'w'");

    SEXP Xwork, ywork, coef, residuals, effects, qraux, pivot, rank_, ans;
    PROTECT(Xwork = allocMatrix(REALSXP, n, p));
    PROTECT(ywork = allocVector(REALSXP, n));
    PROTECT(coef = allocVector(REALSXP, p));
    PROTECT(residuals = allocVector(REALSXP, n));
    PROTECT(effects = allocVector(REALSXP, n));
    PROTECT(qraux = allocVector(REALSXP, p));
    PROTECT(pivot = allocVector(INTSXP, p));
    int *ip = INTEGER(pivot);
    PROTECT(rank_ = allocVector(INTSXP, 1));

    /* copy into the work arrays; with weights, fuse the per-row
     * scaling (X_i * w_i, y_i * w_i) into the copy -- same double
     * multiplications in the same order as X * w %*% ... in R, one
     * pass instead of an extra n x p allocation */
    double *xw = REAL(Xwork), *x = REAL(X), *yw = REAL(ywork), *yv = REAL(y),
           *res = REAL(residuals);
    if (weighted) {
        for (R_xlen_t i = 0; i < (R_xlen_t)n * p; i++)
            xw[i] = x[i] * wp[i % n];
        for (R_xlen_t i = 0; i < (R_xlen_t)n; i++) {
            yw[i] = yv[i] * wp[i];
            res[i] = yw[i];
        }
    } else {
        memcpy(xw, x, (R_xlen_t)n * p * sizeof(double));
        memcpy(yw, yv, (R_xlen_t)n * sizeof(double));
        memcpy(res, yv, (R_xlen_t)n * sizeof(double));
    }
    for (int j = 0; j < p; j++) ip[j] = j + 1;  /* dqrls reads jpvt as input */

    double rtol = tol;
    double *work = (double *) R_alloc(2 * p, sizeof(double));

    /* identical call shape to Cdqrls in src/library/stats/src/lm.c */
    F77_NAME(dqrls)(xw, &n, &p, yw, &ny, &rtol,
                    REAL(coef), res, REAL(effects),
                    INTEGER(rank_), ip, REAL(qraux), work);

    PROTECT(ans = allocVector(VECSXP, 7));
    SET_VECTOR_ELT(ans, 0, coef);
    SET_VECTOR_ELT(ans, 1, effects);
    SET_VECTOR_ELT(ans, 2, residuals);
    SET_VECTOR_ELT(ans, 3, rank_);
    SET_VECTOR_ELT(ans, 4, Xwork);
    SET_VECTOR_ELT(ans, 5, qraux);
    SET_VECTOR_ELT(ans, 6, pivot);
    UNPROTECT(9);
    return ans;
}

static const R_CallMethodDef callMethods[] = {
    {"bglm_qrls", (DL_FUNC) &bglm_qrls, 4},
    {NULL, NULL, 0}
};

void R_init_bayesglm2(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
