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
 * ===================================================================== */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>    /* F77_NAME(dqrls): the lm.fit backend */
#include <string.h>

SEXP bglm_qrls(SEXP X, SEXP y, SEXP Tol)
{
    if (!isMatrix(X) || !isReal(X) || !isReal(y) || !isReal(Tol))
        error("'X', 'y' and 'Tol' must be double");
    int n = nrows(X), p = ncols(X), ny = 1;
    if (length(y) != n) error("'y' must have length nrow(X)");
    if (n < p) error("the design must have at least as many rows as columns");
    if (n == 0 || p == 0) error("empty design");
    double tol = REAL(Tol)[0];
    if (!R_FINITE(tol) || tol <= 0.0) error("'Tol' must be a positive finite number");

    for (R_xlen_t i = 0; i < XLENGTH(X); i++)
        if (!R_FINITE(REAL(X)[i])) error("NA/NaN/Inf in 'X'");
    for (R_xlen_t i = 0; i < XLENGTH(y); i++)
        if (!R_FINITE(REAL(y)[i])) error("NA/NaN/Inf in 'y'");

    SEXP Xwork, ywork, coef, residuals, effects, qraux, pivot, rank_, ans;
    PROTECT(Xwork = allocMatrix(REALSXP, n, p));
    memcpy(REAL(Xwork), REAL(X), (R_xlen_t)n * p * sizeof(double));
    PROTECT(ywork = allocVector(REALSXP, n));
    memcpy(REAL(ywork), REAL(y), (R_xlen_t)n * sizeof(double));
    PROTECT(coef = allocVector(REALSXP, p));
    PROTECT(residuals = allocVector(REALSXP, n));
    memcpy(REAL(residuals), REAL(y), (R_xlen_t)n * sizeof(double));
    PROTECT(effects = allocVector(REALSXP, n));
    PROTECT(qraux = allocVector(REALSXP, p));
    PROTECT(pivot = allocVector(INTSXP, p));
    int *ip = INTEGER(pivot);
    for (int j = 0; j < p; j++) ip[j] = j + 1;  /* dqrls reads jpvt as input */
    PROTECT(rank_ = allocVector(INTSXP, 1));

    double rtol = tol;
    double *work = (double *) R_alloc(2 * p, sizeof(double));

    /* identical call shape to Cdqrls in src/library/stats/src/lm.c */
    F77_NAME(dqrls)(REAL(Xwork), &n, &p, REAL(ywork), &ny, &rtol,
                    REAL(coef), REAL(residuals), REAL(effects),
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
    {"bglm_qrls", (DL_FUNC) &bglm_qrls, 3},
    {NULL, NULL, 0}
};

void R_init_bayesglm2(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
