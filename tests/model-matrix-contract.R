## =====================================================================
## Contract test: model.matrixBayes behaviour pinned by arm's own tests
## and by reverse dependencies (FindIt formula path, qgcomp transformed
## terms, keep.order coefficients ordering).
## =====================================================================

library(bayesglm2)

close <- function(a, b) isTRUE(all.equal(a, b, tolerance = 1e-8))

set.seed(1)
d <- data.frame(y = rnorm(20), x = rnorm(20))

## ---- the formula path: NO intercept column (FindIt contract) --------
X <- bglm_design(y ~ x, data = d)
stopifnot("formula path grew an intercept column" =
            ncol(X) == 1L && identical(colnames(X), "x"))

build <- function(data) {
    f <- as.formula(paste("y ~", "x"))
    bglm_design(f, data = data)
}
Xc <- build(d)
stopifnot("formula-assembled-at-runtime path differs" =
            ncol(Xc) == 1L && identical(colnames(Xc), "x"))

## ---- the terms path: intercept kept ----------------------------------
t <- stats::terms(y ~ x, data = d)
Xt <- bglm_design(t, data = stats::model.frame(y ~ x, data = d))
stopifnot("terms path lost the intercept" =
            ncol(Xt) == 2L && identical(colnames(Xt), c("(Intercept)", "x")))

## ---- no intercept in the formula means no intercept anywhere ---------
Xi <- bglm_design(y ~ x - 1, data = d)
stopifnot("an intercept appeared for a formula without one" = ncol(Xi) == 1L)

## ---- keep.order: columns in written order ----------------------------
## The written order reaches the matrix through the TERMS object (the
## bayesglm() path): model.matrix.default takes its terms from the
## model frame when the frame has a "terms" attribute, so the order is
## whatever the terms object that BUILT the frame had.  A raw formula +
## raw data frame therefore sorts like stats::terms() -- in arm as well
## (pinned by the parity test).
set.seed(8)
n <- 40
dd <- data.frame(y = rnorm(n), x = rnorm(n), z = rnorm(n))
tt <- stats::terms.formula(y ~ x:z + x, data = dd, keep.order = TRUE)
k <- bglm_design(tt, data = stats::model.frame(tt, dd), keep.order = TRUE)
stopifnot("keep.order = TRUE does not follow the written term order" =
            identical(colnames(k), c("(Intercept)", "x:z", "x")))
tt0 <- stats::terms.formula(y ~ x:z + x, data = dd, keep.order = FALSE)
k0 <- bglm_design(tt0, data = stats::model.frame(tt0, dd),
                        keep.order = FALSE)
stopifnot("keep.order = FALSE does not sort alphabetically" =
            identical(colnames(k0), c("(Intercept)", "x", "x:z")))

## ---- drop.baseline = FALSE: full dummies, one column per level -------
## Tested on the TERMS path (intercept present).  On the formula path
## the missing intercept makes even contr.treatment return full dummies
## (model.matrix.default passes contrasts = has.intercept), arm
## included -- pinned by the parity test.
df <- data.frame(y = rnorm(12),
                 g = factor(rep(c("a", "b", "c"), each = 4)),
                 o = factor(rep(c("lo", "md", "hi"), each = 4),
                            ordered = TRUE))
tdf <- stats::terms(y ~ g + o, data = df)
mdf <- stats::model.frame(tdf, df)
Fb <- bglm_design(tdf, data = mdf, drop.baseline = FALSE)
stopifnot("unordered factor did not get one column per level" =
            sum(grepl("^g", colnames(Fb))) == 3L &&
            sum(grepl("^o", colnames(Fb))) == 3L)
Tb <- bglm_design(tdf, data = mdf, drop.baseline = TRUE)
stopifnot("drop.baseline = TRUE should use the session default (n - 1 cols)" =
            sum(grepl("^g", colnames(Tb))) == 2L &&
            sum(grepl("^o", colnames(Tb))) == 2L)
## full dummies reproduce the response exactly when it is a function of
## the cells (every coefficient free: one per level)
y_cells <- c(2, -1, 0.5)[as.integer(df$g)]
stopifnot("full dummy coding does not reproduce a cell-constant response" =
            max(abs(qr.resid(qr(Fb), y_cells))) < 1e-8)

## ---- contrasts.arg: function name and matrix (terms path) ------------
Ca <- bglm_design(tdf, data = mdf, contrasts.arg = list(g = "contr.sum"))
stopifnot("contrasts.arg function name ignored" =
            ncol(Ca) == 6L && sum(Ca[, 2]) == 0 && sum(Ca[, 3]) == 0)
M <- matrix(c(1, 0, 0, 0, 1, 0), 3, 2)
Cm <- bglm_design(tdf, data = mdf, contrasts.arg = list(g = M))
stopifnot("matrix contrasts.arg ignored" =
            identical(unname(Cm[c(1, 5, 9), 2:3, drop = FALSE]), unname(M)))
w <- NULL
withCallingHandlers(
    bglm_design(y ~ g, data = df, contrasts.arg = list(nope = "contr.sum")),
    warning = function(x) { w <<- x$message; invokeRestart("muffleWarning") })
stopifnot("absent contrasts.arg variable did not warn" = grepl("absent", w))

## ---- character columns become factors with a warning -----------------
dc <- data.frame(y = rnorm(12), g = rep(c("a", "b", "c"), each = 4))
w <- NULL
withCallingHandlers(
    Xcc <- bglm_design(y ~ g, data = dc, drop.baseline = TRUE),
    warning = function(x) { w <<- x$message; invokeRestart("muffleWarning") })
stopifnot("character column not converted with the documented warning" =
            grepl("converted to a factor", w) && ncol(Xcc) == 3L)

## ---- xlev: new factor levels must error ------------------------------
lv <- factor(c("a", "b"), levels = c("a", "b"))
stopifnot("xlev did not reject a new level" =
            inherits(try(model.matrixBayes(y ~ g, data = dc,
                                           xlev = list(g = lv)),
                         silent = TRUE), "try-error"))

## ---- transformed terms via a model frame (qgcomp shape) --------------
set.seed(42)
dq <- data.frame(Ya = rnorm(60), psi = rep(0:5, 10))
fq <- as.formula("Ya ~ 1 + poly(psi, degree = 3, raw = TRUE)")
mq <- stats::model.frame(fq, dq)
Xq <- bglm_design(terms(fq), data = mq, keep.order = TRUE)
stopifnot("transformed-term frame path lost columns" =
            ncol(Xq) == 4L && grepl("poly\\(psi", colnames(Xq)[2]))
Xq2 <- bglm_design(fq, data = dq)          # formula path, no intercept
stopifnot("transformed-term formula path lost columns" = ncol(Xq2) == 3L)

## ---- null model: ~1 needs rows, not variables ------------------------
X0 <- bglm_design(y ~ 1, data = d)
stopifnot("null-model formula path wrong" = ncol(X0) == 0L && nrow(X0) == 20L)
t0 <- stats::terms(y ~ 1, data = d)
X0t <- bglm_design(t0, data = stats::model.frame(y ~ 1, data = d))
stopifnot("null-model terms path wrong" =
            ncol(X0t) == 1L && identical(colnames(X0t), "(Intercept)") &&
            nrow(X0t) == 20L)

cat("model-matrix-contract: PASS\n")
## ---- pinned legacy names are the same objects ------------------------
stopifnot("model.matrixBayes is no longer identical to bglm_design" =
            identical(bayesglm2::model.matrixBayes, bayesglm2::bglm_design))
stopifnot("contr.bayes.unordered is no longer identical to contr.treatment.full" =
            identical(bayesglm2::contr.bayes.unordered, bayesglm2::contr.treatment.full))
stopifnot("contr.bayes.ordered is no longer identical to contr.poly.full" =
            identical(bayesglm2::contr.bayes.ordered, bayesglm2::contr.poly.full))

cat("model-matrix-contract: PASS\n")
