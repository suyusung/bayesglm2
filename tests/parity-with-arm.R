## =====================================================================
## Parity test: bayesglm2::bglm_design must agree with
## arm:::model.matrixBayes to the bit on a battery of shapes, until the
## engine is ported and arm's copy is replaced by this one.
## =====================================================================

if (!requireNamespace("arm", quietly = TRUE)) {
    message("arm not installed -- parity check skipped")
    quit(status = 0)
}

library(bayesglm2)

set.seed(42)
n <- 30
d <- data.frame(
    ynum  = rnorm(n),
    x     = rnorm(n),
    z     = rnorm(n),
    g     = factor(sample(c("a", "b", "c"), n, TRUE), levels = c("a", "b", "c")),
    h     = factor(sample(c("p", "q"), n, TRUE)),
    o     = factor(sample(c("lo", "md", "hi"), n, TRUE), ordered = TRUE),
    chr   = sample(c("u", "v"), n, TRUE),
    lg    = sample(c(TRUE, FALSE), n, TRUE),
    psi   = rep(0:5, length.out = n))
d$ybin <- rbinom(n, 1, plogis(d$x))

cases <- list(
    list("numeric",            ynum ~ x),
    list("numeric-no-int",     ynum ~ x - 1),
    list("two-numeric",        ynum ~ x + z),
    list("factor",             ynum ~ g),
    list("factor-ordered",     ynum ~ o),
    list("factor+numeric",     ynum ~ g + x),
    list("interaction",        ynum ~ g:h),
    list("num-interaction",    ynum ~ x:z + x),
    list("logical",            ynum ~ lg),
    list("offset",             ynum ~ x + offset(z)),
    list("transformed-poly",   ynum ~ 1 + poly(psi, degree = 3, raw = TRUE)),
    list("transformed-log",    ynum ~ 1 + log(psi + 1)),
    list("transformed-I",      ynum ~ 1 + I(psi^2) + psi),
    list("binomial",           ybin ~ g + x),
    list("multi",              ynum ~ g + h + o + x:z),
    list("null",               ynum ~ 1))

## our package stores the primary generator names in the "contrasts"
## attribute; arm stores the legacy strings -- compare as equivalent
canon <- function(X) {
    nm <- function(s) if (is.character(s) && length(s) == 1L)
        switch(s,
            "contr.bayes.unordered" = "contr.treatment.full",
            "contr.bayes.ordered"   = "contr.poly.full",
            s)
        else s
    if (!is.null(a <- attr(X, "contrasts"))) attr(X, "contrasts") <- lapply(a, nm)
    X
}


check <- function(label, f, drop.baseline, keep.order, ...) {
    mine <- bayesglm2::bglm_design(f, data = d,
        drop.baseline = drop.baseline, keep.order = keep.order, ...)
    theirs <- arm:::model.matrixBayes(f, data = d,
        drop.baseline = drop.baseline, keep.order = keep.order, ...)
    if (!identical(canon(mine), canon(theirs))) {
        cat("MISMATCH:", label, "drop.baseline =", drop.baseline,
            "keep.order =", keep.order, "\n")
        print(all.equal(theirs, mine))
        quit(status = 1)
    }
}

for (cs in cases) for (db in c(TRUE, FALSE)) for (ko in c(TRUE, FALSE)) {
    check(cs[[1]], cs[[2]], drop.baseline = db, keep.order = ko)
}

## terms-path battery: with the intercept present the session-default
## contrasts actually differ from the bayes full-rank ones (on the
## formula path the absent intercept makes even contr.treatment return
## full dummies, so drop.baseline is unobservable there).
for (cs in cases) for (db in c(TRUE, FALSE)) for (ko in c(TRUE, FALSE)) {
    tt <- stats::terms(cs[[2]], data = d)
    mff <- stats::model.frame(tt, d)
    mine <- bayesglm2::bglm_design(tt, data = mff,
        drop.baseline = db, keep.order = ko)
    theirs <- arm:::model.matrixBayes(tt, data = mff,
        drop.baseline = db, keep.order = ko)
    if (!identical(canon(mine), canon(theirs))) {
        cat("MISMATCH (terms path):", cs[[1]], "drop.baseline =", db,
            "keep.order =", ko, "\n")
        print(all.equal(theirs, mine))
        quit(status = 1)
    }
}

## contrasts.arg battery
check("contr.sum", ynum ~ g, TRUE, TRUE, contrasts.arg = list(g = "contr.sum"))
check("contr.mat", ynum ~ g, FALSE, TRUE,
      contrasts.arg = list(g = matrix(c(1, 0, 0, 0, 1, 0), 3, 2)))
check("contr.absent-warn", ynum ~ g, TRUE, FALSE,
      contrasts.arg = list(nope = "contr.sum"))

## model-frame input (the transformed-terms path predict() uses)
for (db in c(TRUE, FALSE)) {
    f <- ynum ~ 1 + poly(psi, degree = 3, raw = TRUE)
    mf <- stats::model.frame(f, d)
    mine <- bayesglm2::bglm_design(stats::terms(f), data = mf,
                                         drop.baseline = db, keep.order = TRUE)
    theirs <- arm:::model.matrixBayes(stats::terms(f), data = mf,
                                      drop.baseline = db, keep.order = TRUE)
    stopifnot("frame-input parity failed" =
                identical(canon(mine), canon(theirs)))
}

## character column
check("character", ynum ~ chr, TRUE, TRUE)
check("character-full", ynum ~ chr, FALSE, TRUE)

cat("parity-with-arm: PASS\n")
