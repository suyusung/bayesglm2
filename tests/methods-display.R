## =====================================================================
## display() / coef() / se.coef() for bglm objects: internal
## consistency, agreement with summary/vcov, and parity with
## arm::display.bayesglm on the binomial (dispersion == 1) case where
## the SE convention coincides.
## =====================================================================
.libPaths(c("~/R-libraries", .libPaths()))
suppressMessages(library(bayesglm2))

fail <- 0
chk <- function(label, ok) {
    cat(sprintf("%-58s %s\n", label, if (isTRUE(ok)) "OK" else "MISMATCH"))
    if (!isTRUE(ok)) fail <<- fail + 1
}
cmp <- function(a, b, tol = 1e-12) isTRUE(all.equal(a, b, tol = tol))

set.seed(42); n <- 200
d <- data.frame(x1 = rnorm(n), x2 = rnorm(n),
                g = factor(sample(letters[1:4], n, TRUE)))
d$yb <- rbinom(n, 1, plogis(0.7 * d$x1 - 0.4 * d$x2))
d$yg <- 2 + 1.5 * d$x1 - d$x2 + rnorm(n)

fit.b <- bglm(yb ~ x1 + x2 + g, data = d, family = binomial)
fit.g <- bglm(yg ~ x1 + x2 + g, data = d)

## ---- se.coef: sqrt(diag(vcov(.))) and summary column 2 -------------
chk("se.coef == sqrt(diag(vcov)) [binomial]",
    cmp(se.coef(fit.b), sqrt(diag(vcov(fit.b)))))
chk("se.coef == sqrt(diag(vcov)) [gaussian]",
    cmp(se.coef(fit.g), sqrt(diag(vcov(fit.g)))))
chk("se.coef == summary SE column",
    cmp(as.numeric(se.coef(fit.b)), unname(summary(fit.b)$coefficients[, 2])))
chk("se.coef has coefficient names",
    identical(names(se.coef(fit.b)), names(coef(fit.b))))

## ---- coef: the fit's own coefficients, plain numeric vector --------
chk("coef is named numeric of right length",
    is.numeric(coef(fit.b)) && length(coef(fit.b)) == ncol(model.matrix(fit.b)))

## ---- display: returned list ---------------------------------------
out.b <- suppressWarnings(capture.output(dsp <- display(fit.b)))
chk("display returns coef/se/n invisibly",
    identical(dsp$coef, coef(fit.b)) &&
        cmp(as.numeric(dsp$se), as.numeric(se.coef(fit.b))) &&
        dsp$n == n && dsp$k == fit.b$rank)
chk("display deviances match the fit",
    cmp(dsp$deviance, fit.b$deviance) &&
        cmp(dsp$null.deviance, fit.b$null.deviance))
chk("binomial footer has no overdispersion line",
    !any(grepl("overdispersion", out.b)))
chk("display prints the call and the table",
    any(grepl("^bglm\\(", out.b)) && any(grepl("coef.est", out.b)) &&
        any(grepl("---", out.b)))

out.g <- suppressWarnings(capture.output(dsp.g <- display(fit.g)))
chk("gaussian footer shows overdispersion + sigma.hat",
    any(grepl("overdispersion parameter =", out.g)) &&
        any(grepl("residual sd is sqrt\\(overdispersion\\) =", out.g)))
chk("sigma.hat == sqrt(deviance/df)",
    cmp(dsp.g$sigma.hat, sqrt(fit.g$deviance / fit.g$df.residual)))

## detail = TRUE adds z/t and p
invisible(suppressWarnings(capture.output(dsp.d <- display(fit.b, detail = TRUE))))
chk("detail z.value == coef/se",
    cmp(as.numeric(dsp.d$z.value), as.numeric(coef(fit.b) / se.coef(fit.b))))
chk("detail p.value == 2*pnorm(|z|, lower.tail=FALSE)",
    cmp(as.numeric(dsp.d$p.value),
        unname(2 * pnorm(abs(coef(fit.b) / se.coef(fit.b)),
                         lower.tail = FALSE))))

## ---- fround: no "-0" ------------------------------------------------
chk("fround snaps tiny negatives to clean zero",
    identical(format(fround(-0.001, 2)), " 0") ||
        !grepl("-", format(fround(-0.001, 2))))

## ---- parity with arm's display.bayesglm (binomial: dispersion == 1)
if (requireNamespace("arm", quietly = TRUE)) {
    suppressMessages(library(arm))
    a.b <- arm::bayesglm(yb ~ x1 + x2 + g, data = d, family = binomial)
    invisible(suppressWarnings(capture.output(ad <- arm::display(a.b))))
    chk("coef matches arm::bayesglm", cmp(as.numeric(coef(fit.b)), unname(coef(a.b))))
    chk("display se matches arm display.bayesglm se",
        cmp(as.numeric(dsp$se), unname(ad$se), tol = 1e-10))
    chk("display n/k match arm", dsp$n == ad$n && dsp$k == ad$k)
    chk("display deviances match arm",
        cmp(dsp$deviance, ad$deviance) && cmp(dsp$null.deviance, ad$null.deviance))
    ## arm's S4 se.coef generic masks ours once arm is attached; call
    ## ours by namespace.  arm's se.coef on bayesglm dispatches to
    ## se.coef.glm
    chk("se.coef matches arm::se.coef on the bayesglm fit",
        cmp(as.numeric(bayesglm2::se.coef(fit.b)),
            unname(as.numeric(arm::se.coef(a.b))), tol = 1e-10))
} else {
    message("arm not installed -- display parity checks skipped")
}

if (fail > 0) quit(status = fail)
cat("methods-display: all checks passed\n")
