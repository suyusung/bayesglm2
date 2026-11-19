## Engine parity smoke: bglm vs arm::bayesglm across a config matrix,
## and flat-prior bglm vs glm (bit-parity), plus S3 methods.
.libPaths(c("~/R-libraries", .libPaths()))
suppressMessages(library(bayesglm2)); suppressMessages(library(arm))

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
d$yp <- rpois(n, exp(0.3 + 0.5 * d$x1))
d$yg <- 2 + 1.5 * d$x1 - d$x2 + rnorm(n)
wts <- runif(n, 0.5, 1.5); offs <- rnorm(n)

pconfigs <- list(
    list(),                                              # arm defaults
    list(prior.scale = 1, prior.df = 5),
    list(prior.scale = 2.5, prior.df = 1,
         prior.mean.for.intercept = 1, prior.scale.for.intercept = 10),
    list(prior.scale = 10, prior.df = 7, scaled = FALSE))

for (cfg in pconfigs) {
    tag <- paste0("cfg", match(list(cfg), pconfigs))
    for (fam in c("gaussian", "binomial", "poisson")) {
        resp <- switch(fam, gaussian = "yg", binomial = "yb", poisson = "yp")
        F0 <- as.formula(paste(resp, "~ x1 + x2"))
        m1 <- do.call(bglm, c(list(F0, get(fam), d), cfg))
        m2 <- do.call(arm::bayesglm, c(list(F0, get(fam), d), cfg))
        chk(paste0(tag, " ", fam, " coef"),     cmp(coef(m1), coef(m2)))
        chk(paste0(tag, " ", fam, " deviance"), cmp(m1$deviance, m2$deviance))
        chk(paste0(tag, " ", fam, " iter"),     identical(m1$iter, m2$iter))
        chk(paste0(tag, " ", fam, " dispersion"), cmp(m1$dispersion, m2$dispersion))
        chk(paste0(tag, " ", fam, " nulldev"),  cmp(m1$null.deviance, m2$null.deviance))
        chk(paste0(tag, " ", fam, " aic"),      cmp(m1$aic, m2$aic))
        chk(paste0(tag, " ", fam, " linpred"),  cmp(predict(m1, type = "link"),
                                                     predict(m2, type = "link")))
        chk(paste0(tag, " ", fam, " response"), cmp(predict(m1, type = "response"),
                                                     predict(m2, type = "response")))
        chk(paste0(tag, " ", fam, " se"),       cmp(unname(sqrt(diag(vcov(m1)))),
                                                     as.vector(if (is.matrix(se.coef(m2))) se.coef(m2)[, "se"] else se.coef(m2))))
    }
    # weights + offset + factor (rank-deficient) path
    F1 <- as.formula("yg ~ x1 + g")
    m1 <- do.call(bglm, c(list(F1, gaussian, d, weights = wts, offset = offs), cfg))
    m2 <- do.call(arm::bayesglm, c(list(F1, gaussian, d, weights = wts, offset = offs), cfg))
    chk(paste0(tag, " w+off+factor coef"), cmp(coef(m1), coef(m2)))
    chk(paste0(tag, " w+off+factor vcov"), cmp(vcov(m1), vcov(m2)))
}

## flat prior vs glm (these must be bit-identical)
flat <- list(prior.scale = Inf, prior.df = Inf,
             prior.scale.for.intercept = Inf, prior.df.for.intercept = Inf)
for (fam in c("gaussian", "binomial", "poisson")) {
    resp <- switch(fam, gaussian = "yg", binomial = "yb", poisson = "yp")
    F0 <- as.formula(paste(resp, "~ x1 + x2 + g"))
    fb <- do.call(bglm, c(list(F0, get(fam), d, weights = wts, offset = offs), flat))
    fg <- do.call(glm,  list(F0, get(fam), d, weights = wts, offset = offs))
    fbtol <- if (fam == "gaussian") 1e-10 else 0   # EM augmentation round-off
    chk(paste0("flat ", fam, " coef == glm"),
        isTRUE(all.equal(unname(coef(fb)), unname(coef(fg)), tol = fbtol)))
    chk(paste0("flat ", fam, " logLik == glm"),
        isTRUE(all.equal(as.numeric(logLik(fb)), as.numeric(logLik(fg)), tol = 1e-10)))
    chk(paste0("flat ", fam, " AIC == glm"),
        isTRUE(all.equal(AIC(fb), AIC(fg), tol = 1e-10)))
    chk(paste0("flat ", fam, " fitted == glm"),
        isTRUE(all.equal(unname(fitted(fb)), unname(fitted(fg)), tol = fbtol)))
    ## Working residuals: arm's engine stores them mixed-generation
    ## (z from the pre-update iteration against the post-update eta), so
    ## the new engine does NOT bug-compare residuals with arm.  Pin the
    ## mathematically correct definition instead: (y - mu) / mu'(eta)
    ## evaluated at the FINAL fitted model, NA-padded by na.action.
    mu <- fitted(fb)
    mueta <- get(fam)()$mu.eta(predict(fb, type = "link"))
    yv <- switch(fam, gaussian = d$yg, binomial = d$yb, poisson = d$yp)
    chk(paste0("flat ", fam, " residuals == (y-mu)/mu'(eta)"),
        isTRUE(all.equal(unname(residuals(fb)), unname((yv - mu) / mueta), tol = 1e-12)))
    chk(paste0("flat ", fam, " vcov == glm"),
        isTRUE(all.equal(vcov(fb), vcov(fg), tol = 1e-10)))
    nd <- head(d, 10)
    chk(paste0("flat ", fam, " predict newdata"),
        isTRUE(all.equal(predict(fb, newdata = nd, type = "link"),
                         predict(fg, newdata = nd, type = "link"), tol = 0)))
}

cat(if (fail == 0) "engine-parity: PASS\n" else sprintf("engine-parity: %d FAILURES\n", fail))
if (fail > 0) quit(status = 1)
