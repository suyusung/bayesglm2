## ===================================================================
## display / coef / se.coef for bglm objects.
##
## Standalone S3 generics: bayesglm2 must not depend on arm's
## (S4) display/se.coef generics.  The output format mirrors
## arm 1.15.4's display.bayesglm: call, coefficient table rounded to
## `digits` (coef.est / coef.se, plus z-or-t value and p.value with
## detail = TRUE), then the n/k + deviance footer and, when the
## dispersion is not 1, the overdispersion parameter (and residual sd
## for gaussian).
##
## se.coef.bglm mirrors arm's se.coef.glm (which is what se.coef() on
## an arm::bayesglm dispatches to, via its "glm" class): the SEs come
## from vcov(), i.e. the package convention dispersion (vcov.glm's
## deviance/df rule, see bglm.dispersion.for.se), NOT the EM
## dispersion the engine stores internally.
## ===================================================================

display <- function(object, ...) UseMethod("display")

se.coef <- function(object, ...) UseMethod("se.coef")

se.coef.bglm <- function(object, ...) sqrt(diag(vcov(object)))

## arm's fround: round, then snap values that rounded away to nothing
## to a clean 0 (kills "-0.00" and "-0")
fround <- function(x, digits = 2) {
    x <- round(x, digits)
    x[abs(x) < 10^(-digits)] <- 0
    x + 0
}

## column-rounded matrix print (arm's pfround)
pfround <- function(x, digits = 2) {
    if (!is.matrix(x)) x <- as.matrix(x)
    print(format(fround(x, digits)), quote = FALSE)
    invisible(x)
}

display.bglm <- function(object, digits = 2, detail = FALSE, ...) {
    out <- NULL
    out$call <- object$call
    summ <- summary(object, ...)
    coefs <- summ$coefficients
    if (isTRUE(detail)) {
        coef <- coefs
        out$z.value <- coef[, 3]
        out$p.value <- coef[, 4]
    }
    else {
        coef <- coefs[, 1:2, drop = FALSE]
    }
    dimnames(coef)[[2]][1:2] <- c("coef.est", "coef.se")
    out$coef <- coef[, "coef.est"]
    out$se <- coef[, "coef.se"]
    out$n <- object$df.residual + object$rank
    out$k <- object$rank
    out$deviance <- summ$deviance
    out$null.deviance <- summ$null.deviance
    print(out$call)
    pfround(coef, digits)
    cat("---\n")
    cat(paste("n = ", out$n, ", k = ", out$k,
        "\nresidual deviance = ", fround(out$deviance, 1),
        ", null deviance = ", fround(out$null.deviance, 1),
        " (difference = ", fround(out$null.deviance - out$deviance, 1), ")",
        "\n", sep = ""))
    out$dispersion <- summ$dispersion
    if (out$dispersion != 1) {
        out$overdispersion.parameter <- out$dispersion
        cat(paste("overdispersion parameter = ",
            fround(out$dispersion, 1), "\n", sep = ""))
        if (object$family$family == "gaussian") {
            out$sigma.hat <- sqrt(out$dispersion)
            cat(paste("residual sd is sqrt(overdispersion) = ",
                fround(out$sigma.hat, digits), "\n", sep = ""))
        }
    }
    invisible(out)
}
