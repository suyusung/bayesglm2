## ===================================================================
## S3 methods for the bglm fitted object.
##
## The fit is a plain list (class "bglm") whose components carry every
## piece the methods need; nothing here inherits from or delegates to
## glm or lm objects.  fitted()/residuals() need no methods: the stats
## defaults already apply na.action padding through $na.action.
## ===================================================================

coef.bglm <- function(object, ...) object$coefficients

## summary.glm's dispersion rule, shared by summary/vcov so both agree
## with what vcov.glm does on an arm::bayesglm object.  NOTE this is
## deliberately NOT object$dispersion (the EM estimate the engine
## stores, which logLik and the engine's own convergence use): for the
## gaussian family summary.glm recomputes deviance / df.residual.
bglm.dispersion.for.se <- function(object) {
    if (object$family$family %in% c("binomial", "poisson")) 1
    else object$deviance / object$df.residual
}

vcov.bglm <- function(object, ...) {
    p <- object$rank
    nvars <- NCOL(object$R)
    V <- matrix(NA_real_, nvars, nvars,
        dimnames = list(colnames(object$R), colnames(object$R)))
    if (p > 0L) {
        Vp <- chol2inv(object$R[seq_len(p), seq_len(p), drop = FALSE]) *
            bglm.dispersion.for.se(object)
        piv <- object$qr$pivot[seq_len(p)]
        V[piv, piv] <- Vp
    }
    V
}

## the saturated log-likelihood pieces, family by family (glm's own
## decomposition: ll = -dev/2 + saturated); parity against
## stats::logLik.glm is pinned for flat priors in the test suite.
logLik.bglm <- function(object, ...) {
    ## exact mirror of stats:::logLik.glm: p - aic/2 with nobs/df attrs.
    fam <- object$family$family
    p <- object$rank
    disp <- object$family$dispersion
    p <- if (!is.null(disp)) {
        if (is.na(disp)) p + 1L else p
    } else if (fam %in% c("gaussian", "Gamma", "inverse.gaussian")) {
        p + 1L
    }
    val <- p - object$aic / 2
    attr(val, "nobs") <- sum(!is.na(object$residuals))
    attr(val, "df") <- p
    class(val) <- "logLik"
    val
}

AIC.bglm <- function(object, ...) object$aic

nobs.bglm <- function(object, ...) NROW(object$y)

df.residual.bglm <- function(object, ...) object$df.residual

model.matrix.bglm <- function(object, ...) {
    if (!is.null(object[["x"]])) {   # [[ : $ would partial-match xlevels
        return(object[["x"]])
    }
    bglm_design(object$terms, data = model.frame(object),
        contrasts.arg = object$contrasts,
        keep.order = object$keep.order, drop.baseline = object$drop.baseline)
}

model.frame.bglm <- function(formula, ...) stats::model.frame.default(formula, ...)

print.bglm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
    cat("Call:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    print(format(coef(x), digits = digits), quote = FALSE)
    cat("\nDegrees of freedom:", x$df.residual, "residual,", x$df.null,
        "null\nDeviance:", format(x$deviance, digits = digits),
        "on", x$iter, "iterations\n")
    invisible(x)
}

summary.bglm <- function(object, dispersion = NULL, ...) {
    p <- object$rank
    beta <- coef(object)
    nvars <- length(beta)
    se <- rep(NA_real_, nvars)
    if (p > 0L) {
        Vp <- chol2inv(object$R[seq_len(p), seq_len(p), drop = FALSE])
        piv <- object$qr$pivot[seq_len(p)]
        se[piv] <- sqrt(diag(Vp))
    }
    disp <- if (is.null(dispersion)) bglm.dispersion.for.se(object) else dispersion
    se <- se * sqrt(disp)
    is.gauss <- object$family$family == "gaussian"
    stat <- beta / se
    nm <- if (is.gauss) c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
         else c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    prob <- if (is.gauss) {
        2 * pt(abs(stat), df = object$df.residual, lower.tail = FALSE)
    } else {
        2 * pnorm(abs(stat), lower.tail = FALSE)
    }
    coefficients <- cbind(beta, se, stat, prob)
    dimnames(coefficients) <- list(names(beta), nm)
    ans <- list(call = object$call, coefficients = coefficients,
        dispersion = disp, deviance = object$deviance,
        df.residual = object$df.residual, df.null = object$df.null,
        null.deviance = object$null.deviance, converged = object$converged,
        boundary = object$boundary, family = object$family, rank = p)
    class(ans) <- "summary.bglm"
    ans
}

print.summary.bglm <- function(x, digits = max(3L, getOption("digits") - 3L),
    ...) {
    cat("Call:\n")
    print(x$call)
    cat("\nFamily:", x$family$family, "(", x$family$link, ")\n")
    cat("\nCoefficients:\n")
    print(format(x$coefficients, digits = digits), quote = FALSE)
    cat("\nDispersion:", format(x$dispersion, digits = digits),
        "\nNull deviance:", format(x$null.deviance, digits = digits), "on",
        x$df.null, "df\nResidual deviance:",
        format(x$deviance, digits = digits), "on", x$df.residual, "df\n")
    if (!x$converged) cat("\nWarning: algorithm did not converge\n")
    if (x$boundary) cat("\nWarning: algorithm stopped at boundary value\n")
    invisible(x)
}

predict.bglm <- function (object, newdata = NULL,
    type = c("link", "response", "terms"), se.fit = FALSE,
    dispersion = NULL, terms = NULL, na.action = na.pass, ...)
{
    type <- match.arg(type)
    na.act <- object$na.action
    object$na.action <- NULL
    if (!se.fit) {
        if (missing(newdata)) {
            pred <- switch(type, link = object$linear.predictors,
                response = object$fitted.values,
                terms = predictLM(object, se.fit = se.fit, scale = 1,
                    type = "terms", terms = terms))
            if (!is.null(na.act))
                pred <- napredict(na.act, pred)
        }
        else {
            pred <- predictLM(object, newdata, se.fit, scale = 1,
                type = ifelse(type == "link", "response", type),
                terms = terms, na.action = na.action)
            switch(type, response = {
                pred <- family(object)$linkinv(pred)
            }, link = , terms = )
        }
    }
    else {
        if (is.null(dispersion) || dispersion == 0)
            dispersion <- summary(object, dispersion = dispersion)$dispersion
        residual.scale <- as.vector(sqrt(dispersion))
        pred <- predictLM(object, newdata, se.fit, scale = residual.scale,
            type = ifelse(type == "link", "response", type),
            terms = terms, na.action = na.action)
        fit <- pred$fit
        se.fit <- pred$se.fit
        switch(type, response = {
            se.fit <- se.fit * abs(family(object)$mu.eta(fit))
            fit <- family(object)$linkinv(fit)
        }, link = , terms = )
        if (missing(newdata) && !is.null(na.act)) {
            fit <- napredict(na.act, fit)
            se.fit <- napredict(na.act, se.fit)
        }
        pred <- list(fit = fit, se.fit = se.fit,
            residual.scale = residual.scale)
    }
    pred
}

## term contributions and fit standard errors, from the fit's own qr
## (arm's predictLM, minus the confidence/prediction-interval
## machinery predict.bglm never requests)
predictLM <- function (object, newdata, se.fit = FALSE, scale = NULL,
    df = Inf, type = c("response", "terms"), terms = NULL,
    na.action = na.pass, ...)
{
    tt <- terms(object)
    keep.order <- object$keep.order
    drop.baseline <- object$drop.baseline
    if (missing(newdata) || is.null(newdata)) {
        mm <- X <- model.matrix(object)
        mmDone <- TRUE
        offset <- object$offset
    }
    else {
        Terms <- delete.response(tt)
        m <- model.frame(Terms, newdata, na.action = na.action,
            xlev = object$xlevels)
        if (!is.null(cl <- attr(Terms, "dataClasses")))
            .checkMFClasses(cl, m)
        X <- bglm_design(Terms, m, contrasts.arg = object$contrasts,
            keep.order = keep.order, drop.baseline = drop.baseline)
        offset <- rep(0, nrow(X))
        if (!is.null(off.num <- attr(tt, "offset")))
            for (i in off.num) offset <- offset + eval(attr(tt,
                "variables")[[i + 1]], newdata)
        if (!is.null(object$call$offset))
            offset <- offset + eval(object$call$offset, newdata)
        mmDone <- FALSE
    }
    n <- length(object$residuals)
    p <- object$rank
    p1 <- seq_len(p)
    piv <- if (p)
        object$qr$pivot[p1]
    if (p < ncol(X) && !(missing(newdata) || is.null(newdata)))
        warning("prediction from a rank-deficient fit may be misleading")
    beta <- object$coefficients
    predictor <- drop(X[, piv, drop = FALSE] %*% beta[piv])
    if (!is.null(offset))
        predictor <- predictor + offset
    type <- match.arg(type)
    if (se.fit) {
        res.var <- if (is.null(scale)) {
            r <- object$residuals
            w <- object$weights
            rss <- sum(if (is.null(w)) r^2 else r^2 * w)
            df <- object$df.residual
            rss/df
        }
        else scale^2
        if (type != "terms") {
            if (p > 0) {
                XRinv <- if (missing(newdata) && is.null(w))
                    qr.Q(getQr(object))[, p1, drop = FALSE]
                else X[, piv] %*% qr.solve(qr.R(getQr(object))[p1, p1])
                ip <- drop(XRinv^2 %*% rep(res.var, p))
            }
            else ip <- rep(0, n)
        }
    }
    if (type == "terms") {
        if (!mmDone) {
            mm <- bglm_design(object$terms, data = model.frame(object),
                keep.order = keep.order, drop.baseline = drop.baseline)
            mmDone <- TRUE
        }
        aa <- attr(mm, "assign")
        ll <- attr(tt, "term.labels")
        hasintercept <- attr(tt, "intercept") > 0L
        if (hasintercept)
            ll <- c("(Intercept)", ll)
        aaa <- factor(aa, labels = ll)
        asgn <- split(order(aa), aaa)
        if (hasintercept) {
            asgn$"(Intercept)" <- NULL
            if (!mmDone) {
                mm <- bglm_design(object$terms, data = model.frame(object),
                    keep.order = keep.order, drop.baseline = drop.baseline)
                mmDone <- TRUE
            }
            avx <- colMeans(mm)
            termsconst <- sum(avx[piv] * beta[piv])
        }
        nterms <- length(asgn)
        if (nterms > 0) {
            predictor <- matrix(ncol = nterms, nrow = NROW(X))
            dimnames(predictor) <- list(rownames(X), names(asgn))
            if (se.fit) {
                ip <- matrix(ncol = nterms, nrow = NROW(X))
                dimnames(ip) <- list(rownames(X), names(asgn))
                Rinv <- qr.solve(qr.R(getQr(object))[p1, p1])
            }
            if (hasintercept)
                X <- sweep(X, 2L, avx, check.margin = FALSE)
            unpiv <- rep.int(0L, NCOL(X))
            unpiv[piv] <- p1
            for (i in seq.int(1L, nterms, length.out = nterms)) {
                iipiv <- asgn[[i]]
                ii <- unpiv[iipiv]
                iipiv[ii == 0L] <- 0L
                predictor[, i] <- if (any(iipiv > 0L))
                    X[, iipiv, drop = FALSE] %*% beta[iipiv]
                else 0
                if (se.fit)
                    ip[, i] <- if (any(iipiv > 0L))
                        as.matrix(X[, iipiv, drop = FALSE] %*% Rinv[ii,
                            , drop = FALSE])^2 %*% rep.int(res.var, p)
                    else 0
            }
            if (!is.null(terms)) {
                predictor <- predictor[, terms, drop = FALSE]
                if (se.fit)
                    ip <- ip[, terms, drop = FALSE]
            }
        }
        else {
            predictor <- ip <- matrix(0, n, 0L)
        }
        attr(predictor, "constant") <- if (hasintercept)
            termsconst
        else 0
    }
    if (se.fit) {
        se <- sqrt(ip)
        if (type == "terms" && !is.null(terms) && !se.fit)
            se <- se[, terms, drop = FALSE]
    }
    if (missing(newdata) && !is.null(na.act <- object$na.action)) {
        predictor <- napredict(na.act, predictor)
        if (se.fit)
            se <- napredict(na.act, se)
    }
    if (se.fit)
        list(fit = predictor, se.fit = se, df = df,
            residual.scale = sqrt(res.var))
    else predictor
}
