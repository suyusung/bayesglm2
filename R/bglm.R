## ===================================================================
## bglm -- Bayesian generalized linear model by penalized IRLS,
## written from scratch: no glm(), no glm.fit(), no lm.fit().  The
## linear algebra goes through the package's own pivoted-QR kernel
## (src/bglm_qr.c); the statistical structure is the scaled-t MAP
## estimation of Gelman & Hill (arm's bayesglm):
##
##   beta_j ~ scaled-t(prior.mean_j, prior.scale_j, prior.df_j),
##   inference by MAP = iteratively reweighted least squares on the
##   EM scale-mixture representation: each prior coefficient enters
##   as one pseudo-observation row of the design (data augmentation),
##   whose weight is updated every iteration from the current fit
##   (the t as a scale mixture of normals).
##
## Convergence is the same pair of relative criteria the old engine
## used: deviance and (for non-binomial/poisson families) the
## dispersion, both with the 0.1 + |.| floor.
## ===================================================================

bglm <- function (formula, family = gaussian, data, weights, subset,
    na.action, start = NULL, etastart, mustart, offset, control = list(...),
    model = TRUE, method = "bglm.fit", x = FALSE, y = TRUE, contrasts = NULL,
    drop.unused.levels = TRUE,
    prior.mean = 0, prior.scale = NULL, prior.df = 1,
    prior.mean.for.intercept = 0,
    prior.scale.for.intercept = NULL,
    prior.df.for.intercept = 1,
    min.prior.scale = 1e-12,
    scaled = TRUE,
    keep.order = TRUE, drop.baseline = TRUE,
    maxit = 100, print.unnormalized.log.posterior = FALSE,
    Warning = TRUE, ...)
{
    call <- match.call()
    if (is.character(family)) {
        family <- get(family, mode = "function", envir = parent.frame())
    }
    if (is.function(family)) {
        family <- family()
    }
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
    if (missing(data)) {
        data <- environment(formula)
    }
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action",
        "etastart", "mustart", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- drop.unused.levels
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    if (identical(method, "model.frame")) {
        return(mf)
    }
    if (!is.character(method) && !is.function(method)) {
        stop("invalid 'method' argument")
    }
    if (missing(maxit) && !is.null(control$maxit)) {
        maxit <- control$maxit
    }
    if (identical(method, "bglm.fit")) {
        control <- do.call("glm.control", control)
    }
    control$maxit <- maxit

    mt <- attr(mf, "terms")
    # The design matrix is built with keep.order, so the terms stored in
    # the fit must be in that order too: everything that rebuilds the
    # matrix from them -- model.matrix(fit), predict(se.fit = TRUE),
    # predict(type = "terms") -- otherwise pairs columns with the wrong
    # coefficients.  Building the *model frame* from a terms object
    # instead breaks formulas whose terms are transformations, e.g.
    # poly(psi, ...): the frame carries the evaluated column, not the
    # variable the terms refer to.
    mt <- terms(as.formula(formula), data = data, keep.order = keep.order)
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1L) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm)) {
            names(Y) <- nm
        }
    }
    X <- if (!is.empty.model(mt)) {
        bglm_design(object = mt, data = mf, contrasts.arg = contrasts,
            keep.order = keep.order, drop.baseline = drop.baseline)
    } else {
        matrix(, NROW(Y), 0L)
    }
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights)) {
        stop("'weights' must be a numeric vector")
    }
    if (!is.null(weights) && any(weights < 0)) {
        stop("negative weights not allowed")
    }
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
        if (length(offset) != NROW(Y))
            stop(gettextf("number of offsets is %d should equal %d (number of
            observations)",
                length(offset), NROW(Y)), domain = NA)
    }
    mustart <- model.extract(mf, "mustart")
    etastart <- model.extract(mf, "etastart")
    fit <- bglm.fit(x = X, y = Y, weights = weights, start = start,
        etastart = etastart, mustart = mustart, offset = offset,
        family = family, control = control,
        intercept = attr(mt, "intercept") > 0L,
        prior.mean = prior.mean,
        prior.scale = prior.scale,
        prior.df = prior.df,
        prior.mean.for.intercept = prior.mean.for.intercept,
        prior.scale.for.intercept = prior.scale.for.intercept,
        prior.df.for.intercept = prior.df.for.intercept,
        min.prior.scale = min.prior.scale,
        print.unnormalized.log.posterior = print.unnormalized.log.posterior,
        scaled = scaled, Warning = Warning)
    if (length(offset) && attr(mt, "intercept") > 0L) {
        # the null deviance is the offset-only model: same data, an
        # intercept alone, the intercept's prior
        fit2 <- bglm.fit(x = X[, "(Intercept)", drop = FALSE], y = Y,
            weights = weights, offset = offset, family = family,
            control = control, intercept = TRUE,
            prior.mean.for.intercept = prior.mean.for.intercept,
            prior.scale.for.intercept = prior.scale.for.intercept,
            prior.df.for.intercept = prior.df.for.intercept,
            min.prior.scale = min.prior.scale,
            print.unnormalized.log.posterior = print.unnormalized.log.posterior,
            scaled = scaled, Warning = Warning)
        if (!fit2$converged && Warning) {
            warning("fitting to calculate the null deviance did not converge --
            increase 'maxit'?")
        }
        fit$null.deviance <- fit2$deviance
    }
    if (model) {
        fit$model <- mf
    }
    fit$na.action <- attr(mf, "na.action")
    if (x) {
        fit$x <- X
    }
    if (!y) {
        fit$y <- NULL
    }
    fit <- c(fit, list(call = call, formula = formula, terms = mt,
        data = data, offset = offset, control = control, method = method,
        contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, mf)),
        keep.order = keep.order, drop.baseline = drop.baseline)
    class(fit) <- "bglm"
    return(fit)
}


## lm.fit's shape, from our own pivoted-QR kernel.  Same tolerance
## (1e-7) as lm.fit so rank decisions -- which columns are declared
## aliased -- match it exactly (pinned by tests/qr-kernel-vs-lm-fit.R).
bglm_qr_solve <- function(X, y, tol = 1e-7, w = NULL) {
    X <- as.matrix(X); storage.mode(X) <- "double"
    y <- as.vector(y); storage.mode(y) <- "double"
    if (!is.null(w)) {
        w <- as.vector(w); storage.mode(w) <- "double"
    }
    r <- .Call(C_bglm_qrls, X, y, tol, w)  # coef, effects, residuals, rank, qr, qraux, pivot
    rank <- r[[4]]
    coef <- r[[1]]; nms <- colnames(X)
    ## lm.fit (R >= 4.3) reports NA, not 0, for aliased columns:
    if (rank < ncol(X)) coef[r[[7]][(rank + 1):ncol(X)]] <- NA_real_
    names(coef) <- nms
    colnames(r[[5]]) <- nms
    list(coefficients = coef,
         effects = r[[2]],
         residuals = r[[3]],
         rank = rank,
         qr = structure(list(qr = r[[5]], qraux = r[[6]],
                             pivot = r[[7]], tol = tol),
                        class = "qr"))
}


bglm.fit <- function (x, y, weights = rep(1, nobs), start = NULL,
    etastart = NULL, mustart = NULL, offset = rep(0, nobs),
    family = gaussian(), control = list(), intercept = TRUE,
    prior.mean = 0, prior.scale = NULL, prior.df = 1,
    prior.mean.for.intercept = 0,
    prior.scale.for.intercept = NULL,
    prior.df.for.intercept = 1,
    min.prior.scale = 1e-12, scaled = TRUE,
    print.unnormalized.log.posterior = FALSE,
    Warning = TRUE)
{
    control <- do.call("glm.control", control)

    x <- as.matrix(x)
    xnames <- dimnames(x)[[2L]]
    ynames <- if (is.matrix(y)) {
        rownames(y)
    } else {
        names(y)
    }
    conv <- FALSE
    nobs <- NROW(y)
    nvars <- NCOL(x)

    #=============================
    #  initialize priors
    #=============================
    if (is.null(prior.scale)) {
        prior.scale <- 2.5
        if (family$link == "probit") {
            prior.scale <- prior.scale * 1.6
        }
    }
    if (is.null(prior.scale.for.intercept)) {
        prior.scale.for.intercept <- 10
        if (family$link == "probit") {
            prior.scale.for.intercept <- prior.scale.for.intercept * 1.6
        }
    }
    if (intercept) {
        nvars <- nvars - 1
    }
    if (length(prior.mean) == 1L) {
        prior.mean <- rep(prior.mean, nvars)
    } else if (length(prior.mean) != nvars) {
        stop("invalid length for prior.mean")
    }
    if (length(prior.scale) == 1L) {
        prior.scale <- rep(prior.scale, nvars)
    } else if (length(prior.scale) != nvars) {
        stop("invalid length for prior.scale")
    }
    if (length(prior.df) == 1L) {
        prior.df <- rep(prior.df, nvars)
    } else if (length(prior.df) != nvars) {
        stop("invalid length for prior.df")
    }
    if (intercept) {
        prior.mean <- c(prior.mean.for.intercept, prior.mean)
        prior.scale <- c(prior.scale.for.intercept, prior.scale)
        prior.df <- c(prior.df.for.intercept, prior.df)
    }

    # A weight of 2 has to mean the same as a duplicated row, and the
    # prior is scaled by the spread of the predictors and of the
    # response, so those spreads are computed with the weights.
    wts <- if (is.null(weights)) rep.int(1, nobs) else weights
    swts <- sum(wts)
    wtd.sd <- function(v) {
        if (anyNA(v)) {
            keep <- !is.na(v)
            v <- v[keep]
            w <- wts[keep]
            sw <- sum(w)
        } else {
            w <- wts
            sw <- swts
        }
        if (length(v) < 2L || sw <= 1) {
            return(sd(v))
        }
        m <- sum(w * v)/sw
        sqrt(sum(w * (v - m)^2)/(sw - 1))
    }
    if (scaled) {
        if (family$family == "gaussian") {
            prior.scale <- prior.scale * 2 * wtd.sd(y)
        }
        if (nvars == 0 && intercept) {
            nvars <- 1                       # intercept-only: loop over it
        } else if (intercept) {
            nvars <- nvars + 1
        }
        for (j in 1:nvars) {
            x.obs <- x[,j]
            if (anyNA(x.obs)) x.obs <- x.obs[!is.na(x.obs)]
            # the prior only distinguishes < 2, == 2, > 2 distinct
            # values; a capped min/max scan decides that without the
            # O(n) hash table of length(unique(.)) (unique() treats
            # -0 and 0 as one value, same as ==)
            if (length(x.obs) < 2L) {
                num.categories <- length(x.obs)
            } else {
                mn <- min(x.obs)
                mx <- max(x.obs)
                num.categories <- if (mn == mx)
                    1L
                else if (all(x.obs == mn | x.obs == mx))
                    2L
                else
                    3L
            }
            x.scale <- 1
            if (num.categories == 2L) {
                x.scale <- max(x.obs) - min(x.obs)
            } else if (num.categories > 2) {
                x.scale <- 2 * wtd.sd(x[,j])
            }
            prior.scale[j] <- prior.scale[j]/x.scale
            if (prior.scale[j] < min.prior.scale) {
                prior.scale[j] <- min.prior.scale
                warning("prior scale for variable ", j,
                    " set to min.prior.scale = ", min.prior.scale, "\n")
            }
        }
    }

    #===================
    nvars <- NCOL(x)
    EMPTY <- nvars == 0
    if (is.null(weights))
        weights <- rep.int(1, nobs)
    if (is.null(offset))
        offset <- rep.int(0, nobs)
    variance <- family$variance
    linkinv <- family$linkinv
    if (!is.function(variance) || !is.function(linkinv))
        stop("'family' argument seems not to be a valid family object",
            call. = FALSE)
    dev.resids <- family$dev.resids
    aic <- family$aic
    mu.eta <- family$mu.eta
    unless.null <- function(x, if.null) {
        if (is.null(x)) if.null else x
    }
    valideta <- unless.null(family$valideta, function(eta) TRUE)
    validmu <- unless.null(family$validmu, function(mu) TRUE)
    if (is.null(mustart)) {
        eval(family$initialize)              # defines n, adjusts y
    } else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
    }
    if (EMPTY) {
        eta <- rep.int(0, nobs) + offset
        if (!valideta(eta))
            stop("invalid linear predictor values in empty model",
                call. = FALSE)
        mu <- linkinv(eta)
        if (!validmu(mu))
            stop("invalid fitted means in empty model", call. = FALSE)
        dev <- sum(dev.resids(y, mu, weights))
        w <- ((weights * mu.eta(eta)^2)/variance(mu))^0.5
        residuals <- (y - mu)/mu.eta(eta)
        good <- rep_len(TRUE, length(residuals))
        boundary <- conv <- TRUE
        coef <- numeric()
        iter <- 0L
    } else {
        coefold <- NULL
        eta <- if (!is.null(etastart)) {
            etastart
        } else if (!is.null(start)) {
            if (length(start) != nvars) {
                if (length(start) == 1L && start == 0) {
                    start <- rep(0, nvars)
                    offset + as.vector(ifelse((NCOL(x) == 1L), x * start,
                        x %*% start))
                } else {
                    stop(gettextf(
                        "length of 'start' should equal %d and correspond to initial coefs for %s",
                        nvars, paste(deparse(xnames), collapse = ", ")),
                        domain = NA)
                }
            } else {
                coefold <- start
                offset + as.vector(if (NCOL(x) == 1L) x * start else x %*% start)
            }
        } else {
            family$linkfun(mustart)
        }
        mu <- linkinv(eta)
        if (!(validmu(mu) && valideta(eta)))
            stop("cannot find valid starting values: please specify some",
                call. = FALSE)
        devold <- sum(dev.resids(y, mu, weights))
        boundary <- conv <- FALSE

        #======================================
        #   initialize prior.sd
        #======================================
        prior.sd <- prior.scale

        ## the augmented design is constant across EM iterations: the
        ## data block is x, the p prior pseudo-rows are the identity
        ## (the intercept row carrying the weighted predictor means).
        ## Built once; the per-iteration row weighting is fused into
        ## the QR kernel.
        x.star <- rbind(x, diag(NCOL(x)))
        if (intercept && scaled) {
            x.star[nobs + 1,] <- colSums(x * wts)/sum(wts)
        }

        #=====================================
        dispersion <- ifelse((family$family %in% c("poisson", "binomial")),
            1, var(y)/10000)
        dispersionold <- dispersion
        # the rank fix-up after the loop indexes `coef`; if the loop
        # breaks before assigning it (non-finite coefficients at the
        # first iteration) the name must still be a vector.
        coef <- rep(NA_real_, NCOL(x))
        dev <- devold
        for (iter in 1L:control$maxit) {
            good1 <- weights > 0
            vm <- variance(mu)
            varmu <- vm[good1]
            if (anyNA(varmu))
                stop("NAs in V(mu)")
            if (any(varmu == 0))
                stop("0s in V(mu)")
            mu.eta.val <- mu.eta(eta)
            if (any(is.na(mu.eta.val[good1])))
                stop("NAs in d(mu)/d(eta)")
            good <- good1 & (mu.eta.val != 0)
            if (all(!good)) {
                conv <- FALSE
                warning(gettextf("no observations informative at iteration %d",
                    iter), domain = NA)
                break
            }
            z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
            w <- sqrt((weights[good] * mu.eta.val[good]^2)/vm[good])

            #======================
            #  data augmentation
            #======================
            # one prior pseudo-row per coefficient (the scale-mixture
            # representation of the t); the augmented design x.star is
            # constant across iterations (built above the loop, with
            # the intercept pseudo-row at the weighted predictor
            # means), and the per-iteration row weights are applied
            # inside the QR kernel -- same products, no n x p
            # temporaries per iteration.
            z.star <- c(z, prior.mean)
            w.star <- c(w, sqrt(dispersion)/prior.sd)
            Xs <- if (all(good))
                x.star
            else rbind(x[good, , drop = FALSE],
                x.star[nobs + seq_len(NCOL(x)), , drop = FALSE])
            fit <- bglm_qr_solve(Xs, z.star, w = w.star)
            if (any(!is.finite(fit$coefficients))) {
                # The prior's pseudo-rows carry no weight here -- a
                # prior flat enough that the augmented design is
                # numerically singular.  A flat prior means the prior
                # is doing nothing, so the right limit is the
                # maximum-likelihood fit: run THIS engine's loop again
                # with a flat prior (never an external fitter).
                if (Warning) {
                    warning(gettextf("non-finite coefficients at iteration %d",
                        iter), domain = NA)
                }
                fit <- bglm.mle.fit(x, y, weights = weights,
                    offset = offset, family = family, control = control,
                    intercept = intercept)
                conv <- fit$converged
                coef <- fit$coefficients
                eta <- fit$linear.predictors
                mu <- fit$fitted.values
                dev <- fit$deviance
                break
            }
            start[fit$qr$pivot] <- coefs.hat <- fit$coefficients
            fit$qr$qr <- as.matrix(fit$qr$qr)
            V.coefs <- chol2inv(fit$qr$qr[1:NCOL(x.star),
                1:NCOL(x.star), drop = FALSE])
            # EM update of the t scale mixture: the prior sd for the
            # next iteration from the current fit and its variance
            prior.sd <- ifelse(prior.df == Inf, prior.scale,
                sqrt(((coefs.hat - prior.mean)^2 + diag(V.coefs)*dispersion +
                prior.df * prior.scale^2)/(1 + prior.df)))
            start[fit$qr$pivot] <- fit$coefficients
            eta <- drop(x %*% start)
            mu <- linkinv(eta <- eta + offset)
            dev <- sum(dev.resids(y, mu, weights))
            if (!(family$family %in% c("poisson", "binomial"))) {
                x.good <- x[good, , drop = FALSE]
                # averaged over the weighted sample, so that a weight
                # of 2 is the same as a duplicated row here too
                w.good <- weights[good]
                mse.resid <- sum(w.good * (w * (z - x.good %*%
                    coefs.hat))^2)/sum(w.good)
                mse.uncertainty <- sum(w.good *
                    rowSums((x.good %*% V.coefs) * x.good))/sum(w.good) *
                    dispersion
                dispersion <- mse.resid + mse.uncertainty
            }
            if (control$trace)
                cat("Deviance = ", dev, " Iterations - ", iter,
                    "\n", sep = "")
            boundary <- FALSE
            if (!is.finite(dev)) {
                if (is.null(coefold))
                    stop("no valid set of coefficients has been found: please supply starting values",
                        call. = FALSE)
                warning("step size truncated due to divergence",
                    call. = FALSE)
                ii <- 1
                while (!is.finite(dev)) {
                    if (ii > control$maxit)
                        stop("inner loop 1; cannot correct step size",
                            call. = FALSE)
                    ii <- ii + 1
                    start <- (start + coefold)/2
                    eta <- drop(x %*% start)
                    mu <- linkinv(eta <- eta + offset)
                    dev <- sum(dev.resids(y, mu, weights))
                }
                boundary <- TRUE
                if (control$trace)
                    cat("Step halved: new deviance = ", dev, "\n", sep = "")
            }
            if (!(valideta(eta) && validmu(mu))) {
                    if (is.null(coefold))
                        stop("no valid set of coefficients has been found: please supply starting values",
                            call. = FALSE)
                    warning("step size truncated: out of bounds",
                        call. = FALSE)
                    ii <- 1
                    while (!(valideta(eta) && validmu(mu))) {
                        if (ii > control$maxit)
                            stop("inner loop 2; cannot correct step size",
                                call. = FALSE)
                        ii <- ii + 1
                        start <- (start + coefold)/2
                        eta <- drop(x %*% start)
                        mu <- linkinv(eta <- eta + offset)
                    }
                    boundary <- TRUE
                    dev <- sum(dev.resids(y, mu, weights))
                    if (control$trace)
                        cat("Step halved: new deviance = ", dev, "\n", sep = "")
            }
            #===============================
            # print unnormalized log posterior
            #================================
            if (family$family == "binomial" && print.unnormalized.log.posterior) {
                # the prior is a scaled t: location prior.mean, scale
                # prior.scale, df prior.df.  dt()'s third argument is
                # the non-centrality, not the location, and it has no
                # scale argument at all.
                logprior <- sum(dt((coefs.hat - prior.mean)/prior.scale,
                    prior.df, log = TRUE) - log(prior.scale))
                # the link is the family's, not always logit; weights
                # are the binomial totals, so this is the log
                # likelihood glm reports
                xb <- linkinv(drop(x %*% coefs.hat) + offset)
                loglikelihood <- sum(weights * (y * log(xb) +
                    (1 - y) * log(1 - xb)))
                cat("log prior: ", logprior, ", log likelihood: ",
                    loglikelihood, ", unnormalized log posterior: ",
                    loglikelihood + logprior, "\n", sep = "")
            }
            #================================

            if (iter > 1 & abs(dev - devold)/(0.1 + abs(dev)) <
                    control$epsilon & abs(dispersion - dispersionold)/(0.1 +
                    abs(dispersion)) < control$epsilon) {
                    conv <- TRUE
                    coef <- start
                    break
            } else {
                devold <- dev
                dispersionold <- dispersion
                coef <- coefold <- start
            }
        }

        if (!conv && Warning) {
            warning("algorithm did not converge", call. = FALSE)
        }
        if (boundary && Warning) {
            warning("algorithm stopped at boundary value", call. = FALSE)
        }
        eps <- 10 * .Machine$double.eps
        if (family$family == "binomial") {
            if ((any(mu > 1 - eps) || any(mu < eps)) && Warning) {
                warning("fitted probabilities numerically 0 or 1 occurred",
                    call. = FALSE)
            }
        }
        if (family$family == "poisson") {
            if (any(mu < eps) && Warning) {
                warning("fitted rates numerically 0 occurred", call. = FALSE)
            }
        }
        if (fit$rank < nvars) {
            coef[fit$qr$pivot][seq.int(fit$rank + 1, nvars)] <- NA
        }
        xxnames <- xnames[fit$qr$pivot]
        residuals <- rep.int(NA, nobs)
        ## (y - mu) / mu'(eta) at the FINAL fit.  The old engine reused
        ## the IRLS working response z, which is computed before the
        ## coefficient update, so its stored residuals mixed z(iter-1)
        ## with eta(iter) -- an artifact that only cancels for the
        ## gaussian link.  Recompute from the consistent (eta, mu) pair
        ## left by the last update instead.
        residuals[good] <- (y - mu)[good] / mu.eta(eta)[good]
        fit$qr$qr <- as.matrix(fit$qr$qr)
        nr <- min(sum(good), nvars)
        if (nr < nvars) {
            Rmat <- diag(nvars)
            Rmat[1L:nr, 1L:nvars] <- fit$qr$qr[1L:nr, 1L:nvars]
        } else Rmat <- fit$qr$qr[1L:nvars, 1L:nvars]
        Rmat <- as.matrix(Rmat)
        Rmat[row(Rmat) > col(Rmat)] <- 0
        names(coef) <- xnames
        colnames(fit$qr$qr) <- xxnames
        dimnames(Rmat) <- list(xxnames, xxnames)
    }
    names(residuals) <- ynames
    names(mu) <- ynames
    names(eta) <- ynames
    wt <- rep.int(0, nobs)
    wt[good] <- w^2
    names(wt) <- ynames
    names(weights) <- ynames
    names(y) <- ynames
    wtdmu <- if (intercept) {
        sum(weights * y)/sum(weights)
    } else {
        linkinv(offset)
    }
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok - as.integer(intercept)
    rank <- if (EMPTY) {
        0
    } else {
        fit$rank
    }
    resdf <- n.ok - rank
    aic.model <- aic(y, n, mu, weights, dev) + 2 * rank
    list(coefficients = coef,
        residuals = residuals,
        fitted.values = mu,
        effects = if (!EMPTY) fit$effects,
        R = if (!EMPTY) Rmat,
        rank = rank,
        qr = if (!EMPTY) structure(getQr(fit)[c("qr", "rank", "qraux",
            "pivot", "tol")], class = "qr"),
        family = family,
        linear.predictors = eta,
        deviance = dev,
        aic = aic.model,
        null.deviance = nulldev,
        iter = iter,
        weights = wt,
        prior.weights = weights,
        df.residual = resdf,
        df.null = nulldf,
        y = y,
        converged = conv,
        boundary = boundary,
        prior.mean = prior.mean,
        prior.scale = prior.scale,
        prior.df = prior.df,
        prior.sd = prior.sd,
        dispersion = dispersion)
}

getQr <- function (x, ...) {
    if (is.null(r <- x$qr))
        stop("model fit does not have a proper 'qr' component. Rank zero?")
    r
}


## --------------------------------------------------------------------
## bglm.mle.fit -- internal plain maximum-likelihood IRLS, the
## glm.fit() of this engine (no glm/glm.fit anywhere).  Reached only
## when the Bayesian augmented system is numerically singular because
## the prior is flat enough that its pseudo-rows carry no weight
## (bayesglm's rank-deficient + flat-prior path).  Deliberately NOT
## bglm.fit(): calling that recurses forever -- the flat prior does not
## remove the singularity that sent us here.  Same IRLS kernel, plain
## deviance convergence (glm.fit's epsilon criterion with the
## 0.1+|dev| floor), and the lm.fit-style rank handling: the pivoted
## QR solves on the estimable columns, the aliased coefficients stay
## 0 while iterating and are reported NA after the loop, as glm does.
## --------------------------------------------------------------------
bglm.mle.fit <- function(x, y, weights = rep(1, nobs), offset = rep(0, nobs),
    family = gaussian(), control = glm.control(), intercept = TRUE)
{
    nobs <- NROW(y); nvars <- NCOL(x)
    EMPTY <- (nobs == 0L || nvars == 0L)
    dev.resids <- family$dev.resids
    aic <- family$aic
    mu.eta <- family$mu.eta
    linkinv <- family$linkinv
    variance <- family$variance
    unless.null <- function(z, if.null) if (is.null(z)) if.null else z
    valideta <- unless.null(family$valideta, function(eta) TRUE)
    validmu <- unless.null(family$validmu, function(mu) TRUE)
    eval(family$initialize)              # defines n and adjusts y for some families
    if (EMPTY)
        stop("bglm.mle.fit does not take empty models; the engine handles those")
    eta <- family$linkfun(mustart)
    mu <- linkinv(eta)
    if (!(validmu(mu) && valideta(eta)))
        stop("cannot find valid starting values: please specify some",
            call. = FALSE)
    devold <- sum(dev.resids(y, mu, weights))
    conv <- FALSE
    boundary <- FALSE
    coef <- rep(NA_real_, nvars)
    dev <- devold
    for (iter in 1L:control$maxit) {
        good <- (weights > 0)
        varmu <- variance(mu)[good]
        if (anyNA(varmu)) stop("NAs in V(mu)")
        if (any(varmu == 0)) stop("0s in V(mu)")
        mu.eta.val <- mu.eta(eta)
        if (any(is.na(mu.eta.val[good]))) stop("NAs in d(mu)/d(eta)")
        good <- (weights > 0) & (mu.eta.val != 0)
        if (all(!good)) {
            conv <- FALSE
            warning(gettextf("no observations informative at iteration %d",
                iter), domain = NA)
            break
        }
        z <- (eta - offset)[good] + (y - mu)[good] / mu.eta.val[good]
        w <- sqrt((weights[good] * mu.eta.val[good]^2) / varmu)
        fit <- bglm_qr_solve(x[good, , drop = FALSE], z, w = w)
        coefs <- fit$coefficients
        ## aliased columns: lm.fit reports NA; iterate on 0 for them
        coefs[is.na(coefs)] <- 0
        start <- coefs
        eta <- drop(x %*% start)
        mu <- linkinv(eta <- eta + offset)
        dev <- sum(dev.resids(y, mu, weights))
        if (!is.finite(dev)) {
            if (is.null(coefold))
                stop("no valid set of coefficients has been found: please supply starting values",
                    call. = FALSE)
            warning("step size truncated due to divergence", call. = FALSE)
            ii <- 1
            while (!is.finite(dev)) {
                if (ii > control$maxit)
                    stop("inner loop 1; cannot correct step size", call. = FALSE)
                ii <- ii + 1
                start <- (start + coefold) / 2
                eta <- drop(x %*% start)
                mu <- linkinv(eta <- eta + offset)
                dev <- sum(dev.resids(y, mu, weights))
            }
            boundary <- TRUE
        }
        if (!(valideta(eta) && validmu(mu))) {
            if (is.null(coefold))
                stop("no valid set of coefficients has been found: please supply starting values",
                    call. = FALSE)
            warning("step size truncated: out of bounds", call. = FALSE)
            ii <- 1
            while (!(valideta(eta) && validmu(mu))) {
                if (ii > control$maxit)
                    stop("inner loop 2; cannot correct step size", call. = FALSE)
                ii <- ii + 1
                start <- (start + coefold) / 2
                eta <- drop(x %*% start)
                mu <- linkinv(eta <- eta + offset)
            }
            boundary <- TRUE
            dev <- sum(dev.resids(y, mu, weights))
        }
        if (iter > 1 && abs(dev - devold) / (0.1 + abs(dev)) < control$epsilon) {
            conv <- TRUE
            coef <- start
            break
        } else {
            devold <- dev
            coef <- coefold <- start
        }
    }
    if (fit$rank < nvars)
        coef[fit$qr$pivot][seq.int(fit$rank + 1, nvars)] <- NA
    list(coefficients = coef,
         linear.predictors = eta,
         fitted.values = mu,
         deviance = dev,
         converged = conv,
         rank = fit$rank,
         qr = fit$qr,
         iter = iter,
         boundary = boundary)
}

