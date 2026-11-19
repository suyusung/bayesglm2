## ===================================================================
## model.matrixBayes -- the design-matrix layer, rewritten clean.
##
## Builds a model matrix from a formula or a terms object with two
## behaviours arm's bayesglm (and its reverse dependencies) rely on:
##
##   * keep.order    -- columns follow the terms as WRITTEN
##                      (y ~ x:z + x gives x:z before x;
##                      stats::terms() would sort them alphabetically)
##   * drop.baseline -- FALSE builds FULL dummy columns for unordered
##                      factors (contr.treatment.full) and the full
##                      orthogonal-polynomial set for ordered ones
##                      (contr.poly.full), one column per level,
##                      so a coefficient exists for every level.
##
## And one contract that reads like a bug but is pinned by reverse
## dependencies (FindIt's cv.CausalANOVA calls
##   model.matrixBayes(formula, data)
## and prepends its own column of ones, so it must come back WITHOUT
## an intercept column):
##
##   * a FORMULA object carries no "intercept" attribute; the
##     attribute is copied onto the terms, so the formula path
##     produces no intercept column.  A TERMS object carries one and
##     the terms path keeps it.  This is not an accident to fix here
##     -- it is the documented interface of the formula path.
##
## A model frame handed in as `data` is passed to model.matrix.default
## WHOLE, never subsetted to the variables the terms name: a frame
## holds the evaluated column of a transformed term such as
## poly(psi, degree = 3, raw = TRUE), not the raw `psi`, and
## subsetting drops the frame's "terms" attribute, which makes
## model.matrix.default try to rebuild the frame from what is left
## and die with "object 'psi' not found" (qgcomp, issue #54-style
## formulas).
## ===================================================================

bglm_design <- function(object, data = environment(object),
                              contrasts.arg = NULL, xlev = NULL,
                              keep.order = FALSE, drop.baseline = FALSE,
                              ...)
{
    t <- if (missing(data)) {
        terms(object)
    } else {
        stats::terms.formula(object, data = data, keep.order = keep.order)
    }
    ## The FindIt contract, stated as code: NULL for a formula (no
    ## intercept column -- the caller adds its own ones), 1/0 for a
    ## terms object (intercept kept as the terms specify).
    attr(t, "intercept") <- attr(object, "intercept")

    mf <- if (!is.null(attr(data, "terms"))) {
        ## Already a model frame: pass it whole (see header comment).
        data
    } else {
        stats::model.frame(object, data, xlev = xlev)
    }

    int <- attr(t, "response")
    if (length(mf)) {                       # otherwise no rhs terms
        mf <- convert_character_columns(mf)
        isF <- vapply(mf, function(x) is.factor(x) || is.logical(x), NA)
        isF[int] <- FALSE                   # the response gets no contrasts
        isOF <- vapply(mf, is.ordered, NA)
        contr.funs <- default_contrast_funs(drop.baseline)
        for (nn in names(mf)[isF]) {
            if (is.null(attr(mf[[nn]], "contrasts"))) {
                stats::contrasts(mf[[nn]]) <- contr.funs[1L + isOF[nn]]
            }
        }
        mf <- apply_contrasts_arg(mf, contrasts.arg)
    } else {
        ## model.matrix.default wants at least one variable to size from
        isF <- FALSE
        mf <- data.frame(x = rep(0, NROW(mf)))
    }

    ans <- stats::model.matrix.default(object = t, data = mf)
    attr(ans, "contrasts") <- if (any(isF)) {
        lapply(mf[isF], function(x) attr(x, "contrasts"))
    } else {
        NULL
    }
    ans
}

## Full dummies / full orthogonal polynomials, or the session defaults
## when the baseline is dropped (the plain glm-like coding).
default_contrast_funs <- function(drop.baseline) {
    if (drop.baseline) {
        as.character(getOption("contrasts"))
    } else {
        c("contr.treatment.full", "contr.poly.full")
    }
}

## Character columns become factors, with the same warning base has
## always given.  (Kept out of the main flow so that flow reads as the
## four steps it is.)
convert_character_columns <- function(mf) {
    for (nm in names(mf)) {
        if (is.character(mf[[nm]])) {
            mf[[nm]] <- factor(mf[[nm]])
            warning(gettextf("variable '%s' converted to a factor", nm),
                    domain = NA)
        }
    }
    mf
}

## contrasts.arg: a named list, per variable a contrast function name
## or a contrast matrix.  Absent variables warn and are ignored, as in
## model.matrix.default.
apply_contrasts_arg <- function(mf, contrasts.arg) {
    if (is.null(contrasts.arg)) {
        return(mf)
    }
    if (!is.list(contrasts.arg)) {
        stop("invalid 'contrasts.arg' argument")
    }
    namD <- names(mf)
    namC <- names(contrasts.arg)
    if (is.null(namC)) {
        stop("invalid 'contrasts.arg' argument")
    }
    for (nn in namC) {
        ni <- match(nn, namD)
        if (is.na(ni)) {
            warning(gettextf(
                "variable '%s' is absent, its contrast will be ignored", nn),
                domain = NA)
        } else {
            ca <- contrasts.arg[[nn]]
            if (is.matrix(ca)) {
                stats::contrasts(mf[[ni]], ncol(ca)) <- ca
            } else {
                stats::contrasts(mf[[ni]]) <- ca
            }
        }
    }
    mf
}

## ===================================================================
## Pinned legacy name.  FindIt's cv.CausalANOVA calls
## arm:::model.matrixBayes(formula, data) BY NAME, and the no-intercept
## formula-path contract travels with that name; it stays exported as
## an alias of bglm_design().
## ===================================================================

model.matrixBayes <- bglm_design
