## ===================================================================
## The bglm contrast generators.
##
## Full-rank contrast sets: n columns for n levels (drop.baseline =
## FALSE keeps a coefficient for every level).  Used as the defaults
## by model.matrixBayes(); named by string, so they must stay
## exported to be found by model.matrix.default's get().
## ===================================================================

## Ordered factors: the n orthogonal polynomials (constants, linear,
## quadratic, ...), the full set contr.poly() truncates to n - 1.
contr.poly.full <- function(n, scores = 1:n, contrasts = TRUE)
{
    make.poly <- function(n, scores) {
        y <- scores - mean(scores)
        X <- outer(y, seq_len(n) - 1L, "^")
        QR <- qr(X)
        z <- QR$qr
        z <- z * (row(z) == col(z))
        raw <- qr.qy(QR, z)
        Z <- sweep(raw, 2, apply(raw, 2, function(x) sqrt(sum(x^2))), "/")
        colnames(Z) <- paste("^", seq_len(n) - 1L, sep = "")
        Z
    }
    if (is.numeric(n) && length(n) == 1L) {
        levs <- seq_len(n)
    } else {
        levs <- n
        n <- length(levs)
    }
    if (n < 2L) {
        stop(gettextf("contrasts not defined for %d degrees of freedom",
                      n - 1L), domain = NA)
    }
    if (n > 95L) {
        stop(gettextf(
            paste("orthogonal polynomials cannot be represented",
                  "accurately enough for %d degrees of freedom"),
            n - 1L), domain = NA)
    }
    if (length(scores) != n) {
        stop("'scores' argument is of the wrong length")
    }
    if (!is.numeric(scores) || any(duplicated(scores))) {
        stop("'scores' must all be different numbers")
    }
    contr <- make.poly(n, scores)
    if (contrasts) {
        dn <- colnames(contr)
        dn[2:min(4L, n)] <- c(".L", ".Q", ".C")[seq_len(min(3L, n - 1L))]
        colnames(contr) <- dn
        contr[, , drop = FALSE]
    } else {
        contr[, 1L] <- 1
        contr
    }
}

## Unordered factors: the n-by-n identity -- one indicator column per
## level, no baseline dropped.  (base::contr.treat drops a column;
## this is the same dummy coding before the drop.)
contr.treatment.full <- function(n, base = 1, contrasts = TRUE)
{
    if (is.numeric(n) && length(n) == 1L) {
        if (n > 1L) {
            levs <- seq_len(n)
        } else {
            stop("not enough degrees of freedom to define contrasts")
        }
    } else {
        levs <- n
        n <- length(n)
    }
    contr <- array(0, c(n, n), list(levs, levs))
    diag(contr) <- 1
    if (contrasts) {
        if (n < 2L) {
            stop(gettextf("contrasts not defined for %d degrees of freedom",
                          n - 1L), domain = NA)
        }
        if (base < 1L | base > n) {
            stop("baseline group number out of range")
        }
        contr <- contr[, , drop = FALSE]
    }
    contr
}

## ===================================================================
## Pinned legacy names.  arm's model.matrixBayes (and therefore any
## saved arm fit) stores these strings in the "contrasts" attribute,
## and reverse dependencies resolve them by name; they stay exported
## as thin aliases of the primary generators above.
## ===================================================================

contr.bayes.unordered <- contr.treatment.full
contr.bayes.ordered   <- contr.poly.full
