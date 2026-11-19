\name{bglm_design}
\alias{bglm_design}
\alias{model.matrixBayes}
\title{Design Matrix with Full-Rank Contrasts and Written Term Order}
\description{
  Builds a model matrix from a formula or a terms object.  Two features
  go beyond \code{stats::model.matrix}: \code{keep.order} keeps columns
  in the order the terms are WRITTEN (the default sorting is
  alphabetical), and \code{drop.baseline = FALSE} codes factors with
  one column per level (full dummies via
  \code{\link{contr.treatment.full}}, full orthogonal polynomials via
  \code{\link{contr.poly.full}}) so a coefficient exists for every
  level, as a Bayesian model with weakly informative priors wants.

  A formula object carries no \code{"intercept"} attribute, and it is
  that attribute (copied onto the terms) that decides whether an
  intercept column is produced: the FORMULA path therefore returns a
  matrix WITHOUT intercept column (callers such as \code{FindIt}'s
  \code{cv.CausalANOVA} prepend their own), while the TERMS path keeps
  the intercept the terms specify.

  A model frame passed as \code{data} is used whole; it is never
  subsetted to the variables the terms name, so transformed terms such
  as \code{poly(psi, degree = 3, raw = TRUE)} keep working.

  \code{model.matrixBayes} is the pinned legacy name of this function
  (reverse dependencies call it by name); it is identical to
  \code{bglm_design}.
}
\usage{
bglm_design(object, data = environment(object),
            contrasts.arg = NULL, xlev = NULL,
            keep.order = FALSE, drop.baseline = FALSE, \dots)

model.matrixBayes(object, data = environment(object),
                  contrasts.arg = NULL, xlev = NULL,
                  keep.order = FALSE, drop.baseline = FALSE, \dots)
}
\arguments{
  \item{object}{a formula or a terms object (or anything
    \code{\link[stats]{terms}} accepts).}
  \item{data}{a data frame, list or environment containing the
    variables, or a model frame built earlier (passed through whole).}
  \item{contrasts.arg}{a named list, as in
    \code{\link[stats]{model.matrix}}.}
  \item{xlev}{the levels of any factors with fixed levels, as in
    \code{\link[stats]{model.matrix}}, used when \code{data} is not
    already a model frame.}
  \item{keep.order}{logical: keep the terms in written order rather
    than sorting them (the order takes effect through the terms object
    that built the model frame).}
  \item{drop.baseline}{logical: if \code{FALSE} (the default) use the
    full-rank contrast generators \code{\link{contr.treatment.full}}
    and \code{\link{contr.poly.full}}; if \code{TRUE} use the
    session's \code{"contrasts"} option (the baseline is dropped,
    \code{n - 1} columns per factor).}
  \item{\dots}{ignored, for compatibility.}
}
\value{
  A numeric design matrix with \code{"assign"} and \code{"contrasts"}
  attributes, as \code{\link[stats]{model.matrix}} returns.
}
\examples{
d <- data.frame(y = rnorm(24),
                g = factor(rep(c("a","b","c"), each = 8)),
                x = rnorm(24))
bglm_design(y ~ g + x, data = d)          # no intercept (formula path)
X <- bglm_design(terms(y ~ g + x, data = d),
                 data = model.frame(y ~ g + x, d))
dim(X)                                    # intercept kept
bglm_design(terms(y ~ g + x, data = d),
            data = model.frame(y ~ g + x, d),
            drop.baseline = TRUE)          # baseline dropped
}
\keyword{models}
