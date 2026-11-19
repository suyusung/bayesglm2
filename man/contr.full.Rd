\name{contr.full}
\alias{contr.full}
\alias{contr.treatment.full}
\alias{contr.poly.full}
\alias{contr.bayes.unordered}
\alias{contr.bayes.ordered}
\title{Full-Rank Contrast Generators}
\description{
  Contrast generators with \code{n} columns for \code{n} levels, used
  by \code{\link{bglm_design}} when \code{drop.baseline = FALSE}: one
  coefficient per level.  \code{contr.treatment.full} is the identity
  matrix -- one indicator per level, the dummy coding of
  \code{\link[stats]{contr.treatment}} before a column is dropped.
  \code{contr.poly.full} is the \code{n} orthogonal polynomials, the
  full set that \code{\link[stats]{contr.poly}} truncates to
  \code{n - 1}.

  \code{contr.bayes.unordered} and \code{contr.bayes.ordered} are
  pinned legacy aliases (arm stores these strings in the
  \code{"contrasts"} attribute of its model frames); they are
  identical to \code{contr.treatment.full} and
  \code{contr.poly.full}.
}
\usage{
contr.treatment.full(n, base = 1, contrasts = TRUE)
contr.poly.full(n, scores = 1:n, contrasts = TRUE)

contr.bayes.unordered(n, base = 1, contrasts = TRUE)
contr.bayes.ordered(n, scores = 1:n, contrasts = TRUE)
}
\arguments{
  \item{n}{the number of levels, or the levels themselves.}
  \item{base}{the baseline level (kept for compatibility with
    \code{\link[stats]{contr.treatment}}; all \code{n} columns are
    returned either way).}
  \item{scores}{the numeric scores for the ordered case.}
  \item{contrasts}{logical; as in the \code{stats} generators.}
}
\value{
  An \code{n}-column contrast matrix with the levels as row names.
}
\examples{
contr.treatment.full(3)
contr.poly.full(3)
identical(contr.bayes.unordered(3), contr.treatment.full(3))
}
\keyword{models}
