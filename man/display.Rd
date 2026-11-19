\name{display}
\alias{display}
\alias{se.coef}
\alias{display.bglm}
\alias{se.coef.bglm}
\alias{fround}
\alias{pfround}
\title{Coefficient Display and Standard Errors for bglm Fits}
\description{
  \code{display} prints a coefficient table for a fitted object;
  \code{display.bglm} is the method for \code{\link{bglm}} fits.  It
  prints the call, the coefficient table (\code{coef.est} /
  \code{coef.se}, plus test statistic and p-value when
  \code{detail = TRUE}), and an \code{n}/\code{k}/deviance footer,
  following the format of \code{arm::display} for \code{bayesglm}
  objects.  The standard errors follow the package convention
  (\code{vcov}: dispersion \code{deviance/df.residual} for gaussian
  family, 1 for binomial/poisson), \emph{not} the EM dispersion the
  engine stores internally.

  \code{se.coef} returns the vector of standard errors of the
  coefficients; \code{se.coef.bglm} is \code{sqrt(diag(vcov(object)))}.

  \code{fround} rounds to \code{digits} and snaps values that round
  away to nothing to a clean zero (avoiding printed \code{-0});
  \code{pfround} prints a matrix with its columns rounded via
  \code{fround}.
}
\usage{
display(object, \dots)

\method{display}{bglm}(object, digits = 2, detail = FALSE, \dots)

se.coef(object, \dots)

\method{se.coef}{bglm}(object, \dots)

fround(x, digits = 2)

pfround(x, digits = 2)
}
\arguments{
  \item{object}{a fitted object (for the \code{bglm} methods, a fit
    from \code{\link{bglm}}).}
  \item{digits}{number of decimal places for the printed table.}
  \item{detail}{logical; if \code{TRUE} the table also carries the
    test statistic (z or t value) and p-value, and they are returned
    as \code{z.value} and \code{p.value}.}
  \item{x}{numeric vector or matrix (\code{fround}, \code{pfround}).}
  \item{\dots}{further arguments passed to \code{\link{summary}} by
    \code{display.bglm}; currently unused otherwise.}
}
\value{
  \code{display} invisibly returns a list with components \code{call},
  \code{coef}, \code{se}, \code{n}, \code{k}, \code{deviance},
  \code{null.deviance}, \code{dispersion} (and \code{z.value},
  \code{p.value} when \code{detail = TRUE}; \code{sigma.hat} for a
  gaussian fit with dispersion different from 1).  \code{se.coef}
  returns a named numeric vector.  \code{fround} returns the rounded
  numeric; \code{pfround} prints and invisibly returns it.
}
\examples{
set.seed(42)
d <- data.frame(x = rnorm(100))
d$y <- rbinom(100, 1, plogis(d$x))
fit <- bglm(y ~ x, data = d, family = binomial)
display(fit)
display(fit, detail = TRUE)
se.coef(fit)
coef(fit)
}
\seealso{\code{\link{bglm}}, \code{\link{vcov}}}
\note{
  When \code{arm} is attached as well, its S4 \code{display}/
  \code{se.coef} generics mask this package's S3 generics for
  unqualified calls, and they have no method for class \code{bglm}.
  Use \code{bayesglm2::display(fit)} / \code{bayesglm2::se.coef(fit)},
  or the methods directly (\code{display.bglm}, \code{se.coef.bglm}).
}
