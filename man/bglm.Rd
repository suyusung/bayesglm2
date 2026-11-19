\name{bglm}
\alias{bglm}
\alias{bglm.fit}
\title{Bayesian Generalized Linear Models (Standalone Engine)}
\description{
  Fits a generalized linear model with weakly informative scaled-t
  priors on the coefficients, by the same EM-style iteratively
  reweighted data-augmentation loop that \code{arm::bayesglm} uses.
  The engine is standalone: it never calls \code{glm} or
  \code{glm.fit}, and the linear algebra goes through this package's
  own pivoted-QR routine (which agrees with \code{lm.fit} to the last
  bit).  \code{bglm} is the formula interface; \code{bglm.fit} is the
  low-level, matrix interface.
}
\usage{
bglm(formula, family = gaussian, data, weights, subset,
     na.action, start = NULL, etastart, mustart, offset,
     control = list(...), model = TRUE, method = "bglm.fit",
     x = FALSE, y = TRUE, contrasts = NULL,
     drop.unused.levels = TRUE,
     prior.mean = 0, prior.scale = NULL, prior.df = 1,
     prior.mean.for.intercept = 0,
     prior.scale.for.intercept = NULL,
     prior.df.for.intercept = 1,
     min.prior.scale = 1e-12,
     scaled = TRUE,
     keep.order = TRUE, drop.baseline = TRUE,
     maxit = 100, print.unnormalized.log.posterior = FALSE,
     Warning = TRUE, \dots)

bglm.fit(x, y, weights = rep(1, nobs), start = NULL,
         etastart = NULL, mustart = NULL, offset = rep(0, nobs),
         family = gaussian(), control = list(), intercept = TRUE,
         prior.mean = 0, prior.scale = NULL, prior.df = 1,
         prior.mean.for.intercept = 0,
         prior.scale.for.intercept = NULL,
         prior.df.for.intercept = 1,
         min.prior.scale = 1e-12, scaled = TRUE,
         print.unnormalized.log.posterior = FALSE,
         Warning = TRUE)
}
\arguments{
  \item{formula, data, subset, weights, na.action}{model formula and
    the usual \code{\link[stats]{model.frame}} companions, handled
    exactly as in \code{\link[stats]{glm}}.}
  \item{family}{a family object or the name of one; \code{gaussian},
    \code{binomial} and \code{poisson} are exercised by the test
    suite.}
  \item{start, etastart, mustart}{starting values, as in
    \code{\link[stats]{glm}}.}
  \item{offset}{an offset, as in \code{\link[stats]{glm}}.}
  \item{control}{a list of control parameters passed through
    \code{\link[stats]{glm.control}} (\code{epsilon}, \code{maxit},
    \code{trace}).}
  \item{model}{keep the model frame on the fit.}
  \item{method}{\code{"bglm.fit"} (the only supported default name; a
    function may also be given), or \code{"model.frame"} to return the
    frame alone.}
  \item{x, y}{In \code{bglm}: logical flags saying whether to keep the
    design matrix / response on the fit.  In \code{bglm.fit}: the
    design matrix and the response.}
  \item{contrasts}{optional named list of contrast functions, passed
    to \code{\link{bglm_design}}.}
  \item{drop.unused.levels}{drop unused factor levels, as in
    \code{\link[stats]{model.frame}}.}
  \item{prior.mean, prior.scale, prior.df}{location, scale and
    degrees of freedom of the scaled-t prior on the non-intercept
    coefficients.  \code{prior.scale = NULL} uses the arm default
    (2.5, times 1.6 for probit links).  Vector values give
    coefficient-specific priors; \code{Inf} scale (with \code{Inf}
    df) is the flat limit, under which the fit is the maximum-likelihood
    one.}
  \item{prior.mean.for.intercept, prior.scale.for.intercept, prior.df.for.intercept}{the
    prior for the intercept (default scale 10, times 1.6 for probit).}
  \item{min.prior.scale}{floor for the rescaled prior scales when
    \code{scaled = TRUE}; a warning is issued when it binds.}
  \item{scaled}{rescale the prior by the spread of the predictors (and
    of the response for gaussian families), so a weight of 2 acts like
    a duplicated row.}
  \item{keep.order, drop.baseline}{passed to
    \code{\link{bglm_design}}: keep the formula's written term order
    and keep baseline columns; \code{model.matrix(fit)} and
    \code{predict} methods use them to rebuild the design consistently.}
  \item{maxit}{iteration cap (overrides \code{control$maxit} when
    given).}
  \item{print.unnormalized.log.posterior}{for binomial fits, print the
    log prior, log likelihood and their sum at each iteration.}
  \item{Warning}{suppress the engine's warnings (convergence,
    boundary, separation, rank) when \code{FALSE}.}
  \item{intercept}{(for \code{bglm.fit}) whether the first column is
    an intercept; the intercept then gets its own prior.}
  \item{\dots}{further arguments are ignored (\code{bglm}) or passed
    through; \code{predict} methods for class \code{"bglm"} follow
    \code{\link[stats]{glm}}.}
}
\value{
  For \code{bglm}, an object of class \code{"bglm"}: a list with the
  glm-shaped components (\code{coefficients}, \code{residuals}
  (working), \code{fitted.values}, \code{linear.predictors},
  \code{deviance}, \code{null.deviance}, \code{df.residual},
  \code{df.null}, \code{rank}, \code{qr}, \code{family}, \code{aic},
  \code{converged}, \code{iter}, \code{weights}, \code{y}, \code{call},
  \code{terms}, \code{formula}) plus the prior components
  (\code{prior.mean}, \code{prior.scale}, \code{prior.df},
  \code{prior.sd}, \code{dispersion}) and the design-policy components
  (\code{keep.order}, \code{drop.baseline}).  Standard S3 methods are
  provided: \code{coef}, \code{vcov}, \code{logLik}, \code{AIC},
  \code{nobs}, \code{df.residual}, \code{model.matrix},
  \code{model.frame}, \code{predict}, \code{print}, \code{summary}
  (and a \code{print} method for \code{summary}).
  \code{vcov} and \code{summary} follow \code{vcov.glm}, i.e. the
  gaussian dispersion is \code{deviance/df.residual}, not the EM value
  stored in \code{dispersion}.
  For \code{bglm.fit}, the bare list of fit components.
}
\details{
  The MAP under the scaled-t prior is found as an EM scale mixture:
  each coefficient gets one pseudo-observation row whose weight is
  \code{sqrt(dispersion)/prior.sd}, and the prior scale is updated
  each iteration from the current fit (the usual t-as-scale-mixture
  representation).  Convergence uses the old engine's criterion
  (relative deviance change with a \code{0.1 + |.|} floor, plus
  relative dispersion change for estimated-dispersion families).

  When the prior is flat enough that the augmented system is
  numerically singular (e.g. aliased columns with a flat prior), the
  fit falls back to the engine's internal maximum-likelihood IRLS
  (\code{bglm.mle.fit}), the glm.fit of this package; aliased
  coefficients are reported \code{NA}, as \code{glm} does.  Nothing
  here calls \code{\link[stats]{glm}} or \code{glm.fit}.
}
\references{
  The default priors and their scaling follow:

  Gelman, A., Jakulin, A., Pittau, M. L., and Su, Y.-S. (2008). A
  weakly informative default prior distribution for logistic and other
  regression models. \emph{Bayesian Analysis}, 3(4), 1291--1312.  The
  scaled-t (Cauchy at the default df = 1) defaults -- prior scale 2.5
  for coefficients and 10 for the intercept, on predictors scaled by
  two standard deviations -- come from this paper.

  Gelman, A. (2008). Scaling regression inputs by dividing by two
  standard deviations. \emph{Statistics in Medicine}, 27(15),
  2865--2873.  Justifies the \code{scaled = TRUE} prior scaling by
  \code{2 * sd(x)} of the predictors (and the analogous response
  scaling for the gaussian family).

  Gelman, A., and Hill, J. (2007). \emph{Data Analysis Using
  Regression and Multilevel/Hierarchical Models}.  Cambridge
  University Press.  Chapter 9 motivates weakly informative priors for
  GLMs; the engine reproduces the estimator of \code{arm::bayesglm}
  from the companion package \pkg{arm}.
}
\examples{
set.seed(1)
n <- 100
d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
d$y <- rbinom(n, 1, plogis(0.7 * d$x1 - 0.4 * d$x2))
fit <- bglm(y ~ x1 + x2, binomial, d)
summary(fit)
AIC(fit)
}
\seealso{\code{\link{bglm_design}}, \code{\link{contr.treatment.full}}}
