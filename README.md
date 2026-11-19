# bayesglm2

A from-scratch rewrite of `arm::bayesglm`: Bayesian generalized linear
models with weakly informative scaled-t priors, fit by penalized IRLS
(EM on the scale-mixture representation of the t), with **no use of
`glm()`, `glm.fit()`, or `lm.fit()`** — the linear algebra runs through
the package's own pivoted-QR C kernel (`dqrls`, the same netlib routine
behind `lm.fit`, verified bit-for-bit against it).

## Why

- **Correct final-generation working residuals.** `arm::bayesglm`
  stores working residuals mixing IRLS iterations (z from k−1 against
  eta from k); this engine stores `(y − mu)/mu'(eta)` at the final fit.
- **One coherent SE story.** `vcov()`, `summary()`, `se.coef()` and
  `display()` all use the `vcov.glm` dispersion rule.
- **Standalone.** `Imports: stats` only. The design-matrix layer
  (`model.matrixBayes`, `contr.bayes.*`) stands on its own.
- **Faster.** 1.1–1.3× `arm::bayesglm` with bit-identical coefficients
  (`inst/bench/bench-vs-arm.R`).

## Usage

```r
library(bayesglm2)
fit <- bglm(y ~ x1 + x2 + g, data = d, family = binomial)
display(fit)
se.coef(fit)
```

When `arm` is attached too, its S4 `display`/`se.coef` generics mask
these for unqualified calls; use `bayesglm2::display()` /
`bayesglm2::se.coef()`.

## Tests

Full parity suite against `arm::bayesglm` (bit-identical coefficients,
identical rank decisions, 14 acceptance tests ported from arm, and
kernel-vs-`lm.fit` bit parity) in `tests/`; run with `R CMD check`.
Design notes and pitfalls: `DESIGN.md`.
