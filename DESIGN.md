# bayesglm2 design decisions

## No glm / glm.fit
The engine implements penalized IRLS (MAP under scaled-t priors on
standardized coefficients) itself: pivoted QR for aliases and rank
deficiency (aliased columns -> NA, NOT delegation to glm.fit),
step-halving on the Newton/Fisher updates, dispersion update for
non-binomial/poisson families. `family()` objects are used but are
standalone stats machinery, not glm.

## Naming: modern primary API + pinned legacy aliases
Primary names (this package's own API):
* bglm_design()        -- the design-matrix builder (was model.matrixBayes)
* contr.treatment.full -- full-rank dummies (was contr.bayes.unordered)
* contr.poly.full      -- full orthogonal polynomials (was contr.bayes.ordered)
* bglm()               -- the fit function (future), S3 class "bglm"
Pinned legacy names, exported as IDENTICAL aliases (reverse
dependencies call them by string): model.matrixBayes,
contr.bayes.unordered, contr.bayes.ordered. arm stores the legacy
strings in saved objects' "contrasts" attributes, so they must stay
resolvable. The parity test canonicalizes the strings when comparing
against arm.

## Fitted-object API: standard S3 names, own class
The fitted object carries S3 class "bglm". (arm's current bayesglm
returns an S3 object of class c("bayesglm","glm","lm"), NOT the old
S4 that extended glm; the separate "bglm" class keeps dispatch
unambiguous while the two packages coexist.) Porting into arm renames
it to the production class name.

Registered S3 methods: coef, vcov, logLik, AIC, nobs, df.residual,
model.matrix, model.frame, predict, print, summary,
print(summary.bglm). fitted()/residuals() deliberately use the stats
defaults: they read $fitted.values/$residuals and pad via
$na.action, exactly what we want (arm defines no residuals/fitted
methods either).

## Why contrasts() is NOT an accessor here
stats::contrasts is a plain function, not an S3 generic -- there is
nothing to dispatch on. Defining our own `contrasts` generic would
mask stats::contrasts on the search path and capture every
contrasts(factor) call in sessions that load bayesglm2 after stats.
Per-factor contrasts are exposed the way lm does it: as the
"contrasts" attribute of model.matrix(fit).

## Why model.matrixBayes() keeps its own name
* FindIt::cv.CausalANOVA calls arm:::model.matrixBayes(formula, data)
  by name -- reverse-dependency contract.
* Making the formula path an S3 method on model.matrix would override
  stats::model.matrix for formula input globally (every package would
  start receiving no-intercept matrices).
* The FindIt contract itself is pinned: the FORMULA path returns
  WITHOUT an intercept column (the attribute copy is the interface,
  see R/model.matrixBayes.R header); the TERMS path keeps it.

## Contrast generators keep their names
contr.bayes.unordered / contr.bayes.ordered are stored BY STRING in
the "contrasts" attribute and resolved by model.matrix.default's
get(), so they must stay exported under those exact names.

## Engine decisions (implemented in R/bglm.R, R/bglm-methods.R, src/bglm_qr.c)
* Linear algebra: own C routine driving netlib dqrls -- the same
  routine R's Cdqrls (the .Call behind lm.fit) drives, with the same
  argument pattern, pivot init and work size. NOT dqrdc2: the
  dqrdc2_ exported by libR.so is the influence-analysis variant and a
  mismatched hand-declared prototype silently corrupts memory.
  tests/qr-kernel-vs-lm-fit.R pins coefficient/rank equality with
  lm.fit, including R >= 4.3's NA (not 0) for aliased columns.
* Convergence keeps the old engine's exact semantics (relative
  deviance with the 0.1 + |.| floor, plus dispersion for
  estimated-dispersion families) so arm's contract tests pass
  unchanged.
* Flat-prior rank deficiency: the old engine delegated to glm.fit().
  The replacement is bglm.mle.fit, an internal plain-IRLS fitter in
  this package (must NOT be bglm.fit -- the flat prior does not
  remove the singularity, so recursing into bglm.fit loops forever).
  Aliased coefficients reported NA, bit-identical to arm and glm on
  the issue-#54 fixture.
* residuals: arm's engine stores MIXED-GENERATION working residuals
  (z from iteration k-1 against eta from k -- only cancels for the
  gaussian link). The new engine stores (y - mu)/mu'(eta) at the
  FINAL fit. This is an intentional deviation from arm, pinned by
  definition in tests/engine-parity.R rather than bug-compared.
* vcov.bglm / summary.bglm follow vcov.glm, i.e. gaussian dispersion
  = deviance/df.residual -- deliberately NOT the EM $dispersion the
  engine stores (arm's vcov/se.coef do the same through
  summary.glm). logLik.bglm mirrors stats:::logLik.glm exactly:
  rank - aic/2 with nobs/df attributes.
* Pitfall: never use object$name on the fit list where a prefix of
  name is also a field -- b$x partial-matched b$xlevels (empty for
  continuous predictors) and model.matrix.bglm silently returned
  list(). Use object[["name"]].

## Engine acceptance gate (all passing, verified 2026-11-19)
* tests/qr-kernel-vs-lm-fit.R        -- kernel vs lm.fit
* tests/model-matrix-contract.R      -- design layer + alias identity
* tests/parity-with-arm.R            -- 128-combo design parity vs arm:::
* tests/engine-parity.R              -- 133 checks: 4 prior configs x
  3 families x {coef deviance iter dispersion null.deviance aic
  linpred response se} vs arm::bayesglm, weights+offset+factor cases,
  and flat-prior bit-parity vs glm (coef/logLik/AIC/fitted/vcov/
  predict-newdata)
* tests/engine-bayesglm-*.R          -- all 14 of arm's
  tests/bayesglm-*.R acceptance tests, ported to bglm
R CMD check: Status: OK.

## Known base-R quirks the tests pin (all present in arm too)
* model.matrix.default passes contrasts = has.intercept to the
  contrast function: on the no-intercept formula path even
  contr.treatment returns the full identity (full dummies).
* model.matrix.default takes its terms from attr(data, "terms") when
  the frame has one, ignoring the passed object: keep.order
  written-order reaches the matrix only through the terms object that
  BUILT the frame (the bayesglm() path).
