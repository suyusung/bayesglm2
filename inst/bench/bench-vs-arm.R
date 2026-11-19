## =====================================================================
## Speed: arm::bayesglm vs bayesglm2::bglm on identical fits.
## Run manually (it is slow); NOT part of R CMD check:
##   Rscript inst/bench/bench-vs-arm.R
## Reps configurable:  BENCH_REPS=25 Rscript inst/bench/bench-vs-arm.R
## =====================================================================
if (!requireNamespace("arm", quietly = TRUE)) {
    message("arm not installed -- benchmark skipped")
    quit(status = 0)
}

suppressMessages({library(arm); library(bayesglm2)})
set.seed(1)

timeit <- function(f, reps = as.integer(Sys.getenv("BENCH_REPS", "15"))) {
  f()                                  # warm-up
  t <- replicate(reps, {
    gc(FALSE)
    system.time(f())[["elapsed"]]
  })
  t <- t[-1]
  median(t)
}

report <- function(label, fa, fb) {
  coef_a <- coef(fa()); coef_b <- coef(fb())
  maxdiff <- max(abs(coef_a - coef_b))
  ta <- timeit(fa); tb <- timeit(fb)
  cat(sprintf("%-30s arm %8.1f ms | bayesglm2 %8.1f ms | speedup %5.2fx | max|dcoef| %.2e\n",
              label, 1000*ta, 1000*tb, ta/tb, maxdiff))
}

## ---- data -----------------------------------------------------------
n1 <- 2000
X1 <- matrix(rnorm(n1*10), n1, dimnames = list(NULL, paste0("x", 1:10)))
d1 <- data.frame(X1, g = factor(sample(1:6, n1, TRUE)))
d1$y  <- X1 %*% 1:10 + rnorm(n1)
d1$yb <- rbinom(n1, 1, plogis(X1 %*% rep(0.5, 10) - 1))
d1$yc <- rpois(n1, exp(pmin(2, X1 %*% rep(0.3, 10))))

n2 <- 8000;  d2 <- data.frame(f = factor(sample(1:50, n2, TRUE)),
                              x = rnorm(n2), z = rnorm(n2))
d2$y  <- with(d2, 2*sin(as.numeric(f)/8*pi) + 3*x + rnorm(n2))
d2$yb <- with(d2, rbinom(n2, 1, plogis(cos(as.numeric(f)/8*pi) + x)))

n3 <- 20000
X3 <- matrix(rnorm(n3*20), n3, dimnames = list(NULL, paste0("x", 1:20)))
d3 <- data.frame(X3)
d3$y  <- X3 %*% rep(0.5, 20) + rnorm(n3)
d3$yb <- rbinom(n3, 1, plogis(X3 %*% rep(0.4, 20) - 1))

f10 <- as.formula("y ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+g")
f20 <- as.formula("y ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20")
b20 <- as.formula("yb ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20")
b10 <- as.formula("yb ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10")
c10 <- as.formula("yc ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10")

cat("== gaussian ==\n")
report("gauss n=2000 p~10+fac",
       function() arm::bayesglm(f10, data = d1), function() bayesglm2::bglm(f10, data = d1))

cat("== binomial (logit) ==\n")
report("logit n=2000 p~10",
       function() arm::bayesglm(b10, data = d1), function() bayesglm2::bglm(b10, data = d1))
report("logit n=8000 fac50+x+z",
       function() arm::bayesglm(yb ~ f + x + z, data = d2), function() bayesglm2::bglm(yb ~ f + x + z, data = d2))
report("logit n=8000 fac50 only",
       function() arm::bayesglm(yb ~ f, data = d2), function() bayesglm2::bglm(yb ~ f, data = d2))

cat("== poisson ==\n")
report("poisson n=2000 p~10",
       function() arm::bayesglm(c10, data = d1), function() bayesglm2::bglm(c10, data = d1))

cat("== large n ==\n")
report("gauss n=20000 p=20",
       function() arm::bayesglm(f20, data = d3), function() bayesglm2::bglm(f20, data = d3))
report("logit n=20000 p=20",
       function() arm::bayesglm(b20, data = d3), function() bayesglm2::bglm(b20, data = d3))

cat("== flat prior (scale=Inf, df=Inf) ==\n")
report("logit n=2000 p10 flat",
       function() arm::bayesglm(b10, data = d1, prior.scale = Inf, prior.df = Inf),
       function() bayesglm2::bglm(b10, data = d1, prior.scale = Inf, prior.df = Inf))
