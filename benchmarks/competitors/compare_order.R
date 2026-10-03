# Published procedures for a lineage effect on MIC readings, run on the
# replicates of compare_order.py. One row of results per replicate.
#
# Usage: Rscript compare_order.R readings.csv results.csv
#
# The readings are (lo, hi] on the log2 scale with -Inf and Inf at the ends of
# a panel, one row per reading, with the replicate, lineage and laboratory.
# Every procedure receives the readings of lineages with at least two readings,
# the readings the lineage share uses. Procedures that take a single value per
# reading receive the reported MIC: the upper end of the reading, and one
# dilution above the top well for a reading above it.
#
#   anova   one-way analysis of variance of the reported log2 MIC on lineage
#           (the F test of the one-way intraclass correlation), with ICC(1)
#   lmm     linear mixed model of the reported log2 MIC with a random lineage
#           intercept and, with two laboratories, a laboratory fixed effect;
#           likelihood-ratio test of the lineage variance against the mixture
#           of chi-square distributions with 0 and 1 degrees of freedom
#   cens    Gaussian regression for interval-censored log2 MIC (survreg) with
#           lineage and, with two laboratories, laboratory fixed effects;
#           likelihood-ratio test of the lineage effects
#   kw      Kruskal-Wallis test of the reported MIC on lineage
#   ricc    rank intraclass correlation of the reported MIC (rankICC, Tu et
#           al. 2023, observation weights), one-sided z test from its SE
#   ict     K-sample test of generalized Wilcoxon scores for interval-censored
#           data (interval::ictest, scores "wmw", permutational central limit
#           theorem; Fay and Shaw 2010), laboratories pooled, since the
#           package takes no strata
#   strat   the same scores computed within each laboratory, tested with
#           laboratories as blocks (coin::independence_test, quadratic
#           statistic, asymptotic permutation distribution)
.libPaths(c(Sys.getenv("R_LIBS_BENCH", "/root/Rlib"), .libPaths()))
suppressPackageStartupMessages({
  library(survival); library(interval); library(rankICC); library(lme4); library(coin)
})
args <- commandArgs(TRUE)
readings <- read.csv(args[1], stringsAsFactors = FALSE)
readings$lo <- as.numeric(readings$lo)
readings$hi <- as.numeric(readings$hi)

guard <- function(expr) tryCatch(expr, error = function(e) NA_real_)

one <- function(x) {
  x <- x[ave(seq_len(nrow(x)), x$lineage, FUN = length) >= 2, ]
  x$lin <- factor(x$lineage)
  x$labf <- factor(x$lab)
  two <- nlevels(x$labf) > 1
  y <- ifelse(is.finite(x$hi), x$hi, x$lo + 1)
  L <- 2^x$lo
  R <- 2^x$hi
  out <- list(n = nrow(x), lineages = nlevels(x$lin))

  a <- anova(lm(y ~ lin, data = x))
  ni <- as.numeric(table(x$lin)); N <- sum(ni); K <- length(ni)
  n0 <- (N - sum(ni^2) / N) / (K - 1)
  out$anova_p <- a[["Pr(>F)"]][1]
  out$anova_icc <- (a[["Mean Sq"]][1] - a[["Mean Sq"]][2]) / (a[["Mean Sq"]][1] + (n0 - 1) * a[["Mean Sq"]][2])

  lmm <- guard({
    f0 <- if (two) lm(y ~ labf, data = x) else lm(y ~ 1, data = x)
    f1 <- suppressMessages(if (two) lmer(y ~ labf + (1 | lin), data = x, REML = FALSE)
                           else lmer(y ~ 1 + (1 | lin), data = x, REML = FALSE))
    lrt <- max(0, 2 * (as.numeric(logLik(f1)) - as.numeric(logLik(f0))))
    vc <- as.data.frame(VarCorr(f1))$vcov
    c(if (lrt > 0) 0.5 * pchisq(lrt, 1, lower.tail = FALSE) else 1, vc[1] / sum(vc))
  })
  out$lmm_p <- lmm[1]; out$lmm_icc <- if (length(lmm) > 1) lmm[2] else NA_real_

  out$cens_p <- guard({
    s <- Surv(ifelse(is.finite(x$lo), x$lo, NA), ifelse(is.finite(x$hi), x$hi, NA), type = "interval2")
    X0 <- if (two) model.matrix(~ labf, x) else model.matrix(~ 1, x)
    X1 <- if (two) model.matrix(~ labf + lin, x) else model.matrix(~ lin, x)
    q <- qr(X1); X1 <- X1[, q$pivot[seq_len(q$rank)], drop = FALSE]
    m0 <- survreg(s ~ X0 - 1, dist = "gaussian")
    m1 <- survreg(s ~ X1 - 1, dist = "gaussian")
    pchisq(2 * (m1$loglik[2] - m0$loglik[2]), q$rank - ncol(X0), lower.tail = FALSE)
  })

  out$kw_p <- kruskal.test(y, x$lin)$p.value

  r <- guard(rankICC(y, x$lin, weights = "obs"))
  out$ricc <- unname(r[1]); out$ricc_low <- unname(r[3]); out$ricc_high <- unname(r[4])
  out$ricc_p <- unname(pnorm(r[1] / r[2], lower.tail = FALSE))

  out$ict_p <- guard(ictest(L, R, x$lin, scores = "wmw", method = "pclt")$p.value)

  out$strat_p <- guard({
    sc <- numeric(nrow(x))
    for (b in levels(x$labf)) {
      i <- x$labf == b
      sc[i] <- wlr_trafo(L[i], R[i], scores = "wmw")
    }
    it <- if (two) independence_test(sc ~ lin | labf, data = data.frame(sc = sc, lin = x$lin, labf = x$labf),
                                     teststat = "quadratic", distribution = "asymptotic")
          else independence_test(sc ~ lin, data = data.frame(sc = sc, lin = x$lin),
                                 teststat = "quadratic", distribution = "asymptotic")
    as.numeric(pvalue(it))
  })
  unlist(out)
}

keys <- unique(readings[, c("design", "replicate")])
rows <- lapply(seq_len(nrow(keys)), function(k) {
  x <- readings[readings$design == keys$design[k] & readings$replicate == keys$replicate[k], ]
  started <- proc.time()[["elapsed"]]
  v <- one(x)
  data.frame(design = keys$design[k], replicate = keys$replicate[k], t(v),
             seconds = proc.time()[["elapsed"]] - started)
})
write.csv(do.call(rbind, rows), args[2], row.names = FALSE)
