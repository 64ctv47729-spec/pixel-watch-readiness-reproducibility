#!/usr/bin/env Rscript
# Aim 1: within-/between-person association of Google Readiness with same-day
# ordinal wellness. Primary model: cumulative-link mixed model (ordinal::clmm).
# Prespecified in ANALYSIS_DECISION_LOG.md before any fitting.

suppressPackageStartupMessages({
  library(ordinal)
  library(lme4)
})
set.seed(20260901)

root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
ds <- read.csv(file.path(root, "output", "analytic_dataset.csv"), stringsAsFactors = FALSE)
outdir <- file.path(root, "output", "models"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## ---- outcome coding (prespecified merges: extreme category n<5 in full data)
lv <- list(
  Fatigue = c("Very Tired", "Tired", "Normal", "Fresh", "Very Fresh"),
  Stress  = c("Stressed", "Normal", "Relaxed", "Very Relaxed"),      # Very Stressed (n=1) merged down
  Body    = c("Sore/Tight", "Normal", "Feeling Good", "Feeling Great"), # Very Sore/Tight (n=1) merged
  Sleep   = c("Restless Sleep", "Sleep Onset/Early Waking", "Restful", "Very Restful")
)
ds$Stress[ds$Stress == "Very Stressed"] <- "Stressed"
ds$Body[ds$Body == "Very Sore/Tight"] <- "Sore/Tight"

prep <- function(d, centering = c("episode", "participant")) {
  centering <- match.arg(centering)
  d <- d[!is.na(d$Google_Readiness), ]
  d$r10 <- d$Google_Readiness / 10
  key <- if (centering == "episode") interaction(d$Participant_ID, d$episode) else factor(d$Participant_ID)
  mu <- ave(d$r10, key, FUN = mean)
  d$r_within <- d$r10 - mu
  d$r_between <- mu - mean(mu)
  d$day_std <- as.numeric(scale(d$protocol_day))
  d$pid <- factor(d$Participant_ID)
  d
}

fit_one <- function(d, oc, label, try_slope = FALSE) {
  dd <- d[!is.na(d[[oc]]), ]
  dd$y <- factor(dd[[oc]], levels = lv[[oc]], ordered = TRUE)
  dd <- dd[!is.na(dd$y), ]
  n <- nrow(dd); npid <- length(unique(dd$pid))
  f0 <- y ~ r_within + r_between + day_std + (1 | pid)
  m <- tryCatch(clmm(f0, data = dd, link = "logit", nAGQ = 7), error = function(e) NULL)
  nagq <- 7
  if (is.null(m)) { m <- tryCatch(clmm(f0, data = dd, link = "logit"), error = function(e) NULL); nagq <- 1 }
  if (is.null(m)) return(NULL)

  slope_note <- "not attempted"
  if (try_slope) {
    ms <- tryCatch(clmm(y ~ r_within + r_between + day_std + (1 + r_within | pid),
                        data = dd, link = "logit"), error = function(e) NULL)
    if (!is.null(ms) && ms$optRes$convergence == 0) {
      vc <- VarCorr(ms)$pid
      degenerate <- any(!is.finite(vc)) || min(diag(vc)) < 1e-8 || abs(attr(vc, "correlation")[1, 2]) > 0.99
      lrt_p <- tryCatch(anova(m, ms)$`Pr(>Chisq)`[2], error = function(e) NA)
      if (!degenerate && !is.na(lrt_p) && lrt_p < 0.05) {
        m <- ms; slope_note <- sprintf("random slope retained (LRT p=%.4f)", lrt_p)
      } else slope_note <- sprintf("random slope not retained (degenerate=%s, LRT p=%s)",
                                   degenerate, format(lrt_p, digits = 3))
    } else slope_note <- "random slope model failed to converge"
  }

  co <- summary(m)$coefficients
  rows <- lapply(c("r_within", "r_between", "day_std"), function(tm) {
    if (!tm %in% rownames(co)) return(NULL)
    b <- co[tm, "Estimate"]; se <- co[tm, "Std. Error"]; p <- co[tm, "Pr(>|z|)"]
    data.frame(analysis = label, outcome = oc, term = tm, n = n, participants = npid,
               beta = b, se = se, OR = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se),
               p = p, convergence = m$optRes$convergence,
               max_grad = max(abs(m$gradient)), nAGQ = nagq,
               re_sd = sqrt(VarCorr(m)$pid[1, 1]), slope_note = slope_note,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

## proportional-odds diagnostic: threshold-specific binary GLMMs
po_diag <- function(d, oc) {
  dd <- d[!is.na(d[[oc]]), ]
  dd$y <- factor(dd[[oc]], levels = lv[[oc]], ordered = TRUE)
  dd <- dd[!is.na(dd$y), ]
  ks <- seq_len(length(lv[[oc]]) - 1)
  out <- lapply(ks, function(k) {
    dd$b <- as.integer(as.integer(dd$y) > k)
    if (min(table(dd$b)) < 8) return(data.frame(outcome = oc, threshold = k, beta = NA, se = NA, note = "split too sparse"))
    g <- tryCatch(glmer(b ~ r_within + r_between + day_std + (1 | pid), data = dd,
                        family = binomial, control = glmerControl(optimizer = "bobyqa")),
                  error = function(e) NULL)
    if (is.null(g)) return(data.frame(outcome = oc, threshold = k, beta = NA, se = NA, note = "failed"))
    s <- summary(g)$coefficients
    data.frame(outcome = oc, threshold = k, beta = s["r_within", 1], se = s["r_within", 2],
               note = ifelse(isTRUE(g@optinfo$conv$opt == 0), "ok", "check"))
  })
  do.call(rbind, out)
}

## linear mixed model on integer scores (sensitivity)
lmm_scores <- function(d, oc, label) {
  dd <- d[!is.na(d[[oc]]), ]
  dd$y <- as.numeric(factor(dd[[oc]], levels = lv[[oc]], ordered = TRUE))
  m <- lmer(y ~ r_within + r_between + day_std + (1 | pid), data = dd, REML = TRUE)
  s <- summary(m)$coefficients
  data.frame(analysis = label, outcome = oc, term = rownames(s), n = nrow(dd),
             participants = length(unique(dd$pid)), beta = s[, 1], se = s[, 2],
             lo = s[, 1] - 1.96 * s[, 2], hi = s[, 1] + 1.96 * s[, 2], stringsAsFactors = FALSE)
}

## marginal (population-averaged) ordinal model, cluster bootstrap
marginal_boot <- function(d, oc, label, B = 1000) {
  dd <- d[!is.na(d[[oc]]), ]
  dd$y <- factor(dd[[oc]], levels = lv[[oc]], ordered = TRUE)
  fitb <- function(dat) {
    m <- tryCatch(clm(y ~ r_within + r_between + day_std, data = dat, link = "logit"),
                  error = function(e) NULL)
    if (is.null(m)) return(c(NA, NA)) else coef(m)[c("r_within", "r_between")]
  }
  est <- fitb(dd)
  ids <- unique(dd$pid)
  bs <- replicate(B, {
    take <- sample(ids, length(ids), replace = TRUE)
    dat <- do.call(rbind, lapply(seq_along(take), function(i) { x <- dd[dd$pid == take[i], ]; x$pid <- paste0("b", i); x }))
    fitb(dat)
  })
  ok <- colSums(is.na(bs)) == 0
  qs <- apply(bs[, ok, drop = FALSE], 1, quantile, c(.025, .975), na.rm = TRUE)
  data.frame(analysis = label, outcome = oc, term = c("r_within", "r_between"),
             n = nrow(dd), participants = length(ids),
             beta = est, OR = exp(est), lo = exp(qs[1, ]), hi = exp(qs[2, ]),
             boot_ok = sum(ok), note = "population-averaged marginal association (combines within- and between-person variation only via its covariate decomposition; cluster-bootstrap CI)",
             stringsAsFactors = FALSE)
}

outcomes <- c("Fatigue", "Stress", "Body", "Sleep")
d_main <- prep(ds)

primary <- do.call(rbind, lapply(outcomes, function(oc) fit_one(d_main, oc, "primary CLMM (episode-centred)", try_slope = TRUE)))
## BH-FDR within families
for (tm in c("r_within", "r_between")) {
  i <- primary$term == tm & primary$analysis == "primary CLMM (episode-centred)"
  primary$q[i] <- p.adjust(primary$p[i], method = "BH")
}
write.csv(primary, file.path(outdir, "aim1_primary.csv"), row.names = FALSE)

## sensitivities
sens <- list()
d_prot <- prep(ds[ds$episode == 1 & ds$protocol_day <= 14, ])
sens$protocol <- do.call(rbind, lapply(outcomes, function(oc) fit_one(d_prot, oc, "sensitivity: episode 1, protocol days 1-14")))
dup <- read.csv(file.path(root, "output", "audit", "duplicate_questionnaire_adjudication.csv"))
dupkeys <- unique(paste(dup$Participant_ID, substr(dup$q_date, 1, 10)))
d_nodup <- prep(ds[!(paste(ds$Participant_ID, substr(ds$q_date, 1, 10)) %in% dupkeys), ])
sens$nodup <- do.call(rbind, lapply(outcomes, function(oc) fit_one(d_nodup, oc, "sensitivity: duplicate-date participant-days excluded")))
d_nozero <- prep(ds[is.na(ds$Google_Readiness) | ds$Google_Readiness > 0, ])
sens$nozero <- do.call(rbind, lapply(outcomes, function(oc) fit_one(d_nozero, oc, "sensitivity: Readiness=0 excluded")))
d_pc <- prep(ds, centering = "participant")
sens$pcentre <- do.call(rbind, lapply(outcomes, function(oc) fit_one(d_pc, oc, "sensitivity: participant-level centering")))
ds4 <- ds; ds4$Fatigue[ds4$Fatigue == "Very Tired"] <- "Tired"
lv$Fatigue <- c("Tired", "Normal", "Fresh", "Very Fresh")
sens$fat4 <- fit_one(prep(ds4), "Fatigue", "sensitivity: Fatigue 4-level (Very Tired merged)")
lv$Fatigue <- c("Very Tired", "Tired", "Normal", "Fresh", "Very Fresh")
write.csv(do.call(rbind, sens), file.path(outdir, "aim1_sensitivities.csv"), row.names = FALSE)

lmm <- do.call(rbind, lapply(outcomes, function(oc) lmm_scores(d_main, oc, "sensitivity: LMM on integer scores")))
write.csv(lmm, file.path(outdir, "aim1_lmm_scores.csv"), row.names = FALSE)

marg <- do.call(rbind, lapply(outcomes, function(oc) marginal_boot(d_main, oc, "sensitivity: marginal ordinal (cluster bootstrap)")))
write.csv(marg, file.path(outdir, "aim1_marginal.csv"), row.names = FALSE)

pod <- do.call(rbind, lapply(outcomes, function(oc) po_diag(d_main, oc)))
write.csv(pod, file.path(outdir, "aim1_po_diagnostics.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(outdir, "aim1_sessionInfo.txt"))
cat("AIM1 DONE\n")
print(primary[primary$term %in% c("r_within", "r_between"), c("outcome", "term", "n", "participants", "OR", "lo", "hi", "p", "q", "convergence", "slope_note")], digits = 3)
