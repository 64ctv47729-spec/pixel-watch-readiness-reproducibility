#!/usr/bin/env Rscript
# Aim 2: within-/between-person association of Google Readiness with
# independently measured Kubios physiology (Kubios Readiness excluded).
# Linear mixed models; participant-cluster bootstrap CIs (primary).

suppressPackageStartupMessages(library(lme4))
set.seed(20260902)
B <- 2000

root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
ds <- read.csv(file.path(root, "output", "analytic_dataset.csv"), stringsAsFactors = FALSE)
outdir <- file.path(root, "output", "models"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

log_outcomes <- c("RMSSD", "SDNN", "SD1", "SD2")
raw_outcomes <- c("PNS_index", "SNS_index", "HR", "Resp_Rate_breaths_per_min")
outcomes <- c(log_outcomes, raw_outcomes)

prep <- function(d) {
  d <- d[!is.na(d$Google_Readiness) & !is.na(d$RMSSD), ]
  d$r10 <- d$Google_Readiness / 10
  key <- interaction(d$Participant_ID, d$episode)
  mu <- ave(d$r10, key, FUN = mean)
  d$r_within <- d$r10 - mu
  d$r_between <- mu - mean(mu)
  d$day_std <- as.numeric(scale(d$protocol_day))
  d$pid <- factor(d$Participant_ID)
  for (oc in log_outcomes) d[[paste0("log_", oc)]] <- log(d[[oc]])
  d
}

fit_wald <- function(dd, yvar) {
  m <- lmer(as.formula(paste(yvar, "~ r_within + r_between + day_std + (1|pid)")),
            data = dd, REML = TRUE,
            control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-4)))
  s <- summary(m)$coefficients
  list(model = m, coefs = s, singular = isSingular(m))
}

analyse <- function(d, label, boot = TRUE) {
  res_list <- parallel::mclapply(seq_along(outcomes), function(i) {
    set.seed(20260902 + i)
    analyse_one(d, label, boot, outcomes[i])
  }, mc.cores = 2)
  do.call(rbind, res_list)
}

analyse_one <- function(d, label, boot, oc) {
  res <- list()
  {
    yvar <- if (oc %in% log_outcomes) paste0("log_", oc) else oc
    dd <- d[!is.na(d[[yvar]]) & is.finite(d[[yvar]]), ]
    n <- nrow(dd); npid <- length(unique(dd$pid))
    fw <- fit_wald(dd, yvar)
    est <- fw$coefs[c("r_within", "r_between"), 1]
    se <- fw$coefs[c("r_within", "r_between"), 2]

    blo <- bhi <- c(NA, NA); nb_ok <- NA
    if (boot) {
      ids <- as.character(unique(dd$pid))
      idx_by <- split(seq_len(nrow(dd)), as.character(dd$pid))
      form <- as.formula(paste(yvar, "~ r_within + r_between + day_std + (1|pid)"))
      ctrl <- lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-4),
                          calc.derivs = FALSE)
      bs <- matrix(NA_real_, 2, B)
      for (b in seq_len(B)) {
        take <- sample(ids, length(ids), replace = TRUE)
        rows <- unlist(idx_by[take], use.names = FALSE)
        dat <- dd[rows, ]
        dat$pid <- factor(rep(seq_along(take), times = lengths(idx_by[take])))
        fb <- tryCatch({
          mm <- lmer(form, data = dat, REML = TRUE, control = ctrl)
          fixef(mm)[c("r_within", "r_between")]
        }, error = function(e) NULL)
        if (!is.null(fb)) bs[, b] <- fb
      }
      ok <- colSums(is.na(bs)) == 0; nb_ok <- sum(ok)
      q <- apply(bs[, ok, drop = FALSE], 1, quantile, c(.025, .975))
      blo <- q[1, ]; bhi <- q[2, ]
    }
    for (i in 1:2) {
      tm <- c("r_within", "r_between")[i]
      is_log <- oc %in% log_outcomes
      res[[paste(oc, tm)]] <- data.frame(
        analysis = label, outcome = oc, term = tm, n = n, participants = npid,
        beta = est[i], se_wald = se[i],
        wald_lo = est[i] - 1.96 * se[i], wald_hi = est[i] + 1.96 * se[i],
        boot_lo = blo[i], boot_hi = bhi[i], boot_ok = nb_ok,
        p_wald = 2 * pnorm(-abs(est[i] / se[i])),
        pct = if (is_log) 100 * (exp(est[i]) - 1) else NA,
        pct_boot_lo = if (is_log) 100 * (exp(blo[i]) - 1) else NA,
        pct_boot_hi = if (is_log) 100 * (exp(bhi[i]) - 1) else NA,
        scale = if (is_log) "log (reported as % difference per 10 pts)" else "raw units",
        singular = fw$singular, stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, res)
}

d_main <- prep(ds)
prim <- analyse(d_main, "primary LMM (episode-centred), cluster bootstrap", boot = TRUE)
## FDR within families using Wald p (CIs remain bootstrap; reconciliation reported)
for (tm in c("r_within", "r_between")) {
  i <- prim$term == tm
  prim$q_wald[i] <- p.adjust(prim$p_wald[i], method = "BH")
}
write.csv(prim, file.path(outdir, "aim2_primary.csv"), row.names = FALSE)

## supplementary context: Kubios Readiness (NOT part of the primary family)
dd <- d_main[!is.na(d_main$Kubios_Readiness), ]
fk <- fit_wald(dd, "Kubios_Readiness")
kb <- data.frame(analysis = "supplementary context only", outcome = "Kubios_Readiness",
                 term = c("r_within", "r_between"),
                 n = nrow(dd), participants = length(unique(dd$pid)),
                 beta = fk$coefs[c("r_within", "r_between"), 1],
                 se_wald = fk$coefs[c("r_within", "r_between"), 2])
kb$wald_lo <- kb$beta - 1.96 * kb$se_wald; kb$wald_hi <- kb$beta + 1.96 * kb$se_wald
write.csv(kb, file.path(outdir, "aim2_kubios_readiness_context.csv"), row.names = FALSE)

## sensitivities
sens <- list(
  analyse(prep(ds[ds$episode == 1 & ds$protocol_day <= 14, ]),
          "sensitivity: episode 1, protocol days 1-14 (same-date verified)", boot = TRUE),
  analyse(prep(ds[is.na(ds$Google_Readiness) | ds$Google_Readiness > 0, ]),
          "sensitivity: Readiness=0 excluded", boot = FALSE),
  analyse(prep(ds[!is.na(ds$Quality) & ds$Quality >= 90, ]),
          "sensitivity: Kubios Quality >= 90", boot = FALSE)
)
write.csv(do.call(rbind, sens), file.path(outdir, "aim2_sensitivities.csv"), row.names = FALSE)

## HRV outcome correlation matrix (within-person deviations)
key <- interaction(d_main$Participant_ID, d_main$episode)
dev <- sapply(c(paste0("log_", log_outcomes), raw_outcomes), function(v)
  d_main[[v]] - ave(d_main[[v]], key, FUN = function(x) mean(x, na.rm = TRUE)))
cm <- cor(dev, use = "pairwise.complete.obs", method = "pearson")
write.csv(round(cm, 3), file.path(outdir, "aim2_outcome_correlations.csv"))

writeLines(capture.output(sessionInfo()), file.path(outdir, "aim2_sessionInfo.txt"))
cat("AIM2 DONE\n")
w <- prim[prim$term == "r_within", c("outcome", "n", "participants", "beta", "pct", "boot_lo", "boot_hi", "pct_boot_lo", "pct_boot_hi", "p_wald", "q_wald", "singular")]
print(w, digits = 3)
