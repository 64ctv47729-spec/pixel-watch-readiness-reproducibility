#!/usr/bin/env python3
"""Independent Python verification of the submitted Aim 1 and Aim 2 models.

The submitted analysis used R ordinal::clmm and lme4::lmer.  This script does
not reuse those fitted objects.  It evaluates the same likelihoods directly:

* cumulative-logit random-intercept models use Gauss-Hermite integration;
* Gaussian random-intercept models use profiled REML and closed-form GLS;
* participant bootstrap resampling is performed with replacement and sampled
  participants are relabelled as independent clusters.

The implementation is deliberately self-contained so that agreement is an
independent numerical check rather than a second invocation of the same code.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from numpy.polynomial.hermite import hermgauss
from scipy.optimize import minimize, minimize_scalar
from scipy.special import expit, logsumexp
from scipy.stats import norm
from statsmodels.tools.numdiff import approx_hess3


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
MOD = OUT / "models"
MOD.mkdir(parents=True, exist_ok=True)

OUTCOME_LEVELS = {
    "Fatigue": ["Very Tired", "Tired", "Normal", "Fresh", "Very Fresh"],
    "Stress": ["Stressed", "Normal", "Relaxed", "Very Relaxed"],
    "Body": ["Sore/Tight", "Normal", "Feeling Good", "Feeling Great"],
    "Sleep": ["Restless Sleep", "Sleep Onset/Early Waking", "Restful", "Very Restful"],
}
LOG_OUTCOMES = ["RMSSD", "SDNN", "SD1", "SD2"]
RAW_OUTCOMES = ["PNS_index", "SNS_index", "HR", "Resp_Rate_breaths_per_min"]


def bh(pvalues: np.ndarray) -> np.ndarray:
    p = np.asarray(pvalues, float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(ranked)
    out[order] = np.minimum(ranked, 1.0)
    return out


def prepare(d: pd.DataFrame, centering: str = "episode") -> pd.DataFrame:
    d = d.loc[d["Google_Readiness"].notna()].copy()
    d["r10"] = d["Google_Readiness"] / 10.0
    keys = ["Participant_ID", "episode"] if centering == "episode" else ["Participant_ID"]
    d["r_mean"] = d.groupby(keys, observed=True)["r10"].transform("mean")
    d["r_within"] = d["r10"] - d["r_mean"]
    d["r_between"] = d["r_mean"] - d["r_mean"].mean()
    d["day_std"] = (d["protocol_day"] - d["protocol_day"].mean()) / d["protocol_day"].std(ddof=1)
    return d


def _thresholds(raw: np.ndarray) -> np.ndarray:
    if len(raw) == 1:
        return raw.copy()
    return np.r_[raw[0], raw[0] + np.cumsum(np.exp(raw[1:]))]


def _threshold_raw(theta: np.ndarray) -> np.ndarray:
    if len(theta) == 1:
        return theta.copy()
    return np.r_[theta[0], np.log(np.maximum(np.diff(theta), 1e-3))]


def fit_clmm(y: np.ndarray, x: np.ndarray, groups: np.ndarray, nodes: int = 30,
             weights: np.ndarray | None = None) -> dict:
    """Cumulative-logit random-intercept maximum likelihood model.

    Optional non-negative case weights form a weighted pseudo-likelihood and
    are used only by the explicitly labelled IPW selection sensitivity.
    """
    y = np.asarray(y, int)
    x = np.asarray(x, float)
    groups = np.asarray(groups)
    weights = np.ones(len(y), float) if weights is None else np.asarray(weights, float)
    levels = int(y.max()) + 1
    p = x.shape[1]
    gh_x, gh_w = hermgauss(nodes)
    log_w = np.log(gh_w) - 0.5 * math.log(math.pi)
    idx = [np.flatnonzero(groups == g) for g in pd.unique(groups)]

    cum = np.array([(y <= k).mean() for k in range(levels - 1)])
    theta0 = np.log(np.clip(cum, 1e-4, 1 - 1e-4) / np.clip(1 - cum, 1e-4, 1))
    start = np.r_[_threshold_raw(theta0), np.zeros(p), math.log(0.6)]

    def nll(par: np.ndarray) -> float:
        theta = _thresholds(par[: levels - 1])
        beta = par[levels - 1 : levels - 1 + p]
        sigma = math.exp(float(np.clip(par[-1], -12, 8)))
        total = 0.0
        for ii in idx:
            eta0 = x[ii] @ beta
            b = math.sqrt(2.0) * sigma * gh_x
            eta = eta0[:, None] + b[None, :]
            upper = np.ones_like(eta)
            lower = np.zeros_like(eta)
            yi = y[ii]
            has_upper = yi < levels - 1
            has_lower = yi > 0
            if has_upper.any():
                upper[has_upper] = expit(theta[yi[has_upper], None] - eta[has_upper])
            if has_lower.any():
                lower[has_lower] = expit(theta[(yi[has_lower] - 1), None] - eta[has_lower])
            prob = np.clip(upper - lower, 1e-14, 1.0)
            total += logsumexp(log_w + (weights[ii, None] * np.log(prob)).sum(axis=0))
        return -float(total)

    bounds = [(None, None)] * (len(start) - 1) + [(-12, 5)]
    opt = minimize(nll, start, method="L-BFGS-B", bounds=bounds,
                   options={"maxiter": 3000, "ftol": 1e-12, "gtol": 1e-8, "maxls": 60})
    # A second pass with BFGS removes small box-optimizer tolerance artefacts.
    opt2 = minimize(nll, opt.x, method="BFGS",
                    options={"maxiter": 3000, "gtol": 1e-7})
    if opt2.fun <= opt.fun + 1e-7:
        par = opt2.x
        success = bool(opt.success or opt2.success)
        message = f"L-BFGS-B: {opt.message}; BFGS: {opt2.message}"
        nit = int(opt.nit + opt2.nit)
    else:
        par = opt.x
        success = bool(opt.success)
        message = str(opt.message)
        nit = int(opt.nit)

    hess = approx_hess3(par, nll, epsilon=2e-4)
    eig = np.linalg.eigvalsh(hess)
    cov = np.linalg.pinv(hess)
    beta_slice = slice(levels - 1, levels - 1 + p)
    beta = par[beta_slice]
    se = np.sqrt(np.maximum(np.diag(cov)[beta_slice], 0))
    # Finite-difference gradient, saved as a convergence diagnostic.
    eps = 1e-5
    grad = np.array([(nll(par + np.eye(1, len(par), j)[0] * eps) -
                      nll(par - np.eye(1, len(par), j)[0] * eps)) / (2 * eps)
                     for j in range(len(par))])
    return {
        "beta": beta, "se": se, "theta": _thresholds(par[: levels - 1]),
        "re_sd": math.exp(par[-1]), "logLik": -nll(par), "success": success,
        "message": message, "iterations": nit, "max_grad": float(np.max(np.abs(grad))),
        "hessian_min_eigenvalue": float(eig.min()), "nodes": nodes,
    }


def fit_random_intercept_reml(y: np.ndarray, x: np.ndarray, groups: np.ndarray) -> dict:
    """Profiled REML fit for a Gaussian random-intercept model."""
    y = np.asarray(y, float)
    x = np.asarray(x, float)
    groups = np.asarray(groups)
    uniq = pd.unique(groups)
    idx = [np.flatnonzero(groups == g) for g in uniq]
    n, p = x.shape

    def components(log_lambda: float):
        lam = math.exp(float(log_lambda))
        xtvix = np.zeros((p, p))
        xtviy = np.zeros(p)
        logdet_v0 = 0.0
        for ii in idx:
            xi = x[ii]
            yi = y[ii]
            m = len(ii)
            a = lam / (1.0 + lam * m)
            vix = xi - a * np.ones((m, 1)) @ xi.sum(axis=0, keepdims=True)
            viy = yi - a * yi.sum()
            xtvix += xi.T @ vix
            xtviy += xi.T @ viy
            logdet_v0 += math.log1p(lam * m)
        beta = np.linalg.solve(xtvix, xtviy)
        q = 0.0
        for ii in idx:
            r = y[ii] - x[ii] @ beta
            m = len(ii)
            a = lam / (1.0 + lam * m)
            q += float(r @ r - a * r.sum() ** 2)
        sigma2 = q / (n - p)
        sign, logdet_x = np.linalg.slogdet(xtvix)
        if sign <= 0 or sigma2 <= 0:
            return np.inf, None
        obj = 0.5 * (logdet_v0 + logdet_x + (n - p) * math.log(sigma2))
        return obj, (beta, xtvix, sigma2, lam)

    opt = minimize_scalar(lambda z: components(z)[0], bounds=(-16, 12), method="bounded",
                          options={"xatol": 1e-10, "maxiter": 500})
    obj, vals = components(opt.x)
    beta, xtvix, sigma2, lam = vals
    cov = sigma2 * np.linalg.inv(xtvix)
    return {
        "beta": beta, "se": np.sqrt(np.diag(cov)), "resid_sd": math.sqrt(sigma2),
        "re_sd": math.sqrt(lam * sigma2), "lambda": lam, "success": bool(opt.success),
        "objective": obj, "message": str(opt.message),
    }


def aim1(ds: pd.DataFrame) -> None:
    ds = ds.copy()
    ds.loc[ds["Stress"] == "Very Stressed", "Stress"] = "Stressed"
    ds.loc[ds["Body"] == "Very Sore/Tight", "Body"] = "Sore/Tight"

    samples = {
        "primary independent CLMM (episode-centred)": prepare(ds),
        "sensitivity: episode 1, protocol days 1-14": prepare(ds.loc[(ds.episode == 1) & (ds.protocol_day <= 14)]),
        "sensitivity: Readiness=0 excluded": prepare(ds.loc[ds.Google_Readiness.isna() | (ds.Google_Readiness > 0)]),
        "sensitivity: participant-level centering": prepare(ds, centering="participant"),
    }
    dup = pd.read_csv(OUT / "audit" / "duplicate_questionnaire_adjudication.csv")
    dup_keys = set(zip(dup.Participant_ID.astype(int), pd.to_datetime(dup.q_date).dt.date))
    keep = [((int(r.Participant_ID), pd.Timestamp(r.q_date).date()) not in dup_keys) for _, r in ds.iterrows()]
    samples["sensitivity: duplicate-date participant-days excluded"] = prepare(ds.loc[keep])

    all_rows = []
    diagnostics = []
    po_rows = []
    for label, dat in samples.items():
        for outcome, levels in OUTCOME_LEVELS.items():
            dd = dat.loc[dat[outcome].notna()].copy()
            codes = pd.Categorical(dd[outcome], categories=levels, ordered=True).codes
            ok = codes >= 0
            dd = dd.loc[ok]
            codes = codes[ok]
            x = dd[["r_within", "r_between", "day_std"]].to_numpy(float)
            fit = fit_clmm(codes, x, dd.Participant_ID.to_numpy(), nodes=30)
            for j, term in enumerate(["r_within", "r_between", "day_std"]):
                b, se = float(fit["beta"][j]), float(fit["se"][j])
                p = 2 * norm.sf(abs(b / se)) if se > 0 else np.nan
                all_rows.append({
                    "analysis": label, "outcome": outcome, "term": term,
                    "n": len(dd), "participants": dd.Participant_ID.nunique(),
                    "beta": b, "se": se, "OR": math.exp(b),
                    "lo": math.exp(b - 1.96 * se), "hi": math.exp(b + 1.96 * se),
                    "p": p, "converged": fit["success"], "re_sd": fit["re_sd"],
                    "quadrature_nodes": fit["nodes"], "max_grad": fit["max_grad"],
                })
            diagnostics.append({
                "analysis": label, "outcome": outcome, "logLik": fit["logLik"],
                "converged": fit["success"], "message": fit["message"],
                "iterations": fit["iterations"], "max_grad": fit["max_grad"],
                "hessian_min_eigenvalue": fit["hessian_min_eigenvalue"],
                "random_intercept_sd": fit["re_sd"],
                "thresholds": ";".join(f"{v:.12g}" for v in fit["theta"]),
            })
            if label.startswith("primary"):
                for cut in range(len(levels) - 1):
                    binary = (codes > cut).astype(int)
                    if np.bincount(binary, minlength=2).min() < 8:
                        po_rows.append({"outcome": outcome, "threshold": cut + 1,
                                        "beta": np.nan, "se": np.nan, "note": "split too sparse"})
                        continue
                    bf = fit_clmm(binary, x, dd.Participant_ID.to_numpy(), nodes=30)
                    po_rows.append({"outcome": outcome, "threshold": cut + 1,
                                    "beta": bf["beta"][0], "se": bf["se"][0],
                                    "note": "ok" if bf["success"] else "check"})

    res = pd.DataFrame(all_rows)
    for label in res.analysis.unique():
        for term in ["r_within", "r_between"]:
            ii = (res.analysis == label) & (res.term == term)
            res.loc[ii, "q"] = bh(res.loc[ii, "p"].to_numpy())
    res.to_csv(MOD / "aim1_independent_all.csv", index=False)
    res.loc[res.analysis.str.startswith("primary")].to_csv(MOD / "aim1_independent_primary.csv", index=False)
    pd.DataFrame(diagnostics).to_csv(MOD / "aim1_independent_diagnostics.csv", index=False)
    pd.DataFrame(po_rows).to_csv(MOD / "aim1_independent_po_diagnostics.csv", index=False)


def bootstrap_lmm(dd: pd.DataFrame, yvar: str, b: int, seed: int) -> tuple[np.ndarray, int]:
    rng = np.random.default_rng(seed)
    ids = pd.unique(dd.Participant_ID)
    blocks = {pid: dd.loc[dd.Participant_ID == pid] for pid in ids}
    vals = []
    for _ in range(b):
        take = rng.choice(ids, len(ids), replace=True)
        chunks = []
        for j, pid in enumerate(take):
            z = blocks[pid].copy()
            z["boot_group"] = j
            chunks.append(z)
        dat = pd.concat(chunks, ignore_index=True)
        x = np.column_stack([np.ones(len(dat)), dat[["r_within", "r_between", "day_std"]].to_numpy(float)])
        try:
            fit = fit_random_intercept_reml(dat[yvar].to_numpy(float), x, dat.boot_group.to_numpy())
            vals.append(fit["beta"][[1, 2]])
        except (ValueError, np.linalg.LinAlgError, FloatingPointError):
            continue
    return np.asarray(vals), len(vals)


def aim2(ds: pd.DataFrame) -> None:
    def prep2(d: pd.DataFrame) -> pd.DataFrame:
        d = d.loc[d.Google_Readiness.notna() & d.RMSSD.notna()].copy()
        d = prepare(d)
        for oc in LOG_OUTCOMES:
            d[f"log_{oc}"] = np.log(d[oc])
        return d

    samples = {
        "primary independent LMM (episode-centred), cluster bootstrap": (prep2(ds), True),
        "sensitivity: episode 1, protocol days 1-14 (same-date verified)":
            (prep2(ds.loc[(ds.episode == 1) & (ds.protocol_day <= 14)]), True),
        "sensitivity: Readiness=0 excluded":
            (prep2(ds.loc[ds.Google_Readiness.isna() | (ds.Google_Readiness > 0)]), False),
        "sensitivity: Kubios Quality >= 90":
            (prep2(ds.loc[ds.Quality.notna() & (ds.Quality >= 90)]), False),
    }
    rows = []
    diagnostics = []
    for label, (dat, do_boot) in samples.items():
        for oi, oc in enumerate(LOG_OUTCOMES + RAW_OUTCOMES, start=1):
            yvar = f"log_{oc}" if oc in LOG_OUTCOMES else oc
            dd = dat.loc[dat[yvar].notna() & np.isfinite(dat[yvar])].copy()
            x = np.column_stack([np.ones(len(dd)), dd[["r_within", "r_between", "day_std"]].to_numpy(float)])
            fit = fit_random_intercept_reml(dd[yvar].to_numpy(float), x, dd.Participant_ID.to_numpy())
            if do_boot:
                bs, n_ok = bootstrap_lmm(dd, yvar, 2000, 20261027 + oi + (0 if label.startswith("primary") else 100))
                lo, hi = np.quantile(bs, [0.025, 0.975], axis=0)
            else:
                n_ok = 0
                lo = hi = np.array([np.nan, np.nan])
            for pos, (j, term) in enumerate([(1, "r_within"), (2, "r_between")]):
                beta, se = float(fit["beta"][j]), float(fit["se"][j])
                p = 2 * norm.sf(abs(beta / se))
                is_log = oc in LOG_OUTCOMES
                rows.append({
                    "analysis": label, "outcome": oc, "term": term, "n": len(dd),
                    "participants": dd.Participant_ID.nunique(), "beta": beta, "se_wald": se,
                    "wald_lo": beta - 1.96 * se, "wald_hi": beta + 1.96 * se,
                    "boot_lo": lo[pos], "boot_hi": hi[pos], "boot_ok": n_ok,
                    "p_wald": p, "pct": 100 * (math.exp(beta) - 1) if is_log else np.nan,
                    "pct_boot_lo": 100 * (math.exp(lo[pos]) - 1) if is_log and do_boot else np.nan,
                    "pct_boot_hi": 100 * (math.exp(hi[pos]) - 1) if is_log and do_boot else np.nan,
                    "scale": "log (reported as % difference per 10 pts)" if is_log else "raw units",
                })
            diagnostics.append({
                "analysis": label, "outcome": oc, "converged": fit["success"],
                "message": fit["message"], "re_sd": fit["re_sd"],
                "resid_sd": fit["resid_sd"], "variance_ratio": fit["lambda"],
                "reml_objective": fit["objective"], "bootstrap_successes": n_ok,
            })
    res = pd.DataFrame(rows)
    for label in res.analysis.unique():
        for term in ["r_within", "r_between"]:
            ii = (res.analysis == label) & (res.term == term)
            res.loc[ii, "q_wald"] = bh(res.loc[ii, "p_wald"].to_numpy())
    res.to_csv(MOD / "aim2_independent_all.csv", index=False)
    res.loc[res.analysis.str.startswith("primary")].to_csv(MOD / "aim2_independent_primary.csv", index=False)
    pd.DataFrame(diagnostics).to_csv(MOD / "aim2_independent_diagnostics.csv", index=False)

    # Independent supplementary Kubios Readiness fit.
    dat = prep2(ds)
    dd = dat.loc[dat.Kubios_Readiness.notna()].copy()
    x = np.column_stack([np.ones(len(dd)), dd[["r_within", "r_between", "day_std"]].to_numpy(float)])
    fit = fit_random_intercept_reml(dd.Kubios_Readiness.to_numpy(float), x, dd.Participant_ID.to_numpy())
    pd.DataFrame({
        "term": ["r_within", "r_between"], "beta": fit["beta"][[1, 2]],
        "se_wald": fit["se"][[1, 2]],
        "wald_lo": fit["beta"][[1, 2]] - 1.96 * fit["se"][[1, 2]],
        "wald_hi": fit["beta"][[1, 2]] + 1.96 * fit["se"][[1, 2]],
        "n": len(dd), "participants": dd.Participant_ID.nunique(),
    }).to_csv(MOD / "aim2_independent_kubios_readiness_context.csv", index=False)


def main() -> None:
    ds = pd.read_csv(OUT / "analytic_dataset.csv", parse_dates=["q_date"])
    aim1(ds)
    aim2(ds)
    summary = {
        "implementation": "independent Python likelihood and REML verification",
        "aim1_quadrature_nodes": 30,
        "aim2_bootstrap_draws_primary": 2000,
        "aim2_bootstrap_draws_protocol_sensitivity": 2000,
        "seed_base": 20261027,
    }
    (MOD / "independent_models_meta.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

