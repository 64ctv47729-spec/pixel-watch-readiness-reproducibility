#!/usr/bin/env python3
"""Generate corrected/independently verified tables and figures."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
TAB = ROOT / "tables"
FIG = ROOT / "figures"
TAB.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)

BLUE = "#2364AA"
ORANGE = "#D97706"
GREY = "#6B7280"


def save(fig, stem: str) -> None:
    fig.savefig(FIG / f"{stem}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / f"{stem}.tiff", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def tables() -> None:
    audit = json.loads((OUT / "audit" / "audit_summary.json").read_text())
    ds = pd.read_csv(OUT / "analytic_dataset.csv")
    wellness = [int(ds[c].notna().sum()) for c in ["Fatigue", "Stress", "Body", "Sleep"]]
    t1 = pd.DataFrame([
        ["Participant identifiers", audit["identifiers"]],
        ["Distinct consented participants", "not independently verifiable without enrolment log"],
        ["Demographics", "not independently reproducible for all 38 identifiers"],
        ["Monitoring episodes", audit["episodes"]],
        ["Participant-day records", audit["rows"]],
        ["Wellness responses present", f"{min(wellness)}–{max(wellness)}"],
        ["Daily Readiness recorded", f"{audit['readiness']['n']} ({audit['readiness']['n']/audit['rows']*100:.1f}%)"],
        ["Daily Readiness, mean (SD)", f"{audit['readiness']['mean']:.1f} ({audit['readiness']['sd']:.1f})"],
        ["Kubios records", audit["kubios_timing"]["n"]],
        ["Kubios 04:00–11:59", f"{audit['kubios_timing']['morning_0400_1159']} ({audit['kubios_timing']['pct_morning']:.1f}%)"],
        ["Questionnaire after Kubios", f"{audit['kubios_timing']['questionnaire_after_kubios']} ({audit['kubios_timing']['questionnaire_after_kubios']/audit['kubios_timing']['n']*100:.1f}%)"],
        ["Kubios-to-questionnaire lag, median (IQR), min",
         f"{audit['kubios_timing']['median_lag_min']:.1f} ({audit['kubios_timing']['q1_lag']:.1f}–{audit['kubios_timing']['q3_lag']:.1f})"],
    ], columns=["Characteristic", "Audited value"])
    t1.to_csv(TAB / "Table1_audited_sample.csv", index=False)

    a1 = pd.read_csv(OUT / "models" / "aim1_independent_primary.csv")
    rows = []
    for oc in ["Fatigue", "Stress", "Body", "Sleep"]:
        w = a1[(a1.outcome == oc) & (a1.term == "r_within")].iloc[0]
        b = a1[(a1.outcome == oc) & (a1.term == "r_between")].iloc[0]
        rows.append([oc, int(w.n), int(w.participants), w.OR, w.lo, w.hi, w.p, w.q,
                     b.OR, b.lo, b.hi, b.p, b.q])
    pd.DataFrame(rows, columns=["Outcome", "n", "participants", "within_OR", "within_lo",
                                "within_hi", "within_p", "within_q", "between_episode_OR",
                                "between_episode_lo", "between_episode_hi", "between_episode_p",
                                "between_episode_q"]).to_csv(TAB / "Table2_aim1_independent.csv", index=False)

    a2 = pd.read_csv(OUT / "models" / "aim2_independent_primary.csv")
    rows = []
    for oc in ["RMSSD", "SDNN", "SD1", "SD2", "PNS_index", "SNS_index", "HR", "Resp_Rate_breaths_per_min"]:
        w = a2[(a2.outcome == oc) & (a2.term == "r_within")].iloc[0]
        b = a2[(a2.outcome == oc) & (a2.term == "r_between")].iloc[0]
        log_scale = oc in {"RMSSD", "SDNN", "SD1", "SD2"}
        rows.append([oc, int(w.n), int(w.participants), w.pct if log_scale else w.beta,
                     w.pct_boot_lo if log_scale else w.boot_lo,
                     w.pct_boot_hi if log_scale else w.boot_hi, w.p_wald, w.q_wald,
                     b.pct if log_scale else b.beta,
                     b.pct_boot_lo if log_scale else b.boot_lo,
                     b.pct_boot_hi if log_scale else b.boot_hi, b.p_wald, b.q_wald,
                     "percent" if log_scale else "native units"])
    pd.DataFrame(rows, columns=["Outcome", "n", "participants", "within_effect", "within_boot_lo",
                                "within_boot_hi", "within_p", "within_q", "between_episode_effect",
                                "between_episode_boot_lo", "between_episode_boot_hi", "between_episode_p",
                                "between_episode_q", "scale"]).to_csv(TAB / "Table3_aim2_independent.csv", index=False)

    a3 = pd.read_csv(ROOT / "results" / "corrected" / "aim3" / "aim3_grouped_cv_corrected.csv")
    a3.to_csv(TAB / "Table4_aim3_corrected.csv", index=False)


def forest_aim1() -> None:
    d = pd.read_csv(OUT / "models" / "aim1_independent_primary.csv")
    d = d[d.term == "r_within"].set_index("outcome").loc[["Fatigue", "Stress", "Body", "Sleep"]]
    y = np.arange(len(d))[::-1]
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.axvline(1, color="black", lw=1, ls="--")
    ax.errorbar(d.OR, y, xerr=[d.OR - d.lo, d.hi - d.OR], fmt="o", color=BLUE,
                ecolor=BLUE, capsize=4, lw=1.6)
    ax.set_yticks(y, ["Fatigue/freshness", "Stress/relaxation", "Body soreness/recovery", "Sleep quality"])
    ax.set_xlabel("Cumulative odds ratio per 10-point higher-than-usual Readiness")
    ax.set_title("Independent Aim 1 reanalysis")
    ax.grid(axis="x", alpha=.2)
    save(fig, "Fig3_aim1_independent")


def forest_aim1_manuscript() -> None:
    d = pd.read_csv(OUT / "models" / "aim1_independent_primary.csv")
    outcomes = ["Fatigue", "Stress", "Body", "Sleep"]
    labels = ["Fatigue/freshness", "Stress/relaxation", "Body soreness/recovery", "Sleep quality"]
    y = np.arange(len(outcomes))[::-1]
    fig, ax = plt.subplots(figsize=(5.6, 3.0))
    for term, colour, offset, legend in [
        ("r_within", BLUE, 0.16, "Within-person"),
        ("r_between", ORANGE, -0.16, "Episode-mean"),
    ]:
        part = d[d.term == term].set_index("outcome").loc[outcomes]
        for i, outcome in enumerate(outcomes):
            row = part.loc[outcome]
            faded = 0.35 if (term == "r_within" and row.q >= 0.05) or term == "r_between" else 1.0
            ax.errorbar(
                row.OR,
                y[i] + offset,
                xerr=[[row.OR - row.lo], [row.hi - row.OR]],
                fmt="o",
                color=colour,
                ecolor=colour,
                alpha=max(faded, 0.45),
                capsize=2.5,
                markersize=5,
                lw=1.4,
                label=legend if i == 0 else None,
            )
    ax.axvline(1, color=GREY, lw=0.8, ls="--")
    ax.set_xscale("log")
    from matplotlib.ticker import FixedLocator, NullFormatter

    ax.xaxis.set_major_locator(FixedLocator([0.7, 1.0, 1.5, 2.0, 3.0]))
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda value, _: f"{value:g}"))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_yticks(y, labels)
    ax.set_xlabel("Cumulative odds ratio per 10 Readiness points (95% CI)")
    ax.legend(frameon=False, loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=8)
    fig.tight_layout()
    save(fig, "Fig3_aim1_manuscript")


def forest_aim2() -> None:
    d = pd.read_csv(OUT / "models" / "aim2_independent_primary.csv")
    d = d[d.term == "r_within"].set_index("outcome")
    left_names = ["RMSSD", "SDNN", "SD1", "SD2"]
    right_names = ["PNS_index", "SNS_index", "HR", "Resp_Rate_breaths_per_min"]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6))
    for ax, names, pct, title in [(axes[0], left_names, True, "Log-HRV outcomes (%)"),
                                  (axes[1], right_names, False, "Native-unit outcomes")]:
        z = d.loc[names]
        est = z.pct if pct else z.beta
        lo = z.pct_boot_lo if pct else z.boot_lo
        hi = z.pct_boot_hi if pct else z.boot_hi
        y = np.arange(len(names))[::-1]
        ax.axvline(0, color="black", lw=1, ls="--")
        ax.errorbar(est, y, xerr=[est - lo, hi - est], fmt="o", color=BLUE,
                    ecolor=BLUE, capsize=4, lw=1.6)
        labels = ["PNS index" if x == "PNS_index" else "SNS index" if x == "SNS_index"
                  else "Heart rate" if x == "HR" else "Respiratory rate" if x.startswith("Resp") else x
                  for x in names]
        ax.set_yticks(y, labels)
        ax.set_title(title)
        ax.grid(axis="x", alpha=.2)
    axes[0].set_xlabel("Percent difference per 10 Readiness points")
    axes[1].set_xlabel("Difference per 10 Readiness points")
    fig.suptitle("Independent Aim 2 reanalysis (2,000-cluster-bootstrap CIs)")
    fig.tight_layout()
    save(fig, "Fig1_aim2_independent")


def scatter() -> None:
    d = pd.read_csv(OUT / "analytic_dataset.csv")
    d = d.dropna(subset=["Google_Readiness", "RMSSD", "HR"]).copy()
    d["logRMSSD"] = np.log(d.RMSSD)
    key = [d.Participant_ID, d.episode]
    for col in ["Google_Readiness", "logRMSSD", "HR"]:
        d[col + "_dev"] = d[col] - d.groupby(["Participant_ID", "episode"])[col].transform("mean")
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.5))
    for ax, ycol, ylabel in [(axes[0], "logRMSSD_dev", "Within-episode log RMSSD deviation"),
                             (axes[1], "HR_dev", "Within-episode heart-rate deviation")]:
        x = d.Google_Readiness_dev.to_numpy(); y = d[ycol].to_numpy()
        ax.scatter(x, y, s=17, alpha=.45, color=BLUE, edgecolor="none")
        coef = np.polyfit(x, y, 1); xx = np.linspace(x.min(), x.max(), 100)
        # Orientation-only OLS line with a participant-cluster bootstrap band.
        # Resampling participants preserves each participant's repeated rows.
        participant_ids = d.Participant_ID.unique()
        rng = np.random.default_rng(20261027 + (0 if ycol == "logRMSSD_dev" else 1))
        bootstrap_lines = np.empty((2000, len(xx)))
        grouped = {pid: d[d.Participant_ID.eq(pid)] for pid in participant_ids}
        for draw in range(2000):
            selected = rng.choice(participant_ids, len(participant_ids), replace=True)
            xb = np.concatenate([grouped[pid].Google_Readiness_dev.to_numpy() for pid in selected])
            yb = np.concatenate([grouped[pid][ycol].to_numpy() for pid in selected])
            bootstrap_lines[draw] = np.polyval(np.polyfit(xb, yb, 1), xx)
        lo, hi = np.quantile(bootstrap_lines, [.025, .975], axis=0)
        ax.fill_between(xx, lo, hi, color=ORANGE, alpha=.18, linewidth=0)
        ax.plot(xx, np.polyval(coef, xx), color=ORANGE, lw=2)
        ax.axhline(0, color=GREY, lw=.7); ax.axvline(0, color=GREY, lw=.7)
        ax.set_xlabel("Within-episode Readiness deviation")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=.15)
    fig.suptitle("Within-episode physiological relationships")
    fig.tight_layout()
    save(fig, "Fig2_within_scatter_regenerated")


def forest_aim3() -> None:
    d = pd.read_csv(ROOT / "results" / "corrected" / "aim3" / "aim3_grouped_cv_corrected.csv")
    d = d.set_index("outcome").loc[["Fatigue", "Stress", "Body", "Sleep"]]
    y = np.arange(len(d))[::-1]
    fig, ax = plt.subplots(figsize=(7.4, 4.5))
    ax.axvline(0, color="black", lw=1, ls="--")
    ax.errorbar(d.rps_diff, y, xerr=[d.rps_diff - d.diff_boot_lo,
                                     d.diff_boot_hi - d.rps_diff],
                fmt="o", color=ORANGE, ecolor=ORANGE, capsize=4, lw=1.6)
    ax.set_yticks(y, ["Fatigue/freshness", "Stress/relaxation", "Body soreness/recovery", "Sleep quality"])
    ax.set_xlabel("Expanded minus base row-weighted RPS (negative favours Readiness)")
    ax.grid(axis="x", alpha=.2)
    save(fig, "Fig4_aim3_corrected")


def main() -> None:
    tables()
    forest_aim1()
    forest_aim1_manuscript()
    forest_aim2()
    scatter()
    forest_aim3()
    print("corrected tables and figures written")


if __name__ == "__main__":
    main()
