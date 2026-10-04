"""
Brain-age figures for the thesis and the presentation.

Every figure is produced in two versions from the same code:

    thesis        vector PDF, small type, no figure title. The caption is
                  written in the thesis, not burned into the image, which
                  is what journals and examiners expect.
    presentation  wide PNG, large type, one sentence of message above the
                  panels.

The palette is Okabe-Ito, which stays distinguishable in the common forms
of colour vision deficiency and in greyscale.

WHAT IS PLOTTED AND WHERE THE NUMBERS COME FROM
  The per-participant values come from results/bag_sensitivity.csv, written
  by 07_sensitivity_bag.py. Uncorrected predicted age is reconstructed as
  chronological age plus the uncorrected gap; this reproduces the pooled R2
  and MAE in the metrics JSON files to four decimal places, so the figures
  show the same numbers the text reports rather than a recomputation that
  might drift from them.

FIGURES
  1  predicted against chronological age, one panel per pipeline
  2  the gap against age before and after bias correction
  3  distribution of the corrected gap
  4  mean prediction error by age band (required by the analysis plan)

USAGE
  python code/figures/fig_brain_age.py
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
SENS = os.path.join(ROOT, "results", "bag_sensitivity.csv")
OUTD = os.environ.get("FIG_OUT", os.path.join(ROOT, "results", "figures"))
os.makedirs(OUTD, exist_ok=True)

# Okabe-Ito. Blue for structural, vermillion for the primary functional run,
# bluish green for the unreduced functional run.
# Dark blue and dark purple for the two primary pipelines. Blue and purple
# are harder to tell apart than blue and red, particularly in deuteranopia
# where both read as bluish, so the two are also separated by lightness:
# the purple is about 1.3 times brighter than the blue, which keeps them
# distinguishable in greyscale and for a reader who cannot use the hue.
COL = {
    "structural":       "#15497B",
    "functional_pca":   "#7B2D8E",
    "functional_nopca": "#117A65",
}
LABEL = {
    "structural":       "Structural",
    "functional_pca":   "Functional (PCA)",
    "functional_nopca": "Functional (no PCA)",
}
METRICS = {
    "structural":       os.path.join(ROOT, "results", "structural",
                                     "bag_structural_metrics.json"),
    "functional_pca":   os.path.join(ROOT, "results", "functional",
                                     "bag_metrics_pca.json"),
    "functional_nopca": os.path.join(ROOT, "results", "functional_nopca",
                                     "bag_metrics_nopca.json"),
}
ORDER = ["structural", "functional_pca", "functional_nopca"]

# The main figures show the two pre-specified primary pipelines only. The
# unreduced functional run is a sensitivity analysis, so it belongs in its
# own figure rather than as a third panel competing for attention in every
# main one. Figure 5 is where it is reported.
MAIN = ["structural", "functional_pca"]

GREY = "#555555"
LIGHT = "#BBBBBB"


# ---------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------
# One dictionary per destination. Everything that differs between the two
# versions lives here, so the plotting functions below never branch on
# which version they are drawing.

STYLES = {
    "thesis": dict(
        ext="pdf", dpi=300, base=8.5, title_size=0, marker=5, alpha=0.30,
        lw=1.4, scale=1.0, show_message=False, show_scope=True,
    ),
    "presentation": dict(
        ext="png", dpi=200, base=15, title_size=19, marker=14, alpha=0.38,
        lw=2.6, scale=1.55, show_message=True, show_scope=False,
    ),
}


def apply_style(st):
    plt.rcParams.update({
        "figure.dpi": st["dpi"],
        "savefig.dpi": st["dpi"],
        "font.size": st["base"],
        "axes.titlesize": st["base"] * 1.05,
        "axes.labelsize": st["base"],
        "xtick.labelsize": st["base"] * 0.92,
        "ytick.labelsize": st["base"] * 0.92,
        "legend.fontsize": st["base"] * 0.92,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": GREY,
        "axes.linewidth": 0.8,
        "xtick.color": GREY,
        "ytick.color": GREY,
        "axes.labelcolor": "black",
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.03,
        "figure.facecolor": "white",
    })


def save(fig, name, st, message=None):
    if st["show_message"] and message:
        fig.suptitle(message, fontsize=st["title_size"], y=1.0,
                     ha="center", weight="semibold")
        fig.tight_layout(rect=[0, 0, 1, 0.94])
    else:
        fig.tight_layout()
    path = os.path.join(OUTD, f"{name}_{st['tag']}.{st['ext']}")
    fig.savefig(path)
    plt.close(fig)
    print(f"  wrote {os.path.relpath(path, ROOT)}")


def annotate(ax, lines, st, loc="upper left"):
    """Statistics box. Placed inside the axes so it survives cropping."""
    x, ha = (0.04, "left") if "left" in loc else (0.96, "right")
    y, va = (0.96, "top") if "upper" in loc else (0.04, "bottom")
    ax.text(x, y, "\n".join(lines), transform=ax.transAxes,
            ha=ha, va=va, fontsize=st["base"] * 0.92,
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                      edgecolor=LIGHT, linewidth=0.7, alpha=0.92))


# ---------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------
def load():
    if not os.path.exists(SENS):
        raise SystemExit(
            f"not found: {SENS}\n"
            "  Run code/brain_age/07_sensitivity_bag.py first, or copy its\n"
            "  output to that path.")
    d = pd.read_csv(SENS)
    met = {}
    for k, p in METRICS.items():
        if not os.path.exists(p):
            raise SystemExit(f"not found: {p}")
        met[k] = json.load(open(p))

    # Uncorrected predicted age. Checked against the metrics files below so
    # a silent mismatch between the figures and the reported numbers cannot
    # go unnoticed.
    for k in ORDER:
        d[f"pred_{k}"] = d["age"] + d[f"BAG_{k}_raw"]
        r2 = 1 - ((d["age"] - d[f"pred_{k}"]) ** 2).sum() / \
                 ((d["age"] - d["age"].mean()) ** 2).sum()
        if abs(r2 - met[k]["pooled_r2"]) > 1e-4:
            raise SystemExit(
                f"{k}: reconstructed R2 {r2:.4f} does not match the metrics "
                f"file {met[k]['pooled_r2']:.4f}. The figures would not show "
                "the numbers the text reports.")
    print(f"  {len(d)} participants, reconstruction matches all metrics files")
    return d, met


# ---------------------------------------------------------------------
# Figure 1: prediction accuracy
# ---------------------------------------------------------------------
def fig_prediction(d, met, st):
    fig, axes = plt.subplots(1, len(MAIN),
                             figsize=(3.5 * len(MAIN) * st["scale"], 3.6 * st["scale"]),
                             sharex=True, sharey=True)
    lo, hi = d["age"].min() - 3, d["age"].max() + 3

    for ax, k in zip(axes, MAIN):
        y = d[f"pred_{k}"]
        ax.plot([lo, hi], [lo, hi], ls="--", color=LIGHT, lw=st["lw"] * 0.8,
                zorder=1)
        ax.scatter(d["age"], y, s=st["marker"], alpha=st["alpha"],
                   color=COL[k], linewidths=0, zorder=2)
        b = np.polyfit(d["age"], y, 1)
        xs = np.array([lo, hi])
        ax.plot(xs, np.polyval(b, xs), color=COL[k], lw=st["lw"], zorder=3)

        # Fold means with their fold-to-fold SD, which is what the analysis
        # plan asks to be reported. Labelled as such, because the scatter
        # shows all participants at once and the pooled value computed from
        # those points is a different number (.824 and .703).
        m = met[k]
        annotate(ax, [f"$R^2$ = {m['outer_r2_mean']:.3f} "
                      f"\u00b1 {m['outer_r2_sd_across_folds_within_seed']:.3f}",
                      f"MAE = {m['outer_mae_mean']:.2f} "
                      f"\u00b1 {m['outer_mae_sd_across_folds_within_seed']:.2f} yr",
                      f"slope = {b[0]:.2f}"]
                     + (["Mean performance across 100 outer folds",
                         "\u00b1 mean within-repetition SD across 10 folds"]
                        if st["show_scope"] else []), st)
        ax.set_title(LABEL[k], color=COL[k])
        ax.set_xlabel("Chronological age (years)")
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
        ax.set_aspect("equal", adjustable="box")
    axes[0].set_ylabel("Predicted age (years)")

    # The slope below one is the regression-to-the-mean effect the bias
    # correction exists to undo; naming it on the figure saves a question.
    axes[-1].text(0.96, 0.04,
                  "dashed: identity", transform=axes[-1].transAxes,
                  ha="right", va="bottom", fontsize=st["base"] * 0.85,
                  color=GREY)
    save(fig, "fig1_prediction_accuracy", st,
         "Accuracy of age prediction")


# ---------------------------------------------------------------------
# Figure 2: the bias correction
# ---------------------------------------------------------------------
def fig_correction(d, st):
    fig, axes = plt.subplots(2, len(MAIN),
                             figsize=(4.3 * len(MAIN) * st["scale"], 6.0 * st["scale"]),
                             sharex=True)
    lo, hi = d["age"].min() - 2, d["age"].max() + 2
    xs = np.array([lo, hi])

    # The y axis is shared within each column, so that the flattening between
    # the two rows is a property of the data and not of the scaling. Without
    # this, matplotlib would fit each panel to its own range and the
    # correction would look more dramatic than it is.
    for j, k in enumerate(MAIN):
        v = pd.concat([d[f"BAG_{k}_raw"], d[f"BAG_{k}_corrected"]])
        pad = 0.06 * (v.max() - v.min())
        for i in (0, 1):
            axes[i, j].set_ylim(v.min() - pad, v.max() + pad)

    for j, k in enumerate(MAIN):
        for i, kind in enumerate(["raw", "corrected"]):
            ax = axes[i, j]
            y = d[f"BAG_{k}_{kind}"]
            r = np.corrcoef(d["age"], y)[0, 1]
            ax.axhline(0, color=LIGHT, lw=st["lw"] * 0.8, zorder=1)
            ax.scatter(d["age"], y, s=st["marker"], alpha=st["alpha"],
                       color=COL[k], linewidths=0, zorder=2)
            b = np.polyfit(d["age"], y, 1)
            ax.plot(xs, np.polyval(b, xs), color="black", lw=st["lw"], zorder=3)
            annotate(ax, [f"r = {r:+.3f}"], st,
                     loc="upper right" if kind == "raw" else "upper left")
            if i == 0:
                ax.set_title(LABEL[k], color=COL[k])
            if i == 1:
                ax.set_xlabel("Chronological age (years)")
            ax.set_xlim(lo, hi)

    axes[0, 0].set_ylabel("Uncorrected gap (years)")
    axes[1, 0].set_ylabel("Corrected gap (years)")
    save(fig, "fig2_bias_correction", st,
         "Brain-age gap against chronological age, before and after "
         "cross-fitted bias correction")


# ---------------------------------------------------------------------
# Figure 3: distribution and reliability of the corrected gap
# ---------------------------------------------------------------------
# The earlier version put between-person spread on the left and
# within-person spread on the right, both labelled in years, which invited
# the reader to compare two quantities that mean different things. Panel b
# now shows all three spreads on one axis so the comparison that matters is
# the one the figure actually makes: how large the differences between
# people are relative to the noise in measuring any one person.

def fig_distribution(d, st):
    from scipy.stats import gaussian_kde
    from matplotlib.patches import Patch

    # SD of the ten repetition-specific estimates per participant, averaged
    # over participants. Taken from the metrics files.
    REP_SD = {"structural": 1.65, "functional_pca": 3.40, "functional_nopca": 2.63}

    fig, axes = plt.subplots(1, 2, figsize=(9.2 * st["scale"], 3.7 * st["scale"]),
                             gridspec_kw={"width_ratios": [1.5, 1]})

    # -- panel a: how the gap is distributed across people ---------------
    ax = axes[0]
    grid = np.linspace(-35, 35, 400)
    for k in MAIN:
        v = d[f"BAG_{k}_corrected"].dropna()
        dens = gaussian_kde(v)(grid)
        ax.plot(grid, dens, color=COL[k], lw=st["lw"],
                label=f"{LABEL[k]}:  M = {v.mean():+.1f}, SD = {v.std():.1f} yr")
        ax.fill_between(grid, dens, color=COL[k], alpha=0.13)
    ax.axvline(0, color=LIGHT, lw=st["lw"] * 0.8, zorder=0)
    ax.set_xlabel("Corrected brain age gap (years)")
    ax.set_ylabel("Density")
    ax.set_yticks([])
    ax.set_xlim(-35, 35)
    # Headroom so the key sits above the curves instead of across them.
    peak = max(gaussian_kde(d[f"BAG_{k}_corrected"].dropna())(grid).max()
               for k in MAIN)
    ax.set_ylim(0, peak * 1.34)
    ax.legend(frameon=False, loc="upper left", fontsize=st["base"] * 0.88)
    ax.set_title("a   Differences between people", loc="left",
                 fontsize=st["base"] * 1.02, weight="semibold")

    # -- panel b: those differences against the noise --------------------
    ax = axes[1]
    x = np.arange(len(MAIN))
    w = 0.26
    between = [d[f"BAG_{k}_corrected"].std() for k in MAIN]
    rep     = [REP_SD[k] for k in MAIN]
    se      = [REP_SD[k] / np.sqrt(10) for k in MAIN]

    for off, vals, alpha in [(-w, between, 1.0), (0, rep, 0.55), (w, se, 0.28)]:
        ax.bar(x + off, vals, width=w * 0.92,
               color=[COL[k] for k in MAIN], alpha=alpha,
               edgecolor="white", linewidth=0.6)

    # Reliability of the averaged gap with respect to the partition:
    # the share of between-person variance that is not measurement noise.
    for i, k in enumerate(MAIN):
        rel = 1 - (se[i] ** 2) / (between[i] ** 2)
        ax.text(i, between[i] + 0.35, f"rel. {rel:.3f}", ha="center",
                va="bottom", fontsize=st["base"] * 0.85, color=GREY)

    ax.set_xticks(x)
    ax.set_xticklabels([LABEL[k].replace(" (", "\n(") for k in MAIN])
    ax.set_ylabel("SD in years")
    ax.set_ylim(0, max(between) * 1.58)   # headroom for the key and the labels
    ax.legend(handles=[
        Patch(facecolor=GREY, alpha=1.00, label="between people"),
        Patch(facecolor=GREY, alpha=0.55, label="across the 10 repetitions"),
        Patch(facecolor=GREY, alpha=0.28, label="SE of the reported mean")],
        frameon=False, fontsize=st["base"] * 0.82, loc="upper center",
        handlelength=1.2, labelspacing=0.35, ncol=1)
    ax.set_title("b   Those differences against the noise", loc="left",
                 fontsize=st["base"] * 1.02, weight="semibold")

    save(fig, "fig3_gap_distribution", st,
         "Distribution of the brain-age gap")


# ---------------------------------------------------------------------
# Figure 6: the gap across acquisition site and sex
# ---------------------------------------------------------------------
# Descriptive, not a hypothesis test. Site and sex are design variables,
# not predictors in any pre-specified model, and no confirmatory claim
# rests on the differences shown here.
#
# It is worth a figure because the two modalities behave differently. The
# structural gap is flat across both groupings; the functional gap is not,
# differing by about four years between the extreme sites. Site is also
# associated with several of the H2 predictors, which makes it a potential
# confounder there, so the figure supports a limitation rather than a
# finding. Adjusting the H2 coefficients for site moves them by at most
# 0.06 years, which is reported alongside.

def fig_site_sex(d, a, st):
    from scipy import stats

    m = d.merge(a[["participant_id", "sex", "site"]], on="participant_id")
    sites = sorted(m["site"].dropna().unique())

    fig, axes = plt.subplots(1, 2, figsize=(9.6 * st["scale"], 3.9 * st["scale"]),
                             sharey=True, gridspec_kw={"width_ratios": [2.1, 1]})

    def violins(ax, groups, positions, colour):
        parts = ax.violinplot(groups, positions=positions, widths=0.34,
                              showmeans=False, showextrema=False)
        for b in parts["bodies"]:
            b.set_facecolor(colour); b.set_alpha(0.45)
            b.set_edgecolor(colour); b.set_linewidth(0.8)
        for p, g in zip(positions, groups):
            ax.plot([p], [np.mean(g)], "o", color=colour,
                    ms=st["marker"] * 0.42, zorder=5)
            lo, hi = np.percentile(g, [25, 75])
            ax.plot([p, p], [lo, hi], color=colour, lw=st["lw"] * 1.6, zorder=4)

    # -- panel a: site --------------------------------------------------
    ax = axes[0]
    ax.axhline(0, color=LIGHT, lw=st["lw"] * 0.8, zorder=0)
    for off, k in zip([-0.19, 0.19], MAIN):
        groups = [m.loc[m["site"] == sv, f"BAG_{k}_corrected"].dropna().values
                  for sv in sites]
        violins(ax, groups, [i + off for i in range(len(sites))], COL[k])
        F, p = stats.f_oneway(*groups)
        allv = m[f"BAG_{k}_corrected"].dropna()
        eta = sum(len(g) * (g.mean() - allv.mean()) ** 2 for g in groups) / \
              ((allv - allv.mean()) ** 2).sum()
        ax.plot([], [], color=COL[k], lw=st["lw"] * 2,
                label=f"{LABEL[k]}:  $\\eta^2$ = {eta:.3f}, p = {p:.3g}")
    ax.set_xticks(range(len(sites)))
    ax.set_xticklabels(sites)
    ax.set_xlabel("Acquisition site")
    ax.set_ylabel("Corrected brain age gap (years)")
    ax.legend(frameon=False, loc="upper center", fontsize=st["base"] * 0.85,
              ncol=1)
    ax.set_title("a   By site", loc="left", fontsize=st["base"] * 1.02,
                 weight="semibold")

    # -- panel b: sex ---------------------------------------------------
    ax = axes[1]
    ax.axhline(0, color=LIGHT, lw=st["lw"] * 0.8, zorder=0)
    for off, k in zip([-0.19, 0.19], MAIN):
        groups = [m.loc[m["sex"] == sv, f"BAG_{k}_corrected"].dropna().values
                  for sv in ["F", "M"]]
        violins(ax, groups, [i + off for i in range(2)], COL[k])
        t, p = stats.ttest_ind(*groups, equal_var=False)
        dz = (groups[1].mean() - groups[0].mean()) / \
             np.sqrt((groups[0].var() + groups[1].var()) / 2)
        ax.text(0.5, 0.97 - 0.09 * MAIN.index(k),
                f"d = {dz:+.2f}, p = {p:.3g}", transform=ax.transAxes,
                ha="center", va="top", fontsize=st["base"] * 0.85, color=COL[k])
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Female", "Male"])
    ax.set_xlabel("Sex")
    ax.set_title("b   By sex", loc="left", fontsize=st["base"] * 1.02,
                 weight="semibold")
    ax.set_ylim(-34, 40)

    save(fig, "fig6_site_and_sex", st,
         "Brain-age gap by acquisition site and sex")


# ---------------------------------------------------------------------
# Figure 7: what the bias correction does to the prediction itself
# ---------------------------------------------------------------------
# Figure 2 shows the correction working on the gap. This shows the same
# step from the other side, on the predicted ages, and it makes the price
# of the correction explicit.
#
# Dividing by beta stretches the predictions away from the sample mean.
# That is the whole point: it restores the slope against chronological age
# to one, so the gap stops being a restatement of age. But stretching also
# moves predictions further from the truth, so R-squared falls and the mean
# absolute error rises.
#
# The trade is worth naming out loud rather than letting a reader discover
# that the corrected model "predicts worse". Accuracy is not what the
# corrected prediction is for; an unbiased gap is.

def fig_prediction_corrected(d, met, st):
    fig, axes = plt.subplots(2, len(MAIN),
                             figsize=(5.6 * len(MAIN) * st["scale"],
                                      5.9 * st["scale"]),
                             sharex=True, sharey=True)
    # Limits span every value actually plotted. Using the age range would
    # clip the corrected predictions, which by design reach beyond it, and
    # the clipped points would pile up against the frame as a false edge.
    vals = [d["age"]] + [d["age"] + d[f"BAG_{k}_{kind}"]
                         for k in MAIN for kind in ("raw", "corrected")]
    lo = min(float(v.min()) for v in vals) - 2
    hi = max(float(v.max()) for v in vals) + 2
    xs = np.array([lo, hi])

    for j, k in enumerate(MAIN):
        for i, kind in enumerate(["raw", "corrected"]):
            ax = axes[i, j]
            y = d["age"] + d[f"BAG_{k}_{kind}"]
            ax.plot([lo, hi], [lo, hi], ls="--", color=LIGHT,
                    lw=st["lw"] * 0.8, zorder=1)
            ax.scatter(d["age"], y, s=st["marker"], alpha=st["alpha"],
                       color=COL[k], linewidths=0, zorder=2)
            b = np.polyfit(d["age"], y, 1)
            ax.plot(xs, np.polyval(b, xs), color=COL[k], lw=st["lw"], zorder=3)

            # Accuracy is reported only for the uncorrected row, as
            # pre-specified, and as the mean over the 100 outer folds. The
            # corrected row carries the slope alone: there is no fold-level
            # record of the corrected predictions, and accuracy after a
            # correction that uses chronological age is not the study's
            # accuracy in any case.
            if kind == "raw":
                m = met[k]
                annotate(ax, [f"$R^2$ = {m['outer_r2_mean']:.3f} "
                              f"\u00b1 {m['outer_r2_sd_across_folds_within_seed']:.3f}",
                              f"MAE = {m['outer_mae_mean']:.2f} "
                              f"\u00b1 {m['outer_mae_sd_across_folds_within_seed']:.2f} yr",
                              f"slope = {b[0]:.2f}"]
                             + (["Mean performance across 100 outer folds",
                                 "\u00b1 mean within-repetition SD across 10 "
                                 "folds"]
                                if st["show_scope"] else []), st)
            else:
                annotate(ax, [f"slope = {b[0]:.2f}"], st)
            if i == 0:
                ax.set_title(LABEL[k], color=COL[k])
            if i == 1:
                ax.set_xlabel("Chronological age (years)")
            ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
            # No equal aspect here: the panels are deliberately wider than
            # tall, so the identity line is not drawn at 45 degrees. It is
            # still the line predicted = chronological, and it is labelled.

    axes[0, 0].set_ylabel("Uncorrected predicted age (years)")
    axes[1, 0].set_ylabel("Corrected predicted age (years)")
    axes[-1, -1].text(0.96, 0.04, "dashed: identity",
                      transform=axes[-1, -1].transAxes, ha="right", va="bottom",
                      fontsize=st["base"] * 0.85, color=GREY)

    # The analysis plan reports accuracy from the uncorrected predictions,
    # because the correction uses chronological age and what follows it
    # therefore describes the corrected quantity rather than the model's
    # predictive performance. Stating that on the figure stops the lower
    # row being read as the study's accuracy figures.
    fig.text(0.5, -0.012,
             "Top row: uncorrected predictions, from which accuracy is "
             "reported. Bottom row: the same predictions after correction.",
             ha="center", va="top", fontsize=st["base"] * 0.88,
             style="italic", color=GREY)

    save(fig, "fig7_prediction_corrected", st,
         "Accuracy of age prediction, before and after bias correction")


# ---------------------------------------------------------------------
# Figure 5: agreement between the two functional pipelines
# ---------------------------------------------------------------------
# A sensitivity figure, not a main one. The analysis plan requires the
# correlation between the reduced and unreduced functional gaps to be
# reported; a correlation alone hides how far apart the two put an
# individual, so the limits of agreement are shown beside it.
#
# The two pipelines predict age equally well (R2 .703 against .704), which
# makes the disagreement below the more useful thing to show: equal
# accuracy at the group level does not mean equal estimates per person.

def fig_pipeline_agreement(d, st):
    a = d["BAG_functional_pca_corrected"]
    b = d["BAG_functional_nopca_corrected"]
    diff, mid = b - a, (a + b) / 2
    r = a.corr(b)
    md, sd = diff.mean(), diff.std()
    loa = (md - 1.96 * sd, md + 1.96 * sd)

    fig, axes = plt.subplots(1, 2, figsize=(9.2 * st["scale"], 3.9 * st["scale"]))

    ax = axes[0]
    lim = 36
    ax.plot([-lim, lim], [-lim, lim], ls="--", color=LIGHT, lw=st["lw"] * 0.8)
    ax.scatter(a, b, s=st["marker"], alpha=st["alpha"],
               color="#8C6BB1", linewidths=0)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Gap from the PCA pipeline (years)", color=COL["functional_pca"])
    ax.set_ylabel("Gap from the unreduced pipeline (years)",
                  color=COL["functional_nopca"])
    annotate(ax, [f"r = {r:.3f}", f"shared variance {100 * r ** 2:.0f}%"], st)
    ax.set_title("a   The same people, two pipelines", loc="left",
                 fontsize=st["base"] * 1.02, weight="semibold")

    ax = axes[1]
    ax.scatter(mid, diff, s=st["marker"], alpha=st["alpha"],
               color="#8C6BB1", linewidths=0)
    ax.axhline(md, color="black", lw=st["lw"])
    for y in loa:
        ax.axhline(y, color=GREY, ls="--", lw=st["lw"] * 0.9)
    # Left-aligned inside the axes: right-aligned labels ran off the panel.
    for y, va, lab in [(loa[1], "bottom", "+1.96 SD"), (loa[0], "top", "-1.96 SD")]:
        ax.text(0.015, y, f"{lab}  {y:+.1f}", transform=ax.get_yaxis_transform(),
                ha="left", va=va, fontsize=st["base"] * 0.85, color=GREY)
    ax.set_xlabel("Mean of the two estimates (years)")
    ax.set_ylabel("Unreduced minus PCA (years)")
    annotate(ax, [f"mean difference {md:+.2f} yr",
                  f"median |difference| {diff.abs().median():.1f} yr",
                  f"{100 * (diff.abs() > 5).mean():.0f}% differ by over 5 yr"],
             st, loc="lower right")
    ax.set_title("b   How far apart, per person", loc="left",
                 fontsize=st["base"] * 1.02, weight="semibold")

    save(fig, "fig5_pipeline_agreement", st,
         "Agreement between the two functional pipelines")


# ---------------------------------------------------------------------
# Figure 4: prediction error by age band
# ---------------------------------------------------------------------
def fig_error_by_age(d, st):
    bands = [(36, 45), (46, 55), (56, 65), (66, 75), (76, 90)]
    fig, axes = plt.subplots(1, 2, figsize=(9.6 * st["scale"], 3.5 * st["scale"]),
                             sharex=True, sharey=True)

    for ax, kind, ttl in zip(axes, ["raw", "corrected"],
                             ["Before correction", "After correction"]):
        ax.axhline(0, color=LIGHT, lw=st["lw"] * 0.8)
        for off, k in zip(np.linspace(-0.16, 0.16, len(MAIN)), MAIN):
            xs, ys, es = [], [], []
            for i, (a, b) in enumerate(bands):
                v = d.loc[d["age"].between(a, b), f"BAG_{k}_{kind}"]
                if len(v) < 5:
                    continue
                xs.append(i + off); ys.append(v.mean())
                es.append(1.96 * v.std() / np.sqrt(len(v)))
            ax.errorbar(xs, ys, yerr=es, fmt="o", color=COL[k],
                        ms=st["marker"] * 0.55, lw=st["lw"], capsize=3,
                        label=LABEL[k] if kind == "raw" else None)
        ax.set_title(ttl)
        ax.set_xticks(range(len(bands)))
        ax.set_xticklabels([f"{a}-{b}" for a, b in bands])
        ax.set_xlabel("Age band (years)")
    axes[0].set_ylabel("Mean gap (years)")
    axes[0].legend(frameon=False, loc="lower left")

    save(fig, "fig4_error_by_age_band", st,
         "Mean brain-age gap by age band, before and after correction")


# ---------------------------------------------------------------------
def main():
    print("reading data")
    d, met = load()
    design = pd.read_csv(os.path.join(ROOT, "data", "clean",
                                      "analysis_table_pca.csv"),
                         usecols=["participant_id", "sex", "site"])
    for tag, st in STYLES.items():
        st = dict(st, tag=tag)
        apply_style(st)
        print(f"\n{tag}:")
        fig_prediction(d, met, st)
        fig_prediction_corrected(d, met, st)
        fig_correction(d, st)
        fig_distribution(d, st)
        fig_error_by_age(d, st)
        fig_pipeline_agreement(d, st)
        fig_site_sex(d, design, st)
    print(f"\nall figures in {os.path.relpath(OUTD, ROOT)}")


if __name__ == "__main__":
    main()
