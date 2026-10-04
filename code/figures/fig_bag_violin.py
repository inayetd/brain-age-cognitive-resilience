"""The corrected brain-age gap as a distribution rather than a list of numbers.

This is the figure the brain-age literature uses to show a gap, and it is the
right one: a gap is only interpretable as a position within a distribution, so
showing the distribution is showing the quantity itself. Reporting an
individual's gap as a number invites the reading that the person's brain is
that many years old, which it is not.

Only the corrected gap is drawn. The uncorrected gap is not a quantity this
study uses anywhere, and showing it here would imply a comparison that is
already made, with the correct statistics, in the bias-correction figure.

Each violin is a kernel density over all participants; every participant is
also drawn as a point, jittered horizontally, so that the tails can be counted
rather than guessed at. The white box is the interquartile range and the heavy
line the median. Participants beyond four SDs of their own modality's mean are
drawn as open rings, because they are the cases a reader will ask about.

OUTPUT  results/figures/fig9_bag_violin_{thesis.pdf,presentation.png}
"""
import json
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
RES = os.path.join(ROOT, "results")
OUTD = os.environ.get("FIG_OUT", os.path.join(RES, "figures"))

COL = {"structural": "#15497B", "functional_pca": "#7B2D8E"}
LABEL = {"structural": "Structural", "functional_pca": "Functional (PCA)"}
MAIN = ["structural", "functional_pca"]
GREY = "#5A6570"
Z_CUT = 4.0

STYLES = {
    "thesis": dict(ext="pdf", dpi=300, base=9.0, marker=4, alpha=0.26,
                   lw=1.5, scale=1.0, title_size=0, show_title=False),
    "presentation": dict(ext="png", dpi=200, base=16, marker=10, alpha=0.34,
                         lw=2.8, scale=1.55, title_size=20, show_title=True),
}


def apply_style(st):
    plt.rcParams.update({
        "figure.dpi": st["dpi"],
        "font.size": st["base"],
        "axes.labelsize": st["base"] * 1.05,
        "xtick.labelsize": st["base"] * 1.15,
        "ytick.labelsize": st["base"],
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.spines.bottom": False,
        "axes.edgecolor": "#444444",
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.04,
        "figure.facecolor": "white",
    })


def main():
    s = pd.read_csv(os.path.join(RES, "bag_sensitivity.csv"))
    met = {
        "structural": json.load(open(os.path.join(
            RES, "structural", "bag_structural_metrics.json"))),
        "functional_pca": json.load(open(os.path.join(
            RES, "functional", "bag_metrics_pca.json"))),
    }
    data = [(k, s[f"BAG_{k}_corrected"].dropna().to_numpy(float))
            for k in MAIN]

    os.makedirs(OUTD, exist_ok=True)
    for tag, st in STYLES.items():
        st = dict(st, tag=tag)
        apply_style(st)
        fig, ax = plt.subplots(figsize=(7.6 * st["scale"], 5.6 * st["scale"]))
        rng = np.random.default_rng(0)

        ax.axhline(0, color=GREY, lw=st["lw"] * 0.7, zorder=0)

        for i, (k, v) in enumerate(data):
            parts = ax.violinplot([v], positions=[i], widths=0.80,
                                  showmeans=False, showmedians=False,
                                  showextrema=False)
            for b in parts["bodies"]:
                b.set_facecolor(COL[k])
                b.set_alpha(0.26)
                b.set_edgecolor(COL[k])
                b.set_linewidth(st["lw"])

            z = (v - v.mean()) / v.std(ddof=1)
            far = np.abs(z) > Z_CUT
            jit = rng.normal(0, 0.062, len(v))
            ax.scatter(i + jit[~far], v[~far], s=st["marker"] ** 2 * 0.24,
                       color=COL[k], alpha=st["alpha"], linewidths=0, zorder=2)
            if far.any():
                ax.scatter(i + jit[far], v[far], s=st["marker"] ** 2 * 1.1,
                           facecolors="none", edgecolors=COL[k],
                           linewidths=st["lw"], zorder=4)

            q1, med, q3 = np.percentile(v, [25, 50, 75])
            ax.add_patch(plt.Rectangle((i - 0.062, q1), 0.124, q3 - q1,
                                       facecolor="white", edgecolor=COL[k],
                                       lw=st["lw"], zorder=3))
            ax.plot([i - 0.115, i + 0.115], [med, med], color=COL[k],
                    lw=st["lw"] * 1.6, zorder=5, solid_capstyle="butt")

            # The three numbers that let a reader place any individual value.
            # A real minus sign, so the annotation matches the axis ticks.
            fmt = lambda x: f"{x:+.1f}".replace("-", "−")
            ax.annotate(f"Between-person SD = {v.std(ddof=1):.2f} yr\n"
                        f"IQR = {fmt(q1)} to {fmt(q3)}\n"
                        f"r with age = "
                        f"{met[k]['bag_age_corr_corrected']:+.2f}",
                        xy=(i, 1.0), xycoords=("data", "axes fraction"),
                        xytext=(0, -6), textcoords="offset points",
                        ha="center", va="top", linespacing=1.5,
                        fontsize=st["base"] * 0.92, color=GREY)

        ax.set_xticks(range(len(data)))
        ax.set_xticklabels([LABEL[k] for k, _ in data])
        for lbl, (k, _) in zip(ax.get_xticklabels(), data):
            lbl.set_color(COL[k])
            lbl.set_weight("semibold")
        ax.tick_params(axis="x", length=0, pad=10)
        ax.set_ylabel("Corrected brain-age gap (years)")
        ax.set_xlim(-0.62, len(data) - 0.38)

        lo = min(v.min() for _, v in data)
        hi = max(v.max() for _, v in data)
        ax.set_ylim(lo - 0.06 * (hi - lo), hi + 0.20 * (hi - lo))

        ax.annotate(f"n = {len(data[0][1]):,}   ·   "
                    f"open rings: |z| > {Z_CUT:.0f} within modality",
                    xy=(0.5, -0.095), xycoords="axes fraction",
                    ha="center", va="top",
                    fontsize=st["base"] * 0.88, color=GREY)

        if st["show_title"]:
            fig.suptitle("People differ substantially in their brain-age gap",
                         fontsize=st["title_size"], y=1.0, ha="center",
                         weight="semibold")
            fig.tight_layout(rect=[0, 0, 1, 0.95])
        else:
            fig.tight_layout()
        p = os.path.join(OUTD, f"fig9_bag_violin_{tag}.{st['ext']}")
        fig.savefig(p)
        plt.close(fig)
        print("wrote " + os.path.relpath(p, ROOT))


if __name__ == "__main__":
    main()
