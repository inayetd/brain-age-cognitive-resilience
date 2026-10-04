"""PCA against no-PCA: does the dimension reduction cost anything, and do the
two pipelines agree about individual people?

These are two separate questions and the answers differ, which is the reason
for the figure. Accuracy is a property of the model and is nearly identical
either way. The brain-age gap is a property of a person, and there the two
pipelines disagree substantially. A reader who sees only the accuracy numbers
would conclude the choice does not matter.

Pairing matters here. Both pipelines were run over the same 10 repeated
partitions with the same family-grouped folds, so fold i of repetition j is the
same held-out people in both. Comparing the fold-level values as pairs removes
the fold-to-fold variation that dominates the marginal SDs and is the only
honest way to compare two pipelines evaluated on a shared split.

The test across repetitions, not folds: the 10 repetitions are independent
re-partitions, the 100 folds are not (10 folds share a partition, and a fold's
training set overlaps every other fold's). A paired t test over 10 repetition
means respects that; one over 100 folds would not.

OUTPUT  results/figures/fig8_pca_vs_nopca_{thesis.pdf,presentation.png}
        results/table_pca_vs_nopca.csv
"""
import json
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
RES = os.path.join(ROOT, "results")
OUTD = os.environ.get("FIG_OUT", os.path.join(RES, "figures"))

PCA_COL = "#7B2D8E"      # the primary functional pipeline
NOPCA_COL = "#117A65"    # the sensitivity pipeline
GREY = "#5A6570"

# 379 parcels in the released connectivity matrices; the upper triangle is
# what enters the model when no reduction is applied.
N_EDGES = 379 * 378 // 2

STYLES = {
    "thesis": dict(ext="pdf", dpi=300, base=8.5, marker=5, alpha=0.30,
                   lw=1.4, scale=1.0, title_size=0, show_title=False),
    "presentation": dict(ext="png", dpi=200, base=15, marker=14, alpha=0.38,
                         lw=2.6, scale=1.55, title_size=19, show_title=True),
}


def apply_style(st):
    plt.rcParams.update({
        "figure.dpi": st["dpi"],
        "font.size": st["base"],
        "axes.titlesize": st["base"] * 1.05,
        "axes.labelsize": st["base"],
        "xtick.labelsize": st["base"] * 0.92,
        "ytick.labelsize": st["base"] * 0.92,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": "#444444",
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.03,
        "figure.facecolor": "white",
    })


def annotate(ax, lines, st, loc="upper left"):
    x, ha = (0.04, "left") if "left" in loc else (0.96, "right")
    y, va = (0.96, "top") if "upper" in loc else (0.04, "bottom")
    ax.text(x, y, "\n".join(lines), transform=ax.transAxes, ha=ha, va=va,
            fontsize=st["base"] * 0.92,
            bbox=dict(boxstyle="round,pad=0.35", fc="white",
                      ec="#BBBBBB", lw=0.8))


def load():
    P = json.load(open(os.path.join(RES, "functional", "bag_metrics_pca.json")))
    N = json.load(open(os.path.join(RES, "functional_nopca",
                                    "bag_metrics_nopca.json")))
    fp = pd.DataFrame(P["folds"])
    fn = pd.DataFrame(N["folds"])
    pair = fp.merge(fn, on=["seed", "fold"], suffixes=("_p", "_n"))
    if len(pair) != 100:
        raise SystemExit(f"expected 100 paired folds, got {len(pair)}")
    s = pd.read_csv(os.path.join(RES, "bag_sensitivity.csv"))
    return P, N, pair, s


def paired_test(pair, col):
    """Paired t over the 10 repetition means. See the module docstring."""
    a = pair.groupby("seed")[col + "_p"].mean()
    b = pair.groupby("seed")[col + "_n"].mean()
    t, p = stats.ttest_rel(a, b)
    d = a - b
    return d.mean(), d.std(ddof=1), t, p


def panel_paired(ax, pair, col, label, st, better_is_high):
    """Each point is one outer fold, the same held-out people in both runs."""
    x = pair[col + "_n"].to_numpy()
    y = pair[col + "_p"].to_numpy()
    lo = min(x.min(), y.min())
    hi = max(x.max(), y.max())
    pad = 0.05 * (hi - lo)
    lo, hi = lo - pad, hi + pad
    ax.plot([lo, hi], [lo, hi], ls="--", color="#BBBBBB", lw=st["lw"] * 0.8,
            zorder=1)
    ax.scatter(x, y, s=st["marker"] ** 2 * 0.5, color=PCA_COL,
               alpha=st["alpha"] + 0.2, linewidths=0, zorder=2)
    md, sd, t, p = paired_test(pair, col)
    share = (y > x).mean() if better_is_high else (y < x).mean()
    annotate(ax, [f"mean paired difference = {md:+.3f}",
                  f"PCA better in {share * 100:.0f}% of folds",
                  f"paired t(9) = {t:.2f}, p = {p:.3f}"], st,
             loc="lower right")
    ax.set_xlabel(f"{label}, no PCA")
    ax.set_ylabel(f"{label}, PCA")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")


def main():
    P, N, pair, s = load()
    a = s["BAG_functional_pca_corrected"].to_numpy(float)
    b = s["BAG_functional_nopca_corrected"].to_numpy(float)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    r = float(np.corrcoef(a, b)[0, 1])
    diff = a - b

    os.makedirs(OUTD, exist_ok=True)
    for tag, st in STYLES.items():
        st = dict(st, tag=tag)
        apply_style(st)
        fig, axes = plt.subplots(1, 3, figsize=(4.4 * st["scale"] * 3,
                                                4.5 * st["scale"]))

        panel_paired(axes[0], pair, "outer_r2", "$R^2$ per outer fold", st,
                     better_is_high=True)
        axes[0].set_title("a   Accuracy: variance explained", loc="left",
                          fontsize=st["base"] * 1.02, weight="semibold")

        panel_paired(axes[1], pair, "outer_mae", "MAE per outer fold (years)",
                     st, better_is_high=False)
        axes[1].set_title("b   Accuracy: absolute error", loc="left",
                          fontsize=st["base"] * 1.02, weight="semibold")

        ax = axes[2]
        lim = max(np.abs(np.r_[a, b])) * 1.05
        ax.plot([-lim, lim], [-lim, lim], ls="--", color="#BBBBBB",
                lw=st["lw"] * 0.8, zorder=1)
        ax.scatter(b, a, s=st["marker"] ** 2 * 0.35, color=NOPCA_COL,
                   alpha=st["alpha"], linewidths=0, zorder=2)
        annotate(ax, [f"r = {r:.3f}",
                      f"median difference = {np.median(np.abs(diff)):.1f} yr",
                      f"{100 * (np.abs(diff) > 5).mean():.0f}% differ by "
                      f"over 5 yr"], st, loc="lower right")
        ax.set_xlabel("Corrected gap, no PCA (years)")
        ax.set_ylabel("Corrected gap, PCA (years)")
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal", adjustable="box")
        ax.set_title("c   Agreement about individual people", loc="left",
                     fontsize=st["base"] * 1.02, weight="semibold")

        if st["show_title"]:
            fig.suptitle("Functional pipeline with and without PCA",
                         fontsize=st["title_size"], y=1.0, ha="center",
                         weight="semibold")
            fig.tight_layout(rect=[0, 0, 1, 0.93])
        else:
            fig.tight_layout()
        path = os.path.join(OUTD, f"fig8_pca_vs_nopca_{tag}.{st['ext']}")
        fig.savefig(path)
        plt.close(fig)
        print("wrote " + os.path.relpath(path, ROOT))

    # ---------------------------------------------------------------- table
    rows = []
    for label, m, gap, feat in [
            ("Functional (PCA)", P, a, "100 components"),
            ("Functional (no PCA)", N, b, f"{N_EDGES:,} edges")]:
        rows.append({
            "Pipeline": label,
            "Features entering the model": feat,
            "R2 (outer fold mean)": round(m["outer_r2_mean"], 3),
            "R2 SD across folds": round(
                m["outer_r2_sd_across_folds_within_seed"], 3),
            "MAE years (outer fold mean)": round(m["outer_mae_mean"], 2),
            "MAE SD across folds": round(
                m["outer_mae_sd_across_folds_within_seed"], 2),
            "gap-age r corrected": round(m["bag_age_corr_corrected"], 3),
            "BAG SD years": round(float(gap.std(ddof=1)), 2),
            "BAG within-participant SD": round(
                m["bag_sd_within_participant_mean"], 2),
        })
    tbl = pd.DataFrame(rows)
    print()
    print(tbl.to_string(index=False))

    md_r2, sd_r2, t_r2, p_r2 = paired_test(pair, "outer_r2")
    md_mae, sd_mae, t_mae, p_mae = paired_test(pair, "outer_mae")
    print(f"\npaired over 10 repetitions:")
    print(f"  dR2  = {md_r2:+.4f} (SD {sd_r2:.4f})  t(9) = {t_r2:.2f}  "
          f"p = {p_r2:.3f}")
    print(f"  dMAE = {md_mae:+.3f} (SD {sd_mae:.3f})  t(9) = {t_mae:.2f}  "
          f"p = {p_mae:.3f}")
    print(f"  agreement between the corrected gaps: r = {r:.3f}, "
          f"median difference {np.median(np.abs(diff)):.2f} yr, "
          f"{100 * (np.abs(diff) > 5).mean():.0f}% over 5 yr")

    p = os.path.join(RES, "table_pca_vs_nopca.csv")
    tbl.to_csv(p, index=False)
    print("\nwrote " + os.path.relpath(p, ROOT))


if __name__ == "__main__":
    main()
