"""H2: lifestyle and psychosocial correlates of the brain-age gap.

The figure has to carry a null, which is harder than carrying a result.
Eleven of the twelve tests do not survive correction, and the honest job
of the plot is to make that legible rather than to let the eye hunt for
the one interval that clears zero.

So: every test is drawn, in the order the hypothesis states them, with
the survivors filled and the rest hollow. Nothing is sorted by p value -
ranking a confirmatory family by its own results turns a test into a
league table.

The covariates are not here. Age, sex and education are in every model as
adjusters, not results.

SCALE
  The gap is the outcome in H2 and stays in years, so the coefficient
  reads as years of brain-age gap per standard deviation of the exposure.
  This is NOT comparable to the standardised betas of H1 and H3, and the
  axis label says what it is.

WHERE THE q VALUES COME FROM
  08_multiplicity.R, when it has been run for this tag. If it has not,
  the script computes Benjamini-Hochberg itself and says so in the
  console, because a figure that silently invents its own correction is
  worse than one that refuses to draw.

INPUT   results/sem/h2_single_estimates_<TAG>.csv
        results/sem/multiplicity_<TAG>.csv        (preferred, for q)
        results/sem/h2_single_variance_<TAG>.csv  (for delta R2)
OUTPUT  results/figures/fig12_h2_{thesis.pdf,presentation.png}
"""
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
TAG = os.environ.get("FUNCTIONAL_TAG", "pca_it50")

COL = {"structural": "#15497B", "functional": "#7B2D8E"}
LABEL = {"structural": "Structural", "functional": "Functional (PCA)"}
MODS = ["structural", "functional"]
GREY = "#5A6570"
Q_LEVEL = 0.05

# The order the hypothesis states them in. Not sorted by result.
ORDER = [("income",   "Household income"),
         ("sleep",    "Sleep quality"),
         ("activity", "Physical activity\n(high vs low)"),
         ("pwb",      "Psychological\nwell-being"),
         ("soc",      "Social\nconnectedness"),
         ("stress",   "Perceived stress")]

STYLES = {
    "thesis": dict(ext="pdf", dpi=300, base=9.0, lw=1.5, ms=6, scale=1.0,
                   title_size=0, show_title=False),
    "presentation": dict(ext="png", dpi=200, base=15, lw=2.7, ms=11,
                         scale=1.55, title_size=20, show_title=True),
}
if os.environ.get("FIG_NO_TITLE"):
    for _st in STYLES.values():
        _st["show_title"] = False


def apply_style(st):
    plt.rcParams.update({
        "figure.dpi": st["dpi"],
        "font.size": st["base"],
        "axes.labelsize": st["base"],
        "xtick.labelsize": st["base"] * 0.95,
        "ytick.labelsize": st["base"],
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.spines.left": False,
        "axes.edgecolor": "#444444",
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.04,
        "figure.facecolor": "white",
    })


def bh(p):
    """Benjamini-Hochberg, the step-up form, identical to p.adjust(.,'BH')."""
    p = np.asarray(p, float)
    m = len(p)
    o = np.argsort(p)
    q = np.empty(m)
    q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


def load():
    f = os.path.join(RES, "sem", f"h2_single_estimates_{TAG}.csv")
    if not os.path.exists(f):
        raise SystemExit(f"not found: {f}\n  Run 05_h2_single.R first.")
    e = pd.read_csv(f)
    e = e[(e.is_test == True) & (e.type == "unstandardized_years")].copy()
    if len(e) != 12:
        print(f"  WARNING: expected 12 test rows, found {len(e)}")

    mp = os.path.join(RES, "sem", f"multiplicity_{TAG}.csv")
    if os.path.exists(mp):
        m = pd.read_csv(mp)
        m = m[m.family == "H2"][["modality", "label", "p_adjusted"]]
        e = e.merge(m.rename(columns={"label": "exposure",
                                      "p_adjusted": "q"}),
                    on=["modality", "exposure"], how="left")
        if e["q"].isna().any():
            raise SystemExit("multiplicity file does not cover every H2 row")
        print(f"  q values read from {os.path.basename(mp)}")
    else:
        e["q"] = bh(e.pvalue.values)
        print("  NOTE: multiplicity_%s.csv not found; q computed here.\n"
              "        Run 08_multiplicity.R and redraw so that the figure\n"
              "        and the decision table come from one place." % TAG)

    v = os.path.join(RES, "sem", f"h2_single_variance_{TAG}.csv")
    if os.path.exists(v):
        e = e.merge(pd.read_csv(v)[["modality", "exposure", "delta_r2"]],
                    on=["modality", "exposure"], how="left")
    else:
        e["delta_r2"] = np.nan
    return e


def main():
    e = load()
    os.makedirs(OUTD, exist_ok=True)

    for tag, st in STYLES.items():
        st = dict(st, tag=tag)
        apply_style(st)
        fig, ax = plt.subplots(figsize=(8.6 * st["scale"], 6.2 * st["scale"]))

        off = 0.19
        yt, yl = [], []
        for i, (key, lab) in enumerate(ORDER):
            yt.append(-i)
            yl.append(lab)
            for j, m in enumerate(MODS):
                r = e[(e.modality == m) & (e.exposure == key)]
                if not len(r):
                    continue
                r = r.iloc[0]
                y = -i + (off if j == 0 else -off)
                sig = r["q"] < Q_LEVEL

                ax.plot([r["ci.lower"], r["ci.upper"]], [y, y], color=COL[m],
                        lw=st["lw"] * (1.6 if sig else 1.1),
                        alpha=1.0 if sig else 0.55,
                        solid_capstyle="round", zorder=2)
                ax.plot([r["est"]], [y], "o",
                        color=COL[m] if sig else "white",
                        ms=st["ms"] * (1.0 if sig else 0.82),
                        markeredgecolor=COL[m],
                        markeredgewidth=st["lw"] * (0.6 if sig else 1.1),
                        zorder=3)
                if sig:
                    # "q = .000" claims a precision the number does not have.
                    qt = ("q < .001" if r["q"] < .001
                          else f"q = {r['q']:.3f}".replace("0.", "."))
                    ax.text(r["ci.upper"] + 0.07, y, qt,
                            va="center", ha="left",
                            fontsize=st["base"] * 0.9, color=COL[m],
                            weight="semibold")

        ax.axvline(0, color="#333333", lw=st["lw"], zorder=1)
        ax.set_yticks(yt)
        ax.set_yticklabels(yl)
        ax.tick_params(axis="y", length=0, pad=8)
        ax.set_xlabel("Years of brain-age gap per SD of the exposure  "
                      "(95% CI)")
        ax.set_ylim(-len(ORDER) + 0.45, 0.55)

        # Shade the one row that survives, rather than relying on the
        # marker alone to carry it.
        for i, (key, _) in enumerate(ORDER):
            if (e[e.exposure == key]["q"] < Q_LEVEL).any():
                ax.axhspan(-i - 0.42, -i + 0.42, color="#000000",
                           alpha=0.04, zorder=0)

        ax.annotate("filled: survives Benjamini-Hochberg at "
                    f"q = {Q_LEVEL:.2f} across the twelve tests".replace("0.", "."),
                    xy=(0.5, -0.145), xycoords="axes fraction",
                    ha="center", va="top",
                    fontsize=st["base"] * 0.88, color=GREY)

        for j, m in enumerate(MODS):
            fig.text(0.02 + j * 0.15, -0.055, LABEL[m], color=COL[m],
                     fontsize=st["base"], weight="semibold")

        if st["show_title"]:
            fig.suptitle("Household income was the only exposure associated "
                         "with the brain-age gap",
                         fontsize=st["title_size"], y=1.0, ha="center",
                         weight="semibold")
            fig.tight_layout(rect=[0, 0, 1, 0.94])
        else:
            fig.tight_layout()

        p = os.path.join(OUTD, f"fig12_h2_{tag}.{st['ext']}")
        fig.savefig(p)
        plt.close(fig)
        print("wrote " + os.path.relpath(p, ROOT))

    print("\nthe twelve tests, in hypothesis order:")
    for key, lab in ORDER:
        for m in MODS:
            r = e[(e.modality == m) & (e.exposure == key)]
            if not len(r):
                continue
            r = r.iloc[0]
            mark = "*" if r["q"] < Q_LEVEL else " "
            print(f" {mark} {m:11s} {key:9s} b = {r['est']:+.3f} yr  "
                  f"[{r['ci.lower']:+.2f}, {r['ci.upper']:+.2f}]  "
                  f"p = {r['pvalue']:.4f}  q = {r['q']:.4f}  "
                  f"dR2 = {r['delta_r2']:.4f}")


if __name__ == "__main__":
    main()
